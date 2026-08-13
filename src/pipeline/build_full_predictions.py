"""
Extends predictions_with_confidence.parquet from 141 curated genes
to all 19,176 scored genes. Reads core_score and flags in gene-batches
to keep memory low, writes parquet via pyarrow in streaming fashion.
"""
import duckdb
import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path

OUTPUTS = Path("src/pipeline/outputs")

# core_score.parquet is 28.4M rows / 225MB and is NOT clustered by ensg_id --
# all 28 row groups span the full ENSG range, so pyarrow `filters=` cannot prune
# and the previous `pd.read_parquet(...).query("ensg_id in @batch")` materialised
# the entire table into pandas once per batch (96 full reads). DuckDB streams the
# scan and applies the join without holding the table, which matters on a machine
# where this competes with the rest of the pipeline for RAM.
_con = duckdb.connect()

# Load lightweight tables in full
regime_curated = pd.read_parquet(OUTPUTS / "gene_regime.parquet")

# Expand regimes to COSMIC-mapped genes not in the curated 141
gl = pd.read_parquet("reference/gene_lookup.parquet")
gl["ensg_id"] = gl["ensg_id"].astype("string").str.split(".").str[0].str.lower()  # canonical case
role_to_class = {"oncogene": "abundance_tracking", "tsg": "loss_of_function", "both": "activation_driven"}
prior_rows = []
curated_ids = set(regime_curated["ensg_id"])
for _, gr in gl.iterrows():
    if gr["ensg_id"] in curated_ids:
        continue
    cls = role_to_class.get(gr["gene_role"])
    if cls:
        prior_rows.append({"ensg_id": gr["ensg_id"], "class": cls, "regime_source": "prior"})
regime_prior = pd.DataFrame(prior_rows)
regime = pd.concat([regime_curated, regime_prior], ignore_index=True)
regime_map = regime.set_index("ensg_id")[["class","regime_source"]].to_dict("index")
print(f"Regime table: {len(regime_curated)} measured + {len(regime_prior)} prior = {len(regime)} total")
tsg_genes = set(gl[gl["gene_role"] == "tsg"]["ensg_id"])

# Chronos-based validation check for unknown-class genes only (read-only diagnostic,
# NOT a scoring layer -- see build_chronos_validation.py for why: blending Chronos
# into core_score directly was tried and rejected, BCL2's hit@20 collapsed to 0).
# For genes with no curated regime and no COSMIC role, this tells us whether the
# abundance-based ranking actually tracks real CRISPR essentiality for that gene.
chronos_val = pd.read_parquet(OUTPUTS / "chronos_validation.parquet", columns=["ensg_id", "chronos_check"])
chronos_map = chronos_val.set_index("ensg_id")["chronos_check"].to_dict()
print(f"Chronos validation: {chronos_val['chronos_check'].value_counts().to_dict()}")

# Per-gene dynamic range (02_core_score.ipynb Cell 7c). Descriptive only -- it
# never enters core_score or the sort. It exists because core_score is a
# within-gene percentile, so a gene whose expression barely moves across cell
# lines still yields a full 0-1 spread and a top-10 that reads as confidently as
# a sharply tissue-restricted gene's. signal_spread='narrow' is the pipeline
# saying "this ordering is real but the differences behind it are small".
_disp_path = OUTPUTS / "gene_dispersion.parquet"
if not _disp_path.exists():
    raise SystemExit(
        f"{_disp_path} not found. Run 02_core_score.ipynb (Cell 7c) before this "
        "script -- the rebuild order is 02_core_score -> build_full_predictions.")
dispersion = pd.read_parquet(
    _disp_path, columns=["ensg_id", "dynamic_range", "top_vs_median_fold", "signal_spread"])
dispersion = dispersion.drop_duplicates("ensg_id").set_index("ensg_id")
print(f"Gene dispersion: {len(dispersion):,} genes, "
      f"bands {dispersion['signal_spread'].value_counts().to_dict()}")

# Gene list only from core (no data yet)
gene_list = pd.read_parquet(
    OUTPUTS / "core_score.parquet", columns=["ensg_id"]
)["ensg_id"].unique()
print(f"{len(gene_list):,} genes to process")

OUT_COLS = ["model_id","ensg_id","core_score","n_layers","stratum_rank",
            "class","regime_source","rank_basis","has_driver_alteration",
            "has_cna_alteration","confidence","confidence_reason","gene_role",
            "signal_quality","chronos_check",
            "gene_dynamic_range","gene_top_vs_median_fold","signal_spread"]

writer = None
BATCH = 2000  # genes per batch (was 200; each batch cost a full 225MB read)

# flags_with_driver is only 886k rows -- load once, slice per batch, rather than
# re-reading and re-filtering it inside the loop.
flags_all = pd.read_parquet(
    OUTPUTS / "flags_with_driver.parquet",
    columns=["model_id","ensg_id","has_driver_alteration","has_cna_alteration"]
)
flags_all["model_id"] = flags_all["model_id"].str.lower()
flags_all = flags_all.set_index("ensg_id", drop=False).sort_index()

CORE_PATH = str((OUTPUTS / "core_score.parquet").as_posix())

for batch_start in range(0, len(gene_list), BATCH):
    batch_genes = gene_list[batch_start:batch_start+BATCH]

    # Streamed, memory-bounded read of just this batch's genes
    _con.register("batch_genes", pd.DataFrame({"ensg_id": batch_genes}))
    core_b = _con.execute(f"""
        SELECT c.model_id, c.ensg_id, c.core_score, c.n_layers, c.stratum_rank
        FROM read_parquet('{CORE_PATH}') c
        JOIN batch_genes USING (ensg_id)
    """).df()
    _con.unregister("batch_genes")
    core_b["model_id"] = core_b["model_id"].str.lower()

    # Invert scores for TSGs: low abundance = high relevance
    is_tsg = core_b["ensg_id"].isin(tsg_genes)
    core_b.loc[is_tsg, "core_score"]   = 1.0 - core_b.loc[is_tsg, "core_score"]
    core_b.loc[is_tsg, "stratum_rank"] = 1.0 - core_b.loc[is_tsg, "stratum_rank"]

    present = flags_all.index.intersection(pd.Index(batch_genes).unique())
    flags_b = flags_all.loc[present].reset_index(drop=True)

    c = core_b.merge(flags_b, on=["model_id","ensg_id"], how="left")
    assert len(c) == len(core_b), (
        f"fan-out on flags merge: {len(core_b):,} -> {len(c):,} rows "
        f"(duplicate (model_id, ensg_id) in flags_with_driver)")
    c["has_driver_alteration"] = c["has_driver_alteration"].fillna(False)
    c["has_cna_alteration"]    = c["has_cna_alteration"].fillna(False)

    c["class"]         = c["ensg_id"].map(lambda g: regime_map.get(g, {}).get("class", "unknown"))
    c["regime_source"] = c["ensg_id"].map(lambda g: regime_map.get(g, {}).get("regime_source", "unknown"))

    # Chronos essentiality check -- only meaningful for unknown-class genes (curated
    # and COSMIC-prior genes already have their own validation basis). "untested"
    # means no CRISPR essentiality coverage exists for this gene at all.
    c["chronos_check"] = np.where(
        c["class"] == "unknown",
        c["ensg_id"].map(chronos_map).fillna("untested"),
        "n/a"
    )

    # rank_basis must describe the sort that explain_pair.sort_gene_rows actually
    # performs — see DRIVER_GATED_CLASSES there. 'unknown' genes have no census
    # role and no measured regime, so their driver flag only breaks ties; calling
    # that 'unknown_dual' implied a gate that is no longer applied.
    c["rank_basis"] = np.where(
        c["class"] == "abundance_tracking", "core_score",
        np.where(c["class"] == "loss_of_function",
            np.where(c["has_driver_alteration"], "lof_flag+inverted", "lof_inverted"),
        np.where(c["class"] == "unknown", "score+driver_tiebreak",
        np.where(c["has_driver_alteration"], "driver_flag+score", "score_only"))))

    tier = np.full(len(c), "moderate", dtype=object)

    # Base high: abundance_tracking, empirically validated, 2-layer expression+proteomics
    tier = np.where(
        (c.regime_source == "measured") & (c["class"] == "abundance_tracking") & (c.n_layers == 2),
        "high", tier)

    tier = np.where(
        (c["class"] == "activation_driven") & (c.rank_basis == "score_only"),
        "low", tier)
    tier = np.where((c.regime_source == "prior") & ~c["has_cna_alteration"], "prior", tier)
    tier = np.where(c.regime_source == "unknown", "unknown", tier)

    # CNA upgrades applied LAST so they override prior/unknown/low for confirmed dependencies
    # Deletion in TSG/loss_of_function = direct mechanistic evidence
    tier = np.where(
        (c["class"] == "loss_of_function") & c["has_cna_alteration"],
        "high", tier)
    # Amplification in oncogene/activation_driven = direct mechanistic evidence
    tier = np.where(
        (c["class"].isin(["activation_driven", "abundance_tracking"])) & c["has_cna_alteration"],
        "high", tier)
    # Prior-only genes (unknown class, e.g. oncogene not in curated 141) with CNA → moderate
    # (CNA confirms the hypothesis but no empirical validation exists)
    tier = np.where(
        (c.regime_source == "prior") & c["has_cna_alteration"] &
        ~c["class"].isin(["loss_of_function", "activation_driven", "abundance_tracking"]),
        "moderate", tier)

    # Chronos-validated upgrade: unknown-class genes where the abundance ranking is
    # independently confirmed against real CRISPR essentiality data (see
    # build_chronos_validation.py). Does NOT touch "inverted" or "none" genes --
    # those stay "unknown", the chronos_check column carries that information instead.
    tier = np.where(
        (c["class"] == "unknown") & (c["chronos_check"] == "validated"),
        "moderate", tier)
    c["confidence"] = tier

    c["confidence_reason"] = np.where(
        (c.confidence == "high") & c["has_cna_alteration"] & (c["class"] == "loss_of_function"),
        "High: loss_of_function gene with confirmed deletion (CNA) in this line — "
        "direct mechanistic evidence of TSG dependency.",
        np.where(
            (c.confidence == "high") & c["has_cna_alteration"],
            "High: " + c["class"] + " gene with confirmed amplification (CNA) in this line — "
            "direct mechanistic evidence of oncogene dependency.",
        np.where(
            c.confidence == "high",
            "High: abundance_tracking gene (empirically validated), scored on "
                + c.n_layers.astype(str) + " layers.",
        np.where(
            (c.confidence == "moderate") & (c["chronos_check"] == "validated"),
            "Moderate: unknown-class gene whose abundance ranking independently "
            "correlates with real CRISPR essentiality data across cell lines — "
            "not GDSC-validated, but confirmed against a different, direct "
            "functional ground truth.",
        np.where(
            c.confidence == "low",
            "Low: activation-driven gene dependency; ranked on score alone "
            "with no supporting driver alteration in this line.",
        np.where(
            (c.confidence == "unknown") & (c["chronos_check"] == "inverted"),
            "Unknown: no GDSC validation data for this gene, AND its abundance "
            "ranking correlates NEGATIVELY with real CRISPR essentiality — "
            "ranking direction for this gene is likely backwards, treat with "
            "particular caution.",
        np.where(
            (c.confidence == "unknown") & c["chronos_check"].isin(["weak_positive", "weak_negative"]),
            "Unknown: no GDSC validation data for this gene; tested against CRISPR "
            "essentiality data and found a statistically significant but very small "
            "relationship (|rho| < 0.30) — direction recorded, too weak to treat as "
            "validation either way.",
        np.where(
            (c.confidence == "unknown") & (c["chronos_check"] == "none"),
            "Unknown: no GDSC validation data for this gene; tested against "
            "CRISPR essentiality data and found no significant relationship — "
            "a checked unknown, not an assumed one.",
        np.where(
            c.confidence == "unknown",
            "Unknown: no validation data for this gene; scoring regime "
            "unclassified, ranking shown but not regime-weighted.",
        np.where(
            c.confidence == "prior",
            "Prior: regime assigned from COSMIC gene role (" + c["class"] + "); "
            "not empirically validated in DepMap — treat as hypothesis.",
        np.where(
            c.rank_basis == "driver_flag+score",
            "Moderate: activation-driven gene; ranked on driver alteration "
            "present in this line, with score as support.",
            "Moderate: " + c["class"] + " gene, ranked on "
                + c.rank_basis + ", " + c.n_layers.astype(str) + "-layer coverage."
        )))))))))))


    role_lookup = gl.set_index("ensg_id")["gene_role"]
    c["gene_role"] = c["ensg_id"].map(role_lookup).fillna("unknown")

    # Gene-level dispersion, broadcast per row (same value for every cell line of
    # a gene). 'unmeasured' where the gene has no RNA layer to measure spread on.
    c["gene_dynamic_range"]      = c["ensg_id"].map(dispersion["dynamic_range"])
    c["gene_top_vs_median_fold"] = c["ensg_id"].map(dispersion["top_vs_median_fold"])
    c["signal_spread"]           = c["ensg_id"].map(dispersion["signal_spread"]).fillna("unmeasured")

    c["signal_quality"] = np.select(
        [c["confidence"] == "high", c["regime_source"] == "measured",
         c["chronos_check"] == "validated", c["regime_source"] == "prior"],
        ["evidence-based", "evidence-based", "evidence-based", "prior-informed"],
        default="abundance-only",
    )

    table = pa.Table.from_pandas(c[OUT_COLS], preserve_index=False)
    if writer is None:
        writer = pq.ParquetWriter(OUTPUTS / "predictions_with_confidence.parquet", table.schema)
    writer.write_table(table)

    done = min(batch_start + BATCH, len(gene_list))
    print(f"  {done:>6}/{len(gene_list):,} genes written ({len(c):,} rows)", flush=True)

if writer:
    writer.close()

print("Verifying...")
pred = pd.read_parquet(OUTPUTS / "predictions_with_confidence.parquet")
print(f"Shape: {pred.shape}  |  genes: {pred.ensg_id.nunique():,}")
print("Confidence dist:", pred.confidence.value_counts().to_dict())
shank3 = pred[pred.ensg_id == "ensg00000251322"]
print(f"SHANK3 rows: {len(shank3)}")
if not shank3.empty:
    print(shank3[["model_id","core_score","class","confidence","rank_basis"]].head(3).to_string(index=False))
print("Done.")
