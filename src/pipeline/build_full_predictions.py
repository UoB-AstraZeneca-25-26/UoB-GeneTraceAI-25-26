"""
Extends predictions_with_confidence.parquet from 141 curated genes
to all 19,176 scored genes. Reads core_score and flags in gene-batches
to keep memory low, writes parquet via pyarrow in streaming fashion.
"""
import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path

OUTPUTS = Path("src/pipeline/outputs")

# Load lightweight tables in full
regime_curated = pd.read_parquet(OUTPUTS / "gene_regime.parquet")

# Expand regimes to COSMIC-mapped genes not in the curated 141
gl = pd.read_parquet("reference/gene_lookup.parquet")
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

# Gene list only from core (no data yet)
gene_list = pd.read_parquet(
    OUTPUTS / "core_score.parquet", columns=["ensg_id"]
)["ensg_id"].unique()
print(f"{len(gene_list):,} genes to process")

OUT_COLS = ["model_id","ensg_id","core_score","n_layers","stratum_rank",
            "class","regime_source","rank_basis","has_driver_alteration",
            "has_cna_alteration","confidence","confidence_reason","gene_role",
            "signal_quality","chronos_check"]

writer = None
BATCH = 200  # genes per batch

for batch_start in range(0, len(gene_list), BATCH):
    batch_genes = gene_list[batch_start:batch_start+BATCH]

    # Read only this batch's rows from both parquets
    core_b = pd.read_parquet(OUTPUTS / "core_score.parquet").query("ensg_id in @batch_genes")
    core_b["model_id"] = core_b["model_id"].str.upper()

    # Invert scores for TSGs: low abundance = high relevance
    is_tsg = core_b["ensg_id"].isin(tsg_genes)
    core_b.loc[is_tsg, "core_score"]   = 1.0 - core_b.loc[is_tsg, "core_score"]
    core_b.loc[is_tsg, "stratum_rank"] = 1.0 - core_b.loc[is_tsg, "stratum_rank"]

    flags_b = pd.read_parquet(
        OUTPUTS / "flags_with_driver.parquet",
        columns=["model_id","ensg_id","has_driver_alteration","has_cna_alteration"]
    ).query("ensg_id in @batch_genes")
    flags_b["model_id"] = flags_b["model_id"].str.upper()

    c = core_b.merge(flags_b, on=["model_id","ensg_id"], how="left")
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

    c["rank_basis"] = np.where(
        c["class"] == "abundance_tracking", "core_score",
        np.where(c["class"] == "loss_of_function",
            np.where(c["has_driver_alteration"], "lof_flag+inverted", "lof_inverted"),
        np.where(c["class"] == "unknown", "unknown_dual",
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
        ))))))))))


    role_lookup = gl.set_index("ensg_id")["gene_role"]
    c["gene_role"] = c["ensg_id"].map(role_lookup).fillna("unknown")

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
    if (batch_start // BATCH) % 5 == 0:
        print(f"  {done:>6}/{len(gene_list):,} genes written", flush=True)

if writer:
    writer.close()

print("Verifying...")
pred = pd.read_parquet(OUTPUTS / "predictions_with_confidence.parquet")
print(f"Shape: {pred.shape}  |  genes: {pred.ensg_id.nunique():,}")
print("Confidence dist:", pred.confidence.value_counts().to_dict())
shank3 = pred[pred.ensg_id == "ENSG00000251322"]
print(f"SHANK3 rows: {len(shank3)}")
if not shank3.empty:
    print(shank3[["model_id","core_score","class","confidence","rank_basis"]].head(3).to_string(index=False))
print("Done.")
