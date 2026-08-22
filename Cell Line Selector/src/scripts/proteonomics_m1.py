"""
02_Proteinomics/protein.py
---------------------------
Merged pipeline: platform tiering + lineage-conditioned protein z-scores.

STAGE 1 (platform_tier.py logic)
  Computes per-protein Spearman rho between ProCAN and CCLE on their shared
  cell lines, then assigns each protein to a tier:
    consistent  — rho >= 0.5
    cautious    — 0.3 <= rho < 0.5
    conflicting — rho < 0.3   (use ProCAN only downstream)

STAGE 2 (protein_scorer.py logic)
  Lineage-conditioned protein z-scores from two platforms:
    - ProCAN  DIA-MS  (procan_proteomics: wide, model_id pre-joined)
    - CCLE/Gygi TMT   (proteomics: wide, model_id column)
  Platform combination rules (from Stage 1 tiers):
    consistent or cautious tier  -> Stouffer mean of both z-scores
    conflicting tier             -> ProCAN only
  UniProt → ENSG mapping via procan_protein_map + protein_map.

  FIX vs original protein_scorer.py: uniprot values are lowercased before
  the .map(u2e) lookup on BOTH platforms. u2e's keys are built lowercase
  (_build_uniprot_ensg does str(...).lower()), but the wide tables' UniProt
  column headers are uppercase accessions (e.g. P04637). Mapping uppercase
  strings against a lowercase-keyed dict silently returned NaN for most
  rows, which were then dropped by dropna(subset=["gene_id"]) -- so the two
  platforms ended up on almost disjoint gene_id sets and effectively never
  combined (n_sources stayed at 1). Lowercasing both sides before the map
  fixes the join key mismatch.

Both tables in DuckDB are WIDE (cell lines as rows, UniProt IDs as columns).
procan_proteomics also has model_id (ACH) pre-joined at the end of the table.
Non-protein identifier columns are: gdsc_model_name, sanger_model_id,
model_id, matched_via, n_model_id, is_ambiguous.

Reads from:  celllineselector.db
Writes to:   final_pipeline/outputs/protein_platform_tier.parquet
             final_pipeline/outputs/bulk_prot_z.parquet
Schemas:
  protein_platform_tier.parquet -> uniprot, platform_rho, platform_tier
  bulk_prot_z.parquet           -> gene_id (ENSG uppercase), model_id (ACH lowercase), z_t (f32), n_sources (int8)
"""
import sys, time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import PROT_Z, PROT_TIER
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
import common as C

TIER_CUTOFFS = {"consistent": 0.5, "cautious": 0.3}
MIN_SHARED   = 30

_PROCAN_ID_COLS = {"gdsc_model_name", "sanger_model_id", "model_id",
                   "matched_via", "n_model_id", "is_ambiguous"}


# ─────────────────────────────────────────────────────────────────────────
# Stage 1 — platform tier
# ─────────────────────────────────────────────────────────────────────────
def _assign_tier(rho: float) -> str:
    if rho >= TIER_CUTOFFS["consistent"]:
        return "consistent"
    if rho >= TIER_CUTOFFS["cautious"]:
        return "cautious"
    return "conflicting"


def _wide_to_long_procan(con) -> pd.DataFrame:
    df = con.execute("SELECT * FROM main.procan_proteomics").df()
    prot_cols = [c for c in df.columns if c not in _PROCAN_ID_COLS]
    df = df[["model_id"] + prot_cols].set_index("model_id")
    try:
        long = df.stack(future_stack=True).reset_index()
    except TypeError:
        long = df.stack().reset_index()
    long.columns = ["model_id", "uniprot", "val_procan"]
    return long.dropna(subset=["val_procan"])


def _wide_to_long_ccle(con) -> pd.DataFrame:
    df = con.execute("SELECT * FROM main.proteomics").df()
    prot_cols = [c for c in df.columns if c != "model_id"]
    df = df[["model_id"] + prot_cols].set_index("model_id")
    try:
        long = df.stack(future_stack=True).reset_index()
    except TypeError:
        long = df.stack().reset_index()
    long.columns = ["model_id", "uniprot", "val_ccle"]
    return long.dropna(subset=["val_ccle"])


def compute_platform_tier() -> pd.DataFrame:
    """Exact logic of platform_tier.py's run(). Writes PROT_TIER and
    returns tier_df so Stage 2 can consume it directly, in-memory."""
    t0 = time.time()
    print("=" * 70)
    print("02_Proteinomics — protein platform tier (ProCAN vs CCLE)")
    print("=" * 70)

    con = C.connect()
    procan_long = _wide_to_long_procan(con)
    print(f"ProCAN long: {len(procan_long):,} rows  |  "
          f"proteins: {procan_long.uniprot.nunique():,}  |  "
          f"lines: {procan_long.model_id.nunique():,}")

    ccle_long = _wide_to_long_ccle(con)
    print(f"CCLE long:   {len(ccle_long):,} rows  |  "
          f"proteins: {ccle_long.uniprot.nunique():,}  |  "
          f"lines: {ccle_long.model_id.nunique():,}")
    con.close()

    merged = procan_long.merge(ccle_long, on=["model_id", "uniprot"], how="inner")
    eligible = merged.groupby("uniprot").size()
    eligible = eligible[eligible >= MIN_SHARED].index
    merged = merged[merged.uniprot.isin(eligible)]
    print(f"\nShared pairs: {len(merged):,}  |  "
          f"eligible proteins (>={MIN_SHARED} shared lines): {len(eligible):,}")

    results = []
    for uniprot, grp in merged.groupby("uniprot"):
        rho, _ = spearmanr(grp["val_procan"], grp["val_ccle"])
        if np.isfinite(rho):
            results.append({"uniprot": uniprot, "platform_rho": rho,
                            "platform_tier": _assign_tier(rho)})

    tier_df = pd.DataFrame(results)
    tier_df.to_parquet(PROT_TIER, index=False)

    counts = tier_df["platform_tier"].value_counts()
    print(f"\nTier distribution ({len(tier_df):,} proteins):")
    for tier, n in counts.items():
        print(f"  {tier:<15s} {n:>5,}  ({n/len(tier_df)*100:.1f}%)")
    print(f"Median rho: {tier_df.platform_rho.median():.3f}")
    print(f"Written: {PROT_TIER}  ({time.time()-t0:.0f}s)")

    return tier_df


# ─────────────────────────────────────────────────────────────────────────
# Stage 2 — protein z-scorer
# ─────────────────────────────────────────────────────────────────────────
def _build_uniprot_ensg(con) -> dict[str, str]:
    """Map lowercase uniprot → uppercase ENSG."""
    gene_df = con.execute("SELECT gene_id, hugo_symbol FROM gene").df()
    sym_map = dict(zip(gene_df["hugo_symbol"].str.lower(), gene_df["gene_id"].str.upper()))

    ppm = con.execute("SELECT uniprot_id, gene_symbol FROM procan_protein_map").df()
    pm  = con.execute("SELECT uniprot_id, gene_symbol FROM protein_map").df()

    out = {}
    for df in [pm, ppm]:           # ProCAN last → wins on collision
        for _, row in df.iterrows():
            ensg = sym_map.get(str(row["gene_symbol"]).lower())
            if ensg:
                out[str(row["uniprot_id"]).lower()] = ensg
    return out


def _wide_to_long(df: pd.DataFrame, id_col: str, meta_cols: set) -> pd.DataFrame:
    prot_cols = [c for c in df.columns if c not in meta_cols and c != id_col]
    sub = df[[id_col] + prot_cols].set_index(id_col)
    try:
        long = sub.stack(future_stack=True).reset_index()
    except TypeError:
        long = sub.stack().reset_index()
    long.columns = [id_col, "uniprot", "value"]
    return long.dropna(subset=["value"])


def _zscore_platform(df: pd.DataFrame, lineage_map: pd.Series, src: str) -> pd.DataFrame:
    wide = df.pivot_table(index="model_id", columns="gene_id", values="value", aggfunc="first")
    wide["lineage"] = lineage_map.reindex(wide.index)
    wide = wide.dropna(subset=["lineage"])
    rows = []
    for _, grp in wide.groupby("lineage"):
        gene_cols = [c for c in grp.columns if c != "lineage"]
        X  = grp[gene_cols].values.astype(float)
        Z  = C.robust_z_matrix(X)
        zw = pd.DataFrame(Z, index=grp.index.tolist(), columns=gene_cols)
        try:
            zl = zw.stack(future_stack=True).dropna().reset_index()
        except TypeError:
            zl = zw.stack().dropna().reset_index()
        zl.columns = ["model_id", "gene_id", "z"]
        rows.append(zl)
    if not rows:
        return pd.DataFrame(columns=["model_id", "gene_id", "z", "source"])
    out = pd.concat(rows, ignore_index=True)
    out["source"] = src
    return out


def compute_protein_scores(tier_df: pd.DataFrame) -> pd.DataFrame:
    """Exact logic of protein_scorer.py's run(), taking tier_df straight
    from Stage 1 in-memory, with the uniprot case-mismatch fix applied."""
    t0 = time.time()
    print("=" * 70)
    print("02_Proteinomics — protein z-scorer (ProCAN + CCLE)")
    print("=" * 70)

    conflicting = set(tier_df.loc[tier_df.platform_tier == "conflicting", "uniprot"].str.lower())

    con = C.connect()
    u2e = _build_uniprot_ensg(con)
    lineage_map = C.load_lineage(con)

    procan_df = con.execute("SELECT * FROM main.procan_proteomics").df()
    ccle_df   = con.execute("SELECT * FROM main.proteomics").df()
    con.close()

    procan_long = _wide_to_long(procan_df, "model_id", _PROCAN_ID_COLS)
    procan_long["gene_id"] = procan_long["uniprot"].str.lower().map(u2e)   # fixed: lowercase before map
    procan_long = procan_long.dropna(subset=["gene_id"])

    ccle_long = _wide_to_long(ccle_df, "model_id", {"model_id"})
    ccle_long["gene_id"] = ccle_long["uniprot"].str.lower().map(u2e)       # fixed: lowercase before map
    ccle_long = ccle_long.dropna(subset=["gene_id"])

    print(f"ProCAN: {procan_long.gene_id.nunique():,} genes  |  "
          f"CCLE: {ccle_long.gene_id.nunique():,} genes")

    shared_genes = set(procan_long.gene_id) & set(ccle_long.gene_id)
    print(f"Genes shared by both platforms after mapping: {len(shared_genes):,}")

    conflicting_gene_ids = {u2e[u] for u in conflicting if u in u2e}
    ccle_long = ccle_long[~ccle_long["gene_id"].isin(conflicting_gene_ids)]

    z_procan = _zscore_platform(procan_long, lineage_map, "procan")
    z_ccle   = _zscore_platform(ccle_long,   lineage_map, "ccle")
    all_z = pd.concat([z_procan, z_ccle], ignore_index=True)

    src_n = all_z.groupby(["gene_id","source"])["model_id"].nunique().reset_index(name="src_n")
    all_z = all_z.merge(src_n, on=["gene_id","source"])
    all_z["w"]  = np.sqrt(all_z["src_n"])
    all_z["wz"] = all_z["w"] * all_z["z"]
    all_z["w2"] = all_z["w"] ** 2

    agg = all_z.groupby(["gene_id","model_id"]).agg(
        wz_sum=("wz","sum"), w2_sum=("w2","sum"), n_sources=("source","nunique")
    ).reset_index()

    denom  = np.sqrt(agg["w2_sum"].values)
    safe   = denom > 0
    z_raw  = np.where(safe, agg["wz_sum"].values / np.where(safe, denom, 1.0), np.nan)
    agg["z_t"]       = z_raw.astype("float32")
    agg["n_sources"] = agg["n_sources"].astype("int8")

    result = agg[["gene_id","model_id","z_t","n_sources"]]
    result.to_parquet(PROT_Z, index=False)

    n_dual = (agg["n_sources"] == 2).sum()
    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s  |  rows: {len(result):,}  ->  {PROT_Z}")
    print(f"  genes: {result.gene_id.nunique():,}  |  lines: {result.model_id.nunique():,}  |  "
          f"rows with both platforms (n_sources=2): {n_dual:,}")

    return result


# ─────────────────────────────────────────────────────────────────────────
def run():
    tier_df = compute_platform_tier()
    print()
    compute_protein_scores(tier_df)


if __name__ == "__main__":
    run()