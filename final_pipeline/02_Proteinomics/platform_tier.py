"""
02_Proteinomics/platform_tier.py
---------------------------------
Computes per-protein Spearman rho between ProCAN and CCLE on their shared
cell lines, then assigns each protein to a tier:
  consistent  — rho >= 0.5
  cautious    — 0.3 <= rho < 0.5
  conflicting — rho < 0.3   (use ProCAN only downstream)

Both tables in DuckDB are WIDE (cell lines as rows, UniProt IDs as columns).
procan_proteomics also has model_id (ACH) pre-joined at the end of the table.
Non-protein identifier columns are: gdsc_model_name, sanger_model_id,
model_id, matched_via, n_model_id, is_ambiguous.

Reads from:  celllineselector.db
Writes to:   final_pipeline/outputs/protein_platform_tier.parquet
Schema: uniprot, platform_rho, platform_tier
"""
import sys, time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import PROT_TIER
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
import common as C

TIER_CUTOFFS = {"consistent": 0.5, "cautious": 0.3}
MIN_SHARED   = 30

_PROCAN_ID_COLS = {"gdsc_model_name", "sanger_model_id", "model_id",
                   "matched_via", "n_model_id", "is_ambiguous"}


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


def run():
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


if __name__ == "__main__":
    run()
