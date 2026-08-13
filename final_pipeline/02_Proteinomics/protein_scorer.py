"""
02_Proteinomics/protein_scorer.py
-----------------------------------
Lineage-conditioned protein z-scores from two platforms:
  - ProCAN  DIA-MS  (procan_proteomics: wide, model_id pre-joined)
  - CCLE/Gygi TMT   (proteomics: wide, model_id column)

Non-protein columns in procan_proteomics:
  gdsc_model_name, sanger_model_id, model_id, matched_via, n_model_id, is_ambiguous

Platform combination rules (from platform_tier.py):
  consistent or cautious tier  -> Stouffer mean of both z-scores (shrinkage for 1-source)
  conflicting tier             -> ProCAN only

UniProt → ENSG mapping via procan_protein_map + protein_map.

Reads from:  celllineselector.db
             final_pipeline/outputs/protein_platform_tier.parquet
Writes to:   final_pipeline/outputs/bulk_prot_z.parquet
Schema: gene_id (ENSG uppercase), model_id (ACH lowercase), z_t (f32), n_sources (int8)
"""
import sys, time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import PROT_Z, PROT_TIER
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
import common as C

N_SOURCES    = 2
_PROCAN_META = {"gdsc_model_name", "sanger_model_id", "model_id",
                "matched_via", "n_model_id", "is_ambiguous"}


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


def run():
    t0 = time.time()
    print("=" * 70)
    print("02_Proteinomics — protein z-scorer (ProCAN + CCLE)")
    print("=" * 70)

    if not PROT_TIER.exists():
        raise SystemExit(f"Run platform_tier.py first — {PROT_TIER} missing.")

    tier_df = pd.read_parquet(PROT_TIER)
    conflicting = set(tier_df.loc[tier_df.platform_tier == "conflicting", "uniprot"].str.lower())

    con = C.connect()
    u2e = _build_uniprot_ensg(con)
    lineage_map = C.load_lineage(con)

    procan_df = con.execute("SELECT * FROM main.procan_proteomics").df()
    ccle_df   = con.execute("SELECT * FROM main.proteomics").df()
    con.close()

    procan_long = _wide_to_long(procan_df, "model_id", _PROCAN_META)
    procan_long["gene_id"] = procan_long["uniprot"].map(u2e)
    procan_long = procan_long.dropna(subset=["gene_id"])

    ccle_long = _wide_to_long(ccle_df, "model_id", {"model_id"})
    ccle_long["gene_id"] = ccle_long["uniprot"].map(u2e)
    ccle_long = ccle_long.dropna(subset=["gene_id"])

    print(f"ProCAN: {procan_long.gene_id.nunique():,} genes  |  "
          f"CCLE: {ccle_long.gene_id.nunique():,} genes")

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
    shrink = np.where(agg["n_sources"] < N_SOURCES,
                      np.sqrt(N_SOURCES / agg["n_sources"]), 1.0)
    agg["z_t"]       = (z_raw / shrink).astype("float32")
    agg["n_sources"] = agg["n_sources"].astype("int8")

    result = agg[["gene_id","model_id","z_t","n_sources"]]
    result.to_parquet(PROT_Z, index=False)

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s  |  rows: {len(result):,}  ->  {PROT_Z}")
    print(f"  genes: {result.gene_id.nunique():,}  |  lines: {result.model_id.nunique():,}")


if __name__ == "__main__":
    run()
