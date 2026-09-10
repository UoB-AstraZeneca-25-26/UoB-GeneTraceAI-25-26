"""
expression_fusion/06_bulk_protein_scorer.py
--------------------------------------------
Bulk-computes lineage-conditioned protein z-scores for every gene using two
platforms (ProCAN DIA-MS + CCLE/Gygi TMT), then Stouffer-combines them into a
single prot_z_t per (gene, model_id).

OUTPUT
------
expression_fusion/outputs/bulk_prot_z.parquet
    columns: gene_id (lowercase ensg), model_id (lowercase ach-),
             prot_z_t (float32), n_sources (int8)
    one row per (gene, model_id) with a valid combined z-score.

HOW IT FEEDS INTO 02_core_score.ipynb
--------------------------------------
See expression_fusion/patch_core_score_orthogonal.py, which adds three cells:
  Cell 7h — loads bulk_prot_z.parquet → prot_z_wide
  Cell 7i — per-gene rho estimation + orthogonal combination → core_score

SCALE DECISIONS
---------------
- ProCAN  DIA-MS log2 intensity (median ~3.5, range -5 to 15). Already log scale.
- CCLE    TMT log-ratio centred at 0 (median ~0, range -11 to +5). Already log scale.
Neither source needs a log(x+1) transform. The within-lineage z-score removes
the ~3.5-unit offset between the two scales automatically (subtracts each
source's lineage median independently).

PLATFORM CORRELATION
--------------------
Per test_run_proteomics_platform_overlap.py:
  median per-protein Spearman rho = 0.373
  Kish n_eff = 2 / (1 + 0.373) = 1.45
The Stouffer shrinkage for a single-source (model_id, gene) pair:
  sqrt(N_SOURCES / 1) = sqrt(2) ~= 1.41  (score discounted by ~30%)
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

import common as C

N_SOURCES = 2  # ProCAN + CCLE


# ========================================================= UniProt -> ENSG map
def _build_uniprot_ensg(con) -> dict[str, str]:
    """
    Returns lowercase_uniprot -> lowercase_ensg mapping.
    Uses procan_protein_map and protein_map (CCLE) joined to the gene table
    via hugo_symbol. ProCAN entries take priority on collisions (DIA-MS preferred
    per the tier analysis in Cell 7f of 02_core_score.ipynb).
    """
    gene_df = con.execute("SELECT gene_id, hugo_symbol FROM gene").df()
    gene_df["hugo_lower"] = gene_df["hugo_symbol"].str.lower()
    sym_to_ensg = dict(zip(gene_df["hugo_lower"], gene_df["gene_id"]))

    ppm = con.execute("SELECT uniprot_id, gene_symbol FROM procan_protein_map").df()
    ppm["uniprot_id"]   = ppm["uniprot_id"].str.lower()
    ppm["gene_id"]      = ppm["gene_symbol"].str.lower().map(sym_to_ensg)
    ppm                 = ppm.dropna(subset=["gene_id"])

    pm = con.execute("SELECT uniprot_id, gene_symbol FROM protein_map").df()
    pm["uniprot_id"]    = pm["uniprot_id"].str.lower()
    pm["gene_id"]       = pm["gene_symbol"].str.lower().map(sym_to_ensg)
    pm                  = pm.dropna(subset=["gene_id"])

    # ProCAN first so it wins on collision
    combined = pd.concat([ppm[["uniprot_id", "gene_id"]],
                          pm[["uniprot_id", "gene_id"]]], ignore_index=True)
    combined = combined.drop_duplicates("uniprot_id", keep="first")
    mapping  = dict(zip(combined["uniprot_id"], combined["gene_id"]))

    n_ppm = ppm["uniprot_id"].nunique()
    n_pm  = pm["uniprot_id"].nunique()
    print(f"UniProt->ENSG: {n_ppm} ProCAN + {n_pm} CCLE -> {len(mapping)} unique entries")
    return mapping


# ============================================================== source loaders
def _load_procan(con, u2e: dict[str, str]) -> pd.DataFrame:
    """
    Loads ProCAN from procan_proteomics DB table.
    Returns wide DataFrame: model_id (index) x ensg_id (columns), log2 intensity.
    Rows with is_ambiguous=TRUE are excluded (same rule as GEO/HPA loaders).
    """
    print("Loading ProCAN from DB...", flush=True)
    t0 = time.time()

    META = {"gdsc_model_name", "sanger_model_id", "model_id",
            "matched_via", "n_model_id", "is_ambiguous"}

    pc = con.execute(
        "SELECT * FROM procan_proteomics "
        "WHERE model_id IS NOT NULL AND is_ambiguous = FALSE"
    ).df()
    data_cols = [c for c in pc.columns if c not in META]

    pc = pc.set_index("model_id")[data_cols]
    pc.index = pc.index.str.lower()
    pc = pc.groupby(level=0).mean()          # dedupe: take mean across duplicate model_ids

    # Map UniProt columns -> ENSG
    ensg_rename = {c: u2e[c] for c in pc.columns if c in u2e}
    pc = pc.rename(columns=ensg_rename)
    keep = [c for c in pc.columns if c in ensg_rename.values()]
    pc = pc[keep]
    pc = pc.T.groupby(level=0).mean().T      # average duplicate ENSG mappings

    print(f"  ProCAN: {pc.shape[0]} lines x {pc.shape[1]} genes  ({time.time()-t0:.0f}s)")
    return pc


def _load_ccle(con, u2e: dict[str, str]) -> pd.DataFrame:
    """
    Loads CCLE/Gygi TMT proteomics from the proteomics DB table.
    Returns wide DataFrame: model_id (index) x ensg_id (columns), log-ratio.
    """
    print("Loading CCLE proteomics from DB...", flush=True)
    t0 = time.time()

    prot = con.execute("SELECT * FROM proteomics").df()
    key_col = "model_id" if "model_id" in prot.columns else "depmap_id"
    prot = prot.set_index(key_col)
    prot.index = prot.index.str.lower()
    prot = prot.groupby(level=0).mean()

    ensg_rename = {c: u2e[c] for c in prot.columns if c in u2e}
    prot = prot.rename(columns=ensg_rename)
    keep = [c for c in prot.columns if c in ensg_rename.values()]
    prot = prot[keep]
    prot = prot.T.groupby(level=0).mean().T

    print(f"  CCLE:   {prot.shape[0]} lines x {prot.shape[1]} genes  ({time.time()-t0:.0f}s)")
    return prot


# ================================================================== main
def main():
    C.banner("Bulk Protein z-scorer: ProCAN + CCLE/Gygi -> prot_z_t")

    con          = C.connect()
    lineage_map  = C.load_lineage(con)
    print(f"Lineage map: {len(lineage_map):,} model_ids, "
          f"{lineage_map.nunique()} lineages\n")

    u2e          = _build_uniprot_ensg(con)
    procan_wide  = _load_procan(con, u2e)
    ccle_wide    = _load_ccle(con, u2e)

    all_genes = sorted(set(procan_wide.columns) | set(ccle_wide.columns))
    print(f"\nGene union: {len(all_genes):,} genes")

    # ---- per-source lineage z-scores ----
    print("\nComputing ProCAN z-scores per lineage...", flush=True)
    t0 = time.time()
    procan_genes = [g for g in all_genes if g in procan_wide.columns]
    procan_z     = C.score_source_lineage(procan_wide, procan_genes, lineage_map, "procan")
    print(f"  {len(procan_z):,} rows  ({time.time()-t0:.0f}s)")

    print("Computing CCLE z-scores per lineage...", flush=True)
    t0 = time.time()
    ccle_genes = [g for g in all_genes if g in ccle_wide.columns]
    ccle_z     = C.score_source_lineage(ccle_wide, ccle_genes, lineage_map, "ccle_prot")
    print(f"  {len(ccle_z):,} rows  ({time.time()-t0:.0f}s)")

    # ---- Stouffer combine ----
    print("\nStouffer combination...", flush=True)
    t0       = time.time()
    all_z    = pd.concat([procan_z, ccle_z], ignore_index=True)
    combined = C.stouffer_combine(all_z, n_sources_total=N_SOURCES)
    combined = combined.rename(columns={"z_t": "prot_z_t"})
    print(f"  {len(combined):,} (gene, model_id) pairs  ({time.time()-t0:.0f}s)")

    # ---- write ----
    out = C.OUT / "bulk_prot_z.parquet"
    combined.to_parquet(out, index=False)

    print(f"\nWrote {len(combined):,} rows -> {out}")
    print(f"  genes:     {combined.gene_id.nunique():,}")
    print(f"  model_ids: {combined.model_id.nunique():,}")
    print(f"  n_sources  1={int((combined.n_sources==1).sum()):,} "
          f"2={int((combined.n_sources==2).sum()):,}")
    print(f"  prot_z_t range: [{combined.prot_z_t.min():.2f}, "
          f"{combined.prot_z_t.max():.2f}]")

    con.close()


if __name__ == "__main__":
    main()
