"""
02_Proteinomics/export_protein_by_source.py
----------------------------------------------
Standalone export, NOT part of run_all.py: recovers the per-source protein
z-scores (ProCAN / CCLE) that 06_protein_score_M2.py's compute_from_sources()
already computes internally as its `all_z` step, before combining them into
one z_score -- that combination is what gets persisted to protein_z_m2_run_all
/ bulk_prot_z.parquet today; the two inputs behind it are discarded.

Does not modify 06_protein_score_M2.py -- imports its existing, already-tested
functions (load_sources, collapse_isoforms, zscore_platform) and runs the same
first three steps of compute_from_sources() by hand, using the exact same
finalized thresholds (results/m2_thresholds.json) run_protein_m2_stage.py
uses, so results are reproducible with the live pipeline. Stops before the
platform-agreement combination step (not needed for per-source values).

Writes to: final_pipeline/outputs/bulk_prot_z_by_source.parquet
Schema:    gene_id, model_id, z_procan, z_ccle (either may be null -- that
           source didn't measure this gene/line)
"""
import json
import sys
from pathlib import Path

import duckdb
import pandas as pd

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from config import OUT

DB_PATH = PROJECT_ROOT / "outputs" / "celllineselector.db"
THRESHOLDS_PATH = HERE / "results" / "m2_thresholds.json"
OUT_PATH = OUT / "bulk_prot_z_by_source.parquet"


def run():
    import importlib
    m2 = importlib.import_module("06_protein_score_M2")

    thr = json.loads(THRESHOLDS_PATH.read_text())

    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        protein_map = {**m2.build_symbol_to_ensg(con), **m2.build_uniprot_to_ensg(con)}
        lineage_map = m2.load_lineage_map(con)
        raw_longs = m2.load_sources(con, [m2.CCLE_SOURCE, m2.PROCAN_SOURCE], protein_map)
    finally:
        con.close()

    if not raw_longs:
        raise SystemExit("No protein source tables could be loaded.")

    # Same as compute_from_sources() steps 1 + 3: isoform collapse, then
    # lineage-conditioned robust z, per source -- stops before the platform
    # combination step, since that's exactly the part that discards the
    # per-source breakdown we want here.
    blocks = []
    for name, raw in raw_longs.items():
        collapsed, _qc = m2.collapse_isoforms(
            raw, name, min_overlap=thr["isoform_min_overlap"], corr_threshold=thr["isoform_corr"]
        )
        if collapsed.empty:
            continue
        z = m2.zscore_platform(collapsed, lineage_map, name)
        if not z.empty:
            blocks.append(z)

    if not blocks:
        raise SystemExit("No protein rows could be scored.")

    all_z = pd.concat(blocks, ignore_index=True)
    all_z["gene_id"] = all_z["gene_id"].astype(str).str.upper()
    all_z["model_id"] = all_z["model_id"].astype(str).str.upper()

    wide = all_z.pivot_table(index=["gene_id", "model_id"], columns="source",
                             values="z_score", aggfunc="first").reset_index()
    wide = wide.rename(columns={m2.PROCAN_SOURCE: "z_procan", m2.CCLE_SOURCE: "z_ccle"})
    for col in ["z_procan", "z_ccle"]:
        if col not in wide.columns:
            wide[col] = pd.NA
    wide = wide[["gene_id", "model_id", "z_procan", "z_ccle"]]

    wide.to_parquet(OUT_PATH, index=False)
    print(f"Exported {len(wide):,} rows -> {OUT_PATH}")
    print(f"  genes: {wide.gene_id.nunique():,}  |  lines: {wide.model_id.nunique():,}")
    print(f"  z_procan: {wide['z_procan'].notna().sum():,} non-null")
    print(f"  z_ccle: {wide['z_ccle'].notna().sum():,} non-null")


if __name__ == "__main__":
    run()
