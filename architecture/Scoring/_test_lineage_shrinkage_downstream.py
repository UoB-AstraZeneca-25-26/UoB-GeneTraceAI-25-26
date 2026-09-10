"""
Scoring/_test_lineage_shrinkage_downstream.py  (TEST HARNESS -- isolated, no
live files touched.)
-----------------------------------------------------------------------------
Runs 02_Proteinomics/06_protein_score_M2_v3_lineage_shrinkage.py's output
through the CURRENT LIVE core_score.py (protein-residual shrinkage fix +
regime-3 substitution both already in place), by monkeypatching PROT_Z,
CORE_SCORE, GENE_DISP to alternate paths before calling run() -- same
technique used to isolate the protein-scorer test itself. Live
bulk_prot_z.parquet, core_score.parquet, gene_dispersion.parquet are never
read from/written to by this script's protein input; RNA input is the live,
unmodified bulk_rna_z.parquet.
"""
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from config import OUT

PROT_SHRUNK_SRC = PROJECT_ROOT / "02_Proteinomics" / "results" / "protein_z_v3_lineage_shrunk.parquet"
PROT_Z_BRIDGED  = OUT / "_test_bulk_prot_z_lineage_shrunk.parquet"
CORE_SCORE_TEST = OUT / "_test_core_score_lineage_shrunk.parquet"
GENE_DISP_TEST  = OUT / "_test_gene_dispersion_lineage_shrunk.parquet"


def bridge_schema():
    """protein_z_v3_lineage_shrunk.parquet (gene_id, model_id, lineage,
    z_score, z_raw, platform_confidence, n_sources, n_isoforms), model_id
    UPPERCASE -> bulk_prot_z.parquet's schema (gene_id, model_id, z_t,
    n_sources), model_id lowercase -- identical bridge to what
    run_protein_m2_stage.py's _export_to_bulk_prot_z() does for the live file."""
    df = pd.read_parquet(PROT_SHRUNK_SRC, columns=["gene_id", "model_id", "z_score", "n_sources"])
    df = df.rename(columns={"z_score": "z_t"})
    df["model_id"] = df["model_id"].str.lower()
    df.to_parquet(PROT_Z_BRIDGED, index=False)
    print(f"Bridged: {len(df):,} rows -> {PROT_Z_BRIDGED}")
    return df


def run_core_score():
    import importlib
    cs = importlib.import_module("core_score")
    cs.PROT_Z = PROT_Z_BRIDGED
    cs.CORE_SCORE = CORE_SCORE_TEST
    cs.GENE_DISP = GENE_DISP_TEST
    cs.run()


if __name__ == "__main__":
    bridge_schema()
    run_core_score()
