"""
Scoring/_test_huber_lineage_shrinkage.py  (TEST HARNESS -- isolated, no live
files touched.)
-----------------------------------------------------------------------------
Step 3 of the Huber-regression cowork task: runs the lineage-shrinkage-fixed
protein output (already bridged to bulk_prot_z.parquet's schema at
outputs/_test_bulk_prot_z_lineage_shrunk.parquet) through the CORRECTED
Huber-based core_score_v5_huber.py (update_scale=False fix applied), via
monkeypatched PROT_Z/OUT_SCORE/OUT_DISP -- same technique as the OLS version
of this test. No live files touched.
"""
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from config import OUT

PROT_Z_BRIDGED  = OUT / "_test_bulk_prot_z_lineage_shrunk.parquet"   # already built
OUT_SCORE_TEST  = OUT / "_test_core_score_huber_lineage_shrunk.parquet"
OUT_DISP_TEST   = OUT / "_test_gene_dispersion_huber_lineage_shrunk.parquet"


def run():
    import importlib
    cs5 = importlib.import_module("core_score_v5_huber")
    cs5.PROT_Z = PROT_Z_BRIDGED
    cs5.OUT_SCORE = OUT_SCORE_TEST
    cs5.OUT_DISP = OUT_DISP_TEST
    cs5.run()


if __name__ == "__main__":
    run()
