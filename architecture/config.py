"""
final_pipeline/config.py
------------------------
Single source of truth for all paths. Import this in every stage script.
"""
from pathlib import Path

# ── repo root ──────────────────────────────────────────────────────────────
REPO = Path(__file__).resolve().parent.parent      # UoB-GeneTraceAI-25-26/
PIPELINE = Path(__file__).resolve().parent          # final_pipeline/

# ── read-only shared resources ─────────────────────────────────────────────
# Preferred location is inside final_pipeline/ (self-contained bundle); the
# tracked copies live at the repo root, so fall back there when absent.
def _shared(name: str) -> Path:
    local = PIPELINE / name
    return local if local.exists() else REPO / name


DB         = PIPELINE / "outputs/celllineselector.db"
REF        = _shared("reference")
GENE_LKP   = REF / "gene_lookup.parquet"
CELL_LKP   = REF / "cell_line_lookup.parquet"
CLEANED    = _shared("cleaned_track_data")
VALIDATION = PIPELINE / "validation/prepared"
# Not used by pipeline stages (only by harmonisation notebooks):
DATA_CLEAN = REPO / "data/parquet"
COSMIC_CNA = REPO / "data/COSMIC/CellLinesProject_CompleteCNA_v104_GRCh37.tsv"
GDSC_MODEL = REPO / "data/GDSC/model_list_20260709.csv"

# ── final_pipeline outputs (written by each stage) ─────────────────────────
OUT = PIPELINE / "outputs"

# stage outputs — import these instead of hard-coding paths
RNA_Z          = OUT / "bulk_rna_z.parquet"
PROT_Z         = OUT / "bulk_prot_z.parquet"
PROT_TIER      = OUT / "protein_platform_tier.parquet"
CNA_FLAGS      = OUT / "cna_flags.parquet"
MUTATIONS_COLL = OUT / "mutations_collapsed.parquet"
MUTATIONS_SCR  = OUT / "mutations_scores.parquet"
FUSIONS_GENE   = OUT / "fusions_gene_level.parquet"
FUSIONS_SCR    = OUT / "fusions_scores.parquet"
CORE_SCORE     = OUT / "core_score.parquet"
GENE_DISP      = OUT / "gene_dispersion.parquet"
FLAGS_DRIVER   = OUT / "flags_with_driver.parquet"
PREDICTIONS    = OUT / "predictions_with_confidence.parquet"
EVIDENCE_LDG   = OUT / "evidence_ledger.parquet"
EVAL_SUMMARY   = OUT / "stage4_eval_summary.json"

OUT.mkdir(exist_ok=True)
