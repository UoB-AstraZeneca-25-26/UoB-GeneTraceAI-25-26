"""
run_gene_tests.py
------------------
Runs all 8 validation genes through the CLI and prints results.

Genes tested:
  FGFR2  — no proteomics
  CD86   — has proteomics
  ERBB2  — all 3 sources (best represented)
  KLK4   — no proteomics, few GEO cell lines
  GAPDH  — housekeeping; should NOT produce many extreme high scores
  ASGR1  — transcriptomics + proteomics, little GEO
  MUC1   — proteomics data
  CD3E   — T-cell marker; few GEO/no proteomics, interesting biological case
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / "Ranking"))

from Ranking.rank_cell_lines import rank

GENES = [
    ("FGFR2",  "Expected: RNA-only scores; no protein layer"),
    ("CD86",   "Expected: RNA + protein combined scores"),
    ("ERBB2",  "Expected: all 3 sources; strong variation across lineages"),
    ("KLK4",   "Expected: mostly DepMap RNA only; few GEO lines"),
    ("GAPDH",  "Expected: housekeeping — scores compressed near median, few extremes"),
    ("ASGR1",  "Expected: transcriptomics + proteomics; sparse GEO"),
    ("MUC1",   "Expected: proteomics present; tissue-specific expression"),
    ("CD3E",   "Expected: T-cell marker — very few lines express highly; sparse data"),
]

print("=" * 70)
print("GeneTraceAI — Gene Validation Panel")
print("=" * 70)

for gene, note in GENES:
    print(f"\n{'=' * 70}")
    print(f"  {gene}  |  {note}")
    print("=" * 70)
    try:
        df = rank(gene)
        if df is not None and not df.empty:
            tier_counts = df.groupby("confidence_tier").size()
            print(f"\n  Tier breakdown: {dict(tier_counts)}")
            print(f"  Score range: {df.core_score.min():.3f} - {df.core_score.max():.3f}  "
                  f"mean={df.core_score.mean():.3f}  std={df.core_score.std():.3f}")
            n_lay = df.get("n_layers", None)
            if n_lay is not None:
                print(f"  Layers present: {df.n_layers.value_counts().to_dict()}")
    except Exception as e:
        print(f"  ERROR: {e}")

print(f"\n{'=' * 70}")
print("Gene validation complete.")
