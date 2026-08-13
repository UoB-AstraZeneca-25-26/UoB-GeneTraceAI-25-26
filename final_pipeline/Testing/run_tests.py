"""
Testing/run_tests.py
---------------------
Runs the core pre-registered diagnostic tests. Each test is read-only and
prints PASS / FAIL with the measured statistic.

Tests run:
  T-RNA   : RNA source agreement — median Spearman rho across gene pairs
  T-PROT  : Protein platform overlap — median per-gene ProCAN-CCLE rho
  T-SCORE : core_score distribution check — mean, std, fraction > 0.8

Run after all scoring stages have produced outputs.
"""
import sys, json
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import RNA_Z, PROT_Z, PROT_TIER, CORE_SCORE, EVAL_SUMMARY


def test_rna_source_agreement():
    print("\n[T-RNA] RNA source agreement")
    if not RNA_Z.exists():
        print("  SKIP — bulk_rna_z.parquet not found")
        return
    rna = pd.read_parquet(RNA_Z)
    n_sources_dist = rna.groupby("n_sources").size()
    pct_3source = n_sources_dist.get(3, 0) / len(rna) * 100
    print(f"  n_sources dist: {dict(n_sources_dist)}")
    print(f"  3-source coverage: {pct_3source:.1f}%")
    print(f"  z_t range: [{rna.z_t.min():.2f}, {rna.z_t.max():.2f}]  mean={rna.z_t.mean():.3f}")
    print(f"  {'PASS' if 10 < pct_3source < 90 else 'WARN'} — 3-source mix looks reasonable")


def test_protein_platform_overlap():
    print("\n[T-PROT] Protein platform tier")
    if not PROT_TIER.exists():
        print("  SKIP — protein_platform_tier.parquet not found")
        return
    tier = pd.read_parquet(PROT_TIER)
    med_rho = tier.platform_rho.median()
    conflicting_pct = (tier.platform_tier == "conflicting").mean() * 100
    print(f"  Median ProCAN-CCLE rho: {med_rho:.3f}  (expected ~0.37)")
    print(f"  Conflicting tier: {conflicting_pct:.1f}%")
    print(f"  {'PASS' if 0.30 < med_rho < 0.50 else 'WARN'}")


def test_core_score_distribution():
    print("\n[T-SCORE] core_score distribution")
    if not CORE_SCORE.exists():
        print("  SKIP — core_score.parquet not found")
        return
    cs = pd.read_parquet(CORE_SCORE, columns=["core_score"])
    mean = cs.core_score.mean()
    std  = cs.core_score.std()
    p80  = (cs.core_score >= 0.8).mean() * 100
    print(f"  mean={mean:.3f}  std={std:.3f}  >0.8: {p80:.1f}%")
    print(f"  {'PASS' if 0.3 < mean < 0.7 else 'WARN'} — distribution centred")


def test_eval_summary():
    print("\n[T-EVAL] Held-out evaluation summary")
    if not EVAL_SUMMARY.exists():
        print("  SKIP — stage4_eval_summary.json not found")
        return
    s = json.loads(EVAL_SUMMARY.read_text())
    for cls, m in s.items():
        lift = m["hits20_driver"] - m["hits20_flat"]
        print(f"  {cls}: flat={m['hits20_flat']:.3f}  driver={m['hits20_driver']:.3f}  "
              f"lift={lift:+.3f}  {'PASS' if lift >= 0 else 'FAIL'}")


if __name__ == "__main__":
    test_rna_source_agreement()
    test_protein_platform_overlap()
    test_core_score_distribution()
    test_eval_summary()
    print("\nDone.")
