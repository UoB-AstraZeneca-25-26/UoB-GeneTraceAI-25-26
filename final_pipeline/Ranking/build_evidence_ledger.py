"""
Ranking/build_evidence_ledger.py
----------------------------------
Bulk-builds the evidence ledger for all (model_id, ensg_id) pairs in
predictions_with_confidence.parquet. A projection of existing columns — no
new modelling. Shared derivation logic lives in evidence_ledger.py.

Reads from:  final_pipeline/outputs/predictions_with_confidence.parquet
             final_pipeline/outputs/{core_score, cna_flags}.parquet
             reference/{gene_lookup, cell_line_lookup}.parquet
Writes to:   final_pipeline/outputs/evidence_ledger.parquet
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import PREDICTIONS, CORE_SCORE, CNA_FLAGS, GENE_LKP, CELL_LKP, EVIDENCE_LDG
sys.path.insert(0, str(Path(__file__).resolve().parent))
from evidence_ledger import build_ledger_row


def run():
    print("=" * 70)
    print("Ranking — building evidence ledger (all pairs)")
    print("=" * 70)

    for p in [PREDICTIONS, CORE_SCORE, CNA_FLAGS]:
        if not p.exists():
            raise SystemExit(f"Missing input: {p}")

    pred  = pd.read_parquet(PREDICTIONS)
    core  = pd.read_parquet(CORE_SCORE)
    cna   = pd.read_parquet(CNA_FLAGS)
    genes = pd.read_parquet(GENE_LKP)
    lines = pd.read_parquet(CELL_LKP)

    ledger = build_ledger_row(pred, core, cna, genes, lines)
    ledger.to_parquet(EVIDENCE_LDG, index=False)
    print(f"Evidence ledger: {len(ledger):,} rows  ->  {EVIDENCE_LDG}")


if __name__ == "__main__":
    run()
