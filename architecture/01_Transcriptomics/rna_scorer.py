"""
01_Transcriptomics/rna_scorer.py
---------------------------------
Computes lineage-conditioned z-scores for all genes across three RNA sources:
  - DepMap expression (Broad, CCLE v21Q4)
  - HPA RNA (Human Protein Atlas, nTPM)
  - GEO expression (NCBI GEO, log2 MAS5)

Reads from:  celllineselector.db  (Stage-0 warehouse)
Writes to:   final_pipeline/outputs/bulk_rna_z.parquet

Schema:
  gene_id   (str)   — lowercase ENSG
  model_id  (str)   — lowercase ACH-
  z_t       (f32)   — Stouffer-combined, lineage-conditioned z-score
  n_sources (int8)  — number of RNA sources that contributed (1-3)
"""
import sys, time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import RNA_Z
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
import common as C

N_SOURCES = 3
BATCH     = 500   # genes per DuckDB query


def run():
    t0 = time.time()
    print("=" * 70)
    print("01_Transcriptomics — RNA z-scorer (3 sources)")
    print("=" * 70)

    con         = C.connect()
    lineage_map = C.load_lineage(con)
    gene_univ   = C.gene_universe(con)
    genes       = gene_univ["gene_id"].tolist()
    print(f"Gene universe: {len(genes):,}  |  Cell lines: {len(lineage_map):,}")

    all_z = []
    for start in range(0, len(genes), BATCH):
        batch = genes[start:start + BATCH]
        dm = C.fetch_depmap(con, batch)
        hp = C.fetch_hpa(con, batch)
        ge = C.fetch_geo(con, batch)
        for src_df, name in [(dm, "depmap_expr"), (hp, "hpa_rna"), (ge, "geo_expr")]:
            if src_df.empty:
                continue

            def _to_wide(long_df):
                return (long_df.pivot(index="model_id", columns="gene_id", values="value")
                               .reindex(columns=batch))

            z_long = C.score_source_lineage(_to_wide(src_df), batch, lineage_map, name)
            all_z.append(z_long)

        if (start // BATCH + 1) % 10 == 0:
            pct = (start + BATCH) / len(genes) * 100
            print(f"  {pct:.0f}%  {time.time()-t0:.0f}s", flush=True)

    con.close()
    combined = pd.concat(all_z, ignore_index=True)
    result   = C.stouffer_combine(combined, N_SOURCES)
    result["gene_id"] = result["gene_id"].str.upper()   # normalise to canonical ENSG case

    result.to_parquet(RNA_Z, index=False)
    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s  |  rows: {len(result):,}  ->  {RNA_Z}")
    n = result.groupby("n_sources").size()
    print(f"  n_sources distribution: {dict(n)}")


if __name__ == "__main__":
    run()
