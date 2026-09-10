"""
run_full_calibration_cached.py
--------------------------------
Step 4 of the calibration-caching task: runs the FULL 480-combination
grid search via transcriptomics_cached.py, now that the 5-combination
spot check (verify_cache_correctness.py) passed cleanly -- every floor,
every rho, every z-score sample matched the original exactly.

Mirrors 05_transcriptomics_stats_layer.py's own calibration invocation
(step_build_axes + step_params), but imports transcriptomics_cached
explicitly rather than transcriptomics -- the orchestrator's default
import is untouched pending the Step 5 go/no-go decision on whether to
adopt the cached version permanently.

Writes chosen_params.json and grid_search.csv to
results/transcriptomics_cached/ (separate from the real
results/transcriptomics/ output dir) so there's no ambiguity about
which run produced which file.

Usage:
    python run_full_calibration_cached.py --db ../outputs/celllineselector.db
"""
from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

import transcriptomics_cached as tx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--out-dir", default="results/transcriptomics_cached")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect(args.db, read_only=True)
    try:
        t0 = time.time()
        axes = tx.build_method_axes(con)
        lineage_s = tx.load_lineage_map(con)
        colmaps = {ax.name: (None if ax.source == "hpa_rna"
                             else tx.ensg_columns(con, ax.table)) for ax in axes}
        genes_all = tx.gene_universe(con, axes, colmaps)
        models = sorted(set().union(*[set(ax.models) for ax in axes]))
        midx = {m: i for i, m in enumerate(models)}
        print(f"setup: {len(genes_all)} genes, {len(models)} cell lines, "
              f"{len(axes)} axes ({time.time()-t0:.1f}s)")

        rng = np.random.default_rng(args.seed)
        t0 = time.time()
        p, grid, genes_used = tx.calibrate_params(con, axes, colmaps, genes_all, models,
                                                   midx, lineage_s, rng=rng)
        elapsed = time.time() - t0
        print(f"\nfull calibration: {elapsed:.1f}s ({elapsed/60:.1f} min) "
              f"for {len(grid)} scored combinations")

        if len(grid):
            grid.to_csv(out / "grid_search.csv", index=False)
        (out / "chosen_params.json").write_text(json.dumps(asdict(p), indent=2))
        print(f"\nCHOSEN (1-SE rule): {json.dumps(asdict(p), indent=2)}")
        print(f"written to {out}")
    finally:
        con.close()


if __name__ == "__main__":
    main()
