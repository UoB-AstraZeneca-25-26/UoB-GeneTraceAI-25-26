#!/usr/bin/env python3
"""
05_transcriptomics_stat_layer.py

Orchestrator for the transcriptomics statistical layer. All logic lives in
transcriptomics.py; this file only sequences it.

Pipeline
    1. parse geo_meta/*.txt   -> geo_platform  (method = GPL x processing)
    2. build method axes      -> hpa_rna, depmap, N x geo
    3. sensitivity analysis   -> censoring cutoffs, n0, floor quantile
    4. genome-wide z          -> lineage-stratified and global, per cell line
    5. diagnostics + plots

Output table (~20k genes x ~2k cell lines, gene-major):
    ensg | model_id | lineage | z_lineage | z_global | k_src | k_methods
         | det_frac | status | scale_source | z_hpa_rna | z_depmap | z_geo

Usage
    python 05_transcriptomics_stat_layer.py --discover
    python 05_transcriptomics_stat_layer.py --audit-geo
    python 05_transcriptomics_stat_layer.py --build-geo-platform
    python 05_transcriptomics_stat_layer.py --calibrate-only
    python 05_transcriptomics_stat_layer.py \
        --params results/transcriptomics/chosen_params.json \
        --write-table transcriptomics_z
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent   # final_pipeline/
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))  # 01_Transcriptomics/

import transcriptomics as tx
from logging_utils import get_logger

DB_PATH = PROJECT_ROOT / "db" / "celllineselector.duckdb"
RAW_DIR = PROJECT_ROOT / "geo_meta"


def setup_logging(level=logging.INFO):
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-7s %(name)s | %(message)s",
        datefmt="%H:%M:%S",
    )


log = get_logger("05_transcriptomics")


# ----------------------------------------------------------------------

def resolve_project_path(value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (PROJECT_ROOT / path).resolve()


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=str(DB_PATH))
    ap.add_argument("--meta-dir", default=str(RAW_DIR))
    ap.add_argument("--out-dir", default=str(PROJECT_ROOT / "results" / "transcriptomics"))
    ap.add_argument("--write-table", default="transcriptomics_z",
                    help="DuckDB output table for the z-scores")
    ap.add_argument("--genes", default=None, help="comma-separated ENSG subset")
    ap.add_argument("--min-k", type=int, default=1,
                    help="drop rows backed by fewer than this many sources")
    ap.add_argument("--gene-chunk", type=int,
                    default=min(tx.GENE_CHUNK, tx.LOW_MEMORY_GENE_CHUNK),
                    help="Rows of genes per batch; lower values keep CPU/RAM usage down on laptops")
    ap.add_argument("--params", default=None,
                    help="chosen_params.json, to skip the sensitivity analysis")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--discover", action="store_true")
    ap.add_argument("--build-geo-platform", action="store_true")
    ap.add_argument("--audit-geo", action="store_true",
                    help="report GEO method coverage without writing")
    ap.add_argument("--calibrate-only", action="store_true")
    return ap.parse_args()


def step_discover(con):
    log.info("=== schema")
    print(tx.describe_tables(con).to_string(index=False))


def step_geo_platform(con, meta_dir: Path, out: Path, audit_only: bool):
    log.info("=== GEO metadata -> method assignment")
    df, failures = tx.build_geo_platform_frame(meta_dir)
    if len(failures):
        failures.to_csv(out / "geo_meta_parse_failures.csv", index=False)
        log.warning("%d parse failures -> geo_meta_parse_failures.csv",
                    len(failures))
    summary = tx.audit_geo_platform(con, df)
    summary.to_csv(out / "geo_methods.csv", index=False)
    if not audit_only:
        tx.write_geo_platform(con, df)


def step_build_axes(con, args):
    log.info("=== axes")
    axes = tx.build_method_axes(con)
    lineage_s = tx.load_lineage_map(con)

    colmaps = {ax.name: (None if ax.source == "hpa_rna"
                         else tx.ensg_columns(con, ax.table)) for ax in axes}

    if args.genes:
        genes_all = tx.filter_valid_ensg_ids(
            [g.strip() for g in args.genes.split(",") if g.strip()],
            tx.gene_roster_ensg_ids(con),
        )
    else:
        genes_all = tx.gene_universe(con, axes, colmaps)
    log.info("%d genes on the union axis", len(genes_all))

    models = sorted(set().union(*[set(ax.models) for ax in axes]))
    midx = {m: i for i, m in enumerate(models)}
    lineage_vec = lineage_s.reindex(models).fillna("unknown")
    log.info("%d cell lines; top lineages:\n%s", len(models),
             lineage_vec.value_counts().head(12).to_string())
    return axes, colmaps, genes_all, models, midx, lineage_s, lineage_vec


def step_params(con, args, axes, colmaps, genes_all, models, midx, lineage_s, out):
    if args.params:
        p = tx.Params(**json.loads(Path(args.params).read_text()))
        log.info("using supplied parameters: %s", asdict(p))
        return p, pd.DataFrame()
    log.info("=== sensitivity analysis")
    rng = np.random.default_rng(args.seed)
    p, grid, _ = tx.calibrate_params(con, axes, colmaps, genes_all, models,
                                     midx, lineage_s, rng=rng)
    if len(grid):
        grid.to_csv(out / "grid_search.csv", index=False)
    (out / "chosen_params.json").write_text(json.dumps(asdict(p), indent=2))
    return p, grid


def step_profile(con, axes, colmaps, genes_all, p, seed):
    log.info("=== scale detection and detection floors")
    rng = np.random.default_rng(seed + 1)
    probe = list(rng.choice(np.asarray(genes_all, dtype=object),
                            size=min(500, len(genes_all)), replace=False))
    for ax in axes:
        tx.profile_method(con, ax, colmaps.get(ax.name), probe, p.floor_q)


def step_score(con, args, axes, colmaps, genes_all, models, midx, lineage_s,
               lineage_vec, p, out):
    log.info("=== scoring %d genes x %d cell lines in chunks of %d -> %s",
             len(genes_all), len(models), args.gene_chunk, args.write_table)

    first, diag_parts, n_rows = True, [], 0
    lin_np = lineage_vec.to_numpy()
    for ci in range(0, len(genes_all), args.gene_chunk):
        genes = genes_all[ci:ci + args.gene_chunk]

        ps_lin = tx.process_chunk(con, axes, colmaps, genes, models, midx,
                                  lineage_s, p, stratify=True)
        z_lin, k_src, src_lin, extra, order = tx.combine_sources(ps_lin)

        ps_glb = tx.process_chunk(con, axes, colmaps, genes, models, midx,
                                  lineage_s, p, stratify=False)
        z_glb, _, _, _, _ = tx.combine_sources(ps_glb)

        df = tx.chunk_to_frame(models, lin_np, genes, z_lin, z_glb, k_src,
                               extra, order, src_lin, args.min_k)
        if df.empty:
            continue

        tx.write_chunk(con, df, args.write_table, first)
        first = False
        n_rows += len(df)
        diag_parts.append(tx.summarise_chunk(df))
        log.info("  genes %6d-%6d  rows=%9d  cumulative=%d",
                 ci, ci + len(genes), len(df), n_rows)

    diag = (pd.concat(diag_parts)
              .groupby(["lineage", "status", "scale_source", "k_src"],
                       as_index=False).n.sum()) if diag_parts else pd.DataFrame()
    if len(diag):
        diag.to_csv(out / "diagnostics.csv", index=False)
    return diag, n_rows


def report(con, p, diag, n_rows, table):
    print("\n--- parameters (1-SE rule on held-out concordance) ---")
    print(json.dumps(asdict(p), indent=2))
    if len(diag):
        print("\n--- status ---")
        print(diag.groupby("status").n.sum().to_string())
        print("\n--- scale source ---")
        print(diag.groupby("scale_source").n.sum().to_string())
    if n_rows:
        print("\n--- source coverage ---")
        print(tx.verify_output(con, table).to_string(index=False))
    print(f"\n{n_rows:,} rows -> table {table}")


# ----------------------------------------------------------------------

def main():
    args = parse_args()
    setup_logging()

    args.db = str(resolve_project_path(args.db))
    args.meta_dir = str(resolve_project_path(args.meta_dir))
    args.out_dir = str(resolve_project_path(args.out_dir))

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # only pure-reporting modes can open read-only
    read_only = bool(args.discover or
                     (args.audit_geo and not args.build_geo_platform) or
                     args.calibrate_only)
    con = duckdb.connect(args.db, read_only=read_only)

    try:
        if args.discover:
            step_discover(con)
            return

        if args.build_geo_platform or args.audit_geo:
            step_geo_platform(con, Path(args.meta_dir), out, args.audit_geo)
            if args.audit_geo and not args.build_geo_platform:
                return

        axes, colmaps, genes_all, models, midx, lineage_s, lineage_vec = \
            step_build_axes(con, args)

        p, grid = step_params(con, args, axes, colmaps, genes_all, models,
                              midx, lineage_s, out)
        step_profile(con, axes, colmaps, genes_all, p, args.seed)

        if args.calibrate_only:
            tx.make_plots(grid, p, axes, pd.DataFrame(), out)
            print("\n--- parameters ---")
            print(json.dumps(asdict(p), indent=2))
            print(f"\nwritten to {out / 'chosen_params.json'}")
            return

        diag, n_rows = step_score(con, args, axes, colmaps, genes_all, models,
                                  midx, lineage_s, lineage_vec, p, out)
        tx.make_plots(grid, p, axes, diag, out)
        report(con, p, diag, n_rows, args.write_table)
    finally:
        con.close()


if __name__ == "__main__":
    main()