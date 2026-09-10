#!/usr/bin/env python3
"""
05_transcriptomics_stat_layer.py

Orchestrator for the transcriptomics statistical layer. All logic lives in
``src/scripts/transcriptomics.py``; this file only sequences it, so each
``step_*`` function below is a thin wrapper that pulls arguments together,
calls into ``tx``, and logs or writes what comes back.

Pipeline
--------
1. parse ``geo_meta/*.txt``   -> ``geo_platform`` (method = GPL x processing)
2. build method axes          -> ``hpa_rna``, ``depmap``, N x ``geo``
3. sensitivity analysis       -> censoring cutoffs, ``n0``, floor quantile
4. genome-wide z              -> lineage-stratified and global, per cell line
5. diagnostics + plots

Output table (~20k genes x ~2k cell lines, gene-major)
------------------------------------------------------
    ensg | model_id | lineage | z_lineage | z_global | k_src | k_methods
         | det_frac | status | scale_source | z_hpa_rna | z_depmap | z_geo

Design notes
------------
* Scoring runs in gene chunks rather than all at once, appending each
  chunk to the output table, so peak memory stays flat regardless of how
  many genes are on the union axis.
* Every gene chunk is processed twice — once lineage-stratified, once
  global — giving the paired ``z_lineage`` and ``z_global`` columns.
* Read-only DuckDB connections are used for the pure-reporting modes, so
  a discovery or audit run cannot modify the database.
* Parameters can be supplied via ``--params`` to skip the sensitivity
  analysis entirely, which makes a scoring run reproducible from a saved
  ``chosen_params.json``.

Usage
-----
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

PROJECT_ROOT = Path(__file__).resolve().parent
for candidate in [PROJECT_ROOT, *PROJECT_ROOT.parents]:
    if (candidate / "src").exists():
        PROJECT_ROOT = candidate
        break
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.scripts import transcriptomics as tx
from src.scripts.logging_utils import get_logger

DB_PATH = PROJECT_ROOT / "db" / "celllineselector.duckdb"
RAW_DIR = PROJECT_ROOT / "geo_meta"


def setup_logging(level=logging.INFO):
    """
    Configure root logging for the run.

    Sets a compact console format — time, level, logger name, message —
    with second-resolution timestamps, suitable for following a long
    chunked scoring run in a terminal.

    Parameters
    ----------
    level : int, optional
        Standard :mod:`logging` level. Defaults to ``logging.INFO``.

    Returns
    -------
    None
        Configures the root logger as a side effect.

    Notes
    -----
    Uses :func:`logging.basicConfig`, which is a no-op if the root logger
    already has handlers attached.
    """
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-7s %(name)s | %(message)s",
        datefmt="%H:%M:%S",
    )


log = get_logger("05_transcriptomics")


# ----------------------------------------------------------------------

def resolve_project_path(value: str | Path) -> Path:
    """
    Resolve a path argument against the project root.

    Absolute paths are returned unchanged; relative ones are interpreted
    relative to ``PROJECT_ROOT`` rather than the current working
    directory, so the CLI defaults behave the same whether the script is
    run from the repository top or from ``pipelines/``.

    Parameters
    ----------
    value : str or Path
        Path as supplied on the command line.

    Returns
    -------
    Path
        Absolute path. Relative inputs are additionally resolved, so
        ``..`` segments and symlinks are collapsed.
    """
    path = Path(value)
    if path.is_absolute():
        return path
    return (PROJECT_ROOT / path).resolve()


def parse_args():
    """
    Define and parse the command-line interface.

    The parser exposes four mutually compatible mode flags — ``--discover``,
    ``--audit-geo``, ``--build-geo-platform`` and ``--calibrate-only`` —
    plus the paths, scoring parameters and subsetting options used by a
    full run. The module docstring is reused as the help description.

    Returns
    -------
    argparse.Namespace
        Parsed arguments. Notable fields:

        * ``db``, ``meta_dir``, ``out_dir`` (str) — paths, resolved
          against the project root later in :func:`main`;
        * ``write_table`` (str) — DuckDB output table for the z-scores;
        * ``genes`` (str or None) — comma-separated ENSG subset;
        * ``min_k`` (int) — drop rows backed by fewer than this many
          sources;
        * ``gene_chunk`` (int) — genes per scoring batch; the default is
          the smaller of ``tx.GENE_CHUNK`` and ``tx.LOW_MEMORY_GENE_CHUNK``
          so laptops keep CPU and RAM usage down;
        * ``params`` (str or None) — path to ``chosen_params.json``, to
          skip the sensitivity analysis;
        * ``seed`` (int) — RNG seed for calibration and profiling.
    """
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
    """
    Print the schema of the tables this layer reads.

    Inspection mode for checking what is actually in the database before
    committing to a run — column names, types and shapes of the source
    tables.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection, typically read-only in this mode.

    Returns
    -------
    None
        The schema table is printed to stdout.
    """
    log.info("=== schema")
    print(tx.describe_tables(con).to_string(index=False))


def step_geo_platform(con, meta_dir: Path, out: Path, audit_only: bool):
    """
    Parse GEO metadata into method assignments and audit the coverage.

    Reads every ``*.txt`` under ``meta_dir``, derives a method per sample
    as GPL crossed with processing type, audits how that assignment covers
    the GEO samples already in the database, and — unless auditing only —
    writes the resulting ``geo_platform`` table.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection. Must be writable unless ``audit_only`` is True.
    meta_dir : Path
        Directory of raw GEO metadata text files.
    out : Path
        Output directory for the CSV reports.
    audit_only : bool
        When True, report coverage without writing ``geo_platform``.

    Returns
    -------
    None
        Writes ``geo_methods.csv`` to ``out`` always, and
        ``geo_meta_parse_failures.csv`` when any file failed to parse.

    Notes
    -----
    Parse failures are reported and logged as a warning, not raised — a
    handful of unparseable metadata files should not block the rest of the
    assignment. Check the failures CSV to see what was dropped.
    """
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
    """
    Build the method axes and the gene and cell-line universes to score over.

    Constructs one axis per measurement method (``hpa_rna``, ``depmap``,
    and one per GEO method), loads the lineage map, and resolves the
    ENSG column list for each axis. ``hpa_rna`` is long-format rather than
    a gene-column matrix, so it gets no column map.

    The gene universe is either the ``--genes`` subset — validated against
    ``gene_roster`` so an invalid ENSG cannot silently produce an empty
    result — or the union across all axes. The cell-line universe is the
    union of models across axes, with lineages defaulting to ``unknown``
    where the map has no entry.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    args : argparse.Namespace
        Parsed arguments; ``genes`` is the only field read.

    Returns
    -------
    axes : list
        Method axis objects, each carrying ``name``, ``source``, ``table``
        and ``models``.
    colmaps : dict
        Axis name to its ENSG column map, or None for ``hpa_rna``.
    genes_all : list of str
        ENSG IDs on the union axis, in scoring order.
    models : list of str
        Sorted union of ``model_id`` values across all axes.
    midx : dict
        ``model_id`` to its integer position in ``models``.
    lineage_s : pandas.Series
        Raw lineage map as loaded, indexed by ``model_id``.
    lineage_vec : pandas.Series
        Lineage per entry of ``models``, aligned and gap-filled with
        ``"unknown"``.

    Notes
    -----
    ``midx`` and ``lineage_vec`` exist so the scoring loop can address
    cell lines positionally in NumPy arrays rather than by label.
    """
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
    """
    Obtain the scoring parameters, either from file or by calibration.

    If ``--params`` was supplied, the saved JSON is loaded and returned
    as-is, skipping calibration entirely — this is what makes a scoring
    run reproducible. Otherwise a sensitivity analysis runs a grid search
    over censoring cutoffs, ``n0`` and the floor quantile, selecting by the
    1-SE rule on held-out concordance.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    args : argparse.Namespace
        Parsed arguments; ``params`` and ``seed`` are read.
    axes, colmaps, genes_all, models, midx, lineage_s
        Universe and axis objects from :func:`step_build_axes`.
    out : Path
        Output directory for ``grid_search.csv`` and
        ``chosen_params.json``.

    Returns
    -------
    p : tx.Params
        The chosen parameters.
    grid : pandas.DataFrame
        Full grid-search results. Empty when parameters were supplied via
        ``--params``, since no search ran.

    Notes
    -----
    The chosen parameters are always written to
    ``out/chosen_params.json`` when calibration runs, so a later run can
    reproduce this one with ``--params``. Nothing is written in the
    supplied-parameters branch.
    """
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
    """
    Detect each axis's value scale and set its detection floor.

    Profiles every axis against a random probe of up to 500 genes, which
    is enough to identify whether values are log or linear and to place
    the detection floor at quantile ``p.floor_q`` without reading the full
    matrix.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    axes : list
        Method axis objects from :func:`step_build_axes`.
    colmaps : dict
        Axis name to ENSG column map.
    genes_all : list of str
        Gene universe to sample the probe from.
    p : tx.Params
        Chosen parameters; ``floor_q`` is read.
    seed : int
        Base RNG seed. Offset by 1 so the probe differs from the draw used
        during calibration.

    Returns
    -------
    None
        Each axis is updated in place with its detected scale and floor.

    Notes
    -----
    Must run before :func:`step_score`, which relies on the detected
    scale and floors when combining sources. The resulting scale choice
    surfaces per row as ``scale_source`` in the output table.
    """
    log.info("=== scale detection and detection floors")
    rng = np.random.default_rng(seed + 1)
    probe = list(rng.choice(np.asarray(genes_all, dtype=object),
                            size=min(500, len(genes_all)), replace=False))
    for ax in axes:
        tx.profile_method(con, ax, colmaps.get(ax.name), probe, p.floor_q)


def step_score(con, args, axes, colmaps, genes_all, models, midx, lineage_s,
               lineage_vec, p, out):
    """
    Score every gene against every cell line, in chunks, and write the results.

    Iterates the gene universe in batches of ``args.gene_chunk``. Each
    batch is processed twice — once lineage-stratified and once globally —
    and the two z-score sets are assembled into one frame alongside the
    per-source columns, source counts and status flags. Rows backed by
    fewer than ``args.min_k`` sources are dropped.

    Each chunk is appended to the output table as it completes, so memory
    stays flat across the full ~20k gene x ~2k cell line run and progress
    survives as rows in the database rather than only in memory.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open, writable connection.
    args : argparse.Namespace
        Parsed arguments; ``gene_chunk``, ``min_k`` and ``write_table``
        are read.
    axes, colmaps, genes_all, models, midx, lineage_s, lineage_vec
        Universe and axis objects from :func:`step_build_axes`.
    p : tx.Params
        Chosen parameters, with scales and floors already set by
        :func:`step_profile`.
    out : Path
        Output directory for ``diagnostics.csv``.

    Returns
    -------
    diag : pandas.DataFrame
        Per-chunk summaries aggregated into counts by ``lineage``,
        ``status``, ``scale_source`` and ``k_src``. Empty if no chunk
        produced rows.
    n_rows : int
        Total rows written to the output table.

    Notes
    -----
    The first chunk creates or replaces the output table (``first=True``);
    every later chunk appends. An empty chunk is skipped without
    consuming that flag, so the table is still created correctly if the
    first few chunks yield nothing.
    """
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
    """
    Print the end-of-run summary.

    Reports the chosen parameters, the status and scale-source breakdowns
    from the diagnostics, source coverage read back from the written
    table, and the total row count.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection, used to verify the written output.
    p : tx.Params
        Parameters the run used.
    diag : pandas.DataFrame
        Aggregated diagnostics from :func:`step_score`.
    n_rows : int
        Total rows written.
    table : str
        Name of the output table.

    Returns
    -------
    None
        Output is printed to stdout.

    Notes
    -----
    Each block is conditional: the diagnostics breakdowns are skipped when
    ``diag`` is empty, and coverage is only read back when rows were
    actually written.
    """
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
    """
    Parse arguments, open the database, and run the requested mode.

    Resolves all path arguments against the project root, creates the
    output directory, and opens DuckDB read-only for the pure-reporting
    modes (``--discover``, ``--audit-geo`` without
    ``--build-geo-platform``, ``--calibrate-only``) so they cannot modify
    the database.

    Modes, in the order they are checked:

    * ``--discover`` — print the schema and stop.
    * ``--build-geo-platform`` / ``--audit-geo`` — assign GEO methods,
      writing ``geo_platform`` unless auditing only. Auditing alone stops
      here; combined with ``--build-geo-platform`` the run continues.
    * ``--calibrate-only`` — build axes, calibrate, profile, emit plots
      and ``chosen_params.json``, then stop before scoring.
    * default — the full run: axes, parameters, profiling, chunked
      scoring, plots and report.

    Returns
    -------
    None
        The z-score table, CSV reports, plots and log output are written
        as side effects.

    Notes
    -----
    The connection is closed in a ``finally`` block, so it is released
    even when a step raises or an early-return mode exits.
    """
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