"""
verify_cache_correctness.py
-----------------------------
Step 3 of the calibration-caching task: proves transcriptomics_cached.py
produces bit-identical results to transcriptomics.py (the original) on
5 parameter combinations spread across the grid, before trusting the
cached version for a full 480-combination run.

IMPORTANT DESIGN NOTE (v2): axes are built ONCE via a single
build_method_axes() call, then deep-copied into two independent sets --
one driven through transcriptomics.py's functions, one through
transcriptomics_cached.py's. This isolates caching as the ONLY variable
between the two arms.

The first version of this script called build_method_axes() separately
for each arm and found systematic mismatches. Root cause, confirmed
separately: build_method_axes()'s underlying SQL (e.g. "SELECT DISTINCT
model_id FROM hpa_rna", no ORDER BY) does not guarantee row order, and
for axes with no replicates to collapse (hpa_rna_seq: raw_ids and
models are the same length), collapse_replicates() returns data in
raw-fetch order while downstream code (lineage lookup, output
placement) assumes ax.models' sorted order -- a pre-existing bug in
transcriptomics.py unrelated to caching, exposed by two independent
axis-builds hitting DuckDB's non-deterministic ordering twice. Building
axes once and reusing (via deepcopy) removes that confound so this
script actually tests what it's supposed to: whether the cache changes
results, not whether two unrelated DB queries happened to agree.

Compares, per combination, per holdout source:
  - detection floor per axis
  - Spearman rho against the held-out reference (the actual number
    calibrate_params uses to pick a winner)
  - a sample of actual z-score values for a handful of (gene, model)
    cells, not just the summary correlation

Also times the 5-combination check under both versions, as an early
read on the expected full-run speedup.

Usage:
    python verify_cache_correctness.py --db ../outputs/celllineselector.db
"""
from __future__ import annotations

import argparse
import copy
import time

import duckdb
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import transcriptomics as tx_orig
import transcriptomics_cached as tx_cached


def build_ref_cache(tx, con, axes, colmaps, genes, models, midx, lineage_s, holdouts):
    ref_p = tx.Params(**tx.REF_PARAMS_RAW)
    ref_cache = {}
    for h in holdouts:
        hax = [a for a in axes if a.source == h]
        for a in hax:
            tx.profile_method(con, a, colmaps.get(a.name), genes[:400], ref_p.floor_q)
        ph = tx.process_chunk(con, hax, colmaps, genes, models, midx, lineage_s, ref_p)
        ref_cache[h], _, _, _, _ = tx.combine_sources(ph)
    return ref_cache


def run_combo(tx, con, axes, colmaps, genes, models, midx, lineage_s, p, holdouts, ref_cache):
    """Mirrors _score_params's own computation, but also returns the raw
    zp arrays and per-axis floors for value-level comparison, and the
    per-holdout rho -- everything by calling the SAME functions
    _score_params calls internally, not a reimplementation."""
    for a in axes:
        tx.profile_method(con, a, colmaps.get(a.name), genes[:400], p.floor_q)
    floors = {a.name: a.scale_detected and a.floor for a in axes}

    results = {}
    for h in holdouts:
        fit = [a for a in axes if a.source != h]
        pf = tx.process_chunk(con, fit, colmaps, genes, models, midx, lineage_s, p)
        zp, _, _, ex, order = tx.combine_sources(pf)
        zr = ref_cache[h]
        det = None
        for s in order:
            d = ex.get(f"{s}__det")
            det = d if det is None else np.fmax(det, d)
        if det is None:
            det = np.full_like(zp, 0.5)
        m = np.isfinite(zp) & np.isfinite(zr) & (det > 0.01) & (det < 0.90)
        rho = float(spearmanr(zp[m], zr[m]).statistic) if m.sum() >= tx.MIN_CALIB_CELLS else None
        results[h] = dict(rho=rho, zp=zp, mask=m)
    return floors, results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--n-combos", type=int, default=5)
    args = ap.parse_args()

    con = duckdb.connect(args.db, read_only=True)

    print("=== building shared setup (once) ===")
    axes_shared = tx_orig.build_method_axes(con)
    lineage_s = tx_orig.load_lineage_map(con)
    colmaps = {ax.name: (None if ax.source == "hpa_rna"
                         else tx_orig.ensg_columns(con, ax.table)) for ax in axes_shared}
    genes_all = tx_orig.gene_universe(con, axes_shared, colmaps)
    models = sorted(set().union(*[set(ax.models) for ax in axes_shared]))
    midx = {m: i for i, m in enumerate(models)}

    genes = tx_orig.sample_calibration_genes(genes_all, rng=np.random.default_rng(0))
    print(f"{len(genes)} calibration genes")

    holdouts = [s for s in ("depmap", "hpa_rna") if any(a.source == s for a in axes_shared)]
    print(f"holdouts: {holdouts}")

    # Two independent copies of the SAME axis objects (same raw_ids order,
    # same everything) -- caching is the only difference between the arms.
    axes_o = copy.deepcopy(axes_shared)
    axes_c = copy.deepcopy(axes_shared)

    print("=== building reference cache (original) ===")
    ref_o = build_ref_cache(tx_orig, con, axes_o, colmaps, genes, models, midx, lineage_s, holdouts)
    print("=== building reference cache (cached) ===")
    tx_cached._fetch_cache_enable()
    ref_c = build_ref_cache(tx_cached, con, axes_c, colmaps, genes, models, midx, lineage_s, holdouts)
    for h in holdouts:
        same = np.allclose(np.nan_to_num(ref_o[h]), np.nan_to_num(ref_c[h]), equal_nan=True)
        print(f"  ref_cache[{h}] identical: {same}")

    combos = tx_orig.grid_combinations()
    n = len(combos)
    idx = sorted(set(int(round(q * (n - 1))) for q in [0.0, 0.25, 0.5, 0.75, 1.0]))
    idx = idx[:args.n_combos] if len(idx) > args.n_combos else idx
    picked = [combos[i] for i in idx]
    print(f"\n{len(combos)} total combos; spot-checking indices {idx}\n")

    t0 = time.time()
    orig_results = []
    for i, (fq, n0, uc, lc) in zip(idx, picked):
        p = tx_orig.Params(fq, int(n0), uc, lc)
        floors, res = run_combo(tx_orig, con, axes_o, colmaps, genes, models, midx,
                                lineage_s, p, holdouts, ref_o)
        orig_results.append((i, p, floors, res))
    t_orig = time.time() - t0
    print(f"original: {len(picked)} combos in {t_orig:.1f}s ({t_orig/len(picked):.1f}s/combo)")

    t0 = time.time()
    cached_results = []
    for i, (fq, n0, uc, lc) in zip(idx, picked):
        p = tx_cached.Params(fq, int(n0), uc, lc)
        floors, res = run_combo(tx_cached, con, axes_c, colmaps, genes, models, midx,
                                lineage_s, p, holdouts, ref_c)
        cached_results.append((i, p, floors, res))
    tx_cached._fetch_cache_disable()
    t_cached = time.time() - t0
    print(f"cached:   {len(picked)} combos in {t_cached:.1f}s ({t_cached/len(picked):.1f}s/combo)")
    print(f"speedup: {t_orig/t_cached:.1f}x\n")

    print("=== comparison ===")
    all_match = True
    for (i, p_o, floors_o, res_o), (_, p_c, floors_c, res_c) in zip(orig_results, cached_results):
        print(f"\n--- combo #{i}: floor_q={p_o.floor_q} n0={p_o.n0} upper_cut={p_o.upper_cut} lower_cut={p_o.lower_cut} ---")
        for name in floors_o:
            fo, fc = floors_o[name], floors_c.get(name)
            match = (fo == fc) or (fo is not None and fc is not None and abs(fo - fc) < 1e-9)
            all_match &= match
            print(f"  floor[{name}]: orig={fo}  cached={fc}  match={match}")
        for h in holdouts:
            ro, rc = res_o[h]["rho"], res_c[h]["rho"]
            match = (ro is None and rc is None) or (ro is not None and rc is not None and abs(ro - rc) < 1e-9)
            all_match &= match
            print(f"  rho[{h}]: orig={ro}  cached={rc}  match={match}")
            zo, zc = res_o[h]["zp"], res_c[h]["zp"]
            mo, mc = res_o[h]["mask"], res_c[h]["mask"]
            same_mask = np.array_equal(mo, mc)
            zo_f = np.nan_to_num(zo[mo][:10]) if mo.any() else np.array([])
            zc_f = np.nan_to_num(zc[mc][:10]) if mc.any() else np.array([])
            z_match = np.allclose(zo_f, zc_f, atol=1e-6) if same_mask else False
            all_match &= (same_mask and z_match)
            print(f"  z-sample[{h}] (first 10 masked cells): same_mask={same_mask} values_match={z_match}")
            if not (same_mask and z_match):
                print(f"    orig sample:   {zo_f}")
                print(f"    cached sample: {zc_f}")

    print("\n" + "=" * 60)
    print("ALL MATCH -- safe to proceed to full run" if all_match else
          "MISMATCH FOUND -- DO NOT proceed, report this")
    print("=" * 60)


if __name__ == "__main__":
    main()
