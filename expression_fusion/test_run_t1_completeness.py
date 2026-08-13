"""
expression_fusion/test_run_t1_completeness.py
---------------------------------------------
T1: Completeness / fame artefact  — PRE-REGISTRATION GATE TEST

Pre-registration: docs/EXPR_FUSION_PREREGISTRATION.md §4, T1
RNG seeds: BASE_SEED + replicate_index  (BASE_SEED = 20260813, per §3)
Replicates: B = 200

N1 NULL
-------
For each replicate, for each (gene, source) group, permute the z values
within that group. This preserves the missingness mask, per-(gene,source)
coverage, and therefore n_sources per (gene,model_id) pair. Destroys biology.
Then Stouffer-combine the permuted z and compute T1a–T1d.

FAST PATH
---------
Stouffer weights (w = sqrt(n_lines per source)) don't change across replicates.
Precompute them once. Per replicate: permute z, multiply by pre-built w, reduce
via numpy.bincount. Each replicate: O(n_rows), ~1s for 3M rows.

GENE SAMPLE
-----------
N_GENE_SAMPLE = 2000 random genes (seed 42). Artefact is gene-agnostic
(determined by Stouffer formula, not by specific gene biology).

PASS criteria (pre-registered):
  R_sd  in [0.90, 1.10]
  R_abs in [0.90, 1.10]
  T1c   <= 1.10
  T1d   <= 1.10
"""
from __future__ import annotations

import sys, time, json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C

OUT    = HERE / "outputs"
RESULT = OUT / "t1_completeness_result.json"

BASE_SEED       = 20260813
B               = 200
N_GENE_SAMPLE   = 2_000
GENE_SEED       = 42
N_SOURCES_TOTAL = 3
TOP_N           = 1_000
PASS_SD         = (0.90, 1.10)
PASS_ABS        = (0.90, 1.10)
PASS_T1C        = 1.10
PASS_T1D        = 1.10


def _compute_stats(z_t: np.ndarray, n_src: np.ndarray) -> dict:
    m3 = n_src == 3;  m1 = n_src == 1
    if m3.sum() < 10 or m1.sum() < 10:
        return {k: np.nan for k in ("R_sd", "R_abs", "T1c", "T1d")}
    sd3 = np.std(z_t[m3]);  sd1 = np.std(z_t[m1])
    R_sd = float(sd3 / sd1) if sd1 > 0 else np.nan
    ab3  = np.mean(np.abs(z_t[m3]));  ab1 = np.mean(np.abs(z_t[m1]))
    R_abs = float(ab3 / ab1) if ab1 > 0 else np.nan
    base = m3.mean()
    if base == 0 or len(z_t) < TOP_N:
        return {"R_sd": R_sd, "R_abs": R_abs, "T1c": np.nan, "T1d": np.nan}
    top = np.argpartition(z_t, -TOP_N)[-TOP_N:]
    bot = np.argpartition(z_t,  TOP_N)[:TOP_N]
    return {"R_sd": R_sd, "R_abs": R_abs,
            "T1c": float((n_src[top] == 3).mean() / base),
            "T1d": float((n_src[bot] == 3).mean() / base)}


def main():
    t0 = time.time()
    print("=" * 72)
    print("T1 — Completeness artefact  (N1 permutation, B=200)")
    print("Pre-registration: docs/EXPR_FUSION_PREREGISTRATION.md §4")
    print(f"Gene sample: N={N_GENE_SAMPLE}, seed={GENE_SEED}")
    print("=" * 72)

    # ── 1. Load source-level z-scores for gene sample ────────────────────
    print("\n[1] Loading from DB ...", flush=True)
    con         = C.connect()
    lineage_map = C.load_lineage(con)
    gene_univ   = C.gene_universe(con)

    rng_s   = np.random.default_rng(GENE_SEED)
    sidx    = rng_s.choice(len(gene_univ), size=min(N_GENE_SAMPLE, len(gene_univ)),
                           replace=False)
    genes   = gene_univ.iloc[sorted(sidx)]["gene_id"].tolist()
    print(f"  Universe: {len(gene_univ):,}  Sample: {len(genes):,}")

    dm = C.fetch_depmap(con, genes, verbose=True)
    hp = C.fetch_hpa(con,  genes, verbose=True)
    ge = C.fetch_geo(con,  genes, verbose=True)
    con.close()

    # Pivot to wide, then score_source_lineage
    def _to_wide(long_df, val_col="value"):
        return (long_df.pivot(index="model_id", columns="gene_id", values=val_col)
                       .reindex(columns=genes))

    dm_z = C.score_source_lineage(_to_wide(dm), genes, lineage_map, "depmap_expr")
    hp_z = C.score_source_lineage(_to_wide(hp), genes, lineage_map, "hpa_rna")
    if not ge.empty:
        ge_z = C.score_source_lineage(_to_wide(ge), genes, lineage_map, "geo_expr")
    else:
        ge_z = pd.DataFrame(columns=["gene_id", "model_id", "source", "z"])

    all_z = pd.concat([dm_z, hp_z, ge_z], ignore_index=True).dropna(subset=["z"])
    print(f"  Source-level z rows: {len(all_z):,}  ({time.time()-t0:.0f}s)")

    # ── 2. Precompute Stouffer weights (constant across replicates) ───────
    print("\n[2] Precomputing Stouffer weights ...", flush=True)
    # src_n = number of model_ids with valid z for each (gene, source)
    src_n = (all_z.groupby(["gene_id", "source"])["model_id"]
             .nunique().reset_index(name="src_n"))
    wdf = all_z.merge(src_n, on=["gene_id", "source"])
    wdf["w"] = np.sqrt(wdf["src_n"])
    wdf = wdf.sort_values(["gene_id", "source", "model_id"]).reset_index(drop=True)

    # group indices for (gene, source) — for permutation
    gs_key = (wdf["gene_id"] + "||" + wdf["source"]).values
    _, gs_inv, gs_cnt = np.unique(gs_key, return_inverse=True, return_counts=True)
    gs_offsets = np.concatenate([[0], np.cumsum(gs_cnt)])

    # group indices for (gene, model_id) — for Stouffer aggregation
    # wdf has one row per (gene_id, model_id, source) — so bincount gives n_sources per pair
    gm_key = (wdf["gene_id"] + "\x00" + wdf["model_id"]).values
    unique_pairs, gm_inv = np.unique(gm_key, return_inverse=True)
    n_pairs = len(unique_pairs)

    w_arr  = wdf["w"].values.astype(float)
    w2_arr = w_arr ** 2

    # w2_sum per (gene, model_id): constant across replicates
    w2_sum_arr = np.bincount(gm_inv, weights=w2_arr, minlength=n_pairs)
    # n_sources per pair: one row per source in wdf, so bincount = n_sources
    n_src_per_pair = np.bincount(gm_inv, minlength=n_pairs).astype(int)

    safe_n     = np.where(n_src_per_pair > 0, n_src_per_pair, 1)
    shrink_arr = np.where(n_src_per_pair < N_SOURCES_TOTAL,
                          np.sqrt(N_SOURCES_TOTAL / safe_n), 1.0)
    denom_arr  = np.sqrt(w2_sum_arr)

    z_base = wdf["z"].values.astype(float)
    print(f"  {n_pairs:,} (gene, model_id) pairs, "
          f"{len(gs_offsets)-1} (gene, source) groups  ({time.time()-t0:.0f}s)")

    # ── 3. Real-data statistics ───────────────────────────────────────────
    wz_real   = w_arr * z_base
    wz_sum_r  = np.bincount(gm_inv, weights=wz_real,  minlength=n_pairs)
    z_raw_r   = np.where(denom_arr > 0, wz_sum_r / np.where(denom_arr > 0, denom_arr, 1.0), np.nan)
    z_t_real  = z_raw_r / shrink_arr
    real_stats = _compute_stats(z_t_real, n_src_per_pair)
    print("\n[3] Real-data statistics (biological signal present):")
    from collections import Counter
    cnt = Counter(n_src_per_pair.tolist())
    print(f"  n_sources: {dict(sorted(cnt.items()))}")
    for k, v in real_stats.items():
        print(f"  {k} = {v:.4f}")

    # ── 4. N1 replicates (fast path) ──────────────────────────────────────
    print(f"\n[4] Running {B} N1 replicates (fast path) ...", flush=True)
    perm_stats: dict[str, list] = {"R_sd": [], "R_abs": [], "T1c": [], "T1d": []}

    for rep in range(B):
        rng_rep = np.random.default_rng(BASE_SEED + rep)

        # Permute z within each (gene, source) group (vectorised)
        z_perm = z_base.copy()
        for g in range(len(gs_offsets) - 1):
            lo, hi = gs_offsets[g], gs_offsets[g + 1]
            if hi - lo > 1:
                z_perm[lo:hi] = rng_rep.permutation(z_perm[lo:hi])

        # Stouffer combine via numpy bincount (no pandas needed)
        wz_perm  = w_arr * z_perm
        wz_sum_p = np.bincount(gm_inv, weights=wz_perm, minlength=n_pairs)
        z_raw_p  = np.where(denom_arr > 0,
                            wz_sum_p / np.where(denom_arr > 0, denom_arr, 1.0),
                            np.nan)
        z_t_perm = z_raw_p / shrink_arr

        s = _compute_stats(z_t_perm, n_src_per_pair)
        for k in perm_stats:
            perm_stats[k].append(s[k])

        if (rep + 1) % 20 == 0:
            print(f"  rep {rep+1:3d}/{B}  "
                  f"R_sd={np.nanmedian(perm_stats['R_sd']):.3f}  "
                  f"R_abs={np.nanmedian(perm_stats['R_abs']):.3f}  "
                  f"T1c={np.nanmedian(perm_stats['T1c']):.3f}  "
                  f"({time.time()-t0:.0f}s)", flush=True)

    # ── 5. Aggregate and report ───────────────────────────────────────────
    print("\n-- N1 replicate summary (pre-registered) --")
    results: dict = {}
    all_pass = True

    for key, values in perm_stats.items():
        arr = np.array(values, dtype=float)
        med = float(np.nanmedian(arr))
        p5  = float(np.nanpercentile(arr,  5))
        p95 = float(np.nanpercentile(arr, 95))
        if   key == "R_sd":  passed = PASS_SD[0]  <= med <= PASS_SD[1];  band = f"[{PASS_SD[0]}, {PASS_SD[1]}]"
        elif key == "R_abs": passed = PASS_ABS[0] <= med <= PASS_ABS[1]; band = f"[{PASS_ABS[0]}, {PASS_ABS[1]}]"
        elif key == "T1c":   passed = med <= PASS_T1C;                    band = f"<= {PASS_T1C}"
        else:                passed = med <= PASS_T1D;                    band = f"<= {PASS_T1D}"
        verdict  = "PASS" if passed else "FAIL"
        all_pass = all_pass and passed
        results[key] = {"median": med, "p5": p5, "p95": p95,
                        "verdict": verdict, "values": arr.tolist()}
        print(f"  {key:6s}  median={med:7.4f}  95%CI=[{p5:.4f}, {p95:.4f}]  "
              f"band={band:18s}  {verdict}")

    overall = "PASS" if all_pass else "FAIL"
    print(f"\n{'='*72}")
    print(f"T1 VERDICT: {overall}")
    if overall == "FAIL":
        print("KILL SWITCH K1 TRIGGERED.")
        print("Stouffer shrinkage design is withdrawn.")
        print("See EXPR_FUSION_PREREGISTRATION.md §4 for permitted replacements.")
    print(f"Total elapsed: {time.time()-t0:.0f}s")
    print(f"{'='*72}")
    print("\nAdditional (required): Spearman(n_sources, study_freq) needs "
          "external study-frequency index; prior measurement rho ~ 0.60.")

    output = {
        "test": "T1",
        "preregistration": "docs/EXPR_FUSION_PREREGISTRATION.md §4",
        "B": B, "base_seed": BASE_SEED,
        "gene_sample_n": len(genes), "gene_sample_seed": GENE_SEED,
        "overall_verdict": overall, "kill_switch_K1": overall == "FAIL",
        "real_data_stats": real_stats,
        "n1_replicates": {k: {kk: vv for kk, vv in v.items() if kk != "values"}
                          for k, v in results.items()},
    }
    RESULT.write_text(json.dumps(output, indent=2))
    print(f"\nResult saved -> {RESULT}")


if __name__ == "__main__":
    main()
