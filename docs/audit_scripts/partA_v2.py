"""PART A RE-RUN with the gene-complete CN layer wired in.

Judged against docs/RANKING_PRESPEC.md, written before this ran.
Paired with the previous run: identical 300 genes (seed 42), identical lines,
identical rna_gate / protein / alteration / geo_corrob checks reused from cache.
ONLY the cna_gate input changes.
"""
import os, json, numpy as np, pandas as pd
from scipy import stats
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
W = os.path.dirname(os.path.abspath(__file__))
import rank_checks as RC

DEL_PRIMARY = 0.5
DEL_SENS = (0.25, 0.75)

# ---- reuse the previous run's checks so the comparison is paired ---------------
C_old = {c: pd.read_parquet(os.path.join(W, f"partA_check_{c}.parquet")) for c in RC.CHECKS}
genes = list(C_old["rna_gate"].columns)
idx = C_old["rna_gate"].index
print(f"paired re-run: {len(genes)} genes, {len(idx):,} lines")

# ---- the new CN check ---------------------------------------------------------
REL = pd.read_parquet("src/pipeline/outputs/cn_gene_complete.parquet")
have = [g for g in genes if g in REL.columns]
print(f"gene-complete CN covers {len(have)}/{len(genes)} sampled genes, "
      f"{REL.shape[0]:,} lines")
R = REL.reindex(index=idx, columns=genes)


def cn_check(thr):
    """assessed where the CN matrix covers the pair; FAILS where deleted."""
    return pd.DataFrame(np.where(R.notna(), (R >= thr).astype(float), np.nan),
                        index=idx, columns=genes)


cov = pd.read_parquet("src/pipeline/outputs/coverage_matrix_enriched.parquet").set_index("model_id")
frac_expr = C_old["rna_gate"].mean(axis=0)
valid = frac_expr >= RC.SILENT_FRAC
gs_valid = [g for g in genes if valid[g]]
retained = C_old["rna_gate"] == 1.0
print(f"valid genes: {len(gs_valid)}/{len(genes)};  retained pairs "
      f"{int(retained.to_numpy().sum()):,} ({100*retained.to_numpy().mean():.1f}%)")


def coverage_count(cn_assessed_lines):
    """gene-independent evidence-type availability per line, 0-5, with the CN
    availability taken from the NEW layer."""
    cc = pd.Series(0.0, index=idx)
    for chk, cols in RC.COVERAGE_COLS.items():
        if chk == "cna_gate":
            cc += pd.Series(idx.isin(cn_assessed_lines), index=idx).astype(float)
        else:
            cc += cov.reindex(idx)[cols].any(axis=1).fillna(False).astype(float)
    return cc


def evaluate(C, label, gs, tag):
    """The four pre-specified criteria on the PRIMARY statistic:
    raw count, >=2 assessed checks, gate excluded, within the retained set."""
    others = [c for c in RC.CHECKS if c != "rna_gate"]
    npass = sum(C[c].reindex(index=idx, columns=genes).fillna(0.0) for c in others)
    nass = sum(C[c].reindex(index=idx, columns=genes).notna().astype(float) for c in others)
    M = npass.where(retained & (nass >= 2))[gs]

    cn_lines = set(R.dropna(how="all").index)
    cc = coverage_count(cn_lines)

    a = M.to_numpy().ravel()
    b = np.repeat(cc.to_numpy()[:, None], len(gs), axis=1).ravel()
    m = np.isfinite(a) & np.isfinite(b)
    A1 = stats.spearmanr(b[m], a[m]).statistic if m.sum() > 100 else np.nan

    lg, nd = [], []
    for g in gs:
        y = M[g].dropna()
        if len(y) < 50:
            continue
        vc = y.value_counts()
        lg.append(vc.max() / len(y)); nd.append(len(vc))
    A2a = float(np.median(lg)) if lg else np.nan
    A2b = float(np.median(nd)) if nd else np.nan

    rng = np.random.default_rng(0); rr = []
    for i, j in rng.choice(len(gs), size=(400, 2)):
        if i == j:
            continue
        x, y = M[gs[i]], M[gs[j]]
        o = x.notna() & y.notna()
        if o.sum() > 50:
            v = stats.spearmanr(x[o], y[o]).statistic
            if np.isfinite(v):
                rr.append(v)
    A3 = float(np.median(rr)) if rr else np.nan

    kept = float((retained & (nass >= 2))[gs].to_numpy().mean())
    return {"tag": tag, "set": label, "n_genes": len(gs), "A1_coverage_rho": A1,
            "A2a_largest_tie": A2a, "A2b_distinct_values": A2b,
            "A3_between_gene_rho": A3, "frac_pairs_usable": kept,
            "median_lines_per_gene": float(np.median([len(M[g].dropna()) for g in gs]))}


def verdict(r):
    if not np.isfinite(r["A1_coverage_rho"]):
        return "INCONCLUSIVE (not computable)"
    adopt = (abs(r["A1_coverage_rho"]) <= 0.20 and r["A2a_largest_tie"] <= 0.50
             and r["A2b_distinct_values"] >= 4 and abs(r["A3_between_gene_rho"]) <= 0.25)
    aband = (abs(r["A1_coverage_rho"]) >= 0.40 or r["A2a_largest_tie"] >= 0.65
             or r["A2b_distinct_values"] <= 2 or abs(r["A3_between_gene_rho"]) >= 0.40)
    return "ADOPT" if adopt else ("ABANDON" if aband else "INCONCLUSIVE")


rows = []
# BEFORE: previous 5-check configuration, same sample
rows.append({**evaluate(C_old, "valid genes", gs_valid, "BEFORE (COSMIC CNA)"),
             "verdict": None})
rows.append({**evaluate(C_old, "all genes", genes, "BEFORE (COSMIC CNA)"), "verdict": None})
# AFTER at each CN threshold
for thr in (DEL_PRIMARY, *DEL_SENS):
    C_new = dict(C_old); C_new["cna_gate"] = cn_check(thr)
    tag = f"AFTER (rel<{thr})" + ("  <-- PRIMARY" if thr == DEL_PRIMARY else "")
    rows.append({**evaluate(C_new, "valid genes", gs_valid, tag), "verdict": None})
    rows.append({**evaluate(C_new, "all genes", genes, tag), "verdict": None})
T = pd.DataFrame(rows)
for i in T.index:
    T.loc[i, "verdict"] = verdict(T.loc[i])

pd.set_option("display.width", 220)
print("\n" + "=" * 118)
print("PART A RE-RUN -- judged against docs/RANKING_PRESPEC.md")
print("=" * 118)
print(T[["tag", "set", "n_genes", "A1_coverage_rho", "A2a_largest_tie",
         "A2b_distinct_values", "A3_between_gene_rho", "frac_pairs_usable",
         "median_lines_per_gene", "verdict"]].round(4).to_string(index=False))

print("\nPre-specified thresholds:")
print("  ADOPT   : |A1| <= 0.20  AND  A2a <= 0.50  AND  A2b >= 4  AND  |A3| <= 0.25")
print("  ABANDON : |A1| >= 0.40  OR   A2a >= 0.65  OR   A2b <= 2  OR   |A3| >= 0.40")

prim = T[(T.tag.str.startswith(f"AFTER (rel<{DEL_PRIMARY})")) & (T.set == "valid genes")].iloc[0]
print(f"\nDECISION (primary: valid genes, rel<{DEL_PRIMARY}, >=2 assessed, within retained set)")
print(f"  A1 coverage rho     = {prim.A1_coverage_rho:+.4f}   "
      f"(adopt <= 0.20, abandon >= 0.40)")
print(f"  A2a largest tie     = {prim.A2a_largest_tie:.4f}   (adopt <= 0.50, abandon >= 0.65)")
print(f"  A2b distinct values = {prim.A2b_distinct_values:.0f}        (adopt >= 4, abandon <= 2)")
print(f"  A3 between-gene rho = {prim.A3_between_gene_rho:+.4f}   "
      f"(adopt <= 0.25, abandon >= 0.40)")
print(f"  ==> {prim.verdict}")

v = set(T[(T.tag.str.startswith("AFTER")) & (T.set == "valid genes")].verdict)
print(f"\n  verdicts across the three CN thresholds: {v}")
print("  (prespec: if these differ, the result is INCONCLUSIVE regardless)")

# ---- pre-registered prediction 1: assessed fraction on retained pairs ----------
print("\n" + "=" * 118)
print("Pre-registered prediction check: cna_gate assessed-fraction on retained pairs")
print("=" * 118)
old_ass = float((C_old["cna_gate"].notna() & retained).to_numpy().sum()
                / retained.to_numpy().sum())
new_ass = float((cn_check(DEL_PRIMARY).notna() & retained).to_numpy().sum()
                / retained.to_numpy().sum())
print(f"  BEFORE (COSMIC)         : {100*old_ass:.2f}% of retained pairs assessed")
print(f"  AFTER  (gene-complete)  : {100*new_ass:.2f}%   (predicted >55%)")
print(f"  prediction {'HELD' if new_ass > 0.55 else 'FAILED'}")

print("\n  per-check state on retained pairs, AFTER:")
C_new = dict(C_old); C_new["cna_gate"] = cn_check(DEL_PRIMARY)
n = int(retained.to_numpy().sum())
st = []
for chk in RC.CHECKS:
    if chk == "rna_gate":
        continue
    x = C_new[chk].where(retained)
    st.append({"check": chk,
               "passed": round(float((x == 1).to_numpy().sum()) / n, 4),
               "failed": round(float((x == 0).to_numpy().sum()) / n, 4),
               "not_assessed": round(float((x.isna() & retained).to_numpy().sum()) / n, 4)})
print(pd.DataFrame(st).to_string(index=False))

others = [c for c in RC.CHECKS if c != "rna_gate"]
nass_new = sum(C_new[c].notna().astype(float) for c in others).where(retained)
b = nass_new.to_numpy().ravel(); b = b[np.isfinite(b)]
print("\n  number of non-gate checks assessed per retained pair, AFTER:")
for k in range(5):
    print(f"    {k}: {100*(b == k).mean():6.2f}%")

T.to_json(os.path.join(W, "partA_v2_results.json"), orient="records", indent=1)
print("\nwritten partA_v2_results.json")
