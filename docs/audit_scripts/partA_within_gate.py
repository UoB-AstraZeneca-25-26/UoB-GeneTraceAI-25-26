"""Part A, the fair test: the count is proposed as LEVEL 2, applied WITHIN the
gate-retained set (Level 1). Measure it there, not over all lines.

Retained := passed the RNA abundance gate. That is the only check that is 100%
assessed and materially discriminating, so it defines Level 1 here.
"""
import os, json, numpy as np, pandas as pd
from scipy import stats
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
W = os.path.dirname(os.path.abspath(__file__))
import rank_checks as RC

npass = pd.read_parquet(os.path.join(W, "partA_npass.parquet"))
nass = pd.read_parquet(os.path.join(W, "partA_nassessed.parquet"))
genes = list(npass.columns)
C = {c: pd.read_parquet(os.path.join(W, f"partA_check_{c}.parquet")) for c in RC.CHECKS}

cov = pd.read_parquet("src/pipeline/outputs/coverage_matrix_enriched.parquet").set_index("model_id")
cc = pd.Series(0.0, index=npass.index)
for chk, cols in RC.COVERAGE_COLS.items():
    cc += cov.reindex(npass.index)[cols].any(axis=1).fillna(False).astype(float)

retained = C["rna_gate"] == 1.0            # Level 1
frac_expr = C["rna_gate"].mean(axis=0)
valid = frac_expr >= RC.SILENT_FRAC
gsets = {"all genes": genes, "valid genes only": [g for g in genes if valid[g]]}

# Level 2 candidates, computed on the RETAINED set only, EXCLUDING the gate check
# itself (a check cannot confirm the tier it defined).
sub_pass = npass - C["rna_gate"].fillna(0.0)          # drop rna_gate from the count
sub_ass = nass - C["rna_gate"].notna().astype(float)
raw2 = sub_pass.where(retained)
frac2 = (sub_pass / sub_ass.replace(0, np.nan)).where(retained)

print("=" * 78)
print("LEVEL 2 measured WITHIN the gate-retained set (gate check excluded)")
print("=" * 78)
print(f"retained pairs: {int(retained.to_numpy().sum()):,} of {retained.size:,} "
      f"({100*retained.to_numpy().mean():.1f}%)")
print(f"median retained lines per gene: {retained.sum(axis=0).median():.0f}")

out = {}
for label, gs in gsets.items():
    print(f"\n[{label}]  n_genes={len(gs)}")
    row = {}
    for name, M in [("raw count (4 checks)", raw2[gs]), ("fraction assessable", frac2[gs])]:
        # A1 confound
        a = M.to_numpy().ravel()
        b = np.repeat(cc.to_numpy()[:, None], len(gs), axis=1).ravel()
        m = np.isfinite(a) & np.isfinite(b)
        rho_cov = stats.spearmanr(b[m], a[m]).statistic
        # A2 ties
        lg, nd, nret = [], [], []
        for g in gs:
            y = M[g].dropna()
            if len(y) < 50:
                continue
            vc = y.value_counts()
            lg.append(vc.max() / len(y)); nd.append(len(vc)); nret.append(len(y))
        # A3 between-gene agreement
        rng = np.random.default_rng(0); rr = []
        for i, j in rng.choice(len(gs), size=(300, 2)):
            if i == j:
                continue
            x, y = M[gs[i]], M[gs[j]]
            o = x.notna() & y.notna()
            if o.sum() > 50:
                v = stats.spearmanr(x[o], y[o]).statistic
                if np.isfinite(v):
                    rr.append(v)
        print(f"  {name:<22} coverage rho={rho_cov:+.4f}   "
              f"largest tie={np.median(lg):.3f}   distinct={np.median(nd):.0f}   "
              f"between-gene rho={np.median(rr):+.4f}")
        row[name] = {"coverage_rho": round(float(rho_cov), 4),
                     "largest_tie_frac": round(float(np.median(lg)), 4),
                     "distinct_values": float(np.median(nd)),
                     "between_gene_rho": round(float(np.median(rr)), 4),
                     "median_retained_lines": float(np.median(nret))}
    out[label] = row

# what does the level-2 count actually consist of, on retained pairs?
print("\n" + "=" * 78)
print("What is left to count, on retained pairs?")
print("=" * 78)
rows = []
for chk in RC.CHECKS:
    if chk == "rna_gate":
        continue
    v = C[chk].where(retained)
    n = int(retained.to_numpy().sum())
    rows.append({"check": chk,
                 "passed": round(float((v == 1).to_numpy().sum()) / n, 4),
                 "failed": round(float((v == 0).to_numpy().sum()) / n, 4),
                 "not_assessed": round(float((v.isna() & retained).to_numpy().sum()) / n, 4)})
print(pd.DataFrame(rows).to_string(index=False))

print("\ndistribution of the level-2 raw count on retained pairs:")
a = raw2.to_numpy().ravel(); a = a[np.isfinite(a)]
for k in range(5):
    print(f"  {k}: {int((a==k).sum()):>9,}  ({100*(a==k).mean():5.2f}%)")

print("\ndistribution of n_assessed (excluding the gate) on retained pairs:")
b = sub_ass.where(retained).to_numpy().ravel(); b = b[np.isfinite(b)]
for k in range(5):
    print(f"  {k}: {int((b==k).sum()):>9,}  ({100*(b==k).mean():5.2f}%)")

json.dump(out, open(os.path.join(W, "partA_within_gate_results.json"), "w"), indent=1)
print("\nwritten partA_within_gate_results.json")
