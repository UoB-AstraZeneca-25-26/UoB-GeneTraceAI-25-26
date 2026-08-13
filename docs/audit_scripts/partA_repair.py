"""Part A, continued: does the obvious repair rescue the count?

B3 offers three ways to handle unmeasured checks. Test each against the A1
confound before declaring the design dead:
  (i)   raw count                     n_pass
  (ii)  fraction of assessable        n_pass / n_assessed
  (iii) minimum-assessed gate         n_pass restricted to n_assessed >= m
Also asks the sharper question: is the count just the RNA gate wearing four
other labels?
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
res = json.load(open(os.path.join(W, "partA_results.json")))

con = RC.connect()
cov = pd.read_parquet("src/pipeline/outputs/coverage_matrix_enriched.parquet").set_index("model_id")
cc = pd.Series(0.0, index=npass.index)
for chk, cols in RC.COVERAGE_COLS.items():
    cc += cov.reindex(npass.index)[cols].any(axis=1).fillna(False).astype(float)

# validity flag (recompute cheaply from the rna_gate check: it is 100% assessed)
frac_expr = C["rna_gate"].mean(axis=0)
valid = frac_expr >= RC.SILENT_FRAC
gsets = {"all genes": genes, "valid genes only": [g for g in genes if valid[g]]}

frac = (npass / nass.replace(0, np.nan))

print("=" * 78)
print("Does the repair remove the coverage confound?")
print("=" * 78)
out = {}
for label, gs in gsets.items():
    print(f"\n[{label}]  n_genes={len(gs)}")
    row = {}
    for name, M in [("(i)   raw count            n_pass", npass[gs]),
                    ("(ii)  fraction assessable  n_pass/n_assessed", frac[gs])]:
        a = M.to_numpy().ravel()
        b = np.repeat(cc.to_numpy()[:, None], len(gs), axis=1).ravel()
        m = np.isfinite(a) & np.isfinite(b)
        r = stats.spearmanr(b[m], a[m]).statistic
        # per-gene
        pg = []
        for g in gs:
            y = M[g]; o = y.notna() & cc.notna()
            if o.sum() > 30:
                rr = stats.spearmanr(cc[o], y[o]).statistic
                if np.isfinite(rr):
                    pg.append(rr)
        print(f"  {name:<48} pooled rho={r:+.4f}   per-gene median={np.median(pg):+.4f}")
        row[name.strip()] = {"pooled": round(float(r), 4),
                             "per_gene_median": round(float(np.median(pg)), 4)}
    # (iii) minimum-assessed gate
    for m_min in (2, 3, 4):
        keep = nass[gs] >= m_min
        M = npass[gs].where(keep)
        a = M.to_numpy().ravel()
        b = np.repeat(cc.to_numpy()[:, None], len(gs), axis=1).ravel()
        mk = np.isfinite(a) & np.isfinite(b)
        r = stats.spearmanr(b[mk], a[mk]).statistic
        retained = float(keep.to_numpy().mean())
        print(f"  (iii) min_assessed >= {m_min}   pooled rho={r:+.4f}   "
              f"retains {100*retained:.1f}% of pairs")
        row[f"min_assessed_{m_min}"] = {"pooled": round(float(r), 4),
                                        "frac_pairs_retained": round(retained, 4)}
    out[label] = row

print("\n" + "=" * 78)
print("Is the count just the RNA gate wearing four other labels?")
print("=" * 78)
for label, gs in gsets.items():
    a = npass[gs].to_numpy().ravel()
    r = C["rna_gate"][gs].to_numpy().ravel()
    f = frac[gs].to_numpy().ravel()
    m = np.isfinite(a) & np.isfinite(r)
    print(f"\n[{label}]")
    print(f"  Spearman(n_pass, rna_gate alone)          = "
          f"{stats.spearmanr(r[m], a[m]).statistic:+.4f}")
    mf = np.isfinite(f) & np.isfinite(r)
    print(f"  Spearman(fraction, rna_gate alone)        = "
          f"{stats.spearmanr(r[mf], f[mf]).statistic:+.4f}")
    # how much does the count add beyond the RNA gate? tie structure within gate state
    sub = pd.DataFrame({"rna": r[m], "n": a[m]})
    print("  distribution of n_pass within each rna_gate state:")
    print(sub.groupby("rna")["n"].agg(["size", "mean", "std"]).round(3).to_string())

print("\n" + "=" * 78)
print("Tie structure of the REPAIRED (fraction) score")
print("=" * 78)
for label, gs in gsets.items():
    lg, nd = [], []
    for g in gs:
        y = frac[g].dropna()
        if len(y) < 50:
            continue
        vc = y.value_counts()
        lg.append(vc.max() / len(y)); nd.append(len(vc))
    print(f"  [{label}] largest tied group: median {np.median(lg):.3f}  "
          f"IQR {np.percentile(lg,25):.3f}-{np.percentile(lg,75):.3f};  "
          f"distinct values median {np.median(nd):.0f}")
    out[label]["fraction_largest_tie"] = round(float(np.median(lg)), 4)
    out[label]["fraction_distinct_values"] = float(np.median(nd))

print("\n" + "=" * 78)
print("Between-gene ordering agreement for the REPAIRED score")
print("=" * 78)
for label, gs in gsets.items():
    rng = np.random.default_rng(0)
    rr = []
    for i, j in rng.choice(len(gs), size=(300, 2)):
        if i == j:
            continue
        a, b = frac[gs[i]], frac[gs[j]]
        o = a.notna() & b.notna()
        if o.sum() > 50:
            x = stats.spearmanr(a[o], b[o]).statistic
            if np.isfinite(x):
                rr.append(x)
    print(f"  [{label}] Spearman between two random genes' orderings: "
          f"median {np.median(rr):+.4f}  IQR {np.percentile(rr,25):+.4f}-{np.percentile(rr,75):+.4f}")
    out[label]["fraction_between_gene_rho"] = round(float(np.median(rr)), 4)

json.dump(out, open(os.path.join(W, "partA_repair_results.json"), "w"), indent=1)
print("\nwritten partA_repair_results.json")
