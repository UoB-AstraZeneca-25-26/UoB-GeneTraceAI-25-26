"""PART A -- can a confirmation count discriminate at all?

A1 coverage confound (decisive)
A2 tie structure
A3 between-line vs between-(gene,line) variance

Standing rule 3: every measurement runs twice -- ALL genes, and VALID genes only
(valid = passes the silence guard, expressed in >=20% of the panel).
"""
import os, json, numpy as np, pandas as pd
from scipy import stats
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
W = os.path.dirname(os.path.abspath(__file__))
import rank_checks as RC

con = RC.connect()
genes, sym = RC.sample_genes(con, n=300, seed=42)
print(f"sampled {len(genes)} protein-coding genes (seed 42)")
B = RC.build(con, genes, sym)
C, idx = B["checks"], B["lines"]
print(f"lines: {len(idx):,}")
npass, nass = RC.confirmation_count(C, genes, idx)
valid = B["valid_gene"]
print(f"genes passing the silence guard: {int(valid.sum())} of {len(genes)} "
      f"({100*valid.mean():.1f}%)")

# per-check state fractions (feeds B3 as well)
print("\n=== per-check state fractions over all (gene, line) pairs ===")
rows = []
for chk in RC.CHECKS:
    v = C[chk].reindex(index=idx, columns=genes)
    n = v.size
    rows.append({"check": chk,
                 "passed": round(float((v == 1).sum().sum()) / n, 4),
                 "failed": round(float((v == 0).sum().sum()) / n, 4),
                 "not_assessed": round(float(v.isna().sum().sum()) / n, 4)})
S = pd.DataFrame(rows)
print(S.to_string(index=False))

# ------------------------------------------------------------------ A1
print("\n" + "=" * 78)
print("A1 -- THE COVERAGE CONFOUND (decisive)")
print("=" * 78)
cc = B["cov_count"]
print("gene-independent evidence-type availability per line (0-5):")
print(cc.value_counts().sort_index().to_string())

lin = pd.read_parquet("src/pipeline/outputs/gdsc_models.parquet",
                      columns=["model_id", "tissue"]).dropna(subset=["model_id"])
lin = lin.drop_duplicates("model_id").set_index("model_id")["tissue"].reindex(idx)

a1 = {}
for label, gset in [("all genes", genes), ("valid genes only", [g for g in genes if valid[g]])]:
    npc = npass[gset]
    # per-gene correlation between coverage_count and confirmation_count
    rho = []
    for g in gset:
        y = npc[g]
        ok = y.notna() & cc.notna()
        if ok.sum() > 30:
            rho.append(stats.spearmanr(cc[ok], y[ok]).statistic)
    rho = np.array([r for r in rho if np.isfinite(r)])
    # pooled over all pairs
    flat_y = npc.to_numpy().ravel()
    flat_c = np.repeat(cc.to_numpy()[:, None], len(gset), axis=1).ravel()
    m = np.isfinite(flat_y) & np.isfinite(flat_c)
    pooled = stats.spearmanr(flat_c[m], flat_y[m]).statistic
    # within lineage
    wl = []
    for t, sub in lin.dropna().groupby(lin.dropna()):
        li = sub.index
        if len(li) < 30:
            continue
        yy = npc.loc[li].to_numpy().ravel()
        ccx = np.repeat(cc.loc[li].to_numpy()[:, None], len(gset), axis=1).ravel()
        mm = np.isfinite(yy) & np.isfinite(ccx)
        if mm.sum() > 100:
            r = stats.spearmanr(ccx[mm], yy[mm]).statistic
            if np.isfinite(r):
                wl.append(r)
    a1[label] = {"n_genes": len(gset), "pooled_spearman": round(float(pooled), 4),
                 "per_gene_median": round(float(np.median(rho)), 4),
                 "per_gene_q25": round(float(np.percentile(rho, 25)), 4),
                 "per_gene_q75": round(float(np.percentile(rho, 75)), 4),
                 "within_lineage_median": round(float(np.median(wl)), 4) if wl else None,
                 "n_lineages": len(wl)}
    print(f"\n  [{label}]  n_genes={len(gset)}")
    print(f"    pooled Spearman(coverage_count, confirmation_count) = {pooled:+.4f}")
    print(f"    per-gene Spearman: median {np.median(rho):+.4f}  "
          f"IQR {np.percentile(rho,25):+.4f} to {np.percentile(rho,75):+.4f}")
    if wl:
        print(f"    within-lineage pooled Spearman: median {np.median(wl):+.4f} "
              f"over {len(wl)} lineages")

# does coverage correlate with things that make this actively misleading?
print("\n  -- is coverage itself confounded? --")
rna_mean = B["RNA"].mean(axis=1)
ok = cc.notna() & rna_mean.notna()
r_rna = stats.spearmanr(cc[ok], rna_mean[ok])
print(f"    Spearman(coverage_count, mean RNA per line) = {r_rna.statistic:+.4f} "
      f"p={r_rna.pvalue:.2e}  n={int(ok.sum())}")
kw = [cc[lin == t].dropna() for t in lin.dropna().unique() if (lin == t).sum() >= 15]
H = stats.kruskal(*[k for k in kw if len(k) > 2])
print(f"    coverage_count vs lineage: Kruskal-Wallis H={H.statistic:.1f} p={H.pvalue:.3e} "
      f"({len(kw)} lineages)")
byl = pd.DataFrame({"cov": cc, "tis": lin}).dropna().groupby("tis")["cov"].agg(["size", "mean"])
byl = byl[byl["size"] >= 15].sort_values("mean")
print(f"    per-lineage mean coverage: {byl['mean'].min():.2f} ({byl.index[0]}) to "
      f"{byl['mean'].max():.2f} ({byl.index[-1]})")
reg = pd.read_parquet("src/pipeline/outputs/gene_regime.parquet", columns=["ensg_id"])
cur = set(reg.ensg_id.str.lower())
print(f"    curated-gene-set membership: {len(set(genes)&cur)} of {len(genes)} sampled "
      f"genes are in the 141-gene curated set")

# ------------------------------------------------------------------ A2
print("\n" + "=" * 78)
print("A2 -- THE TIE STRUCTURE")
print("=" * 78)
a2 = {}
for label, gset in [("all genes", genes), ("valid genes only", [g for g in genes if valid[g]])]:
    npc = npass[gset]
    lg, nd, tot = [], [], []
    for g in gset:
        y = npc[g].dropna()
        if len(y) < 50:
            continue
        vc = y.value_counts()
        lg.append(vc.max() / len(y)); nd.append(len(vc)); tot.append(len(y))
    a2[label] = {"n_genes": len(lg),
                 "median_largest_tie_frac": round(float(np.median(lg)), 4),
                 "q25": round(float(np.percentile(lg, 25)), 4),
                 "q75": round(float(np.percentile(lg, 75)), 4),
                 "median_distinct_values": float(np.median(nd)),
                 "median_lines_per_gene": float(np.median(tot))}
    print(f"\n  [{label}]  genes={len(lg)}")
    print(f"    largest tied group as a fraction of lines: median "
          f"{np.median(lg):.3f}  IQR {np.percentile(lg,25):.3f}-{np.percentile(lg,75):.3f}")
    print(f"    distinct count values achieved: median {np.median(nd):.0f} of 6 possible")
    print(f"    lines per gene: median {np.median(tot):.0f}")

print("\n  distribution of confirmation counts, pooled over all (gene, line):")
allv = npass[genes].to_numpy().ravel()
allv = allv[np.isfinite(allv)]
for k in range(6):
    print(f"    count = {k}: {int((allv==k).sum()):>9,}  ({100*(allv==k).mean():5.2f}%)")

# ------------------------------------------------------------------ A3
print("\n" + "=" * 78)
print("A3 -- BETWEEN-LINE vs BETWEEN-(GENE,LINE) VARIANCE")
print("=" * 78)
a3 = {}
for label, gset in [("all genes", genes), ("valid genes only", [g for g in genes if valid[g]])]:
    M = npass[gset]
    v = M.to_numpy(float)
    ok = np.isfinite(v)
    grand = v[ok].mean()
    line_mean = np.nanmean(np.where(ok, v, np.nan), axis=1)
    n_per_line = ok.sum(axis=1)
    good = n_per_line > 10
    ss_between = float((n_per_line[good] * (line_mean[good] - grand) ** 2).sum())
    ss_total = float(((v[ok] - grand) ** 2).sum())
    icc = ss_between / ss_total if ss_total else np.nan
    # rank agreement of per-gene orderings against the line-mean ordering
    ref = pd.Series(line_mean, index=M.index)
    cors = []
    for g in gset:
        y = M[g]
        o = y.notna() & ref.notna()
        if o.sum() > 50:
            r = stats.spearmanr(y[o], ref[o]).statistic
            if np.isfinite(r):
                cors.append(r)
    a3[label] = {"variance_between_lines_frac": round(icc, 4),
                 "median_rho_to_line_mean": round(float(np.median(cors)), 4),
                 "q25": round(float(np.percentile(cors, 25)), 4),
                 "q75": round(float(np.percentile(cors, 75)), 4),
                 "n_genes": len(cors)}
    print(f"\n  [{label}]")
    print(f"    fraction of variance BETWEEN LINES: {icc:.4f}")
    print(f"    per-gene ordering vs the line-mean ordering: Spearman median "
          f"{np.median(cors):+.4f}  IQR {np.percentile(cors,25):+.4f}-{np.percentile(cors,75):+.4f}")

# gene-to-gene ordering agreement: would the user get the same answer every time?
print("\n  -- would two different genes return the same ordering? --")
for label, gset in [("all genes", genes), ("valid genes only", [g for g in genes if valid[g]])]:
    rng = np.random.default_rng(0)
    pairs = rng.choice(len(gset), size=(300, 2))
    rr = []
    for i, j in pairs:
        if i == j:
            continue
        a, b = npass[gset[i]], npass[gset[j]]
        o = a.notna() & b.notna()
        if o.sum() > 50:
            r = stats.spearmanr(a[o], b[o]).statistic
            if np.isfinite(r):
                rr.append(r)
    print(f"    [{label}] Spearman between two random genes' orderings: "
          f"median {np.median(rr):+.4f}  IQR {np.percentile(rr,25):+.4f}-"
          f"{np.percentile(rr,75):+.4f}  (n={len(rr)} pairs)")
    a3[label]["median_rho_between_genes"] = round(float(np.median(rr)), 4)

json.dump({"per_check_states": S.to_dict("records"), "A1": a1, "A2": a2, "A3": a3,
           "cov_count_dist": {int(k): int(v) for k, v in cc.value_counts().items()},
           "n_genes": len(genes), "n_valid_genes": int(valid.sum()), "n_lines": len(idx),
           "cov_vs_meanRNA_spearman": float(r_rna.statistic)},
          open(os.path.join(W, "partA_results.json"), "w"), indent=1)
npass.to_parquet(os.path.join(W, "partA_npass.parquet"))
nass.to_parquet(os.path.join(W, "partA_nassessed.parquet"))
for chk in RC.CHECKS:
    C[chk].reindex(index=idx, columns=genes).to_parquet(os.path.join(W, f"partA_check_{chk}.parquet"))
print("\nwritten partA_results.json")
