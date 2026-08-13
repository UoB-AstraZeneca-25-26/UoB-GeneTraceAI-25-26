"""B1.3 / B1.4 -- the two metric definitions, computed separately, on the
113 tuning genes and the 28 held-out genes, with hypergeometric chance
baselines and a paired bootstrap over genes.

DEF-A  recall@k      = |top-k inter sensitive| / |sensitive|      (what ALL code computes)
DEF-B  any-hit@k     = 1[ |top-k inter sensitive| >= 1 ]          (what the .docx reports)

Chance baselines, per gene, with N scored lines and K sensitive lines:
  E[recall@k]   = k/N
  P[any-hit@k]  = 1 - C(N-K, k)/C(N, k)      (hypergeometric)
"""
import os, json, numpy as np, pandas as pd
from scipy.special import gammaln
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")

OUT = "src/pipeline/outputs"
core = pd.read_parquet(f"{OUT}/core_score.parquet", columns=["model_id", "ensg_id", "core_score"])
regime = pd.read_parquet(f"{OUT}/gene_regime.parquet")
gdsc = pd.read_parquet("validation/prepared/gdsc_scored_ready.parquet",
                       columns=["model_id", "target_ensg", "sensitive"])
core["model_id"] = core["model_id"].str.lower()
gdsc["model_id"] = gdsc["model_id"].str.lower()
gdsc["target_ensg"] = gdsc["target_ensg"].astype("string").str.lower()

# EXACT split reproduction from 06_held_out_eval.ipynb cell 3
curated = sorted(set(regime.ensg_id) & set(gdsc.target_ensg))
rng = np.random.default_rng(42)
test_genes = set(rng.choice(curated, size=max(1, len(curated) // 5), replace=False))
train_genes = set(g for g in curated if g not in test_genes)
print(f"curated={len(curated)}  train={len(train_genes)}  test={len(test_genes)}")
shipped = set(json.load(open(f"{OUT}/stage4_eval_summary.json"))["test_genes"])
print(f"split reproduces shipped stage4 test set: {test_genes == shipped}")

sens = (gdsc[gdsc.sensitive][["target_ensg", "model_id"]]
        .drop_duplicates().rename(columns={"target_ensg": "ensg_id"}).assign(is_sensitive=True))
c = core[core.ensg_id.isin(curated)].merge(sens, on=["ensg_id", "model_id"], how="left")
c["is_sensitive"] = c["is_sensitive"].fillna(False)
c["rank"] = c.groupby("ensg_id")["core_score"].rank(ascending=False, method="first")


def logC(n, k):
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


rows = []
for g, d in c.groupby("ensg_id"):
    N = len(d); K = int(d.is_sensitive.sum())
    if K == 0:
        continue
    r = {"ensg_id": g, "N": N, "K": K, "split": "test" if g in test_genes else "train"}
    for k in (1, 5, 20):
        hits = int(d.loc[d["rank"] <= k, "is_sensitive"].sum())
        r[f"recall@{k}"] = hits / K
        r[f"anyhit@{k}"] = float(hits >= 1)
        r[f"base_recall@{k}"] = min(k, N) / N
        r[f"base_anyhit@{k}"] = 1.0 - (np.exp(logC(N - K, k) - logC(N, k)) if N - K >= k else 0.0)
    rows.append(r)
per = pd.DataFrame(rows)
per.to_parquet(os.path.join(os.path.dirname(__file__), "b14_per_gene.parquet"), index=False)


def boot(vals, base, n=10000, seed=42):
    """Paired bootstrap over genes of mean(metric - chance baseline)."""
    rg = np.random.default_rng(seed)
    d = np.asarray(vals) - np.asarray(base)
    idx = rg.integers(0, len(d), size=(n, len(d)))
    m = d[idx].mean(axis=1)
    return np.percentile(m, [2.5, 97.5])


print("\n" + "=" * 96)
print(f"{'split':6s} {'n':>4s} {'metric':12s} {'observed':>9s} {'chance':>9s} {'excess':>9s}  {'95% CI on excess (paired boot)'}")
print("=" * 96)
for split in ("train", "test"):
    s = per[per.split == split]
    for k in (1, 5, 20):
        for name in ("recall", "anyhit"):
            o = s[f"{name}@{k}"].mean(); b = s[f"base_{name}@{k}"].mean()
            lo, hi = boot(s[f"{name}@{k}"], s[f"base_{name}@{k}"])
            star = "" if (lo <= 0 <= hi) else "  *"
            print(f"{split:6s} {len(s):4d} {name+'@'+str(k):12s} {o:9.4f} {b:9.4f} {o-b:+9.4f}  [{lo:+.4f}, {hi:+.4f}]{star}")
    print("-" * 96)
print("* = 95% CI on (observed - chance) excludes zero")
print(f"\nmedian N (scored lines/gene) = {per.N.median():.0f}; median K (sensitive/gene) = {per.K.median():.0f}")
print("\n.docx claims: top-1 23.4%  top-5 65.2%  top-20 93.6% over 'the 141 curated genes'")
allg = per
for k in (1, 5, 20):
    print(f"  all {len(allg)} curated genes, anyhit@{k} = {allg[f'anyhit@{k}'].mean()*100:.1f}%"
          f"   recall@{k} = {allg[f'recall@{k}'].mean()*100:.1f}%")
