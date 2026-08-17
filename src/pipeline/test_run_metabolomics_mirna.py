"""
test_run_metabolomics_mirna.py
--------------------------------
Read-only diagnostic. Does NOT touch core_score.parquet or any production
output -- answers one question: if we added a curated metabolomics enzyme
layer and a curated miRNA-inhibitor-of-expression layer, would known genes
score better or worse?

Two separate tests, because the two candidate signals have different
available ground truth (see findings below):

  A) miRNA-as-expression-inhibitor, tested on BCL2 + MET only.
     These are the only two of the curated-mapping genes that have GDSC
     drug-sensitivity ground truth in validation/prepared/curated_validation_pairs.parquet
     (the same table build_full_predictions.py's validation harness uses).
     Metric: hit@20, same definition as validate_chronos.py.

  B) Metabolomics-as-enzyme-activity, tested on the 9 curated enzyme genes
     against Chronos CRISPR essentiality (validation/prepared/chronos_long.parquet),
     because NONE of them have GDSC drug coverage -- same substitute ground
     truth build_chronos_validation.py uses for genes without a curated drug.
     Metric: Spearman rho between the metabolomics proxy and real essentiality,
     same validated/inverted/none classification as build_chronos_validation.py.

Curated mappings (both literature-sourced, see comments below) -- built
because neither metabolomics nor miRNA has a native gene axis (as noted in
02_core_score.ipynb's own scope statement).
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, t as tdist
import sys
sys.path.insert(0, "src/pipeline")
from scoring_variants import combine

REF = Path("reference")
CLEANED = Path("cleaned_track_data")
DATA_CLEAN = Path("data/parquet/data_clean")
OUT = Path("src/pipeline/outputs")
VAL = Path("validation/prepared")

RHO_THRESHOLD = 0.1
P_THRESHOLD = 0.05

# ── Curated gene -> metabolite enzyme map ────────────────────────────────
# Only genes in the 141-curated set that are genuine metabolic enzymes with
# a substrate/product pair present in the 227-metabolite CCLE panel.
# proxy sign convention: proxy = product_pct - substrate_pct (percentile
# ranks across cell lines) -- higher proxy = more flux through that enzyme.
METABOLITE_MAP = {
    # NAD+ salvage pathway, rate-limiting step (Revollo et al. 2007)
    "NAMPT": {"substrate": "niacinamide", "product": "nad"},
    # 2-oxoglutarate-dependent dioxygenase (prolyl hydroxylase), HIF regulation
    # (Kaelin & Ratcliffe 2008) -- alpha-KG consumed, succinate produced
    "EGLN1": {"substrate": "alpha-ketoglutarate", "product": "succinate/methylmalonate"},
    # JmjC-domain histone demethylases -- same Fe(II)/2-OG oxygenase family as EGLN1
    "KDM3A": {"substrate": "alpha-ketoglutarate", "product": "succinate/methylmalonate"},
    "KDM4A": {"substrate": "alpha-ketoglutarate", "product": "succinate/methylmalonate"},
    "KDM4C": {"substrate": "alpha-ketoglutarate", "product": "succinate/methylmalonate"},
    "KDM4E": {"substrate": "alpha-ketoglutarate", "product": "succinate/methylmalonate"},
    "KDM6B": {"substrate": "alpha-ketoglutarate", "product": "succinate/methylmalonate"},
    # PARylation consumes NAD+, releases nicotinamide (opposite direction to NAMPT)
    "PARP1": {"substrate": "nad", "product": "niacinamide"},
    "PARP2": {"substrate": "nad", "product": "niacinamide"},
}

# ── Curated gene -> validated miRNA regulator map ────────────────────────
# Restricted to well-cited, experimentally validated pairs (not predicted).
MIRNA_MAP = {
    "BCL2": ["hsa-mir-34a", "hsa-mir-15a", "hsa-mir-16"],   # Cimmino 2005 PNAS; He 2007 Nature
    "MET":  ["hsa-mir-1"],                                    # Nasser 2008 JBC (rhabdomyosarcoma/HCC)
}

GENE_SYMBOL_TO_ENSG = {
    "NAMPT": "ensg00000105835", "EGLN1": "ensg00000135766",
    "KDM3A": "ensg00000115548", "KDM4A": "ensg00000066135",
    "KDM4C": "ensg00000107077", "KDM4E": "ensg00000235268",
    "KDM6B": "ensg00000132510", "PARP1": "ensg00000143799",
    "PARP2": "ensg00000129484", "BCL2": "ensg00000171791",
    "MET": "ensg00000105976",
}


def percentile_rank(df):
    return df.rank(axis=0, pct=True, na_option="keep")


def dedup_mean(df):
    if df.index.duplicated(keep=False).sum():
        df = df.groupby(level=0).mean()
    return df


def spearman_with_p(x, y):
    rho, p = spearmanr(x, y)
    return rho, p


def classify(rho, p):
    if rho > RHO_THRESHOLD and p < P_THRESHOLD:
        return "validated"
    if rho < -RHO_THRESHOLD and p < P_THRESHOLD:
        return "inverted"
    return "none"


def hit_at_k(scores_df, ensg_id, sensitive_model_ids, k=20, score_col="core_score"):
    gene_rows = scores_df[scores_df.ensg_id == ensg_id].copy()
    if gene_rows.empty or len(sensitive_model_ids) == 0:
        return np.nan
    gene_rows = gene_rows.sort_values(score_col, ascending=False).reset_index(drop=True)
    top_k = set(gene_rows.head(k)["model_id"])
    return len(top_k & sensitive_model_ids) / len(sensitive_model_ids)


print("=" * 70)
print("PART A -- miRNA-as-expression-inhibitor (BCL2, MET)")
print("=" * 70)

genes_a = list(MIRNA_MAP.keys())
ensg_a = [GENE_SYMBOL_TO_ENSG[g] for g in genes_a]

# -- rebuild expression, restricted to the 2 target genes (mirrors 02_core_score.ipynb) --
prof = pd.read_parquet(REF / "depmap_profiles.parquet")
rna_prof = prof[prof["datatype"] == "rna"][["profileid", "modelid"]].copy()
rna_prof["model_id"] = rna_prof["modelid"].str.lower()
rna_prof = rna_prof.drop(columns=["modelid"])

expr_raw = pd.read_parquet(DATA_CLEAN / "depmap_expr_clean.parquet")
expr_raw.index.name = "profileid"
expr_raw.columns = expr_raw.columns.str.lower()
keep_cols = [c for c in expr_raw.columns if c in ensg_a]
expr_raw = expr_raw[keep_cols]
expr = expr_raw.reset_index().merge(rna_prof, on="profileid", how="inner")
expr = expr.drop(columns=["profileid"]).set_index("model_id")
expr = dedup_mean(expr)
print(f"expr (raw log2 TPM) restricted shape: {expr.shape}")

# -- rebuild proteomics, restricted to the 2 target genes --
gene_lookup = pd.read_parquet(REF / "gene_lookup.parquet")
universe = gene_lookup[(gene_lookup["biotype"] == "protein_coding") &
                        (gene_lookup["hgnc_status"] == "Approved")][["ensg_id", "uniprot_ids"]].copy()
universe["ensg_id"] = universe["ensg_id"].astype("string").str.split(".").str[0].str.lower()
uniprot_to_ensg = (universe.dropna(subset=["uniprot_ids"])
                   .assign(uniprot_id=lambda x: x["uniprot_ids"].str.strip().str.lower())
                   .drop_duplicates(subset=["uniprot_id"])[["uniprot_id", "ensg_id"]])

prot_raw = pd.read_parquet(CLEANED / "proteomics.parquet")
key_col = "model_id" if "model_id" in prot_raw.columns else "depmap_id"
prot_long = prot_raw.melt(id_vars=[key_col], var_name="uniprot_id", value_name="log_ratio").dropna(subset=["log_ratio"])
prot_long["model_id"] = prot_long[key_col].str.lower()
prot_mapped = prot_long.merge(uniprot_to_ensg, on="uniprot_id", how="inner")
prot_mapped = prot_mapped[prot_mapped["ensg_id"].isin(ensg_a)]
prot = prot_mapped.pivot_table(index="model_id", columns="ensg_id", values="log_ratio", aggfunc="mean")
print(f"prot restricted shape: {prot.shape}")

E_pct_baseline = percentile_rank(expr)
P_pct = percentile_rank(prot)

# -- baseline candidate score, matched recipe (percentile + weighted_sum 0.5/0.5) --
all_models = sorted(set(E_pct_baseline.index) | set(P_pct.index))
all_genes = sorted(set(E_pct_baseline.columns) | set(P_pct.columns))
E_al = E_pct_baseline.reindex(index=all_models, columns=all_genes)
P_al = P_pct.reindex(index=all_models, columns=all_genes)
core_baseline, n_layers = combine(E_al, P_al, "weighted_sum")
core_baseline_long = core_baseline.stack().reset_index()
core_baseline_long.columns = ["model_id", "ensg_id", "core_score"]

# -- miRNA inhibitor: per gene, per model, percentile of curated miRNA expression --
mirna = pd.read_parquet(CLEANED / "mirna_model_level.parquet")
mirna["model_id"] = mirna["model_id"].str.lower()

ALPHA = 1.0  # full dampening at max miRNA percentile
E_eff = expr.copy()
for gene, mirs in MIRNA_MAP.items():
    ensg = GENE_SYMBOL_TO_ENSG[gene]
    if ensg not in E_eff.columns:
        continue
    sub = mirna[mirna["mirna_symbol"].isin(mirs)]
    burden = sub.groupby("model_id")["expression_value"].mean()
    burden_pct = burden.rank(pct=True)
    burden_pct = burden_pct.reindex(E_eff.index)
    dampener = 1.0 - ALPHA * burden_pct.fillna(0.0)
    E_eff[ensg] = E_eff[ensg] * dampener
    print(f"  {gene}: miRNA burden coverage {burden_pct.notna().sum()}/{len(E_eff)} cell lines, "
          f"mean dampener {dampener.mean():.3f}")

E_pct_inhibited = percentile_rank(E_eff)
E_al2 = E_pct_inhibited.reindex(index=all_models, columns=all_genes)
core_candidate, _ = combine(E_al2, P_al, "weighted_sum")
core_candidate_long = core_candidate.stack().reset_index()
core_candidate_long.columns = ["model_id", "ensg_id", "core_score"]

# -- production reference (exact current pipeline output, for sanity anchoring) --
prod = pd.read_parquet(OUT / "core_score.parquet")
prod["model_id"] = prod["model_id"].str.lower()
prod = prod[prod["ensg_id"].isin(ensg_a)]

curated = pd.read_parquet(VAL / "curated_validation_pairs.parquet")
curated["model_id"] = curated["model_id"].str.lower()
curated["ensg_id"] = curated["ensg_id"].str.lower()

print(f"\n{'Gene':<8}{'N sens':<8}{'Prod hit@20':<14}{'Baseline hit@20':<18}{'miRNA-cand hit@20':<20}{'Delta':<8}")
results_a = {}
for gene, ensg in zip(genes_a, ensg_a):
    cur_gene = curated[curated.ensg_id == ensg]
    sens = set(cur_gene[cur_gene.ground_truth < 0]["model_id"])
    h_prod = hit_at_k(prod, ensg, sens, score_col="core_score")
    h_base = hit_at_k(core_baseline_long, ensg, sens, score_col="core_score")
    h_cand = hit_at_k(core_candidate_long, ensg, sens, score_col="core_score")
    delta = h_cand - h_base
    results_a[gene] = dict(n_sensitive=len(sens), hit_prod=h_prod, hit_baseline=h_base,
                            hit_candidate=h_cand, delta=delta)
    print(f"{gene:<8}{len(sens):<8}{h_prod:<14.4f}{h_base:<18.4f}{h_cand:<20.4f}{delta:+.4f}")

# hit@20 has a very low ceiling here (n_sensitive >> k=20 for both genes),
# so it barely discriminates. Continuous Spearman rho against raw GDSC
# ground_truth (NOT sign-flipped -- same convention as validate_chronos.py's
# own bootstrap_rho check, so deltas are directly comparable to the project's
# established baseline numbers) is the more sensitive comparison.
print(f"\n{'Gene':<8}{'N pairs':<10}{'Baseline rho':<16}{'miRNA-cand rho':<18}{'Delta rho'}")
for gene, ensg in zip(genes_a, ensg_a):
    cur_gene = curated[curated.ensg_id == ensg][["model_id", "ground_truth"]]
    m_base = core_baseline_long[core_baseline_long.ensg_id == ensg].merge(cur_gene, on="model_id", how="inner").dropna()
    m_cand = core_candidate_long[core_candidate_long.ensg_id == ensg].merge(cur_gene, on="model_id", how="inner").dropna()
    rho_base, p_base = spearmanr(m_base["core_score"], m_base["ground_truth"])
    rho_cand, p_cand = spearmanr(m_cand["core_score"], m_cand["ground_truth"])
    results_a[gene]["rho_baseline"] = rho_base
    results_a[gene]["rho_candidate"] = rho_cand
    results_a[gene]["rho_delta"] = rho_cand - rho_base
    results_a[gene]["rho_p_baseline"] = p_base
    results_a[gene]["rho_p_candidate"] = p_cand
    print(f"{gene:<8}{len(m_base):<10}{rho_base:<16.4f}{rho_cand:<18.4f}{rho_cand - rho_base:+.4f}")

print("\n" + "=" * 70)
print("PART B -- metabolomics-as-enzyme-activity vs Chronos essentiality")
print("=" * 70)

metab = pd.read_parquet(CLEANED / "metabolomics.parquet")
<<<<<<< HEAD
metab["model_id"] = metab["model_id"].str.lower()
=======
metab["model_id"] = metab["model_id"].str.upper()
>>>>>>> 20fd4e72bc1fd7b405f50772ccc06f338c5a0bb4
metab = metab.set_index("model_id")

metab_proxy = pd.DataFrame(index=metab.index)
for gene, pair in METABOLITE_MAP.items():
    sub_pct = metab[pair["substrate"]].rank(pct=True)
    prod_pct = metab[pair["product"]].rank(pct=True)
    metab_proxy[GENE_SYMBOL_TO_ENSG[gene]] = prod_pct - sub_pct

metab_proxy_pct = (metab_proxy.rank(pct=True) + 1) / 2  # re-centre to [0,1]-ish for readability only

ch = pd.read_parquet(VAL / "chronos_long.parquet")
ib = pd.read_parquet(VAL / "id_bridge.parquet")
ch = ch.merge(ib, on="sanger_model_id", how="inner")
ch["model_id"] = ch["model_id"].str.lower()
ch["ensg_id"] = ch["ensg_id"].str.lower()
ch["chronos_pct"] = 1.0 - ch.groupby("ensg_id")["essentiality"].rank(pct=True, method="average")

print(f"\n{'Gene':<8}{'N pairs':<10}{'rho':<10}{'p-value':<12}{'classification'}")
results_b = {}
for gene, pair in METABOLITE_MAP.items():
    ensg = GENE_SYMBOL_TO_ENSG[gene]
    proxy = metab_proxy[ensg].reset_index().rename(columns={ensg: "proxy"})
    ch_gene = ch[ch.ensg_id == ensg][["model_id", "chronos_pct"]]
    merged = proxy.merge(ch_gene, on="model_id", how="inner").dropna()
    if len(merged) < 30:
        print(f"{gene:<8}{len(merged):<10} too few paired cell lines, skipped")
        results_b[gene] = dict(n=len(merged), rho=np.nan, p=np.nan, classification="insufficient_n")
        continue
    rho, p = spearman_with_p(merged["proxy"], merged["chronos_pct"])
    cls = classify(rho, p)
    results_b[gene] = dict(n=len(merged), rho=rho, p=p, classification=cls)
    print(f"{gene:<8}{len(merged):<10}{rho:<10.4f}{p:<12.4g}{cls}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("\nA) miRNA inhibitor on expression (positive delta = candidate beats baseline):")
for gene, r in results_a.items():
    verdict = "IMPROVES" if r["rho_delta"] > 0.01 else ("WORSENS" if r["rho_delta"] < -0.01 else "NO CHANGE")
    print(f"  {gene}: {verdict}  (rho {r['rho_baseline']:.4f} -> {r['rho_candidate']:.4f}, "
          f"delta {r['rho_delta']:+.4f} | hit@20 {r['hit_baseline']:.4f} -> {r['hit_candidate']:.4f}, "
          f"n_sensitive={r['n_sensitive']})")

print("\nB) metabolomics enzyme-activity proxy vs real CRISPR essentiality:")
for gene, r in results_b.items():
    print(f"  {gene}: {r['classification']}  (rho={r['rho']}, n={r['n']})" if not np.isnan(r.get("rho", np.nan))
          else f"  {gene}: {r['classification']} (n={r['n']})")

import json as _json
out_path = OUT / "test_run_metabolomics_mirna_results.json"
with open(out_path, "w") as f:
    _json.dump({"mirna_inhibitor_test": results_a, "metabolomics_test": results_b}, f, indent=2, default=float)
print(f"\nWritten: {out_path}")
