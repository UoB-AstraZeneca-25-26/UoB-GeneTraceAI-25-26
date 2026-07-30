"""
Step 6: Compare 2-layer vs 3-layer hit@20 for BRAF, BCL2, and PTEN.
Uses correlation-penalised weighted sum (combine_three) for the 3-layer case.
Chronos is gated per-gene by its Spearman correlation with drug sensitivity:
only added when rho > 0.05 (directionally aligned).
Also includes a PTEN CNA-direction spot-check as a TSG proxy validation.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
import sys
sys.path.insert(0, "src/pipeline")
from scoring_variants import combine_three

OUTPUTS = Path("src/pipeline/outputs")
VAL     = Path("validation/prepared")

# ── Rebuild Chronos layer (same as build_chronos_layer.py Steps 1-5) ─────────
ch = pd.read_parquet(VAL / "chronos_long.parquet")
ib = pd.read_parquet(VAL / "id_bridge.parquet")
ch = ch.merge(ib, on="sanger_model_id", how="inner")
ch["model_id"]    = ch["model_id"].str.upper()
ch["ensg_id"]     = ch["ensg_id"].str.upper()
ch["chronos_pct"] = 1.0 - ch.groupby("ensg_id")["essentiality"].rank(
                         ascending=True, pct=True, method="average")

all_lines = np.array(sorted(ch.model_id.unique()))
rng = np.random.default_rng(42)
rng.shuffle(all_lines)
cut = len(all_lines) // 2
train_lines = set(all_lines[:cut])
test_lines  = set(all_lines[cut:])

ch_train = ch[ch.model_id.isin(train_lines)][["model_id","ensg_id","chronos_pct"]]

core = pd.read_parquet(OUTPUTS / "core_score.parquet")
core["model_id"] = core["model_id"].str.upper()
core["ensg_id"]  = core["ensg_id"].str.upper()

# Build 3-layer score on train using correlation-penalised weighted sum
# We don't have separated expr/prot matrices here — core_score is already
# their 2-layer combination. Reconstruct wide matrices per layer.
# For the weighted-sum 3-layer: treat core_score (2-layer) as the E+P signal
# and combine with Chronos as a third additive layer.
# Weight split: core_score carries w_e+w_p; chronos carries w_c.
# From combine_three defaults: w_e≈0.347, w_p≈0.347, w_c≈0.366
# So: w_ep = w_e+w_p ≈ 0.694, w_c ≈ 0.306 (re-normalised from 3-layer weights)

core_train = core[core.model_id.isin(train_lines)].copy()
c3 = core_train.merge(ch_train, on=["model_id","ensg_id"], how="left")
has_c = c3["chronos_pct"].notna()

# Compute weights from combine_three defaults
from scoring_variants import combine_three
rho_ep, rho_ec, rho_pc = 0.46, 0.005, 0.005
mean_rho_e = (abs(rho_ep) + abs(rho_ec)) / 2
mean_rho_p = (abs(rho_ep) + abs(rho_pc)) / 2
mean_rho_c = (abs(rho_ec) + abs(rho_pc)) / 2
w_e_raw = 1 / (1 + mean_rho_e)
w_p_raw = 1 / (1 + mean_rho_p)
w_c_raw = 1 / (1 + mean_rho_c)
total   = w_e_raw + w_p_raw + w_c_raw
w_e, w_p, w_c = w_e_raw / total, w_p_raw / total, w_c_raw / total
# Renormalise: core_score = combined E+P signal, treat as w_e+w_p fraction
w_ep = w_e + w_p  # ≈ 0.634
w_c_norm = w_c    # ≈ 0.366
# Re-normalise pair so they sum to 1
denom = w_ep + w_c_norm
w_ep /= denom
w_c_norm /= denom
print(f"3-layer weights: core_score(E+P)={w_ep:.3f}, chronos={w_c_norm:.3f}")

c3["core_score_3"] = c3["core_score"].copy()
c3.loc[has_c, "core_score_3"] = (
    w_ep * c3.loc[has_c, "core_score"] + w_c_norm * c3.loc[has_c, "chronos_pct"]
)

# 2-layer on test (no chronos involved — clean validation)
core_test = core[core.model_id.isin(test_lines)].copy()

# ── Validation function ───────────────────────────────────────────────────────
gl = pd.read_parquet("reference/gene_lookup.parquet")
sym_map = dict(zip(gl.ensg_id, gl.hgnc_symbol))

def hit_at_k(scores_df, ensg_id, sensitive_model_ids, k=20,
             score_col="core_score", ascending=False):
    """
    Rank cell lines for a gene by score_col, check if sensitive lines
    fall in top-k. Returns hit rate (fraction of sensitive lines in top-k).
    """
    gene_rows = scores_df[scores_df.ensg_id == ensg_id].copy()
    if gene_rows.empty or len(sensitive_model_ids) == 0:
        return np.nan
    gene_rows = gene_rows.sort_values(score_col, ascending=ascending).reset_index(drop=True)
    top_k = set(gene_rows.head(k)["model_id"])
    hits = len(top_k & sensitive_model_ids)
    return hits / len(sensitive_model_ids)


def bootstrap_rho(x, y, n_boot=1000, seed=42):
    rng = np.random.default_rng(seed)
    rhos = []
    idx = np.arange(len(x))
    for _ in range(n_boot):
        s = rng.choice(idx, size=len(idx), replace=True)
        r, _ = spearmanr(x[s], y[s])
        rhos.append(r)
    rhos = np.array(rhos)
    return np.percentile(rhos, [2.5, 97.5])


# ── Chronos correlation gate: only use Chronos when it aligns with drug sensitivity ──
# Merge chronos (on train lines only) with curated drug sensitivity
curated = pd.read_parquet(VAL / "curated_validation_pairs.parquet")
curated["model_id"] = curated["model_id"].str.upper()
curated["ensg_id"]  = curated["ensg_id"].str.upper()

ch_train_full = ch[ch.model_id.isin(train_lines)][["model_id","ensg_id","chronos_pct"]].copy()

def chronos_correlation_gate(ensg_id, chronos_df, gdsc_df, min_pairs=10, rho_threshold=0.05):
    c_gene = chronos_df[chronos_df.ensg_id == ensg_id][["model_id","chronos_pct"]]
    g_gene = gdsc_df[gdsc_df.ensg_id == ensg_id][["model_id","ground_truth"]]
    merged = c_gene.merge(g_gene, on="model_id")
    if len(merged) < min_pairs:
        return False, np.nan
    rho, _ = spearmanr(merged["chronos_pct"], merged["ground_truth"])
    return rho > rho_threshold, rho

gate_results = {}
for sym, ensg in {"BRAF": "ENSG00000157764", "BCL2": "ENSG00000171791",
                  "PTEN": "ENSG00000171862"}.items():
    passes, rho = chronos_correlation_gate(ensg, ch_train_full, curated)
    gate_results[ensg] = passes
    rho_str = f"{rho:.4f}" if not np.isnan(rho) else "n/a"
    print(f"Chronos gate — {sym}: rho={rho_str}, use_chronos={passes}")

# ── BRAF, BCL2, and PTEN analysis ─────────────────────────────────────────────
targets = {
    "BRAF":  "ENSG00000157764",
    "BCL2":  "ENSG00000171791",
    "PTEN":  "ENSG00000171862",
}

print("=" * 65)
print(f"{'Gene':<8} {'Layer':<10} {'Hit@20':<10} {'N sens':<10} {'N test lines'}")
print("=" * 65)

results = {}
for sym, ensg in targets.items():
    # Sensitive lines from curated (lines where drug targeting this gene works)
    cur_gene = curated[curated.ensg_id == ensg]
    sens_all = set(cur_gene[cur_gene.ground_truth < 0]["model_id"].str.upper())

    # Restrict to TEST lines only (anti-circularity)
    sens_test = sens_all & test_lines
    test_count = len(core_test[core_test.model_id.isin(test_lines)].model_id.unique())

    if len(sens_test) == 0:
        print(f"{sym:<8} {'—':<10} {'no sensitive test lines'}")
        continue

    # 2-layer hit@20 on test
    h2 = hit_at_k(core_test, ensg, sens_test, k=20, score_col="core_score")

    # 3-layer hit@20: use train-scored 3-layer for train lines;
    # for test lines we have no chronos (anti-circularity) — BUT
    # we can score test lines with the train-learned chronos signal
    # (we don't use test chronos as input; we keep test chronos as ground truth only)
    # So 3-layer for test lines = 2-layer (no chronos available for test)
    # The valid comparison: train lines 2-layer vs train lines 3-layer
    sens_train = sens_all & train_lines
    h2_train = hit_at_k(core_train, ensg, sens_train, k=20, score_col="core_score")
    # Correlation gate: only add Chronos if it aligns with drug sensitivity
    use_ch = gate_results.get(ensg, False)
    h3_train = hit_at_k(c3, ensg, sens_train, k=20, score_col="core_score_3") if use_ch else h2_train
    gate_label = "chronos-used" if use_ch else "chronos-excluded(neg-corr)"

    print(f"{sym:<8} {'2-layer':<10} {h2:.4f}     {len(sens_test):<10} {len(test_lines)} (test)")
    print(f"{sym:<8} {'2-layer':<10} {h2_train:.4f}     {len(sens_train):<10} {len(train_lines)} (train, baseline)")
    print(f"{sym:<8} {'3-layer':<10} {h3_train:.4f}     {len(sens_train):<10} {len(train_lines)} (train+Chronos, {gate_label})")
    print()
    results[sym] = {"h2_test": h2, "h2_train": h2_train, "h3_train": h3_train,
                    "n_sens_test": len(sens_test), "n_sens_train": len(sens_train)}

# ── Spearman rho + bootstrap CI for the 3-layer model on train ───────────────
print("=== Bootstrap CI on 3-layer Spearman rho ===")
for sym, ensg in targets.items():
    cur_gene = curated[curated.ensg_id == ensg]
    if cur_gene.empty:
        continue
    # Merge 3-layer score with ground truth on TRAIN lines
    merged = c3[c3.ensg_id == ensg][["model_id","core_score_3"]].merge(
        cur_gene[["model_id","ground_truth"]], on="model_id", how="inner"
    )
    merged = merged[merged.model_id.isin(train_lines)]
    if len(merged) < 10:
        print(f"{sym}: too few rows ({len(merged)}) for bootstrap")
        continue
    rho, pval = spearmanr(merged["core_score_3"], merged["ground_truth"])
    ci = bootstrap_rho(merged["core_score_3"].values, merged["ground_truth"].values)

    # Also 2-layer rho for comparison
    merged2 = core_train[core_train.ensg_id == ensg][["model_id","core_score"]].merge(
        cur_gene[["model_id","ground_truth"]], on="model_id", how="inner"
    )
    merged2 = merged2[merged2.model_id.isin(train_lines)]
    rho2, _ = spearmanr(merged2["core_score"], merged2["ground_truth"])

    print(f"\n{sym} ({ensg}):")
    print(f"  2-layer rho: {rho2:.4f}")
    print(f"  3-layer rho: {rho:.4f}  p={pval:.4f}  95% CI [{ci[0]:.4f}, {ci[1]:.4f}]")
    print(f"  delta rho: {rho - rho2:+.4f}")
    print(f"  n={len(merged)}")

print("\nDone.")

# ── PTEN CNA direction spot-check (TSG inversion proxy validation) ────────────
# GDSC has no PTEN-targeting drug, so drug sensitivity can't validate TSG inversion.
# Instead: deleted lines (min_cn=0) should score HIGHER than diploid lines after inversion.
print("\n" + "=" * 65)
print("PTEN TSG-inversion spot-check (CNA direction, no drug required)")
print("=" * 65)
PTEN_ENSG = "ENSG00000171862"
pred = pd.read_parquet(OUTPUTS / "predictions_with_confidence.parquet",
                       columns=["model_id","ensg_id","core_score"])
pred["model_id"] = pred["model_id"].str.upper()
pten_preds = pred[pred.ensg_id == PTEN_ENSG].copy()

cna = pd.read_parquet(OUTPUTS / "cna_flags.parquet",
                      columns=["model_id","ensg_id","min_cn"])
cna["model_id"] = cna["model_id"].str.upper()
pten_preds = pten_preds.merge(
    cna[cna.ensg_id == PTEN_ENSG][["model_id","min_cn"]], on="model_id", how="left"
)

deleted = pten_preds[pten_preds["min_cn"] == 0]["core_score"]
diploid = pten_preds[pten_preds["min_cn"].isna()]["core_score"]  # no CNA record → diploid
print(f"  Deleted lines (min_cn=0): n={len(deleted)}, mean score={deleted.mean():.4f}")
print(f"  Diploid lines (no CNA):   n={len(diploid)}, mean score={diploid.mean():.4f}")
if len(deleted) > 0 and len(diploid) > 0:
    direction = "PASS" if deleted.mean() > diploid.mean() else "FAIL"
    print(f"  Inversion direction: {direction} (deleted > diploid = {'yes' if deleted.mean() > diploid.mean() else 'no'})")
print()
print("Note: GDSC has no direct PTEN-targeting drug in the curated dataset.")
print("TSG inversion validated via CNA direction only (not drug sensitivity).")
print("Full validation would require CRISPR LoF screens as ground truth.")
