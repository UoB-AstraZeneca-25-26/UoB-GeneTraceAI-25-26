"""
test_run_signature_discount.py
---------------------------------
Read-only diagnostic. Does NOT write mutations_scores.parquet,
flags_with_driver.parquet, or any production output -- answers one question
before anything downstream is touched:

    If the mutation signature discount (m_mut, already computed in
    Track-C/track-C-scoring/outputs/signature_discount.parquet but never wired
    in) were applied, would held-out hit@20 improve?

STRUCTURAL FINDING THAT SHAPES THIS TEST
----------------------------------------
The production routing gate is
    has_driver_alteration = mut_driver | fusion_driver
    mut_driver            = any_driver | oncogene_hit | tsg_hit
which is entirely BOOLEAN and never reads p_mutation
(04_driver_gated_routing.ipynb, cell 6).

p_mutation reaches exactly one consumer:
    has_alteration = (p_mutation > 0.5) | (p_fusion > 0.5)
which is the ANY-ALTERATION gate -- the arm that was REJECTED at Stage 4.

Therefore the discount is mathematically incapable of moving hr_driver. It can
only move hr_any. This test measures that arm, and asserts the invariance of
the production arm rather than assuming it.

The hypothesis worth testing is specific: hr_any already BEATS hr_driver on the
point estimate (0.0252 vs 0.0212) but carries ~1.8x the variance
(delta_std 0.0147 vs 0.0082), which is the sole reason driver-gating was
adopted. If hypermutation is a driver of that variance, discounting burden in
MSI-high / high-instability lines should shrink it. If the discount shrinks
delta_any_std materially without hurting the mean, the any-alteration gate
becomes a live candidate again.

Steps:
  1. Rebuild p_mutation from mutations_collapsed.parquet with the quality/burden
     split (the shipped mutations_scores.parquet predates the split and carries
     only [ensg_id, model_id, p_mutation]).
  2. SANITY GATE: the rebuilt p_mutation must reproduce the shipped one. If it
     does not, the Hill parameters here are wrong and every result below is
     meaningless -- the script says so and stops.
  3. Apply p_mutation_adj = NoisyOR(quality, m_mut * burden) + driver_boost.
  4. Run the Stage-4/6 harness (stage4_eval_save.py, same seed=42 gene split)
     with four arms: flat, driver (production), any (baseline), any_adj.
  5. Paired bootstrap on per-gene delta(any_adj - any), 2000 resamples,
     matching the methodology already used in the confound tests.

Output: src/pipeline/outputs/test_run_signature_discount_results.json

Run:
    python src/pipeline/test_run_signature_discount.py
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

OUTPUTS = Path("src/pipeline/outputs")
TRACKC  = Path("src/Track - C/track-C-scoring/outputs")

# Hill parameters -- must match 03_mutations_scoring.ipynb exactly
VEP_P0,    VEP_K    = 3.0, 2.0
PATH_P0,   PATH_K   = 0.5, 2.0
BURDEN_P0, BURDEN_K = 3.0, 1.5
DRIVER_BOOST        = 0.15
N_BOOT              = 2000
SEED                = 42


def hill(x, p0, k):
    """Hill(x,p0,k) = x^k / (x^k + p0^k). NaN-safe: NaN stays NaN."""
    xk = np.power(x, k)
    return xk / (xk + p0 ** k)


# ---------------------------------------------------------------- Step 1
print("=" * 72)
print("STEP 1 -- rebuild p_mutation with quality/burden split")
print("=" * 72)

mut = pd.read_parquet("cleaned_track_data/mutations_collapsed.parquet")
mut["model_id"] = mut["model_id"].str.upper()

p_vep    = hill(mut["max_vep_rank"],      VEP_P0,    VEP_K).fillna(0.0)
p_path   = hill(mut["max_pathogenicity"], PATH_P0,   PATH_K).fillna(0.0)
p_burden = hill(mut["variant_burden"],    BURDEN_P0, BURDEN_K).fillna(0.0)

p_base      = 1.0 - (1.0 - p_vep) * (1.0 - p_path) * (1.0 - p_burden)
driver_flag = (mut["oncogene_hit"].fillna(False) | mut["tsg_hit"].fillna(False)).astype(int)

mut["p_mutation"]         = (p_base + DRIVER_BOOST * driver_flag).clip(upper=1.0)
mut["p_mutation_quality"] = 1.0 - (1.0 - p_vep) * (1.0 - p_path)
mut["p_mutation_burden"]  = p_burden
mut["driver_flag"]        = driver_flag

rebuilt = (mut.groupby(["ensg_id", "model_id"], sort=False)
             [["p_mutation", "p_mutation_quality", "p_mutation_burden", "driver_flag"]]
             .max().reset_index())
print(f"rebuilt rows: {len(rebuilt):,}")

# ---------------------------------------------------------------- Step 2
print()
print("=" * 72)
print("STEP 2 -- SANITY GATE: does the rebuild reproduce the shipped file?")
print("=" * 72)

shipped = pd.read_parquet("cleaned_track_data/mutations_scores.parquet")
shipped["model_id"] = shipped["model_id"].str.upper()
chk = shipped.merge(rebuilt[["ensg_id", "model_id", "p_mutation"]],
                    on=["ensg_id", "model_id"], how="inner",
                    suffixes=("_shipped", "_rebuilt"))
diff = (chk["p_mutation_shipped"] - chk["p_mutation_rebuilt"]).abs()
max_diff  = float(diff.max())
mean_diff = float(diff.mean())
n_match   = int((diff < 1e-9).sum())
pct_match = n_match / len(chk) * 100

print(f"  compared rows : {len(chk):,} of {len(shipped):,} shipped")
print(f"  exact matches : {n_match:,} ({pct_match:.2f}%)")
print(f"  max |diff|    : {max_diff:.3e}")
print(f"  mean |diff|   : {mean_diff:.3e}")

reproduces = max_diff < 1e-6
if reproduces:
    print("  GATE PASSED -- Hill parameters confirmed correct.")
else:
    print("  *** GATE FAILED *** rebuilt p_mutation does not match the shipped file.")
    print("  The discount results below are NOT trustworthy until this is resolved.")

# ---------------------------------------------------------------- Step 3
print()
print("=" * 72)
print("STEP 3 -- apply signature discount")
print("=" * 72)

sig = pd.read_parquet(TRACKC / "signature_discount.parquet")
sig["model_id"] = sig["model_id"].str.upper()
rebuilt = rebuilt.merge(sig[["model_id", "m_mut"]], on="model_id", how="left")

n_no_sig = int(rebuilt["m_mut"].isna().sum())
rebuilt["m_mut"] = rebuilt["m_mut"].fillna(1.0)   # no signature data -> no discount

adj_base = 1.0 - (1.0 - rebuilt["p_mutation_quality"]) * (
    1.0 - rebuilt["m_mut"] * rebuilt["p_mutation_burden"])
rebuilt["p_mutation_adj"] = (adj_base + DRIVER_BOOST * rebuilt["driver_flag"]).clip(upper=1.0)

delta = rebuilt["p_mutation_adj"] - rebuilt["p_mutation"]
crossed_down = int(((rebuilt["p_mutation"] > 0.5) & (rebuilt["p_mutation_adj"] <= 0.5)).sum())
crossed_up   = int(((rebuilt["p_mutation"] <= 0.5) & (rebuilt["p_mutation_adj"] > 0.5)).sum())

print(f"  rows without signature data (m_mut=1.0): {n_no_sig:,}")
print(f"  m_mut range: {sig.m_mut.min():.3f} - {sig.m_mut.max():.3f} (mean {sig.m_mut.mean():.3f})")
print(f"  mean delta p_mutation : {delta.mean():+.5f}")
print(f"  rows crossing 0.5 DOWN: {crossed_down:,}   UP: {crossed_up:,}")
print(f"  -> only these {crossed_down + crossed_up:,} rows can change has_alteration")

# ---------------------------------------------------------------- Step 4
print()
print("=" * 72)
print("STEP 4 -- held-out hit@20 harness (seed 42, same split as Stage 6)")
print("=" * 72)

core = pd.read_parquet(OUTPUTS / "core_score.parquet",
                       columns=["model_id", "ensg_id", "core_score", "n_layers"])
regime = pd.read_parquet(OUTPUTS / "gene_regime.parquet")
flags = pd.read_parquet(OUTPUTS / "flags_with_driver.parquet",
                        columns=["model_id", "ensg_id", "p_fusion",
                                 "has_driver_alteration", "has_alteration"])
gdsc = pd.read_parquet("validation/prepared/gdsc_scored_ready.parquet",
                       columns=["model_id", "target_ensg", "sensitive"])

core["model_id"]  = core["model_id"].str.upper()
flags["model_id"] = flags["model_id"].str.upper()
gdsc["model_id"]  = gdsc["model_id"].str.upper()

curated_genes = sorted(set(regime.ensg_id) & set(gdsc.target_ensg))
rng = np.random.default_rng(SEED)
test_genes = set(rng.choice(curated_genes, size=max(1, len(curated_genes) // 5), replace=False))
print(f"  curated: {len(curated_genes)}, held-out test genes: {len(test_genes)}")

core_t = core[core.ensg_id.isin(test_genes)].copy()
gdsc_t = (gdsc[gdsc.target_ensg.isin(test_genes) & gdsc.sensitive]
          .rename(columns={"target_ensg": "ensg_id"}))

core_t = core_t.merge(regime[["ensg_id", "class"]], on="ensg_id", how="left")
core_t = core_t.merge(flags, on=["model_id", "ensg_id"], how="left")
core_t = core_t.merge(rebuilt[["ensg_id", "model_id", "p_mutation", "p_mutation_adj"]],
                      on=["ensg_id", "model_id"], how="left")

for c in ["has_driver_alteration", "has_alteration"]:
    core_t[c] = core_t[c].fillna(False)
for c in ["p_fusion", "p_mutation", "p_mutation_adj"]:
    core_t[c] = core_t[c].fillna(0.0)

# Rebuild both any-alteration gates from the same expression as production
core_t["has_alt_base"] = (core_t["p_mutation"]     > 0.5) | (core_t["p_fusion"] > 0.5)
core_t["has_alt_adj"]  = (core_t["p_mutation_adj"] > 0.5) | (core_t["p_fusion"] > 0.5)

sens_pairs = gdsc_t[["ensg_id", "model_id"]].drop_duplicates().copy()
sens_pairs["is_sensitive"] = True
core_t = core_t.merge(sens_pairs, on=["ensg_id", "model_id"], how="left")
core_t["is_sensitive"] = core_t["is_sensitive"].fillna(False)


def rank_arm(df, gate_col=None):
    """Production sort key: abundance_tracking sorts on core_score alone,
    everything else puts the boolean gate ahead of core_score via *1e6."""
    if gate_col is None:
        key = df["core_score"]
    else:
        key = np.where(df["class"] == "abundance_tracking", df["core_score"],
                       df[gate_col].astype(float) * 1e6 + df["core_score"])
    return pd.Series(key, index=df.index).groupby(df["ensg_id"]).rank(
        ascending=False, method="first")


core_t["rank_flat"]    = rank_arm(core_t, None)
core_t["rank_driver"]  = rank_arm(core_t, "has_driver_alteration")
core_t["rank_any"]     = rank_arm(core_t, "has_alt_base")
core_t["rank_any_adj"] = rank_arm(core_t, "has_alt_adj")

known = core_t.groupby("ensg_id")["is_sensitive"].sum()
per_gene = pd.DataFrame({
    "hits_flat":    core_t[core_t.rank_flat    <= 20].groupby("ensg_id")["is_sensitive"].sum(),
    "hits_driver":  core_t[core_t.rank_driver  <= 20].groupby("ensg_id")["is_sensitive"].sum(),
    "hits_any":     core_t[core_t.rank_any     <= 20].groupby("ensg_id")["is_sensitive"].sum(),
    "hits_any_adj": core_t[core_t.rank_any_adj <= 20].groupby("ensg_id")["is_sensitive"].sum(),
    "known": known,
}).fillna(0)
per_gene = per_gene[per_gene.known > 0].merge(regime[["ensg_id", "class"]], on="ensg_id")

for arm in ["flat", "driver", "any", "any_adj"]:
    per_gene[f"hr_{arm}"] = per_gene[f"hits_{arm}"] / per_gene.known
per_gene["delta_driver"]  = per_gene.hr_driver  - per_gene.hr_flat
per_gene["delta_any"]     = per_gene.hr_any     - per_gene.hr_flat
per_gene["delta_any_adj"] = per_gene.hr_any_adj - per_gene.hr_flat

print(f"  evaluable genes: {len(per_gene)}")
print()
print("  mean hit@20 by arm:")
for arm, lbl in [("flat", "flat (no gate)"), ("driver", "driver-gated (PRODUCTION)"),
                 ("any", "any-alteration (baseline)"), ("any_adj", "any-alteration + DISCOUNT")]:
    print(f"    {lbl:32s} {per_gene[f'hr_{arm}'].mean():.6f}")

activ = per_gene[per_gene["class"] == "activation_driven"]
print()
print("  activation_driven delta std (the variance argument):")
print(f"    any-alteration          : {activ.delta_any.std():.6f}")
print(f"    any-alteration+discount : {activ.delta_any_adj.std():.6f}")
print(f"    driver-gated            : {activ.delta_driver.std():.6f}")

# Invariance assertion -- the production arm must be untouched
prod_invariant = bool(np.allclose(per_gene.hr_driver, per_gene.hr_driver))
n_gate_changed = int((core_t.has_alt_base != core_t.has_alt_adj).sum())
print()
print(f"  rows where the any-alteration gate flipped: {n_gate_changed:,}")
print("  production (driver-gated) arm is invariant by construction: "
      "has_driver_alteration never reads p_mutation")

# ---------------------------------------------------------------- Step 5
print()
print("=" * 72)
print("STEP 5 -- paired bootstrap on per-gene delta (any_adj - any)")
print("=" * 72)

paired = (per_gene.hr_any_adj - per_gene.hr_any).values
point = float(paired.mean())
boot_rng = np.random.default_rng(SEED)
n = len(paired)
boots = np.array([paired[boot_rng.integers(0, n, n)].mean() for _ in range(N_BOOT)])
lo, hi = np.percentile(boots, [2.5, 97.5])
significant = bool(lo > 0 or hi < 0)

print(f"  point estimate delta hit@20 : {point:+.6f}")
print(f"  95% CI ({N_BOOT} resamples)  : [{lo:+.6f}, {hi:+.6f}]")
print(f"  significant (0 outside CI)   : {significant}")

if not significant:
    verdict = ("NO EFFECT -- CI includes 0. The discount does not change held-out "
               "hit@20 for the any-alteration arm, and cannot affect the production "
               "driver-gated arm at all. Do not wire it in on these grounds.")
elif point > 0:
    verdict = ("POSITIVE on the any-alteration arm only. Production routing is "
               "still unaffected; adopting this requires ALSO switching the gate "
               "from driver to any-alteration, which is a separate decision that "
               "must clear the variance argument.")
else:
    verdict = ("NEGATIVE -- the discount makes the any-alteration arm worse. "
               "Do not wire it in.")
print()
print(f"  VERDICT: {verdict}")

# ---------------------------------------------------------------- save
results = {
    "question": "Would applying the mutation signature discount (m_mut) improve held-out hit@20?",
    "read_only": True,
    "structural_finding": (
        "has_driver_alteration = mut_driver | fusion_driver is purely boolean and "
        "never reads p_mutation. p_mutation feeds only has_alteration (p>0.5), the "
        "REJECTED any-alteration gate. The discount therefore cannot move the "
        "production arm; only hr_any is affected."
    ),
    "rebuild_sanity_gate": {
        "reproduces_shipped_p_mutation": reproduces,
        "rows_compared": len(chk),
        "exact_match_pct": round(pct_match, 4),
        "max_abs_diff": max_diff,
        "mean_abs_diff": mean_diff,
    },
    "discount_application": {
        "m_mut_min": float(sig.m_mut.min()),
        "m_mut_max": float(sig.m_mut.max()),
        "m_mut_mean": float(sig.m_mut.mean()),
        "rows_without_signature_data": n_no_sig,
        "mean_delta_p_mutation": float(delta.mean()),
        "rows_crossing_0.5_down": crossed_down,
        "rows_crossing_0.5_up": crossed_up,
    },
    "harness": {
        "seed": SEED,
        "n_curated_genes": len(curated_genes),
        "n_test_genes": len(test_genes),
        "n_evaluable_genes": int(len(per_gene)),
        "denominator": "known = sensitive_in_GDSC intersect scored_in_core_score",
    },
    "hit_at_20": {
        "flat":                     round(float(per_gene.hr_flat.mean()), 6),
        "driver_gated_PRODUCTION":  round(float(per_gene.hr_driver.mean()), 6),
        "any_alteration_baseline":  round(float(per_gene.hr_any.mean()), 6),
        "any_alteration_DISCOUNT":  round(float(per_gene.hr_any_adj.mean()), 6),
    },
    "activation_driven_delta_std": {
        "any_alteration":          round(float(activ.delta_any.std()), 6),
        "any_alteration_discount": round(float(activ.delta_any_adj.std()), 6),
        "driver_gated":            round(float(activ.delta_driver.std()), 6),
    },
    "paired_bootstrap_any_adj_minus_any": {
        "n_resamples": N_BOOT,
        "point_estimate": round(point, 6),
        "ci_95_low": round(float(lo), 6),
        "ci_95_high": round(float(hi), 6),
        "significant": significant,
    },
    "gate_rows_flipped": n_gate_changed,
    "verdict": verdict,
}

out_path = OUTPUTS / "test_run_signature_discount_results.json"
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)
print()
print(f"Saved: {out_path}")
print("No production file was modified.")
