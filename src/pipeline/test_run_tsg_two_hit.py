"""
test_run_tsg_two_hit.py
--------------------------
Read-only diagnostic. Writes no production file.

QUESTION
--------
Knudson's two-hit hypothesis: a tumour suppressor generally needs BOTH alleles
inactivated to lose function. A TSG with a damaging mutation AND loss of the
remaining wild-type allele is genuinely inactivated; either alone is much weaker
evidence. Does gating on two-hit status rank cell lines better than the current
production gate for TSGs?

TWO-HIT DEFINITION (per model_id, ensg_id -- all per-locus, no genome-wide scalars)
    hit 1 : tsg_hit                    -- a damaging mutation in a TSG
    hit 2 : max_vaf >= 0.75            -- high variant allele fraction. If both
                                          alleles were present with one mutated,
                                          VAF ~ 0.5; VAF near 1.0 implies the
                                          wild-type copy is gone (LOH proxy)
         OR multi_hit_high_impact      -- two independent damaging mutations
         OR min_cn < 1.5               -- deletion at this locus (CNA second hit)

    two_hit = hit1 AND (any of hit2)

NOTE ON lohfraction: signatures_model_level.parquet carries `lohfraction`, but it
is a GENOME-WIDE per-cell-line scalar, not per-locus. Using it here would repeat
the gene-agnostic-scalar mistake that got metabolomics and miRNA rejected as
lineage confounds (see DECISIONS.md). It is deliberately NOT used. max_vaf is the
correct per-locus LOH proxy and is 100% populated on TSG rows.

*** POWER LIMITATION -- READ BEFORE INTERPRETING ***
Only 14 TSG/both genes have any GDSC sensitive pairs, and the curated 141-gene
regime table contains ZERO loss_of_function genes (95 activation_driven,
46 abundance_tracking). The seed-42 held-out split therefore contains no TSG
genes at all, so this CANNOT be run as a held-out test. It runs on all available
TSG genes and is exploratory, not confirmatory. n=14 is small; treat any result
as hypothesis-generating.

DIRECTION
---------
For a TSG, LOW abundance is the actionable signal (loss of the suppressor).
core_score.parquet is NOT inverted -- inversion happens later in
build_full_predictions.py. So this script tests both directions explicitly,
which also checks whether the TSG inversion rule is itself correct.

Output: src/pipeline/outputs/test_run_tsg_two_hit_results.json

Run:
    python src/pipeline/test_run_tsg_two_hit.py
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

OUTPUTS = Path("src/pipeline/outputs")
VAF_LOH_THRESHOLD = 0.75
DEL_THRESHOLD     = 1.5
TOP_K             = 20
N_BOOT            = 2000
SEED              = 42

# ---------------------------------------------------------------- Step 1
print("=" * 72)
print("STEP 1 -- build the two-hit flag")
print("=" * 72)

gl = pd.read_parquet("reference/gene_lookup.parquet",
                     columns=["ensg_id", "hgnc_symbol", "gene_role"])
gl["ensg_id"] = gl["ensg_id"].str.split(".").str[0]
tsg_genes = set(gl[gl.gene_role.isin(["tsg", "both"])].ensg_id)
sym = gl.drop_duplicates("ensg_id").set_index("ensg_id")["hgnc_symbol"].to_dict()
print(f"  TSG/both genes in lookup: {len(tsg_genes):,}")

mut = pd.read_parquet("cleaned_track_data/mutations_collapsed.parquet")
mut["model_id"] = mut["model_id"].str.upper()
mut = mut[mut.ensg_id.isin(tsg_genes)].copy()

cna = pd.read_parquet(OUTPUTS / "cna_flags.parquet",
                      columns=["model_id", "ensg_id", "min_cn"])
cna["model_id"] = cna["model_id"].str.upper()

mut = mut.merge(cna, on=["model_id", "ensg_id"], how="left")

mut["hit1"]        = mut["tsg_hit"].fillna(False).astype(bool)
mut["h2_vaf"]      = mut["max_vaf"].fillna(0) >= VAF_LOH_THRESHOLD
mut["h2_multihit"] = mut["multi_hit_high_impact"].fillna(False).astype(bool)
mut["h2_del"]      = mut["min_cn"].fillna(99) < DEL_THRESHOLD
mut["hit2"]        = mut["h2_vaf"] | mut["h2_multihit"] | mut["h2_del"]
mut["two_hit"]     = mut["hit1"] & mut["hit2"]

print(f"  TSG (model,gene) rows      : {len(mut):,}")
print(f"    hit1 (tsg_hit)           : {int(mut.hit1.sum()):,}")
print(f"    hit2 via VAF>={VAF_LOH_THRESHOLD}        : {int(mut.h2_vaf.sum()):,}")
print(f"    hit2 via multi-hit       : {int(mut.h2_multihit.sum()):,}")
print(f"    hit2 via deletion        : {int(mut.h2_del.sum()):,}")
print(f"    TWO-HIT                  : {int(mut.two_hit.sum()):,}"
      f"  ({mut.two_hit.mean()*100:.1f}% of TSG rows)")

# ---------------------------------------------------------------- Step 2
print()
print("=" * 72)
print("STEP 2 -- assemble evaluation set")
print("=" * 72)

core = pd.read_parquet(OUTPUTS / "core_score.parquet",
                       columns=["model_id", "ensg_id", "core_score"])
core["model_id"] = core["model_id"].str.upper()
core = core[core.ensg_id.isin(tsg_genes)].copy()

flags = pd.read_parquet(OUTPUTS / "flags_with_driver.parquet",
                        columns=["model_id", "ensg_id", "has_driver_alteration"])
flags["model_id"] = flags["model_id"].str.upper()

gdsc = pd.read_parquet("validation/prepared/gdsc_scored_ready.parquet",
                       columns=["model_id", "target_ensg", "sensitive"])
gdsc["model_id"] = gdsc["model_id"].str.upper()
sens = (gdsc[gdsc.sensitive & gdsc.target_ensg.isin(tsg_genes)]
        .rename(columns={"target_ensg": "ensg_id"})[["ensg_id", "model_id"]]
        .drop_duplicates())
sens["is_sensitive"] = True

eval_genes = sorted(set(sens.ensg_id) & set(core.ensg_id))
print(f"  TSG genes with GDSC sensitive pairs AND a core_score: {len(eval_genes)}")
print(f"  genes: {[sym.get(g, g) for g in eval_genes]}")

df = core[core.ensg_id.isin(eval_genes)].copy()
df = df.merge(flags, on=["model_id", "ensg_id"], how="left")
df = df.merge(mut[["model_id", "ensg_id", "two_hit", "hit1"]],
              on=["model_id", "ensg_id"], how="left")
df = df.merge(sens, on=["model_id", "ensg_id"], how="left")
for c in ["has_driver_alteration", "two_hit", "hit1", "is_sensitive"]:
    df[c] = df[c].fillna(False).astype(bool)

print(f"  evaluation rows: {len(df):,}")

# ---------------------------------------------------------------- Step 3
print()
print("=" * 72)
print("STEP 3 -- hit@20 by arm")
print("=" * 72)


def hit_rate(frame, gate_col=None, invert=False):
    """Per-gene hit@20. invert=True means LOW core_score ranks first (TSG logic)."""
    f = frame.copy()
    base = (1.0 - f["core_score"]) if invert else f["core_score"]
    key = base if gate_col is None else f[gate_col].astype(float) * 1e6 + base
    f["_k"] = key
    f["_r"] = f.groupby("ensg_id")["_k"].rank(ascending=False, method="first")
    known = f.groupby("ensg_id")["is_sensitive"].sum()
    hits = f[f._r <= TOP_K].groupby("ensg_id")["is_sensitive"].sum()
    out = pd.DataFrame({"known": known, "hits": hits}).fillna(0)
    out = out[out.known > 0]
    return (out.hits / out.known).rename(gate_col or "flat")


arms = {
    "flat_high_abundance":      hit_rate(df, None, invert=False),
    "flat_low_abundance_TSG":   hit_rate(df, None, invert=True),
    "production_driver_gate":   hit_rate(df, "has_driver_alteration", invert=True),
    "hit1_only_tsg_mutation":   hit_rate(df, "hit1", invert=True),
    "TWO_HIT_gate":             hit_rate(df, "two_hit", invert=True),
}
res = pd.DataFrame(arms)
print(f"  evaluable genes: {len(res)}")
print()
for name in arms:
    print(f"    {name:26s} {res[name].mean():.6f}")

print()
print("  per-gene detail (hit@20):")
det = res.copy()
det.index = [sym.get(g, g) for g in det.index]
print(det.round(4).to_string())

# ---------------------------------------------------------------- Step 4
print()
print("=" * 72)
print("STEP 4 -- paired bootstrap, TWO_HIT vs production")
print("=" * 72)

paired = (res["TWO_HIT_gate"] - res["production_driver_gate"]).values
point = float(paired.mean())
rng = np.random.default_rng(SEED)
n = len(paired)
boots = np.array([paired[rng.integers(0, n, n)].mean() for _ in range(N_BOOT)])
lo, hi = np.percentile(boots, [2.5, 97.5])
sig = bool(lo > 0 or hi < 0)
print(f"  delta hit@20 (two_hit - production) : {point:+.6f}")
print(f"  95% CI ({N_BOOT} resamples)          : [{lo:+.6f}, {hi:+.6f}]")
print(f"  significant                          : {sig}")

inv_delta = float(res["flat_low_abundance_TSG"].mean() - res["flat_high_abundance"].mean())
print()
print(f"  TSG inversion check (low - high abundance): {inv_delta:+.6f}")
print(f"    -> inversion is {'CORRECT' if inv_delta > 0 else 'NOT SUPPORTED'} on these {len(res)} genes")

if not sig:
    verdict = (f"NO EFFECT -- CI includes 0 on n={len(res)} genes. Two-hit gating does not "
               "measurably beat the production driver gate. Note this is exploratory, not "
               "held-out: the seed-42 split contains zero TSG genes.")
elif point > 0:
    verdict = (f"POSITIVE on n={len(res)} genes (exploratory, not held-out). Two-hit gating "
               "beats the production driver gate. Needs confirmation on a proper held-out "
               "set, which does not currently exist for TSGs.")
else:
    verdict = f"NEGATIVE -- two-hit gating is worse than the production gate on n={len(res)} genes."
print()
print(f"  VERDICT: {verdict}")

results = {
    "question": "Does Knudson two-hit gating rank TSG cell lines better than the production driver gate?",
    "read_only": True,
    "POWER_LIMITATION": (
        "Only 14 TSG genes have GDSC sensitive pairs; the curated 141-gene regime table "
        "contains ZERO loss_of_function genes (95 activation_driven, 46 abundance_tracking), "
        "so the seed-42 held-out split has no TSG genes. This is exploratory, not confirmatory."
    ),
    "two_hit_definition": {
        "hit1": "tsg_hit",
        "hit2": f"max_vaf >= {VAF_LOH_THRESHOLD} OR multi_hit_high_impact OR min_cn < {DEL_THRESHOLD}",
        "lohfraction_excluded_because": "genome-wide per-cell-line scalar, not per-locus",
    },
    "flag_counts": {
        "tsg_rows": int(len(mut)), "hit1": int(mut.hit1.sum()),
        "hit2_vaf": int(mut.h2_vaf.sum()), "hit2_multihit": int(mut.h2_multihit.sum()),
        "hit2_deletion": int(mut.h2_del.sum()), "two_hit": int(mut.two_hit.sum()),
    },
    "n_eval_genes": int(len(res)),
    "eval_genes": [sym.get(g, g) for g in res.index],
    "hit_at_20_mean": {k: round(float(res[k].mean()), 6) for k in arms},
    "per_gene": {sym.get(g, g): {k: round(float(res.loc[g, k]), 4) for k in arms}
                 for g in res.index},
    "paired_bootstrap_twohit_minus_production": {
        "n_resamples": N_BOOT, "point_estimate": round(point, 6),
        "ci_95_low": round(float(lo), 6), "ci_95_high": round(float(hi), 6),
        "significant": sig,
    },
    "tsg_inversion_check": {
        "low_minus_high_abundance": round(inv_delta, 6),
        "inversion_supported": bool(inv_delta > 0),
    },
    "verdict": verdict,
}

out = OUTPUTS / "test_run_tsg_two_hit_results.json"
with open(out, "w") as f:
    json.dump(results, f, indent=2)
print()
print(f"Saved: {out}")
print("No production file was modified.")
