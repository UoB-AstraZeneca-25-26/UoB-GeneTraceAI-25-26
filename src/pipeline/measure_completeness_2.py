"""
measure_completeness_2.py
--------------------------
Q1, Q4, Q5, Q7 -- follow-up measurements after the 0-vs-unknown bug fix in
rank_convergent.py (see its evidence() docstring for what was wrong and why).

Q1  measured-negative (0, but layer WAS measured) vs true-unknown (0, layer NOT
    measured), as a fraction of evidence-vector cells, for representative
    queries. Now measurable because evidence() exposes `<layer>__measured`.
Q4  does the completeness/study-frequency correlation (rho=0.60, prior script)
    hold WITHIN lineage, or is it a cross-lineage artefact?
Q5  is the mutation-layer correlation (rho=0.72) driven by per-line sequencing
    depth, or by how many lines per lineage entered the mutation pipeline?
Q7  within the confirmed subset (n_confirm >= 1), what is the largest tied
    band? Does ranking add anything even there?

    python src/pipeline/measure_completeness_2.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import multi_gene_kleene as K   # noqa: E402
import rank_convergent as RC    # noqa: E402

OUT = ROOT / "src" / "pipeline" / "outputs"
RNG = np.random.default_rng(42)
res = {}

print("=" * 78)
print("FOLLOW-UP MEASUREMENTS -- Q1, Q4, Q5, Q7")
print("=" * 78)

E = K.expression()
fe = K.gene_validity()
CN = RC.cn_matrix()
lin = K.lineage()

pool = [g for g in E.columns if fe.get(g, 0) >= K.SILENT_FRAC and g in CN.columns]
genes = list(RNG.choice(pool, 15, replace=False))
for s in ["TP53", "PTEN", "RB1", "CDKN2A", "BRAF", "EGFR", "KRAS", "ERBB2"]:
    e, _ = K.resolve_gene(s)
    if e in pool and e not in genes:
        genes.append(e)

# ---------------------------------------------------------------- Q1
print("\nQ1 -- MEASURED-NEGATIVE (0, measured) vs TRUE-UNKNOWN (0, not measured)")
print("  Now distinguishable after the evidence() fix.")
print()
all_ev = []
for g in genes:
    terms = K.parse_terms(not_expressed=[g])
    q = K.query(terms, k_alternatives=0)
    surv = q["buckets"]["matches"]
    if len(surv) < 20:
        continue
    V = RC.evidence(g, surv, "loss")
    V["gene"] = g
    all_ev.append(V)
A = pd.concat(all_ev)
layers = RC.LOSS_LAYERS
print(f"  {'layer':<22}{'confirm +1':>12}{'measured-0':>13}{'unknown-0':>12}{'contradict -1':>15}")
res["q1"] = {}
for l in layers:
    v, m = A[l], A[f"{l}__measured"]
    n = len(A)
    confirm = int((v == 1).sum())
    contradict = int((v == -1).sum())
    meas0 = int(((v == 0) & m).sum())
    unk0 = int(((v == 0) & ~m).sum())
    print(f"  {l:<22}{confirm:>11,} {meas0:>12,} {unk0:>11,} {contradict:>14,}")
    res["q1"][l] = {"n": n, "confirm": confirm, "measured_zero": meas0,
                    "unknown_zero": unk0, "contradict": contradict,
                    "measured_zero_frac_of_zeros": meas0 / max(1, meas0 + unk0)}
print()
print("  'measured-0' = layer measured, did not confirm (a real negative).")
print("  'unknown-0'  = layer never measured. Both used to be indistinguishable.")
tm = res["q1"]["truncating_mutation"]
print(f"\n  truncating_mutation: of {tm['measured_zero']+tm['unknown_zero']:,} zeros, "
      f"{tm['measured_zero_frac_of_zeros']:.1%} were actually measured-negative.")

# ---------------------------------------------------------------- Q4
print("\nQ4 -- DOES THE FAME CORRELATION HOLD WITHIN LINEAGE?")
prof = pd.read_parquet(ROOT / "reference" / "depmap_profiles.parquet")
prof["model_id"] = prof.modelid.astype(str).str.lower()
study_freq = prof.groupby("model_id").datatype.nunique()

import measure_completeness as MC  # reuses completeness_for
all_c = []
for g in genes:
    terms = K.parse_terms(not_expressed=[g])
    q = K.query(terms, k_alternatives=0)
    surv = q["buckets"]["matches"]
    if len(surv) < 20:
        continue
    D = MC.completeness_for(g, surv)
    all_c.append(D)
C = pd.concat(all_c).reset_index().drop_duplicates("model_id").set_index("model_id")
C["study_freq"] = study_freq.reindex(C.index)
C["lineage"] = lin.reindex(C.index).lineage
C = C.dropna(subset=["study_freq", "lineage"])

pooled_rho, pooled_p = spearmanr(C.n_measured, C.study_freq)
print(f"  pooled (all lineages)   : rho={pooled_rho:.4f}  p={pooled_p:.3g}  n={len(C):,}")

within = []
for lg, grp in C.groupby("lineage"):
    if len(grp) < 30:
        continue
    r, p = spearmanr(grp.n_measured, grp.study_freq)
    within.append({"lineage": lg, "n": len(grp), "rho": r, "p": p})
W = pd.DataFrame(within).sort_values("rho")
print(f"\n  within-lineage (n>=30 lines), {len(W)} lineages tested:")
print(f"  {'lineage':<26}{'n':>6}{'rho':>9}{'p':>12}")
for _, r in W.iterrows():
    print(f"  {r['lineage']:<26}{int(r['n']):>6}{r['rho']:>9.3f}{r['p']:>12.3g}")
print()
print(f"  median within-lineage rho: {W.rho.median():.4f}")
print(f"  pooled rho               : {pooled_rho:.4f}")
if W.rho.median() < pooled_rho * 0.6:
    q4v = "MOSTLY_CROSS_LINEAGE"
    print("  -> the correlation is substantially a CROSS-lineage effect: some")
    print("     lineages are more completely measured than others as a block.")
    print("     A researcher querying within one lineage sees much less fame")
    print("     bias than the pooled rho suggests.")
else:
    q4v = "HOLDS_WITHIN_LINEAGE"
    print("  -> the correlation persists within lineage too -- not purely a")
    print("     between-lineage artefact. Fame bias operates even inside a")
    print("     single cancer type.")
res["q4"] = {"pooled_rho": float(pooled_rho), "pooled_p": float(pooled_p),
            "n_lineages_tested": int(len(W)),
            "within_lineage_rho_median": float(W.rho.median()),
            "per_lineage": W.to_dict("records"), "verdict": q4v}

# ---------------------------------------------------------------- Q5
print("\nQ5 -- MUTATION rho=0.72: SEQUENCING DEPTH OR PANEL INCLUSION?")
print("  Proxy for depth: WES profile count per line (>1 = replicate WES runs).")
print("  Proxy for panel inclusion: lineage size in the mutation table.")
wes = prof[prof.datatype == "wes"].groupby("model_id").size()
mut = pd.read_parquet(ROOT / "cleaned_track_data" / "mutations_collapsed.parquet",
                      columns=["model_id"])
mut["model_id"] = mut.model_id.astype(str).str.lower()
seq_lines = set(mut.model_id.unique())

D5 = C.copy()
D5["wes_profiles"] = wes.reindex(D5.index).fillna(0)
D5["in_mutation_table"] = D5.index.isin(seq_lines).astype(int)
lineage_n = D5.groupby("lineage").size()
D5["lineage_size"] = D5.lineage.map(lineage_n)

r_depth, p_depth = spearmanr(D5.mutation.astype(int), D5.wes_profiles)
r_lineage, p_lineage = spearmanr(D5.mutation.astype(int), D5.lineage_size)
print(f"\n  mutation-measured vs WES replicate count : rho={r_depth:.4f}  p={p_depth:.3g}")
print(f"  mutation-measured vs lineage panel size  : rho={r_lineage:.4f}  p={p_lineage:.3g}")
# WES profile count is near-binary (0 or 1 for almost every line), so this
# mostly degenerates to "was WES run at all" -- report that explicitly.
print(f"\n  WES profile count distribution: {D5.wes_profiles.value_counts().to_dict()}")
if D5.wes_profiles.max() <= 1:
    print("  -> no line in this set has more than one WES profile, so 'depth' in")
    print("     the sense of replicate coverage cannot be distinguished from")
    print("     simple inclusion here. The available data speaks to INCLUSION")
    print("     (was this line ever sequenced), not to depth of sequencing.")
    q5v = "CANNOT_SEPARATE_DEPTH_FROM_INCLUSION_WITH_AVAILABLE_DATA"
else:
    q5v = "MEASURED"
res["q5"] = {"rho_vs_wes_replicates": float(r_depth), "p_vs_wes": float(p_depth),
            "rho_vs_lineage_size": float(r_lineage), "p_vs_lineage": float(p_lineage),
            "wes_profile_distribution": {str(k): int(v) for k, v in
                                         D5.wes_profiles.value_counts().items()},
            "verdict": q5v}
print(f"\n  VERDICT: {q5v}")
if q5v.startswith("CANNOT"):
    print("  The remedy this implies: the fame proxy is really an ASCERTAINMENT")
    print("  proxy (which lines got included in a sequencing panel at all), not a")
    print("  technical-depth proxy. Fixing it means expanding panel coverage to")
    print("  under-profiled lines, not re-sequencing existing ones more deeply.")

# ---------------------------------------------------------------- Q7
print("\nQ7 -- WITHIN THE CONFIRMED SUBSET, DOES RANKING ADD ANYTHING?")
rows7 = []
for g in genes:
    for direction, kwarg in (("loss", {"not_expressed": [g]}),
                             ("gain", {"expressed": [g]})):
        terms = K.parse_terms(**kwarg)
        q = K.query(terms, k_alternatives=0)
        surv = q["buckets"]["matches"]
        if len(surv) < 20:
            continue
        R = RC.rank(g, surv, direction)
        conf = R[R.n_confirm >= 1]
        if len(conf) < 5:
            continue
        biggest = int(conf.band_size[conf.band.isin(
            conf.band.value_counts().head(1).index)].max())
        # recompute banding fraction within JUST the confirmed subset
        big_band_size = conf.groupby("band").size().max()
        rows7.append({"gene": g, "direction": direction,
                      "n_survivors": len(surv), "n_confirmed": len(conf),
                      "pct_confirmed": len(conf) / len(surv),
                      "largest_band_in_confirmed": int(big_band_size),
                      "largest_band_frac_in_confirmed": big_band_size / len(conf)})
R7 = pd.DataFrame(rows7)
print(f"  {'gene':<18}{'dir':<6}{'survivors':>10}{'confirmed':>11}{'%conf':>8}"
      f"{'big band':>10}{'%band':>8}")
for _, r in R7.sort_values("largest_band_frac_in_confirmed", ascending=False).head(14).iterrows():
    print(f"  {r['gene'][:16]:<18}{r['direction']:<6}{r['n_survivors']:>10,}"
          f"{r['n_confirmed']:>11,}{r['pct_confirmed']:>8.1%}"
          f"{r['largest_band_in_confirmed']:>10,}{r['largest_band_frac_in_confirmed']:>8.1%}")
print()
print(f"  median %% of survivors with >=1 confirming layer: "
      f"{R7.pct_confirmed.median():.1%}")
print(f"  median largest tied band WITHIN the confirmed subset: "
      f"{R7.largest_band_frac_in_confirmed.median():.1%}")
res["q7"] = {
    "median_pct_confirmed": float(R7.pct_confirmed.median()),
    "median_largest_band_frac_in_confirmed": float(R7.largest_band_frac_in_confirmed.median()),
    "n_queries": int(len(R7)),
}
if res["q7"]["median_largest_band_frac_in_confirmed"] >= 0.60:
    print("\n  Also near-total ties, even restricted to the confirmed subset.")
    print("  The honest output is a FLAG, not a ranking: 'N lines qualify; K")
    print("  have additional corroborating evidence', with no order claimed")
    print("  within that K.")
    q7v = "NO_ORDER_EVEN_WITHIN_CONFIRMED"
else:
    print("\n  The confirmed subset DOES discriminate meaningfully -- ranking")
    print("  within it is defensible, unlike ranking the full survivor set.")
    q7v = "ORDER_MEANINGFUL_WITHIN_CONFIRMED"
res["q7"]["verdict"] = q7v

p = OUT / "completeness_measurement_2.json"
p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {p}")
