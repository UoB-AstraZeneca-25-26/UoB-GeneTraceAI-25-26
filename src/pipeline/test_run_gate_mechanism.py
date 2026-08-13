"""
test_run_gate_mechanism.py
--------------------------
Read-only diagnostic. Writes no production file.

Two questions about the BOTTOM of the abundance score, both raised by
test_run_abundance_gate.py --label chronos.

PART A -- IS THE FLOOR RESULT WEAK, OR AT THE LABEL'S NOISE FLOOR?
------------------------------------------------------------------
The floor test found dependency rate 0.0139 among (gene, line) pairs where the
gene sits at its expression floor, against a base of 0.0169 -- reported as
"12.5% depletion, a measured weak effect".

That reading assumes 0.0139 is residual biology. It probably is not. Knocking
out a gene that is not transcribed cannot impair fitness, so to first order that
1.4% is the FALSE POSITIVE RATE OF THE LABEL, not a real signal the gate failed
to remove. If so, the gate is not achieving 12.5% depletion against a true
baseline -- it is achieving near-complete depletion and hitting a floor imposed
by screen noise.

Estimated with negative-control gene families that are not expressed in cancer
cell lines and are not essential anywhere: olfactory receptors (OR[0-9]*),
keratin-associated proteins (KRTAP*), taste receptors (TAS2R*), vomeronasal
receptors (VN1R*). These are the standard non-expressed controls in CRISPR
screen QC. Their rate of `essentiality <= DEP_THRESHOLD` across all lines is a
clean estimate of the label's FPR, with no reference to the abundance score.

    control FPR ~ 0.014  -> the floor result is AT the noise ceiling and should
                            be reported as complete depletion within label noise
    control FPR << 0.014 -> 1.4% at the floor is real residual signal and the
                            original "weak effect" reading was right

PART B -- WHY THE MEAN BEATS RELIABILITY-WEIGHTED STOUFFER AT THE BOTTOM
------------------------------------------------------------------------
Proposed mechanism: protein detection is CATEGORICAL evidence of expression. A
protein that was quantified cannot have zero transcript, so a floor RNA reading
on a line with detected protein is more likely a measurement artefact than true
absence. The two rules treat that case very differently:

    floor RNA (p=0.0007), mid protein (p=0.50)
      mean      (0.0007 + 0.50)/2                    = 0.25   mid-low
      Stouffer  (0.987*-3.2 + 0.381*0)/1.21 = -2.61  -> 0.005  very bottom

w_P = 0.381 means a protein reading barely resists a floor RNA reading. The mean
happens to respect the detection evidence; reliability-weighted Stouffer
overrides it.

Test: among RNA-floor lines, split by whether protein was detected and compare
dependency rates. CONTROL: run the identical split on a mid-RNA stratum. If
protein-detected lines are more dependent everywhere, that is just a property of
which lines get profiled (they are bigger, faster-growing, better-screened). The
mechanism is confirmed only if the gap is specific to the floor.

If confirmed, the implication is a RULE, not a weight: a line with detected
protein cannot be gated out on RNA alone. That connects to the MNAR work --
detection is informative -- and predicts equal-weight Stouffer will recover the
gate-region deficit.

SIGN CONVENTION
---------------
chronos `essentiality`: more negative = more essential. Positive (dependent)
= essentiality <= DEP_THRESHOLD. Verified in-script against core-essential genes.

Outputs:
    src/pipeline/outputs/test_run_gate_mechanism_results.json

Run:
    python src/pipeline/test_run_gate_mechanism.py
"""
import json
import re
import warnings
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import fisher_exact, norm, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

OUTPUTS = Path("src/pipeline/outputs")
REF_DIR = Path("reference")
DATA_CLEAN = Path("data/parquet/data_clean")
CLEANED_TRACK = Path("cleaned_track_data")

DEP_THRESHOLD = -0.5
N_GENES = 3000
SEED = 42
results = {}

# ============================================================ shared load
print("=" * 78)
print("STEP 1 -- chronos labels and gene symbols")
print("=" * 78)

gl = pd.read_parquet(REF_DIR / "gene_lookup.parquet",
                     columns=["ensg_id", "hgnc_symbol", "biotype",
                              "hgnc_status", "uniprot_ids"])
gl["ensg_id"] = gl.ensg_id.astype("string").str.split(".").str[0].str.lower()
sym_of = dict(zip(gl.ensg_id, gl.hgnc_symbol.astype(str).str.upper()))
uni = gl[(gl.biotype == "protein_coding") & (gl.hgnc_status == "Approved")]
valid_ensg = set(uni.ensg_id)
u2e = {}
for e, u in zip(uni.ensg_id, uni.uniprot_ids):
    if isinstance(u, str):
        for acc in u.split("|"):
            acc = acc.strip().lower()
            if acc:
                u2e.setdefault(acc, e)

ch = pd.read_parquet("validation/prepared/chronos_long.parquet")
ch["ensg_id"] = ch.ensg_id.astype(str).str.split(".").str[0].str.lower()
core = [e for e, s in sym_of.items() if s in {"RPL13A", "RPS6", "EIF4A3", "PSMB2"}]
rest = [e for e, s in sym_of.items() if s in {"ALB", "KLK4", "CD3E", "MYOD1"}]
med = ch.groupby("ensg_id").essentiality.median()
m_core, m_rest = med.reindex(core).median(), med.reindex(rest).median()
print(f"  sign check: core-essential {m_core:+.2f} vs tissue-restricted {m_rest:+.2f}")
if not (m_core < m_rest and m_core < 0):
    raise SystemExit("chronos sign convention check FAILED. Aborting.")

ib = pd.read_parquet("validation/prepared/id_bridge.parquet")
ch = ch.merge(ib, on="sanger_model_id", how="inner")
ch["model_id"] = ch.model_id.str.lower()
ch["dep"] = ch.essentiality <= DEP_THRESHOLD
print(f"  chronos: {ch.ensg_id.nunique():,} genes x {ch.model_id.nunique():,} lines")
print(f"  global dependency prevalence: {ch.dep.mean():.4f}")

# ============================================================ PART A
print()
print("=" * 78)
print("PART A -- NEGATIVE-CONTROL ESTIMATE OF THE LABEL'S FALSE POSITIVE RATE")
print("=" * 78)

# NOTE ON THE THRESHOLD, found while running this: on this file's scale
# DEP_THRESHOLD = -0.5 sits at roughly the 12th percentile of ALL (gene, line)
# pairs. The binary label is therefore closer to "bottom ~12% of fitness scores"
# than to an absolute dependency call. That does not affect the decile profile or
# the partial AUROC, which are within-gene relative comparisons, but it does mean
# the "12.2% prevalence" figure is a property of the cut, not of biology.
pct_of_thresh = float((ch.essentiality <= DEP_THRESHOLD).mean())
print(f"  NOTE: threshold {DEP_THRESHOLD} sits at the {100*pct_of_thresh:.1f}th "
      f"percentile of all pairs -- this label is a quantile, not an absolute call")

# PARALOGY WARNING. Multi-targeting guides in near-identical gene families produce
# APPARENT fitness effects: one guide cuts many loci, the cell suffers, and the
# gene scores essential without being essential. KRTAP and VN1R are exactly such
# families, so they are reported but excluded from the FPR estimate.
FAMILIES = {
    # name: (matcher, usable_as_control)
    "olfactory receptors (OR*)": (lambda s: bool(re.fullmatch(r"OR\d+[A-Z]\d*", s)), True),
    "taste receptors (TAS2R*)": (lambda s: s.startswith("TAS2R"), True),
    "keratin-associated (KRTAP*)": (lambda s: s.startswith("KRTAP"), False),
    "vomeronasal (VN1R*)": (lambda s: s.startswith("VN1R"), False),
}
in_chr = set(ch.ensg_id.unique())
print()
print(f"  {'family':<30} {'genes':>6} {'in screen':>10} {'pairs':>9} "
      f"{'dep rate':>10}  {'usable':>7}")
partA, clean, absent = {}, set(), []
for name, (fn, usable) in FAMILIES.items():
    s_all = {e for e, sy in sym_of.items() if sy and fn(sy)}
    s = s_all & in_chr
    sub = ch[ch.ensg_id.isin(s)]
    rate = float(sub.dep.mean()) if len(sub) else float("nan")
    print(f"  {name:<30} {len(s_all):>6} {len(s):>10} {len(sub):>9,} "
          f"{rate:>10.4f}  {'yes' if usable else 'PARALOG':>7}")
    partA[name] = {"genes_in_lookup": len(s_all), "genes_in_screen": len(s),
                   "pairs": int(len(sub)), "dep_rate": rate, "usable": usable}
    if len(s) == 0:
        absent.append(name)
    elif usable:
        clean |= s

for name in absent:
    print(f"  !! {name}: 0 of its genes are in the screen at all -- the library")
    print(f"     excludes them, so the strongest available control is unavailable.")

FLOOR_RATE = 0.0139        # measured in test_run_abundance_gate.py --label chronos
FLOOR_BASE = 0.0169
sub_clean = ch[ch.ensg_id.isin(clean)]
ctrl_fpr = float(sub_clean.dep.mean()) if len(sub_clean) else float("nan")
n_clean_genes = len(clean)

print()
print(f"  genome-wide prevalence                        : {ch.dep.mean():.4f}")
print(f"  paralogy-contaminated families (KRTAP, VN1R)  : "
      f"{float(ch[ch.ensg_id.isin({e for e, sy in sym_of.items() if sy.startswith(('KRTAP','VN1R'))} & in_chr)].dep.mean()):.4f}"
      f"   <- inflated, not an FPR")
print(f"  USABLE negative-control FPR ({n_clean_genes} genes)        : {ctrl_fpr:.4f}")
print(f"  dependency rate at the RNA floor              : {FLOOR_RATE:.4f}")
print()
if not np.isfinite(ctrl_fpr) or n_clean_genes < 10:
    verdict_a = "inconclusive__no_usable_negative_control"
    print(f"  -> INCONCLUSIVE. Only {n_clean_genes} usable control genes remain after")
    print(f"     excluding paralogous families, and the olfactory receptors that")
    print(f"     would normally carry this test are not in the screen library.")
    print(f"     The label's FPR cannot be estimated reliably from gene families here.")
elif FLOOR_RATE <= ctrl_fpr:
    verdict_a = "floor_result_is_at_or_below_the_label_noise_floor"
    print(f"  -> The floor rate ({FLOOR_RATE:.4f}) is AT OR BELOW the cleanest")
    print(f"     available estimate of the label's false-positive rate ({ctrl_fpr:.4f}).")
    print(f"     Depletion at the floor is complete within label noise, and")
    print(f"     '12.5% depletion vs base' understates what the gate achieves.")
    print(f"     CAVEAT: {n_clean_genes} control genes only, and the estimate")
    print(f"     ranges {min(v['dep_rate'] for v in partA.values() if v['pairs']):.4f}"
          f"-{max(v['dep_rate'] for v in partA.values() if v['pairs']):.4f} across")
    print(f"     families, so it is an order-of-magnitude statement, not a number.")
else:
    verdict_a = "floor_result_is_real_residual_signal"
    print(f"  -> Control FPR ({ctrl_fpr:.4f}) is below the floor rate "
          f"({FLOOR_RATE:.4f}),")
    print(f"     so part of the floor rate is not label noise and the 'weak effect'")
    print(f"     reading stands.")
results["part_a"] = {"families": partA, "usable_control_fpr": ctrl_fpr,
                     "n_usable_control_genes": n_clean_genes,
                     "threshold_percentile": pct_of_thresh,
                     "genome_prevalence": float(ch.dep.mean()),
                     "floor_rate_measured": FLOOR_RATE,
                     "verdict": verdict_a}

# ============================================================ PART B
print()
print("=" * 78)
print("PART B -- DOES PROTEIN DETECTION RESCUE FLOOR-RNA LINES?")
print("=" * 78)

expr_names = pq.ParquetFile(DATA_CLEAN / "depmap_expr_clean.parquet").schema_arrow.names
expr_col = {c.split(".")[0].lower(): c for c in expr_names if c != "index"}
prot_names = pq.ParquetFile(CLEANED_TRACK / "proteomics.parquet").schema_arrow.names
key_col = "model_id" if "model_id" in prot_names else "depmap_id"
ccle_by_ensg = {}
for c in prot_names:
    if c == key_col:
        continue
    a = c.strip().lower()
    e = u2e.get(a) or u2e.get(a.split("-")[0])
    if e in valid_ensg:
        ccle_by_ensg.setdefault(e, []).append(c)
con = duckdb.connect(str(OUTPUTS / "celllineselector.db"), read_only=True)
pc_names = [r[0] for r in con.execute("DESCRIBE procan_proteomics").fetchall()]
PC_META = {"gdsc_model_name", "sanger_model_id", "model_id",
           "matched_via", "n_model_id", "is_ambiguous"}
procan_by_ensg = {}
for c in pc_names:
    if c in PC_META:
        continue
    e = u2e.get(c.strip().lower())
    if e in valid_ensg:
        procan_by_ensg.setdefault(e, []).append(c)

rng0 = np.random.default_rng(SEED)
cands = sorted(set(ch.ensg_id.unique()) & set(expr_col)
               & (set(ccle_by_ensg) | set(procan_by_ensg)))
sweep = sorted(rng0.choice(cands, size=min(N_GENES, len(cands)), replace=False))
print(f"  genes with chronos + RNA + protein columns: {len(cands):,}  "
      f"(sampled {len(sweep):,})")

prof = pd.read_parquet(REF_DIR / "depmap_profiles.parquet")
rp = prof[prof.datatype == "rna"][["profileid", "modelid"]].copy()
rp["model_id"] = rp.modelid.str.lower()
rp["profileid"] = rp.profileid.astype(str)
E = pq.read_table(DATA_CLEAN / "depmap_expr_clean.parquet",
                  columns=["index"] + [expr_col[g] for g in sweep]).to_pandas()
if "index" not in E.columns:
    E = E.reset_index()
E["profileid"] = E["index"].astype(str)
E = E.drop(columns=["index"]).merge(rp[["profileid", "model_id"]],
                                    on="profileid", how="inner")
E = E.drop(columns=["profileid"]).groupby("model_id").mean()
E.columns = [c.split(".")[0].lower() for c in E.columns]

need_ccle = sorted({a for g in sweep for a in ccle_by_ensg.get(g, [])})
C = pq.read_table(CLEANED_TRACK / "proteomics.parquet",
                  columns=[key_col] + need_ccle).to_pandas()
C.index = C[key_col].str.lower().values
C = C.drop(columns=[key_col]).astype(float).groupby(level=0).mean()
need_pc = sorted({a for g in sweep for a in procan_by_ensg.get(g, [])})
q = ", ".join(f'"{a}"' for a in need_pc)
Pm = con.execute(f"SELECT model_id, {q} FROM procan_proteomics "
                 "WHERE model_id IS NOT NULL").df()
Pm.index = Pm.model_id.str.lower().values
Pm = Pm.drop(columns=["model_id"]).astype(float).groupby(level=0).mean()
con.close()

dep_by_gene = {g: s.set_index("model_id")["dep"]
               for g, s in ch[ch.ensg_id.isin(sweep)][
                   ["ensg_id", "model_id", "dep"]].groupby("ensg_id")}

# Pooled 2x2 counts, floor stratum and mid-RNA control stratum
tot = {"floor": {"det": [0, 0], "nodet": [0, 0]},
       "mid": {"det": [0, 0], "nodet": [0, 0]}}
per_gene = []
for g in sweep:
    if g not in dep_by_gene or g not in E.columns:
        continue
    e = E[g].dropna()
    dep = dep_by_gene[g]
    lines = [m for m in e.index if m in dep.index]
    if len(lines) < 100:
        continue
    e = e.reindex(lines)
    y = dep.reindex(lines).astype(bool)

    det = pd.Series(False, index=lines)
    accs = [a for a in ccle_by_ensg.get(g, []) if a in C.columns]
    if accs:
        s = C[accs].mean(axis=1).dropna()
        det.loc[[m for m in lines if m in s.index]] = True
    accs = [a for a in procan_by_ensg.get(g, []) if a in Pm.columns]
    if accs:
        s = Pm[accs].mean(axis=1).dropna()
        det.loc[[m for m in lines if m in s.index]] = True

    floor = e == e.min()
    r = e.rank(pct=True, method="min")
    mid = (r >= 0.45) & (r <= 0.55)          # control stratum, same profiling mix
    for nm, mask in [("floor", floor), ("mid", mid)]:
        for dk, dmask in [("det", det), ("nodet", ~det)]:
            sel = mask & dmask
            tot[nm][dk][0] += int(y[sel].sum())
            tot[nm][dk][1] += int(sel.sum())
    if floor.sum() >= 20 and (floor & det).sum() >= 5 and (floor & ~det).sum() >= 5:
        per_gene.append({
            "ensg_id": g,
            "floor_det": float(y[floor & det].mean()),
            "floor_nodet": float(y[floor & ~det].mean()),
        })

print()
print(f"  {'stratum':<22} {'protein':<10} {'lines':>9} {'dependent':>10} {'rate':>9}")
partB = {}
for nm in ("floor", "mid"):
    for dk, disp in (("det", "detected"), ("nodet", "not detected")):
        k, n = tot[nm][dk]
        rate = k / n if n else np.nan
        print(f"  {nm:<22} {disp:<10} {n:>9,} {k:>10,} {rate:>9.4f}")
        partB[f"{nm}_{dk}"] = {"lines": n, "dependent": k, "rate": float(rate)}
    k1, n1 = tot[nm]["det"]
    k0, n0 = tot[nm]["nodet"]
    orr, pf = fisher_exact([[k1, n1 - k1], [k0, n0 - k0]])
    ratio = (k1 / n1) / (k0 / n0) if n0 and k0 else np.nan
    print(f"  {'':<22} {'-> ratio det/nodet':<10} {ratio:>28.3f}   "
          f"OR {orr:.3f}  Fisher p {pf:.3g}")
    partB[f"{nm}_ratio"] = float(ratio)
    partB[f"{nm}_fisher_p"] = float(pf)
    print()

floor_ratio, mid_ratio = partB["floor_ratio"], partB["mid_ratio"]
specific = floor_ratio > mid_ratio * 1.15 and partB["floor_fisher_p"] < 0.05
print(f"  floor stratum ratio {floor_ratio:.3f}  vs  mid-RNA control ratio "
      f"{mid_ratio:.3f}")
if specific:
    verdict_b = "detection_rescues_floor_lines__mechanism_confirmed"
    print(f"  -> The gap is SPECIFIC TO THE FLOOR. Protein detection carries")
    print(f"     information that RNA-floor does not, exactly where the gate acts.")
    print(f"     Implication is a RULE: a line with detected protein should not be")
    print(f"     gated out on RNA alone.")
else:
    verdict_b = "no_floor_specific_effect__profiling_confound_only"
    print(f"  -> NOT specific to the floor. Protein-detected lines differ in the")
    print(f"     mid stratum too, so this is a property of which lines get")
    print(f"     profiled, not evidence that detection rescues floor readings.")
    print(f"     The proposed mechanism is NOT supported.")

if per_gene:
    PG = pd.DataFrame(per_gene)
    d = (PG.floor_det - PG.floor_nodet).dropna()
    print(f"\n  per-gene check ({len(PG)} genes with >=5 lines in both cells):")
    print(f"    median (rate_detected - rate_not_detected) = {d.median():+.4f}   "
          f"p = {wilcoxon(d).pvalue:.3g}   positive in {100*(d > 0).mean():.1f}%")
    partB["per_gene_median_diff"] = float(d.median())
    partB["per_gene_p"] = float(wilcoxon(d).pvalue)
    partB["per_gene_n"] = int(len(PG))

results["part_b"] = {**partB, "verdict": verdict_b}

print()
print("=" * 78)
print("VERDICT")
print("=" * 78)
print(f"  A: {verdict_a}")
print(f"  B: {verdict_b}")

out = OUTPUTS / "test_run_gate_mechanism_results.json"
out.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
