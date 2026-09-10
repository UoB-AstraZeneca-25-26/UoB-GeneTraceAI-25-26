"""
combined_score_b0.py
--------------------
TASK B0 -- the precondition, and the kill switch.

Three questions, in order:
  1. How many gene-drug pairs exist where the drug targets the gene and GDSC2
     covers both?
  2. Reproduce the flat MARGINAL on exactly that pair set -- not on the
     141-target set the gate audit used.
  3. How many pairs survive three-way coverage: expression for A, alteration for
     B, drug response for the pair?

KILL SWITCH: fewer than 20 usable pairs -> stop, report, ship Task A.

WHAT THIS SCRIPT WILL NOT DO
----------------------------
It will not choose the modifier genes B. Per B1 the A-B list is a USER INPUT with
citations, pre-registered before any measurement. Searching for modifiers that
work is selection by outcome -- the `.head(30)` error at larger scale -- and it
would invalidate the exercise regardless of what came out.

So the B-side is reported as a CAPACITY: for each covered A, how many candidate
modifier genes could in principle support an interaction test at this coverage.
That is a statement about what the data could support, not a selection.

    python src/pipeline/combined_score_b0.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import multi_gene_kleene as K  # noqa: E402

GDSC = ROOT / "validation" / "prepared" / "gdsc_scored_ready.parquet"
OUT = ROOT / "src" / "pipeline" / "outputs"
MIN_PAIRS = 20
MIN_LINES_PER_PAIR = 100
RNG = np.random.default_rng(42)
res = {"kill_switch_threshold": MIN_PAIRS}

print("=" * 78)
print("TASK B0 -- PRECONDITION FOR A COMBINED MULTI-GENE SCORE")
print("=" * 78)

E = K.expression()
fe = K.gene_validity()
sequenced, alt, _mod = K.alterations()

gd = pd.read_parquet(GDSC, columns=["model_id", "drug_id", "drug_name",
                                    "target_ensg", "target_symbol", "AUC",
                                    "sensitive"])
gd["model_id"] = gd.model_id.astype(str).str.lower()
gd["target_ensg"] = (gd.target_ensg.astype("string").str.split(".").str[0]
                     .str.lower())

# ---------------------------------------------------------------- 1. pairs
print("\n1. GENE-DRUG PAIRS")
pairs = (gd.groupby(["target_ensg", "drug_id", "drug_name", "target_symbol"])
         .agg(n_lines=("model_id", "nunique"),
              n_sensitive=("sensitive", "sum"))
         .reset_index())
print(f"  raw (target gene, drug) pairs in GDSC2 : {len(pairs):,}")
print(f"  distinct target genes                  : {pairs.target_ensg.nunique()}")
print(f"  distinct drugs                         : {pairs.drug_id.nunique()}")

expr_lines = set(E.index)
gd_cov = gd[gd.model_id.isin(expr_lines)]
pairs_cov = (gd_cov.groupby(["target_ensg", "drug_id", "drug_name", "target_symbol"])
             .agg(n_lines=("model_id", "nunique"),
                  n_sensitive=("sensitive", "sum"))
             .reset_index())
pairs_cov = pairs_cov[(pairs_cov.n_lines >= MIN_LINES_PER_PAIR)
                      & (pairs_cov.n_sensitive >= 10)]
pairs_cov = pairs_cov[pairs_cov.target_ensg.isin(E.columns)]
print(f"\n  after requiring RNA coverage for A, >={MIN_LINES_PER_PAIR} lines and "
      f">=10 sensitive:")
print(f"    usable (A, drug) pairs               : {len(pairs_cov):,}")
print(f"    distinct target genes                : {pairs_cov.target_ensg.nunique()}")
valid_A = pairs_cov[pairs_cov.target_ensg.map(lambda g: fe.get(g, 0)) >= K.SILENT_FRAC]
print(f"    ... whose A passes the validity guard: {len(valid_A):,} "
      f"({valid_A.target_ensg.nunique()} genes)")
res["pairs"] = {
    "raw": int(len(pairs)), "raw_genes": int(pairs.target_ensg.nunique()),
    "raw_drugs": int(pairs.drug_id.nunique()),
    "with_rna_and_min_lines": int(len(pairs_cov)),
    "with_valid_A": int(len(valid_A)),
    "valid_A_genes": int(valid_A.target_ensg.nunique())}

# ---------------------------------------------------------------- 2. marginal
print("\n2. THE MARGINAL, REPRODUCED ON THIS EXACT PAIR SET")
print("   Does expression of the target gene predict sensitivity to a drug")
print("   against that target? Per-gene prevalence-matched binarisation, as in")
print("   gate_audit/07 -- threshold transfer across assays is invalid.")
print()

QS = (0.05, 0.10, 0.25)
rows = []
for _, r in valid_A.iterrows():
    g, did = r.target_ensg, r.drug_id
    sub = gd_cov[(gd_cov.target_ensg == g) & (gd_cov.drug_id == did)]
    lines = sorted(set(sub.model_id) & expr_lines)
    if len(lines) < MIN_LINES_PER_PAIR:
        continue
    x = E[g].reindex(lines).values
    auc = sub.set_index("model_id")["AUC"].reindex(lines).values
    ok = ~np.isnan(x) & ~np.isnan(auc)
    x, auc, n = x[ok], auc[ok], int(ok.sum())
    if n < MIN_LINES_PER_PAIR:
        continue
    rec = {"target_ensg": g, "symbol": r.target_symbol, "drug_id": int(did),
           "drug": r.drug_name, "n_lines": n}
    for q in QS:
        k = max(10, int(round(n * q)))
        if n - k < 10:
            continue
        cut = np.partition(auc, k - 1)[k - 1]
        pos = auc <= cut                       # most sensitive q fraction
        base = pos.mean()
        # decile-1 ratio on the expression score: does the LOW-expression decile
        # deplete for sensitivity, as an abundance gate would predict?
        order = np.argsort(x, kind="mergesort")
        d1 = order[:max(5, n // 10)]
        rec[f"d1_ratio_q{q}"] = float(pos[d1].mean() / base) if base > 0 else np.nan
        rec[f"rho_q{q}"] = float(np.corrcoef(x, auc.astype(float))[0, 1])
    rows.append(rec)

M = pd.DataFrame(rows)
print(f"  pairs evaluated: {len(M):,}   median lines/pair {M.n_lines.median():.0f}")
print()
# Standing rule: cluster bootstrap over CELL LINES, never over pairs. A pair
# bootstrap would treat 294 pairs built from ~640 shared lines as 294
# independent units, which is the same pseudo-replication fault the gate audit
# found (design effect 2.4x) and the similarity work corrected for.
per_pair = []
for _, r in valid_A.iterrows():
    g, did = r.target_ensg, r.drug_id
    sub = gd_cov[(gd_cov.target_ensg == g) & (gd_cov.drug_id == did)]
    lines = sorted(set(sub.model_id) & expr_lines)
    if len(lines) < MIN_LINES_PER_PAIR:
        continue
    x = E[g].reindex(lines).values
    auc = sub.set_index("model_id")["AUC"].reindex(lines).values
    ok = ~np.isnan(x) & ~np.isnan(auc)
    if ok.sum() < MIN_LINES_PER_PAIR:
        continue
    per_pair.append((np.array(lines)[ok], x[ok], auc[ok]))

all_lines = sorted({m for ls, _, _ in per_pair for m in ls})
line_ix = {m: i for i, m in enumerate(all_lines)}


def d1_median(q, weights=None):
    """Median per-pair decile-1 ratio. `weights[line]` = draw multiplicity."""
    vals = []
    for ls, x, auc in per_pair:
        if weights is not None:
            w = weights[[line_ix[m] for m in ls]]
            if w.sum() < MIN_LINES_PER_PAIR:
                continue
            x2 = np.repeat(x, w)
            a2 = np.repeat(auc, w)
        else:
            x2, a2 = x, auc
        n = len(x2)
        k = max(10, int(round(n * q)))
        if n - k < 10:
            continue
        cut = np.partition(a2, k - 1)[k - 1]
        pos = a2 <= cut
        base = pos.mean()
        if base <= 0:
            continue
        order = np.argsort(x2, kind="mergesort")
        d1 = order[:max(5, n // 10)]
        vals.append(pos[d1].mean() / base)
    return float(np.median(vals)) if vals else np.nan


print(f"  {'q (target prevalence)':<24} {'pairs':>6} {'decile-1 ratio':>16} "
      f"{'95% CI (LINE cluster boot)':>28} {'frac<1':>8}")
res["marginal"] = {}
res["bootstrap_unit"] = "cell line (cluster bootstrap), not pair"
NB = 600
for q in QS:
    c = f"d1_ratio_q{q}"
    if c not in M:
        continue
    v = M[c].dropna().values
    obs = d1_median(q)
    b = []
    for _ in range(NB):
        w = np.bincount(RNG.integers(0, len(all_lines), len(all_lines)),
                        minlength=len(all_lines))
        val = d1_median(q, weights=w)
        if np.isfinite(val):
            b.append(val)
    lo, hi = np.percentile(b, [2.5, 97.5])
    print(f"  q = {q:<20.2f} {len(v):>6} {obs:>16.4f} "
          f"{f'[{lo:.4f}, {hi:.4f}]':>28} {np.mean(v < 1):>8.2f}")
    res["marginal"][f"{q:.2f}"] = {
        "n_pairs": int(len(v)), "d1_ratio_median": float(obs),
        "ci": [float(lo), float(hi)], "frac_below_1": float(np.mean(v < 1)),
        "n_boot": len(b)}
print()
print("  A working abundance gate predicts decile-1 ratio << 1 (low expression")
print("  depleted of sensitivity). An interval spanning 1.0 is a flat marginal.")

# ---------------------------------------------------------------- 3. three-way
print("\n3. THREE-WAY COVERAGE CAPACITY")
print("   For each usable (A, drug) pair, how many candidate modifier genes B")
print("   could support an interaction test at this coverage? This is CAPACITY,")
print("   not a selection -- B1 requires the A-B list to be a pre-registered")
print("   user input, and choosing B by outcome would invalidate the test.")
print()

alt_counts = {g: len(s) for g, s in alt.items()}
cap = []
for _, r in M.iterrows():
    g, did = r.target_ensg, r.drug_id
    sub = gd_cov[(gd_cov.target_ensg == g) & (gd_cov.drug_id == did)]
    lines = sorted(set(sub.model_id) & expr_lines & sequenced)
    n3 = len(lines)
    lineset = set(lines)
    n_B = sum(1 for gb, s in alt.items()
              if len(s & lineset) >= 10 and len(lineset - s) >= 10)
    cap.append({"symbol": r.symbol, "drug": r.drug, "n_lines_3way": n3,
                "n_candidate_B": n_B})
C = pd.DataFrame(cap)
usable3 = C[C.n_lines_3way >= MIN_LINES_PER_PAIR]
print(f"  pairs with >= {MIN_LINES_PER_PAIR} lines having expression + alteration "
      f"+ drug: {len(usable3):,}")
print(f"  median lines per pair after three-way intersection: "
      f"{usable3.n_lines_3way.median():.0f}")
print(f"  median candidate modifier genes per pair (>=10 altered and >=10 "
      f"unaltered): {usable3.n_candidate_B.median():.0f}")
res["three_way"] = {
    "usable_pairs": int(len(usable3)),
    "median_lines": float(usable3.n_lines_3way.median()) if len(usable3) else 0,
    "median_candidate_B": float(usable3.n_candidate_B.median()) if len(usable3) else 0,
}

M.to_parquet(OUT / "combined_score_b0_pairs.parquet", index=False)

# ---------------------------------------------------------------- verdict
print("\n" + "=" * 78)
print("VERDICT")
print("=" * 78)
n_usable = len(usable3)
print(f"  usable pairs after three-way coverage : {n_usable}")
print(f"  kill-switch threshold                 : {MIN_PAIRS}")
flat = all(v["ci"][0] <= 1.0 <= v["ci"][1] for v in res["marginal"].values())
print(f"  marginal on this pair set             : "
      f"{'FLAT (all intervals span 1.0)' if flat else 'NOT flat'}")
print()
if n_usable < MIN_PAIRS:
    verdict = "KILL_SWITCH_FIRED__insufficient_pairs"
    print("  KILL SWITCH FIRED. Too few usable pairs to validate an interaction.")
    print("  Stop. Ship Task A.")
else:
    verdict = "PAIR_COUNT_SUFFICIENT__blocked_on_preregistered_modifier_list"
    print("  Pair count clears the kill switch.")
    print("  BUT Stage 1 cannot start: B1 requires the A-B modifier list as a")
    print("  pre-registered user input with citations. It has not been supplied,")
    print("  and deriving it here would be selection by outcome.")
    print("  -> B0 continue, B1 BLOCKED pending that list.")
res["verdict"] = verdict
res["marginal_flat"] = bool(flat)

p = OUT / "combined_score_b0_results.json"
p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {p}")
