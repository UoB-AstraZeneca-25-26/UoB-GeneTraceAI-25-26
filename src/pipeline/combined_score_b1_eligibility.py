"""
combined_score_b1_eligibility.py
--------------------------------
TASK B1 -- eligibility of the pre-registered A-B list.

Reads the list from docs/COMBINED_SCORE_PREREGISTRATION.md (committed BEFORE this
runs) and applies the five inclusion criteria fixed there.

USES NO DRUG-RESPONSE DATA except to establish that an (A, drug) pair exists and
has enough lines. No AUC value is read, no sensitivity is computed. Eligibility
cannot therefore be contaminated by outcome.

KILL SWITCH: fewer than 10 eligible ESTABLISHED pairs -> stop, ship Task A.

    python src/pipeline/combined_score_b1_eligibility.py
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
MIN_LINES = 100
MIN_ALT = 10
KILL = 10

# Verbatim from docs/COMBINED_SCORE_PREREGISTRATION.md section 1.
ESTABLISHED = [
    ("EGFR", "KRAS"), ("BRAF", "MAP2K1"), ("BRCA1", "BRCA2"),
    ("TP53", "MDM2"), ("PTEN", "PIK3CA"), ("ERBB2", "ERBB3"),
    ("ALK", "EML4"), ("MYC", "MAX"), ("RB1", "CDK6"), ("VHL", "HIF1A"),
]
CONTRAST = [
    ("GAPDH", "KRAS"), ("ACTB", "EGFR"), ("BRCA1", "ALK"), ("TP53", "EML4"),
    ("PTEN", "MYC"), ("VHL", "BRAF"), ("RB1", "VHL"), ("ERBB2", "TP53"),
    ("MAX", "PTEN"), ("CDK6", "VHL"),
]

print("=" * 78)
print("TASK B1 -- ELIGIBILITY OF THE PRE-REGISTERED A-B LIST")
print("=" * 78)
print("  no drug-response VALUE is read here; only whether a pair exists")

E = K.expression()
fe = K.gene_validity()
sequenced, alt, _mod = K.alterations()

gd = pd.read_parquet(GDSC, columns=["model_id", "drug_id", "drug_name",
                                    "target_ensg", "target_symbol"])
gd["model_id"] = gd.model_id.astype(str).str.lower()
gd["target_ensg"] = (gd.target_ensg.astype("string").str.split(".").str[0]
                     .str.lower())
expr_lines = set(E.index)
gd = gd[gd.model_id.isin(expr_lines) & gd.target_ensg.isin(E.columns)]

alt_counts = {g: len(s) for g, s in alt.items()}
varying = {g for g, n in alt_counts.items() if n >= MIN_ALT}
print(f"\n  varying modifier genes (>= {MIN_ALT} altered lines): {len(varying):,}")


def check(a_sym, b_sym):
    r = {"A": a_sym, "B": b_sym}
    a_ensg, _ = K.resolve_gene(a_sym)
    b_ensg, _ = K.resolve_gene(b_sym)
    r["A_ensg"], r["B_ensg"] = a_ensg, b_ensg
    if a_ensg is None or b_ensg is None:
        r["eligible"] = False
        r["fail"] = "gene symbol did not resolve"
        return r
    # 1. A is a GDSC2 target
    sub_all = gd[gd.target_ensg == a_ensg]
    r["A_is_gdsc_target"] = bool(len(sub_all))
    if not len(sub_all):
        r["eligible"] = False
        r["fail"] = "1: A is not a GDSC2 drug target"
        return r
    # 2. A passes the validity guard
    r["A_frac_expressed"] = float(fe.get(a_ensg, 0.0))
    if r["A_frac_expressed"] < K.SILENT_FRAC:
        r["eligible"] = False
        r["fail"] = f"2: A fails validity guard ({r['A_frac_expressed']:.1%})"
        return r
    # 3. B is in the varying set
    r["B_altered_lines_panel"] = int(alt_counts.get(b_ensg, 0))
    if b_ensg not in varying:
        r["eligible"] = False
        r["fail"] = (f"3: B has only {r['B_altered_lines_panel']} altered lines "
                     f"(< {MIN_ALT})")
        return r
    # 4 + 5, per (A, drug) pair
    b_alt = alt.get(b_ensg, set())
    best = None
    for did, s in sub_all.groupby("drug_id"):
        lines = sorted(set(s.model_id) & expr_lines & sequenced)
        if len(lines) < MIN_LINES:
            continue
        ls = set(lines)
        n_alt = len(ls & b_alt)
        n_un = len(ls - b_alt)
        if n_alt < MIN_ALT or n_un < MIN_ALT:
            continue
        cand = {"drug_id": int(did),
                "drug": s.drug_name.iloc[0],
                "n_lines": len(lines), "n_B_altered": n_alt,
                "n_B_unaltered": n_un}
        if best is None or cand["n_lines"] > best["n_lines"]:
            best = cand
    n_drugs_ok = 0
    for did, s in sub_all.groupby("drug_id"):
        lines = sorted(set(s.model_id) & expr_lines & sequenced)
        ls = set(lines)
        if (len(lines) >= MIN_LINES and len(ls & b_alt) >= MIN_ALT
                and len(ls - b_alt) >= MIN_ALT):
            n_drugs_ok += 1
    r["n_eligible_drugs"] = n_drugs_ok
    if best is None:
        r["eligible"] = False
        r["fail"] = "4/5: no (A, drug) pair with >=100 lines and >=10 in each B arm"
        return r
    r.update(best)
    r["eligible"] = True
    r["fail"] = None
    return r


rows = []
for grp, pairs in (("established", ESTABLISHED), ("contrast", CONTRAST)):
    for a, b in pairs:
        r = check(a, b)
        r["set"] = grp
        rows.append(r)
R = pd.DataFrame(rows)

for grp in ("established", "contrast"):
    G = R[R.set == grp]
    print(f"\n{grp.upper()} SET")
    print(f"  {'A':<8} {'B':<8} {'elig':>5} {'drugs':>6} {'lines':>6} "
          f"{'B alt':>6} {'B unalt':>8}  detail")
    for _, r in G.iterrows():
        if r["eligible"]:
            print(f"  {r['A']:<8} {r['B']:<8} {'YES':>5} "
                  f"{r.get('n_eligible_drugs', 0):>6} {int(r['n_lines']):>6} "
                  f"{int(r['n_B_altered']):>6} {int(r['n_B_unaltered']):>8}  "
                  f"{r.get('drug')}")
        else:
            print(f"  {r['A']:<8} {r['B']:<8} {'no':>5} {'':>6} {'':>6} "
                  f"{'':>6} {'':>8}  {r['fail']}")

n_est = int(((R.set == "established") & R.eligible).sum())
n_con = int(((R.set == "contrast") & R.eligible).sum())

print("\n" + "=" * 78)
print("VERDICT")
print("=" * 78)
print(f"  eligible ESTABLISHED pairs : {n_est} of {len(ESTABLISHED)}")
print(f"  eligible CONTRAST pairs    : {n_con} of {len(CONTRAST)}")
print(f"  kill switch                : < {KILL} established -> stop")
print()
fails = R[(R.set == "established") & (~R.eligible)]
if len(fails):
    print("  why established pairs dropped:")
    for _, r in fails.iterrows():
        print(f"    {r['A']:<8} - {r['B']:<8}  {r['fail']}")
    print()
if n_est < KILL:
    verdict = "KILL_SWITCH_FIRED__too_few_eligible_established_pairs"
    print(f"  KILL SWITCH FIRED. {n_est} eligible established pairs, need {KILL}.")
    print("  An interaction cannot be validated on this many. Stop, report, and")
    print("  ship Task A.")
else:
    verdict = "ELIGIBLE__proceed_to_stage_1"
    print("  Proceed to Stage 1 under the pre-registered adoption rule.")
print(f"\n  VERDICT: {verdict}")

R.to_parquet(OUT / "combined_score_b1_eligibility.parquet", index=False)
p = OUT / "combined_score_b1_eligibility.json"
p.write_text(json.dumps({"n_established_eligible": n_est,
                         "n_contrast_eligible": n_con,
                         "kill_switch": KILL, "verdict": verdict,
                         "rows": rows}, indent=2, default=str), encoding="utf-8")
print(f"\nwrote {p}")
