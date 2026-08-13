"""
test_kleene_ground_truth.py
---------------------------
Single-gene ground-truth test of the Kleene layer against named cell lines with
known molecular status.

This is the "boring test with a known answer" applied to `multi_gene_kleene.py`.
Every row is a published, textbook fact about a specific line. The layer either
reproduces it or it does not, and a disagreement is either a data-coverage
finding or a bug.

EXPECTATIONS BY ANNOTATION TYPE -- and what the layer actually claims
---------------------------------------------------------------------
    annotation          expected                         in scope?
    wild-type           altered(G) == FALSE              yes
    mutation            altered(G) == TRUE               yes
    fusion              altered(G) == TRUE               yes  (in-frame fusions)
    overexpressed       expressed(G) == TRUE             yes
    amplification       expressed(G) == TRUE             yes for expression
                        altered(G) == TRUE               NO -- see below
    deletion / null     expressed(G) == FALSE            PARTIAL -- see below

TWO KNOWN SCOPE LIMITS, DECLARED BEFORE MEASURING
--------------------------------------------------
1. `altered()` is built from `any_driver` (point mutations, indels) plus
   in-frame fusions. It does NOT read copy number. So an AMPLIFIED line is not
   expected to return altered == TRUE, and a homozygous DELETION is not expected
   to return altered == TRUE either. Those rows test the EXPRESSION term instead,
   and the altered() result is recorded as out-of-scope rather than as a failure.

2. "null" is usually a PROTEIN-level statement. A line can be protein-null and
   still transcribe the gene -- PTEN-null lines frequently retain mRNA, and a
   nonsense-mutant transcript may escape decay. So `expressed(G) == FALSE` is the
   expectation for null lines but a disagreement there is evidence about the
   annotation's level, not necessarily about the layer.

Both limits are stated here so the pass rate is read against what the layer
claims, not against what one might assume it claims.

    python src/pipeline/test_kleene_ground_truth.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import multi_gene_kleene as K  # noqa: E402

OUT = ROOT / "src" / "pipeline" / "outputs"

# (gene, line name, annotation, category)
GT = [
    ("TP53", "MCF7", "wild-type", "wildtype"),
    ("TP53", "HCT116", "wild-type", "wildtype"),
    ("TP53", "SAOS2", "null", "deletion_null"),
    ("TP53", "H1299", "null", "deletion_null"),
    ("EGFR", "A431", "overexpressed", "overexpression"),
    ("EGFR", "PC9", "exon 19 deletion", "mutation"),
    ("EGFR", "HCC827", "exon 19 deletion", "mutation"),
    ("EGFR", "H1975", "T790M", "mutation"),
    ("KRAS", "A549", "G12S", "mutation"),
    ("KRAS", "HCT116", "G13D", "mutation"),
    ("KRAS", "SW480", "G12V", "mutation"),
    ("KRAS", "PANC1", "G12D", "mutation"),
    ("BRAF", "A375", "V600E", "mutation"),
    ("BRAF", "SKMEL28", "V600E", "mutation"),
    ("BRAF", "COLO205", "V600E", "mutation"),
    ("PTEN", "U87MG", "null/mutant", "deletion_null"),
    ("PTEN", "LNCAP", "null", "deletion_null"),
    ("PTEN", "MDAMB468", "null", "deletion_null"),
    ("ERBB2", "SKBR3", "amplified", "amplification"),
    ("ERBB2", "BT474", "amplified", "amplification"),
    ("ERBB2", "AU565", "amplified", "amplification"),
    ("ERBB2", "NCIN87", "amplified", "amplification"),
    ("BRCA1", "HCC1937", "mutant", "mutation"),
    ("BRCA1", "MDAMB436", "mutant", "mutation"),
    ("BRCA2", "CAPAN1", "mutant", "mutation"),
    ("RB1", "SAOS2", "null", "deletion_null"),
    ("RB1", "WERIRB1", "null", "deletion_null"),
    ("MYC", "HL60", "amplified", "amplification"),
    ("MYC", "RAJI", "translocated", "fusion"),
    ("MYC", "P4936", "tet-inducible model", "mutation"),
    ("VHL", "786O", "mutant/null", "mutation"),
    ("VHL", "RCC4", "mutant", "mutation"),
    ("VHL", "A498", "mutant", "mutation"),
    ("ALK", "NCIH2228", "EML4-ALK fusion", "fusion"),
    ("ALK", "NCIH3122", "EML4-ALK fusion", "fusion"),
    ("ALK", "KARPAS299", "NPM-ALK fusion", "fusion"),
    ("PIK3CA", "MCF7", "E545K", "mutation"),
    ("PIK3CA", "T47D", "H1047R", "mutation"),
    ("CDKN2A", "U2OS", "deleted", "deletion_null"),
    ("CDKN2A", "MIAPACA2", "deleted", "deletion_null"),
    ("APC", "SW480", "mutant", "mutation"),
    ("APC", "CACO2", "mutant", "mutation"),
    ("NRAS", "SKMEL2", "Q61R", "mutation"),
    ("NRAS", "IPC298", "Q61L", "mutation"),
    ("MET", "EBC1", "amplified", "amplification"),
    ("MET", "HS746T", "amplified / exon 14 skipping", "amplification"),
    # CORRECTED 2026-08-10: supplied as "wild-type, hormone-sensitive".
    # LNCaP carries AR T877A (T878A) in the ligand-binding domain -- the
    # classic promiscuous-AR model. The layer flagged this before the
    # reference was corrected; see MULTIGENE_GROUND_TRUTH.md.
    ("AR", "LNCAP", "T878A mutant [corrected]", "mutation"),
    ("AR", "22RV1", "splice variant", "mutation"),
    ("AR", "VCAP", "amplified", "amplification"),
    ("STK11", "A549", "mutant", "mutation"),
    ("STK11", "NCIH460", "mutant", "mutation"),
    ("CTNNB1", "HEPG2", "mutant", "mutation"),
    ("CTNNB1", "DLD1", "mutant", "mutation"),
]

# What each category asserts, and about which term.
EXPECT = {
    "wildtype":       {"altered": "FALSE"},
    "mutation":       {"altered": "TRUE"},
    "fusion":         {"altered": "TRUE"},
    "overexpression": {"expressed": "TRUE"},
    "amplification":  {"expressed": "TRUE"},        # altered() has no CN input
    "deletion_null":  {"expressed": "FALSE"},       # protein-level caveat applies
}
OUT_OF_SCOPE_NOTE = {
    "amplification": "altered() reads mutations and fusions, not copy number",
    "deletion_null": "'null' is usually protein-level; mRNA may persist",
}


def norm(s):
    return "".join(ch for ch in str(s).upper() if ch.isalnum())


# DepMap names carry vendor prefixes and clone suffixes that a plain normalised
# match misses: LNCaP is "LNCaP clone FGC", H1299 is "NCI-H1299", 786-O is
# "786O_KIDNEY". These aliases are generated mechanically, not curated per line.
_PREFIXES = ("NCI", "ACC", "CCL")
_SUFFIXES = ("CLONEFGC", "CLONE", "CELLS")


def _aliases(v: str):
    """Every reasonable normalised form of one catalogue name."""
    out = set()
    n = norm(v)
    if not n:
        return out
    out.add(n)
    out.add(norm(v.split("_")[0]))            # strip CCLE tissue suffix
    for pre in _PREFIXES:                      # NCIH1299 -> H1299
        if n.startswith(pre) and len(n) > len(pre) + 1:
            out.add(n[len(pre):])
    for suf in _SUFFIXES:                      # LNCAPCLONEFGC -> LNCAP
        if n.endswith(suf) and len(n) > len(suf) + 1:
            out.add(n[: -len(suf)])
    return {a for a in out if len(a) >= 3}


def build_resolver():
    """cell line name -> model_id, from every name column available."""
    m = {}
    si = pd.read_parquet(ROOT / "cleaned_track_data" / "sample_info.parquet")
    cols = [c for c in si.columns
            if "name" in c.lower() and si[c].dtype == object]
    for _, r in si.iterrows():
        mid = str(r["model_id"]).lower()
        for c in cols:
            v = r.get(c)
            if isinstance(v, str) and v.strip():
                for a in _aliases(v):
                    m.setdefault(a, mid)
    cl = pd.read_parquet(ROOT / "reference" / "cell_line_lookup.parquet",
                         columns=["model_id", "cell_line_name", "canonical_name"])
    for _, r in cl.iterrows():
        mid = str(r["model_id"]).lower()
        for c in ("cell_line_name", "canonical_name"):
            v = r.get(c)
            if isinstance(v, str) and v.strip():
                for a in _aliases(v):
                    m.setdefault(a, mid)
    return m


print("=" * 78)
print("KLEENE LAYER -- SINGLE-GENE GROUND-TRUTH TEST")
print("=" * 78)

resolver = build_resolver()
E = K.expression()
sequenced, alt, _mod = K.alterations()
print(f"  name index: {len(resolver):,} aliases")
print(f"  ground-truth rows: {len(GT)}   genes: {len(set(g for g, *_ in GT))}")

rows = []
for gene, line, ann, cat in GT:
    mid = resolver.get(norm(line))
    r = {"gene": gene, "line": line, "annotation": ann, "category": cat,
         "model_id": mid}
    if mid is None:
        r["status"] = "LINE_NOT_FOUND"
        rows.append(r)
        continue
    ensg, sym = K.resolve_gene(gene)
    r["ensg"] = ensg
    if ensg is None:
        r["status"] = "GENE_NOT_FOUND"
        rows.append(r)
        continue
    te = K.Term("expressed", gene)
    ta = K.Term("altered", gene)
    st_e, rs_e, v_e = te.evaluate([mid])
    st_a, rs_a, _ = ta.evaluate([mid])
    r["expressed_state"] = st_e.iloc[0]
    r["expressed_reason"] = rs_e.iloc[0]
    r["log2tpm"] = float(v_e.iloc[0]) if pd.notna(v_e.iloc[0]) else None
    r["altered_state"] = st_a.iloc[0]
    r["altered_reason"] = rs_a.iloc[0]

    exp = EXPECT[cat]
    term = "altered" if "altered" in exp else "expressed"
    want = exp[term]
    got = r[f"{term}_state"]
    r["tested_term"] = term
    r["expected"] = want
    r["observed"] = got
    if got == "UNKNOWN":
        r["status"] = "UNKNOWN"
    elif got == want:
        r["status"] = "PASS"
    else:
        r["status"] = "FAIL"
    rows.append(r)

R = pd.DataFrame(rows)

# ---------------------------------------------------------------- sensitivity
# `any_driver` misses activating IN-FRAME INDELS. VEP scores an inframe deletion
# as MODERATE impact (max_vep_rank 2), and only 0.5% of rank-2 rows carry
# any_driver, against 100% of rank-3 (HIGH). EGFR exon-19 deletion and CTNNB1
# exon-3 deletion are exactly this class. Test a widened definition on the same
# ground truth rather than silently changing the layer.
mut = pd.read_parquet(ROOT / "cleaned_track_data" / "mutations_collapsed.parquet",
                      columns=["ensg_id", "model_id", "any_driver",
                               "oncogene_hit", "tsg_hit", "max_vep_rank"])
mut["model_id"] = mut.model_id.astype(str).str.lower()
mut["ensg_id"] = mut.ensg_id.astype(str).str.split(".").str[0].str.lower()
wide = mut[(mut.any_driver.fillna(False).astype(bool))
           | (mut.oncogene_hit.fillna(False).astype(bool))
           | (mut.tsg_hit.fillna(False).astype(bool))]
wide_set = wide.groupby("ensg_id")["model_id"].apply(set).to_dict()
# The widened rule must keep the FUSION contribution too, or the comparison is
# not like-for-like -- a first version omitted it and appeared to turn the three
# ALK fusion calls FALSE, which was an artefact of the check, not of the rule.
_fus = pd.read_parquet(ROOT / "cleaned_track_data" / "fusions_gene_level.parquet",
                       columns=["ensg_id", "model_id", "any_in_frame"])
_fus["model_id"] = _fus.model_id.astype(str).str.lower()
_fus["ensg_id"] = _fus.ensg_id.astype(str).str.split(".").str[0].str.lower()
for _g, _s in (_fus[_fus.any_in_frame.fillna(False).astype(bool)]
               .groupby("ensg_id")["model_id"].apply(set).items()):
    wide_set[_g] = wide_set.get(_g, set()) | _s

alt_rows = []
for _, r in R.iterrows():
    if r.get("tested_term") != "altered" or pd.isna(r.get("model_id")):
        continue
    got_wide = "TRUE" if r["model_id"] in wide_set.get(r["ensg"], set()) else "FALSE"
    if r["observed"] == "UNKNOWN":
        got_wide = "UNKNOWN"
    alt_rows.append({"gene": r["gene"], "line": r["line"],
                     "any_driver": r["observed"], "widened": got_wide,
                     "expected": r["expected"],
                     "changed": got_wide != r["observed"]})
A = pd.DataFrame(alt_rows)

R.to_parquet(OUT / "kleene_ground_truth.parquet", index=False)

# ------------------------------------------------------------ report
print("\n" + "-" * 78)
print(f"  {'gene':<8} {'line':<11} {'annotation':<22} {'term':<10} "
      f"{'want':<6} {'got':<8} {'log2TPM':>8}  status")
print("-" * 78)
for _, r in R.iterrows():
    if r["status"] in ("LINE_NOT_FOUND", "GENE_NOT_FOUND"):
        print(f"  {r['gene']:<8} {r['line']:<11} {r['annotation']:<22} "
              f"{'--':<10} {'--':<6} {'--':<8} {'--':>8}  {r['status']}")
        continue
    lt = f"{r['log2tpm']:.2f}" if r.get("log2tpm") is not None else "--"
    mark = {"PASS": "PASS", "FAIL": "**FAIL**", "UNKNOWN": "unknown"}[r["status"]]
    print(f"  {r['gene']:<8} {r['line']:<11} {str(r['annotation'])[:22]:<22} "
          f"{r['tested_term']:<10} {r['expected']:<6} {r['observed']:<8} "
          f"{lt:>8}  {mark}")

res = {}
print("\n" + "=" * 78)
print("SUMMARY")
print("=" * 78)
found = R[~R.status.isin(["LINE_NOT_FOUND", "GENE_NOT_FOUND"])]
print(f"  rows                    : {len(R)}")
print(f"  lines resolved to DepMap: {len(found)} "
      f"({int((R.status == 'LINE_NOT_FOUND').sum())} not found)")
if len(found):
    n_pass = int((found.status == "PASS").sum())
    n_fail = int((found.status == "FAIL").sum())
    n_unk = int((found.status == "UNKNOWN").sum())
    print(f"  PASS {n_pass}   FAIL {n_fail}   UNKNOWN {n_unk}")
    resolved = n_pass + n_fail
    if resolved:
        print(f"  agreement among RESOLVED calls: {n_pass}/{resolved} "
              f"({n_pass/resolved:.1%})")
    res["overall"] = {"rows": len(R), "resolved_lines": len(found),
                      "pass": n_pass, "fail": n_fail, "unknown": n_unk,
                      "agreement": n_pass / resolved if resolved else None}

print(f"\n  BY CATEGORY")
print(f"  {'category':<16} {'n':>4} {'pass':>5} {'fail':>5} {'unk':>5} "
      f"{'agreement':>10}  scope note")
res["by_category"] = {}
for cat, g in found.groupby("category"):
    p_ = int((g.status == "PASS").sum())
    f_ = int((g.status == "FAIL").sum())
    u_ = int((g.status == "UNKNOWN").sum())
    ag = p_ / (p_ + f_) if (p_ + f_) else float("nan")
    print(f"  {cat:<16} {len(g):>4} {p_:>5} {f_:>5} {u_:>5} {ag:>9.0%}  "
          f"{OUT_OF_SCOPE_NOTE.get(cat, '')}")
    res["by_category"][cat] = {"n": len(g), "pass": p_, "fail": f_,
                               "unknown": u_, "agreement": float(ag)}

bad = found[found.status == "FAIL"]
if len(bad):
    print(f"\n  DISAGREEMENTS ({len(bad)})")
    for _, r in bad.iterrows():
        lt = f"{r['log2tpm']:.2f}" if r.get("log2tpm") is not None else "n/a"
        print(f"    {r['gene']:<8} {r['line']:<11} {r['annotation']:<24} "
              f"expected {r['expected']}, got {r['observed']}  "
              f"(log2TPM {lt})")

unk = found[found.status == "UNKNOWN"]
if len(unk):
    print(f"\n  UNKNOWN ({len(unk)}) -- the layer declining to call")
    for _, r in unk.iterrows():
        # altered() now has TWO distinct UNKNOWN reasons -- never-sequenced, and
        # the Fix 1 MODERATE-only abstention. Hardcoding "line never sequenced"
        # attributed every altered abstention to the wrong cause and made Fix 1
        # look inert in the report.
        why = (r["expressed_reason"] if r["tested_term"] == "expressed"
               else r.get("altered_reason") or "line never sequenced")
        print(f"    {r['gene']:<8} {r['line']:<11} {r['tested_term']:<10} {why}")

missing = R[R.status == "LINE_NOT_FOUND"]
if len(missing):
    print(f"\n  LINES NOT IN DEPMAP ({len(missing)})")
    print("    " + ", ".join(sorted(set(missing.line))))

print("\n" + "=" * 78)
print("SENSITIVITY -- would a widened altered() definition help?")
print("=" * 78)
print("  widened = any_driver OR oncogene_hit OR tsg_hit")
if len(A):
    ch = A[A.changed]
    n_ok_now = int((A.any_driver == A.expected).sum())
    n_ok_wide = int((A.widened == A.expected).sum())
    print(f"  altered() rows tested       : {len(A)}")
    print(f"  correct with any_driver     : {n_ok_now}")
    print(f"  correct with widened rule   : {n_ok_wide}")
    print(f"  calls changed               : {len(ch)}")
    for _, r in ch.iterrows():
        print(f"    {r['gene']:<8} {r['line']:<11} {r['any_driver']} -> "
              f"{r['widened']}  (expected {r['expected']})")
    if n_ok_wide <= n_ok_now:
        print("\n  The widened rule does NOT recover the misses. The failures are")
        print("  not a flag-threshold problem -- EGFR exon-19 rows exist but are")
        print("  scored MODERATE impact, and the CTNNB1/MYC lesions have no row at")
        print("  all. altered() is left unchanged.")
    res["widened_definition"] = {"correct_any_driver": n_ok_now,
                                 "correct_widened": n_ok_wide,
                                 "n_changed": int(len(ch))}

p = OUT / "kleene_ground_truth_results.json"
p.write_text(json.dumps({**res, "rows": rows}, indent=2, default=str),
             encoding="utf-8")
print(f"\nwrote {p}")
