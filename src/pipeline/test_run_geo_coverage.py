"""
test_run_geo_coverage.py
------------------------
Read-only diagnostic. Writes no production file.

QUESTION
--------
GEO could play three roles in this project:

  Role 1  a second opinion on DepMap's expression values for lines both cover
  Role 2  an independent expression source to re-test the gate on
          (test_run_geo_gate_replication.py)
  Role 3  EXTRA CELL LINES -- coverage DepMap does not have

This script settles Role 3, and only Role 3, by counting. If nearly every GEO
cell line is already in the DepMap expression cohort, GEO adds no coverage and
Role 3 closes itself; the discussion section then only has to defend Roles 1
and 2. If a meaningful number of lines are GEO-only, that is a genuine
extension of the addressable cell-line space and has to be reported as such.

The count that matters is not "GEO lines" but GEO lines that are RESOLVABLE --
mapped to a stable identity the rest of the pipeline can join on. A GSM whose
cell line cannot be resolved to a Cellosaurus accession adds nothing usable, so
the census walks the resolution chain explicitly and reports the attrition at
each step rather than a single headline number.

RESOLUTION CHAIN
----------------
    GSM sample --> cellosaurus_id (CVCL)  [geo_info_clean, already harmonised]
               --> DepMap model_id        [sample_info_clean.rrid, which IS CVCL]

A GEO line is "DepMap-known" if its CVCL resolves to a DepMap model at all, and
separately "DepMap-expressed" if that model actually has an RNA profile. These
differ, and the second is the one Role 3 turns on: a line DepMap has metadata
for but no expression for is still a coverage gain on the expression axis.

Lines whose CVCL resolves to no DepMap model are reported as GEO-only. That is
a lower bound on novelty, not an upper one -- some are likely non-cancer or
non-human samples rather than useful new models -- so the census also breaks
GEO-only lines down by whether Cellosaurus knows them at all.

Outputs:
    src/pipeline/outputs/test_run_geo_coverage_results.json
    src/pipeline/outputs/test_run_geo_coverage_per_series.parquet

Run:
    python src/pipeline/test_run_geo_coverage.py
"""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

warnings.filterwarnings("ignore", category=RuntimeWarning)

OUTPUTS = Path("src/pipeline/outputs")
REF_DIR = Path("reference")
DATA_CLEAN = Path("data/parquet/data_clean")
CHR_PATH = Path("data/DepMap_Chronos/GeneFitnessEffect_Chronos_Achilles.hdf5")

results = {}

# ============================================================ Step 1
print("=" * 78)
print("STEP 1 -- the GEO sample table and its resolution chain")
print("=" * 78)

info = pd.read_parquet(DATA_CLEAN / "geo_info_clean.parquet")
info["geo_accession"] = info.geo_accession.astype(str).str.lower()
info["cvcl"] = info.cellosaurus_id.astype(str).str.lower()
info.loc[~info.cvcl.str.startswith("cvcl_"), "cvcl"] = np.nan
info["gse_id"] = info.gse_id.astype("string").str.lower()

# Only GSMs that actually have a column in the expression matrix are usable.
expr_cols = set(pq.ParquetFile(DATA_CLEAN / "geo_expr_clean.parquet")
                .schema_arrow.names) - {"gene"}
info["has_expr"] = info.geo_accession.isin(expr_cols)

si = pd.read_parquet(DATA_CLEAN / "sample_info_clean.parquet",
                     columns=["depmap_id", "rrid", "stripped_cell_line_name",
                              "lineage", "primary_disease"])
si["rrid"] = si.rrid.astype(str).str.lower()
si["depmap_id"] = si.depmap_id.astype(str).str.lower()
cvcl2model = dict(zip(si.rrid, si.depmap_id))
info["model_id"] = info.cvcl.map(cvcl2model)

prof = pd.read_parquet(REF_DIR / "depmap_profiles.parquet")
rna_lines = set(prof[prof.datatype == "rna"].modelid.astype(str).str.lower())

# The DepMap expression cohort as the pipeline actually sees it: models with an
# RNA profile that survives into the cleaned expression matrix.
_ep = pq.read_table(DATA_CLEAN / "depmap_expr_clean.parquet",
                    columns=["index"]).to_pandas()
if "index" not in _ep.columns:          # arrow restores it as the frame index
    _ep = _ep.reset_index()
expr_prof = set(_ep["index"].astype(str))
rna_prof = prof[prof.datatype == "rna"].copy()
rna_prof["profileid"] = rna_prof.profileid.astype(str)
depmap_expressed = set(rna_prof[rna_prof.profileid.isin(expr_prof)]
                       .modelid.astype(str).str.lower())

print(f"  GEO GSM records                          : {len(info):,}")
print(f"  ... with a column in geo_expr_clean      : {int(info.has_expr.sum()):,}")
print(f"  ... resolved to a Cellosaurus accession  : "
      f"{int(info.cvcl.notna().sum()):,}")
print(f"  ... resolved onward to a DepMap model    : "
      f"{int(info.model_id.notna().sum()):,}")
print()
print(f"  DepMap RNA profiles (any)                : {len(rna_lines):,} models")
print(f"  DepMap expression cohort (in the matrix) : {len(depmap_expressed):,} models")

usable = info[info.has_expr].copy()
geo_lines_cvcl = set(usable.cvcl.dropna())
geo_lines_model = set(usable.model_id.dropna())
print()
print(f"  distinct GEO cell lines (by CVCL, expression-backed) : "
      f"{len(geo_lines_cvcl):,}")
print(f"  ... resolving to a DepMap model                      : "
      f"{len(geo_lines_model):,}")

results["chain"] = {
    "gsm_records": int(len(info)),
    "gsm_with_expression": int(info.has_expr.sum()),
    "gsm_with_cvcl": int(info.cvcl.notna().sum()),
    "gsm_with_depmap_model": int(info.model_id.notna().sum()),
    "depmap_rna_profiles_models": len(rna_lines),
    "depmap_expression_cohort": len(depmap_expressed),
    "geo_distinct_lines_cvcl": len(geo_lines_cvcl),
    "geo_distinct_lines_depmap_resolved": len(geo_lines_model),
}

# ============================================================ Step 2
print()
print("=" * 78)
print("STEP 2 -- WHAT DOES GEO ADD?  lines GEO has that DepMap expression does not")
print("=" * 78)

# A GEO line is novel on the expression axis when it is not in the DepMap
# expression cohort. Two flavours, because they answer different questions.
resolved = usable.dropna(subset=["cvcl"]).copy()
resolved["in_depmap_expr"] = resolved.model_id.isin(depmap_expressed)
resolved["in_depmap_any"] = resolved.model_id.notna()

by_line = (resolved.groupby("cvcl")
           .agg(model_id=("model_id", "first"),
                in_depmap_expr=("in_depmap_expr", "any"),
                in_depmap_any=("in_depmap_any", "any"),
                n_gsm=("geo_accession", "size"),
                name=("cellline", "first")))

n_total = len(by_line)
n_overlap = int(by_line.in_depmap_expr.sum())
n_new_expr = int((~by_line.in_depmap_expr).sum())
n_known_no_expr = int((by_line.in_depmap_any & ~by_line.in_depmap_expr).sum())
n_unknown = int((~by_line.in_depmap_any).sum())

print(f"  GEO expression-backed cell lines (CVCL)      : {n_total:,}")
print(f"    already in the DepMap expression cohort    : {n_overlap:,}  "
      f"({100*n_overlap/n_total:.1f}%)")
print(f"    NOT in the DepMap expression cohort        : {n_new_expr:,}  "
      f"({100*n_new_expr/n_total:.1f}%)")
print(f"      of which DepMap knows the model but has")
print(f"      no expression for it                     : {n_known_no_expr:,}")
print(f"      of which DepMap has no model at all      : {n_unknown:,}")
print()
print(f"  -> GEO would grow the expression-covered cell-line space from")
print(f"     {len(depmap_expressed):,} to {len(depmap_expressed) + n_new_expr:,} "
      f"(+{100*n_new_expr/len(depmap_expressed):.1f}%)")

results["coverage_gain"] = {
    "geo_lines_total": n_total,
    "overlap_with_depmap_expression": n_overlap,
    "geo_only_vs_depmap_expression": n_new_expr,
    "geo_only_depmap_model_exists": n_known_no_expr,
    "geo_only_no_depmap_model": n_unknown,
    "overlap_fraction": float(n_overlap / n_total),
    "depmap_expression_cohort": len(depmap_expressed),
    "combined_cohort": len(depmap_expressed) + n_new_expr,
    "pct_growth": float(100 * n_new_expr / len(depmap_expressed)),
}

# ============================================================ Step 3
print()
print("=" * 78)
print("STEP 3 -- PER SERIES.  n samples, n lines, and what each series adds")
print("=" * 78)
print("  'new' = line not in the DepMap expression cohort")
print("  'chronos' = line has a Broad Chronos label, i.e. usable for the gate test")

chronos_lines = set()
if CHR_PATH.exists():
    import h5py
    with h5py.File(CHR_PATH, "r") as fh:
        chronos_lines = {x.decode().lower() for x in fh["dim_0"][:]}
    print(f"  (Chronos Achilles: {len(chronos_lines):,} lines)")
else:
    print(f"  (WARNING: {CHR_PATH} absent -- chronos columns will read 0)")

resolved["in_chronos"] = resolved.model_id.isin(chronos_lines)
resolved["series"] = resolved.gse_id.fillna("(no gse_id)")

rows = []
for gse, sub in resolved.groupby("series"):
    lines = sub.groupby("cvcl").agg(
        in_expr=("in_depmap_expr", "any"), in_chr=("in_chronos", "any"))
    rows.append({
        "gse_id": gse,
        "n_gsm": int(len(sub)),
        "n_lines": int(len(lines)),
        "n_lines_in_depmap_expr": int(lines.in_expr.sum()),
        "n_lines_new": int((~lines.in_expr).sum()),
        "pct_new": float(100 * (~lines.in_expr).mean()),
        "n_lines_chronos": int(lines.in_chr.sum()),
        "n_lines_chronos_and_depmap_expr": int((lines.in_chr & lines.in_expr).sum()),
    })
S = pd.DataFrame(rows).sort_values("n_lines", ascending=False).reset_index(drop=True)
S.to_parquet(OUTPUTS / "test_run_geo_coverage_per_series.parquet")

print()
print(f"  {'series':<13} {'n_gsm':>6} {'lines':>6} {'in DepMap':>10} {'new':>6} "
      f"{'%new':>6} {'chronos':>8} {'gate-usable':>12}")
for _, r in S.iterrows():
    print(f"  {r.gse_id:<13} {r.n_gsm:>6} {r.n_lines:>6} "
          f"{r.n_lines_in_depmap_expr:>10} {r.n_lines_new:>6} "
          f"{r.pct_new:>6.1f} {r.n_lines_chronos:>8} "
          f"{r.n_lines_chronos_and_depmap_expr:>12}")

results["per_series"] = S.to_dict(orient="records")

# The union across series is not the sum -- lines recur. Report it explicitly so
# nobody adds the 'new' column up.
union_new = set(by_line[~by_line.in_depmap_expr].index)
print()
print(f"  NOTE: the 'new' column does not sum -- {S.n_lines_new.sum()} across")
print(f"  series but {len(union_new)} distinct lines, because series overlap.")
results["coverage_gain"]["new_lines_sum_over_series"] = int(S.n_lines_new.sum())
results["coverage_gain"]["new_lines_distinct"] = int(len(union_new))

# ============================================================ Step 4
print()
print("=" * 78)
print("STEP 4 -- ARE THE NEW LINES USABLE?")
print("=" * 78)
print("  A new line only helps if the pipeline can do something with it. The")
print("  binding constraint is a LABEL: without a Chronos or GDSC label the line")
print("  can be scored but never validated.")

new_lines = by_line[~by_line.in_depmap_expr]
new_with_model = new_lines[new_lines.model_id.notna()]
n_new_chronos = int(new_with_model.model_id.isin(chronos_lines).sum())
print()
print(f"  GEO-only lines (vs DepMap expression)      : {len(new_lines):,}")
print(f"    with a DepMap model_id                   : {len(new_with_model):,}")
print(f"    with a Chronos label                     : {n_new_chronos:,}")
print(f"    with neither                             : "
      f"{len(new_lines) - len(new_with_model):,}")

# how many GSMs back each new line -- a line seen once on one platform is
# weaker evidence than one seen repeatedly
if len(new_lines):
    print()
    print(f"  GSMs per GEO-only line: median {new_lines.n_gsm.median():.0f}, "
          f"{(new_lines.n_gsm == 1).mean()*100:.0f}% seen exactly once")

results["new_line_usability"] = {
    "geo_only_lines": int(len(new_lines)),
    "with_depmap_model_id": int(len(new_with_model)),
    "with_chronos_label": n_new_chronos,
    "with_neither": int(len(new_lines) - len(new_with_model)),
}

# ============================================================ Step 5
print()
print("=" * 78)
print("STEP 5 -- VERDICT ON ROLE 3")
print("=" * 78)
frac_new = n_new_expr / n_total
if frac_new < 0.05:
    v = "role3_closed__overlap_is_near_total"
    msg = ("GEO is a second opinion, not extra coverage. Role 3 needs no "
           "further discussion beyond one sentence recording the overlap.")
elif n_new_chronos < 30:
    v = "role3_marginal__new_lines_exist_but_are_unlabelled"
    msg = ("GEO does add lines, but almost none carry a dependency label, so "
           "they can be scored and never checked. Role 3 is a coverage claim "
           "only, not a validation claim.")
else:
    v = "role3_live__genuine_labelled_coverage_gain"
    msg = ("GEO adds labelled lines DepMap expression does not cover. Role 3 "
           "is worth arguing and the gain should be quantified in the writeup.")
results["verdict"] = v
print(f"  {100*frac_new:.1f}% of GEO's expression-backed lines are outside the")
print(f"  DepMap expression cohort ({n_new_expr:,} lines), of which "
      f"{n_new_chronos:,} carry a Chronos label.")
print()
print(f"  VERDICT: {v}")
print(f"  {msg}")

out = OUTPUTS / "test_run_geo_coverage_results.json"
out.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
print(f"wrote {OUTPUTS / 'test_run_geo_coverage_per_series.parquet'}")
