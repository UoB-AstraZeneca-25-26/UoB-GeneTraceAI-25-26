"""
================================================================================
WORKING MODEL — all four tracks feeding the gene x cell-line matrix
================================================================================
A presentable, runnable walkthrough of how every track's cleaned output becomes
the model's input. Synthetic stubs, but shaped with the REAL key columns and
grains from the team's pre-cleaning report. Swap stubs for the real parquet
spokes when they're keyed and the same code runs.

The point this makes for the team: the 14 datasets do NOT behave the same way.
There are four contribution patterns, and the matrix is assembled by handling
each correctly — not by merging everything into one table.

  PATTERN 1  per-pair evidence    gene + line -> one (gene, line) flag
             depmap_expr, mutations, fusions, proteomics
  PATTERN 2  cell-line evidence   line only   -> flag applies to ALL genes
             signatures, metabolomics                of that cell line
  PATTERN 3  metadata / context   not coverage -> held as annotation
             geo_info (tissue), hpa_desc (gene notes)
  PATTERN 4  own entity (Opt. B)  kept separate -> NOT folded into the matrix
             mirna                                  (its own node type later)

This is the DATA FLOW only. The scoring function is Phase 3 — not in here.
================================================================================
"""

import pandas as pd
pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 40)

# ===========================================================================
# 1. EACH TRACK'S CLEANED OUTPUT  (synthetic; real key columns + grain)
# ===========================================================================

# ---- Track A : identity spine + the profile->model bridge ----
sample_info = pd.DataFrame({
    "model_id":       ["ach-000001","ach-000002","ach-000003","ach-000004"],
    "cell_line_name": ["a549","hela","mcf7","k562"],
})
profiles = pd.DataFrame({                       # mutations are keyed on profile
    "profileid": ["pr-aaa","pr-bbb","pr-ccc","pr-ddd"],
    "model_id":  ["ach-000001","ach-000002","ach-000003","ach-000004"],
})

# ---- Track B : expression (gene x line) + gene annotation (gene only) ----
depmap_expr = pd.DataFrame({
    "ensg_id":  ["ensg001","ensg001","ensg002","ensg002","ensg003"],
    "model_id": ["ach-000001","ach-000002","ach-000001","ach-000003","ach-000002"],
    "tpm":      [42.0, 5.1, 88.0, 0.3, 12.0],
})
hpa_desc = pd.DataFrame({                        # gene-level only, no cell line
    "ensg_id":          ["ensg001","ensg002","ensg003"],
    "gene_description": ["kinase","tumour suppressor","membrane receptor"],
})

# ---- Track C : the four alteration datasets, each with its own quirk ----
mutations = pd.DataFrame({                       # per-variant (1,066,869 rows); real cols: ensemblgeneid, profileid, vepimpact, hotspot, af
    "ensemblgeneid": ["ensg001","ensg002","ensg002","ensg003"],
    "profileid":     ["pr-aaa","pr-aaa","pr-ccc","pr-bbb"],
    "vepimpact":     ["high","moderate","high","low"],
})
fusions = pd.DataFrame({                          # per gene-pair; real key: modelid, canonicalfusionname, ffpm
    "gene_a":   ["ensg001","ensg003"],
    "gene_b":   ["ensg009","ensg002"],
    "model_id": ["ach-000002","ach-000001"],      # real col name: modelid (no underscore)
})
signatures = pd.DataFrame({                       # per cell line only; real cols: modelid, msiscore, lohfraction, wgd, ploidy
    "model_id":      ["ach-000001","ach-000002"],  # real col name: modelid
    "sbs_signature": ["sbs1","sbs5"],
})
mirna = pd.DataFrame({                            # per miRNA x line, OWN ENTITY
    "mirbase_id": ["hsa-mir-21","hsa-mir-21","hsa-mir-155"],
    "model_id":   ["ach-000001","ach-000002","ach-000003"],
    "expression": [120.0, 45.0, 88.0],
})

# ---- Track D : protein + metabolic + metadata ----
proteomics = pd.DataFrame({                       # WIDE in reality (375 rows x 12,559 uniprot cols); stub is post-melt long form
    "ensg_id":   ["ensg001","ensg002"],            # real row key: depmap_id; uniprot->ensg is an extra hop
    "model_id":  ["ach-000001","ach-000003"],
    "abundance": [3.4, 1.1],
})
metabolomics = pd.DataFrame({                     # per cell line, MODEL ONLY
    "model_id":           ["ach-000001","ach-000004"],
    "metabolite_profile": ["m-high","m-low"],
})
geo_info = pd.DataFrame({                          # per cell line, METADATA only
    "model_id": ["ach-000002","ach-000003"],
    "tissue":   ["cervix","breast"],
})

print("=" * 78)
print("STEP 1 — each track's REAL cleaned table, at its own grain")
print("=" * 78)
# Numbers and column lists from the team's Pre-Cleaning Report (Table Report xlsx).
# LONG  = already one row per record; keys ready to attach.
# WIDE  = matrix (genes or cell lines as columns); must be melted to long form
#         before per-pair keys can be attached — an extra reshape step.
real = [
    ("A", "sample_info",      "1,840",       "29",     "LONG", "depmap_id (=model_id), cell_line_name, rrid, ccle_name"),
    ("A", "cellosaurus",      "152,231",     "17",     "LONG", "cellosaurus_accession (CVCL), synonyms, cross-references"),
    ("A", "depmap_profiles",  "3,830",       "5",      "LONG", "profileid -> modelid  (the bridge)"),
    ("B", "depmap_expr",      "1,495",       "53,961", "WIDE", "rows=model, cols=ensg ids; melt to (ensg_id, model_id, tpm)"),
    ("B", "hpa_rna",          "24,315,372",  "6",      "LONG", "gene, cell line, tpm/ptpm/ntpm"),
    ("B", "geo_expr",         "19,914",      "3,268",  "WIDE", "rows=gene, cols=gsm; melt to (gene, gsm, value)"),
    ("B", "hpa_desc",         "1,206",       "7",      "LONG", "cell line, disease, cellosaurus id  (annotation)"),
    ("C", "mutations",        "1,066,869",   "70",     "LONG", "ensemblgeneid, profileid, vepimpact, hotspot, af"),
    ("C", "fusions",          "184,237",     "30",     "LONG", "modelid, canonicalfusionname, ffpm, confidence"),
    ("C", "signatures",       "3,021",       "12",     "LONG", "modelid, msiscore, lohfraction, wgd, ploidy  (no gene)"),
    ("C", "mirna",            "734",         "956",    "WIDE", "rows=miRNA, cols=cell line; own entity (Option B)"),
    ("D", "proteomics",       "375",         "12,559", "WIDE", "rows=model (depmap_id), cols=uniprot; melt to long"),
    ("D", "metabolomics",     "928",         "227",    "WIDE", "rows=model, cols=metabolite; cell-line level (no gene)"),
    ("D", "geo_info",         "3,267",       "23",     "LONG", "geo_accession, characteristics_ch1  (metadata)"),
]
print(f"  {'trk':<4}{'table':<17}{'rows':>13}{'cols':>8}  {'form':<6} key columns / note")
print("  " + "-" * 74)
for trk, name, rows, cols, form, note in real:
    print(f"  {trk:<4}{name:<17}{rows:>13}{cols:>8}  {form:<6} {note}")
print()
print("  WIDE tables (depmap_expr, geo_expr, proteomics, metabolomics, mirna)")
print("  must be melted to long form before per-pair keys can be attached.")
print("  Same architecture as LONG tables — one extra reshape step.")
print()
print("  (Steps 2–5 below use small synthetic stubs to show the join mechanics")
print("   on real column names. The counts above are the verified real ones.)")
print()

# ===========================================================================
# 2. REDUCE EACH TO ITS CONTRIBUTION  (four patterns, handled differently)
# ===========================================================================

def gene_line(df, gcol, mcol):
    out = df[[gcol, mcol]].dropna().drop_duplicates()
    out.columns = ["ensg_id","model_id"]
    return out

# PATTERN 1 — per-pair evidence
mut_keyed = mutations.merge(profiles, on="profileid", how="left", validate="m:1")
fus_long = pd.concat([                            # explode pair -> two gene rows
    fusions.rename(columns={"gene_a":"ensg_id"})[["ensg_id","model_id"]],
    fusions.rename(columns={"gene_b":"ensg_id"})[["ensg_id","model_id"]],
]).drop_duplicates()

expr_p = gene_line(depmap_expr, "ensg_id", "model_id")
mut_p  = gene_line(mut_keyed,   "ensemblgeneid", "model_id")
prot_p = gene_line(proteomics,  "ensg_id", "model_id")

universe = pd.concat([expr_p, mut_p, fus_long, prot_p]).drop_duplicates()

matrix = (universe
    .merge(expr_p.assign(has_expr=1),         on=["ensg_id","model_id"], how="left")
    .merge(mut_p.assign(has_mutation=1),      on=["ensg_id","model_id"], how="left")
    .merge(fus_long.assign(has_fusion=1),     on=["ensg_id","model_id"], how="left")
    .merge(prot_p.assign(has_proteomics=1),   on=["ensg_id","model_id"], how="left")
    .fillna(0))

# PATTERN 2 — cell-line evidence: flag joins on model_id, applies to ALL genes
sig_flag = signatures[["model_id"]].drop_duplicates().assign(has_signature=1)
met_flag = metabolomics[["model_id"]].drop_duplicates().assign(has_metabolomics=1)
matrix = (matrix
    .merge(sig_flag, on="model_id", how="left")
    .merge(met_flag, on="model_id", how="left")
    .fillna(0))

pair_flags = ["has_expr","has_mutation","has_fusion","has_proteomics"]
line_flags = ["has_signature","has_metabolomics"]
for c in pair_flags + line_flags:
    matrix[c] = matrix[c].astype(int)
matrix["coverage"] = matrix[pair_flags + line_flags].sum(axis=1)

print("=" * 72)
print("STEP 2 — how each pattern contributes")
print("=" * 72)
print("  PATTERN 1  per-pair    expr, mutation, fusion, proteomics")
print("             -> a flag on the exact (gene, cell line) pair")
print("  PATTERN 2  cell-line   signature, metabolomics")
print("             -> a flag on the cell line, shared by all its genes")
print("  PATTERN 3  metadata    geo_info (tissue), hpa_desc (gene notes)")
print("             -> held as context, does NOT count as coverage")
print("  PATTERN 4  own entity  mirna  -> kept as its own table (Option B)")
print()

# ===========================================================================
# 3. THE MATRIX  (per-pair flags | cell-line flags | coverage)
# ===========================================================================
ordered = ["ensg_id","model_id"] + pair_flags + line_flags + ["coverage"]
print("=" * 72)
print("STEP 3 — the model's input, derived from all four tracks")
print("=" * 72)
print(matrix[ordered].sort_values(["ensg_id","model_id"]).to_string(index=False))
print()
print("  left block  = per-pair evidence (gene-specific)")
print("  right block = cell-line evidence (same for every gene of that line)")
print("  NOTE: whether cell-line flags should count toward coverage the same as")
print("  per-pair flags is a DECISION for the model owner — flagged, not assumed.")
print()

# ===========================================================================
# 4. HELD SEPARATELY  (not folded into the gene x line matrix)
# ===========================================================================
print("=" * 72)
print("STEP 4 — what is deliberately NOT in the matrix")
print("=" * 72)
print("mirna — kept as its OWN entity (Option B), for the knowledge-graph phase:")
print(mirna.to_string(index=False))
print()
print("hpa_desc — gene-level annotation, joins by gene as enrichment, not coverage:")
print(hpa_desc.to_string(index=False))
print()
print("geo_info — cell-line metadata (tissue), context for results, not evidence:")
print(geo_info.to_string(index=False))
print()

# ===========================================================================
# 5. REVERSE-ENGINEERING POINT
# ===========================================================================
print("=" * 72)
print("STEP 5 — the matrix is DERIVED, so it is re-derivable")
print("=" * 72)
print("Every table above stayed independent at its own grain. The matrix was")
print("built by joins, not by merging. If the model changes what it needs, we")
print("re-run this script with different columns — no track is re-cleaned.")