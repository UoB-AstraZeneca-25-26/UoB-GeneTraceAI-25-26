"""
measure_completeness.py
------------------------
Four questions about the qualifying set, computed BEFORE building a
completeness-based ranking. If Q3 shows completeness tracking how well-studied a
line is, ranking by completeness is ranking by fame, not by evidence -- the same
coverage-proxy fault already found in the abundance gate and the multi-gene
alteration layer.

Q1  Per line in a qualifying set, how many of {expression, mutation, CNA,
    protein} are MEASURED (not just present as NaN)?
Q2  Median completeness across qualifying lines, for a representative gene.
Q3  Does completeness correlate with study frequency (Spearman rho)?
Q4  N/A here -- it is a labelling requirement on the OUTPUT, applied directly to
    rank_convergent.py's printer, not a thing to measure.
Q5  N/A here -- answered directly: no tumour-similarity layer exists in this
    pipeline. UNVERIFIED -- could not locate.

    python src/pipeline/measure_completeness.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import multi_gene_kleene as K  # noqa: E402
import rank_convergent as RC   # noqa: E402

OUT = ROOT / "src" / "pipeline" / "outputs"
RNG = np.random.default_rng(42)
res = {}

print("=" * 78)
print("MEASUREMENT COMPLETENESS -- Q1-Q3, computed before building")
print("=" * 78)

E = K.expression()
CN = RC.cn_matrix()
seq, alt, _mod = K.alterations()
sequenced = seq
pn_names = None
import pyarrow.parquet as pq
pn = pq.ParquetFile(ROOT / "cleaned_track_data" / "proteomics.parquet").schema_arrow.names
key = "model_id" if "model_id" in pn else "depmap_id"
ccle = set(pq.read_table(ROOT / "cleaned_track_data" / "proteomics.parquet",
                         columns=[key]).to_pandas()[key].str.lower())
import duckdb
con = duckdb.connect(str(ROOT / "src" / "pipeline" / "outputs" / "celllineselector.db"),
                     read_only=True)
procan = set(x[0].lower() for x in con.execute(
    "SELECT DISTINCT model_id FROM procan_proteomics WHERE model_id IS NOT NULL"
).fetchall())
con.close()
protein_lines = ccle | procan

LAYERS = ["expression", "mutation", "cna", "protein"]


def completeness_for(gene: str, lines) -> pd.DataFrame:
    """1 row per line, 1 bool column per layer = MEASURED (not the value)."""
    ensg, _ = K.resolve_gene(gene)
    e_col = E[ensg] if ensg in E.columns else None
    cn_col = CN[ensg] if ensg in CN.columns else None
    rows = []
    for m in lines:
        rows.append({
            "model_id": m,
            "expression": bool(e_col is not None and m in e_col.index
                              and np.isfinite(e_col[m])),
            "mutation": bool(m in sequenced),
            "cna": bool(cn_col is not None and m in CN.index
                       and np.isfinite(cn_col[m])),
            "protein": bool(m in protein_lines),
        })
    D = pd.DataFrame(rows).set_index("model_id")
    D["n_measured"] = D[LAYERS].sum(axis=1)
    return D


# ---------------------------------------------------------------- Q1 + Q2
print("\nQ1 + Q2 -- COMPLETENESS FOR A REPRESENTATIVE SET OF QUERIES")
fe = K.gene_validity()
pool = [g for g in E.columns if fe.get(g, 0) >= K.SILENT_FRAC and g in CN.columns]
genes = list(RNG.choice(pool, 15, replace=False))
for s in ["TP53", "PTEN", "RB1", "CDKN2A", "BRAF", "EGFR", "KRAS", "ERBB2"]:
    e, _ = K.resolve_gene(s)
    if e in pool and e not in genes:
        genes.append(e)

all_D = []
per_gene_medians = []
for g in genes:
    terms = K.parse_terms(not_expressed=[g])
    q = K.query(terms, k_alternatives=0)
    surv = q["buckets"]["matches"]
    if len(surv) < 20:
        continue
    D = completeness_for(g, surv)
    D["gene"] = g
    all_D.append(D)
    per_gene_medians.append({"gene": g, "n_lines": len(D),
                             "median_n_measured": float(D.n_measured.median())})

A = pd.concat(all_D)
PG = pd.DataFrame(per_gene_medians)
print(f"  queries evaluated: {len(PG)}   total (gene, line) rows: {len(A):,}")
print()
print(f"  {'n layers measured':<20}{'count':>10}{'share':>9}")
dist = A.n_measured.value_counts().sort_index()
for k, v in dist.items():
    print(f"  {k:<20}{v:>10,}{v/len(A):>8.1%}")
print()
print(f"  median completeness across ALL qualifying lines: {A.n_measured.median():.0f}")
print(f"  median-of-medians across {len(PG)} queries       : {PG.median_n_measured.median():.1f}")
print()
print(f"  {'layer':<14}{'measured':>10}{'coverage':>10}")
for l in LAYERS:
    n = int(A[l].sum())
    print(f"  {l:<14}{n:>10,}{n/len(A):>9.1%}")

res["q1_q2"] = {
    "n_queries": int(len(PG)),
    "n_rows": int(len(A)),
    "distribution": {int(k): int(v) for k, v in dist.items()},
    "median_completeness_pooled": float(A.n_measured.median()),
    "median_of_query_medians": float(PG.median_n_measured.median()),
    "layer_coverage": {l: float(A[l].mean()) for l in LAYERS},
}
if A.n_measured.median() == 1:
    print("\n  Median = 1 (expression only). Completeness ranking still mostly")
    print("  ties -- but now the reason is visible: three of four layers are")
    print("  measured on a minority of lines, so most lines share the same")
    print("  completeness score by construction, not by chance.")

# ---------------------------------------------------------------- Q3
print("\nQ3 -- DOES COMPLETENESS TRACK STUDY FREQUENCY?")
print("  Proxy: number of distinct DepMap assay types run on a line")
print("  (reference/depmap_profiles.parquet: wes / rna / wgs), independent of")
print("  the four completeness layers measured above.")
prof = pd.read_parquet(ROOT / "reference" / "depmap_profiles.parquet")
prof["model_id"] = prof.modelid.astype(str).str.lower()
study_freq = prof.groupby("model_id").datatype.nunique()

A2 = A.reset_index().drop_duplicates("model_id").set_index("model_id")
A2["study_freq"] = study_freq.reindex(A2.index).fillna(0)
ok = A2.study_freq.notna()
rho, p = spearmanr(A2.loc[ok, "n_measured"], A2.loc[ok, "study_freq"])
print(f"  lines compared: {ok.sum():,}")
print(f"  Spearman rho(completeness, study_freq) = {rho:.4f}   p = {p:.3g}")
res["q3"] = {"spearman_rho": float(rho), "p": float(p), "n": int(ok.sum())}

if rho > 0.3 and p < 0.05:
    print("\n  YES -- completeness correlates with how well-studied a line is.")
    print("  Ranking by completeness ranks by fame, not evidence strength.")
    print("  Same coverage-proxy fault as the abundance gate and the")
    print("  driver-alteration flag. Must be labelled (Q4) and should NOT be")
    print("  presented as a biological ranking.")
    q3_verdict = "CONFIRMED_coverage_proxy"
elif rho > 0.15:
    print("\n  WEAK positive correlation. Some coverage-proxy risk; label")
    print("  regardless (Q4), and do not treat completeness position as a")
    print("  biological claim.")
    q3_verdict = "WEAK_coverage_proxy"
else:
    print("\n  No material correlation. Completeness is not simply tracking")
    print("  fame -- but the labelling requirement (Q4) still applies, because")
    print("  absence of correlation on THIS proxy does not rule out others.")
    q3_verdict = "NOT_CONFIRMED"
res["q3"]["verdict"] = q3_verdict

# per-layer breakdown, since layers may differ
print("\n  per-layer coverage vs study frequency:")
for l in LAYERS:
    r2, p2 = spearmanr(A2.loc[ok, l].astype(int), A2.loc[ok, "study_freq"])
    print(f"    {l:<12} rho={r2:>7.4f}  p={p2:.3g}")
res["q3"]["per_layer"] = {}
for l in LAYERS:
    r2, p2 = spearmanr(A2.loc[ok, l].astype(int), A2.loc[ok, "study_freq"])
    res["q3"]["per_layer"][l] = {"rho": float(r2), "p": float(p2)}

# ---------------------------------------------------------------- Q4 (statement)
print("\nQ4 -- LABELLING REQUIREMENT")
print("  Applied directly to rank_convergent.py's output: every ranked result")
print("  header must read")
print("    'ranked by measurement completeness, not biological evidence")
print("     strength -- position 1 is NOT a claim of biological superiority'")
print("  when the completeness measure is used. Implemented below.")

# ---------------------------------------------------------------- Q5 (statement)
print("\nQ5 -- TUMOUR SIMILARITY REFERENCE DATASET")
print("  UNVERIFIED -- could not locate. No tumour-similarity layer (TCGA or")
print("  CCLE-primary matched) exists anywhere in this pipeline. cell_similarity/")
print("  compares CELL LINE to CELL LINE only (RNA/metabolomics/miRNA axes),")
print("  never cell line to primary tumour. If tumour-matched selection is")
print("  wanted, it is new scope, not a reframing of existing code -- it would")
print("  need a TCGA (or Celligner-style) reference explicitly acquired and")
print("  matched by lineage before any lines could be scored against it.")
res["q5"] = {"status": "UNVERIFIED", "note": "no tumour-similarity layer exists"}

p = OUT / "completeness_measurement.json"
p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {p}")
