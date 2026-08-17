"""
05_biological_context_audit.py

MEASUREMENT-ONLY. No transformations, no aggregation, no writes to feature tables.
Answers one question per source, before any expression aggregation is implemented:

    For each (model_id, gene) pair, how many observations exist, and if >1,
    what experimental context (sample_id / study / treatment / condition /
    timepoint / replicate) distinguishes them?

This is the audit that decides whether you have a real aggregation problem
(different biological conditions that must NOT be averaged) or a non-problem
(one baseline value per line, or simply many genes per line).

Point it at whatever DepMap / GEO / HPA tables you are actually about to
aggregate. It is deliberately schema-tolerant: it inventories columns,
classifies them, and picks the right audit automatically.

Usage:
    python 05_biological_context_audit.py depmap  /path/OmicsExpression...csv
    python 05_biological_context_audit.py geo      /path/geo_expr.parquet  --meta /path/geo_gsm_metadata.parquet
    python 05_biological_context_audit.py hpa      /path/hpa_rna.parquet
"""

import argparse
import sys
import pandas as pd

# --- column classification vocabulary -------------------------------------
ID_TOKENS      = ("model_id", "ach-", "depmap_id", "cell_line", "cellline", "cell.line")
GENE_TOKENS    = ("gene", "symbol", "hgnc", "ensg", "entrez")
VALUE_TOKENS   = ("expr", "tpm", "log2", "logp1", "intensity", "value", "count", "nx")
SAMPLE_TOKENS  = ("sample_id", "gsm", "gse", "srr", "run", "experiment", "study", "series", "sample")
CONTEXT_TOKENS = ("treatment", "treated", "drug", "compound", "dose",
                  "condition", "perturb", "stimul", "agent",
                  "time", "timepoint", "hour", "day", "hr",
                  "replicate", "rep", "batch", "passage")


def classify_columns(cols):
    buckets = {"id": [], "gene": [], "value": [], "sample": [], "context": [], "other": []}
    for c in cols:
        lc = str(c).lower()
        if any(t in lc for t in CONTEXT_TOKENS):   buckets["context"].append(c)
        elif any(t in lc for t in SAMPLE_TOKENS):  buckets["sample"].append(c)
        elif any(t in lc for t in ID_TOKENS):      buckets["id"].append(c)
        elif any(t in lc for t in GENE_TOKENS):    buckets["gene"].append(c)
        elif any(t in lc for t in VALUE_TOKENS):   buckets["value"].append(c)
        else:                                       buckets["other"].append(c)
    return buckets


def looks_like_wide_matrix(df):
    """Wide expression matrix = index/first col is line ids, remaining cols are genes.
    In this layout each (line, gene) cell is unique BY CONSTRUCTION -> no within-line
    replicate aggregation is possible or needed. This is the DepMap case."""
    first = df.columns[0]
    first_vals = df[first].astype(str).head(50)
    id_like_rows = first_vals.str.contains("ACH-", case=False, na=False).mean()
    # many numeric-ish columns beyond the first = genes across the top
    numeric_cols = df.select_dtypes("number").shape[1]
    return id_like_rows > 0.5 and numeric_cols > 50


def audit_wide(df, source):
    print(f"\n[{source}] Detected WIDE matrix (lines x genes).")
    print(f"  rows (lines): {df.shape[0]:,}   columns (genes+meta): {df.shape[1]:,}")
    print("  -> Exactly one value per (line, gene) by construction.")
    print("  -> 'Multiple rows for a line' in a melted view = multiple GENES, not")
    print("     repeated measurements of the same gene. NO within-line aggregation needed.")
    dup_lines = df.iloc[:, 0].duplicated().sum()
    print(f"  duplicate line ids in first column: {dup_lines} "
          f"({'clean' if dup_lines == 0 else 'INVESTIGATE — a line appears on >1 row'})")


def audit_long(df, source, meta=None):
    b = classify_columns(df.columns)
    print(f"\n[{source}] Detected LONG table.")
    print(f"  id cols:      {b['id']}")
    print(f"  gene cols:    {b['gene']}")
    print(f"  value cols:   {b['value']}")
    print(f"  sample cols:  {b['sample']}")
    print(f"  context cols: {b['context']}")

    if not b["id"] or not b["gene"]:
        print("  !! Could not identify id and gene columns automatically.")
        print("     Rename them or edit ID_TOKENS/GENE_TOKENS, then re-run.")
        return

    id_col, gene_col = b["id"][0], b["gene"][0]
    grp = df.groupby([id_col, gene_col], observed=True).size()
    dup = grp[grp > 1]

    print(f"\n  (model_id, gene) pairs total: {len(grp):,}")
    print(f"  pairs with >1 observation:    {len(dup):,} "
          f"({100*len(dup)/max(len(grp),1):.1f}%)")

    if dup.empty:
        print("  -> No within-source duplication. One value per (line, gene).")
        print("     NO aggregation problem on this source.")
        return

    print(f"  max observations for one pair: {dup.max()}")
    print("  -> Duplicates exist. Now checking WHAT distinguishes them...")

    # Take a sample of duplicated pairs and report which candidate columns vary.
    varying = {c: 0 for c in b["sample"] + b["context"]}
    sample_pairs = dup.sample(min(500, len(dup)), random_state=0).index
    key = df.set_index([id_col, gene_col])
    for pair in sample_pairs:
        block = key.loc[[pair]]
        for c in varying:
            if c in block.columns and block[c].nunique(dropna=False) > 1:
                varying[c] += 1

    n = len(sample_pairs)
    print(f"\n  Of {n} sampled duplicated pairs, columns that VARY within the pair:")
    for c, k in sorted(varying.items(), key=lambda x: -x[1]):
        pct = 100 * k / n
        tag = ""
        if any(t in c.lower() for t in CONTEXT_TOKENS):
            tag = "  <-- EXPERIMENTAL CONTEXT: do NOT blindly average across these"
        print(f"    {c:<28} varies in {pct:5.1f}% of pairs{tag}")

    ctx_varies = any(varying[c] > 0 for c in b["context"])
    if ctx_varies:
        print("\n  VERDICT: duplicates span different experimental contexts.")
        print("  -> Separate by context first; aggregate only comparable (e.g. baseline) samples.")
    elif b["sample"]:
        print("\n  VERDICT: duplicates differ only by sample id, no context columns vary.")
        print("  -> Looks like replicates of one condition. Safe to aggregate (retain variability).")
    else:
        print("\n  VERDICT: duplicates present but no sample/context columns to explain them.")
        print("  -> Trace sample ids back to source metadata before deciding.")


def audit_geo_metadata(meta):
    """For GEO specifically: does the GSM-level metadata encode perturbation at all,
    or only line identity? This is what tells you baseline-study vs multi-study grab-bag."""
    b = classify_columns(meta.columns)
    print("\n[GEO metadata] GSM-level context scan.")
    print(f"  sample cols:  {b['sample']}")
    print(f"  context cols: {b['context']}")
    # Also scan free-text characteristics fields for perturbation tokens.
    text_cols = [c for c in meta.columns
                 if any(t in str(c).lower()
                        for t in ("characteristics", "source_name", "title", "description"))]
    print(f"  free-text cols scanned: {text_cols}")
    hits = {}
    for c in text_cols:
        joined = meta[c].astype(str).str.lower()
        for tok in ("treat", "drug", "dose", "µm", "um ", "vehicle", "control",
                    "24h", "48h", "hour", " day", "stimul", "sirna", "shrna",
                    "knockdown", "overexpress"):
            n = joined.str.contains(tok, na=False).sum()
            if n:
                hits[(c, tok)] = n
    if hits:
        print("  Perturbation-like tokens found in free-text metadata:")
        for (c, tok), n in sorted(hits.items(), key=lambda x: -x[1]):
            print(f"    {c}: '{tok}' in {n} samples")
        print("  -> Multi-context GEO source. Filter to baseline before per-line aggregation.")
    else:
        print("  No perturbation-like tokens found. Metadata appears to encode LINE IDENTITY only.")
        print("  -> Consistent with a single baseline cell-line deposit (e.g. GSE36133 / CCLE).")
        print("     One baseline profile per line; the treatment/timepoint audit does not apply.")


def load(path):
    if path.endswith((".parquet", ".pq")):
        return pd.read_parquet(path)
    return pd.read_csv(path, low_memory=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source", choices=["depmap", "geo", "hpa"])
    ap.add_argument("path")
    ap.add_argument("--meta", help="GEO GSM-level metadata table (optional but recommended for geo)")
    a = ap.parse_args()

    df = load(a.path)
    print("=" * 72)
    print(f"AUDIT: {a.source}  <-  {a.path}")
    print(f"shape: {df.shape}")
    print("=" * 72)

    if looks_like_wide_matrix(df):
        audit_wide(df, a.source)
    else:
        audit_long(df, a.source)

    if a.source == "geo" and a.meta:
        audit_geo_metadata(load(a.meta))

    print("\nDone. This module only measures; it changes nothing.")


if __name__ == "__main__":
    sys.exit(main())
