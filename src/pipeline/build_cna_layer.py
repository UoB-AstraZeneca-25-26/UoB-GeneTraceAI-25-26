"""
build_cna_layer.py
------------------
Builds CNA confidence modifier flags from COSMIC CellLines CNA data.
Outputs: src/pipeline/outputs/cna_flags.parquet

Role in pipeline:
    CNA is a CONFIDENCE MODIFIER, not a ranking layer.
    Same pattern as has_driver_alteration in flags_with_driver.parquet.
    Gene-class gating:
        oncogene / activation_driven  → amplification (CN > 2.5) raises confidence
        tsg / loss_of_function        → deletion (CN < 1.5) raises confidence
        unknown_dual                  → CNA surfaced as context only, no confidence change

[TRACK-C][CNA] Initial build
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
CNA_FILE   = "data/COSMIC/CellLinesProject_CompleteCNA_v104_GRCh37.tsv"
MODEL_LIST = "data/GDSC/model_list_20260709.csv"
GENE_LOOKUP = "reference/gene_lookup.parquet"
GENE_REGIME = "src/pipeline/outputs/gene_regime.parquet"
OUT_FILE    = "src/pipeline/outputs/cna_flags.parquet"

# ── Thresholds ─────────────────────────────────────────────────────────────
AMP_THRESHOLD = 2.5
DEL_THRESHOLD = 1.5

# ── 1. Load CNA ────────────────────────────────────────────────────────────
print("Loading CNA file...")
cna = pd.read_csv(CNA_FILE, sep="\t", usecols=[
    "GENE_SYMBOL", "SAMPLE_NAME", "TOTAL_CN", "MUT_TYPE"
])
print(f"  Raw rows: {len(cna):,}  |  genes: {cna.GENE_SYMBOL.nunique():,}  |  lines: {cna.SAMPLE_NAME.nunique():,}")
print(f"  TOTAL_CN range: min={cna.TOTAL_CN.min()}  max={cna.TOTAL_CN.max()}  mean={cna.TOTAL_CN.mean():.2f}")
print(f"  MUT_TYPE values: {cna.MUT_TYPE.value_counts().to_dict()}")

# ── 2. Bridge SAMPLE_NAME → ACH model_id ───────────────────────────────────
print("\nBridging SAMPLE_NAME → ACH...")
model = pd.read_csv(MODEL_LIST, usecols=["model_name", "BROAD_ID"])
model = model.dropna(subset=["BROAD_ID"])
model["model_id"]    = model["BROAD_ID"].str.lower()
model["sample_name"] = model["model_name"].str.strip()
# Drop duplicate model_names (e.g. MS-1 maps to two ACH IDs — keep first)
dupes = model.sample_name.duplicated(keep=False).sum()
if dupes:
    print(f"  Warning: {dupes} duplicate model_name rows — keeping first ACH per name")
model = model.drop_duplicates(subset="sample_name", keep="first")

before = len(cna)
cna = cna.merge(
    model[["sample_name", "model_id"]],
    left_on="SAMPLE_NAME", right_on="sample_name", how="left"
)
if len(cna) != before:
    raise RuntimeError(f"Left join still changed row count: {before} → {len(cna)}")
mapped = cna["model_id"].notna().sum()
print(f"  Bridged: {mapped:,}/{before:,} ({mapped/before*100:.1f}%)")
cna = cna[cna["model_id"].notna()].copy()

# ── 3. Join GENE_SYMBOL → ensg_id ──────────────────────────────────────────
print("\nJoining gene symbols → ENSG...")
gl = pd.read_parquet(GENE_LOOKUP, columns=["ensg_id", "hgnc_symbol", "gene_role"])
gl["ensg_id"] = gl["ensg_id"].str.split(".").str[0].str.lower()   # canonical case
gl["hgnc_symbol"] = gl["hgnc_symbol"].astype("string").str.lower()

# COSMIC writes GENE_SYMBOL uppercase; gene_lookup is normalised to lowercase
# above. Both sides must be lowered or the join matches nothing -- and because
# this is a LEFT join, that failure produces 0 mapped rows rather than an error.
cna["GENE_SYMBOL"] = cna["GENE_SYMBOL"].astype("string").str.strip().str.lower()

before = len(cna)
cna = cna.merge(gl, left_on="GENE_SYMBOL", right_on="hgnc_symbol", how="left")
assert len(cna) == before, f"Left join changed row count: {before} → {len(cna)}"
mapped = cna["ensg_id"].notna().sum()
assert mapped > 0, ("gene symbol join matched 0 rows -- GENE_SYMBOL/hgnc_symbol "
                    "case mismatch")
print(f"  Gene mapped: {mapped:,}/{before:,} ({mapped/before*100:.1f}%)")
cna = cna[cna["ensg_id"].notna()].copy()

# ── 4. Load gene regime (column is 'class', not 'regime_class') ────────────
print("\nLoading gene regime...")
regime = pd.read_parquet(GENE_REGIME, columns=["ensg_id", "class"])
regime = regime.rename(columns={"class": "regime_class"})
regime["ensg_id"] = regime["ensg_id"].str.split(".").str[0].str.lower()
cna = cna.merge(regime, on="ensg_id", how="left")
cna["regime_class"] = cna["regime_class"].fillna("unknown_dual")
cna["gene_role"]    = cna["gene_role"].fillna("unknown")

# ── 5. Apply CN thresholds and gene-class gating ───────────────────────────
print("\nApplying CN thresholds and gene-class gates...")
cna["is_amplification"] = cna["TOTAL_CN"] > AMP_THRESHOLD
cna["is_deletion"]      = cna["TOTAL_CN"] < DEL_THRESHOLD

cna["has_cna_alteration"] = (
    (cna["is_amplification"] & cna["gene_role"].isin(["oncogene", "both"])) |
    (cna["is_deletion"]      & cna["gene_role"].isin(["tsg", "both"])) |
    (cna["is_amplification"] & cna["regime_class"].isin(["activation_driven", "abundance_tracking"])) |
    (cna["is_deletion"]      & cna["regime_class"] == "loss_of_function")
)

cna["has_cna_context"] = (
    (cna["is_amplification"] | cna["is_deletion"]) &
    ~cna["has_cna_alteration"]
)

print(f"  is_amplification: {cna.is_amplification.sum():,}")
print(f"  is_deletion:      {cna.is_deletion.sum():,}")
print(f"  has_cna_alteration (role-gated): {cna.has_cna_alteration.sum():,}")
print(f"  has_cna_context (unknown):       {cna.has_cna_context.sum():,}")

# ── 6. Collapse to one row per (model_id, ensg_id) ─────────────────────────
print("\nCollapsing to one row per (model_id, ensg_id)...")
flags = (
    cna.groupby(["model_id", "ensg_id"])
    .agg(
        has_cna_alteration=("has_cna_alteration", "any"),
        has_cna_context   =("has_cna_context",    "any"),
        max_cn            =("TOTAL_CN",            "max"),
        min_cn            =("TOTAL_CN",            "min"),
        cna_type          =("MUT_TYPE", lambda x: "both" if x.nunique() > 1 else x.iloc[0])
    )
    .reset_index()
)
print(f"  Output pairs: {len(flags):,}")
print(f"  has_cna_alteration=True: {flags.has_cna_alteration.sum():,}")
print(f"  has_cna_context=True:    {flags.has_cna_context.sum():,}")
print(f"  Unique cell lines: {flags.model_id.nunique():,}")
print(f"  Unique genes:      {flags.ensg_id.nunique():,}")

# ── 7. Spot-check key genes ────────────────────────────────────────────────
print("\nSpot checks:")
gl_sym = gl.set_index("hgnc_symbol")["ensg_id"].to_dict()
for symbol in ["erbb2", "myc", "pten", "rb1", "braf", "kras", "tp53"]:
    ensg = gl_sym.get(symbol, "not found")
    if ensg == "not found":
        print(f"  {symbol}: not in gene_lookup")
        continue
    rows = flags[flags.ensg_id == ensg]
    if rows.empty:
        print(f"  {symbol} ({ensg}): no CNA data")
        continue
    alt = rows.has_cna_alteration.sum()
    ctx = rows.has_cna_context.sum()
    max_cn = rows.max_cn.max()
    min_cn = rows.min_cn.min()
    print(f"  {symbol} ({ensg}): {len(rows)} lines | alteration={alt} | context={ctx} | CN range [{min_cn}–{max_cn}]")

# ── 8. Write output ────────────────────────────────────────────────────────
Path(OUT_FILE).parent.mkdir(parents=True, exist_ok=True)
flags.to_parquet(OUT_FILE, index=False)
print(f"\nWritten: {OUT_FILE}  ({len(flags):,} rows)")
print("Done.")
