"""
T7 — D24: CNA attrition waterfall.
Traces how many (gene, line) pairs survive each filter in cna_layer.py.
Run from repo root: python diagnostics/T7_D24_cna_waterfall.py
"""
import sys
from pathlib import Path
import pandas as pd
import duckdb

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import DB, GENE_LKP as GENE_LOOKUP, CNA_FLAGS

print("=" * 72)
print("T7 — D24: CNA ATTRITION WATERFALL")
print("=" * 72)

gl = pd.read_parquet(GENE_LOOKUP, columns=["ensg_id","gene_role"])
gl["ensg_id"] = gl["ensg_id"].str.lower()
n_genes = len(gl)
print(f"[COUNT] Gene lookup rows: {n_genes:,}")
print(f"[COUNT] Role breakdown: {gl.gene_role.value_counts().to_dict()}")

con = duckdb.connect(str(DB), read_only=True)

# Step 1: all raw COSMIC CNA rows for our genes
raw = con.execute("""
    SELECT lower(gene_id) as ensg_id, lower(model_id) as model_id,
           cna_call, total_cn
    FROM main.cosmic_cna
    WHERE lower(gene_id) IN (SELECT lower(ensg_id) FROM parquet_scan(?)
                              WHERE gene_role != 'unknown')
""", [str(GENE_LOOKUP)]).df()
con.close()

print(f"\n[COUNT] Step 1 — raw COSMIC rows for screened genes: {len(raw):,}")
print(f"  distinct genes: {raw.ensg_id.nunique():,}")
print(f"  distinct models: {raw.model_id.nunique():,}")

cna_call_dist = raw.cna_call.value_counts().head(10).to_dict()
print(f"  cna_call values (top 10): {cna_call_dist}")

# Step 2: join gene_role
merged = raw.merge(gl, on="ensg_id", how="inner")
print(f"\n[COUNT] Step 2 — after gene_role join: {len(merged):,}")

# Step 3: amplification filter
merged["is_amp"] = merged["cna_call"].str.lower().str.contains("amplif", na=False)
merged["is_del"] = merged["cna_call"].str.lower().str.contains("delet|loss|homozyg", na=False)
amp_rows = merged["is_amp"].sum()
del_rows = merged["is_del"].sum()
print(f"\n[COUNT] Step 3 — cna_call contains 'amplif': {amp_rows:,}")
print(f"[COUNT] Step 3 — cna_call contains deletion keyword: {del_rows:,}")

# Step 4: direction + role gate
onc_amp = (merged["is_amp"]) & (merged["gene_role"].isin(["oncogene","both"]))
tsg_del = (merged["is_del"]) & (merged["gene_role"].isin(["tsg","both"]))
has_alt = onc_amp | tsg_del
n_pass = has_alt.sum()
print(f"\n[COUNT] Step 4 — oncogene AND amp: {onc_amp.sum():,}")
print(f"[COUNT] Step 4 — TSG AND del: {tsg_del.sum():,}")
print(f"[COUNT] Step 4 — total has_cna_alteration: {n_pass:,}")
print(f"[MEASURED] Attrition: {len(raw):,} → {len(merged):,} → {n_pass:,}")
print(f"[MEASURED] Pass rate: {100*n_pass/len(raw):.2f}%")

# Load the saved flags and compare
flags = pd.read_parquet(CNA_FLAGS)
flags_true = flags["has_cna_alteration"].sum()
print(f"\n[COUNT] Saved cna_flags.parquet has_cna_alteration=TRUE: {flags_true:,}")
print(f"[MEASURED] Expected ~2,925 from prior run; difference: {abs(flags_true - n_pass):,}")

# Breakdown by role
for role in ["oncogene","tsg","both"]:
    sub = merged[merged["gene_role"]==role]
    if role == "oncogene":
        n = (sub["is_amp"]).sum()
    elif role == "tsg":
        n = (sub["is_del"]).sum()
    else:
        n = (sub["is_amp"] | sub["is_del"]).sum()
    print(f"[MEASURED] {role}: {len(sub):,} rows → {n:,} pass ({100*n/max(len(sub),1):.1f}%)")

# ── Ploidy correlation (remediation item 5) ──────────────────────────────────
# r(mean_total_CN, alteration_rate) — does higher ploidy explain CNA call rate?
print("\n--- PLOIDY CORRELATION (item 5) ---")
raw_num = raw.copy()
raw_num["total_cn"] = pd.to_numeric(raw_num["total_cn"], errors="coerce")
raw_num = raw_num.dropna(subset=["total_cn"])

# Per-gene: mean_total_CN and alteration_rate (fraction of alt rows with amp/del in driver)
gene_cn = raw_num.groupby("ensg_id").agg(
    mean_cn=("total_cn","mean"),
    n_rows=("total_cn","count")
).reset_index()

driver_genes = gl[gl["gene_role"].isin(["oncogene","tsg","both"])]["ensg_id"].tolist()
merged_drv = merged[merged["ensg_id"].isin(driver_genes)].copy()
gene_alt = merged_drv.groupby("ensg_id").agg(
    n_total=("model_id","count"),
    n_alt=("is_amp","sum")  # amp or del
).reset_index()
gene_alt["n_alt"] = merged_drv.groupby("ensg_id")["is_del"].sum().values + gene_alt["n_alt"].values
gene_alt["alt_rate"] = gene_alt["n_alt"] / gene_alt["n_total"]

ploidy_df = gene_cn.merge(gene_alt, on="ensg_id", how="inner")
print(f"[COUNT] Genes with both CN and driver alt data: {len(ploidy_df):,}")

if len(ploidy_df) > 5:
    from scipy import stats
    r_sp, p_sp = stats.spearmanr(ploidy_df["mean_cn"], ploidy_df["alt_rate"])
    r_pe, p_pe = stats.pearsonr(ploidy_df["mean_cn"], ploidy_df["alt_rate"])
    print(f"[MEASURED] Spearman r(mean_total_CN, driver_alt_rate) = {r_sp:.4f}  p={p_sp:.4e}")
    print(f"[MEASURED] Pearson  r(mean_total_CN, driver_alt_rate) = {r_pe:.4f}  p={p_pe:.4e}")
    print(f"[MEASURED] mean CN by cna_call:")
    for call, grp in raw_num.groupby("cna_call"):
        print(f"  {call}: mean_CN={grp.total_cn.mean():.2f}  n={len(grp):,}")

print("\n" + "=" * 72)
print("VERDICT: Waterfall printed above. No pre-declared numeric threshold.")
print(f"[DOCUMENTED] Ploidy correlation above addresses remediation item 5.")
print("=" * 72)
