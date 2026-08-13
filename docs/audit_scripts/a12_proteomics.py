"""A1.2 + C6: characterise BOTH proteomics matrices on disk.
  (a) Nusinow/CCLE Gygi TMT  -> data/parquet/4_Harmonized_MS_CCLE_Gygi_subsetted.parquet
  (b) ProCan DIA             -> data/proteomics_procan/Protein_matrix_averaged_20250211.tsv
                                (pipeline copy: src/pipeline/outputs/procan_proteomics.parquet)
Questions: ratios or intensities? histones present? cross-sample normalised? depth?
"""
import os, re, numpy as np, pandas as pd
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
pd.set_option("display.width", 200)

HIST = re.compile(r"^(H2A|H2B|H3|H4|HIST\d|H1-\d|H2AC|H2BC|H3C|H4C|H3-\d|H2AZ|H2AX|MACROH2A)", re.I)


def describe(name, df, idcol=None):
    num = df.select_dtypes(include=[np.number])
    v = num.to_numpy(dtype=float)
    flat = v[np.isfinite(v)]
    print(f"\n{'='*78}\n{name}\n{'='*78}")
    print(f"shape                : {df.shape[0]:,} rows x {df.shape[1]:,} cols "
          f"({num.shape[1]:,} numeric)")
    print(f"finite values        : {len(flat):,}   missing: "
          f"{100*(1-len(flat)/max(v.size,1)):.1f}%")
    print(f"range                : min={flat.min():.3f}  p1={np.percentile(flat,1):.3f}  "
          f"median={np.median(flat):.3f}  p99={np.percentile(flat,99):.3f}  max={flat.max():.3f}")
    print(f"mean={flat.mean():.4f}  sd={flat.std():.4f}   frac<0 = {100*(flat<0).mean():.1f}%")
    return num


# ---------------------------------------------------------------- (a) TMT
tmt = pd.read_parquet("data/parquet/4_Harmonized_MS_CCLE_Gygi_subsetted.parquet")
print("TMT first 4 col names:", list(tmt.columns[:4]))
num_t = describe("(a) Nusinow / CCLE Gygi  TMT", tmt)
print("\nINTERPRETATION: values centred on ~0 spanning negative and positive with unit-ish")
print("sd are LOG-RATIOS to a bridge channel, not intensities.")
# per-sample centring evidence: TMT matrix is rows=lines, cols=proteins
rowmed = num_t.median(axis=1); colmed = num_t.median(axis=0)
print(f"per-ROW (cell line) median: median={rowmed.median():.4f} sd={rowmed.std():.4f} "
      f"range=[{rowmed.min():.3f},{rowmed.max():.3f}]")
print(f"per-COL (protein)  median: median={colmed.median():.4f} sd={colmed.std():.4f}")
print("-> per-row medians pinned near a constant == cross-sample median centring applied.")

tmt_hist = [c for c in tmt.columns if HIST.match(str(c).split(" ")[0].split("_")[0])]
print(f"\nhistone-matching columns in TMT: {len(tmt_hist)}")
print("  e.g.", tmt_hist[:12])

# ---------------------------------------------------------------- (b) ProCan DIA
pro = pd.read_parquet("src/pipeline/outputs/procan_proteomics.parquet")
print("\nProCan first 4 col names:", list(pro.columns[:4]))
num_p = describe("(b) ProCan DIA (pipeline copy)", pro)
rowmed_p = num_p.median(axis=1); colmed_p = num_p.median(axis=0)
print(f"per-ROW (cell line) median: median={rowmed_p.median():.4f} sd={rowmed_p.std():.4f} "
      f"range=[{rowmed_p.min():.3f},{rowmed_p.max():.3f}]")
print(f"per-COL (protein)  median: median={colmed_p.median():.4f} sd={colmed_p.std():.4f}")

pro_hist = [c for c in pro.columns if HIST.match(str(c).split(";")[0])]
print(f"\nhistone-matching columns in ProCan: {len(pro_hist)}")
print("  ", sorted(pro_hist)[:40])

# depth per sample
depth_p = num_p.notna().sum(axis=1)
print(f"\nProCan proteins quantified per cell line: min={depth_p.min():,} "
      f"median={depth_p.median():,.0f} max={depth_p.max():,}")
depth_t = num_t.notna().sum(axis=1)
print(f"TMT    proteins quantified per cell line: min={depth_t.min():,} "
      f"median={depth_t.median():,.0f} max={depth_t.max():,}")

# histone share of total signal per sample -- the ruler's key quantity
if pro_hist:
    hs = num_p[[c for c in pro_hist if c in num_p.columns]]
    if hs.shape[1]:
        lin = np.power(2.0, num_p)          # if values are log2 intensities
        tot = lin.sum(axis=1, skipna=True)
        hsum = np.power(2.0, hs).sum(axis=1, skipna=True)
        ratio = (hsum / tot).replace([np.inf, -np.inf], np.nan).dropna()
        print(f"\nhistone / total intensity ratio (assuming log2 intensities):")
        print(f"  n={len(ratio)}  median={ratio.median():.5f}  "
              f"IQR={ratio.quantile(.25):.5f}-{ratio.quantile(.75):.5f}  "
              f"CV={ratio.std()/ratio.mean():.3f}")
