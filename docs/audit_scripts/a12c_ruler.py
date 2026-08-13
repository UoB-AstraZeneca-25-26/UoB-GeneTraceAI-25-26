"""C6 Q2/Q3 finish: normalisation state of the raw ProCan TSV, and whether the
histone/total ratio has stabilised at the available depth."""
import os, re, numpy as np, pandas as pd
from scipy.stats import spearmanr
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")

raw = pd.read_csv("data/proteomics_procan/Protein_matrix_averaged_20250211.tsv",
                  sep="\t", low_memory=False)
print("RAW ProCan TSV shape:", raw.shape)
print("first 3 cols:", list(raw.columns[:3]), "| dtypes:", raw.dtypes.iloc[:3].tolist())
print("head of the two ID columns:\n", raw.iloc[:3, :2].to_string())

data = raw.iloc[:, 2:].apply(pd.to_numeric, errors="coerce")
arr = data.to_numpy(dtype=float); fl = arr[np.isfinite(arr)]
print(f"\nvalues: min={fl.min():.3f} p1={np.percentile(fl,1):.3f} median={np.median(fl):.3f} "
      f"p99={np.percentile(fl,99):.3f} max={fl.max():.3f}")
print(f"  frac<0 = {100*(fl<0).mean():.2f}%   missing = {100*(1-len(fl)/arr.size):.1f}%")
print("  -> overwhelmingly positive, bounded ~[0,16] => LOG2 INTENSITIES, not log-ratios.")

rowmed = data.median(axis=1)
print(f"\nper-CELL-LINE median log2 intensity: median={rowmed.median():.3f} sd={rowmed.std():.3f}")
print(f"  range=[{rowmed.min():.3f}, {rowmed.max():.3f}]  "
      f"IQR={rowmed.quantile(.25):.3f}-{rowmed.quantile(.75):.3f}")
print(f"  spread p1-p99 = {rowmed.quantile(.99)-rowmed.quantile(.01):.3f} log2 units "
      f"= {2**(rowmed.quantile(.99)-rowmed.quantile(.01)):.2f}x in linear intensity")
print("  A hard median-scaled matrix would have sd ~0 here. It does not.")

# --------- histone ratio stability vs depth ---------
hg = pd.read_csv("data/gene_with_protein_product.txt", sep="\t", low_memory=False)
acc2sym = {}
for sym, uids in zip(hg["symbol"], hg["uniprot_ids"]):
    if isinstance(uids, str):
        for a in uids.split("|"):
            if a.strip():
                acc2sym[a.strip().lower()] = sym
HIST = re.compile(r"^(H1-\d|H2A|H2B|H3-\d|H3C\d|H4C\d|H2AC\d|H2BC\d|H2AZ|H2AX|MACROH2A)", re.I)
hist_cols = [c for c in data.columns if HIST.match(str(acc2sym.get(str(c).split("-")[0].lower(), "")))]
print(f"\nhistone columns in raw TSV: {len(hist_cols)}")

lin = np.power(2.0, data)
tot = lin.sum(axis=1, min_count=1)
hs = lin[hist_cols].sum(axis=1, min_count=1)
ratio = (hs / tot).replace([np.inf, -np.inf], np.nan)
depth = data.notna().sum(axis=1)
j = pd.DataFrame({"ratio": ratio, "depth": depth}).dropna()
print(f"\nhistone/total linear-intensity ratio: n={len(j)}")
print(f"  median={j.ratio.median():.5f}  IQR={j.ratio.quantile(.25):.5f}-{j.ratio.quantile(.75):.5f}"
      f"  CV={j.ratio.std()/j.ratio.mean():.3f}")
r, p = spearmanr(j.ratio, j.depth)
print(f"  Spearman(ratio, proteins quantified) = {r:.3f}  p={p:.3e}")
print(f"  proteins/line: min={depth.min():,} median={depth.median():,.0f} max={depth.max():,}")
print("\n  Published ruler guidance: histone/total stabilises from ~12,000 PEPTIDES.")
print("  This matrix is protein-level only; peptide counts are not in the file, so")
print("  depth adequacy cannot be assessed on the published criterion directly.")

# binned to show trend
j["bin"] = pd.qcut(j.depth, 5, duplicates="drop")
print("\n  ratio by depth quintile:")
print(j.groupby("bin", observed=True).agg(n=("ratio", "size"), depth_med=("depth", "median"),
                                          ratio_med=("ratio", "median")).round(5).to_string())
