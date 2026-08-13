"""C6 Q1/Q2/Q3: histone presence, normalisation, depth -- done on SYMBOLS, by
mapping UniProt accessions properly instead of regexing the accession string."""
import os, re, numpy as np, pandas as pd
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")

HIST = re.compile(r"^(H1-\d|H2A[A-Z]?\d*|H2B[A-Z]?\d*|H3-\d|H3C\d|H4C\d|H2AC\d|H2BC\d|"
                  r"H2AZ\d|H2AX|H2AJ|H3\d?$|H4$|MACROH2A\d|HIST\d)", re.I)

# ---- symbol source: HGNC protein-product file (accession -> symbol) ----
hg = pd.read_csv("data/gene_with_protein_product.txt", sep="\t", low_memory=False)
acc2sym = {}
for sym, uid, uids in zip(hg["symbol"], hg.get("uniprot_ids", ""), hg.get("uniprot_ids", "")):
    if isinstance(uids, str):
        for a in uids.split("|"):
            a = a.strip().lower()
            if a:
                acc2sym[a] = sym
print(f"HGNC accession->symbol map: {len(acc2sym):,} accessions")
hgnc_hist = sorted({s for s in hg["symbol"] if HIST.match(str(s))})
print(f"HGNC histone symbols matched by pattern: {len(hgnc_hist)}  e.g. {hgnc_hist[:12]}")

# ---------------- ProCan (columns are lowercase uniprot accessions) -------------
pro = pd.read_parquet("src/pipeline/outputs/procan_proteomics.parquet")
meta = [c for c in pro.columns if not re.fullmatch(r"[a-z0-9\-]+", str(c))] + \
       ["gdsc_model_name", "sanger_model_id"]
acc_cols = [c for c in pro.columns if c not in set(meta)]
pro_sym = {c: acc2sym.get(str(c).split("-")[0], None) for c in acc_cols}
mapped = {c: s for c, s in pro_sym.items() if s}
print(f"\nProCan: {len(acc_cols):,} accession columns, {len(mapped):,} mapped to an HGNC symbol")
pro_hist = sorted([c for c, s in mapped.items() if HIST.match(s)])
print(f"ProCan HISTONE columns: {len(pro_hist)}")
for c in pro_hist:
    n = pro[c].notna().sum()
    print(f"   {c:12s} -> {mapped[c]:10s}  quantified in {n:4d}/{len(pro)} lines "
          f"({100*n/len(pro):5.1f}%)  median={pro[c].median():.3f}")

# ---------------- TMT (columns are 'ACC (SYMBOL)') -----------------------------
tmt = pd.read_parquet("data/parquet/4_Harmonized_MS_CCLE_Gygi_subsetted.parquet")
tsym = {}
for c in tmt.columns[1:]:
    m = re.match(r"^(\S+)\s*\((.+)\)$", str(c))
    tsym[c] = m.group(2) if m else acc2sym.get(str(c).split("-")[0].lower())
tmt_hist = sorted([c for c, s in tsym.items() if s and HIST.match(s)])
print(f"\nTMT HISTONE columns: {len(tmt_hist)}")
for c in tmt_hist[:25]:
    n = tmt[c].notna().sum()
    print(f"   {c:28s} -> {tsym[c]:10s} quantified in {n:4d}/{len(tmt)} ({100*n/len(tmt):5.1f}%)")

# ---------------- C6 Q2: normalisation state of the RAW ProCan TSV -------------
raw = pd.read_csv("data/proteomics_procan/Protein_matrix_averaged_20250211.tsv",
                  sep="\t", nrows=6000, low_memory=False)
print(f"\nRAW ProCan TSV (first {len(raw):,} rows): shape {raw.shape}")
print("first cols:", list(raw.columns[:6]))
rn = raw.select_dtypes(include=[np.number])
arr = rn.to_numpy(dtype=float); fl = arr[np.isfinite(arr)]
print(f"  values: min={fl.min():.3f} median={np.median(fl):.3f} max={fl.max():.3f} "
      f"frac<0={100*(fl<0).mean():.1f}%  missing={100*(1-len(fl)/arr.size):.1f}%")
colmed = rn.median(axis=0)
print(f"  per-COLUMN median across sample columns: median={colmed.median():.3f} "
      f"sd={colmed.std():.3f} range=[{colmed.min():.3f},{colmed.max():.3f}]")
print("  -> a WIDE spread of per-sample medians means between-sample intensity")
print("     differences survive; a pinned spread means median-scaling was applied.")

# ---------------- C6 Q3: histone fraction stability ---------------------------
if pro_hist:
    lin = np.power(2.0, pro[acc_cols].astype(float))
    tot = lin.sum(axis=1, min_count=1)
    hs = np.power(2.0, pro[pro_hist].astype(float)).sum(axis=1, min_count=1)
    ratio = (hs / tot).replace([np.inf, -np.inf], np.nan).dropna()
    depth = pro[acc_cols].notna().sum(axis=1)
    print(f"\nhistone/total intensity ratio: n={len(ratio)} median={ratio.median():.5f} "
          f"IQR={ratio.quantile(.25):.5f}-{ratio.quantile(.75):.5f} CV={ratio.std()/ratio.mean():.3f}")
    from scipy.stats import spearmanr
    j = pd.concat([ratio.rename("r"), depth.rename("d")], axis=1).dropna()
    r, p = spearmanr(j.r, j.d)
    print(f"Spearman(histone ratio, proteins quantified) = {r:.3f} p={p:.2e} n={len(j)}")
    print("  -> a strong depth dependence means the ratio has NOT stabilised at this depth.")
