import pandas as pd
import numpy as np

live = pd.read_parquet("results/protein_z_v1_live_snapshot.parquet")
shrunk = pd.read_parquet("results/protein_z_v3_lineage_shrunk.parquet")

print("="*70)
print("STRUCTURAL CHECK")
print("="*70)
print(f"live:   {len(live):,} rows  genes={live.gene_id.nunique():,}  lines={live.model_id.nunique():,}  "
      f"nulls(z_score)={live.z_score.isna().sum():,}")
print(f"shrunk: {len(shrunk):,} rows  genes={shrunk.gene_id.nunique():,}  lines={shrunk.model_id.nunique():,}  "
      f"nulls(z_score)={shrunk.z_score.isna().sum():,}")

m = live[["gene_id","model_id","lineage","z_score","z_raw","n_sources"]].merge(
    shrunk[["gene_id","model_id","z_score"]], on=["gene_id","model_id"],
    suffixes=("_live","_shrunk"), how="outer", indicator=True)
print(f"\nmerge indicator: {m._merge.value_counts().to_dict()}")
print(f"null pattern match (both null or both non-null): "
      f"{(m.z_score_live.isna() == m.z_score_shrunk.isna()).mean()*100:.4f}%")

print("\n" + "="*70)
print("DISTRIBUTION COMPARISON")
print("="*70)
def stats(s, label):
    s = s.dropna()
    print(f"{label}: n={len(s):,}  SD={s.std():.4f}  "
          f"|z|>=5: {(s.abs()>=5).mean()*100:.4f}%  |z|>=3: {(s.abs()>=3).mean()*100:.4f}%")

stats(live.z_score, "live  ")
stats(shrunk.z_score, "shrunk")

print("\n" + "="*70)
print("SMALL vs LARGE LINEAGE: mean|delta z|")
print("="*70)
m2 = m.dropna(subset=["z_score_live","z_score_shrunk"]).copy()
m2["delta"] = (m2.z_score_shrunk - m2.z_score_live).abs()
lineage_size = live.groupby("lineage")["model_id"].nunique().rename("n_lines_in_lineage")
m2 = m2.merge(lineage_size, left_on="lineage", right_index=True, how="left")

median_size = lineage_size.median()
small = m2[m2.n_lines_in_lineage <= median_size]
large = m2[m2.n_lines_in_lineage > median_size]
print(f"lineage size median (n cell lines): {median_size:.0f}")
print(f"small lineages (<= median): n_rows={len(small):,}  mean|delta z|={small.delta.mean():.4f}")
print(f"large lineages (>  median): n_rows={len(large):,}  mean|delta z|={large.delta.mean():.4f}")

print("\nBy lineage-size decile:")
m2["size_decile"] = pd.qcut(m2.n_lines_in_lineage, 10, duplicates="drop")
print(m2.groupby("size_decile", observed=True).agg(
    n_lines=("n_lines_in_lineage","first"), mean_abs_delta=("delta","mean"), n_rows=("delta","size")
).to_string())

print("\n" + "="*70)
print("WORKED EXAMPLES -- small lineage, largest deltas")
print("="*70)
small_sorted = small.sort_values("delta", ascending=False).head(8)
print(small_sorted[["gene_id","model_id","lineage","n_lines_in_lineage","z_score_live","z_score_shrunk","delta"]].to_string(index=False))

small_sorted.to_parquet("results/_worked_examples_small_lineage.parquet")
m2.to_parquet("results/_delta_by_row.parquet")
print("\nSaved delta comparison.")
