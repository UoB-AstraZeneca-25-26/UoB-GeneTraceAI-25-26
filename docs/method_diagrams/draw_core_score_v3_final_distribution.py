import duckdb
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CORE = "C:/Disertation/UoB-GeneTraceAI-25-26/final_pipeline/outputs/core_score.parquet"  # = core_score_v3, now live

con = duckdb.connect()

overall = con.execute(f"""
    SELECT count(*) n, avg(core_score) mean, stddev_samp(core_score) sd, median(core_score) med,
           avg(CASE WHEN core_score<0.01 THEN 1.0 ELSE 0 END) lo,
           avg(CASE WHEN core_score>0.99 THEN 1.0 ELSE 0 END) hi,
           avg(CASE WHEN core_score>=0.499 AND core_score<=0.501 THEN 1.0 ELSE 0 END) mid
    FROM read_parquet('{CORE}') WHERE core_score IS NOT NULL
""").fetchdf().iloc[0]
print(f"n={int(overall.n):,}  mean={overall['mean']:.4f}  median={overall.med:.4f}  sd={overall.sd:.4f}")
print(f"<0.01: {overall.lo*100:.3f}%  >0.99: {overall.hi*100:.3f}%  "
      f"combined: {(overall.lo+overall.hi)*100:.3f}%  ~0.5 cluster: {overall.mid*100:.3f}%")

by_layer = con.execute(f"""
    SELECT n_layers, count(*) n, avg(core_score) mean, stddev_samp(core_score) sd
    FROM read_parquet('{CORE}') WHERE core_score IS NOT NULL GROUP BY n_layers ORDER BY n_layers
""").fetchdf()
print("\nBy n_layers:")
print(by_layer.to_string(index=False))

samp = con.execute(f"SELECT core_score, n_layers FROM read_parquet('{CORE}') WHERE core_score IS NOT NULL USING SAMPLE 400000 ROWS").fetchdf()

BG, BOX_BG, TEXT, GRID = "#0a0d11", "#121820", "#dbe2e9", "#2a333d"
plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": BOX_BG, "axes.edgecolor": GRID,
    "text.color": TEXT, "axes.labelcolor": TEXT, "xtick.color": TEXT,
    "ytick.color": TEXT, "font.family": "DejaVu Sans", "grid.color": GRID,
})
fig, axes = plt.subplots(1, 3, figsize=(15.5, 5))

axes[0].hist(samp.core_score, bins=200, range=(0,1), color="#c792ea", density=True)
axes[0].set_title(f"core_score (LIVE, v3 fix)\nn={int(overall.n):,}  mean={overall['mean']:.3f} sd={overall.sd:.3f}")
axes[0].set_xlabel("core_score")

s2 = samp[samp.n_layers==2].core_score
axes[1].hist(s2, bins=200, range=(0,1), color="#72d9be", density=True)
axes[1].set_title(f"n_layers=2 (RNA + protein)\nmean={s2.mean():.3f} sd={s2.std():.3f}")
axes[1].set_xlabel("core_score")

s1 = samp[samp.n_layers==1].core_score
axes[2].hist(s1, bins=200, range=(0,1), color="#e8ab68", density=True)
axes[2].set_title(f"n_layers=1 (RNA-only fallback)\nmean={s1.mean():.3f} sd={s1.std():.3f}")
axes[2].set_xlabel("core_score")

fig.suptitle(f"core_score distribution — LIVE (protein shrinkage+clip promoted)\n"
             f"edge spike <0.01 or >0.99: {(overall.lo+overall.hi)*100:.3f}%   "
             f"~0.5 cluster: {overall.mid*100:.3f}%",
             fontsize=13.5, color=TEXT, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.88])
fig.savefig("core_score_v3_final_distribution.png", dpi=170, facecolor=BG)
print("\nwrote core_score_v3_final_distribution.png")
