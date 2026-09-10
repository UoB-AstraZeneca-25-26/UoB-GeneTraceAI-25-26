"""
gate_audit/05_region_sweep.py
-----------------------------
AUDIT 5 of 6.  "The gate region was almost certainly chosen post-hoc."

THE CHALLENGE
-------------
FPR in [0.8, 1.0] as the gate region was identified AFTER full AUROC came back
flat at 0.52. That is post-hoc region selection on the discovery dataset, and it
is a fair challenge. Three free parameters sit inside the headline:

    the region boundary       0.8
    the decile granularity    10
    the NPV exclusion point   5%

The project's own standing rule says sweep free parameters. This does that. A
result stable across the sweep is a far stronger claim than one boundary.

THE DEFENCE, WHICH MUST BE STATED RATHER THAN LEFT TO BE FOUND
--------------------------------------------------------------
Post-hoc selection is answered by replication on data not used to choose the
region. The region was fixed BEFORE the independent screen replication
(`src/pipeline/test_run_screen_replication.py`: Broad 0.9950, Sanger 0.9960,
difference p = 0.19). Audit 6 re-runs that on the real Chronos release for both
screens. Say this explicitly in the write-up; do not let a panel find it.

WHAT THIS SCRIPT SWEEPS
-----------------------
  A. Region boundary. pAUC over FPR in [lo, 1.0] for lo in
     {0.70, 0.75, 0.80, 0.85, 0.90}, plus the symmetric ranker region for
     contrast. Stability across lo is the claim.
  B. Full 2-D grid over [lo, hi], so the reader can see whether 0.8 sits on a
     plateau or on a spike. A spike is a red flag; a plateau is the answer.
  C. Bin granularity. Bottom-bin depletion ratio for n_bins in
     {4, 5, 8, 10, 16, 20}. The decile is a convention; the depletion should not
     depend on it beyond the mechanical fact that a narrower bottom bin is more
     extreme.
  D. NPV exclusion point. Uplift over base NPV at cuts
     {1, 2, 5, 10, 15, 20, 25}%.

Every sweep is run on ALL genes and on VALID genes (audit 2), so post-hoc
sensitivity and validity sensitivity are visible in one table.

Read-only. Writes only into gate_audit/outputs/.

Run:
    python gate_audit/05_region_sweep.py
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import common as C

panel = C.read_panel()
G = pd.read_parquet(C.OUT / "genes.parquet")
valid = set(G.loc[G.is_valid, "ensg_id"])

by_gene = {g: (d.Z.values.astype(float), d.dep.values.astype(bool))
           for g, d in panel.groupby("ensg_id", sort=False)}
genes_all = sorted(by_gene)
genes_val = [g for g in genes_all if g in valid]
SETS = [("all", genes_all), ("valid", genes_val)]

res = {"prespecification_note": (
    "FPR region [0.8,1.0] was identified after full AUROC returned 0.52 on the "
    "discovery data, i.e. post-hoc. The defence is replication on data not used "
    "to select it (audit 6, Broad vs Sanger on the real Chronos release). This "
    "sweep tests whether the choice sits on a plateau or a spike.")}

C.banner("AUDIT 5 -- SWEEPING THE FREE PARAMETERS IN THE HEADLINE")
print(f"  genes: all {len(genes_all):,}   valid {len(genes_val):,}")
print(f"  NOTE: the region was chosen post-hoc. This sweep does not undo that.")
print(f"        It shows whether the result depends on the choice.")

# ------------------------------------------------------------ A
C.banner("5A -- REGION BOUNDARY SWEEP,  pAUC over FPR in [lo, 1.0]")
print(f"  {'lo':>6} {'ALL pAUC':>10} {'ALL p':>11} {'VALID pAUC':>12} "
      f"{'VALID p':>11} {'>chance':>9}")
res["region_boundary"] = {}
for lo in (0.70, 0.75, 0.80, 0.85, 0.90):
    row = {}
    cells = []
    for name, gs in SETS:
        v = np.array([C.partial_auc(*by_gene[g], lo, 1.0) for g in gs])
        v = v[np.isfinite(v)]
        p = wilcoxon(v - 0.5).pvalue if len(v) > 5 else np.nan
        row[name] = {"pauc": float(np.median(v)), "p": float(p),
                     "frac_above_chance": float((v > 0.5).mean()),
                     "n": int(len(v))}
        cells.append((np.median(v), p, (v > 0.5).mean()))
    print(f"  {lo:>6.2f} {cells[0][0]:>10.4f} {cells[0][1]:>11.3g} "
          f"{cells[1][0]:>12.4f} {cells[1][1]:>11.3g} {cells[1][2]:>9.3f}")
    res["region_boundary"][f"{lo:.2f}"] = row

pav = [res["region_boundary"][f"{lo:.2f}"]["valid"]["pauc"]
       for lo in (0.70, 0.75, 0.80, 0.85, 0.90)]
print()
print(f"  spread across boundaries (valid genes): "
      f"{min(pav):.4f} .. {max(pav):.4f}   range {max(pav)-min(pav):.4f}")
res["region_boundary_spread_valid"] = float(max(pav) - min(pav))

print()
print("  CONTRAST -- the symmetric ranker region, FPR in [0.0, hi]:")
print(f"  {'hi':>6} {'ALL pAUC':>10} {'VALID pAUC':>12}")
res["ranker_region"] = {}
for hi in (0.10, 0.20, 0.30):
    cells = []
    for name, gs in SETS:
        v = np.array([C.partial_auc(*by_gene[g], 0.0, hi) for g in gs])
        v = v[np.isfinite(v)]
        cells.append(float(np.median(v)))
    print(f"  {hi:>6.2f} {cells[0]:>10.4f} {cells[1]:>12.4f}")
    res["ranker_region"][f"{hi:.2f}"] = {"all": cells[0], "valid": cells[1]}

# ------------------------------------------------------------ B
C.banner("5B -- FULL [lo, hi] GRID  (valid genes).  Plateau or spike?")
los = [0.60, 0.70, 0.75, 0.80, 0.85, 0.90]
his = [0.90, 0.95, 1.00]
grid = {}
print(f"  {'lo\\hi':>7}" + "".join(f"{h:>10.2f}" for h in his))
for lo in los:
    cells = []
    for hi in his:
        if hi <= lo:
            cells.append(np.nan)
            continue
        v = np.array([C.partial_auc(*by_gene[g], lo, hi) for g in genes_val])
        v = v[np.isfinite(v)]
        cells.append(float(np.median(v)))
        grid[f"{lo:.2f}_{hi:.2f}"] = cells[-1]
    print(f"  {lo:>7.2f}" + "".join(
        "       n/a" if not np.isfinite(c) else f"{c:>10.4f}" for c in cells))
res["region_grid_valid"] = grid
fin = [v for v in grid.values() if np.isfinite(v)]
print()
print(f"  grid range {min(fin):.4f} .. {max(fin):.4f}   "
      f"({'PLATEAU' if max(fin)-min(fin) < 0.05 else 'NOT flat'})")
res["region_grid_range"] = float(max(fin) - min(fin))

# ------------------------------------------------------------ C
C.banner("5C -- BIN GRANULARITY.  Bottom-bin depletion ratio.")
print("  The decile is a convention. A narrower bottom bin is mechanically more")
print("  extreme; what matters is that the sign and the story do not change.")
print()
print(f"  {'bins':>6} {'bottom bin %':>13} {'ALL ratio':>11} {'VALID ratio':>13} "
      f"{'VALID p':>11}")
res["bin_granularity"] = {}
for nb in (4, 5, 8, 10, 16, 20):
    cells = []
    for name, gs in SETS:
        vals = []
        for g in gs:
            s, p = by_gene[g]
            base = p.mean()
            if base <= 0:
                continue
            q, keep = C.decile_index(s, nb, "ordinal")
            qv, pv = q[keep], p[keep]
            if (qv == 0).sum() < 5:
                continue
            vals.append(pv[qv == 0].mean() / base)
        vals = np.array(vals)
        cells.append((float(np.median(vals)),
                      float(wilcoxon(vals - 1.0).pvalue) if len(vals) > 5 else np.nan,
                      int(len(vals))))
    print(f"  {nb:>6} {100/nb:>12.1f}% {cells[0][0]:>11.4f} "
          f"{cells[1][0]:>13.4f} {cells[1][1]:>11.3g}")
    res["bin_granularity"][str(nb)] = {
        "all_ratio": cells[0][0], "valid_ratio": cells[1][0],
        "valid_p": cells[1][1], "valid_n": cells[1][2]}

# ------------------------------------------------------------ D
C.banner("5D -- NPV EXCLUSION POINT.  Uplift over base NPV.")
print("  base NPV = 1 - prevalence. Uplift is the only honest quantity, because")
print("  NPV is high whenever prevalence is low.")
print()
print(f"  {'cut %':>7} {'ALL uplift':>12} {'VALID uplift':>14} {'VALID NPV':>11} "
      f"{'positive in':>12}")
res["npv_cut"] = {}
for cut in (1, 2, 5, 10, 15, 20, 25):
    cells = []
    for name, gs in SETS:
        up, npv = [], []
        for g in gs:
            s, p = by_gene[g]
            k = max(1, int(round(len(s) * cut / 100)))
            order = np.argsort(s, kind="mergesort")
            n = 1.0 - p[order[:k]].mean()
            npv.append(n)
            up.append(n - (1.0 - p.mean()))
        cells.append((float(np.median(up)), float(np.median(npv)),
                      float(np.mean(np.array(up) > 0))))
    print(f"  {cut:>7} {cells[0][0]:>+12.5f} {cells[1][0]:>+14.5f} "
          f"{cells[1][1]:>11.4f} {100*cells[1][2]:>11.1f}%")
    res["npv_cut"][str(cut)] = {"all_uplift": cells[0][0],
                                "valid_uplift": cells[1][0],
                                "valid_npv": cells[1][1],
                                "valid_frac_positive": cells[1][2]}

# ------------------------------------------------------------ verdict
C.banner("VERDICT")
spread = res["region_boundary_spread_valid"]
grng = res["region_grid_range"]
stable = spread < 0.02 and grng < 0.05
print(f"  pAUC across boundaries 0.70-0.90 (valid) : range {spread:.4f}")
print(f"  pAUC across the full [lo,hi] grid        : range {grng:.4f}")
print()
if stable:
    verdict = "stable_across_the_sweep__report_the_range_not_the_point"
    print("  The result is a PLATEAU, not a spike. 0.8 is not doing any work, and")
    print("  that is the strong form of the claim. Report the range across")
    print("  boundaries rather than the single boundary, and state the")
    print("  pre-specification defence (replication in audit 6) explicitly.")
else:
    verdict = "sensitive_to_the_region__the_post_hoc_challenge_lands"
    print("  The result moves materially with the boundary. The post-hoc")
    print("  selection challenge lands, and the single-boundary number cannot")
    print("  carry the claim on its own.")
print(f"\n  VERDICT: {verdict}")
res["verdict"] = verdict

out = C.OUT / "05_region_sweep_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
