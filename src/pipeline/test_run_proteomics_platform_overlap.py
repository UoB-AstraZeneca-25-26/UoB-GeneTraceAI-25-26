"""
test_run_proteomics_platform_overlap.py
------------------------------------------
Read-only diagnostic. Writes no production file.

QUESTION
--------
ProCan-DepMapSanger (Goncalves et al., Cancer Cell 2022; 948 lines x 8,453
proteins, DIA-MS) would raise proteomics coverage from 375 to 773 of the 1,485
scored cell lines, taking 2-layer rows from 9.6% to ~21%. Before merging it with
the existing CCLE/Gygi proteomics (375 lines, TMT), measure the batch effect:

    On cell lines measured by BOTH platforms, do the two agree per protein?

Why it matters: Stage 2 percentile-ranks each protein across cell lines. If the
two platforms disagree, a naive merge lets a protein look "high" merely because
it was measured on the more sensitive platform, and that artifact is baked
straight into core_score.

MATCHING
--------
Both datasets are keyed by UniProt accession -- ProCan in row 1 of the matrix,
CCLE as (lowercased) column names. Matching on UniProt avoids a symbol round-trip.
ProCan SIDM ids bridge to ACH via model_list_20260709.csv (measured 99.3%).

INTERPRETATION GUIDE (decided before running)
    median rho > 0.5  -> platforms broadly agree; merging after per-dataset
                         percentile ranking is safe
    0.3 - 0.5         -> partial agreement; merge only with per-dataset
                         normalisation, and report the caveat
    < 0.3             -> platforms disagree; do NOT merge. Use ProCan alone
                         (689 lines, single platform) instead.

Output: src/pipeline/outputs/test_run_proteomics_platform_overlap_results.json
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr

OUTPUTS = Path("src/pipeline/outputs")
PROCAN  = "data/proteomics_procan/Protein_matrix_averaged_20250211.tsv"

print("=" * 72)
print("STEP 1 -- load both proteomics datasets")
print("=" * 72)

meta = pd.read_csv(PROCAN, sep="\t", nrows=2, header=None, low_memory=False)
uniprots = [str(u).strip().lower() for u in meta.iloc[0, 2:].tolist()]
symbols  = [str(s).strip() for s in meta.iloc[1, 2:].tolist()]

pc = pd.read_csv(PROCAN, sep="\t", skiprows=3, header=None, low_memory=False)
pc_lines = pc.iloc[:, 1].astype(str).str.strip()          # SIDM
pc_vals = pc.iloc[:, 2:].astype("float32")
pc_vals.columns = uniprots
pc_vals.index = pc_lines
print(f"  ProCan : {pc_vals.shape[0]:,} lines x {pc_vals.shape[1]:,} proteins")

ccle = pd.read_parquet("cleaned_track_data/proteomics.parquet")
ccle = ccle.set_index(ccle.columns[0])
ccle.index = ccle.index.astype(str).str.lower()
ccle.columns = [str(c).strip().lower() for c in ccle.columns]
print(f"  CCLE   : {ccle.shape[0]:,} lines x {ccle.shape[1]:,} proteins")

print()
print("=" * 72)
print("STEP 2 -- bridge ProCan SIDM -> ACH and find shared lines")
print("=" * 72)

ml = pd.read_csv("data/GDSC/model_list_20260709.csv", low_memory=False).dropna(subset=["BROAD_ID"])
sidm2ach = dict(zip(ml.model_id, ml.BROAD_ID.str.lower()))
pc_vals.index = [sidm2ach.get(s) for s in pc_vals.index]
pc_vals = pc_vals[[i is not None for i in pc_vals.index]]
pc_vals = pc_vals[~pc_vals.index.duplicated(keep="first")]
print(f"  ProCan bridged to ACH: {pc_vals.shape[0]:,} lines")

shared_lines = sorted(set(pc_vals.index) & set(ccle.index))
shared_prot  = sorted(set(pc_vals.columns) & set(ccle.columns))
print(f"  lines measured by BOTH platforms : {len(shared_lines):,}")
print(f"  proteins present in BOTH         : {len(shared_prot):,}")

if len(shared_lines) < 20 or len(shared_prot) < 100:
    raise SystemExit("Too little overlap to assess a batch effect.")

A = pc_vals.loc[shared_lines, shared_prot].to_numpy(dtype="float64")
B = ccle.loc[shared_lines, shared_prot].to_numpy(dtype="float64")

print()
print("=" * 72)
print("STEP 3 -- per-protein agreement across shared lines")
print("=" * 72)

rhos, ns, prot_names = [], [], []
for j in range(A.shape[1]):
    m = ~np.isnan(A[:, j]) & ~np.isnan(B[:, j])
    n = int(m.sum())
    if n < 20:
        continue
    a, b = A[m, j], B[m, j]
    if np.all(a == a[0]) or np.all(b == b[0]):
        continue
    r, _ = spearmanr(a, b)
    if not np.isnan(r):
        rhos.append(float(r)); ns.append(n); prot_names.append(shared_prot[j])

rhos = np.array(rhos)
print(f"  proteins assessed (n>=20 shared lines): {len(rhos):,}")
print(f"  median rho : {np.median(rhos):+.4f}")
print(f"  mean   rho : {rhos.mean():+.4f}")
for q in [10, 25, 50, 75, 90]:
    print(f"    p{q:<2d} : {np.percentile(rhos, q):+.4f}")
print(f"  proteins with rho > 0.5 : {(rhos > 0.5).sum():,}  ({(rhos>0.5).mean()*100:.1f}%)")
print(f"  proteins with rho > 0.3 : {(rhos > 0.3).sum():,}  ({(rhos>0.3).mean()*100:.1f}%)")
print(f"  proteins with rho < 0   : {(rhos < 0).sum():,}  ({(rhos<0).mean()*100:.1f}%)")

med = float(np.median(rhos))
if med > 0.5:
    verdict = ("PLATFORMS AGREE (median rho {:.3f}). Merge is safe provided each dataset is "
               "percentile-ranked separately before combination, which Stage 2 already does "
               "per gene.".format(med))
    action = "merge_with_per_dataset_ranking"
elif med > 0.3:
    verdict = ("PARTIAL AGREEMENT (median rho {:.3f}). Merge only with per-dataset "
               "normalisation, and report the batch effect as a caveat.".format(med))
    action = "merge_with_caution"
else:
    verdict = ("PLATFORMS DISAGREE (median rho {:.3f}). Do NOT merge -- use ProCan alone "
               "(689 scored lines, single platform) rather than mixing.".format(med))
    action = "use_procan_alone"
print()
print(f"  VERDICT: {verdict}")

print()
print("=" * 72)
print("STEP 4 -- per-protein platform conflict tiers")
print("=" * 72)
print("  consistent  (rho >= 0.5) : safe to merge after per-dataset percentile ranking")
print("  cautious    (0.3 - 0.5)  : merge with caveat; flag in output")
print("  conflicting (rho < 0.3)  : use only one platform's score for this protein")
print()

rho_arr = rhos   # already a numpy array after Step 3
tier_labels = np.where(rho_arr >= 0.5, "consistent",
              np.where(rho_arr >= 0.3, "cautious", "conflicting"))

tier_df = pd.DataFrame({
    "uniprot":       prot_names,
    "platform_rho":  rho_arr,
    "platform_tier": tier_labels,
})
n_consistent  = int((tier_labels == "consistent").sum())
n_cautious    = int((tier_labels == "cautious").sum())
n_conflicting = int((tier_labels == "conflicting").sum())
print(f"  consistent  (rho >= 0.5) : {n_consistent:,}  ({n_consistent/len(rhos)*100:.1f}%)")
print(f"  cautious    (0.3 - 0.5)  : {n_cautious:,}  ({n_cautious/len(rhos)*100:.1f}%)")
print(f"  conflicting (rho <  0.3) : {n_conflicting:,}  ({n_conflicting/len(rhos)*100:.1f}%)")
tier_out = OUTPUTS / "protein_platform_tier.parquet"
tier_df.to_parquet(tier_out, index=False)
print(f"  Saved: {tier_out}")

print()
print("=" * 72)
print("STEP 5 -- per-cell-line cross-protein agreement")
print("=" * 72)
print("  For each of the 291 shared lines: Spearman rho across all shared proteins.")
print("  Lines with rho < 0.3 are suspect 2-layer entries; flag or exclude.")
print()

prot_to_rho = dict(zip(prot_names, rhos))
line_records = []
for line in shared_lines:
    a = pc_vals.loc[line, shared_prot].to_numpy(dtype="float64")
    b = ccle.loc[line, shared_prot].to_numpy(dtype="float64")
    mask = ~np.isnan(a) & ~np.isnan(b)
    if mask.sum() < 50:
        continue
    av, bv = a[mask], b[mask]
    if np.all(av == av[0]) or np.all(bv == bv[0]):
        continue
    r, _ = spearmanr(av, bv)
    if not np.isnan(r):
        line_records.append({
            "model_id":          line,
            "cross_protein_rho": float(r),
            "n_proteins_used":   int(mask.sum()),
        })

line_df  = pd.DataFrame(line_records)
line_arr = line_df["cross_protein_rho"].to_numpy()
n_trusted = int((line_arr >= 0.3).sum())
n_suspect = int((line_arr <  0.3).sum())
print(f"  lines assessed (>= 50 shared proteins) : {len(line_arr):,}")
print(f"  median cross-protein rho               : {np.median(line_arr):+.4f}")
for q in [10, 25, 50, 75, 90]:
    print(f"    p{q:<2d} : {np.percentile(line_arr, q):+.4f}")
print(f"  lines rho >= 0.3 (trusted 2-layer) : {n_trusted:,}  ({n_trusted/len(line_arr)*100:.1f}%)")
print(f"  lines rho <  0.3 (suspect 2-layer) : {n_suspect:,}  ({n_suspect/len(line_arr)*100:.1f}%)")
line_out = OUTPUTS / "cell_line_platform_agreement.parquet"
line_df.to_parquet(line_out, index=False)
print(f"  Saved: {line_out}")

print()
print("=" * 72)
print("STEP 6 -- abundance-stratified per-protein rho")
print("=" * 72)
print("  Tests whether the 0.373 median is pulled down by low-abundance proteins")
print("  that DIA-MS detects but TMT does not reliably quantify.")
print()

prot_median_abund = pc_vals[shared_prot].median(axis=0)
quartile_edges = np.nanpercentile(prot_median_abund.values, [0, 25, 50, 75, 100])
quartile_names = ["Q1_low", "Q2", "Q3", "Q4_high"]
strat_rows = []
for i, (lo, hi, qname) in enumerate(zip(quartile_edges[:-1], quartile_edges[1:], quartile_names)):
    if i < 3:
        in_q = prot_median_abund[(prot_median_abund >= lo) & (prot_median_abund < hi)].index
    else:
        in_q = prot_median_abund[(prot_median_abund >= lo) & (prot_median_abund <= hi)].index
    q_rhos = np.array([prot_to_rho[p] for p in in_q if p in prot_to_rho])
    if len(q_rhos) == 0:
        continue
    med_q   = float(np.median(q_rhos))
    pct_con = float((q_rhos >= 0.5).mean() * 100)
    strat_rows.append({
        "quartile":          qname,
        "abundance_range":   f"{lo:.2f} – {hi:.2f}",
        "n_proteins":        int(len(q_rhos)),
        "median_rho":        round(med_q, 4),
        "pct_consistent":    round(pct_con, 1),
    })
    print(f"  {qname:10s}  abundance {lo:+.2f}–{hi:+.2f} : "
          f"n={len(q_rhos):,}  median_rho={med_q:+.4f}  pct_consistent={pct_con:.1f}%")

results = {
    "question": "Do ProCan (DIA-MS) and CCLE/Gygi (TMT) proteomics agree, i.e. is a merge safe?",
    "read_only": True,
    "coverage": {
        "procan_lines":      int(pc_vals.shape[0]),
        "procan_proteins":   int(pc_vals.shape[1]),
        "ccle_lines":        int(ccle.shape[0]),
        "ccle_proteins":     int(ccle.shape[1]),
        "shared_lines":      len(shared_lines),
        "shared_proteins":   len(shared_prot),
        "proteins_assessed": int(len(rhos)),
    },
    "per_protein_rho": {
        "median":         round(med, 6),
        "mean":           round(float(rho_arr.mean()), 6),
        "p10":            round(float(np.percentile(rho_arr, 10)), 6),
        "p25":            round(float(np.percentile(rho_arr, 25)), 6),
        "p75":            round(float(np.percentile(rho_arr, 75)), 6),
        "p90":            round(float(np.percentile(rho_arr, 90)), 6),
        "pct_above_0.5":  round(float((rho_arr > 0.5).mean() * 100), 2),
        "pct_above_0.3":  round(float((rho_arr > 0.3).mean() * 100), 2),
        "pct_negative":   round(float((rho_arr < 0).mean()   * 100), 2),
    },
    "platform_conflict_tiers": {
        "n_consistent":   n_consistent,
        "n_cautious":     n_cautious,
        "n_conflicting":  n_conflicting,
        "pct_consistent": round(n_consistent  / len(rhos) * 100, 1),
        "pct_cautious":   round(n_cautious    / len(rhos) * 100, 1),
        "pct_conflicting":round(n_conflicting / len(rhos) * 100, 1),
        "sidecar_file":   str(tier_out),
    },
    "per_line_agreement": {
        "n_lines_assessed":    len(line_arr),
        "median_cross_prot_rho": round(float(np.median(line_arr)), 4),
        "n_trusted_lines":     n_trusted,
        "n_suspect_lines":     n_suspect,
        "pct_trusted":         round(n_trusted / len(line_arr) * 100, 1),
        "sidecar_file":        str(line_out),
    },
    "abundance_stratified_rho": strat_rows,
    "recommended_action": action,
    "verdict": verdict,
}
out = OUTPUTS / "test_run_proteomics_platform_overlap_results.json"
with open(out, "w") as f:
    json.dump(results, f, indent=2)
print()
print(f"Saved: {out}")
print(f"Saved: {tier_out}")
print(f"Saved: {line_out}")
print("No production file was modified.")
