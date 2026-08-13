"""
build_cn_layer_complete.py
--------------------------
Gene-complete copy-number layer, built from the DepMap 20Q2 release matrix that
had been sitting unused in data/DepMap_Chronos/.

WHY THIS EXISTS
The shipped CNA layer (build_cna_layer.py -> cna_flags.parquet) is derived from
COSMIC CellLinesProject CNA, which covers 997 of 2,145 panel models and only
0.31% of scored (gene, line) pairs. Measured in docs/RANKING_FEASIBILITY.md, the
resulting cna_gate is NOT ASSESSED for 99.72% of gate-retained pairs -- it is
inert, and it was the binding constraint on building any confirmation count.

This file does not modify or replace build_cna_layer.py. It is additive: a
second, denser CN source written to its own artefact, so the two can be compared
before anything downstream is repointed.

THE SOURCE
data/DepMap_Chronos/CCLEGeneCopyNumber20Q2.hdf5
    /data   908 cell lines x 27,639 genes, float64, ZERO missing values
    /dim_0  'ACH-000001' style Broad ids
    /dim_1  'SYMBOL (ENTREZ)' style gene labels

SCALE -- established empirically, not assumed:
    values are log2(relative_CN + 1), neutral = 1.0
      median 1.0019  =>  implied relative CN 1.0026
      p5 0.6816, p95 1.3277
    RELATIVE means relative to each cell line's OWN ploidy. Confirmed by
    Spearman(per-line mean value, GDSC ploidy_wgs/wes) = +0.115 with per-line
    mean sd = 0.020 -- i.e. essentially constant across a panel whose ploidy
    ranges 1.52 to 5.40. A raw absolute-copy matrix would track ploidy strongly.

    CONSEQUENCE: no external ploidy correction is applied, and the "30% of the
    panel has no ploidy estimate" limitation (MATH_REFERENCE.md gap 8) does not
    bite for this layer. That gap is closed here for CN calls, though it remains
    open for the COSMIC-derived thresholds in build_cna_layer.py.

THRESHOLDS
    Fixed a priori in docs/RANKING_PRESPEC.md before any downstream measurement,
    and deliberately NOT selected on any outcome:
        deletion       relative CN < 0.5     (primary; ~one-copy loss from diploid)
        deletion       relative CN < 0.25    (sensitivity, deep loss)
        deletion       relative CN < 0.75    (sensitivity, shallow)
        amplification  relative CN > 2.0     (reported, not gated on)

OUTPUT
    src/pipeline/outputs/cn_gene_complete.parquet
        wide, index = model_id (lowercase), columns = ensg_id, values = relative CN
    src/pipeline/outputs/cn_gene_complete_meta.json
        provenance, scale evidence, thresholds, coverage against the panel
"""
import json
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

SRC = Path("data/DepMap_Chronos/CCLEGeneCopyNumber20Q2.hdf5")
REF = Path("reference")
OUTPUTS = Path("src/pipeline/outputs")

DEL_PRIMARY = 0.5
DEL_SENSITIVITY = (0.25, 0.75)
AMP_THRESHOLD = 2.0

t0 = time.time()

with h5py.File(SRC, "r") as fh:
    D = fh["data"][:]
    lines_raw = [x.decode() for x in fh["dim_0"][:]]
    genes_raw = [x.decode() for x in fh["dim_1"][:]]
print(f"[{time.time()-t0:.1f}s] source: {D.shape[0]:,} lines x {D.shape[1]:,} genes, "
      f"{100*np.isnan(D).mean():.2f}% missing")

# ---- scale evidence, recomputed every run rather than trusted -----------------
flat = D[np.isfinite(D)]
scale_evidence = {
    "median_value": float(np.median(flat)),
    "p5": float(np.percentile(flat, 5)),
    "p95": float(np.percentile(flat, 95)),
    "implied_median_relative_cn": float(2 ** np.median(flat) - 1),
    "interpretation": "log2(relative_CN + 1), neutral = 1.0",
}
assert 0.9 < scale_evidence["implied_median_relative_cn"] < 1.1, (
    "median relative CN is not ~1.0 -- the scale assumption is wrong, refusing to build")
print(f"[{time.time()-t0:.1f}s] scale check passed: implied median relative CN = "
      f"{scale_evidence['implied_median_relative_cn']:.4f}")

# ---- gene symbols -> ensg -----------------------------------------------------
gl = pd.read_parquet(REF / "gene_lookup.parquet", columns=["ensg_id", "hgnc_symbol"])
gl["ensg_id"] = gl.ensg_id.astype(str).str.split(".").str[0].str.lower()
sym2ensg = dict(zip(gl.hgnc_symbol.astype(str).str.upper(), gl.ensg_id))

keep_idx, keep_ensg = [], []
for i, label in enumerate(genes_raw):
    e = sym2ensg.get(label.split(" (")[0].strip().upper())
    if e:
        keep_idx.append(i)
        keep_ensg.append(e)
print(f"[{time.time()-t0:.1f}s] {len(keep_idx):,}/{len(genes_raw):,} gene labels mapped to ENSG")

W = pd.DataFrame(D[:, keep_idx], index=[m.lower() for m in lines_raw], columns=keep_ensg)
# duplicate symbol -> ensg collapses: mean of the copy-number values
n_dup = len(keep_ensg) - len(set(keep_ensg))
W = W.T.groupby(level=0).mean().T
print(f"[{time.time()-t0:.1f}s] collapsed {n_dup} duplicate ENSG columns -> "
      f"{W.shape[1]:,} unique genes")

# ---- log2(rel+1) -> relative CN ----------------------------------------------
REL = np.power(2.0, W) - 1.0
REL = REL.clip(lower=0.0).astype("float32")

# ---- coverage against the panel ----------------------------------------------
cov = pd.read_parquet(OUTPUTS / "coverage_matrix_enriched.parquet", columns=["model_id"])
panel = set(cov.model_id)
covered = set(REL.index) & panel
print(f"[{time.time()-t0:.1f}s] panel coverage: {len(covered):,} of {len(panel):,} models "
      f"({100*len(covered)/len(panel):.1f}%)")

old = pd.read_parquet(OUTPUTS / "cna_flags.parquet", columns=["model_id", "ensg_id"])
print(f"           existing COSMIC cna_flags: {old.model_id.nunique():,} models, "
      f"{old.ensg_id.nunique():,} genes, {len(old):,} pairs")
print(f"           this layer:                {REL.shape[0]:,} models, "
      f"{REL.shape[1]:,} genes, {REL.size:,} pairs")

rates = {}
for thr in (DEL_PRIMARY, *DEL_SENSITIVITY):
    rates[f"deletion_rate_at_rel_lt_{thr}"] = float((REL < thr).to_numpy().mean())
rates[f"amplification_rate_at_rel_gt_{AMP_THRESHOLD}"] = float(
    (REL > AMP_THRESHOLD).to_numpy().mean())
for k, v in rates.items():
    print(f"           {k}: {100*v:.2f}%")

OUTPUTS.mkdir(parents=True, exist_ok=True)
REL.to_parquet(OUTPUTS / "cn_gene_complete.parquet")

meta = {
    "source_file": str(SRC),
    "source_shape": {"lines": int(D.shape[0]), "genes": int(D.shape[1])},
    "missing_fraction_in_source": float(np.isnan(D).mean()),
    "scale_evidence": scale_evidence,
    "ploidy_normalised": True,
    "ploidy_evidence": ("Spearman(per-line mean value, GDSC ploidy) = +0.115 with "
                        "per-line mean sd 0.020 over a panel of ploidy 1.52-5.40; "
                        "no external ploidy correction applied"),
    "genes_mapped_to_ensg": int(len(keep_idx)),
    "genes_in_source": int(len(genes_raw)),
    "unique_ensg_after_collapse": int(REL.shape[1]),
    "duplicate_ensg_collapsed": int(n_dup),
    "panel_models_covered": int(len(covered)),
    "panel_models_total": int(len(panel)),
    "thresholds": {"deletion_primary_relative_cn": DEL_PRIMARY,
                   "deletion_sensitivity": list(DEL_SENSITIVITY),
                   "amplification_relative_cn": AMP_THRESHOLD,
                   "fixed_where": "docs/RANKING_PRESPEC.md, a priori, not outcome-selected"},
    "call_rates": rates,
    "values_stored": "relative copy number (linear), float32",
    "supersedes": None,
    "note": ("Additive. Does not modify build_cna_layer.py or cna_flags.parquet. "
             "Compare before repointing anything downstream."),
}
with open(OUTPUTS / "cn_gene_complete_meta.json", "w") as fh:
    json.dump(meta, fh, indent=2)

print(f"\n[{time.time()-t0:.1f}s] Written: {OUTPUTS/'cn_gene_complete.parquet'} "
      f"({REL.shape[0]:,} x {REL.shape[1]:,})")
print(f"[{time.time()-t0:.1f}s] Written: {OUTPUTS/'cn_gene_complete_meta.json'}")
