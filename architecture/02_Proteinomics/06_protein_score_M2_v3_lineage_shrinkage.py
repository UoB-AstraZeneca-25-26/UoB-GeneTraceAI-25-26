"""
02_Proteinomics/06_protein_score_M2_v3_lineage_shrinkage.py  (TEST HARNESS --
not live, not wired into run_all.py, does not modify or write to the live DB
or bulk_prot_z.parquet.)
-----------------------------------------------------------------------------
Tests partial-pooling shrinkage in zscore_platform()'s per-lineage robust z,
mirroring transcriptomics.py's method_z() design as closely as this
function's structure allows. This is NOT the same mechanism as the
un-promoted 06_protein_score_M2_v2.py (which added shrinkage to the platform
*combination* step for a different reason -- extreme-tail z's). This targets
the lineage-level centering/scaling step instead.

Live v1 (06_protein_score_M2.py) confirmed as the file actually wired into
run_all.py via run_protein_m2_stage.py -- this variant is built against v1's
logic, imported and monkeypatched (only zscore_platform/robust_z_matrix are
replaced; everything else -- isoform collapse, platform agreement, symmetric
combination -- runs the live v1 code unmodified, called through the live
module object). No DB writes happen anywhere in this script.

Shrinkage design:
  scale = sqrt(w * lineage_scale^2 + (1-w) * global_scale^2),  w = n/(n+n0)
    -- "scale" here is robust_z_matrix's own scale (1.4826*MAD, falling back
    to SD when MAD==0/non-finite), computed the same way at both the
    lineage level (per lineage group, matching the live function) and the
    global level (once, panel-wide, this platform/source only, BEFORE the
    per-lineage loop -- computed fresh here, not reusing any other file's
    global stats, since this is source-specific).
  n = n_valid: the per-gene, per-lineage count of non-missing values already
    computed inside robust_z_matrix (finer-grained than transcriptomics.py's
    per-LINEAGE nb, since a lineage can have many cell lines but few with
    non-missing protein data for a specific gene -- using n_valid is the
    natural quantity for this function's existing structure).
  n0 = 25: reused from RNA_STRATUM_N0 / PROT_RESID_N0 by convention, same as
    those two -- NOT independently derived or validated for this context.
    Flagged, same honesty standard as the other two.

Center (median) is NOT blended -- decision, not an oversight. method_z()
mirrors: ctr = lineage_median if nb >= MIN_PEERS else gene_global_median (a
hard cutoff, not method_z's own continuous center blend -- that blend lives
in core_score.py's SEPARATE RNA-stratum step, a different function). Applying
that same hard cutoff here is a structural no-op: this function already NaNs
every cell with n_valid < MIN_PEERS via the existing min_peers gate (kept
unchanged, per the task), so by the time a cell survives to be assigned a
center, n_valid >= MIN_PEERS is already guaranteed -- the "fall back to
global median" branch of method_z's own logic would never fire here. Mirrored
precisely; the different-looking result (center never blended) is a
consequence of protein's pre-existing NaN gate, not a deviation from RNA's
design.
"""
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent
DB_PATH = PROJECT_ROOT / "outputs" / "celllineselector.db"
THRESHOLDS_PATH = HERE / "results" / "m2_thresholds.json"

sys.path.insert(0, str(HERE))
import importlib
m2 = importlib.import_module("06_protein_score_M2")

N0 = 25          # reused shrinkage pseudo-count; unvalidated for this context, see docstring
MIN_PEERS = m2.MIN_PEERS  # kept identical to live (5) -- shrinkage is additive to this gate, not a replacement


def robust_z_matrix_shrunk(X: np.ndarray, global_med: np.ndarray, global_scale: np.ndarray,
                           min_peers: int = MIN_PEERS, n0: int = N0) -> np.ndarray:
    """Same as live robust_z_matrix(), plus: scale is shrunk toward the
    platform-global (panel-wide) scale, weighted by w=n_valid/(n_valid+n0).
    Center stays the lineage's own median for every surviving cell -- see
    module docstring for why that mirrors method_z()'s design rather than
    deviating from it."""
    X = np.asarray(X, dtype=float)
    if X.size == 0:
        return X
    n_valid = np.sum(np.isfinite(X), axis=0)
    with np.errstate(all="ignore"):
        med_L = np.nanmedian(X, axis=0, keepdims=True)
        mad_L = np.nanmedian(np.abs(X - med_L), axis=0, keepdims=True)
        std_L = np.nanstd(X, axis=0, keepdims=True)
    scale_L = np.where(np.isclose(mad_L, 0.0) | ~np.isfinite(mad_L), std_L, 1.4826 * mad_L)
    scale_L = np.where(np.isclose(scale_L, 0.0) | ~np.isfinite(scale_L), 1.0, scale_L)

    w = (n_valid / (n_valid + n0))[None, :].astype(float)
    scale = np.sqrt(w * scale_L ** 2 + (1.0 - w) * global_scale ** 2)

    Z = (X - med_L) / scale
    Z[:, n_valid < min_peers] = np.nan
    return Z


def zscore_platform_shrunk(df, lineage_map, source_name, n0: int = N0) -> pd.DataFrame:
    """Same structure as live zscore_platform(), with a platform-global
    (panel-wide, this source only) median/scale computed once before the
    per-lineage loop, then blended in via robust_z_matrix_shrunk()."""
    wide = df.pivot_table(index="model_id", columns="gene_id", values="value", aggfunc="first")
    wide["lineage"] = lineage_map.reindex(wide.index)
    wide = wide.dropna(subset=["lineage"])
    empty = pd.DataFrame(columns=["gene_id", "model_id", "lineage", "z_score", "source"])
    if wide.empty:
        return empty

    gene_cols = [c for c in wide.columns if c != "lineage"]
    if not gene_cols:
        return empty

    Xg = wide[gene_cols].to_numpy(dtype=float)
    with np.errstate(all="ignore"):
        global_med = np.nanmedian(Xg, axis=0, keepdims=True)
        global_mad = np.nanmedian(np.abs(Xg - global_med), axis=0, keepdims=True)
        global_std = np.nanstd(Xg, axis=0, keepdims=True)
    global_scale = np.where(np.isclose(global_mad, 0.0) | ~np.isfinite(global_mad),
                            global_std, 1.4826 * global_mad)
    global_scale = np.where(np.isclose(global_scale, 0.0) | ~np.isfinite(global_scale),
                            1.0, global_scale)

    rows = []
    for lineage, grp in wide.groupby("lineage", sort=True):
        gcols = [c for c in grp.columns if c != "lineage"]
        if not gcols:
            continue
        Z = robust_z_matrix_shrunk(grp[gcols].to_numpy(dtype=float), global_med, global_scale)
        zdf = pd.DataFrame(Z, index=grp.index, columns=gcols).stack().reset_index()
        zdf.columns = ["model_id", "gene_id", "z_score"]
        zdf["lineage"] = lineage
        zdf["source"] = source_name
        rows.append(zdf)

    if not rows:
        return empty
    out = pd.concat(rows, ignore_index=True)
    return out[["gene_id", "model_id", "lineage", "z_score", "source"]]


def run():
    import json
    thresholds = json.loads(THRESHOLDS_PATH.read_text())

    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        print("=" * 70)
        print("Protein scorer v3 -- lineage-level partial-pooling shrinkage test")
        print("(isolated: live v1 module, zscore_platform monkeypatched, no DB writes)")
        print("=" * 70)

        # LIVE v1 result, unmodified, for direct comparison
        result_live, _, _ = m2.compute_protein_z(con, thresholds=thresholds)

        # SHRUNK: monkeypatch only zscore_platform (compute_from_sources looks
        # it up via the module's own global namespace at call time)
        m2.zscore_platform = zscore_platform_shrunk
        result_shrunk, _, _ = m2.compute_protein_z(con, thresholds=thresholds)
    finally:
        con.close()

    print(f"\nlive:   {len(result_live):,} rows  genes={result_live.gene_id.nunique():,}  "
          f"lines={result_live.model_id.nunique():,}")
    print(f"shrunk: {len(result_shrunk):,} rows  genes={result_shrunk.gene_id.nunique():,}  "
          f"lines={result_shrunk.model_id.nunique():,}")

    result_live.to_parquet(HERE / "results" / "protein_z_v1_live_snapshot.parquet", index=False)
    result_shrunk.to_parquet(HERE / "results" / "protein_z_v3_lineage_shrunk.parquet", index=False)
    print(f"\nSaved -> results/protein_z_v1_live_snapshot.parquet, "
          f"results/protein_z_v3_lineage_shrunk.parquet")


if __name__ == "__main__":
    run()
