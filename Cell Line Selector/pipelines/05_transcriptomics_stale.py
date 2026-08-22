"""
transcriptomics.py

Integrates depmap_expr, hpa_rna, and geo_expr onto a common log2
scale, produces publication-quality cluster density figures, and
computes robust z-scores for every protein-coding gene x cell line
into a queryable DuckDB store.

GENE SET SCOPE -- two distinct lists, deliberately:
  all_genes     -- EVERY protein-coding gene in gene_roster. Used for
                   the integrated matrix and the z-score database. A
                   gene measured by only one source is still validly
                   scored against that source's own strata; it simply
                   has fewer sources in its Stouffer combination.
  cluster_genes -- the 3-way intersection. Used ONLY for embedding
                   and density figures, which need one shared
                   coordinate space and therefore every gene present
                   in every source.

Z-SCORE DESIGN:
  Baselines are built within TECHNICAL strata (source | method),
  never within biological ones. Grouping by lineage to build a
  baseline would define away the signal of interest -- a gene high
  across all breast lines would score ~0 for every breast line.
  Lineage enters as a reporting dimension and a secondary diagnostic.

  Sources combine by STOUFFER'S METHOD. DepMap, HPA and GEO are
  independent experiments run by different institutions, so their
  measurement errors are unrelated and cross-institution agreement is
  genuine corroboration. Plain averaging (kept as z_mean) shrinks
  toward zero as sources are added, understating that agreement.

  The MAD scale factor is DERIVED from the normal quantile function
  (see mad_consistency_constant) and empirically calibrated against
  the data (estimate_empirical_mad_scale), not hardcoded.

ASSUMPTIONS NEEDING VERIFICATION AGAINST YOUR SCHEMA:
  - hpa_rna's expression VALUE column name (auto-detected)
  - GEO methods other than RMA/gcRMA treated as LINEAR and
    log2-transformed -- verify with summarize_geo_scale_by_method
  - depmap and hpa are RNA-seq; GEO's RMA/gcRMA/MAS5 are Affymetrix
    microarray pipelines
  - standardize_per_source is a lightweight partial batch correction,
    not a substitute for ComBat

Dependencies: GEOparse, scikit-learn, scipy, duckdb; optionally umap-learn.
    pip install GEOparse scikit-learn scipy duckdb umap-learn
"""

import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

try:
    from sklearn.decomposition import PCA  # type: ignore[import-not-found]
    from sklearn.preprocessing import StandardScaler  # type: ignore[import-not-found]
    from sklearn.cluster import KMeans  # type: ignore[import-not-found]
    from sklearn.metrics import adjusted_rand_score  # type: ignore[import-not-found]
except ImportError:
    PCA = None
    StandardScaler = None
    KMeans = None
    adjusted_rand_score = None

try:
    import umap  # type: ignore
except ImportError:
    umap = None

try:
    from scipy.stats import gaussian_kde, norm  # type: ignore[import-not-found]
except ImportError:
    gaussian_kde = None
    norm = None

try:
    import duckdb  # type: ignore[import-not-found]
except ImportError:
    duckdb = None

from src.scripts.eda_functions import style_axis, _safe_savefig, cols_of, numeric_summary, cached_csv
from src.scripts.logging_utils import get_logger, log_error

try:
    import GEOparse  # type: ignore[import-not-found]
except ImportError:
    GEOparse = None

logger = get_logger("transcriptomics")

MAX_WORKERS = min(8, os.cpu_count() or 4)


# ============================================================
# 0. Figure style (thesis presentation)
# ============================================================

# Colourblind-safe qualitative palette (Okabe-Ito), chosen over tab10
# because it stays distinguishable in greyscale print and under the
# common forms of colour vision deficiency -- both matter for a
# printed thesis.
OKABE_ITO = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#F0E442", "#000000",
]

# Low-contrast sequential ramp for density contours, so the scatter
# layer on top stays readable.
DENSITY_CMAP = LinearSegmentedColormap.from_list(
    "density_blue", ["#FFFFFF", "#E8EEF4", "#C3D4E3", "#8FAFC9", "#5B84A8"]
)


def apply_thesis_style() -> None:
    """Restrained typography, hairline spines, print-grade DPI."""
    plt.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 9,
        "axes.titlesize": 10.5,
        "axes.titleweight": "semibold",
        "axes.titlepad": 10,
        "axes.labelsize": 9,
        "axes.labelcolor": "#333333",
        "axes.edgecolor": "#BBBBBB",
        "axes.linewidth": 0.8,
        "axes.grid": True,
        "grid.color": "#E8E8E8",
        "grid.linewidth": 0.6,
        "xtick.color": "#555555",
        "ytick.color": "#555555",
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "legend.title_fontsize": 8.5,
    })


apply_thesis_style()


def _clean_axis(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#BBBBBB")
        ax.spines[side].set_linewidth(0.8)
    ax.grid(True, alpha=0.5, zorder=0)
    ax.set_axisbelow(True)


# ============================================================
# 1. Caching helpers
# ============================================================

def cached_parquet(path: Path, compute_fn, force: bool = False) -> pd.DataFrame:
    path = Path(path)
    if not force and path.exists():
        return pd.read_parquet(path)
    df = compute_fn()
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    return df


def cached_gene_list(path: Path, compute_fn, force: bool = False) -> list:
    path = Path(path)
    if not force and path.exists():
        return pd.read_csv(path)["gene_id"].tolist()
    genes = compute_fn()
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene_id": genes}).to_csv(path, index=False)
    return genes


# ============================================================
# 2. GEO per-GSM method metadata
# ============================================================

def harvest_geo_methods(gsms: list, destdir: str = "./geo_meta", verbose: bool = True) -> pd.DataFrame:
    """
    Fetch per-GSM metadata from NCBI and classify the normalisation
    method. GEO archives; it does not compute -- each GSM carries
    whatever its submitter uploaded, so the method must be read per
    sample rather than assumed for a whole series.

    Checks destdir for an already-downloaded file per GSM first, so
    reruns only hit the network for genuinely new accessions.
    """
    if GEOparse is None:
        raise ImportError("GEOparse is required. Install it with: pip install GEOparse")

    import glob
    os.makedirs(destdir, exist_ok=True)

    rows = []
    n_cached = 0
    for i, g in enumerate(gsms, 1):
        try:
            existing = glob.glob(os.path.join(destdir, f"{g.upper()}*"))
            if existing:
                s = GEOparse.get_GEO(filepath=existing[0], silent=True)
                n_cached += 1
            else:
                s = GEOparse.get_GEO(geo=g.upper(), destdir=destdir, silent=True)
            m = s.metadata
            rows.append({
                "gsm": g.lower(),
                "platform": m.get("platform_id", [None])[0],
                "processing": (m.get("data_processing", [""])[0] or "")[:200],
                "value_def": next((v for k, v in s.columns.description.items()
                                    if k.upper() == "VALUE"), None),
                "has_abs_call": "ABS_CALL" in s.table.columns,
                "series": "; ".join(m.get("series_id", [])),
                "title": (m.get("title", [""])[0] or "")[:120],
            })
        except Exception as e:
            rows.append({"gsm": g.lower(), "platform": None,
                         "processing": f"FETCH FAILED: {e}"})
        if verbose and i % 200 == 0:
            print(f"  {i:,} / {len(gsms):,}  ({n_cached:,} from local cache)")

    meta = pd.DataFrame(rows)

    # gcRMA MUST be tested before RMA -- the string "gcRMA" contains "RMA"
    meta["method"] = np.select(
        [meta.processing.str.contains("gcrma", case=False, na=False),
         meta.processing.str.contains("rma|robust multi", case=False, na=False),
         meta.processing.str.contains("mas5|mas 5|gcos", case=False, na=False),
         meta.processing.str.contains("quantile", case=False, na=False)],
        ["gcRMA", "RMA", "MAS5", "quantile"], default="other")

    if verbose:
        print(f"\n{len(meta):,} GSMs ({n_cached:,} loaded locally)\n")
        print(meta.platform.value_counts().to_string())
        print()
        print(meta.method.value_counts().to_string())
        failed = meta.processing.str.startswith("FETCH FAILED").sum()
        if failed:
            print(f"FAILED to fetch: {failed:,}")
    return meta


def load_or_harvest_geo_methods(con, cache_path, force: bool = False) -> pd.Series:
    cache_path = Path(cache_path)
    if not force and cache_path.exists():
        meta = pd.read_csv(cache_path)
        print(f"  loaded cached GEO metadata: {len(meta):,} GSMs")
    else:
        gsms = con.execute(
            "SELECT DISTINCT geo_accession FROM geo_info WHERE geo_accession IS NOT NULL"
        ).df()["geo_accession"].tolist()
        meta = harvest_geo_methods(gsms)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        meta.to_csv(cache_path, index=False)
        print(f"  harvested and cached {len(meta):,} GSMs")

    gsm_method = meta.set_index(meta["gsm"].str.lower())["method"]
    print(gsm_method.value_counts().to_string())
    return gsm_method


def register_gsm_method(con, gsm_method: pd.Series) -> None:
    df = gsm_method.rename("method").rename_axis("gsm").reset_index()
    con.register("gsm_method_lookup", df)


# ============================================================
# 3. Scale verification
# ============================================================

def row_wise_scale_stats(con, table: str, id_col: str) -> pd.DataFrame:
    """
    Per-ROW min/mean/max via DuckDB UNPIVOT + GROUP BY -- one query
    regardless of column count. Scale verification is inherently
    per-sample here, since GEO's scale varies sample to sample.
    """
    value_cols = [c for c in cols_of(con, table) if c != id_col]
    if not value_cols:
        return pd.DataFrame(columns=[id_col, "row_min", "row_mean", "row_max"])

    query = f'''
        SELECT "{id_col}" AS id, min(value) AS row_min, avg(value) AS row_mean, max(value) AS row_max
        FROM (UNPIVOT "{table}" ON COLUMNS(* EXCLUDE ("{id_col}")) INTO NAME gene VALUE value)
        WHERE value IS NOT NULL
        GROUP BY "{id_col}"
    '''
    return con.execute(query).fetchdf().rename(columns={"id": id_col})


def summarize_source_scale(con, table: str, id_col: str, label: str) -> dict:
    stats = row_wise_scale_stats(con, table, id_col)
    if stats.empty:
        return {"source": label, "n_rows": 0}
    return {
        "source": label,
        "n_rows": len(stats),
        "median_row_mean": stats["row_mean"].median(),
        "median_row_max": stats["row_max"].median(),
        "overall_min": stats["row_min"].min(),
        "overall_max": stats["row_max"].max(),
    }


def summarize_geo_scale_by_method(con, gsm_method: pd.Series, table: str = "geo_expr",
                                   id_col: str = "gsm_id") -> pd.DataFrame:
    """
    Scale fingerprint PER GEO METHOD -- the check that matters, since
    GEO scale is submitter-defined. RMA/gcRMA should show small
    maxima (log2 range); MAS5/other much larger (linear). A group
    that doesn't match invalidates the log2/linear split for it.
    """
    stats = row_wise_scale_stats(con, table, id_col)
    stats["method"] = stats[id_col].str.lower().map(gsm_method).fillna("unknown")

    return stats.groupby("method").agg(
        n_samples=("row_max", "size"),
        median_row_mean=("row_mean", "median"),
        median_row_max=("row_max", "median"),
        min_row_max=("row_max", "min"),
        max_row_max=("row_max", "max"),
    ).reset_index()


# ============================================================
# 4. Scale harmonisation -- everything to log2(x+1)
# ============================================================

GEO_ALREADY_LOG_METHODS = {"RMA", "gcRMA"}
_HPA_VALUE_COL_CANDIDATES = ["ntpm", "ptpm", "tpm", "value", "expression"]


def _find_value_col(available_cols: list, candidates: list) -> Optional[str]:
    return next((c for c in candidates if c in available_cols), None)


def harmonized_hpa_query(con, table: str = "hpa_rna", value_col: Optional[str] = None) -> str:
    available = cols_of(con, table)
    value_col = value_col or _find_value_col(available, _HPA_VALUE_COL_CANDIDATES)
    if value_col is None:
        raise ValueError(
            f"Couldn't auto-detect hpa_rna's expression value column. "
            f"Available columns: {available}. Pass value_col explicitly."
        )
    return f'SELECT *, log2("{value_col}" + 1) AS log2_value FROM "{table}"'


def harmonized_geo_long_query(table: str = "geo_expr", id_col: str = "gsm_id") -> str:
    """geo_expr in LONG format with the per-method log2 transform applied."""
    already_log_list = ", ".join(f"'{m}'" for m in GEO_ALREADY_LOG_METHODS)
    return f'''
        WITH long AS (
            UNPIVOT "{table}" ON COLUMNS(* EXCLUDE ("{id_col}")) INTO NAME gene VALUE value
        )
        SELECT l."{id_col}", l.gene, l.value,
               CASE WHEN m.method IN ({already_log_list}) THEN l.value
                    ELSE log2(greatest(l.value, 0) + 1) END AS log2_value
        FROM long l
        LEFT JOIN gsm_method_lookup m ON lower(l."{id_col}") = lower(m.gsm)
        WHERE l.value IS NOT NULL
    '''


# ============================================================
# 5. Cell-line lineage
# ============================================================

def build_all_cell_lines_with_lineage(tables: dict, cell_line_roster: pd.DataFrame) -> pd.DataFrame:
    """
    Every cell line from cell_line_roster, joined with lineage from
    sample_info's sample_collection_site. A model_id with conflicting
    lineage values keeps all of them semicolon-joined, so the
    conflict stays visible rather than being silently resolved.
    """
    sample_info = tables.get("sample_info")
    if sample_info is None:
        raise ValueError("sample_info table not found in `tables`")
    if "sample_collection_site" not in sample_info.columns:
        raise ValueError(
            f"Expected 'sample_collection_site' in sample_info. "
            f"Available: {sample_info.columns.tolist()}"
        )
    if "model_id" not in sample_info.columns:
        raise ValueError("sample_info has no model_id column")

    lineage_lookup = (
        sample_info[["model_id", "sample_collection_site"]]
        .rename(columns={"sample_collection_site": "lineage"})
        .dropna(subset=["model_id", "lineage"])
        .drop_duplicates()
        .groupby("model_id")["lineage"]
        .apply(lambda s: "; ".join(sorted(set(s))))
        .reset_index()
    )

    result = cell_line_roster.merge(lineage_lookup, on="model_id", how="left")

    if "cell_line_name(s)" in result.columns:
        result["cell_line_name"] = (
            result["cell_line_name(s)"].fillna("").str.split("; ").str[0].replace("", np.nan)
        )
    return result


# ============================================================
# 6. Gene selection -- two distinct lists
# ============================================================

def select_all_protein_coding_genes(gene_roster: pd.DataFrame) -> list:
    """
    EVERY protein-coding gene_id in gene_roster (already filtered to
    protein_coding upstream). This is the scoring gene set -- z-scores
    do NOT require presence in all three sources, since each gene is
    scored independently within its own source|method strata.
    """
    genes = sorted(set(gene_roster["gene_id"].dropna()))
    print(f"  {len(genes):,} protein-coding genes in gene_roster")
    return genes


def select_common_genes(con, gene_roster: pd.DataFrame, top_n: Optional[int] = None) -> list:
    """
    The 3-way intersection: genes covered in depmap_expr AND geo_expr
    (via gene_roster's coverage columns) AND hpa_rna (queried
    directly). Used ONLY for clustering -- a shared embedding needs
    every gene present in every source.
    """
    depmap_geo = set(gene_roster[
        (gene_roster["profile_id(s)"].fillna("") != "") &
        (gene_roster["geo_gsm_id(s)"].fillna("") != "")
    ]["gene_id"])
    print(f"  covered in depmap_expr AND geo_expr: {len(depmap_geo):,}")

    hpa_covered = set(con.execute(
        "SELECT DISTINCT gene_id FROM hpa_rna WHERE gene_id IS NOT NULL"
    ).fetchdf()["gene_id"])
    print(f"  covered in hpa_rna: {len(hpa_covered):,}")

    covered = sorted(depmap_geo & hpa_covered)
    print(f"  covered in ALL THREE sources: {len(covered):,}")

    if not covered:
        raise ValueError("No protein-coding genes covered in all three sources.")
    if top_n is None:
        return covered

    summary = numeric_summary(con, "depmap_expr", cols=covered)
    top = summary.assign(var=lambda d: d["std"] ** 2).sort_values("var", ascending=False).head(top_n)
    return top["column"].tolist()


# ============================================================
# 7. Integrated matrix -- all samples, all protein-coding genes
# ============================================================

def _fetch_depmap_wide(con, genes: list, chunk: int = 300,
                        max_workers: int = MAX_WORKERS) -> pd.DataFrame:
    """
    Fetch depmap_expr's selected gene columns in batches, run
    CONCURRENTLY (each on its own con.cursor(), safe for parallel use
    on one DuckDB database; DuckDB releases the GIL during execute()).

    Genes absent from depmap_expr are skipped in SQL and added back as
    all-NaN columns -- with the full protein-coding set, many genes
    exist in hpa/geo but not depmap, and SELECTing a missing column
    would otherwise raise.
    """
    available = set(cols_of(con, "depmap_expr"))
    present = [g for g in genes if g in available]
    missing = [g for g in genes if g not in available]
    if missing:
        print(f"    {len(missing):,} genes absent from depmap_expr -- carried as NaN")

    if not present:
        base = con.execute('SELECT profile_id FROM depmap_expr').fetchdf()
        for g in genes:
            base[g] = np.nan
        return base

    batches = [present[i:i + chunk] for i in range(0, len(present), chunk)]

    def _fetch(batch):
        cur = con.cursor()
        cols_sql = ", ".join(f'"{g}"' for g in batch)
        df = cur.execute(f'SELECT profile_id, {cols_sql} FROM depmap_expr').fetchdf()
        cur.close()
        return df.set_index("profile_id")

    frames = []
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(_fetch, b) for b in batches]
        for fut in as_completed(futures):
            frames.append(fut.result())

    combined = pd.concat(frames, axis=1).reset_index()
    for g in missing:
        combined[g] = np.nan
    return combined


def build_integrated_matrix(con, genes: list, gsm_method: pd.Series,
                             all_cell_lines: pd.DataFrame,
                             hpa_value_col: Optional[str] = None) -> pd.DataFrame:
    """
    One integrated matrix, one row per sample. EVERY sample from
    EVERY source is included -- a cell line unique to one source is
    still kept, with model_id/lineage NaN where no match exists.
    """
    gene_list_sql = ", ".join(f"'{g}'" for g in genes)

    print("    fetching depmap...")
    depmap = _fetch_depmap_wide(con, genes)
    depmap_profiles = con.execute('SELECT profile_id, model_id FROM depmap_profiles').fetchdf()
    depmap = depmap.merge(depmap_profiles, on="profile_id", how="left")
    depmap_long = depmap.melt(id_vars=["profile_id", "model_id"],
                               var_name="gene", value_name="log2_value")
    depmap_long["sample_id"] = depmap_long["profile_id"]
    depmap_long["source"] = "depmap"

    print("    fetching geo...")
    register_gsm_method(con, gsm_method)
    geo_query = harmonized_geo_long_query()
    geo_long = con.execute(f'''
        SELECT gsm_id, gene, log2_value FROM ({geo_query}) WHERE gene IN ({gene_list_sql})
    ''').fetchdf()
    geo_info = con.execute(
        'SELECT geo_accession, model_id FROM geo_info WHERE geo_accession IS NOT NULL'
    ).fetchdf()
    geo_long = geo_long.merge(geo_info, left_on="gsm_id", right_on="geo_accession", how="left")
    geo_long["sample_id"] = geo_long["gsm_id"]
    geo_long["source"] = "geo"

    print("    fetching hpa...")
    hpa_query = harmonized_hpa_query(con, value_col=hpa_value_col)
    hpa_long = con.execute(f'''
        SELECT gene_id AS gene, model_id, log2_value FROM ({hpa_query})
        WHERE gene_id IN ({gene_list_sql})
    ''').fetchdf()
    hpa_long["sample_id"] = hpa_long["model_id"]
    hpa_long["source"] = "hpa"

    cols = ["source", "sample_id", "model_id", "gene", "log2_value"]
    combined_long = pd.concat(
        [depmap_long[cols], geo_long[cols], hpa_long[cols]], ignore_index=True)

    print("    pivoting to wide...")
    wide = combined_long.pivot_table(
        index=["source", "sample_id", "model_id"], columns="gene", values="log2_value"
    ).reset_index()
    wide.columns.name = None

    lineage_lookup = (
        all_cell_lines[["model_id", "lineage"]]
        .dropna(subset=["model_id"]).drop_duplicates(subset=["model_id"])
    )
    return wide.merge(lineage_lookup, on="model_id", how="left")


# ============================================================
# 8. Batch correction + full sample retention
# ============================================================

def standardize_per_source(matrix: pd.DataFrame, gene_cols: list) -> pd.DataFrame:
    """
    Z-scores each gene WITHIN each source so no source's raw scale
    dominates the embedding. VISUALISATION ONLY -- the z-score
    pipeline works from the raw matrix, since scoring
    already-standardised values is double normalisation.
    """
    out = matrix.copy()
    for _, idx in matrix.groupby("source").groups.items():
        sub = out.loc[idx, gene_cols]
        std = sub.std().replace(0, np.nan)
        out.loc[idx, gene_cols] = (sub - sub.mean()) / std
    return out


def prepare_for_clustering(matrix: pd.DataFrame, gene_cols: list) -> pd.DataFrame:
    """
    Keeps EVERY sample regardless of missingness. Missing values are
    filled with 0 -- after standardize_per_source that is the
    source's own per-gene mean, a neutral fill.

    Trade-off: a sparse sample is pulled toward its source centroid
    for every gene it lacks, understating how different it truly is.
    That is the accepted cost of dropping nothing.
    """
    out = matrix.copy()
    numeric = out[gene_cols].apply(pd.to_numeric, errors="coerce")
    n_missing = numeric.isna().sum(axis=1)
    if (n_missing > 0).any():
        sparse = n_missing[n_missing > 0]
        print(f"    filled missing values for {len(sparse):,} / {len(out):,} samples "
              f"(median {sparse.median():.0f} of {len(gene_cols):,} genes)")
    out[gene_cols] = numeric.fillna(0)
    return out


# ============================================================
# 9. Sample labelling
# ============================================================

_SOURCE_TECHNOLOGY = {"depmap": "rna_seq", "hpa": "rna_seq"}


def attach_method_label(matrix: pd.DataFrame, gsm_method: pd.Series) -> pd.DataFrame:
    """
    Adds `method`: depmap/hpa -> 'rna_seq'; geo -> that GSM's own
    normalisation method, or 'geo_unknown' if its metadata fetch
    failed.

    Mixed granularity is deliberate: depmap and hpa each run one
    known pipeline so naming it adds nothing, whereas GEO's variation
    is precisely what's worth seeing.
    """
    out = matrix.copy()
    out["method"] = out["source"].map(_SOURCE_TECHNOLOGY)
    is_geo = out["source"] == "geo"
    out.loc[is_geo, "method"] = (
        out.loc[is_geo, "sample_id"].astype(str).str.lower().map(gsm_method).fillna("geo_unknown")
    )
    out["method"] = out["method"].fillna("unknown")
    return out


def attach_stratum_label(matrix: pd.DataFrame, gsm_method: pd.Series) -> pd.DataFrame:
    """
    Adds `stratum` = source | method -- the technical unit within
    which values are directly comparable. depmap and hpa stay
    SEPARATE despite both being RNA-seq: their normalisation
    pipelines differ, so pooling them would leave a technical offset
    inside one stratum.
    """
    out = attach_method_label(matrix, gsm_method)
    out["stratum"] = out["source"].astype(str) + " | " + out["method"].astype(str)
    return out


def common_cell_line_model_ids(matrix: pd.DataFrame) -> set:
    """model_ids present in all three sources."""
    ids_by_source = {
        s: set(matrix.loc[matrix["source"] == s, "model_id"].dropna())
        for s in matrix["source"].unique()
    }
    if len(ids_by_source) < 3:
        return set()
    return set.intersection(*ids_by_source.values())


# ============================================================
# 10. Density figures
# ============================================================

def _reduce_to_2d(numeric: pd.DataFrame, use_umap: bool):
    scaled = StandardScaler().fit_transform(numeric)
    if use_umap and umap is not None:
        return umap.UMAP(n_components=2, random_state=0).fit_transform(scaled), "UMAP"
    return PCA(n_components=2, random_state=0).fit_transform(scaled), "PC"


def _density_panel(ax, embedding: np.ndarray, labels: pd.Series, title: str,
                    axis_label: str, legend_title: str, max_legend: int = 12) -> None:
    """
    Filled KDE contours underneath, categorical scatter on top.
    Groups are drawn largest-first so small groups aren't buried, and
    the legend is capped so a many-category panel stays readable.
    """
    labels = labels.fillna("unknown").astype(str)

    if gaussian_kde is not None and len(embedding) > 20:
        try:
            kde = gaussian_kde(embedding.T)
            pad_x = (embedding[:, 0].max() - embedding[:, 0].min()) * 0.05
            pad_y = (embedding[:, 1].max() - embedding[:, 1].min()) * 0.05
            xx, yy = np.mgrid[
                embedding[:, 0].min() - pad_x:embedding[:, 0].max() + pad_x:180j,
                embedding[:, 1].min() - pad_y:embedding[:, 1].max() + pad_y:180j,
            ]
            zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
            ax.contourf(xx, yy, zz, levels=12, cmap=DENSITY_CMAP, alpha=0.85, zorder=1)
            ax.contour(xx, yy, zz, levels=6, colors="#7A96AD", linewidths=0.4, alpha=0.5, zorder=2)
        except Exception:
            pass  # near-degenerate clouds make KDE singular; scatter alone still informative

    counts = labels.value_counts()
    for i, cat in enumerate(counts.index.tolist()):
        mask = (labels == cat).values
        ax.scatter(embedding[mask, 0], embedding[mask, 1],
                   s=11, alpha=0.75, linewidths=0,
                   color=OKABE_ITO[i % len(OKABE_ITO)],
                   label=f"{cat}  (n={counts[cat]:,})", zorder=3 + i)

    handles, legend_labels = ax.get_legend_handles_labels()
    if len(handles) > max_legend:
        n_hidden = len(counts) - max_legend
        handles, legend_labels = handles[:max_legend], legend_labels[:max_legend]
        legend_labels[-1] += f"  (+{n_hidden} more)"
    ax.legend(handles, legend_labels, title=legend_title, loc="best",
              markerscale=1.8, handletextpad=0.4, labelspacing=0.45)

    ax.set_xlabel(f"{axis_label} 1")
    ax.set_ylabel(f"{axis_label} 2")
    ax.set_title(title)
    _clean_axis(ax)


def plot_source_and_method_density(matrix: pd.DataFrame, gsm_method: pd.Series, gene_cols: list,
                                    out_path: Path, use_umap: bool = True) -> Path:
    """
    Two-panel density figure over ALL cell lines, nothing dropped:

      (a) coloured by SOURCE  -- depmap / geo / hpa
      (b) coloured by METHOD  -- rna_seq / RMA / gcRMA / MAS5 / ...

    Both panels share ONE embedding, so any difference between them
    is purely about which grouping explains the structure -- a fair
    comparison rather than two separate projections.
    """
    if StandardScaler is None or PCA is None:
        raise ImportError("scikit-learn is required (pip install scikit-learn)")

    labeled = attach_method_label(matrix, gsm_method)
    print("  samples per source:", labeled["source"].value_counts().to_dict())
    print("  samples per method:", labeled["method"].value_counts().to_dict())

    filled = prepare_for_clustering(labeled, gene_cols)
    embedding, axis_label = _reduce_to_2d(filled[gene_cols], use_umap)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2))
    _density_panel(axes[0], embedding, filled["source"].reset_index(drop=True),
                    "(a) Grouped by data source", axis_label, "Source")
    _density_panel(axes[1], embedding, filled["method"].reset_index(drop=True),
                    "(b) Grouped by measurement method", axis_label, "Method")

    fig.suptitle(
        f"Transcriptomic cluster density across {len(filled):,} cell line samples\n"
        f"{len(gene_cols):,} genes common to DepMap, HPA and GEO",
        fontsize=11.5, y=1.005,
    )
    fig.tight_layout()
    _safe_savefig(fig, out_path)
    return out_path


def plot_lineage_density(matrix: pd.DataFrame, gene_cols: list, out_path: Path,
                          use_umap: bool = True, top_lineages: int = 8) -> Optional[Path]:
    """
    Single-panel density coloured by lineage, restricted to the
    top_lineages most frequent tissues (rest pooled as 'other') --
    with 30+ lineages every colour becomes indistinguishable and the
    figure stops communicating.
    """
    if "lineage" not in matrix.columns or matrix["lineage"].isna().all():
        print("    (no lineage data -- skipping lineage figure)")
        return None
    if StandardScaler is None or PCA is None:
        raise ImportError("scikit-learn is required (pip install scikit-learn)")

    filled = prepare_for_clustering(matrix, gene_cols)
    lineage = filled["lineage"].fillna("unknown").reset_index(drop=True)
    keep = lineage.value_counts().head(top_lineages).index
    lineage_binned = lineage.where(lineage.isin(keep), "other")

    embedding, axis_label = _reduce_to_2d(filled[gene_cols], use_umap)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    _density_panel(ax, embedding, lineage_binned,
                    f"Cluster density by tissue lineage (top {top_lineages})",
                    axis_label, "Lineage", max_legend=top_lineages + 1)
    fig.suptitle(f"{len(filled):,} samples, {len(gene_cols):,} genes", fontsize=10, y=0.995)
    fig.tight_layout()
    _safe_savefig(fig, out_path)
    return out_path


def score_clustering_concordance(matrix: pd.DataFrame, gsm_method: pd.Series, gene_cols: list,
                                  out_path: Path, n_components: int = 30) -> pd.DataFrame:
    """
    KMeans scored by Adjusted Rand Index against source, method and
    lineage. High ARI for source/method and low for lineage means
    technical variation still dominates biology, and the
    harmonisation is not sufficient for cross-source comparison.
    """
    if StandardScaler is None or PCA is None or KMeans is None:
        raise ImportError("scikit-learn is required (pip install scikit-learn)")

    labeled = attach_method_label(matrix, gsm_method)
    filled = prepare_for_clustering(labeled, gene_cols)

    scaled = StandardScaler().fit_transform(filled[gene_cols])
    n_comp = max(2, min(n_components, scaled.shape[0] - 1, scaled.shape[1]))
    pcs = PCA(n_components=n_comp, random_state=0).fit_transform(scaled)

    rows = [{"metric": "n_samples", "value": len(filled)},
            {"metric": "n_pcs_used", "value": n_comp}]

    for col in ["source", "method", "lineage"]:
        if col not in filled.columns:
            continue
        groups = filled[col].fillna("unknown").reset_index(drop=True)
        k = max(2, min(groups.nunique(), len(groups) - 1))
        labels = KMeans(n_clusters=k, random_state=0, n_init=10).fit_predict(pcs)
        rows.append({"metric": f"n_groups_{col}", "value": groups.nunique()})
        rows.append({"metric": f"ari_vs_{col}", "value": adjusted_rand_score(groups, labels)})

    result = pd.DataFrame(rows)
    result.to_csv(out_path, index=False)
    return result


# ============================================================
# 11. Robust z-scores -- derived scale factor
# ============================================================

def mad_consistency_constant() -> float:
    """
    The factor converting MAD to a standard-deviation-equivalent
    scale, DERIVED rather than hardcoded.

    MAD estimates median(|x - median(x)|). For normally distributed
    data that equals 0.6745*sigma, because 0.6745 is the 0.75
    quantile of the standard normal. So sigma ~= MAD / Phi^-1(0.75),
    and the constant is 1/Phi^-1(0.75).

    Computing it from the quantile function makes the normality
    assumption explicit and auditable -- the number is a consequence
    of a stated model, not a magic literal.
    """
    if norm is not None:
        return float(1.0 / norm.ppf(0.75))
    return 1.4826022185056018  # closed form, for scipy-less environments


MAD_SCALE = mad_consistency_constant()


def estimate_empirical_mad_scale(matrix: pd.DataFrame, gene_cols: list,
                                  sample_genes: int = 500, seed: int = 0) -> dict:
    """
    Empirical check on the derived constant: for a random gene
    sample, the observed SD/MAD ratio across cell lines.

    Under exact normality this centres on MAD_SCALE (~1.4826). A
    substantially higher median indicates heavier tails than normal
    -- expected for expression data, and precisely why MAD is used
    instead of SD. Report this to justify the robust estimator rather
    than asserting it.
    """
    rng = np.random.default_rng(seed)
    chosen = rng.choice(gene_cols, size=min(sample_genes, len(gene_cols)), replace=False)

    ratios = []
    for gene in chosen:
        vals = pd.to_numeric(matrix[gene], errors="coerce").dropna()
        if len(vals) < 20:
            continue
        mad = (vals - vals.median()).abs().median()
        sd = vals.std()
        if mad > 0 and np.isfinite(sd):
            ratios.append(sd / mad)

    if not ratios:
        return {"derived_constant": MAD_SCALE, "empirical_median_ratio": np.nan,
                "empirical_q1": np.nan, "empirical_q3": np.nan, "n_genes_used": 0}

    r = pd.Series(ratios)
    return {
        "derived_constant": MAD_SCALE,
        "empirical_median_ratio": float(r.median()),
        "empirical_q1": float(r.quantile(0.25)),
        "empirical_q3": float(r.quantile(0.75)),
        "n_genes_used": len(r),
    }


# TPM < 1 is below reliable detection: at that level a genuinely
# silent gene and a failed measurement are indistinguishable, so the
# value carries no usable information. On log2(TPM+1), TPM=1 -> 1.0.
TPM_DETECTION_FLOOR = 1.0
LOG2_DETECTION_FLOOR = float(np.log2(TPM_DETECTION_FLOOR + 1))

# Only these are genuinely TPM-derived. GEO microarray values are
# intensities on an unrelated scale where 1.0 means nothing, so a TPM
# floor there would be a category error.
TPM_SCALE_SOURCES = {"depmap", "hpa"}


def build_sample_metadata(matrix: pd.DataFrame, gsm_method: pd.Series) -> pd.DataFrame:
    """One row per sample: identity and stratum labels, no expression values."""
    labeled = attach_stratum_label(matrix, gsm_method)
    keep = ["source", "method", "stratum", "sample_id", "model_id"]
    if "lineage" in labeled.columns:
        keep.append("lineage")
    meta = labeled[keep].copy()
    if "lineage" in meta.columns:
        meta["lineage"] = meta["lineage"].fillna("unknown")
    return meta.reset_index(drop=True)


def load_expression_long(con, matrix: pd.DataFrame, genes: list, sample_meta: pd.DataFrame,
                          table_name: str = "expr_long", chunk: int = 500,
                          max_workers: int = MAX_WORKERS) -> int:
    """
    Melts the wide matrix into long (sample, gene, value) form in
    DuckDB, in column chunks -- an 18k x 5k melt in pandas would
    materialise ~90M rows at once.

    Melting runs CONCURRENTLY across threads (pure pandas work,
    which releases the GIL for the underlying numpy operations),
    while inserts are serialised through a lock since they mutate one
    table.

    NULL log2_value rows are KEPT, so an unmeasured (gene, sample)
    pair is an explicit row rather than an absence to be inferred.
    """
    import threading

    con.execute(f'DROP TABLE IF EXISTS "{table_name}"')
    con.register("sample_meta_src", sample_meta)
    con.execute("CREATE OR REPLACE TABLE sample_meta AS SELECT * FROM sample_meta_src")

    has_lineage = "lineage" in sample_meta.columns
    tpm_sources = ", ".join(f"'{s}'" for s in TPM_SCALE_SOURCES)
    batches = [genes[i:i + chunk] for i in range(0, len(genes), chunk)]

    con.execute(f'''
        CREATE TABLE "{table_name}" (
            sample_id VARCHAR, gene_id VARCHAR, log2_value DOUBLE,
            source VARCHAR, method VARCHAR, stratum VARCHAR, model_id VARCHAR,
            lineage VARCHAR, below_detection BOOLEAN, not_measured BOOLEAN
        )
    ''')

    insert_lock = threading.Lock()
    counter = {"rows": 0, "done": 0}

    def _melt_and_insert(batch):
        long_chunk = matrix[["source", "sample_id"] + batch].melt(
            id_vars=["source", "sample_id"], var_name="gene_id", value_name="log2_value")
        if long_chunk.empty:
            return 0

        with insert_lock:
            cur = con.cursor()
            cur.register("long_chunk_src", long_chunk)
            cur.execute(f'''
                INSERT INTO "{table_name}"
                SELECT l.sample_id, l.gene_id, l.log2_value,
                       m.source, m.method, m.stratum, m.model_id,
                       {"m.lineage" if has_lineage else "'unknown'"},
                       (m.source IN ({tpm_sources}) AND l.log2_value IS NOT NULL
                        AND l.log2_value < {LOG2_DETECTION_FLOOR}),
                       (l.log2_value IS NULL)
                FROM long_chunk_src l
                JOIN sample_meta m ON l.sample_id = m.sample_id AND l.source = m.source
            ''')
            cur.unregister("long_chunk_src")
            cur.close()
            counter["rows"] += len(long_chunk)
            counter["done"] += 1
            print(f"    chunk {counter['done']}/{len(batches)}  ({counter['rows']:,} rows)")
        return len(long_chunk)

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        list(as_completed([ex.submit(_melt_and_insert, b) for b in batches]))

    con.execute(f'CREATE INDEX IF NOT EXISTS idx_{table_name}_gene ON "{table_name}"(gene_id)')
    return counter["rows"]


def compute_all_gene_zscores(con, expr_table: str = "expr_long",
                              out_table: str = "gene_sample_zscores",
                              min_group_size: int = 5,
                              mad_scale: float = MAD_SCALE) -> int:
    """
    Robust z-scores for every gene x sample, in DuckDB via window
    functions -- one pass for medians, one for MAD (which needs the
    median first), independent of gene count. Already parallel across
    cores internally.

      z_stratum -- vs samples in the same source|method stratum.
                   PRIMARY: strips technical artefact, preserves the
                   biological differences between lineages.
      z_lineage -- vs samples of the same lineage. SECONDARY
                   diagnostic, NOT a ranking signal for
                   lineage-specific expression, since grouping by
                   lineage defines that signal away.

    Below-detection and unmeasured rows are excluded from both
    baselines and score NULL -- "could not measure" is a different
    claim from "measured as low". Groups under min_group_size score
    NULL rather than estimating dispersion from too few points.
    """
    con.execute(f'DROP TABLE IF EXISTS "{out_table}"')
    con.execute(f'''
        CREATE TABLE "{out_table}" AS
        WITH usable AS (
            SELECT *, CASE WHEN below_detection OR not_measured THEN NULL
                           ELSE log2_value END AS usable_value
            FROM "{expr_table}"
        ),
        medians AS (
            SELECT *,
                median(usable_value) OVER w_strat AS med_stratum,
                count(usable_value)  OVER w_strat AS n_stratum,
                median(usable_value) OVER w_lin   AS med_lineage,
                count(usable_value)  OVER w_lin   AS n_lineage
            FROM usable
            WINDOW w_strat AS (PARTITION BY gene_id, stratum),
                   w_lin   AS (PARTITION BY gene_id, lineage)
        ),
        mads AS (
            SELECT *,
                median(CASE WHEN usable_value IS NULL THEN NULL
                            ELSE abs(usable_value - med_stratum) END)
                    OVER (PARTITION BY gene_id, stratum) AS mad_stratum,
                median(CASE WHEN usable_value IS NULL THEN NULL
                            ELSE abs(usable_value - med_lineage) END)
                    OVER (PARTITION BY gene_id, lineage) AS mad_lineage
            FROM medians
        )
        SELECT gene_id, sample_id, model_id, source, method, stratum, lineage,
               log2_value, below_detection, not_measured,
               CASE WHEN usable_value IS NULL OR n_stratum < {min_group_size}
                         OR mad_stratum IS NULL OR mad_stratum = 0
                    THEN NULL
                    ELSE (log2_value - med_stratum) / ({mad_scale} * mad_stratum)
               END AS z_stratum,
               CASE WHEN usable_value IS NULL OR n_lineage < {min_group_size}
                         OR mad_lineage IS NULL OR mad_lineage = 0
                    THEN NULL
                    ELSE (log2_value - med_lineage) / ({mad_scale} * mad_lineage)
               END AS z_lineage
        FROM mads
    ''')
    con.execute(f'CREATE INDEX IF NOT EXISTS idx_{out_table}_gene ON "{out_table}"(gene_id)')
    con.execute(f'CREATE INDEX IF NOT EXISTS idx_{out_table}_model ON "{out_table}"(model_id)')
    return con.execute(f'SELECT count(*) FROM "{out_table}"').fetchone()[0]


def collapse_all_to_cell_line(con, zscore_table: str = "gene_sample_zscores",
                               out_table: str = "gene_cell_line_scores") -> int:
    """
    Collapses per-sample scores to one row per (gene, cell line):

      1. median z_stratum within (gene, model_id, source) -- so a cell
         line with 40 GEO samples doesn't outvote one with a single
         DepMap sample.
      2. combine across sources by STOUFFER'S METHOD,
         sum(z)/sqrt(n_sources). DepMap, HPA and GEO are independent
         experiments from different institutions, so agreement across
         them is corroborating evidence; plain averaging (kept as
         z_mean) shrinks toward zero as sources are added.

    EVERY (gene, cell line) pair with any sample is kept. A NULL
    z_stouffer is diagnosable via n_not_measured /
    n_below_detection / n_scoreable_total rather than silently absent.
    """
    con.execute(f'DROP TABLE IF EXISTS "{out_table}"')
    con.execute(f'''
        CREATE TABLE "{out_table}" AS
        WITH per_source AS (
            SELECT gene_id, model_id, lineage, source,
                   median(z_stratum) AS z_source,
                   count(*) AS n_samples,
                   count(z_stratum) AS n_scoreable,
                   count(*) FILTER (WHERE below_detection) AS n_below_detection,
                   count(*) FILTER (WHERE not_measured) AS n_not_measured
            FROM "{zscore_table}"
            WHERE model_id IS NOT NULL
            GROUP BY gene_id, model_id, lineage, source
        )
        SELECT gene_id, model_id, any_value(lineage) AS lineage,
               sum(z_source) / sqrt(count(z_source)) AS z_stouffer,
               avg(z_source) AS z_mean,
               count(z_source) AS n_sources_scored,
               count(DISTINCT source) AS n_sources_present,
               max(z_source) - min(z_source) AS source_spread,
               sum(n_samples) AS n_samples_total,
               sum(n_scoreable) AS n_scoreable_total,
               sum(n_below_detection) AS n_below_detection,
               sum(n_not_measured) AS n_not_measured,
               max(CASE WHEN source = 'depmap' THEN z_source END) AS z_depmap,
               max(CASE WHEN source = 'geo'    THEN z_source END) AS z_geo,
               max(CASE WHEN source = 'hpa'    THEN z_source END) AS z_hpa
        FROM per_source
        GROUP BY gene_id, model_id
    ''')
    con.execute(f'CREATE INDEX IF NOT EXISTS idx_{out_table}_gene ON "{out_table}"(gene_id)')
    con.execute(f'CREATE INDEX IF NOT EXISTS idx_{out_table}_model ON "{out_table}"(model_id)')
    return con.execute(f'SELECT count(*) FROM "{out_table}"').fetchone()[0]


def source_concordance(con, table: str = "gene_cell_line_scores", min_pairs: int = 50) -> pd.DataFrame:
    """
    Pearson correlation between each pair of sources' per-cell-line
    scores, across all genes.

    The empirical test of whether pooling is valid at all. Strong
    positive correlation means the institutions agree and z_stouffer
    is meaningful. Weak correlation means they measure materially
    different things, and per-source scores should be reported
    separately rather than combined.
    """
    pairs = [("z_depmap", "z_geo"), ("z_depmap", "z_hpa"), ("z_geo", "z_hpa")]
    rows = []
    for a, b in pairs:
        res = con.execute(f'''
            SELECT corr({a}, {b}) AS r,
                   count(*) FILTER (WHERE {a} IS NOT NULL AND {b} IS NOT NULL) AS n
            FROM "{table}"
        ''').fetchone()
        rows.append({
            "source_a": a.replace("z_", ""),
            "source_b": b.replace("z_", ""),
            "pearson_r": res[0] if res[1] and res[1] >= min_pairs else np.nan,
            "n_gene_cell_line_pairs": res[1],
        })
    return pd.DataFrame(rows)


def build_all_gene_scores(matrix: pd.DataFrame, genes: list, gsm_method: pd.Series,
                           db_path: Path, force: bool = False) -> dict:
    """
    All-gene scoring into a persistent DuckDB store: expr_long,
    gene_sample_zscores, gene_cell_line_scores.

    Uses the RAW matrix, not the per-source-standardised one --
    standardize_per_source has already z-scored that, and robust
    scoring on top would be double normalisation.
    """
    if duckdb is None:
        raise ImportError("duckdb is required (pip install duckdb)")

    db_path = Path(db_path)
    if db_path.exists() and not force:
        con = duckdb.connect(str(db_path))
        existing = {r[0] for r in con.execute("SHOW TABLES").fetchall()}
        if {"expr_long", "gene_sample_zscores", "gene_cell_line_scores"} <= existing:
            print(f"  using existing score database at {db_path}")
            return {"con": con, "rebuilt": False}
        con.close()

    db_path.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(db_path))
    con.execute(f"PRAGMA threads={MAX_WORKERS}")

    sample_meta = build_sample_metadata(matrix, gsm_method)
    print(f"  melting {len(genes):,} genes into long format ({MAX_WORKERS} workers)...")
    n_long = load_expression_long(con, matrix, genes, sample_meta)
    print(f"    {n_long:,} (gene, sample) rows")

    print(f"  computing robust z-scores (MAD scale factor = {MAD_SCALE:.6f})...")
    n_z = compute_all_gene_zscores(con)
    print(f"    {n_z:,} scored rows")

    print("  collapsing to per-(gene, cell line) scores via Stouffer's method...")
    n_cl = collapse_all_to_cell_line(con)
    print(f"    {n_cl:,} (gene, cell line) rows")

    return {"con": con, "rebuilt": True}


# ============================================================
# 12. Runner
# ============================================================

def run_transcriptomics_cluster_pipeline(con, tables: dict, out_dir: Path,
                                          top_n_genes: Optional[int] = None,
                                          force: bool = False) -> dict:
    """
    GEO metadata -> scale verification -> lineage -> gene selection
    (all protein-coding for scoring, 3-way intersection for
    clustering) -> integrated matrix -> density figures -> clustering
    concordance -> all-gene robust z-score database.

    No cell line is ever dropped for missingness.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading GEO per-GSM methods...")
    gsm_method = load_or_harvest_geo_methods(con, out_dir / "geo_gsm_methods.csv", force=force)

    print("\nVerifying scale per source...")
    scale_summary = cached_csv(
        out_dir / "scale_summary_by_source.csv",
        lambda: pd.DataFrame([
            summarize_source_scale(con, "depmap_expr", "profile_id", "depmap_expr"),
            summarize_source_scale(con, "geo_expr", "gsm_id", "geo_expr (all methods)"),
        ]),
        force=force,
    )
    print(scale_summary.to_string(index=False))

    print("\nVerifying scale per GEO method...")
    geo_scale = cached_csv(
        out_dir / "scale_summary_by_geo_method.csv",
        lambda: summarize_geo_scale_by_method(con, gsm_method),
        force=force,
    )
    print(geo_scale.to_string(index=False))

    print("\nBuilding cell line lineage table...")
    cell_line_roster = tables.get("cell_line_roster")
    if cell_line_roster is None:
        raise ValueError("cell_line_roster not found -- run 02_build_rosters.py first")
    all_cell_lines = cached_csv(
        out_dir / "all_cell_lines_with_lineage.csv",
        lambda: build_all_cell_lines_with_lineage(tables, cell_line_roster),
        force=force,
    )
    print(f"  {len(all_cell_lines):,} cell lines, "
          f"{all_cell_lines['model_id'].notna().sum():,} with a model_id")

    print("\nSelecting genes...")
    gene_roster = tables.get("gene_roster")
    if gene_roster is None:
        raise ValueError("gene_roster not found -- run 02_build_rosters.py first")

    all_genes = cached_gene_list(
        out_dir / "all_protein_coding_genes.csv",
        lambda: select_all_protein_coding_genes(gene_roster),
        force=force,
    )
    cluster_genes = cached_gene_list(
        out_dir / "cluster_genes_three_way.csv",
        lambda: select_common_genes(con, gene_roster, top_n=top_n_genes),
        force=force,
    )
    print(f"  scoring set:    {len(all_genes):,} protein-coding genes")
    print(f"  clustering set: {len(cluster_genes):,} genes (3-way intersection)")

    print(f"\nBuilding integrated matrix over all {len(all_genes):,} genes...")
    matrix = cached_parquet(
        out_dir / "integrated_matrix.parquet",
        lambda: build_integrated_matrix(con, all_genes, gsm_method, all_cell_lines),
        force=force,
    )
    print(f"  {len(matrix):,} samples x {len(all_genes):,} genes")
    print("  rows per source:", matrix["source"].value_counts().to_dict())

    # Genes actually present as columns (a roster gene absent from all
    # three sources never becomes a column).
    matrix_genes = [g for g in all_genes if g in matrix.columns]
    cluster_genes = [g for g in cluster_genes if g in matrix.columns]

    print("\nCalibrating MAD scale factor against the data...")
    mad_check = cached_csv(
        out_dir / "mad_scale_calibration.csv",
        lambda: pd.DataFrame([estimate_empirical_mad_scale(matrix, matrix_genes)]),
        force=force,
    )
    print(mad_check.to_string(index=False))

    print("\nApplying per-source standardisation (visualisation only)...")
    matrix_std = cached_parquet(
        out_dir / "integrated_matrix_standardized.parquet",
        lambda: standardize_per_source(matrix, cluster_genes),
        force=force,
    )

    print("\nPlotting source/method density figure...")
    try:
        p = plot_source_and_method_density(
            matrix_std, gsm_method, cluster_genes, out_dir / "fig_density_source_method.png")
        print(f"  saved -> {p}")
    except Exception as e:
        log_error(logger, step="plot_source_and_method_density", error=e)
        print("  [FAILED] density figure errored — see log")

    print("\nPlotting lineage density figure...")
    try:
        p = plot_lineage_density(matrix_std, cluster_genes, out_dir / "fig_density_lineage.png")
        if p:
            print(f"  saved -> {p}")
    except Exception as e:
        log_error(logger, step="plot_lineage_density", error=e)
        print("  [FAILED] lineage figure errored — see log")

    print("\nScoring clustering concordance (ARI vs source/method/lineage)...")
    try:
        ari = score_clustering_concordance(
            matrix_std, gsm_method, cluster_genes, out_dir / "clustering_concordance.csv")
        print(ari.to_string(index=False))
    except Exception as e:
        log_error(logger, step="score_clustering_concordance", error=e)
        ari = None

    print(f"\nBuilding robust z-score database over all {len(matrix_genes):,} genes...")
    score_db = None
    try:
        score_db = build_all_gene_scores(
            matrix, matrix_genes, gsm_method, out_dir / "gene_scores.duckdb", force=force)
        conc = source_concordance(score_db["con"])
        conc.to_csv(out_dir / "source_concordance.csv", index=False)
        print("\n  cross-institution concordance:")
        print(conc.to_string(index=False))
    except Exception as e:
        log_error(logger, step="build_all_gene_scores", error=e)
        print("  [FAILED] score database build errored — see log")

    print(f"\nDone. Output saved to {out_dir}")
    return {"gsm_method": gsm_method, "all_cell_lines": all_cell_lines,
            "matrix": matrix, "matrix_standardized": matrix_std,
            "all_genes": matrix_genes, "cluster_genes": cluster_genes,
            "clustering_concordance": ari, "score_db": score_db}