"""
transcriptomics.py

Integrates depmap_expr, hpa_rna, and geo_expr onto a common log2
scale and produces 2D cluster density maps spanning all three
transcriptomic sources.

Why depmap/hpa looked "invisible" in earlier plots: they WERE in the
data -- but GEO pools many independent studies on different
platforms/years, so its cross-study variance dwarfs depmap's/hpa's
internal variance, collapsing them to a single point on the same
axes. standardize_per_source() is the fix: z-score each gene WITHIN
each source before combining.

Gene selection: select_common_genes requires a gene to be covered in
ALL THREE sources (depmap_expr, geo_expr, AND hpa_rna).

Sample inclusion: EVERY sample from EVERY source is kept regardless
of missingness -- prepare_for_clustering fills missing gene values
with 0 (that source's own per-gene mean, post-standardization)
rather than dropping samples.

ASSUMPTIONS NEEDING VERIFICATION AGAINST YOUR SCHEMA:
  - hpa_rna's expression VALUE column name (auto-detected)
  - GEO methods other than RMA/gcRMA treated as LINEAR and
    log2-transformed -- verify with summarize_geo_scale_by_method
  - Technology mapping (see _SOURCE_TECHNOLOGY / _MICROARRAY_METHODS):
    depmap and hpa are RNA-seq; GEO's RMA/gcRMA/MAS5 are Affymetrix
    MICROARRAY pipelines. GEO samples that don't match land in
    'geo_unclassified' rather than being guessed at.
  - standardize_per_source is a LIGHTWEIGHT partial batch correction,
    not a substitute for ComBat.

Caching: every expensive step is written to out_dir and re-loaded on
a rerun -- see cached_csv/cached_parquet/cached_gene_list and force=.

Dependencies: GEOparse, scikit-learn, optionally umap-learn.
    pip install GEOparse scikit-learn umap-learn
"""

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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
    from scipy.stats import gaussian_kde  # type: ignore[import-not-found]
except ImportError:
    gaussian_kde = None

from src.scripts.eda_functions import style_axis, _safe_savefig, cols_of, numeric_summary, cached_csv
from src.scripts.logging_utils import get_logger, log_error

try:
    import GEOparse  # type: ignore[import-not-found]
except ImportError:
    GEOparse = None

logger = get_logger("transcriptomics")


# ============================================================
# 0. Caching helpers
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
# 1. GEO per-GSM method metadata
# ============================================================

def harvest_geo_methods(gsms: list, destdir: str = "./geo_meta", verbose: bool = True) -> pd.DataFrame:
    if GEOparse is None:
        raise ImportError("GEOparse is required for harvest_geo_methods(). Install it with: pip install GEOparse")

    rows = []
    for i, g in enumerate(gsms, 1):
        try:
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
            print(f"  {i:,} / {len(gsms):,}")

    meta = pd.DataFrame(rows)

    # gcRMA MUST be tested before RMA -- the string "gcRMA" contains "RMA"
    meta["method"] = np.select(
        [meta.processing.str.contains("gcrma", case=False, na=False),
         meta.processing.str.contains("rma|robust multi", case=False, na=False),
         meta.processing.str.contains("mas5|mas 5|gcos", case=False, na=False),
         meta.processing.str.contains("quantile", case=False, na=False)],
        ["gcRMA", "RMA", "MAS5", "quantile"], default="other")

    if verbose:
        print(f"\n{len(meta):,} GSMs\n")
        print(meta.platform.value_counts().to_string())
        print()
        print(meta.method.value_counts().to_string())
        print(f"\nwith ABS_CALL (Affymetrix present/absent): {meta.has_abs_call.sum():,}")
        failed = meta.processing.str.startswith("FETCH FAILED").sum()
        if failed:
            print(f"FAILED to fetch: {failed:,}")
    return meta


def load_or_harvest_geo_methods(con, cache_path, force: bool = False) -> pd.Series:
    cache_path = Path(cache_path)
    if not force and cache_path.exists():
        meta = pd.read_csv(cache_path)
        print(f"loaded cached GEO metadata: {len(meta):,} GSMs")
    else:
        gsms = con.execute(
            "SELECT DISTINCT geo_accession FROM geo_info WHERE geo_accession IS NOT NULL"
        ).df()["geo_accession"].tolist()
        meta = harvest_geo_methods(gsms)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        meta.to_csv(cache_path, index=False)
        print(f"harvested and cached {len(meta):,} GSMs to {cache_path}")

    gsm_method = meta.set_index(meta["gsm"].str.lower())["method"]
    print(gsm_method.value_counts().to_string())
    return gsm_method


def register_gsm_method(con, gsm_method: pd.Series) -> None:
    df = gsm_method.rename("method").rename_axis("gsm").reset_index()
    con.register("gsm_method_lookup", df)


# ============================================================
# 2. Scale verification
# ============================================================

def row_wise_scale_stats(con, table: str, id_col: str) -> pd.DataFrame:
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
    stats = row_wise_scale_stats(con, table, id_col)
    stats["method"] = stats[id_col].str.lower().map(gsm_method).fillna("unknown")

    summary = stats.groupby("method").agg(
        n_samples=("row_max", "size"),
        median_row_mean=("row_mean", "median"),
        median_row_max=("row_max", "median"),
        min_row_max=("row_max", "min"),
        max_row_max=("row_max", "max"),
    ).reset_index()
    return summary


# ============================================================
# 3. Scale harmonization -- everything to log2(x+1)
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
# 4. Cell-line lineage
# ============================================================

def build_all_cell_lines_with_lineage(tables: dict, cell_line_roster: pd.DataFrame) -> pd.DataFrame:
    sample_info = tables.get("sample_info")
    if sample_info is None:
        raise ValueError("sample_info table not found in `tables`")
    if "sample_collection_site" not in sample_info.columns:
        raise ValueError(
            f"Expected a 'sample_collection_site' column in sample_info. "
            f"Available columns: {sample_info.columns.tolist()}"
        )
    if "model_id" not in sample_info.columns:
        raise ValueError("sample_info has no model_id column -- can't join lineage by model_id")

    lineage_lookup = (
        sample_info[["model_id", "sample_collection_site"]]
        .rename(columns={"sample_collection_site": "lineage"})
        .dropna(subset=["model_id", "lineage"])
        .drop_duplicates()
    )
    lineage_lookup = (
        lineage_lookup.groupby("model_id")["lineage"]
        .apply(lambda s: "; ".join(sorted(set(s))))
        .reset_index()
    )

    result = cell_line_roster.merge(lineage_lookup, on="model_id", how="left")

    if "cell_line_name(s)" in result.columns:
        result["cell_line_name"] = (
            result["cell_line_name(s)"]
            .fillna("")
            .str.split("; ")
            .str[0]
            .replace("", np.nan)
        )

    return result


# ============================================================
# 5. Gene selection -- ALL THREE sources required
# ============================================================

def select_common_genes(con, gene_roster: pd.DataFrame, top_n: Optional[int] = None) -> list:
    """
    Genes covered in ALL THREE sources: depmap_expr (via
    gene_roster's profile_id(s)), geo_expr (via geo_gsm_id(s)), AND
    hpa_rna (queried directly, since gene_roster doesn't track hpa
    coverage as a column). Restricted to protein-coding gene_ids.

    top_n=None (default) returns every 3-way-covered gene.
    """
    depmap_geo_covered = set(gene_roster[
        (gene_roster["profile_id(s)"].fillna("") != "") &
        (gene_roster["geo_gsm_id(s)"].fillna("") != "")
    ]["gene_id"])
    print(f"  genes covered in depmap_expr AND geo_expr: {len(depmap_geo_covered):,}")

    hpa_covered = set(con.execute("SELECT DISTINCT gene_id FROM hpa_rna WHERE gene_id IS NOT NULL")
                       .fetchdf()["gene_id"])
    print(f"  genes covered in hpa_rna: {len(hpa_covered):,}")

    covered = sorted(depmap_geo_covered & hpa_covered)
    print(f"  genes covered in ALL THREE sources: {len(covered):,}")

    if not covered:
        raise ValueError("No protein-coding genes found with coverage in depmap_expr, geo_expr, AND hpa_rna.")

    if top_n is None:
        return covered

    summary = numeric_summary(con, "depmap_expr", cols=covered)
    top = summary.assign(var=lambda d: d["std"] ** 2).sort_values("var", ascending=False).head(top_n)
    return top["column"].tolist()


# ============================================================
# 6. Integrated matrix -- ALL samples from ALL sources
# ============================================================

def _fetch_depmap_wide(con, genes: list, chunk: int = 300) -> pd.DataFrame:
    """Fetch depmap_expr's profile_id + selected gene columns in batches, joined on profile_id."""
    frames = []
    for i in range(0, len(genes), chunk):
        batch = genes[i:i + chunk]
        cols_sql = ", ".join(f'"{g}"' for g in batch)
        df = con.execute(f'SELECT profile_id, {cols_sql} FROM depmap_expr').fetchdf()
        frames.append(df.set_index("profile_id"))
    combined = pd.concat(frames, axis=1)
    return combined.reset_index()


def build_integrated_matrix(con, genes: list, gsm_method: pd.Series,
                             all_cell_lines: pd.DataFrame,
                             hpa_value_col: Optional[str] = None) -> pd.DataFrame:
    """
    Build one integrated expression matrix. EVERY sample from EVERY
    source is included (no cell-line-overlap filtering).
    """
    gene_list_sql = ", ".join(f"'{g}'" for g in genes)

    depmap = _fetch_depmap_wide(con, genes)
    depmap_profiles = con.execute('SELECT profile_id, model_id FROM depmap_profiles').fetchdf()
    depmap = depmap.merge(depmap_profiles, on="profile_id", how="left")
    depmap_long = depmap.melt(id_vars=["profile_id", "model_id"], var_name="gene", value_name="log2_value")
    depmap_long["sample_id"] = depmap_long["profile_id"]
    depmap_long["source"] = "depmap"

    register_gsm_method(con, gsm_method)
    geo_query = harmonized_geo_long_query()
    geo_long = con.execute(f'''
        SELECT gsm_id, gene, log2_value FROM ({geo_query}) WHERE gene IN ({gene_list_sql})
    ''').fetchdf()
    geo_info = con.execute('SELECT geo_accession, model_id FROM geo_info WHERE geo_accession IS NOT NULL').fetchdf()
    geo_long = geo_long.merge(geo_info, left_on="gsm_id", right_on="geo_accession", how="left")
    geo_long["sample_id"] = geo_long["gsm_id"]
    geo_long["source"] = "geo"

    hpa_query = harmonized_hpa_query(con, value_col=hpa_value_col)
    hpa_long = con.execute(f'''
        SELECT gene_id AS gene, model_id, log2_value FROM ({hpa_query}) WHERE gene_id IN ({gene_list_sql})
    ''').fetchdf()
    hpa_long["sample_id"] = hpa_long["model_id"]
    hpa_long["source"] = "hpa"

    combined_long = pd.concat([
        depmap_long[["source", "sample_id", "model_id", "gene", "log2_value"]],
        geo_long[["source", "sample_id", "model_id", "gene", "log2_value"]],
        hpa_long[["source", "sample_id", "model_id", "gene", "log2_value"]],
    ], ignore_index=True)

    wide = combined_long.pivot_table(index=["source", "sample_id", "model_id"], columns="gene", values="log2_value")
    wide = wide.reset_index()
    wide.columns.name = None

    lineage_lookup = (
        all_cell_lines[["model_id", "lineage"]]
        .dropna(subset=["model_id"])
        .drop_duplicates(subset=["model_id"])
    )
    wide = wide.merge(lineage_lookup, on="model_id", how="left")
    return wide


# ============================================================
# 6b. Per-source batch correction + full sample retention
# ============================================================

def standardize_per_source(matrix: pd.DataFrame, gene_cols: list) -> pd.DataFrame:
    """Z-scores each gene column WITHIN each source group. See module docstring."""
    out = matrix.copy()
    for _, idx in matrix.groupby("source").groups.items():
        sub = out.loc[idx, gene_cols]
        mean = sub.mean()
        std = sub.std().replace(0, np.nan)
        out.loc[idx, gene_cols] = (sub - mean) / std
    return out


def prepare_for_clustering(matrix: pd.DataFrame, gene_cols: list) -> pd.DataFrame:
    """
    Keeps EVERY sample regardless of missingness -- no threshold,
    nothing dropped. Missing values are filled with 0, which (since
    this runs AFTER standardize_per_source) is that source's own
    per-gene mean rather than an arbitrary guess.

    Trade-off: a very sparse sample gets pulled toward its source's
    centroid for every gene it's missing, understating how different
    it really is. That's the cost of keeping every sample.
    """
    out = matrix.copy()
    numeric = out[gene_cols].apply(pd.to_numeric, errors="coerce")
    n_missing = numeric.isna().sum(axis=1)
    if (n_missing > 0).any():
        sparse = n_missing[n_missing > 0]
        print(f"    prepare_for_clustering: filling missing values for "
              f"{len(sparse):,} / {len(out):,} samples "
              f"(median {sparse.median():.0f} of {len(gene_cols)} genes missing per sparse sample)")
    out[gene_cols] = numeric.fillna(0)
    return out


# ============================================================
# 7. Labeling: source and assay technology
# ============================================================

# depmap (OmicsExpressionTPMLogp1) and hpa (nTPM) are both RNA-seq.
# GEO is mixed -- technology inferred from per-GSM normalisation
# method, since RMA/gcRMA/MAS5 are Affymetrix MICROARRAY pipelines.
_SOURCE_TECHNOLOGY = {"depmap": "rna_seq", "hpa": "rna_seq"}
_MICROARRAY_METHODS = {"RMA", "gcRMA", "MAS5"}


def attach_method_label(matrix: pd.DataFrame, gsm_method: pd.Series) -> pd.DataFrame:
    """
    Adds a 'method' column: GEO samples get their actual normalisation
    method from gsm_method; depmap/hpa use their source name.
    """
    out = matrix.copy()
    out["method"] = out["source"]
    is_geo = out["source"] == "geo"
    out.loc[is_geo, "method"] = (
        out.loc[is_geo, "sample_id"].astype(str).str.lower().map(gsm_method).fillna("unknown")
    )
    return out


def attach_technology_label(matrix: pd.DataFrame, gsm_method: pd.Series) -> pd.DataFrame:
    """
    Adds a 'technology' column describing HOW each sample was
    measured, as distinct from WHICH source it came from:

      depmap -> rna_seq
      hpa    -> rna_seq
      geo    -> 'microarray' if its method is RMA/gcRMA/MAS5
                (Affymetrix pipelines), else 'geo_unclassified' --
                'quantile'/'other' alone don't identify the platform,
                and a wrong guess would silently mislabel the panel.

    VERIFY against scale_summary_by_geo_method.csv. If your GEO set
    includes genuine RNA-seq submissions they'll currently land in
    'geo_unclassified' rather than 'rna_seq'.
    """
    out = matrix.copy()
    out["technology"] = out["source"].map(_SOURCE_TECHNOLOGY)

    is_geo = out["source"] == "geo"
    geo_methods = out.loc[is_geo, "sample_id"].astype(str).str.lower().map(gsm_method)
    out.loc[is_geo, "technology"] = np.where(
        geo_methods.isin(_MICROARRAY_METHODS), "microarray", "geo_unclassified"
    )

    out["technology"] = out["technology"].fillna("unknown")
    return out


def common_cell_line_model_ids(matrix: pd.DataFrame) -> set:
    """model_ids present in ALL THREE sources (depmap AND geo AND hpa)."""
    sources = matrix["source"].unique().tolist()
    ids_by_source = {
        s: set(matrix.loc[matrix["source"] == s, "model_id"].dropna())
        for s in sources
    }
    if len(ids_by_source) < 3:
        return set()
    return set.intersection(*ids_by_source.values())


# ============================================================
# 8. Plots
# ============================================================

def _reduce_to_2d(numeric: pd.DataFrame, use_umap: bool):
    scaled = StandardScaler().fit_transform(numeric)
    if use_umap and umap is not None:
        return umap.UMAP(n_components=2, random_state=0).fit_transform(scaled), "UMAP"
    return PCA(n_components=2, random_state=0).fit_transform(scaled), "PCA"


def _plot_one_facet(ax, df: pd.DataFrame, gene_cols: list, color_col: str, title: str, use_umap: bool) -> None:
    """Draws one facet onto an existing ax. Never drops rows -- prepare_for_clustering fills instead."""
    if df.empty:
        ax.text(0.5, 0.5, "no samples", ha="center", va="center", color="#999", fontsize=10, transform=ax.transAxes)
        ax.set_title(title, fontsize=9.5)
        ax.set_xticks([]); ax.set_yticks([])
        return

    filled = prepare_for_clustering(df, gene_cols)
    numeric = filled[gene_cols]
    meta = filled[[color_col]].reset_index(drop=True)

    if len(numeric) < 3:
        ax.text(0.5, 0.5, f"only {len(numeric)} sample(s)", ha="center", va="center",
                color="#999", fontsize=10, transform=ax.transAxes)
        ax.set_title(title, fontsize=9.5)
        ax.set_xticks([]); ax.set_yticks([])
        return

    embedding, method_label = _reduce_to_2d(numeric, use_umap)

    if gaussian_kde is not None and len(embedding) > 10:
        try:
            kde = gaussian_kde(embedding.T)
            xmin, xmax = embedding[:, 0].min(), embedding[:, 0].max()
            ymin, ymax = embedding[:, 1].min(), embedding[:, 1].max()
            xx, yy = np.mgrid[xmin:xmax:150j, ymin:ymax:150j]
            zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
            ax.contourf(xx, yy, zz, levels=15, cmap="Blues", alpha=0.5, zorder=1)
        except Exception:
            pass  # near-duplicate points can make KDE singular -- keep scatter, skip contour

    categories = sorted(meta[color_col].fillna("unknown").unique())
    cmap = plt.get_cmap("tab10" if len(categories) <= 10 else "tab20")
    for i, cat in enumerate(categories):
        mask = (meta[color_col].fillna("unknown") == cat).values
        ax.scatter(embedding[mask, 0], embedding[mask, 1], s=10, alpha=0.7,
                   color=cmap(i % cmap.N), label=str(cat), zorder=3)

    ax.legend(fontsize=6.5, markerscale=1.5, loc="best", frameon=False)
    ax.set_xlabel(f"{method_label} 1", fontsize=8)
    ax.set_ylabel(f"{method_label} 2", fontsize=8)
    ax.set_title(f"{title}\n({len(numeric):,} samples)", fontsize=9.5)
    style_axis(ax, grid_axis="both")


def plot_source_and_technology_facet(matrix: pd.DataFrame, gsm_method: pd.Series, gene_cols: list,
                                      out_path: Path, use_umap: bool = True) -> Path:
    """
    Two-panel cluster density figure over ALL cell lines (nothing
    dropped):

      Left:  colored by SOURCE      (depmap / geo / hpa)
      Right: colored by TECHNOLOGY  (rna_seq / microarray / ...)

    Both panels use the same samples -- only the coloring differs, so
    differences between them are purely about which grouping explains
    the structure better.
    """
    if StandardScaler is None or PCA is None:
        raise ImportError("scikit-learn is required for clustering (pip install scikit-learn)")

    labeled = attach_technology_label(matrix, gsm_method)

    print("  samples per source:")
    print(labeled["source"].value_counts().to_string())
    print("  samples per technology:")
    print(labeled["technology"].value_counts().to_string())

    fig, axes = plt.subplots(1, 2, figsize=(17, 7))

    _plot_one_facet(axes[0], labeled, gene_cols, "source",
                     "Cell line expression — colored by SOURCE", use_umap)
    _plot_one_facet(axes[1], labeled, gene_cols, "technology",
                     "Cell line expression — colored by TECHNOLOGY", use_umap)

    fig.suptitle(f"Transcriptomics cluster density — {len(gene_cols):,} genes, "
                 f"{len(labeled):,} samples (no cell lines dropped)",
                 fontsize=13, y=1.02)
    fig.tight_layout()

    _safe_savefig(fig, out_path)
    return out_path


# ============================================================
# 9. Clustering + concordance scoring
# ============================================================

def score_source_and_technology(matrix: pd.DataFrame, gsm_method: pd.Series, gene_cols: list,
                                 out_path: Path, n_components: int = 30) -> pd.DataFrame:
    """
    KMeans over all samples, scored against BOTH groupings via
    Adjusted Rand Index (ARI, -1 to 1, 0 = random). Higher ARI means
    that label explains the clustering better -- i.e. whether the
    structure tracks the data SOURCE or the assay TECHNOLOGY.
    """
    if StandardScaler is None or PCA is None or KMeans is None:
        raise ImportError("scikit-learn is required (pip install scikit-learn)")

    labeled = attach_technology_label(matrix, gsm_method)
    filled = prepare_for_clustering(labeled, gene_cols)
    numeric = filled[gene_cols]

    scaled = StandardScaler().fit_transform(numeric)
    n_comp = max(2, min(n_components, scaled.shape[0] - 1, scaled.shape[1]))
    pcs = PCA(n_components=n_comp, random_state=0).fit_transform(scaled)

    rows = [{"metric": "n_samples", "value": len(filled)}, {"metric": "n_pcs_used", "value": n_comp}]

    for group_col in ["source", "technology"]:
        groups = filled[group_col].fillna("unknown").reset_index(drop=True)
        n_clusters = max(2, min(groups.nunique(), len(groups) - 1))
        labels = KMeans(n_clusters=n_clusters, random_state=0, n_init=10).fit_predict(pcs)
        rows.append({"metric": f"n_clusters_for_{group_col}", "value": n_clusters})
        rows.append({"metric": f"ari_vs_{group_col}", "value": adjusted_rand_score(groups, labels)})

    result = pd.DataFrame(rows)
    result.to_csv(out_path, index=False)
    return result


# ============================================================
# 10. Runner
# ============================================================

def run_transcriptomics_cluster_pipeline(con, tables: dict, out_dir: Path,
                                          top_n_genes: Optional[int] = None, force: bool = False) -> dict:
    """
    Full pipeline: GEO method harvesting -> scale verification ->
    lineage attachment -> gene selection (3-way common genes) ->
    integrated matrix (all samples, all sources) -> per-source
    standardization -> source/technology cluster density facet ->
    KMeans concordance scoring. No sample is ever dropped for
    missingness (see prepare_for_clustering).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading/harvesting GEO per-GSM methods...")
    gsm_method = load_or_harvest_geo_methods(con, cache_path=out_dir / "geo_gsm_methods.csv", force=force)

    print("\nVerifying scale per source...")
    scale_summary = cached_csv(
        out_dir / "scale_summary_by_source.csv",
        lambda: pd.DataFrame([
            summarize_source_scale(con, "depmap_expr", "profile_id", "depmap_expr"),
            summarize_source_scale(con, "geo_expr", "gsm_id", "geo_expr (all methods combined)"),
        ]),
        force=force,
    )
    print(scale_summary.to_string(index=False))

    print("\nVerifying scale per GEO method...")
    geo_scale_by_method = cached_csv(
        out_dir / "scale_summary_by_geo_method.csv",
        lambda: summarize_geo_scale_by_method(con, gsm_method),
        force=force,
    )
    print(geo_scale_by_method.to_string(index=False))

    print("\nBuilding all-cell-lines-with-lineage table...")
    cell_line_roster = tables.get("cell_line_roster")
    if cell_line_roster is None:
        raise ValueError("cell_line_roster not found in tables -- run 02_build_rosters.py first")
    all_cell_lines = cached_csv(
        out_dir / "all_cell_lines_with_lineage.csv",
        lambda: build_all_cell_lines_with_lineage(tables, cell_line_roster),
        force=force,
    )
    print(f"  {len(all_cell_lines):,} cell lines, {all_cell_lines['model_id'].notna().sum():,} with a matched model_id")

    print("\nSelecting genes common to depmap_expr, geo_expr, AND hpa_rna...")
    gene_roster = tables.get("gene_roster")
    if gene_roster is None:
        raise ValueError("gene_roster not found in tables -- run 02_build_rosters.py first")
    genes = cached_gene_list(
        out_dir / "selected_genes.csv",
        lambda: select_common_genes(con, gene_roster, top_n=top_n_genes),
        force=force,
    )
    print(f"  {len(genes)} genes selected"
          f"{' (ALL 3-way-covered protein-coding genes)' if top_n_genes is None else f' (top {top_n_genes} by variance)'}")

    print("\nBuilding integrated matrix (ALL samples from ALL sources)...")
    matrix = cached_parquet(
        out_dir / "integrated_matrix.parquet",
        lambda: build_integrated_matrix(con, genes, gsm_method, all_cell_lines),
        force=force,
    )
    print(f"  {len(matrix):,} samples x {len(genes)} genes")
    print("  rows per source:")
    print(matrix["source"].value_counts().to_string())

    print("\nApplying per-source standardization (mitigates cross-study batch effects)...")
    matrix_std = cached_parquet(
        out_dir / "integrated_matrix_standardized.parquet",
        lambda: standardize_per_source(matrix, genes),
        force=force,
    )

    print("\nPlotting source + technology cluster density facet...")
    try:
        facet_path = plot_source_and_technology_facet(
            matrix_std, gsm_method, genes, out_dir / "cluster_facet_source_technology.png")
        print(f"  saved -> {facet_path}")
    except Exception as e:
        log_error(logger, step="plot_source_and_technology_facet", error=e)
        print("  [FAILED] plot_source_and_technology_facet errored — see log")

    print("\nScoring clusters against source and technology...")
    try:
        cluster_scores = score_source_and_technology(
            matrix_std, gsm_method, genes, out_dir / "source_technology_cluster_scores.csv")
        print(cluster_scores.to_string(index=False))
    except Exception as e:
        log_error(logger, step="score_source_and_technology", error=e)
        print("  [SKIPPED] score_source_and_technology failed — see log")
        cluster_scores = None

    print(f"\nDone. Output saved to {out_dir}")
    return {"gsm_method": gsm_method, "all_cell_lines": all_cell_lines,
            "matrix": matrix, "matrix_standardized": matrix_std, "cluster_scores": cluster_scores}