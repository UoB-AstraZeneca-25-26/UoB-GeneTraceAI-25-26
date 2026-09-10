"""
eda_functions.py

Exploratory data analysis for the cleaned pipeline output, backed by
DuckDB rather than pandas -- every statistic is computed with a
DuckDB aggregate query, and only small, bounded samples are ever
pulled into a pandas DataFrame (for plots that genuinely need
row-level data, like missingno's matrix/heatmap or histograms). This
is what lets this run over very wide tables (depmap_expr, geo_expr,
proteomics can have thousands of columns) without materializing them
in memory.

Plots are plain matplotlib (no seaborn) via the shared `style_axis`
helper, matching the house style used elsewhere in this project.

Speed note: missing_value_report and numeric_summary are each
computed ONCE per table in run_eda_on_table and passed into every
downstream plot/test via report=/summary= -- recomputing either from
scratch inside every plot function is what made wide tables like
depmap_expr (~54k columns) slow. Both also batch every statistic
into ONE query per chunk of columns via chunked_multi_agg.

Caching: missing_value_report.csv and numeric_summary.csv under each
table's output dir act as an on-disk cache (see cached_csv) -- a
rerun of the pipeline loads these instead of recomputing, unless
force=True is passed (04_eda.py exposes this as --force).

Gene-matrix restriction: depmap_expr and geo_expr (see
GENE_MATRIX_TABLES) are restricted to gene_roster's protein-coding
gene_ids for every stat and plot in this file (see
get_protein_coding_gene_ids), instead of profiling every gene column.

Layout of this file:
  1. Setup (style, constants)
  2. DuckDB helpers (chunked aggregation, sampling, caching)
  3. Missing value analysis (per-table report + bar/matrix/heatmap/dendrogram)
  4. Missingness randomness testing (Little's MCAR + dependency heuristic)
  5. Distribution / scale analysis (summary stats + scale profiling + plots)
  6. Cross-table summary plots (span ALL tables in one figure)
  7. Per-table and full-pipeline runners
"""

from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

try:
    from scipy import stats  # type: ignore[import-not-found]
except ImportError:
    stats = None

try:
    import missingno as msno  # type: ignore
except ImportError:
    msno = None

from src.scripts.logging_utils import get_logger, log_error

logger = get_logger("eda_functions")

if stats is None:
    log_error(logger, step="import", error=ImportError("scipy not available"), note="some statistical tests disabled")
if msno is None:
    log_error(logger, step="import", error=ImportError("missingno not available"), note="nullity plots disabled")

# ============================================================
# 1. Setup
# ============================================================

plt.rcParams["figure.dpi"] = 110
plt.rcParams["savefig.dpi"] = 200
plt.rcParams["axes.titleweight"] = "bold"
plt.rcParams["axes.titlesize"] = 13
plt.rcParams["font.size"] = 9

MAX_BAR_COLS = 40
MAX_MATRIX_COLS = 60
MAX_MATRIX_ROWS = 500
MAX_HEATMAP_COLS = 40
MAX_DENDROGRAM_COLS = 60
MAX_MCAR_COLS = 30
MAX_MCAR_ROWS = 5000
MAX_DEPENDENCY_FOCUS = 10
MAX_DEPENDENCY_COMPARE = 20
MAX_DIST_COLS = 12
MAX_SCALE_COLS = 20
MAX_CORR_COLS = 40
MAX_GAUSSIAN_COLS = 250
MAX_AGG_CHUNK = 250
MAX_MULTI_AGG_CHUNK = 100
MAX_RAW_FACET_VALUES = 20_000
GAUSSIAN_FACET_TABLES = {
    "depmap_expr", "geo_expr", "hpa_rna", "proteomics", "metabolomics", "mirna",
}
MISSINGNESS_FACET_TABLES = (
    "depmap_expr", "geo_expr", "hpa_rna", "mutations", "fusions", "proteomics",
    "metabolomics", "mirna", "signatures", "sample_info", "cellosaurus", "hpa_desc",
    "depmap_profiles", "geo_info",
)

GENE_MATRIX_TABLES = {"depmap_expr", "geo_expr"}

from concurrent.futures import ThreadPoolExecutor, as_completed
import os

MAX_WORKERS = min(8, os.cpu_count() or 4)


def _safe_savefig(fig, out_path: Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def style_axis(ax, grid_axis: str = "y") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis=grid_axis, alpha=0.3, zorder=0)
    ax.set_axisbelow(True)


def cached_csv(path: Path, compute_fn, force: bool = False) -> pd.DataFrame:
    """
    Load a DataFrame from `path` if it already exists (and force is
    False); otherwise compute it via compute_fn(), save it to `path`,
    and return it. Delete the CSV, or pass force=True, to recompute.
    """
    path = Path(path)
    if not force and path.exists():
        return pd.read_csv(path)
    df = compute_fn()
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return df


# ============================================================
# 2. DuckDB helpers
# ============================================================

_NUMERIC_DUCKDB_TYPES = {
    "TINYINT", "SMALLINT", "INTEGER", "BIGINT", "HUGEINT",
    "UTINYINT", "USMALLINT", "UINTEGER", "UBIGINT",
    "FLOAT", "DOUBLE", "DECIMAL", "REAL",
}


def cols_of(con, table: str) -> list:
    return con.execute(f'DESCRIBE "{table}"').fetchdf()["column_name"].tolist()


def numeric_cols(con, table: str) -> list:
    desc = con.execute(f'DESCRIBE "{table}"').fetchdf()
    base_type = desc["column_type"].str.upper().str.extract(r"^(\w+)")[0]
    return desc.loc[base_type.isin(_NUMERIC_DUCKDB_TYPES), "column_name"].tolist()


def row_count(con, table: str) -> int:
    return con.execute(f'SELECT count(*) FROM "{table}"').fetchone()[0]


def get_protein_coding_gene_ids(con) -> Optional[list]:
    try:
        return con.execute('SELECT DISTINCT gene_id FROM gene_roster').fetchdf()["gene_id"].tolist()
    except Exception as e:
        log_error(logger, step="get_protein_coding_gene_ids", error=e)
        print("    [WARN] gene_roster not available -- gene-matrix tables will use ALL columns, not just protein-coding")
        return None


def _run_single_chunk(con, table: str, batch: list, expr_fn) -> pd.Series:
    cur = con.cursor()
    select = ", ".join(f'{expr_fn(c)} AS "{c}"' for c in batch)
    row = cur.execute(f'SELECT {select} FROM "{table}"').fetchdf().iloc[0]
    cur.close()
    return row


def chunked_agg(con, table: str, cols: list, expr_fn, chunk: int = MAX_AGG_CHUNK,
                 max_workers: int = MAX_WORKERS) -> pd.Series:
    if not cols:
        return pd.Series(dtype=float)

    batches = [cols[i:i + chunk] for i in range(0, len(cols), chunk)]
    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(_run_single_chunk, con, table, b, expr_fn) for b in batches]
        for future in as_completed(futures):
            results.update(future.result().to_dict())
    return pd.Series(results)


def _run_multi_chunk(con, table: str, batch: list, expr_map: dict) -> pd.Series:
    cur = con.cursor()
    select_parts = [
        f'{expr_fn(c)} AS "{stat_name}__{c}"'
        for stat_name, expr_fn in expr_map.items()
        for c in batch
    ]
    row = cur.execute(f'SELECT {", ".join(select_parts)} FROM "{table}"').fetchdf().iloc[0]
    cur.close()
    return row


def chunked_multi_agg(con, table: str, cols: list, expr_map: dict, chunk: int = MAX_MULTI_AGG_CHUNK,
                       progress_label: Optional[str] = None, max_workers: int = MAX_WORKERS) -> dict:
    if not cols:
        return {name: pd.Series(dtype=float) for name in expr_map}

    batches = [cols[i:i + chunk] for i in range(0, len(cols), chunk)]
    results = {name: {} for name in expr_map}
    n_chunks = len(batches)
    completed = 0

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_run_multi_chunk, con, table, b, expr_map): b for b in batches}
        for future in as_completed(futures):
            row = future.result()
            for key, val in row.items():
                stat_name, col = key.split("__", 1)
                results[stat_name][col] = val
            completed += 1
            if progress_label:
                print(f"    {progress_label}: chunk {completed}/{n_chunks}")

    return {name: pd.Series(vals) for name, vals in results.items()}


def fetch_sample(con, table: str, cols: Optional[list] = None, max_rows: Optional[int] = None,
                  seed: int = 0) -> pd.DataFrame:
    col_sql = "*" if not cols else ", ".join(f'"{c}"' for c in cols)
    query = f'SELECT {col_sql} FROM "{table}"'
    if max_rows is not None:
        total = row_count(con, table)
        if total > max_rows:
            query += f" USING SAMPLE {max_rows} ROWS (reservoir, {seed})"
    return con.execute(query).fetchdf()


def raw_measurement_cols(con, table: str) -> list:
    """Return numeric columns that are likely measurements, not identifiers."""
    excluded_tokens = ("id", "index", "start", "end", "position", "chromosome", "coord")
    return [
        col for col in numeric_cols(con, table)
        if not any(token in col.lower() for token in excluded_tokens)
    ]


def fetch_raw_value_sample(con, table: str, cols: list,
                           max_values: int = MAX_RAW_FACET_VALUES,
                           seed: int = 7) -> np.ndarray:
    """Sample pooled values from wide raw tables without materializing them."""
    if not cols:
        return np.array([], dtype=float)

    max_batches = 5
    if len(cols) > max_batches * 100:
        positions = np.linspace(0, len(cols) - 1, max_batches * 100, dtype=int)
        cols = [cols[position] for position in positions]
    rows_per_query = max(1, min(500, max_values // len(cols)))
    values = []
    for start in range(0, len(cols), 100):
        batch = cols[start:start + 100]
        query_cols = ", ".join(f'"{col}"' for col in batch)
        query = f'SELECT {query_cols} FROM "{table}" LIMIT {rows_per_query}'
        sampled = con.execute(query).fetchdf()
        values.append(sampled.to_numpy(dtype=float, na_value=np.nan).ravel())

    pooled = np.concatenate(values) if values else np.array([], dtype=float)
    pooled = pooled[np.isfinite(pooled)]
    if len(pooled) > max_values:
        rng = np.random.default_rng(seed)
        pooled = rng.choice(pooled, size=max_values, replace=False)
    return pooled


# ============================================================
# 3. Missing value analysis
# ============================================================

def missing_value_report(con, table: str, chunk: int = MAX_MULTI_AGG_CHUNK,
                          cols: Optional[list] = None) -> pd.DataFrame:
    cs = cols if cols is not None else cols_of(con, table)
    n = row_count(con, table)
    if n == 0 or not cs:
        return pd.DataFrame(columns=["column", "dtype", "n_missing", "pct_missing", "n_unique"])

    dtypes = con.execute(f'DESCRIBE "{table}"').fetchdf().set_index("column_name")["column_type"]
    agg = chunked_multi_agg(con, table, cs, {
        "n_missing": lambda c: f'sum(CASE WHEN "{c}" IS NULL THEN 1 ELSE 0 END)',
        "n_unique": lambda c: f'count(DISTINCT "{c}")',
    }, chunk=chunk, progress_label=f"{table} missing_value_report" if len(cs) > 1000 else None)

    n_missing = agg["n_missing"].reindex(cs)
    report = pd.DataFrame({
        "column": cs,
        "dtype": [str(dtypes.get(c, "")) for c in cs],
        "n_missing": n_missing.values,
        "pct_missing": n_missing.values / n * 100,
        "n_unique": agg["n_unique"].reindex(cs).values,
    }).sort_values("pct_missing", ascending=False).reset_index(drop=True)
    return report


def plot_missing_bar(con, table: str, out_path: Path, top_n: int = MAX_BAR_COLS,
                      report: Optional[pd.DataFrame] = None) -> Optional[Path]:
    report = report if report is not None else missing_value_report(con, table)
    report = report[report["pct_missing"] > 0].head(top_n)
    if report.empty:
        print(f"    (no missing values -- skipping missing-value bar chart for {table})")
        return None

    fig, ax = plt.subplots(figsize=(9, max(3, 0.28 * len(report))))
    y = np.arange(len(report))[::-1]
    ax.barh(y, report["pct_missing"], color="#B33A3A", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(report["column"], fontsize=8)
    ax.set_xlabel("% missing")
    ax.set_title(f"{table} — missing values by column (top {len(report)})")
    ax.set_xlim(0, 100)
    for yi, v in zip(y, report["pct_missing"]):
        ax.text(v + 1, yi, f"{v:.1f}%", va="center", fontsize=8)
    style_axis(ax, grid_axis="x")

    _safe_savefig(fig, out_path)
    return out_path


def plot_missing_matrix(con, table: str, out_path: Path,
                         max_cols: int = MAX_MATRIX_COLS, max_rows: int = MAX_MATRIX_ROWS,
                         report: Optional[pd.DataFrame] = None) -> Optional[Path]:
    if msno is None:
        print("    [SKIPPED] missingno not available -- missing matrix disabled")
        return None

    report = report if report is not None else missing_value_report(con, table)
    partial = report[(report["pct_missing"] > 0) & (report["pct_missing"] < 100)]
    fully_missing = report[report["pct_missing"] == 100]["column"].tolist()

    if partial.empty:
        if fully_missing:
            print(f"    (only fully-missing columns found ({len(fully_missing)}) -- no partial-missingness pattern for {table})")
        else:
            print(f"    (no missing values -- skipping missing matrix for {table})")
        return None

    top_cols = partial.head(max_cols)["column"].tolist()
    sub = fetch_sample(con, table, cols=top_cols, max_rows=max_rows)

    ax = msno.matrix(sub, figsize=(min(14, 0.25 * len(top_cols) + 4), 6), fontsize=8, sparkline=False)
    subtitle = f"{table} — nullity matrix ({len(top_cols)} cols shown, {len(sub):,} rows sampled)"
    if fully_missing:
        shown = ", ".join(fully_missing[:5]) + ("..." if len(fully_missing) > 5 else "")
        subtitle += f"\n({len(fully_missing)} fully-missing column(s) omitted: {shown})"
    ax.set_title(subtitle, fontsize=11)
    _safe_savefig(ax.figure, out_path)
    return out_path


def plot_missing_heatmap(con, table: str, out_path: Path, max_cols: int = MAX_HEATMAP_COLS,
                          report: Optional[pd.DataFrame] = None) -> Optional[Path]:
    if msno is None:
        print("    [SKIPPED] missingno not available -- nullity heatmap disabled")
        return None

    report = report if report is not None else missing_value_report(con, table)
    partial = report[(report["pct_missing"] > 0) & (report["pct_missing"] < 100)]
    if len(partial) < 2:
        print(f"    (fewer than 2 partially-missing columns -- skipping nullity heatmap for {table})")
        return None

    top_cols = partial.head(max_cols)["column"].tolist()
    sub = fetch_sample(con, table, cols=top_cols, max_rows=MAX_MATRIX_ROWS)
    ax = msno.heatmap(sub, figsize=(min(12, 0.35 * len(top_cols) + 3), min(12, 0.35 * len(top_cols) + 3)), fontsize=8)
    ax.set_title(f"{table} — nullity correlation heatmap ({len(top_cols)} cols)")
    _safe_savefig(ax.figure, out_path)
    return out_path


def plot_missing_dendrogram(con, table: str, out_path: Path, max_cols: int = MAX_DENDROGRAM_COLS,
                             report: Optional[pd.DataFrame] = None) -> Optional[Path]:
    if msno is None:
        print("    [SKIPPED] missingno not available -- dendrogram disabled")
        return None

    report = report if report is not None else missing_value_report(con, table)
    partial = report[report["pct_missing"] > 0]
    if len(partial) < 2:
        print(f"    (fewer than 2 columns with missingness -- skipping dendrogram for {table})")
        return None

    top_cols = partial.head(max_cols)["column"].tolist()
    sub = fetch_sample(con, table, cols=top_cols, max_rows=MAX_MATRIX_ROWS)
    ax = msno.dendrogram(sub, figsize=(min(14, 0.3 * len(top_cols) + 4), 6), fontsize=8)
    ax.set_title(f"{table} — nullity dendrogram ({len(top_cols)} cols)")
    _safe_savefig(ax.figure, out_path)
    return out_path


# ============================================================
# 4. Missingness randomness testing
# ============================================================

def little_mcar_test(con, table: str, max_cols: int = MAX_MCAR_COLS, max_rows: int = MAX_MCAR_ROWS,
                      report: Optional[pd.DataFrame] = None) -> Optional[dict]:
    try:
        from pyampute.exploration.mcar_statistical_tests import MCARTest  # type: ignore[import-not-found]
    except ImportError:
        print("    [SKIPPED] pyampute not available -- Little's MCAR test disabled")
        return None

    numeric = numeric_cols(con, table)
    if not numeric:
        print("    (no numeric columns -- Little's MCAR test not applicable)")
        return None

    report = report if report is not None else missing_value_report(con, table)
    report = report[report["column"].isin(numeric)]
    partial = report[(report["pct_missing"] > 0) & (report["pct_missing"] < 100)]
    if len(partial) < 2:
        print("    (fewer than 2 numeric columns with partial missingness -- Little's MCAR test not applicable)")
        return None

    cols = partial.head(max_cols)["column"].tolist()
    sub = fetch_sample(con, table, cols=cols, max_rows=max_rows)

    try:
        p_value = MCARTest(method="little").little_mcar_test(sub)
    except Exception as e:
        log_error(logger, step="little_mcar_test", error=e, n_cols=len(cols), n_rows=len(sub))
        print("    [FAILED] Little's MCAR test errored (often a singular covariance matrix "
              "on wide/correlated data) — see log")
        return None

    return {
        "p_value": p_value,
        "n_cols_tested": len(cols),
        "n_rows_tested": len(sub),
        "columns_tested": cols,
        "reject_mcar_at_0.05": bool(p_value < 0.05),
    }


def missingness_dependency_check(con, table: str, max_focus: int = MAX_DEPENDENCY_FOCUS,
                                  max_compare: int = MAX_DEPENDENCY_COMPARE, alpha: float = 0.05,
                                  max_rows: int = MAX_MCAR_ROWS,
                                  report: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    if stats is None:
        print("    [SKIPPED] scipy not available -- dependency check disabled")
        return pd.DataFrame()

    report = report if report is not None else missing_value_report(con, table)
    focus_candidates = report[(report["pct_missing"] > 0) & (report["pct_missing"] < 100)]
    focus_cols = focus_candidates.head(max_focus)["column"].tolist()
    if not focus_cols:
        return pd.DataFrame()

    numeric = numeric_cols(con, table)
    compare_cols = [c for c in numeric if c not in focus_cols][:max_compare]
    if not compare_cols:
        return pd.DataFrame()

    sub = fetch_sample(con, table, cols=list(dict.fromkeys(focus_cols + compare_cols)), max_rows=max_rows)

    results = []
    for focus_col in focus_cols:
        is_missing = sub[focus_col].isnull()
        if is_missing.sum() < 5 or (~is_missing).sum() < 5:
            continue

        for compare_col in compare_cols:
            vals = sub[compare_col]
            group_missing = vals[is_missing].dropna()
            group_present = vals[~is_missing].dropna()
            if len(group_missing) < 5 or len(group_present) < 5:
                continue
            try:
                _, p = stats.mannwhitneyu(group_missing, group_present, alternative="two-sided")
            except Exception:
                continue

            results.append({
                "focus_col_missing": focus_col,
                "compare_col": compare_col,
                "p_value": p,
                "significant": bool(p < alpha),
                "n_missing_group": len(group_missing),
                "n_present_group": len(group_present),
            })

    result_df = pd.DataFrame(results)
    if not result_df.empty:
        result_df = result_df.sort_values("p_value")
    return result_df


def classify_missingness_mechanism(mcar_result: Optional[dict], dependency_result: pd.DataFrame) -> str:
    n_significant = int(dependency_result["significant"].sum()) if not dependency_result.empty else 0
    n_tested = len(dependency_result)

    if mcar_result is None and n_tested == 0:
        return "Not tested (insufficient missingness or numeric columns)."

    lines = []
    if mcar_result is not None:
        if mcar_result["reject_mcar_at_0.05"]:
            lines.append(
                f"Little's MCAR test REJECTS MCAR (p={mcar_result['p_value']:.4f} < 0.05) "
                f"on {mcar_result['n_cols_tested']} numeric columns — missingness likely depends "
                f"on other variables (MAR) or on the missing values themselves (MNAR)."
            )
        else:
            lines.append(
                f"Little's MCAR test does NOT reject MCAR (p={mcar_result['p_value']:.4f}) "
                f"on {mcar_result['n_cols_tested']} numeric columns — no detected evidence "
                f"against complete randomness, though this doesn't prove MCAR."
            )
    else:
        lines.append("Little's MCAR test: not run (see console/log for reason).")

    if n_tested > 0:
        lines.append(
            f"Dependency check: {n_significant} / {n_tested} focus-column x compare-column pairs "
            f"showed a significant distribution difference (p<0.05) between missing and present rows."
        )
        if n_significant > 0:
            lines.append("  -> Localized evidence against MCAR for at least some columns; see the dependency report for which ones.")
    else:
        lines.append("Dependency check: not run (insufficient missingness or numeric columns).")

    return " ".join(lines)


# ============================================================
# 5. Distribution / scale analysis
# ============================================================

def numeric_summary(con, table: str, chunk: int = MAX_MULTI_AGG_CHUNK, cols: Optional[list] = None) -> pd.DataFrame:
    cs = cols if cols is not None else numeric_cols(con, table)
    if not cs:
        return pd.DataFrame()

    label = f"{table} numeric_summary" if len(cs) > 1000 else None
    agg = chunked_multi_agg(con, table, cs, {
        "count": lambda c: f'count("{c}")',
        "mean": lambda c: f'avg(CAST("{c}" AS DOUBLE))',
        "std": lambda c: f'stddev_samp(CAST("{c}" AS DOUBLE))',
        "min": lambda c: f'min(CAST("{c}" AS DOUBLE))',
        "q1": lambda c: f'quantile_cont("{c}", 0.25)',
        "q2": lambda c: f'quantile_cont("{c}", 0.5)',
        "q3": lambda c: f'quantile_cont("{c}", 0.75)',
        "max": lambda c: f'max(CAST("{c}" AS DOUBLE))',
        "skew": lambda c: f'skewness(CAST("{c}" AS DOUBLE))',
        "kurt": lambda c: f'kurtosis(CAST("{c}" AS DOUBLE))',
    }, chunk=chunk, progress_label=label)

    q1, q3 = agg["q1"].reindex(cs), agg["q3"].reindex(cs)
    iqr = q3 - q1
    lower, upper = (q1 - 1.5 * iqr), (q3 + 1.5 * iqr)

    def outlier_expr(c):
        lo, hi = lower.get(c), upper.get(c)
        if pd.isna(lo) or pd.isna(hi):
            return "NULL"
        return f'sum(CASE WHEN "{c}" < {lo} OR "{c}" > {hi} THEN 1 ELSE 0 END)'

    outliers = chunked_multi_agg(con, table, cs, {"n_outliers": outlier_expr}, chunk=chunk,
                                  progress_label=f"{table} outliers" if len(cs) > 1000 else None)

    summary = pd.DataFrame({
        "column": cs,
        "count": agg["count"].reindex(cs).values,
        "mean": agg["mean"].reindex(cs).values,
        "std": agg["std"].reindex(cs).values,
        "min": agg["min"].reindex(cs).values,
        "25%": q1.values, "50%": agg["q2"].reindex(cs).values, "75%": q3.values,
        "max": agg["max"].reindex(cs).values,
        "skew": agg["skew"].reindex(cs).values,
        "kurtosis": agg["kurt"].reindex(cs).values,
        "iqr": iqr.values,
        "n_outliers_iqr": outliers["n_outliers"].reindex(cs).values,
    })
    summary["pct_outliers_iqr"] = summary["n_outliers_iqr"] / summary["count"] * 100
    return summary.dropna(subset=["mean", "std"])


def detect_scale_profile(summary: pd.DataFrame) -> pd.DataFrame:
    def tag(row):
        if row["min"] >= 0 and row["max"] <= 1:
            return "binary" if row["min"] == 0 and row["max"] == 1 and row["std"] < 0.5 else "proportion_0_1"
        if abs(row["mean"]) < 0.5 and 0.7 <= row["std"] <= 1.3:
            return "standardized"
        if row["min"] >= -1 and row["max"] < 20 and row["max"] > row["min"]:
            return "log_like_small_range"
        if row["min"] >= 0 and row["max"] > 1000:
            return "large_scale_count"
        return "unknown"

    out = summary.copy()
    out["scale_profile"] = out.apply(tag, axis=1)
    return out[["column", "min", "max", "mean", "std", "scale_profile"]]


def _pick_top_variance_cols(con, table: str, max_cols: int, summary: Optional[pd.DataFrame] = None) -> list:
    summary = summary if summary is not None else numeric_summary(con, table)
    if summary.empty:
        return []
    return (
        summary.assign(var=lambda d: d["std"] ** 2)
        .sort_values("var", ascending=False)
        .head(max_cols)["column"].tolist()
    )


def plot_distribution_grid(con, table: str, out_path: Path, cols: Optional[list] = None,
                            max_cols: int = MAX_DIST_COLS, summary: Optional[pd.DataFrame] = None) -> Optional[Path]:
    cols = cols[:max_cols] if cols is not None else _pick_top_variance_cols(con, table, max_cols, summary=summary)
    if not cols:
        print(f"    (no numeric columns -- skipping distribution grid for {table})")
        return None

    sub = fetch_sample(con, table, cols=cols, max_rows=MAX_MCAR_ROWS)
    ncols = min(4, len(cols))
    nrows = int(np.ceil(len(cols) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for ax, col in zip(axes, cols):
        ax.hist(sub[col].dropna(), bins=30, color="#4C72B0", zorder=3)
        ax.set_title(col, fontsize=9)
        style_axis(ax)
    for ax in axes[len(cols):]:
        ax.axis("off")

    fig.suptitle(f"{table} — distributions ({len(cols)} numeric columns shown)", y=1.02)
    _safe_savefig(fig, out_path)
    return out_path


def plot_boxplot_grid(con, table: str, out_path: Path, cols: Optional[list] = None,
                       max_cols: int = MAX_DIST_COLS, summary: Optional[pd.DataFrame] = None) -> Optional[Path]:
    cols = cols[:max_cols] if cols is not None else _pick_top_variance_cols(con, table, max_cols, summary=summary)
    if not cols:
        print(f"    (no numeric columns -- skipping boxplot grid for {table})")
        return None

    sub = fetch_sample(con, table, cols=cols, max_rows=MAX_MCAR_ROWS)
    ncols = min(4, len(cols))
    nrows = int(np.ceil(len(cols) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for ax, col in zip(axes, cols):
        ax.boxplot(sub[col].dropna(), vert=True, patch_artist=True,
                    boxprops=dict(facecolor="#DD8452"), medianprops=dict(color="black"))
        ax.set_title(col, fontsize=9)
        ax.set_xticks([])
        style_axis(ax)
    for ax in axes[len(cols):]:
        ax.axis("off")

    fig.suptitle(f"{table} — boxplots ({len(cols)} numeric columns shown)", y=1.02)
    _safe_savefig(fig, out_path)
    return out_path


def plot_scale_comparison(con, table: str, out_path: Path, max_cols: int = MAX_SCALE_COLS,
                           summary: Optional[pd.DataFrame] = None) -> Optional[Path]:
    summary = summary if summary is not None else numeric_summary(con, table)
    if summary.empty:
        print(f"    (no numeric columns -- skipping scale comparison for {table})")
        return None

    cols = summary.sort_values("50%").head(max_cols)["column"].tolist()
    sub = fetch_sample(con, table, cols=cols, max_rows=MAX_MCAR_ROWS)

    fig, ax = plt.subplots(figsize=(9, max(3, 0.3 * len(cols))))
    data = [sub[c].dropna() for c in cols]
    ax.boxplot(data, vert=False, patch_artist=True, labels=cols,
               boxprops=dict(facecolor="#4C8577"), medianprops=dict(color="black"))
    ax.set_xscale("symlog")
    ax.set_xlabel("value (symlog scale)")
    ax.set_title(f"{table} — scale comparison across columns ({len(cols)} shown, sorted by median)")
    ax.tick_params(axis="y", labelsize=8)
    style_axis(ax, grid_axis="x")

    _safe_savefig(fig, out_path)
    return out_path


def plot_correlation_heatmap(con, table: str, out_path: Path, max_cols: int = MAX_CORR_COLS,
                              summary: Optional[pd.DataFrame] = None) -> Optional[Path]:
    cols = _pick_top_variance_cols(con, table, max_cols, summary=summary)
    if len(cols) < 2:
        print(f"    (fewer than 2 numeric columns -- skipping correlation heatmap for {table})")
        return None

    sub = fetch_sample(con, table, cols=cols, max_rows=MAX_MCAR_ROWS)
    corr = sub.corr()

    fig, ax = plt.subplots(figsize=(min(14, 0.4 * len(cols) + 3), min(12, 0.4 * len(cols) + 3)))
    im = ax.imshow(corr, cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(cols)))
    ax.set_yticklabels(cols, fontsize=8)
    if len(cols) <= 15:
        for i in range(len(cols)):
            for j in range(len(cols)):
                ax.text(j, i, f"{corr.iloc[i, j]:.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, shrink=0.8)
    ax.set_title(f"{table} — correlation heatmap ({len(cols)} numeric columns shown)")

    _safe_savefig(fig, out_path)
    return out_path


# ============================================================
# 6. Cross-table summary plots -- span ALL tables at once
# ============================================================

def column_moments(con, table: str, chunk: int = MAX_AGG_CHUNK, cols: Optional[list] = None) -> pd.DataFrame:
    cs = cols if cols is not None else numeric_cols(con, table)
    if not cs:
        return pd.DataFrame(columns=["column", "mean", "sd"])
    mean = chunked_agg(con, table, cs, lambda c: f'avg(CAST("{c}" AS DOUBLE))', chunk)
    sd = chunked_agg(con, table, cs, lambda c: f'stddev_samp(CAST("{c}" AS DOUBLE))', chunk)
    d = pd.DataFrame({"column": mean.index, "mean": mean.values, "sd": sd.values})
    return d.replace([np.inf, -np.inf], np.nan).dropna().query("sd > 0")


def plot_gaussian_per_column(con, tables: list, out_path: Path, max_cols: int = MAX_GAUSSIAN_COLS,
                              ncols: int = 3, grid: int = 400, seed: int = 7,
                              table_cols: Optional[dict] = None) -> tuple:
    tables = [table for table in tables if table in GAUSSIAN_FACET_TABLES]
    nr = -(-len(tables) // ncols)
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    })
    fig, axes = plt.subplots(nr, ncols, figsize=(11.69, 8.27), squeeze=False)
    axes = axes.flatten()
    allm = []

    for ax, t in zip(axes, tables):
        m = column_moments(con, t, cols=(table_cols or {}).get(t))
        if m.empty:
            ax.axis("off")
            continue
        m["table"] = t
        allm.append(m)
        d = m.sample(max_cols, random_state=seed) if len(m) > max_cols else m

        lo = np.percentile(d["mean"] - 3 * d["sd"], 1)
        hi = np.percentile(d["mean"] + 3 * d["sd"], 99)
        xs = np.linspace(lo, hi, grid)
        cmap = plt.get_cmap("viridis")
        rank = d["mean"].rank(pct=True).to_numpy()

        for (mu, sd), r in zip(d[["mean", "sd"]].to_numpy(), rank):
            ax.plot(xs, np.exp(-0.5 * ((xs - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi)),
                    lw=0.55, alpha=0.32, color=cmap(r), zorder=3)

        display_name = t.replace("_", " ").title()
        ax.set_title(f"{display_name}   ({len(d):,} of {len(m):,} columns)",
                 fontsize=11.5, fontweight="bold", loc="left", pad=7)
        ax.set_xlabel("Value", fontsize=12)
        ax.set_ylabel("Density", fontsize=12)
        ax.tick_params(axis="both", labelsize=10)
        style_axis(ax, grid_axis="both")

    for ax in axes[len(tables):]:
        ax.axis("off")

    fig.suptitle("Distribution and Scale of Numeric Variables Across Datasets",
                 fontsize=18, fontweight="bold", x=0.5, ha="center", y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.94], h_pad=1.8, w_pad=1.4)

    moments = pd.concat(allm, ignore_index=True) if allm else pd.DataFrame()
    _safe_savefig(fig, out_path)
    return fig, moments


MISS_BINS = [-0.001, 0.0, 10, 25, 50, 75, 95, 100.001]
MISS_LABELS = ["complete", "<10%", "10–25%", "25–50%", "50–75%", "75–95%", "≥95%"]
MISS_COLORS = ["#2E7D5B", "#7FB069", "#C9CE6E", "#EDC26B", "#E39A5C", "#D9704E", "#B33A3A"]


def missing_pct(con, table: str, chunk: int = MAX_AGG_CHUNK, cols: Optional[list] = None) -> pd.Series:
    cs = cols if cols is not None else cols_of(con, table)
    total = row_count(con, table)
    if total == 0 or not cs:
        return pd.Series(dtype=float)
    return chunked_agg(con, table, cs,
                        lambda c: f'100.0*sum(CASE WHEN "{c}" IS NULL THEN 1 ELSE 0 END)/{total}', chunk)


def plot_missingness(con, tables: list, out_path: Path, ncols: int = 4, table_cols: Optional[dict] = None):
    tables = [table for table in MISSINGNESS_FACET_TABLES if table in tables]
    data = {}
    for t in tables:
        s = missing_pct(con, t, cols=(table_cols or {}).get(t))
        data[t] = None if s.empty else (
            pd.cut(s, bins=MISS_BINS, labels=MISS_LABELS).value_counts().reindex(MISS_LABELS, fill_value=0)
        )

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    })
    nrows = -(-len(tables) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(11.69, 8.27), squeeze=False)
    axes = axes.flatten()
    ypos = np.arange(len(MISS_LABELS))[::-1]

    for ax, table in zip(axes, tables):
        counts = data[table]
        if counts is None:
            ax.axis("off")
            continue
        n_total = int(counts.sum())
        ax.barh(ypos, counts.values, color=MISS_COLORS, height=0.7,
                edgecolor="white", linewidth=0.3, zorder=3)
        ax.set_yticks(ypos)
        ax.set_yticklabels(MISS_LABELS, fontsize=8.5)
        ax.set_title(f"{table.replace('_', ' ').title()}  ({n_total:,} columns)",
                     fontsize=10.5, fontweight="bold", loc="left", pad=5)
        ax.set_xlim(0, max(counts.max() * 1.3, 1))
        ax.tick_params(axis="x", labelsize=8)
        style_axis(ax)

    for ax in axes[len(tables):]:
        ax.axis("off")

    fig.legend(handles=[Patch(facecolor=c, label=l) for c, l in zip(MISS_COLORS, MISS_LABELS)],
               loc="lower center", ncol=7, frameon=False, fontsize=10,
               bbox_to_anchor=(0.5, -0.055))
    fig.suptitle("Distribution of Column Missingness Across Selected Tables",
                 fontsize=17, fontweight="bold", x=0.5, ha="center", y=0.985)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.925, bottom=0.34,
                        wspace=0.42, hspace=0.85)

    _safe_savefig(fig, out_path)
    return fig


# ============================================================
# 7. Per-table and full-pipeline runners
# ============================================================

def run_eda_on_table(con, name: str, out_dir: Path, protein_coding_genes: Optional[list] = None,
                      force: bool = False) -> dict:
    table_dir = Path(out_dir) / name
    table_dir.mkdir(parents=True, exist_ok=True)
    n_rows = row_count(con, name)
    all_cols = cols_of(con, name)

    restrict_cols = None
    if protein_coding_genes is not None and name in GENE_MATRIX_TABLES:
        pc_set = set(protein_coding_genes)
        restrict_cols = [c for c in all_cols if c in pc_set]

    n_cols_analyzed = len(restrict_cols) if restrict_cols is not None else len(all_cols)
    header = f"{name}  ({n_rows:,} rows x {len(all_cols):,} cols"
    header += f", {n_cols_analyzed:,} protein-coding analyzed)" if restrict_cols is not None else ")"
    print(f"\n{'=' * 70}\n{header}\n{'=' * 70}")

    summary_row = {
        "dataset": name,
        "rows": n_rows,
        "cols": len(all_cols),
        "cols_analyzed": n_cols_analyzed,
        "overall_pct_missing": None,
        "mcar_p_value": None,
        "mcar_reject_at_0.05": None,
        "missingness_verdict": None,
    }

    def step(step_name, fn, *args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except Exception as e:
            log_error(logger, step=f"run_eda_on_table:{step_name}", error=e, dataset=name)
            print(f"  [SKIPPED] {step_name} failed — see log")
            return None

    report_path = table_dir / "missing_value_report.csv"
    report = step("missing_value_report", cached_csv, report_path,
                  lambda: missing_value_report(con, name, cols=restrict_cols), force)
    if report is not None:
        summary_row["overall_pct_missing"] = (
            report["n_missing"].sum() / (n_rows * n_cols_analyzed) * 100
        ) if n_rows and n_cols_analyzed else 0.0

    step("plot_missing_bar", plot_missing_bar, con, name, table_dir / "missing_bar.png", report=report)
    step("plot_missing_matrix", plot_missing_matrix, con, name, table_dir / "missing_matrix.png", report=report)
    step("plot_missing_heatmap", plot_missing_heatmap, con, name, table_dir / "missing_nullity_heatmap.png", report=report)
    step("plot_missing_dendrogram", plot_missing_dendrogram, con, name, table_dir / "missing_dendrogram.png", report=report)

    mcar_result = step("little_mcar_test", little_mcar_test, con, name, report=report)
    dependency_result = step("missingness_dependency_check", missingness_dependency_check, con, name, report=report)
    if dependency_result is None:
        dependency_result = pd.DataFrame()
    if not dependency_result.empty:
        dependency_result.to_csv(table_dir / "missingness_dependency_check.csv", index=False)

    verdict = step("classify_missingness_mechanism", classify_missingness_mechanism, mcar_result, dependency_result)
    if verdict is not None:
        (table_dir / "missingness_verdict.txt").write_text(verdict)
        print(f"  Missingness verdict: {verdict}")
        summary_row["missingness_verdict"] = verdict
    if mcar_result is not None:
        summary_row["mcar_p_value"] = mcar_result["p_value"]
        summary_row["mcar_reject_at_0.05"] = mcar_result["reject_mcar_at_0.05"]

    summary_path = table_dir / "numeric_summary.csv"
    num_summary = step("numeric_summary", cached_csv, summary_path,
                       lambda: numeric_summary(con, name, cols=restrict_cols), force)
    if num_summary is not None and not num_summary.empty:
        scale_profile = step("detect_scale_profile", detect_scale_profile, num_summary)
        if scale_profile is not None:
            scale_profile.to_csv(table_dir / "scale_profile.csv", index=False)

    step("plot_distribution_grid", plot_distribution_grid, con, name, table_dir / "distributions.png", summary=num_summary)
    step("plot_boxplot_grid", plot_boxplot_grid, con, name, table_dir / "boxplots.png", summary=num_summary)
    step("plot_scale_comparison", plot_scale_comparison, con, name, table_dir / "scale_comparison.png", summary=num_summary)
    step("plot_correlation_heatmap", plot_correlation_heatmap, con, name, table_dir / "correlation_heatmap.png", summary=num_summary)

    print(f"  Saved EDA output to {table_dir}")
    return summary_row


def run_eda_pipeline(con, table_names: list, out_dir: Path, force: bool = False) -> pd.DataFrame:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    protein_coding_genes = get_protein_coding_gene_ids(con)
    if protein_coding_genes is not None:
        print(f"Restricting {', '.join(sorted(GENE_MATRIX_TABLES))} to "
              f"{len(protein_coding_genes):,} protein-coding genes from gene_roster")

    summaries = []
    for name in table_names:
        try:
            summaries.append(run_eda_on_table(con, name, out_dir,
                                               protein_coding_genes=protein_coding_genes, force=force))
        except Exception as e:
            log_error(logger, step="run_eda_pipeline", error=e, dataset=name)
            print(f"[SKIPPED] {name}: EDA failed entirely — see log")

    overview = pd.DataFrame(summaries)
    overview.to_csv(out_dir / "_overview.csv", index=False)

    if not overview.empty:
        fig, ax = plt.subplots(figsize=(9, max(3, 0.35 * len(overview))))
        plot_df = overview.sort_values("overall_pct_missing", ascending=False)
        y = np.arange(len(plot_df))[::-1]
        ax.barh(y, plot_df["overall_pct_missing"], color="#B33A3A", zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(plot_df["dataset"], fontsize=8)
        ax.set_xlabel("overall % missing")
        ax.set_title("Overall missing % by dataset")
        ax.set_xlim(0, 100)
        for yi, v in zip(y, plot_df["overall_pct_missing"]):
            ax.text(v + 1, yi, f"{v:.1f}%", va="center", fontsize=8)
        style_axis(ax, grid_axis="x")
        _safe_savefig(fig, out_dir / "_overview_missing_by_dataset.png")

    table_cols = None
    if protein_coding_genes is not None:
        pc_set = set(protein_coding_genes)
        table_cols = {
            t: [c for c in cols_of(con, t) if c in pc_set]
            for t in GENE_MATRIX_TABLES if t in table_names
        }

    print("\nBuilding cross-table summary plots...")
    try:
        _, moments = plot_gaussian_per_column(con, table_names, out_dir / "_summary_gaussian_per_column.png",
                                               table_cols=table_cols)
        if not moments.empty:
            moments.to_csv(out_dir / "_summary_column_moments.csv", index=False)
    except Exception as e:
        log_error(logger, step="plot_gaussian_per_column", error=e)
        print("  [SKIPPED] plot_gaussian_per_column failed — see log")

    try:
        plot_missingness(con, table_names, out_dir / "_summary_missingness_by_table.png", table_cols=table_cols)
    except Exception as e:
        log_error(logger, step="plot_missingness", error=e)
        print("  [SKIPPED] plot_missingness failed — see log")

    print(f"\nEDA complete. Overview saved to {out_dir / '_overview.csv'}")
    return overview
