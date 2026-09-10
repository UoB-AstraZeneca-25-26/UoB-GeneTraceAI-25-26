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

Further design notes
--------------------
* Optional dependencies degrade, they don't break. scipy, missingno and
  pyampute are each imported defensively; when one is absent the
  analyses that need it print a skip message and return None, and the
  rest of the run proceeds.
* Nothing raises out of a table's run. Every step inside
  :func:`run_eda_on_table` goes through a local wrapper that logs and
  returns None, so one failed plot cannot cost you the other twenty.
* Every ``MAX_*`` constant is a legibility or tractability bound, not a
  statistical choice. Plots cap the columns they show because a
  200-column heatmap is unreadable; the MCAR test caps them because the
  covariance matrix goes singular. Both are stated where they apply.
* All figures are written through :func:`_safe_savefig`, which closes
  the figure — matplotlib holds figures open otherwise, and a run over
  every table would exhaust memory.
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

#: Display caps for the per-table plots. Each bounds how many columns a
#: figure can show before it stops being readable — not a statistical
#: cut. Columns are chosen by the relevant ranking (missingness for the
#: nullity plots, variance for the distribution ones).
MAX_BAR_COLS = 40
MAX_MATRIX_COLS = 60
MAX_MATRIX_ROWS = 500
MAX_HEATMAP_COLS = 40
MAX_DENDROGRAM_COLS = 60

#: Caps for the missingness randomness tests. MAX_MCAR_COLS is a
#: tractability bound rather than a display one: Little's test inverts a
#: covariance matrix, which goes singular on wide correlated data.
MAX_MCAR_COLS = 30
MAX_MCAR_ROWS = 5000
MAX_DEPENDENCY_FOCUS = 10
MAX_DEPENDENCY_COMPARE = 20

MAX_DIST_COLS = 12
MAX_SCALE_COLS = 20
MAX_CORR_COLS = 40
MAX_GAUSSIAN_COLS = 250

#: Columns per DuckDB query. The multi-agg chunk is smaller because each
#: column contributes one expression per statistic, so a chunk of 100
#: columns x 10 statistics is already a 1000-expression SELECT.
MAX_AGG_CHUNK = 250
MAX_MULTI_AGG_CHUNK = 100

MAX_RAW_FACET_VALUES = 20_000

#: Tables included in the cross-table density figure — the ones whose
#: numeric columns are genuine measurements rather than identifiers.
GAUSSIAN_FACET_TABLES = {
    "depmap_expr", "geo_expr", "hpa_rna", "proteomics", "metabolomics", "mirna",
}

#: Tables included in the cross-table missingness figure, in the order
#: they are laid out. A tuple rather than a set, since order matters here.
MISSINGNESS_FACET_TABLES = (
    "depmap_expr", "geo_expr", "hpa_rna", "mutations", "fusions", "proteomics",
    "metabolomics", "mirna", "signatures", "sample_info", "cellosaurus", "hpa_desc",
    "depmap_profiles", "geo_info",
)

#: Tables whose columns are gene IDs, and which are therefore restricted
#: to gene_roster's protein-coding set before anything is computed.
GENE_MATRIX_TABLES = {"depmap_expr", "geo_expr"}

from concurrent.futures import ThreadPoolExecutor, as_completed
import os

MAX_WORKERS = min(8, os.cpu_count() or 4)


def _safe_savefig(fig, out_path: Path) -> None:
    """
    Write a figure to disk and close it.

    Creates the parent directory if needed, tightens the layout, and
    saves with a tight bounding box. Closing matters: matplotlib keeps
    every figure alive until closed, and a full run produces dozens.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to write. Closed afterwards, so the caller must not reuse
        it.
    out_path : Path or str
        Destination file. The extension determines the format.

    Returns
    -------
    None
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def style_axis(ax, grid_axis: str = "y") -> None:
    """
    Apply the project's house style to an axis.

    Drops the top and right spines and puts a faint grid behind the data.
    Applied to every axis in this file, so the figures stay visually
    consistent with the plots produced elsewhere in the project.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to style, modified in place.
    grid_axis : str, optional
        Which axis gets grid lines: ``"y"`` (default), ``"x"`` or
        ``"both"``. Choose the one the values run along.

    Returns
    -------
    None
    """
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis=grid_axis, alpha=0.3, zorder=0)
    ax.set_axisbelow(True)


def cached_csv(path: Path, compute_fn, force: bool = False) -> pd.DataFrame:
    """
    Load a DataFrame from `path` if it already exists (and force is
    False); otherwise compute it via compute_fn(), save it to `path`,
    and return it. Delete the CSV, or pass force=True, to recompute.

    This is what makes a rerun cheap. The two expensive per-table
    statistics — the missing-value report and the numeric summary — each
    cost one pass over every column, which on a 54k-column table is the
    bulk of the runtime.

    Parameters
    ----------
    path : Path or str
        Cache file. Its parent is created on write.
    compute_fn : callable
        Zero-argument callable returning the DataFrame. Only called on a
        miss, so an expensive computation is skipped entirely on a hit.
    force : bool, optional
        Ignore any existing cache and recompute. Default False.

    Returns
    -------
    pandas.DataFrame
        Either the cached frame or the freshly computed one.

    Notes
    -----
    The cache key is the path alone — nothing about the inputs is
    recorded. If the underlying table or the protein-coding gene set
    changes, the stale cache is returned silently; pass ``force=True``
    after any upstream change.

    A cached frame comes back through the CSV round trip, so dtypes are
    re-inferred rather than preserved.
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

#: DuckDB base type names treated as numeric. Matched against the first
#: word of the declared type, so parameterised types like
#: ``DECIMAL(18,3)`` are caught.
_NUMERIC_DUCKDB_TYPES = {
    "TINYINT", "SMALLINT", "INTEGER", "BIGINT", "HUGEINT",
    "UTINYINT", "USMALLINT", "UINTEGER", "UBIGINT",
    "FLOAT", "DOUBLE", "DECIMAL", "REAL",
}


def cols_of(con, table: str) -> list:
    """
    Return every column name in a table.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.

    Returns
    -------
    list of str
        Column names in declaration order.
    """
    return con.execute(f'DESCRIBE "{table}"').fetchdf()["column_name"].tolist()


def numeric_cols(con, table: str) -> list:
    """
    Return the numeric columns of a table, by declared type.

    Selection is on the schema rather than on the values, so a numeric
    column that happens to be all-null is still included — which matters,
    since that is exactly the kind of column the missingness analysis is
    looking for.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.

    Returns
    -------
    list of str
        Names of columns whose base type is in
        :data:`_NUMERIC_DUCKDB_TYPES`.
    """
    desc = con.execute(f'DESCRIBE "{table}"').fetchdf()
    base_type = desc["column_type"].str.upper().str.extract(r"^(\w+)")[0]
    return desc.loc[base_type.isin(_NUMERIC_DUCKDB_TYPES), "column_name"].tolist()


def row_count(con, table: str) -> int:
    """
    Return the number of rows in a table.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.

    Returns
    -------
    int
        Row count. Used as the denominator for every missingness
        percentage, and to decide whether sampling is needed.
    """
    return con.execute(f'SELECT count(*) FROM "{table}"').fetchone()[0]


def get_protein_coding_gene_ids(con) -> Optional[list]:
    """
    Load the protein-coding gene IDs that define the gene-matrix restriction.

    Every gene_id in ``gene_roster`` counts, since the roster is already
    built protein-coding-only upstream.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.

    Returns
    -------
    list of str or None
        Distinct gene IDs, or None if ``gene_roster`` is unavailable.

    Notes
    -----
    None means no restriction: the gene-matrix tables are then profiled
    over all their columns. That is far slower and mixes non-coding genes
    into every statistic, so the warning printed here is worth reading
    rather than scrolling past.
    """
    try:
        return con.execute('SELECT DISTINCT gene_id FROM gene_roster').fetchdf()["gene_id"].tolist()
    except Exception as e:
        log_error(logger, step="get_protein_coding_gene_ids", error=e)
        print("    [WARN] gene_roster not available -- gene-matrix tables will use ALL columns, not just protein-coding")
        return None


def _run_single_chunk(con, table: str, batch: list, expr_fn) -> pd.Series:
    """
    Run one aggregate expression across a batch of columns, in one query.

    Uses its own cursor so the query can run on a worker thread without
    contending for the shared connection.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to derive a cursor from.
    table : str
        Table or view name.
    batch : list of str
        Columns for this query.
    expr_fn : callable
        Takes a column name and returns the SQL aggregate expression for
        it.

    Returns
    -------
    pandas.Series
        One value per column, indexed by column name.
    """
    cur = con.cursor()
    select = ", ".join(f'{expr_fn(c)} AS "{c}"' for c in batch)
    row = cur.execute(f'SELECT {select} FROM "{table}"').fetchdf().iloc[0]
    cur.close()
    return row


def chunked_agg(con, table: str, cols: list, expr_fn, chunk: int = MAX_AGG_CHUNK,
                 max_workers: int = MAX_WORKERS) -> pd.Series:
    """
    Compute one aggregate across many columns, chunked and run in parallel.

    Columns are split into batches, each batch becomes a single query, and
    the batches run concurrently on their own cursors. Chunking is what
    keeps the SQL statement a workable size on a table with tens of
    thousands of columns.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    cols : list of str
        Columns to aggregate. An empty list returns an empty Series.
    expr_fn : callable
        Takes a column name, returns its SQL aggregate expression.
    chunk : int, optional
        Columns per query. Defaults to :data:`MAX_AGG_CHUNK`.
    max_workers : int, optional
        Concurrent queries. Defaults to :data:`MAX_WORKERS`.

    Returns
    -------
    pandas.Series
        One value per column, indexed by column name.

    Notes
    -----
    Results arrive as batches complete, so the index order does not match
    ``cols`` — reindex if order matters. Where several statistics are
    wanted, :func:`chunked_multi_agg` computes them in one pass instead of
    one pass each.
    """
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
    """
    Run every requested aggregate across a batch of columns, in one query.

    Each output column is aliased ``<stat_name>__<column>``, which is how
    the caller unpacks the flat result row back into per-statistic Series.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to derive a cursor from.
    table : str
        Table or view name.
    batch : list of str
        Columns for this query.
    expr_map : dict
        Statistic name to a callable returning its SQL expression for a
        given column.

    Returns
    -------
    pandas.Series
        One value per (statistic, column) pair, indexed by the composite
        alias.

    Notes
    -----
    The query contains ``len(batch) * len(expr_map)`` expressions, which
    is why the multi-agg chunk size is smaller than the single-agg one.
    """
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
    """
    Compute several aggregates across many columns in one pass per chunk.

    The point of batching every statistic into a single query is that a
    wide table is then scanned once per chunk rather than once per
    statistic — the difference between ten passes over 54k columns and
    one.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    cols : list of str
        Columns to aggregate. Empty returns empty Series per statistic.
    expr_map : dict
        Statistic name to a callable returning its SQL expression. Names
        must not contain ``__``, which separates them from column names
        in the result aliases.
    chunk : int, optional
        Columns per query. Defaults to :data:`MAX_MULTI_AGG_CHUNK`.
    progress_label : str, optional
        When given, prints chunk progress — worth passing for tables
        large enough that the wait is noticeable.
    max_workers : int, optional
        Concurrent queries. Defaults to :data:`MAX_WORKERS`.

    Returns
    -------
    dict
        Statistic name to a :class:`pandas.Series` of its values, indexed
        by column name.

    Notes
    -----
    As with :func:`chunked_agg`, results arrive out of order — reindex
    against ``cols`` where order matters, which the callers do.
    """
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
    """
    Pull a bounded sample of rows into pandas.

    This is the only route by which row-level data enters memory. Plots
    that genuinely need rows — histograms, nullity matrices, correlations
    — go through here; everything else is computed in DuckDB.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    cols : list of str, optional
        Columns to fetch. All columns when omitted, which on a wide table
        is exactly what this function exists to avoid — pass a subset.
    max_rows : int, optional
        Row cap. Sampling only kicks in when the table exceeds it, so a
        small table is returned whole.
    seed : int, optional
        Reservoir sampling seed, so a rerun draws the same rows and the
        figures are reproducible. Default 0.

    Returns
    -------
    pandas.DataFrame
        The sampled rows.
    """
    col_sql = "*" if not cols else ", ".join(f'"{c}"' for c in cols)
    query = f'SELECT {col_sql} FROM "{table}"'
    if max_rows is not None:
        total = row_count(con, table)
        if total > max_rows:
            query += f" USING SAMPLE {max_rows} ROWS (reservoir, {seed})"
    return con.execute(query).fetchdf()


def raw_measurement_cols(con, table: str) -> list:
    """Return numeric columns that are likely measurements, not identifiers.

    Filters out numeric columns whose names suggest they are IDs or
    genomic coordinates — profiling the distribution of a chromosome
    position or a row index says nothing useful.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.

    Returns
    -------
    list of str
        Numeric columns whose names contain none of ``id``, ``index``,
        ``start``, ``end``, ``position``, ``chromosome``, ``coord``.

    Notes
    -----
    A substring match, so it is deliberately blunt: a genuine measurement
    with ``id`` anywhere in its name is excluded too.
    """
    excluded_tokens = ("id", "index", "start", "end", "position", "chromosome", "coord")
    return [
        col for col in numeric_cols(con, table)
        if not any(token in col.lower() for token in excluded_tokens)
    ]


def fetch_raw_value_sample(con, table: str, cols: list,
                           max_values: int = MAX_RAW_FACET_VALUES,
                           seed: int = 7) -> np.ndarray:
    """Sample pooled values from wide raw tables without materializing them.

    Reduces a wide table to a flat array of finite values, sampling in
    both directions: columns are thinned to an evenly spaced subset when
    there are too many, and only a few rows are read per query.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    cols : list of str
        Candidate columns. Empty returns an empty array.
    max_values : int, optional
        Cap on returned values. Defaults to
        :data:`MAX_RAW_FACET_VALUES`.
    seed : int, optional
        Seed for the final down-sample. Default 7.

    Returns
    -------
    numpy.ndarray
        Pooled finite values, at most ``max_values`` of them.

    Notes
    -----
    Rows are taken with ``LIMIT``, not a random sample — this reads the
    first N rows of each column batch. Fine for characterising a value
    scale; not a random sample of the table, so don't read anything
    distributional into it.
    """
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
    """
    Per-column missingness and cardinality for one table.

    Computed entirely in DuckDB, with the null count and distinct count
    batched into one query per chunk. This is the expensive statistic on
    a wide table, which is why it is computed once per table, cached to
    disk, and passed into every downstream plot rather than recomputed.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    chunk : int, optional
        Columns per query. Defaults to :data:`MAX_MULTI_AGG_CHUNK`.
    cols : list of str, optional
        Restrict to these columns — how the gene-matrix tables are
        limited to protein-coding genes. All columns when omitted.

    Returns
    -------
    pandas.DataFrame
        Columns ``column``, ``dtype``, ``n_missing``, ``pct_missing``,
        ``n_unique``, sorted by ``pct_missing`` descending. Empty with
        those columns when the table has no rows or no columns.

    Notes
    -----
    Sorting descending is what makes ``.head(n)`` the right selection in
    every plot below: the most-missing columns are the informative ones.
    Progress is printed only past 1000 columns.
    """
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
    """
    Horizontal bar chart of the most-missing columns.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection. Used only if ``report`` is not supplied.
    table : str
        Table or view name, used in the title.
    out_path : Path
        Destination image file.
    top_n : int, optional
        Columns to show. Defaults to :data:`MAX_BAR_COLS`.
    report : pandas.DataFrame, optional
        Precomputed missing-value report. Always pass this — recomputing
        it per plot is what made wide tables slow.

    Returns
    -------
    Path or None
        The output path, or None when no column has any missing values,
        in which case nothing is written.
    """
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
    """
    missingno nullity matrix over the partially-missing columns.

    Shows where the gaps fall row by row, which reveals whether
    missingness clusters — several columns going null together points at
    a merge or a source that dropped out, not at random loss.

    Only *partially* missing columns are plotted. A fully-missing column
    is a solid bar carrying no pattern, so those are excluded and named
    in the subtitle instead.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    out_path : Path
        Destination image file.
    max_cols, max_rows : int, optional
        Display caps. Default to :data:`MAX_MATRIX_COLS` and
        :data:`MAX_MATRIX_ROWS`.
    report : pandas.DataFrame, optional
        Precomputed missing-value report.

    Returns
    -------
    Path or None
        The output path, or None when missingno is unavailable or no
        column is partially missing.
    """
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
    """
    missingno nullity correlation heatmap.

    Quantifies what the matrix shows visually: how strongly one column's
    nullity predicts another's. A strong pair means the two go missing
    together, which usually traces back to a shared source rather than to
    anything about the values.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    out_path : Path
        Destination image file.
    max_cols : int, optional
        Display cap. Defaults to :data:`MAX_HEATMAP_COLS`.
    report : pandas.DataFrame, optional
        Precomputed missing-value report.

    Returns
    -------
    Path or None
        The output path, or None when missingno is unavailable or fewer
        than two columns are partially missing — a correlation needs a
        pair.
    """
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
    """
    missingno nullity dendrogram.

    Clusters columns by how similarly they go missing. Where the heatmap
    shows pairs, this shows groups — a tight cluster is a set of columns
    that share a fate, which is the clearest signal that they share a
    source.

    Unlike the matrix and heatmap, fully-missing columns are included
    here, since they legitimately cluster together.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    out_path : Path
        Destination image file.
    max_cols : int, optional
        Display cap. Defaults to :data:`MAX_DENDROGRAM_COLS`.
    report : pandas.DataFrame, optional
        Precomputed missing-value report.

    Returns
    -------
    Path or None
        The output path, or None when missingno is unavailable or fewer
        than two columns have any missingness.
    """
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
    """
    Little's test for Missing Completely At Random.

    Tests the joint hypothesis that missingness is unrelated to any
    observed value. Rejection means the data are MAR or MNAR, which
    matters downstream: MCAR is the only mechanism under which
    complete-case analysis is unbiased.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    max_cols : int, optional
        Columns tested. Defaults to :data:`MAX_MCAR_COLS` — a
        tractability bound, since the test inverts a covariance matrix.
    max_rows : int, optional
        Rows sampled. Defaults to :data:`MAX_MCAR_ROWS`.
    report : pandas.DataFrame, optional
        Precomputed missing-value report.

    Returns
    -------
    dict or None
        ``p_value``, ``n_cols_tested``, ``n_rows_tested``,
        ``columns_tested``, ``reject_mcar_at_0.05``. None when pyampute
        is unavailable, when fewer than two numeric columns are partially
        missing, or when the test itself errored.

    Notes
    -----
    Failure is expected on wide correlated data — a singular covariance
    matrix — and is logged and returned as None rather than raised.

    Two limits on how far the result travels: only the ``max_cols``
    most-missing numeric columns are tested, so this is not a statement
    about the table as a whole; and not rejecting MCAR is not evidence
    for it, only absence of evidence against.
    """
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
    """
    Per-column heuristic complement to Little's test.

    For each focus column with partial missingness, splits every compare
    column's values by whether the focus column is null and tests the two
    groups with Mann-Whitney U. A significant difference says missingness
    in the focus column tracks the values of the compare column — the
    localised version of what Little's test asks globally.

    Where Little's test gives one verdict for a whole table, this names
    which columns are implicated, which is what makes it actionable.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    max_focus, max_compare : int, optional
        Caps on the two column sets. Default to
        :data:`MAX_DEPENDENCY_FOCUS` and
        :data:`MAX_DEPENDENCY_COMPARE`.
    alpha : float, optional
        Significance level for the ``significant`` flag. Default 0.05.
    max_rows : int, optional
        Rows sampled. Defaults to :data:`MAX_MCAR_ROWS`.
    report : pandas.DataFrame, optional
        Precomputed missing-value report.

    Returns
    -------
    pandas.DataFrame
        Columns ``focus_col_missing``, ``compare_col``, ``p_value``,
        ``significant``, ``n_missing_group``, ``n_present_group``, sorted
        by p-value. Empty when scipy is unavailable or nothing was
        testable.

    Notes
    -----
    Mann-Whitney is used rather than a t-test because it makes no
    normality assumption, which suits the skewed expression
    distributions here.

    No multiple-comparison correction is applied, and up to 200 pairs are
    tested — expect roughly ten false positives at alpha 0.05 by chance
    alone. Read the count as a signal of degree, not as a set of
    individually confirmed findings. Groups smaller than five
    observations are skipped as untestable.
    """
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
    """
    Turn the two missingness tests into a written verdict.

    Combines the global test and the per-column heuristic into prose
    saved alongside the table's plots, so the finding is readable without
    re-deriving it from the CSVs.

    Parameters
    ----------
    mcar_result : dict or None
        Result from :func:`little_mcar_test`. None when it did not run.
    dependency_result : pandas.DataFrame
        Result from :func:`missingness_dependency_check`. May be empty.

    Returns
    -------
    str
        A sentence or three describing what was found, or a note that
        nothing was tested when neither test ran.

    Notes
    -----
    The wording is deliberately hedged in both directions: a
    non-rejection is reported as no detected evidence rather than as
    support for MCAR, and a rejection does not distinguish MAR from MNAR,
    which no test on observed data can.
    """
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
    """
    Full descriptive statistics per numeric column, computed in DuckDB.

    Ten statistics are batched into one query per chunk, then a second
    pass counts IQR outliers — the second pass is unavoidable, since the
    outlier bounds depend on the quartiles from the first.

    Like the missing-value report, this is the expensive per-table
    statistic: computed once, cached to disk, and passed into every plot
    below via ``summary=``.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    chunk : int, optional
        Columns per query. Defaults to :data:`MAX_MULTI_AGG_CHUNK`.
    cols : list of str, optional
        Restrict to these columns. All numeric columns when omitted.

    Returns
    -------
    pandas.DataFrame
        Columns ``column``, ``count``, ``mean``, ``std``, ``min``,
        ``25%``, ``50%``, ``75%``, ``max``, ``skew``, ``kurtosis``,
        ``iqr``, ``n_outliers_iqr``, ``pct_outliers_iqr``. Empty when the
        table has no numeric columns.

    Notes
    -----
    Outliers use the standard Tukey fences, 1.5 x IQR beyond the
    quartiles. On the heavy-tailed expression data here a high outlier
    percentage is expected and is not itself a data-quality finding.

    Rows with a null mean or standard deviation are dropped, so an
    all-null or constant-null numeric column is absent from the result
    even though it exists in the table. Compare against
    :func:`numeric_cols` if a column seems to have vanished.
    """
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
        """SQL count of values outside this column's Tukey fences, or NULL."""
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
    """
    Label each numeric column with the kind of scale it appears to be on.

    A heuristic read of the range, mean and spread, used to spot columns
    that are not on the scale they are assumed to be — a supposedly
    log-transformed matrix showing a large-scale count profile is the
    finding this is for.

    Parameters
    ----------
    summary : pandas.DataFrame
        Output of :func:`numeric_summary`; ``min``, ``max``, ``mean`` and
        ``std`` are read.

    Returns
    -------
    pandas.DataFrame
        Columns ``column``, ``min``, ``max``, ``mean``, ``std``,
        ``scale_profile``, where the profile is one of ``binary``,
        ``proportion_0_1``, ``standardized``, ``log_like_small_range``,
        ``large_scale_count`` or ``unknown``.

    Notes
    -----
    Rules are checked in order and the first match wins, so the labels
    are not independent tests. They describe the summary statistics, not
    the actual transformation history — ``log_like_small_range`` means
    the values look like logs, not that a log was taken.
    """
    def tag(row):
        """Classify one column's scale from its summary statistics."""
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
    """
    Choose the highest-variance numeric columns for a plot.

    Variance is the selection criterion because a near-constant column
    produces an uninformative histogram or a blank correlation row — the
    columns that vary are the ones worth the plot's limited space.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection. Used only if ``summary`` is not supplied.
    table : str
        Table or view name.
    max_cols : int
        How many columns to return.
    summary : pandas.DataFrame, optional
        Precomputed numeric summary.

    Returns
    -------
    list of str
        Column names, highest variance first. Empty when there are no
        numeric columns.

    Notes
    -----
    Variance is compared across columns without standardisation, so on a
    table mixing scales this favours the large-scale columns regardless
    of how much relative variation they carry.
    """
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
    """
    Grid of histograms for the highest-variance numeric columns.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    out_path : Path
        Destination image file.
    cols : list of str, optional
        Explicit columns, truncated to ``max_cols``. Chosen by variance
        when omitted.
    max_cols : int, optional
        Display cap. Defaults to :data:`MAX_DIST_COLS`.
    summary : pandas.DataFrame, optional
        Precomputed numeric summary, used for the variance ranking.

    Returns
    -------
    Path or None
        The output path, or None when there are no numeric columns.

    Notes
    -----
    Drawn from a row sample, not the full table, so the shapes are
    indicative rather than exact. Each subplot is scaled independently —
    use :func:`plot_scale_comparison` to compare columns against each
    other.
    """
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
    """
    Grid of boxplots for the highest-variance numeric columns.

    The counterpart to the histogram grid: same columns, but showing
    quartiles and outliers rather than shape.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    out_path : Path
        Destination image file.
    cols : list of str, optional
        Explicit columns, truncated to ``max_cols``. Chosen by variance
        when omitted.
    max_cols : int, optional
        Display cap. Defaults to :data:`MAX_DIST_COLS`.
    summary : pandas.DataFrame, optional
        Precomputed numeric summary.

    Returns
    -------
    Path or None
        The output path, or None when there are no numeric columns.
    """
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
    """
    Boxplots of several columns on one shared symlog axis.

    Where the boxplot grid scales each column independently, this puts
    them on a common axis so their relative magnitudes are visible. A
    column sitting orders of magnitude away from its neighbours is on a
    different scale and would dominate any unstandardised model built
    over them.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    out_path : Path
        Destination image file.
    max_cols : int, optional
        Display cap. Defaults to :data:`MAX_SCALE_COLS`.
    summary : pandas.DataFrame, optional
        Precomputed numeric summary.

    Returns
    -------
    Path or None
        The output path, or None when there are no numeric columns.

    Notes
    -----
    The axis is symlog, which is logarithmic away from zero but linear
    near it, so negative values plot correctly. Columns are selected by
    lowest median, not by variance — the intent is to show the bottom of
    the range against everything else.
    """
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
    """
    Pearson correlation heatmap over the highest-variance numeric columns.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    out_path : Path
        Destination image file.
    max_cols : int, optional
        Display cap. Defaults to :data:`MAX_CORR_COLS`.
    summary : pandas.DataFrame, optional
        Precomputed numeric summary, used for the variance ranking.

    Returns
    -------
    Path or None
        The output path, or None when fewer than two numeric columns
        exist.

    Notes
    -----
    Computed on a row sample with pandas' pairwise-complete default, so
    different cells may rest on different row subsets where missingness
    differs. Numeric annotations are drawn only at 15 columns or fewer,
    beyond which they are unreadable.
    """
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
    """
    Mean and standard deviation per numeric column.

    The minimal pair of statistics needed to draw a column as a Gaussian
    curve in the cross-table density figure.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    chunk : int, optional
        Columns per query. Defaults to :data:`MAX_AGG_CHUNK`.
    cols : list of str, optional
        Restrict to these columns. All numeric columns when omitted.

    Returns
    -------
    pandas.DataFrame
        Columns ``column``, ``mean``, ``sd``. Non-finite rows and
        zero-variance columns are dropped, since neither can be drawn as
        a curve.
    """
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
    """
    Cross-table figure: one faceted panel per table, one curve per column.

    Each numeric column is drawn as the Gaussian implied by its mean and
    standard deviation. Individual curves are not the point — the shape
    of the pile is. A panel of tightly overlapping curves means the
    columns share a scale; a wide spread means they do not, which is what
    the downstream z-scoring has to handle.

    Restricted to :data:`GAUSSIAN_FACET_TABLES`, the tables whose numeric
    columns are genuine measurements.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    tables : list of str
        Candidate tables, filtered against the facet set.
    out_path : Path
        Destination image file.
    max_cols : int, optional
        Curves per panel, sampled when a table has more. Defaults to
        :data:`MAX_GAUSSIAN_COLS`.
    ncols : int, optional
        Panels per row. Default 3.
    grid : int, optional
        Points per curve. Default 400.
    seed : int, optional
        Seed for the column sample, so the figure is reproducible.
        Default 7.
    table_cols : dict, optional
        Per-table column restriction, used to limit the gene matrices to
        protein-coding genes.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure, already saved and closed.
    moments : pandas.DataFrame
        Every column's mean and sd across all panels, with a ``table``
        column — saved by the caller as the figure's underlying data.

    Notes
    -----
    Curves are coloured by the percentile rank of their mean, so colour
    encodes position along the x-axis rather than any category. The
    x-range is clipped to the 1st and 99th percentile of the +/-3 sd
    bounds, so extreme columns are cropped rather than flattening
    everything else.

    Sets a serif font family globally via ``rcParams``, which persists
    for the rest of the session.
    """
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


#: Bin edges for the missingness bands. The first bin isolates exactly
#: zero — a complete column is categorically different from a nearly
#: complete one — and the last extends past 100 so fully-missing columns
#: fall inside it.
MISS_BINS = [-0.001, 0.0, 10, 25, 50, 75, 95, 100.001]
MISS_LABELS = ["complete", "<10%", "10–25%", "25–50%", "50–75%", "75–95%", "≥95%"]

#: Green-to-red ramp aligned to MISS_LABELS, so severity reads off colour
#: alone across every panel of the cross-table figure.
MISS_COLORS = ["#2E7D5B", "#7FB069", "#C9CE6E", "#EDC26B", "#E39A5C", "#D9704E", "#B33A3A"]


def missing_pct(con, table: str, chunk: int = MAX_AGG_CHUNK, cols: Optional[list] = None) -> pd.Series:
    """
    Percentage of rows missing, per column, computed in DuckDB.

    The lighter counterpart to :func:`missing_value_report` — just the
    percentage, for the cross-table figure that needs nothing else.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Table or view name.
    chunk : int, optional
        Columns per query. Defaults to :data:`MAX_AGG_CHUNK`.
    cols : list of str, optional
        Restrict to these columns. All columns when omitted.

    Returns
    -------
    pandas.Series
        Percent missing, indexed by column name. Empty when the table has
        no rows or no columns.
    """
    cs = cols if cols is not None else cols_of(con, table)
    total = row_count(con, table)
    if total == 0 or not cs:
        return pd.Series(dtype=float)
    return chunked_agg(con, table, cs,
                        lambda c: f'100.0*sum(CASE WHEN "{c}" IS NULL THEN 1 ELSE 0 END)/{total}', chunk)


def plot_missingness(con, tables: list, out_path: Path, ncols: int = 4, table_cols: Optional[dict] = None):
    """
    Cross-table figure: how each table's columns distribute across missingness bands.

    One panel per table, counting columns in each band rather than
    listing them. This is the whole-project view of completeness — which
    tables are largely complete and which are mostly gaps — where the
    per-table bar chart names individual columns.

    Restricted to :data:`MISSINGNESS_FACET_TABLES` and laid out in that
    order, so the figure is stable across runs.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    tables : list of str
        Candidate tables, filtered against the facet set.
    out_path : Path
        Destination image file.
    ncols : int, optional
        Panels per row. Default 4.
    table_cols : dict, optional
        Per-table column restriction for the gene matrices.

    Returns
    -------
    matplotlib.figure.Figure
        The figure, already saved and closed.

    Notes
    -----
    Each panel's x-limit is set from its own maximum, so bar lengths are
    comparable within a panel but not across panels — the column counts
    in the panel titles differ by orders of magnitude. Read the shape of
    the distribution, not the bar lengths.

    Sets a serif font family globally via ``rcParams``, which persists
    for the rest of the session.
    """
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
    """
    Run the full EDA suite over one table.

    Computes the two expensive statistics once — the missing-value report
    and the numeric summary, both cached to disk — then passes them into
    every plot and test rather than letting each recompute. That sharing
    is what makes a wide table tractable.

    In order: the missing-value report and its four plots, the two
    missingness randomness tests and their written verdict, then the
    numeric summary, the scale profile, and four distribution plots.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    name : str
        Table or view name. Also the output subdirectory name.
    out_dir : Path
        Parent output directory; a subdirectory per table is created.
    protein_coding_genes : list of str, optional
        Gene IDs to restrict the gene-matrix tables to. No restriction
        when None.
    force : bool, optional
        Ignore the on-disk caches and recompute. Default False.

    Returns
    -------
    dict
        One overview row: ``dataset``, ``rows``, ``cols``,
        ``cols_analyzed``, ``overall_pct_missing``, ``mcar_p_value``,
        ``mcar_reject_at_0.05``, ``missingness_verdict``. Fields whose
        step failed or was skipped stay None.

    Notes
    -----
    Every step runs through a local wrapper that catches, logs and
    returns None, so no single failure costs the rest of the table's
    output. A run with several ``[SKIPPED]`` lines still produces a valid
    overview row — check the log rather than assuming the output is
    complete.

    ``overall_pct_missing`` is computed over the analysed columns only,
    so for a gene-matrix table it describes the protein-coding subset,
    not the whole table.

    Written to the table's directory: ``missing_value_report.csv``,
    ``numeric_summary.csv``, ``scale_profile.csv``,
    ``missingness_dependency_check.csv``, ``missingness_verdict.txt``,
    and eight PNGs.
    """
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
        """Run one EDA step, logging and swallowing any failure."""
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
    """
    Run the EDA suite over every table, then build the cross-table figures.

    Loads the protein-coding restriction once and applies it to the
    gene-matrix tables throughout, runs :func:`run_eda_on_table` per
    table, assembles the overview, and finishes with the two figures that
    span all tables at once.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection with every table registered.
    table_names : list of str
        Tables to analyse.
    out_dir : Path
        Output root, created if absent.
    force : bool, optional
        Ignore the on-disk caches and recompute. Default False.

    Returns
    -------
    pandas.DataFrame
        One row per successfully analysed table, as returned by
        :func:`run_eda_on_table`. Tables that failed entirely are absent.

    Notes
    -----
    A table that fails outright is logged and skipped, so the run
    continues — the returned overview may therefore be shorter than
    ``table_names``.

    Written to ``out_dir``: ``_overview.csv``,
    ``_overview_missing_by_dataset.png``,
    ``_summary_gaussian_per_column.png``,
    ``_summary_column_moments.csv``, and
    ``_summary_missingness_by_table.png``. The leading underscore keeps
    the cross-table outputs sorted above the per-table directories.

    Each cross-table figure is wrapped in its own try/except, so a
    failure in one still leaves the other and the whole per-table output
    intact.
    """
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