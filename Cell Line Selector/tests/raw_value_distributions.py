"""
05_raw_value_distributions.py

Build one small histogram facet per raw parquet dataset. Values are pooled
across measurement columns, while obvious identifiers and genomic coordinates
are excluded. Each facet is sampled and clipped only for display so that
outliers do not flatten the central distribution.

Run with:
    python pipelines/05_raw_value_distributions.py
"""

import sys
import textwrap
from pathlib import Path

project_root = Path(__file__).resolve().parent
for candidate in [project_root, *project_root.parents]:
    if (candidate / "src").exists():
        project_root = candidate
        break
sys.path.insert(0, str(project_root))

import duckdb  # type: ignore
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

from src.scripts.data_utils import EDA_DIR, PARQUET_RAW
from src.scripts.eda_functions import (
    fetch_raw_value_sample,
    raw_measurement_cols,
    style_axis,
    _safe_savefig,
)


def load_raw_views():
    con = duckdb.connect(database=":memory:")
    table_names = []
    for path in sorted(PARQUET_RAW.glob("*.parquet")):
        name = path.stem
        con.execute(f'CREATE VIEW "{name}" AS SELECT * FROM read_parquet(\'{path.as_posix()}\')')
        table_names.append(name)
    return con, table_names


def plot_raw_value_distribution_facets(con, table_names, out_path: Path):
    facets = []
    for table in table_names:
        cols = raw_measurement_cols(con, table)
        values = fetch_raw_value_sample(con, table, cols)
        if len(values):
            shown_cols = min(len(cols), 5 * 100)
            facets.append((table, shown_cols, values))

    if not facets:
        raise ValueError("No numeric raw measurement values were found")

    ncols = 3
    nrows = int(np.ceil(len(facets) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11.0, 3.0 * nrows), squeeze=False)
    axes = axes.ravel()

    for panel_number, (ax, (table, n_measurement_cols, values)) in enumerate(zip(axes, facets), start=1):
        low, high = np.percentile(values, [0.5, 99.5])
        shown = values if low == high else values[(values >= low) & (values <= high)]
        ax.hist(shown, bins=35, density=True, color="#356A8A", alpha=0.9,
                edgecolor="white", linewidth=0.3)
        title = textwrap.fill(table, width=22)
        ax.set_title(f"({chr(96 + panel_number)}) {title}", fontsize=9,
                     fontweight="bold", loc="left", pad=5)
        ax.text(0.01, 0.91, f"{n_measurement_cols:,} variables; n = {len(values):,}",
                transform=ax.transAxes, fontsize=7.5, color="#444444")
        ax.set_xlabel("Raw value", fontsize=8.5)
        ax.set_ylabel("Density", fontsize=8.5)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.tick_params(axis="both", labelsize=7.5)
        style_axis(ax)

    for ax in axes[len(facets):]:
        ax.axis("off")

    fig.suptitle("Distribution of raw measurement values across datasets",
                 fontsize=13, fontweight="bold", x=0.01, ha="left")
    fig.text(0.01, 0.006,
             "Histograms pool numeric measurement variables within each dataset. Values beyond the 0.5th and 99.5th percentiles are omitted from display.",
             fontsize=8, color="#444444")
    fig.subplots_adjust(left=0.07, right=0.99, top=0.94, bottom=0.045,
                        wspace=0.28, hspace=0.58)
    _safe_savefig(fig, out_path)
    return out_path


def main() -> None:
    if not PARQUET_RAW.exists():
        raise FileNotFoundError(f"Raw parquet directory not found: {PARQUET_RAW}")
    con, table_names = load_raw_views()
    out_path = EDA_DIR / "_raw_value_distribution_facets.png"
    plot_raw_value_distribution_facets(con, table_names, out_path)
    print(f"Saved raw distribution facets to {out_path}")
    con.close()


if __name__ == "__main__":
    main()
