from __future__ import annotations

import argparse
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy.stats import spearmanr  # type: ignore


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DB = ROOT / "db" / "celllineselector.duckdb"

# Method 1 platform-tier rules
TIER_CUTOFFS = {
    "consistent": 0.5,
    "cautious": 0.3,
}
MIN_SHARED = 30

# Current DuckDB source names.
# Method 1 semantics:
#   ProCAN = protein_matrix_averaged_20250211
#   CCLE   = proteomics
PROCAN_SOURCE = "protein_matrix_averaged_20250211"
CCLE_SOURCE = "proteomics"


def canonical_ensg(value) -> str | None:
    if value is None or pd.isna(value):
        return None

    text = str(value).strip().upper()

    if not text:
        return None

    match = re.fullmatch(r"ENSG\d+(?:\.\d+)?", text)
    return text.split(".", 1)[0] if match else None


def normalize_uniprot(value) -> str:
    if value is None or pd.isna(value):
        return ""

    return str(value).strip().lower()


def split_id_list(value) -> list[str]:
    if value is None or pd.isna(value):
        return []

    text = str(value).strip()

    if not text:
        return []

    parts = re.split(r"[;,|/\\]", text)

    out = []

    for p in parts:
        token = str(p).strip().upper()

        if token and token.lower() not in {"nan", "none", "null"}:
            out.append(token)

    return out


def build_uniprot_to_ensg(
    con: duckdb.DuckDBPyConnection,
) -> dict[str, str]:

    mapping: dict[str, str] = {}

    roster = con.execute(
        'SELECT gene_id, "uniprot_id(s)", procan_uniprot_id '
        'FROM gene_roster'
    ).fetchdf()

    for col in ['"uniprot_id(s)"', "procan_uniprot_id"]:

        try:
            subset = roster[
                ["gene_id", col]
            ].dropna(
                subset=["gene_id", col]
            )
        except KeyError:
            continue

        for _, row in subset.iterrows():

            gene = canonical_ensg(row["gene_id"])

            if gene is None:
                continue

            for token in split_id_list(row[col]):

                if token:
                    mapping[token.lower()] = gene

    return mapping


def build_symbol_to_ensg(
    con: duckdb.DuckDBPyConnection,
) -> dict[str, str]:

    mapping: dict[str, str] = {}

    roster = con.execute(
        'SELECT gene_id, "gene_name(s)" FROM gene_roster'
    ).fetchdf()

    if "gene_name(s)" not in roster.columns:
        return mapping

    subset = roster[
        ["gene_id", "gene_name(s)"]
    ].dropna(
        subset=["gene_id", "gene_name(s)"]
    )

    for _, row in subset.iterrows():

        gene = canonical_ensg(row["gene_id"])

        if gene is None:
            continue

        for token in split_id_list(row["gene_name(s)"]):

            if token:
                mapping[token.lower()] = gene

    return mapping


def robust_z_matrix(X: np.ndarray) -> np.ndarray:

    X = np.asarray(X, dtype=float)

    if X.size == 0:
        return X

    med = np.median(
        X,
        axis=0,
        keepdims=True,
    )

    mad = np.median(
        np.abs(X - med),
        axis=0,
        keepdims=True,
    )

    scale = np.where(
        np.isclose(mad, 0.0),
        np.std(X, axis=0, keepdims=True),
        1.4826 * mad,
    )

    scale = np.where(
        np.isclose(scale, 0.0),
        1.0,
        scale,
    )

    return (X - med) / scale


def normalize_name(value) -> str:

    if value is None or pd.isna(value):
        return ""

    return re.sub(
        r"[^a-z0-9]+",
        "",
        str(value).strip().lower(),
    )


def attach_model_ids(
    con: duckdb.DuckDBPyConnection,
    df: pd.DataFrame,
) -> pd.DataFrame:

    if "model_id" in df.columns:
        return df

    if "ccle_name" not in df.columns:
        raise ValueError(
            "The protein matrix table is missing both "
            "model_id and ccle_name."
        )

    meta = con.execute(
        """
        SELECT
            model_id,
            cell_line_name,
            stripped_cell_line_name,
            ccle_id
        FROM sample_info
        WHERE model_id IS NOT NULL
        """
    ).fetchdf()

    if meta.empty:
        raise ValueError(
            "sample_info has no usable cell-line mappings "
            "for model_id lookup."
        )

    lookup: dict[str, str] = {}

    for col in [
        "ccle_id",
        "cell_line_name",
        "stripped_cell_line_name",
    ]:

        if col not in meta.columns:
            continue

        for _, row in meta[
            ["model_id", col]
        ].dropna().iterrows():

            key = normalize_name(row[col])

            if key:
                lookup.setdefault(
                    key,
                    str(row["model_id"]),
                )

    out = df.copy()

    out["model_id"] = out["ccle_name"].map(
        lambda v: lookup.get(
            normalize_name(v),
            None,
        )
    )

    return out.dropna(
        subset=["model_id"]
    ).copy()


def wide_to_long(
    df: pd.DataFrame,
    id_col: str,
    meta_cols: set[str] | None = None,
) -> pd.DataFrame:

    meta_cols = meta_cols or set()

    protein_cols = [
        c
        for c in df.columns
        if c != id_col and c not in meta_cols
    ]

    sub = df[
        [id_col] + protein_cols
    ].set_index(id_col)

    try:
        long = sub.stack(
            future_stack=True
        ).reset_index()
    except TypeError:
        long = sub.stack().reset_index()

    long.columns = [
        id_col,
        "protein",
        "value",
    ]

    return long.dropna(
        subset=["value"]
    ).reset_index(drop=True)


def table_to_long(
    con: duckdb.DuckDBPyConnection,
    df: pd.DataFrame,
    source_name: str,
    protein_map: dict[str, str],
) -> pd.DataFrame:

    work = df.copy()

    if (
        "model_id" not in work.columns
        and "ccle_name" in work.columns
    ):
        work = attach_model_ids(
            con,
            work,
        )

    if "model_id" not in work.columns:
        raise ValueError(
            f"Table '{source_name}' has neither "
            "model_id nor ccle_name."
        )

    long = wide_to_long(
        work,
        "model_id",
    )

    long["protein"] = (
        long["protein"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    long["gene_id"] = long[
        "protein"
    ].map(protein_map)

    long = long.dropna(
        subset=["gene_id"]
    ).copy()

    long["gene_id"] = long[
        "gene_id"
    ].map(canonical_ensg)

    long = long.dropna(
        subset=["gene_id"]
    ).copy()

    long["source"] = source_name

    return long


def load_lineage_map(
    con: duckdb.DuckDBPyConnection,
) -> pd.Series:

    df = con.execute(
        """
        SELECT DISTINCT
            model_id,
            lineage
        FROM sample_info
        WHERE model_id IS NOT NULL
        """
    ).fetchdf()

    if df.empty:
        return pd.Series(dtype=object)

    return (
        df.drop_duplicates(
            subset=["model_id"]
        )
        .set_index("model_id")["lineage"]
        .fillna("unknown")
    )


def zscore_platform(
    df: pd.DataFrame,
    lineage_map: pd.Series,
    source_name: str,
) -> pd.DataFrame:

    wide = df.pivot_table(
        index="model_id",
        columns="gene_id",
        values="value",
        aggfunc="first",
    )

    wide["lineage"] = lineage_map.reindex(
        wide.index
    )

    wide = wide.dropna(
        subset=["lineage"]
    )

    if wide.empty:
        return pd.DataFrame(
            columns=[
                "gene_id",
                "model_id",
                "lineage",
                "z_score",
                "source",
            ]
        )

    rows = []

    for lineage, grp in wide.groupby(
        "lineage",
        sort=True,
    ):

        gene_cols = [
            c
            for c in grp.columns
            if c != "lineage"
        ]

        if not gene_cols:
            continue

        X = grp[
            gene_cols
        ].values.astype(float)

        Z = robust_z_matrix(X)

        zdf = (
            pd.DataFrame(
                Z,
                index=grp.index,
                columns=gene_cols,
            )
            .stack()
            .reset_index()
        )

        zdf.columns = [
            "model_id",
            "gene_id",
            "z_score",
        ]

        zdf["lineage"] = lineage
        zdf["source"] = source_name

        rows.append(zdf)

    if not rows:
        return pd.DataFrame(
            columns=[
                "gene_id",
                "model_id",
                "lineage",
                "z_score",
                "source",
            ]
        )

    out = pd.concat(
        rows,
        ignore_index=True,
    )

    return out[
        [
            "gene_id",
            "model_id",
            "lineage",
            "z_score",
            "source",
        ]
    ]


# ---------------------------------------------------------------------
# Method 1: platform agreement
# ---------------------------------------------------------------------

def assign_tier(rho: float) -> str:

    if rho >= TIER_CUTOFFS["consistent"]:
        return "consistent"

    if rho >= TIER_CUTOFFS["cautious"]:
        return "cautious"

    return "conflicting"


def calculate_platform_tiers(
    procan_long: pd.DataFrame,
    ccle_long: pd.DataFrame,
) -> pd.DataFrame:
    """
    Method 1 platform comparison.

    For each UniProt:
      1. Find ProCAN/CCLE observations on shared model_id.
      2. Require >= MIN_SHARED shared models.
      3. Calculate Spearman rho.
      4. Assign consistent/cautious/conflicting.
    """

    procan = procan_long[
        ["model_id", "protein", "value"]
    ].rename(
        columns={
            "value": "value_procan",
        }
    )

    ccle = ccle_long[
        ["model_id", "protein", "value"]
    ].rename(
        columns={
            "value": "value_ccle",
        }
    )

    # Only observations available in BOTH platforms
    merged = procan.merge(
        ccle,
        on=[
            "model_id",
            "protein",
        ],
        how="inner",
    )

    if merged.empty:
        return pd.DataFrame(
            columns=[
                "protein",
                "platform_rho",
                "platform_tier",
                "n_shared",
            ]
        )

    # Number of shared cell lines per protein
    shared_counts = (
        merged.groupby("protein")
        .size()
        .rename("n_shared")
    )

    eligible = shared_counts[
        shared_counts >= MIN_SHARED
    ].index

    merged = merged[
        merged["protein"].isin(eligible)
    ].copy()

    results = []

    for protein, grp in merged.groupby(
        "protein"
    ):

        if len(grp) < MIN_SHARED:
            continue

        rho, _ = spearmanr(
            grp["value_procan"],
            grp["value_ccle"],
        )

        if not np.isfinite(rho):
            continue

        results.append(
            {
                "protein": protein,
                "platform_rho": float(rho),
                "platform_tier": assign_tier(
                    float(rho)
                ),
                "n_shared": int(len(grp)),
            }
        )

    return pd.DataFrame(results)


def compute_protein_z(
    con: duckdb.DuckDBPyConnection,
    table_names: list[str] | None = None,
) -> pd.DataFrame:

    if table_names is None:
        table_names = [
            CCLE_SOURCE,
            PROCAN_SOURCE,
        ]

    # ---------------------------------------------------------------
    # 1. Build identifier mappings
    # ---------------------------------------------------------------

    uniprot_map = build_uniprot_to_ensg(con)
    symbol_map = build_symbol_to_ensg(con)

    if not (
        uniprot_map
        or symbol_map
    ):
        raise ValueError(
            "No UniProt or gene-symbol -> ENSG mapping "
            "could be built from gene_roster."
        )

    protein_map = {
        **symbol_map,
        **uniprot_map,
    }

    # ---------------------------------------------------------------
    # 2. Load lineage mapping
    # ---------------------------------------------------------------

    lineage_map = load_lineage_map(con)

    # ---------------------------------------------------------------
    # 3. Load and convert both sources
    # ---------------------------------------------------------------

    source_longs = {}

    for table_name in table_names:

        try:
            table = con.execute(
                f'SELECT * FROM "{table_name}"'
            ).fetchdf()

        except Exception:
            continue

        if table.empty:
            continue

        long = table_to_long(
            con,
            table,
            table_name,
            protein_map,
        )

        if not long.empty:
            source_longs[
                table_name
            ] = long

    if not source_longs:
        raise ValueError(
            "No protein rows could be loaded "
            "from the available source tables."
        )

    # ---------------------------------------------------------------
    # 4. PLATFORM TIER — Method 1
    #
    # Compare raw ProCAN vs CCLE values BEFORE z-scoring.
    # ---------------------------------------------------------------

    if (
        PROCAN_SOURCE in source_longs
        and CCLE_SOURCE in source_longs
    ):

        tiers = calculate_platform_tiers(
            source_longs[PROCAN_SOURCE],
            source_longs[CCLE_SOURCE],
        )

    else:

        tiers = pd.DataFrame(
            columns=[
                "protein",
                "platform_rho",
                "platform_tier",
                "n_shared",
            ]
        )

    # Proteins classified as conflicting
    conflicting_proteins = set(
        tiers.loc[
            tiers["platform_tier"]
            == "conflicting",
            "protein",
        ].astype(str)
        .str.lower()
    )

    # Convert conflicting UniProt/symbol IDs to ENSG IDs.
    conflicting_gene_ids = {
        protein_map[p]
        for p in conflicting_proteins
        if p in protein_map
    }

    # ---------------------------------------------------------------
    # 5. REMOVE CCLE for conflicting proteins
    #
    # This is the critical Method 1 behaviour.
    # ---------------------------------------------------------------

    if (
        CCLE_SOURCE in source_longs
        and conflicting_gene_ids
    ):

        source_longs[
            CCLE_SOURCE
        ] = source_longs[
            CCLE_SOURCE
        ][
            ~source_longs[
                CCLE_SOURCE
            ]["gene_id"].isin(
                conflicting_gene_ids
            )
        ].copy()

    # ---------------------------------------------------------------
    # 6. Calculate lineage-conditioned z-scores independently
    # ---------------------------------------------------------------

    blocks = []

    for source_name, long in source_longs.items():

        if long.empty:
            continue

        scored = zscore_platform(
            long,
            lineage_map,
            source_name,
        )

        if not scored.empty:
            blocks.append(scored)

    if not blocks:
        raise ValueError(
            "No protein rows could be scored "
            "from the available source tables."
        )

    all_z = pd.concat(
        blocks,
        ignore_index=True,
    )

    # ---------------------------------------------------------------
    # 7. Source coverage weighting
    # ---------------------------------------------------------------

    source_counts = (
        all_z.groupby(
            [
                "gene_id",
                "source",
            ],
            as_index=False,
        )["model_id"]
        .nunique()
        .rename(
            columns={
                "model_id": "src_n",
            }
        )
    )

    all_z = all_z.merge(
        source_counts,
        on=[
            "gene_id",
            "source",
        ],
        how="left",
    )

    all_z["w"] = np.sqrt(
        all_z["src_n"].astype(float)
    )

    all_z["wz"] = (
        all_z["w"]
        * all_z["z_score"]
    )

    all_z["w2"] = (
        all_z["w"] ** 2
    )

    # ---------------------------------------------------------------
    # 8. Combine source scores
    #
    # At this point conflicting CCLE proteins have already been
    # removed, so this automatically gives:
    #
    #   n_sources = 1 -> one source
    #   n_sources = 2 -> ProCAN + CCLE
    # ---------------------------------------------------------------

    agg = (
        all_z.groupby(
            [
                "gene_id",
                "model_id",
                "lineage",
            ],
            as_index=False,
        )
        .agg(
            wz_sum=("wz", "sum"),
            w2_sum=("w2", "sum"),
            n_sources=(
                "source",
                "nunique",
            ),
        )
    )

    denom = np.sqrt(
        agg[
            "w2_sum"
        ].to_numpy(
            dtype=float
        )
    )

    z_raw = np.where(
        denom > 0,
        agg[
            "wz_sum"
        ].to_numpy(
            dtype=float
        ) / denom,
        0.0,
    )

    result = agg[
        [
            "gene_id",
            "model_id",
            "lineage",
        ]
    ].copy()

    result["z_score"] = (
        z_raw.astype(float)
    )

    result["n_sources"] = (
        agg["n_sources"]
        .astype(int)
    )

    result = result.dropna(
        subset=[
            "gene_id",
            "model_id",
            "lineage",
            "z_score",
        ]
    ).copy()

    result["gene_id"] = (
        result["gene_id"]
        .astype(str)
        .str.upper()
    )

    result["model_id"] = (
        result["model_id"]
        .astype(str)
        .str.upper()
    )

    result["lineage"] = (
        result["lineage"]
        .astype(str)
    )

    result = (
        result
        .sort_values(
            [
                "gene_id",
                "model_id",
            ]
        )
        .reset_index(drop=True)
    )

    return result[
        [
            "gene_id",
            "model_id",
            "lineage",
            "z_score",
            "n_sources",
        ]
    ]


def write_protein_z(
    con: duckdb.DuckDBPyConnection,
    df: pd.DataFrame,
    output_table: str = "protein_z",
) -> None:

    tmp = df[
        [
            "gene_id",
            "model_id",
            "lineage",
            "z_score",
            "n_sources",
        ]
    ].copy()

    con.register(
        "protein_z_tmp",
        tmp,
    )

    con.execute(
        f'CREATE OR REPLACE TABLE "{output_table}" AS '
        "SELECT "
        "gene_id, "
        "model_id, "
        "lineage, "
        "z_score, "
        "n_sources "
        "FROM protein_z_tmp"
    )


def parse_args() -> argparse.Namespace:

    ap = argparse.ArgumentParser(
        description=__doc__
    )

    ap.add_argument(
        "--db",
        type=str,
        default=str(DEFAULT_DB),
        help="Path to the DuckDB file",
    )

    ap.add_argument(
        "--tables",
        nargs="*",
        default=[
            CCLE_SOURCE,
            PROCAN_SOURCE,
        ],
        help=(
            "Protein source tables. "
            "Default: both project protein tables."
        ),
    )

    ap.add_argument(
        "--output-table",
        type=str,
        default="protein_z",
        help=(
            "DuckDB table to write "
            "the final score results to."
        ),
    )

    return ap.parse_args()


def main() -> None:

    args = parse_args()

    con = duckdb.connect(
        args.db,
        read_only=False,
    )

    try:

        table_names = list(
            args.tables
        )

        result = compute_protein_z(
            con,
            table_names=table_names,
        )

        write_protein_z(
            con,
            result,
            output_table=args.output_table,
        )

        # -----------------------------------------------------------
        # Final summary
        # -----------------------------------------------------------

        summary = con.execute(
            f"""
            SELECT
                COUNT(*) AS n_rows,
                COUNT(DISTINCT gene_id) AS n_genes,
                COUNT(DISTINCT model_id) AS n_models,
                COUNT(DISTINCT lineage) AS n_lineages,
                COUNT(*) FILTER (
                    WHERE n_sources = 1
                ) AS n_sources_1,
                COUNT(*) FILTER (
                    WHERE n_sources = 2
                ) AS n_sources_2
            FROM "{args.output_table}"
            """
        ).fetchdf()

        print(
            summary.to_string(
                index=False
            )
        )

        print(
            "\nSource counts:"
        )

        print(
            result[
                "n_sources"
            ]
            .value_counts()
            .sort_index()
            .to_string()
        )

    finally:
        con.close()


if __name__ == "__main__":
    main()