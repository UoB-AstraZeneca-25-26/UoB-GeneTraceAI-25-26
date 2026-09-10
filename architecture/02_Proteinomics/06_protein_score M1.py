from __future__ import annotations

"""
Cell-line protein z-scores from two proteomics sources (ProCAN + CCLE).

Method 1 semantics
------------------
Both platforms are loaded, mapped to ENSG gene IDs and compared *before*
z-scoring. For each protein measured on both platforms with enough shared
cell lines, a Spearman rho decides a tier:

    rho >= 0.5   -> consistent
    rho >= 0.3   -> cautious
    otherwise    -> conflicting

Proteins tiered as conflicting have their **CCLE** rows dropped, so the
merged score for those genes rests on ProCAN alone. This asymmetry — a
standing preference for ProCAN where the platforms disagree — is the
defining behaviour of Method 1.

Surviving rows are then z-scored per platform, conditioned on lineage and
robust to outliers (median / MAD), and combined across platforms with
coverage weights ``w = sqrt(n_models)``.

Output table
------------
    gene_id | model_id | lineage | z_score | n_sources

where ``n_sources`` is 1 (one platform, or CCLE dropped as conflicting)
or 2 (ProCAN + CCLE).

Thresholds
----------
``TIER_CUTOFFS`` and ``MIN_SHARED`` are fixed domain values set by hand.

Usage
-----
    python <this_file>.py
    python <this_file>.py --tables proteomics protein_matrix_averaged_20250211 \
        --output-table protein_z
"""

import argparse
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy.stats import spearmanr  # type: ignore


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DB = ROOT / "db" / "celllineselector.duckdb"

#: Spearman rho cutoffs defining the platform-agreement tiers. A protein at
#: or above ``consistent`` is trusted on both platforms; at or above
#: ``cautious`` it is kept with reservations; below that it is conflicting
#: and its CCLE rows are dropped. Hand-set domain values.
TIER_CUTOFFS = {
    "consistent": 0.5,
    "cautious": 0.3,
}

#: Minimum shared cell lines before a ProCAN/CCLE correlation is trusted
#: enough to tier a protein. Proteins below this are left untiered rather
#: than tiered on thin evidence.
MIN_SHARED = 30

# Current DuckDB source names.
# Method 1 semantics:
#   ProCAN = protein_matrix_averaged_20250211
#   CCLE   = proteomics
PROCAN_SOURCE = "protein_matrix_averaged_20250211"
CCLE_SOURCE = "proteomics"


def canonical_ensg(value) -> str | None:
    """
    Normalise a value to a bare, version-stripped Ensembl gene ID.

    Accepts ``ENSG00000141510`` or ``ENSG00000141510.17`` in any case and
    returns the uppercase ID with any version suffix removed. Anything
    that is not a well-formed ENSG ID returns None, so callers can drop
    unmapped rows with ``dropna``.

    Parameters
    ----------
    value : any
        Candidate identifier. NaN, None and empty strings are tolerated.

    Returns
    -------
    str or None
        e.g. ``"ENSG00000141510"``, or None if the value is missing or
        does not match the ENSG pattern.
    """
    if value is None or pd.isna(value):
        return None

    text = str(value).strip().upper()

    if not text:
        return None

    match = re.fullmatch(r"ENSG\d+(?:\.\d+)?", text)
    return text.split(".", 1)[0] if match else None


def normalize_uniprot(value) -> str:
    """
    Normalise a UniProt accession for use as a lookup key.

    Trims surrounding whitespace and lowercases, matching the key
    convention used by :func:`build_uniprot_to_ensg`.

    Parameters
    ----------
    value : any
        Candidate accession. NaN and None are tolerated.

    Returns
    -------
    str
        Lowercased accession, or an empty string if the value is missing.

    Notes
    -----
    Currently unused — the lookup paths call ``.str.lower()`` inline.
    """
    if value is None or pd.isna(value):
        return ""

    return str(value).strip().lower()


def split_id_list(value) -> list[str]:
    """
    Split a delimited identifier string into a list of clean tokens.

    Roster columns hold multiple identifiers per gene in a single cell,
    with no one consistent delimiter, so semicolons, commas, pipes and
    both slashes are all treated as separators. Tokens are uppercased and
    the string forms of missing values are discarded.

    Parameters
    ----------
    value : any
        Delimited identifier string. NaN, None and empty strings are
        tolerated.

    Returns
    -------
    list of str
        Uppercased tokens, empty if the value held nothing usable.

    Examples
    --------
    >>> split_id_list("P04637; Q9H3D4|nan")
    ['P04637', 'Q9H3D4']
    """
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
    """
    Build a lowercase UniProt accession -> ENSG lookup from ``gene_roster``.

    Reads both accession columns — ``uniprot_id(s)`` and
    ``procan_uniprot_id`` — splits each multi-value cell, and maps every
    token to its canonical gene ID.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the project database.

    Returns
    -------
    dict
        Lowercased accession to canonical ENSG ID. Later columns
        overwrite earlier ones where an accession appears in both, so
        ``procan_uniprot_id`` wins on collision.

    Notes
    -----
    A missing column is skipped via ``KeyError`` rather than raising, so
    the map is built from whatever columns the roster actually has.
    """

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
    """
    Build a lowercase gene-symbol -> ENSG lookup from ``gene_roster``.

    Reads the ``gene_name(s)`` column, splits each multi-value cell, and
    maps every symbol and alias to its canonical gene ID. Used as the
    fallback identifier space for platforms whose columns are gene
    symbols rather than UniProt accessions.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the project database.

    Returns
    -------
    dict
        Lowercased symbol to canonical ENSG ID. Empty if the roster has
        no ``gene_name(s)`` column.
    """

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
    """
    Column-wise robust z-scores using median and MAD.

    Centres each column on its median and scales by ``1.4826 * MAD``, the
    MAD-to-sigma constant for a normal distribution. Median/MAD is used in
    place of mean/SD because protein abundance distributions across cell
    lines carry genuine extreme values that would otherwise inflate the
    denominator and flatten the very signal being looked for.

    Parameters
    ----------
    X : numpy.ndarray
        2-D array, rows are cell lines and columns are genes. Cast to
        float internally.

    Returns
    -------
    numpy.ndarray
        Array of the same shape, z-scored per column. An empty input is
        returned unchanged.

    Notes
    -----
    Two degenerate cases are handled so no column returns infinities: a
    zero MAD falls back to the standard deviation, and a zero scale after
    that fallback becomes 1.0, leaving the column as plain deviations
    from its median.
    """

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
    """
    Reduce a cell-line name to a comparison key.

    Lowercases and strips every non-alphanumeric character, so that
    ``MCF-7``, ``MCF 7`` and ``mcf7`` all collapse to ``mcf7``. This is
    what lets the protein matrix's ``ccle_name`` match ``sample_info``
    entries that punctuate the same line differently.

    Parameters
    ----------
    value : any
        Cell-line name. NaN and None are tolerated.

    Returns
    -------
    str
        Alphanumeric-only lowercase key, empty if the value is missing.
    """

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
    """
    Attach ``model_id`` to a protein table keyed only by ``ccle_name``.

    Builds a name lookup from ``sample_info`` across three identifier
    columns — ``ccle_id``, ``cell_line_name`` and
    ``stripped_cell_line_name`` — normalising each through
    :func:`normalize_name`. The columns are tried in that order and
    ``setdefault`` keeps the first hit, so ``ccle_id`` takes precedence
    where two columns normalise to the same key.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the project database.
    df : pandas.DataFrame
        Protein table. Returned unchanged if it already has ``model_id``.

    Returns
    -------
    pandas.DataFrame
        Copy of ``df`` with a ``model_id`` column, with rows whose name
        matched nothing dropped.

    Raises
    ------
    ValueError
        If ``df`` has neither ``model_id`` nor ``ccle_name``, or if
        ``sample_info`` yields no usable mappings — in both cases nothing
        downstream could be keyed, so failing loudly beats returning an
        empty frame.
    """

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
    """
    Melt a wide protein matrix into long form.

    Treats every column that is neither ``id_col`` nor listed in
    ``meta_cols`` as a protein measurement, stacking them into one row per
    (id, protein) pair. Missing measurements are dropped rather than
    carried as NaN, so the long frame holds observations only.

    Parameters
    ----------
    df : pandas.DataFrame
        Wide table, one row per cell line and one column per protein.
    id_col : str
        Column identifying the cell line, typically ``model_id``.
    meta_cols : set of str, optional
        Non-measurement columns to exclude from the melt. Defaults to
        none.

    Returns
    -------
    pandas.DataFrame
        Long frame with columns ``[id_col, "protein", "value"]`` and a
        reset index.

    Notes
    -----
    ``DataFrame.stack`` is called with ``future_stack=True`` where the
    installed pandas supports it, falling back to the legacy signature on
    ``TypeError`` — this keeps the file working across the pandas 2.x
    stack deprecation.
    """

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
    """
    Convert one wide protein source table into a mapped long frame.

    Ensures the table is keyed by ``model_id`` (attaching it from
    ``ccle_name`` if needed), melts it to long form, lowercases the
    protein identifiers to match the lookup convention, and maps each to
    a canonical ENSG ID. Rows whose protein has no mapping, or whose
    mapped value is not a valid ENSG, are dropped.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection, used only if ``model_id`` must be attached.
    df : pandas.DataFrame
        Wide protein table as read from DuckDB.
    source_name : str
        Platform name, recorded in the ``source`` column and used in
        error messages.
    protein_map : dict
        Lowercased UniProt accession or gene symbol to ENSG ID.

    Returns
    -------
    pandas.DataFrame
        Long frame with columns ``model_id``, ``protein``, ``value``,
        ``gene_id`` and ``source``.

    Raises
    ------
    ValueError
        If the table has neither ``model_id`` nor ``ccle_name``.

    Notes
    -----
    Unmapped proteins are dropped silently. A platform whose column
    identifiers are absent from ``protein_map`` will therefore yield an
    empty frame rather than an error — worth checking row counts per
    source when a platform contributes nothing.
    """

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
    """
    Load the ``model_id`` -> lineage mapping from ``sample_info``.

    Lineage is what the z-scores are conditioned on, so that a protein is
    scored against cell lines of the same tissue of origin rather than
    against the whole panel.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the project database.

    Returns
    -------
    pandas.Series
        Lineage indexed by ``model_id``, with missing lineages filled as
        ``"unknown"``. Empty if ``sample_info`` yields no rows.

    Notes
    -----
    Duplicates are collapsed by keeping the first row per ``model_id``,
    so a model with conflicting lineage labels resolves arbitrarily but
    deterministically.
    """

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
    """
    Compute lineage-conditioned robust z-scores for one platform.

    Pivots the long frame to a cell-line x gene matrix, attaches lineage,
    and z-scores each lineage group independently via
    :func:`robust_z_matrix`. Scoring within lineage means a value is
    judged against comparable cell lines, so lineage-typical abundance
    does not read as gene-specific signal.

    Parameters
    ----------
    df : pandas.DataFrame
        Long frame for one platform, with ``model_id``, ``gene_id`` and
        ``value``.
    lineage_map : pandas.Series
        Lineage indexed by ``model_id``, from :func:`load_lineage_map`.
    source_name : str
        Platform name, recorded in the ``source`` column.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_id``, ``model_id``, ``lineage``, ``z_score``,
        ``source``. An empty frame with those columns is returned when
        nothing survives the lineage join.

    Notes
    -----
    Cell lines with no lineage are dropped rather than pooled into an
    ``unknown`` group at this stage — ``lineage_map`` only fills
    ``unknown`` for models it already knows about. Duplicate
    (model, gene) pairs are resolved by ``aggfunc="first"``.
    """

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
    """
    Map a platform-agreement correlation to its tier label.

    Parameters
    ----------
    rho : float
        Spearman correlation between the two platforms for one protein.

    Returns
    -------
    str
        ``"consistent"``, ``"cautious"`` or ``"conflicting"``, per
        :data:`TIER_CUTOFFS`.

    Notes
    -----
    Only ``"conflicting"`` changes behaviour downstream — it triggers
    removal of the protein's CCLE rows. ``"cautious"`` is recorded but
    carries no consequence in this version.
    """

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

    Correlation is computed on the raw values, before z-scoring, so the
    tier reflects agreement between the platforms themselves rather than
    agreement after each has been normalised within lineage.

    Parameters
    ----------
    procan_long : pandas.DataFrame
        ProCAN long frame with ``model_id``, ``protein``, ``value``.
    ccle_long : pandas.DataFrame
        CCLE long frame with the same columns.

    Returns
    -------
    pandas.DataFrame
        Columns ``protein``, ``platform_rho``, ``platform_tier``,
        ``n_shared`` — one row per protein that cleared
        :data:`MIN_SHARED`. Empty (with those columns) when the platforms
        share no observations.

    Notes
    -----
    The inner merge is on ``model_id`` **and** ``protein``, so only
    proteins that both platforms label with the same identifier string
    can ever be tiered.
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
    """
    Run the full Method 1 protein scoring, end to end.

    Eight steps, in order:

    1. Build the identifier maps. Symbols first, then UniProt accessions
       merged over them, so an accession wins where a token appears in
       both spaces.
    2. Load the lineage map.
    3. Load each source table and convert it to a mapped long frame. A
       table that is missing, empty or wholly unmapped is skipped.
    4. Tier the platforms on raw values, before any z-scoring.
    5. **Drop CCLE rows for conflicting genes.** This is the critical
       Method 1 behaviour — where the platforms disagree, ProCAN is kept
       and CCLE discarded.
    6. Z-score each surviving platform independently, conditioned on
       lineage.
    7. Weight each platform by its coverage, ``w = sqrt(n_models)``, so a
       platform measuring a gene across more cell lines counts for more.
    8. Combine as ``sum(w*z) / sqrt(sum(w^2))``, which preserves the
       z-scale of the inputs. Because step 5 already removed conflicting
       CCLE rows, ``n_sources`` falls out automatically: 2 means both
       platforms agreed enough to contribute, 1 means one platform only.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the project database.
    table_names : list of str, optional
        Source tables to score. Defaults to ``[CCLE_SOURCE,
        PROCAN_SOURCE]``.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_id``, ``model_id``, ``lineage``, ``z_score``,
        ``n_sources``, sorted by gene then model, with IDs uppercased.

    Raises
    ------
    ValueError
        If no identifier mapping could be built from ``gene_roster``, if
        no source table yielded rows, or if nothing could be scored.
    """

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
    """
    Write the scored results to DuckDB, replacing any existing table.

    Registers the frame as a temporary view and materialises it with an
    explicit column list, so the output schema is fixed regardless of any
    extra columns the caller's frame happens to carry.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open, writable connection.
    df : pandas.DataFrame
        Result frame from :func:`compute_protein_z`.
    output_table : str, optional
        Destination table name. Defaults to ``"protein_z"``.

    Returns
    -------
    None
        The table is created or replaced as a side effect.

    Notes
    -----
    ``CREATE OR REPLACE`` means a rerun overwrites the previous results
    outright — there is no append path and no versioning.
    """

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
    """
    Define and parse the command-line interface.

    Returns
    -------
    argparse.Namespace
        Parsed arguments:

        * ``db`` (str) — path to the DuckDB file;
        * ``tables`` (list of str) — protein source tables, defaulting to
          both project protein tables;
        * ``output_table`` (str) — destination table for the scores.

    Notes
    -----
    The parser description is taken from the module docstring, so the
    Method 1 semantics above appear in ``--help``.
    """

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
    """
    Run the Method 1 protein scoring from the command line.

    Opens the database, scores the requested source tables, writes the
    results, then prints a summary read back from the written table —
    row, gene, model and lineage counts, plus how many rows rest on one
    versus two platforms — followed by the source-count distribution.

    Returns
    -------
    None
        The output table and console summary are the side effects.

    Raises
    ------
    ValueError
        Propagated from :func:`compute_protein_z` when no mapping, no
        source rows or no scoreable rows are available.

    Notes
    -----
    The connection is closed in a ``finally`` block, so it is released
    even when scoring raises. The summary is read back from the database
    rather than computed from the in-memory frame, which verifies the
    write actually landed.
    """

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