"""
cleaning_functions.py

Per-dataset cleaning logic for the pipeline.

Every dataset's cleaning function calls `base_clean()` first (shared
lowercase/whitespace/Ensembl-version normalization), then applies
whatever dataset-specific column renaming, splitting, or reshaping
that dataset needs.

`CLEANERS` maps each raw dataset key (matching the parquet filename
produced by 00_data_loading.py) to its cleaning function. `clean_all()`
is the dispatcher `01_data_cleaning.py` used to call directly; the
pipeline script now iterates `CLEANERS` itself for per-dataset error
handling, but `clean_all()` is kept here for standalone/ad-hoc use.

Layout of this file:
  1. Regex patterns
  2. Generic helper functions (used by multiple cleaners)
  3. base_clean() — shared first pass for every dataset
  4. Per-dataset cleaning functions
  5. CLEANERS dispatch dict
  6. clean_all()

Standing conventions
--------------------
* Headers are lowercased, values are not. ``base_clean`` normalises
  column names but leaves cell values in their original capitalisation,
  so case-insensitive matching happens at comparison time and never by
  mutating the data. Every function that compares values does so
  explicitly.
* Rename, don't reorder or drop. Cleaners bring each source's column
  names onto the project's vocabulary — ``model_id``, ``gene_id``,
  ``cvcl_id``, ``cell_line_name`` — which is what lets the harmonisation
  stage join them.
* Species filtering is by accession, never by name. See the
  ``NON_HUMAN_CVCLS`` comment for the measured reason.
* Flag don't drop, with one exception. Ambiguity is preserved rather
  than resolved, except where a schema assumption fails — a missing
  species column drops every row loudly rather than letting non-human
  data through silently.

Module state
------------
``NON_HUMAN_CVCLS`` is populated as a side effect of
:func:`clean_cellosaurus` and read by :func:`filter_non_human_rows`.
This makes call order load-bearing: cellosaurus must be cleaned first or
the blocklist is empty and nothing is filtered. ``01_data_cleaning.py``
enforces that ordering.
"""

import re

import pandas as pd


# ============================================================
# Regex patterns
# ============================================================

# Matches Ensembl IDs with a version suffix, e.g. ENSG00000000003.15
# Case-insensitive since base_clean lowercases values before this runs.
ENSEMBL_VERSION_RE = re.compile(r"^(ENS[A-Za-z]*\d+)\.\d+$", re.IGNORECASE)

# Matches an Ensembl ID inside parentheses, e.g. the '(ensg00000000005)'
# in 'tnmd_(ensg00000000005)'.
ENSEMBL_ID_IN_PARENS_RE = re.compile(r"\(([^)]+)\)")

# Matches 'gene name (ensg_id.version)' cell values, e.g.
# 'dlg1 (ensg00000075711.21)', used to split fusion gene1/gene2 columns.
GENE_NAME_ID_RE = re.compile(r"^(.*?)\s*\(([^)]+)\)\s*$")

# UniProt accession pattern (official format), used to tell which half
# of a proteomics header ('a0av96_(rbm47)') is the UniProt ID and which
# is the gene symbol, regardless of which order the source file uses.
UNIPROT_ACCESSION_RE = re.compile(
    r"^[A-Z][0-9][A-Z0-9]{3}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?$", re.IGNORECASE
)

# Matches 'x_(y)' after base_clean's space-to-underscore pass, e.g.
# 'a0av96_(rbm47)'.
PROTEOMICS_COL_RE = re.compile(r"^(.+?)_\(([^)]+)\)$")


# ============================================================
# Generic helper functions
# ============================================================

def strip_ensembl_version(value):
    """Strip the trailing .N version suffix from an Ensembl ID string.

    Applied to every string value by :func:`base_clean`, so IDs from
    sources that carry versions match IDs from sources that don't.

    Parameters
    ----------
    value : any
        Candidate ID. Non-strings pass through untouched, so this is safe
        to apply across a mixed column.

    Returns
    -------
    any
        The unversioned ID, or the input unchanged when it does not match
        the versioned-Ensembl pattern.

    Notes
    -----
    Anchored, so only a value that is *entirely* a versioned Ensembl ID
    is stripped — an ID embedded in a longer string is left alone, which
    is what the parenthesised-ID helpers below handle instead.
    """
    if isinstance(value, str):
        match = ENSEMBL_VERSION_RE.match(value.strip())
        if match:
            return match.group(1)
    return value


def extract_ensembl_id(column_name: str) -> str:
    """
    Given a column like 'tnmd_(ensg00000000005)' or
    'tnmd_(ensg00000000005.6)', return just the Ensembl ID with any
    version suffix stripped. Columns without a parenthesised ID
    (e.g. a sample/index column) are returned unchanged.

    Parameters
    ----------
    column_name : str
        Header to reduce.

    Returns
    -------
    str
        The bare Ensembl ID, or the original header if it carries no
        parenthesised part.

    Notes
    -----
    Returning the input unchanged on no match is what lets a caller map
    this over every column without special-casing the identifier columns.
    """
    match = ENSEMBL_ID_IN_PARENS_RE.search(column_name)
    if not match:
        return column_name
    return strip_ensembl_version(match.group(1))


def split_gene_name_id(value):
    """
    Split a 'genename (ensgxxxxxxxxxxx.version)' string into
    (gene_name, ensembl_id) with the version suffix stripped from the
    ID. Returns (value, None) for anything that doesn't match the
    pattern (missing values, malformed cells, etc).

    Parameters
    ----------
    value : any
        Cell value to split.

    Returns
    -------
    tuple
        ``(gene_name, ensembl_id)``. A non-matching string returns
        ``(value, None)`` — keeping the name rather than discarding it —
        and a non-string returns ``(None, None)``.

    Notes
    -----
    The two failure modes differ on purpose: a malformed string still has
    a usable name, whereas a NaN has nothing.
    """
    if not isinstance(value, str):
        return None, None
    match = GENE_NAME_ID_RE.match(value)
    if not match:
        return value, None
    name, ens_id = match.groups()
    return name.strip(), strip_ensembl_version(ens_id.strip())


def is_blank_header(col: str) -> bool:
    """
    True for a column header carrying no real information — blank, or
    pandas' 'Unnamed: 0' placeholder for an unlabeled index column
    (after base_clean's lowercase/underscore pass this reads
    'unnamed:_0').

    Parameters
    ----------
    col : str
        Header to test.

    Returns
    -------
    bool
        True for an empty string, ``"index"``, or an ``unnamed``
        placeholder in any of its punctuation variants.

    Notes
    -----
    These headers mark the identifier column that a source file left
    unlabelled — the cleaners rename them to the appropriate ID rather
    than dropping them.
    """
    stripped = str(col).strip()
    return stripped == "" or stripped == "index" or bool(re.match(r"^unnamed:?_?\d*$", stripped))


def split_gene_uniprot(col: str):
    """
    Split a proteomics header into (gene_name, uniprot_id), regardless
    of which side of the parens the raw file puts the UniProt
    accession on. Returns (None, None) for headers that don't match
    the 'x_(y)' pattern (e.g. the model-id column).

    Order is decided by testing each half against the official UniProt
    accession format rather than by trusting a documented convention,
    since the source files are not consistent about it.

    Parameters
    ----------
    col : str
        Header of the form ``x_(y)``.

    Returns
    -------
    tuple
        ``(gene_name, uniprot_id)``, or ``(None, None)`` for a header
        that does not match the pattern.

    Notes
    -----
    When neither half looks like an accession, the documented
    ``genename_(uniprotid)`` order is assumed. The result is then a guess:
    the "accession" may be anything, and downstream lookups keyed on it
    simply won't match.
    """
    match = PROTEOMICS_COL_RE.match(col)
    if not match:
        return None, None
    first, second = match.group(1), match.group(2)
    if UNIPROT_ACCESSION_RE.match(first):
        return second, first  # first is the accession -> (gene_name=second, uniprot_id=first)
    if UNIPROT_ACCESSION_RE.match(second):
        return first, second  # second is the accession -> (gene_name=first, uniprot_id=second)
    # Neither looks like a UniProt accession — fall back to the
    # order described: genename_(uniprotid)
    return first, second


def rename_by_prefix(df: pd.DataFrame, prefix_map: dict) -> pd.DataFrame:
    """
    Rename each column whose cleaned name STARTS WITH a given prefix.

    Used for headers that carry extra punctuation after base_clean
    (e.g. 'Identifier (cell line Name)' -> 'identifier_(cell_line_name)')
    — an exact-string match on the raw header would silently miss
    these, since base_clean's lowercasing/whitespace-to-underscore
    step shifts the exact punctuation. Only the first matching column
    per prefix is renamed; no match = left alone.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to rename. Not modified — a renamed copy is returned.
    prefix_map : dict
        Prefix to the new column name.

    Returns
    -------
    pandas.DataFrame
        Copy with matching columns renamed.

    Notes
    -----
    First match wins, in column order, so a prefix matching several
    columns renames only one of them — and which one depends on the
    source file's column order. Keep prefixes specific enough to be
    unambiguous.
    """
    rename_map = {}
    for prefix, new_name in prefix_map.items():
        col = next((c for c in df.columns if c.startswith(prefix)), None)
        if col is not None:
            rename_map[col] = new_name
    return df.rename(columns=rename_map)


def transpose_with_id_col(df: pd.DataFrame, new_id_name: str) -> pd.DataFrame:
    """
    Transpose a wide table whose first column is a row identifier
    (gene/miRNA ID) and whose remaining columns are per-sample values.

    After transposing, samples become rows and the former column
    headers (sample IDs) become a single identifier column named
    `new_id_name`.

    Value columns are coerced to numeric (pd.to_numeric, errors="coerce").
    Raw data sometimes has non-numeric placeholders (e.g. 'ND', 'NA',
    blank strings) sitting in what should be a purely numeric column;
    left as-is, that column stays 'object' dtype (a silent mix of
    float and str) and fails when writing to parquet with an
    ArrowTypeError. Coercing turns any such placeholder into NaN,
    which is also the more honest representation of "missing" anyway.

    Parameters
    ----------
    df : pandas.DataFrame
        Features-by-samples table, identifier in the first column.
    new_id_name : str
        Name for the identifier column in the transposed result, e.g.
        ``"gsm_id"`` or ``"ccle_id"``.

    Returns
    -------
    pandas.DataFrame
        Samples-by-features table with ``new_id_name`` as its first
        column and every other column numeric.

    Notes
    -----
    Transposing to samples-as-rows is what puts these tables in the same
    orientation as the rest of the project, so they can be keyed by cell
    line like everything else.

    The coercion is silent: a genuinely mis-parsed column becomes all-NaN
    rather than raising. The EDA stage's missing-value report is where
    that would show up.
    """
    id_col = df.columns[0]
    df = df.set_index(id_col).transpose()
    df.columns.name = None
    df.index.name = new_id_name
    df = df.reset_index()

    value_cols = [c for c in df.columns if c != new_id_name]
    df[value_cols] = df[value_cols].apply(pd.to_numeric, errors="coerce")

    return df


def extract_proteomics_gene_map(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build the gene_name <-> uniprot_id lookup implied by the raw
    proteomics headers. Must be called on the raw (pre-clean_proteomics)
    DataFrame — clean_proteomics's output only keeps uniprot_id, so the
    gene_name side has to be captured first.

    Gene symbols are UPPERCASED. The raw headers carry them lowercase
    ('a0av96_(rbm47)'), which is a formatting artefact of the source
    file rather than the symbol itself -- HGNC's canonical form is
    RBM47. mutations/fusions/hpa_rna all use canonical casing, so
    uppercasing here is what makes the gene_name bridge in
    build_gene_roster actually match.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw proteomics table. Only its headers are read.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_name`` (uppercase) and ``uniprot_id``, one row per
        header that yielded both.

    Notes
    -----
    This is the only place the gene-symbol side of the proteomics headers
    is preserved — the proteomics table itself is keyed on UniProt ID
    alone after cleaning. Headers that do not parse are silently skipped,
    so compare the row count against the table's column count if the map
    seems short.
    """
    df = base_clean(df)
    rows = []
    for col in df.columns:
        gene_name, uniprot_id = split_gene_uniprot(col)
        if gene_name and uniprot_id:
            rows.append({
                "gene_name": gene_name.upper(),
                "uniprot_id": uniprot_id,
            })
    return pd.DataFrame(rows)

def extract_procan_gene_map(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract the raw ProCan gene_name <-> UniProt bridge before the metadata
    rows are dropped. Values are NOT lowercased; the metadata-row exclusion
    below compares case-insensitively without mutating the data.

    ProCan carries its identifiers in the first two rows rather than in
    headers: row 0 holds UniProt accessions, row 1 holds gene symbols. So
    the bridge has to be read before :func:`clean_procan_raw_tsv` drops
    those rows.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw ProCan frame, read with ``header=None`` so the identifier
        rows are data.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_name`` (uppercased, to match the HGNC-canonical
        casing used elsewhere) and ``procan_uniprot_ids``. Empty with
        those columns when the frame is unusable.

    Notes
    -----
    Column labels from the metadata block (``symbol``, ``model_name``,
    ``model_id``, ``uniprot_id``) are excluded on both sides, so the
    row-and-column headers do not enter the map as if they were data.
    """
    if df is None or df.empty:
        return pd.DataFrame(columns=["gene_name", "procan_uniprot_ids"])

    raw = df.copy()

    if raw.shape[0] < 2:
        return pd.DataFrame(columns=["gene_name", "procan_uniprot_ids"])

    symbol_row = raw.iloc[1].fillna("").astype(str).str.strip()
    uniprot_row = raw.iloc[0].fillna("").astype(str).str.strip()

    pairs = pd.DataFrame({
        "gene_name": symbol_row.str.upper(),
        "procan_uniprot_ids": uniprot_row,
    })

    META = {"symbol", "model_name", "model_id", "uniprot_id", "nan"}
    pairs = pairs[(pairs["gene_name"] != "") & (pairs["procan_uniprot_ids"] != "")]
    pairs = pairs[~pairs["gene_name"].str.lower().isin(META)]
    pairs = pairs[~pairs["procan_uniprot_ids"].str.lower().isin(META)]
    pairs = pairs.drop_duplicates().reset_index(drop=True)
    return pairs[["gene_name", "procan_uniprot_ids"]]


def add_procan_uniprot_ids_to_proteomics_map(
    proteomics_gene_map: pd.DataFrame,
    raw_procan_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge the raw ProCan gene_name <-> UniProt pairs into the existing
    proteomics_gene_map so the final table has:
      gene_name, uniprot_id, procan_uniprot_ids

    The raw ProCan map is built before the initial metadata rows are stripped,
    which is required for the correct aliasing and avoids the wrong values that
    appear if you try to infer this after the row deletion step.

    Parameters
    ----------
    proteomics_gene_map : pandas.DataFrame
        Output of :func:`extract_proteomics_gene_map`.
    raw_procan_df : pandas.DataFrame or None
        Raw ProCan frame. None or empty is tolerated.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_name``, ``uniprot_id``, ``procan_uniprot_ids``.
        The ProCan column is present but all-null when no ProCan data was
        available — the column always exists, so downstream code can
        check its content rather than its presence.

    Notes
    -----
    A left merge on the uppercased name, so every proteomics row survives
    and only the ProCan side may be null. Both sides are uppercased for
    the join without either being mutated in the output.

    An all-null ``procan_uniprot_ids`` after a successful merge means the
    names did not match, which is a different problem from ProCan being
    absent — ``02_data_harmonisation.py`` warns about exactly this case.
    """
    if proteomics_gene_map is None:
        return pd.DataFrame(columns=["gene_name", "uniprot_id", "procan_uniprot_ids"])
    if raw_procan_df is None or raw_procan_df.empty:
        proteomics_gene_map = proteomics_gene_map.copy()
        proteomics_gene_map["procan_uniprot_ids"] = None
        return proteomics_gene_map

    procan_map = extract_procan_gene_map(raw_procan_df)
    if procan_map.empty:
        proteomics_gene_map = proteomics_gene_map.copy()
        proteomics_gene_map["procan_uniprot_ids"] = None
        return proteomics_gene_map

    # case-insensitive merge: proteomics names are lowercase, procan names are uppercase
    left = proteomics_gene_map.copy()
    left["_key"] = left["gene_name"].str.upper()
    right = procan_map.copy()
    right["_key"] = right["gene_name"].str.upper()
    merged = left.merge(
        right[["_key", "procan_uniprot_ids"]],
        on="_key",
        how="left",
    ).drop(columns=["_key"])
    return merged.reindex(columns=["gene_name", "uniprot_id", "procan_uniprot_ids"])

# ============================================================
# base_clean — shared first pass for every dataset
# ============================================================

def base_clean(df: pd.DataFrame) -> pd.DataFrame:
    """
    Shared cleaning applied to every dataset:
      - lowercase column names, strip/collapse whitespace in them
      - strip/collapse whitespace in string cell values
      - strip Ensembl version suffixes anywhere they appear

    Cell VALUES retain their original capitalization. Only headers
    are lowercased. Case-insensitive matching is done at comparison
    time (see filter_non_human_rows), never by mutating the data.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to clean. Not modified — a cleaned copy is returned.

    Returns
    -------
    pandas.DataFrame
        Copy with normalised headers and string values.

    Notes
    -----
    Also the fallback cleaner for any dataset without a specific one, so
    a newly added source gets the normalisation without needing setup.

    Value cleaning applies only to object-dtype columns, and each is
    processed with three chained ``apply`` calls — the slowest part of
    cleaning on a wide table.

    Stripping Ensembl versions from *values* means a column of versioned
    IDs is silently unversioned. That is what makes IDs comparable across
    sources, but it does discard the version, which is not recoverable
    afterwards.
    """
    df = df.copy()

    df.columns = (
        pd.Index(df.columns)
        .astype(str)
        .str.strip()
        .str.lower()
        .str.replace(r"\s+", "_", regex=True)
    )

    obj_cols = df.select_dtypes(include="object").columns
    for col in obj_cols:
        df[col] = (
            df[col]
            .apply(lambda v: v.strip() if isinstance(v, str) else v)
            .apply(lambda v: re.sub(r"\s+", " ", v) if isinstance(v, str) else v)
            .apply(strip_ensembl_version)
        )

    return df


# ============================================================
# Per-dataset cleaning functions
# ============================================================

def clean_hpa_rna(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the HPA RNA table.
      'gene' -> 'gene_id' (version suffix stripped first)
      'cell_line' -> 'cell_line_name'

    Parameters
    ----------
    df : pandas.DataFrame
        Raw HPA RNA table.

    Returns
    -------
    pandas.DataFrame
        Cleaned copy. Renames are conditional, so a table missing either
        column passes through without error.
    """
    df = base_clean(df)
    if "gene" in df.columns:
        df["gene"] = df["gene"].apply(strip_ensembl_version)
    df = df.rename(columns={
        c: n for c, n in {"gene": "gene_id", "cell_line": "cell_line_name"}.items()
        if c in df.columns
    })
    return df


def clean_hpa_desc(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the HPA cell-line description table.
      'cellosaurus_id' -> 'cvcl_id'
      'cell_line' -> 'cell_line_name'

    Both renames bring this table onto the identifier vocabulary the
    cell-line roster bridges on.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw HPA description table.

    Returns
    -------
    pandas.DataFrame
        Cleaned copy.
    """
    df = base_clean(df)
    df = df.rename(columns={
        c: n for c, n in {"cellosaurus_id": "cvcl_id", "cell_line": "cell_line_name"}.items()
        if c in df.columns
    })
    return df


def clean_depmap_expr(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the DepMap expression matrix.

    Gene columns come in as 'genename_(ensgxxxxxxxxxxx)' or
    'genename_(ensgxxxxxxxxxxx.version)' after base_clean lowercases
    and underscores them. Rename each to just its Ensembl ID.

    The profile ID sometimes arrives as the DataFrame's index rather
    than a real column (e.g. if it was read upstream with
    index_col=0) — reset it into a real column first so it isn't
    silently dropped. Then the column with no real header (blank, or
    pandas' 'Unnamed: 0' placeholder) is renamed to 'profile_id'.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw DepMap expression matrix, tens of thousands of columns wide.

    Returns
    -------
    pandas.DataFrame
        Cleaned copy with ``profile_id`` plus one column per Ensembl gene
        ID.

    Notes
    -----
    The index promotion happens *before* ``base_clean``, so the former
    index values get the same normalisation as every other value.

    Column names come out lowercase, since ``base_clean`` lowercased them
    — the harmonisation stage upper-cases them when joining against the
    roster.
    """
    # Promote the index to a column BEFORE base_clean, so its values
    # also get the standard lowercase/strip cleaning, and it isn't lost.
    if not isinstance(df.index, pd.RangeIndex):
        df = df.reset_index()

    df = base_clean(df)

    def rename_col(col: str) -> str:
        """Map one header to profile_id, a bare Ensembl ID, or itself."""
        stripped = str(col).strip()
        if stripped == "" or re.match(r"^unnamed:?_?\d*$", stripped) or stripped == "index":
            return "profile_id"
        if re.search(r"ensg\d+", stripped):
            return extract_ensembl_id(stripped)
        return col

    df = df.rename(columns={col: rename_col(col) for col in df.columns})
    return df


def clean_depmap_profiles(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the DepMap omics profiles table.
      'profileid' -> 'profile_id'
      'modelid' -> 'model_id'

    This is the bridge from expression profiles to cell lines: it is what
    lets ``depmap_expr``'s ``profile_id`` resolve to a ``model_id``.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw DepMap profiles table.

    Returns
    -------
    pandas.DataFrame
        Cleaned copy.
    """
    df = base_clean(df)
    df = df.rename(columns={
        c: n for c, n in {"profileid": "profile_id", "modelid": "model_id"}.items()
        if c in df.columns
    })
    return df

def clean_geo_expr(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the GEO expression table.
    Genes-x-samples -> transpose to samples-x-genes; id column -> 'gsm_id'.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw GEO expression matrix, genes as rows.

    Returns
    -------
    pandas.DataFrame
        Cleaned copy, samples as rows, keyed by ``gsm_id``.

    Notes
    -----
    Carries no species or CVCL column of its own, so the species filter
    cannot reach it directly — see
    :func:`filter_geo_expr_by_kept_gsms`, which must run afterwards.
    """
    df = base_clean(df)
    return transpose_with_id_col(df, "gsm_id")

def filter_geo_expr_by_kept_gsms(geo_expr: pd.DataFrame, geo_info: pd.DataFrame) -> tuple:
    """
    Restrict geo_expr to GSM samples that survive in cleaned geo_info.

    geo_expr is keyed on gsm_id only -- it carries no species or CVCL
    column, so it can't be filtered directly. geo_info is the authority:
    any GSM dropped there (non-human organism_ch1, or a non-human CVCL)
    must also leave the expression matrix, or the two tables disagree on
    which samples exist.

    Returns (geo_expr, n_dropped).

    Parameters
    ----------
    geo_expr : pandas.DataFrame
        Cleaned GEO expression table, keyed by ``gsm_id``.
    geo_info : pandas.DataFrame
        Cleaned GEO info table, already species-filtered.

    Returns
    -------
    geo_expr : pandas.DataFrame
        The restricted table, or the input unchanged when filtering was
        not possible.
    n_dropped : int
        Samples removed.

    Notes
    -----
    Run this AFTER both tables are cleaned — it depends on ``geo_info``
    having already had its non-human rows removed.

    Accession matching is case-insensitive on stripped values, without
    mutating either table. The accession column is resolved by trying
    ``geo_accession``, ``gsm_id``, ``gsm`` in turn.

    Both non-filterable cases — no ``gsm_id``, no accession column —
    return unchanged with a printed message rather than raising, so a
    schema change shows up as unfiltered data with a warning rather than
    a crash. Worth reading those messages: unfiltered here means
    non-human samples stay in the expression matrix.
    """
    if geo_expr is None or geo_info is None:
        return geo_expr, 0
    if "gsm_id" not in geo_expr.columns:
        print("  geo_expr has no gsm_id column - cannot filter")
        return geo_expr, 0

    acc_col = next(
        (c for c in ["geo_accession", "gsm_id", "gsm"] if c in geo_info.columns), None
    )
    if acc_col is None:
        print(f"  geo_info has no accession column - columns: {list(geo_info.columns)}")
        return geo_expr, 0

    kept = set(geo_info[acc_col].dropna().astype(str).str.strip().str.lower())
    have = geo_expr["gsm_id"].astype(str).str.strip().str.lower()

    mask_keep = have.isin(kept)
    n_dropped = int((~mask_keep).sum())

    if n_dropped == 0:
        print(f"  geo_expr: all {len(geo_expr):,} GSMs present in geo_info - nothing dropped")
        return geo_expr, 0

    dropped_ids = sorted(geo_expr.loc[~mask_keep, "gsm_id"].astype(str).str.strip())[:5]
    out = geo_expr[mask_keep].reset_index(drop=True)
    print(f"  geo_expr: {len(geo_expr):,} -> {len(out):,} GSMs "
          f"(dropped {n_dropped:,} not in cleaned geo_info)")
    print(f"      e.g. {dropped_ids}")
    return out, n_dropped


def clean_fusions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the DepMap fusions table.

    'gene1' and 'gene2' come in as 'genename (ensgxxxxxxxxxxx.version)'
    pairs. Split each into a gene-name column and an Ensembl-ID column,
    then drop the original combined column.
      gene1 -> gene1_name, gene1_ens_id
      gene2 -> gene2_name, gene2_ens_id

    Splitting is what makes fusions usable as a gene-name source in
    ``build_gene_roster``, which needs the name and the ID as separate
    fields.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw fusions table.

    Returns
    -------
    pandas.DataFrame
        Cleaned copy with the split columns in place of the combined
        ones.

    Notes
    -----
    Columns are found by prefix, so a suffixed variant
    (``gene1_something``) still matches — first match wins.

    A cell that does not parse keeps its whole value as the name and gets
    a null ID, so nothing is lost, but the ID column may be sparser than
    the name column.
    """
    df = base_clean(df)

    for prefix in ["gene1", "gene2"]:
        col = next((c for c in df.columns if c.startswith(prefix)), None)
        if col is None:
            continue
        split = df[col].apply(split_gene_name_id)
        df[f"{prefix}_name"] = split.apply(lambda t: t[0])
        df[f"{prefix}_ens_id"] = split.apply(lambda t: t[1])
        df = df.drop(columns=[col])

    return df

# ============================================================
# Non-human blocklist
# ============================================================
#
# CVCL accessions only -- names are deliberately NOT collected.
#
# audit_name_overlap.py measured 580 names shared between human and
# non-human Cellosaurus entries: 30 primary-name-to-primary-name
# (e.g. 'ca', 'ham-1', 'me1', 'mk2', '11a'), 220 human-synonym-to-
# non-human-primary, 134 the reverse, 196 synonym-to-synonym. Many are
# 1-3 characters ('m', 'a3', 'b9'). Matching on name would silently drop
# legitimate human rows from hpa_rna, sample_info and others.
#
# CVCL accessions are unique per cell line by construction, so they
# cannot collide. Tables with no CVCL column are left unfiltered here
# and handled downstream by the cell_line_roster model_id join.

#: Non-human CVCL accessions, in their original casing. Populated as a
#: side effect of clean_cellosaurus and read by filter_non_human_rows —
#: which makes cleaning order load-bearing. Empty until cellosaurus has
#: been cleaned, and filter_non_human_rows filters nothing while it is.
NON_HUMAN_CVCLS = set()   # cvcl_id, original casing


# Matches any indication of human origin in Cellosaurus' 'Species of
# origin' field, e.g. 'NCBI_TaxID=9606; ! Homo sapiens (Human)'.
#
# Multi-species rows use '||' as a separator:
#   'NCBI_TaxID=9606; ! Homo sapiens (Human) || NCBI_TaxID=10116; ! Rattus norvegicus (Rat)'
# A row counts as human if ANY segment indicates human, so hybrids and
# human-derived xenograft lines stay in (flag don't drop). They're
# counted separately so the decision stays visible.
#
# IGNORECASE is required: base_clean preserves original capitalisation,
# so the raw text is 'Homo sapiens (Human)', not lowercased.
HUMAN_ORIGIN_RE = re.compile(
    r"""
      \b9606\b                 # NCBI taxonomy ID for Homo sapiens
    | homo\s*sapiens           # 'Homo sapiens', 'Homosapiens'
    | \bh\.\s*sapiens\b        # 'H. sapiens'
    | \bhuman\b                # 'Human', 'human', 'HUMAN'
    | \bhumans\b
    """,
    re.IGNORECASE | re.VERBOSE,
)


def clean_cellosaurus(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the Cellosaurus cell-line dictionary and restrict to human lines.
      'identifier_(cell_line_name)' -> 'cell_line_name'
      'accession_(cvcl_xxxxx)'      -> 'cvcl_id'

    A row is HUMAN if its 'species_of_origin' text contains any human
    marker (see HUMAN_ORIGIN_RE) anywhere in the field. NON-HUMAN means no
    human spelling or taxid appears at all -- a pure 'Rattus norvegicus
    (Rat)' row. Blank species_of_origin counts as non-human and is
    reported separately.

    Non-human CVCL accessions are captured into NON_HUMAN_CVCLS BEFORE the
    rows are dropped, for use by filter_non_human_rows(). Names and
    synonyms are NOT captured -- see the NON_HUMAN_CVCLS comment.

    If species_of_origin is missing entirely (schema drift upstream), ALL
    rows are dropped and a warning is printed -- a loud failure by design,
    so a renamed column can't silently let non-human rows through.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw Cellosaurus table.

    Returns
    -------
    pandas.DataFrame
        Human rows only, with renamed identifier columns. Empty when the
        species column is missing.

    Notes
    -----
    **Must be cleaned first.** Populating ``NON_HUMAN_CVCLS`` is a side
    effect that every subsequent :func:`filter_non_human_rows` call
    depends on. ``01_data_cleaning.py`` enforces the ordering.

    The multi-species rule is deliberately inclusive: a human-plus-other
    hybrid stays in, since dropping human-derived material would be the
    worse error. Those rows are counted separately so the choice stays
    visible in the run output.

    Blank species is treated as non-human, which is the cautious
    direction but does mean an under-annotated human line is lost. The
    count is printed so the cost is visible.
    """
    global NON_HUMAN_CVCLS

    df = base_clean(df)
    df = rename_by_prefix(df, {
        "identifier": "cell_line_name",
        "accession":  "cvcl_id",
    })

    species_col = next((c for c in df.columns if c.startswith("species_of_origin")), None)
    if species_col is None:
        print(f"  [WARNING] no 'species_of_origin' column -- columns: {list(df.columns)}. "
              "Dropping ALL rows.")
        return df.iloc[0:0]

    species = df[species_col].fillna("").astype(str)

    is_human = species.str.contains(HUMAN_ORIGIN_RE, na=False)
    is_multi = species.str.contains(r"\|\|", regex=True, na=False)
    is_blank = species.str.strip() == ""

    NON_HUMAN_CVCLS = set(
        df.loc[~is_human, "cvcl_id"].dropna().astype(str).str.strip()
    ) - {"", "nan"}

    n_before = len(df)
    df = df[is_human].reset_index(drop=True)

    print(f"  clean_cellosaurus: kept {len(df):,} / {n_before:,} human rows "
          f"(dropped {n_before - len(df):,}).")
    print(f"    multi-species rows kept (human + other): {int((is_human & is_multi).sum()):,}")
    print(f"    blank species_of_origin (dropped as non-human): {int(is_blank.sum()):,}")
    print(f"  Blocklist: {len(NON_HUMAN_CVCLS):,} non-human cvcl_ids.")

    return df


def filter_non_human_rows(name: str, df: pd.DataFrame) -> tuple:
    """
    Drop rows whose CVCL accession is in the non-human blocklist.

    Matches on cvcl_id / cellosaurus_id / rrid ONLY -- never on
    cell_line_name. See the NON_HUMAN_CVCLS comment for the measured
    reason (580 human/non-human name collisions in Cellosaurus).

    RRID values ('RRID:CVCL_1234') have their accession extracted before
    comparison. Matching is case-insensitive; neither the blocklist nor
    the DataFrame is mutated, so stored values keep original casing.

    Tables with no CVCL column return unchanged -- their non-human rows
    are excluded later by the cell_line_roster model_id join.

    Returns (df, n_dropped).

    Parameters
    ----------
    name : str
        Dataset name, for the printed report only.
    df : pandas.DataFrame
        Cleaned table to filter.

    Returns
    -------
    df : pandas.DataFrame
        The filtered table, or the input unchanged.
    n_dropped : int
        Rows removed.

    Notes
    -----
    Reads the module-level ``NON_HUMAN_CVCLS``, so
    :func:`clean_cellosaurus` must have run first. An empty blocklist
    means every table passes through unfiltered — the printed
    ``blocklist empty`` line is the signal that ordering went wrong.

    Every table prints one line whichever branch it takes, so the run
    output shows what happened to all of them, not just the filtered
    ones.
    """
    if not NON_HUMAN_CVCLS:
        print(f"  {name:<40} blocklist empty - skipped")
        return df, 0

    CVCL_COLS = ["cvcl_id", "cellosaurus_id", "rrid"]
    cvcl_cols = [c for c in CVCL_COLS if c in df.columns]

    if not cvcl_cols:
        print(f"  {name:<40} NO cvcl column - deferred to roster join")
        return df, 0

    cvcls_lc = {c.lower() for c in NON_HUMAN_CVCLS}

    mask_drop = pd.Series(False, index=df.index)
    hits = {}

    for col in cvcl_cols:
        vals = df[col].fillna("").astype(str).str.strip().str.lower()
        if col == "rrid":
            # 'RRID:CVCL_1234' -> 'cvcl_1234'; non-CVCL RRIDs become "".
            vals = vals.str.extract(r"(cvcl_[a-z0-9]+)", expand=False).fillna("")
        m = vals.isin(cvcls_lc)
        if m.any():
            hits[col] = int(m.sum())
        mask_drop |= m

    checked = ", ".join(cvcl_cols)
    n_dropped = int(mask_drop.sum())

    if n_dropped == 0:
        print(f"  {name:<40} NOT PRESENT  (checked: {checked})")
        return df, 0

    sample_vals = set()
    for col in hits:
        sample_vals.update(df.loc[mask_drop, col].dropna().astype(str).str.strip().head(5))
    sample = sorted(sample_vals)[:5]

    before = len(df)
    df = df[~mask_drop].reset_index(drop=True)

    print(f"  {name:<40} PRESENT      {before:,} -> {len(df):,}  (dropped {n_dropped:,})")
    print(f"      hits by column: {hits}")
    print(f"      e.g. {sample}")

    return df, n_dropped


def clean_metabolomics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the CCLE metabolomics table.
    'ccle_id' already matches target format, left as-is.
    'depmap_id' -> 'model_id'.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw CCLE metabolomics table.

    Returns
    -------
    pandas.DataFrame
        Cleaned copy.

    Notes
    -----
    Carries both ``model_id`` and ``ccle_id``, which makes it one of the
    direct sources the cell-line roster is built from.
    """
    df = base_clean(df)
    if "depmap_id" in df.columns:
        df = df.rename(columns={"depmap_id": "model_id"})
    return df


def clean_mirna(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the CCLE miRNA table.
    miRNAs-x-samples -> transpose to samples-x-miRNAs; id column -> 'ccle_id'.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw CCLE miRNA table, miRNAs as rows.

    Returns
    -------
    pandas.DataFrame
        Cleaned copy, samples as rows, keyed by ``ccle_id``.

    Notes
    -----
    Keyed by ``ccle_id`` rather than ``model_id``, so it depends on the
    roster's name bridge to reach a canonical cell-line identifier.
    """
    df = base_clean(df)
    return transpose_with_id_col(df, "ccle_id")


def clean_mutations(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the somatic mutations table.
    'ensemblgeneid' -> 'gene_id'.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw somatic mutations table.

    Returns
    -------
    pandas.DataFrame
        Cleaned copy.

    Notes
    -----
    This table defines the gene universe: ``build_gene_roster`` takes its
    protein-coding gene IDs from here, so a change to its biotype column
    changes what the whole pipeline considers a gene.
    """
    df = base_clean(df)
    if "ensemblgeneid" in df.columns:
        df = df.rename(columns={"ensemblgeneid": "gene_id"})
    return df

def clean_signature(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the signatures table.
    Normalise the model identifier to `model_id`.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw signatures table.

    Returns
    -------
    pandas.DataFrame
        Cleaned copy.

    Notes
    -----
    Both source spellings are mapped in one rename, which is safe only
    because a table would not carry both. If one ever did, the second
    would overwrite the first.
    """
    df = base_clean(df)

    df = df.rename(columns={
        "modelid": "model_id",
        "depmap_id": "model_id"})

    return df


def clean_sample_info(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the DepMap sample info table.
      'depmap_id' -> 'model_id'
      'ccle_name' -> 'ccle_id'

    The central cell-line metadata table: it supplies ``model_id``,
    names, and the lineage map that the transcriptomics and proteomics
    layers stratify on.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw DepMap sample info table.

    Returns
    -------
    pandas.DataFrame
        Cleaned copy.

    Notes
    -----
    The ``ccle_name -> ccle_id`` rename brings this onto the identifier
    vocabulary used elsewhere, but note the direction: what other sources
    call a CCLE *name* is stored here under ``ccle_id``.
    """
    df = base_clean(df)
    df = df.rename(columns={
        c: n for c, n in {"depmap_id": "model_id", "ccle_name": "ccle_id"}.items()
        if c in df.columns
    })
    return df


def clean_procan_raw_tsv(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the raw ProCan TSV without dropping the actual protein-value columns.

    Raw structure:
      row 0: uniprot_id + UniProt accessions
      row 1: symbol + gene symbols
      row 2: model_name / model_id metadata row
      row 3+: cell-line measurements

    Gene symbols from row 1 become COLUMN HEADERS, so they are lowercased
    (headers only). `ccle_name` VALUES keep their original capitalization —
    downstream joins must compare case-insensitively.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw ProCan frame, read with ``header=None`` so the three
        identifier rows arrive as data.

    Returns
    -------
    pandas.DataFrame
        Measurements only, with ``ccle_name`` plus one column per gene
        symbol. Returned early and unchanged when the frame is empty, has
        fewer than three rows, or already carries ``ccle_name``.

    Notes
    -----
    Idempotent: a frame that already has ``ccle_name`` is returned as-is,
    which prevents a second pass from treating measurement rows as
    headers — the failure that would otherwise turn a whole matrix into
    nonsense.

    The UniProt row is consumed here without being preserved. Capture it
    first with :func:`extract_procan_gene_map` if the accessions are
    needed, which is what ``01_data_cleaning.py`` does.

    Duplicate gene symbols get a numeric suffix so the frame can be
    written; the ``model_id`` column is dropped, since ProCan's is not
    the project's canonical one.
    """
    df = base_clean(df)
    if df is None or df.empty:
        return df.copy() if isinstance(df, pd.DataFrame) else pd.DataFrame()

    raw = df.copy().reset_index(drop=True)

    # Already cleaned (ccle_name column present) — return as-is to avoid
    # treating measurement rows as gene-symbol headers.
    if "ccle_name" in raw.columns:
        return raw

    if raw.shape[0] < 3:
        return raw

    symbol_row = raw.iloc[1].fillna("").astype(str).str.strip()
    data = raw.iloc[3:].copy().reset_index(drop=True)

    BAD = {"nan", "none", "null", "na", ""}
    gene_idx = [
        j for j in range(2, len(symbol_row))
        if str(symbol_row.iloc[j]).strip()
        and str(symbol_row.iloc[j]).strip() not in BAD
    ]

    keep_cols = [j for j in ([0, 1] + gene_idx) if j < data.shape[1]]
    data = data.iloc[:, keep_cols].copy()

    # Keep gene symbol case as-is from the raw file.
    gene_names = [str(symbol_row.iloc[j]).strip() for j in gene_idx]
    gene_names = [g for g in gene_names if g and g.lower() not in BAD]

    if not gene_names:
        return data

    seen = {}
    unique_gene_names = []
    for gname in gene_names:
        key = gname.lower()
        count = seen.get(key, 0) + 1
        seen[key] = count
        unique_gene_names.append(gname if count == 1 else f"{gname}_{count}")

    data.columns = ["ccle_name", "model_id"] + unique_gene_names
    data = data.drop(columns=["model_id"], errors="ignore")

    data["ccle_name"] = data["ccle_name"].astype(str).str.strip()

    return data

def clean_proteomics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the proteomics matrix.

    'unnamed:_0' (blank/default index header) -> 'model_id'.

    Each gene column comes in as 'a0av96_(rbm47)' — a UniProt
    accession paired with a gene symbol, in either order. Renamed to
    just the UniProt ID, so the table ends up as model_id + one
    column per UniProt ID.

    If two raw headers resolve to the same UniProt ID (e.g. two
    different gene symbols both mapped to it in the source file),
    that collision is instead named by its gene symbol alone (e.g.
    'rbm47'), since the UniProt ID can't disambiguate them. If the
    gene symbol ALSO collides, a numeric suffix is appended as a
    last resort so the file can still be written.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw proteomics matrix.

    Returns
    -------
    pandas.DataFrame
        Cleaned copy, ``model_id`` plus one column per protein.

    Notes
    -----
    The gene-symbol half of each header is discarded here. Call
    :func:`extract_proteomics_gene_map` on the *raw* frame first if the
    symbol-to-accession mapping is needed, which the cleaning stage does.

    The collision fallbacks mean a column name is not guaranteed to be a
    UniProt ID: it may be a gene symbol, or a suffixed variant of either.
    Anything joining on these headers should tolerate that.
    """
    df = base_clean(df)

    parsed = {}  # original col -> (gene_name, uniprot_id)
    for col in df.columns:
        if is_blank_header(col):
            parsed[col] = (None, "model_id")
        else:
            gene_name, uniprot_id = split_gene_uniprot(col)
            parsed[col] = (gene_name, uniprot_id if uniprot_id else col)

    uniprot_counts = pd.Series([v[1] for v in parsed.values()]).value_counts()
    dup_uniprot_ids = set(uniprot_counts[uniprot_counts > 1].index) - {"model_id"}

    rename_map = {}
    seen_final = {}
    for col, (gene_name, uniprot_id) in parsed.items():
        if uniprot_id in dup_uniprot_ids and gene_name:
            new_name = gene_name
        else:
            new_name = uniprot_id

        if new_name in seen_final:
            seen_final[new_name] += 1
            new_name = f"{new_name}_{seen_final[new_name]}"
        else:
            seen_final[new_name] = 1

        rename_map[col] = new_name

    df = df.rename(columns=rename_map)
    return df


def clean_geo_info(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the GEO series info table.
      'cellosaurus_id' -> 'cvcl_id'
      'cellline' -> 'cellline_name'
    'geo_accession' kept as-is.

    Also drops non-human samples using GEO's own `organism_ch1` field
    (e.g. 'Homo sapiens' vs 'Mus musculus'). This is independent of the
    Cellosaurus CVCL blocklist: a GEO sample can have a blank or
    unmapped cvcl_id but still declare its organism, so both filters
    are needed and they catch different rows.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw GEO info table.

    Returns
    -------
    pandas.DataFrame
        Human samples only, with renamed identifier columns. Returned
        unfiltered when no organism column is found.

    Notes
    -----
    This table is the species authority for GEO — see
    :func:`filter_geo_expr_by_kept_gsms`, which propagates its decisions
    to the expression matrix, and must run after this.

    A missing organism column warns and returns unfiltered, unlike
    :func:`clean_cellosaurus`, which drops everything. The asymmetry is
    deliberate: the CVCL blocklist still catches mapped GEO samples,
    whereas Cellosaurus has no second line of defence.

    Renames ``cellline`` to ``cellline_name``, not ``cell_line_name`` —
    so it does not match the roster's name column. Worth checking against
    ``build_cell_line_roster``, which bridges on ``cell_line_name``.
    """
    df = base_clean(df)
    df = df.rename(columns={"cellosaurus_id": "cvcl_id", "cellline": "cellline_name"})

    org_col = next((c for c in df.columns if c.startswith("organism_ch")), None)
    if org_col is None:
        print(f"  [WARNING] clean_geo_info: no 'organism_ch1' column - "
              f"columns: {list(df.columns)}. No species filter applied.")
        return df

    organism = df[org_col].fillna("").astype(str)
    is_human = organism.str.contains(HUMAN_ORIGIN_RE, na=False)

    n_before = len(df)
    non_human_orgs = sorted(set(organism[~is_human].str.strip()) - {""})
    df = df[is_human].reset_index(drop=True)

    print(f"  clean_geo_info: kept {len(df):,} / {n_before:,} human samples "
          f"(dropped {n_before - len(df):,} by {org_col})")
    if non_human_orgs:
        print(f"      organisms dropped: {non_human_orgs[:8]}")

    return df


def clean_model_list(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the model list table.

    This file contains a source identifier column named `model_id` that is not
    the same as the harmonised ACH-style `model_id` generated later from the
    roster. To keep the source identifier distinct, rename it to `sidm_id`
    before the roster-based attachment step adds the canonical `model_id`.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw model list table.

    Returns
    -------
    pandas.DataFrame
        Cleaned copy with ``model_id`` renamed to ``sidm_id``.

    Notes
    -----
    The rename is not cosmetic. ``attach_model_id`` skips any table that
    already has a ``model_id`` column, so leaving the source name in
    place would block the canonical ID from ever being attached — and the
    table would silently join on the wrong identifier.
    """
    df = base_clean(df)
    if "model_id" in df.columns:
        df = df.rename(columns={"model_id": "sidm_id"})
    return df


# ============================================================
# Dispatch table + clean_all
# ============================================================

# Maps each raw dataset key (matching its parquet filename from
# 00_data_loading.py) to the cleaning function that should run on it.
# Datasets not listed here fall back to base_clean in clean_all().
CLEANERS = {
    "hpa_rna": clean_hpa_rna,
    "hpa_desc": clean_hpa_desc,
    "depmap_expr": clean_depmap_expr,
    "depmap_profiles": clean_depmap_profiles,
    "geo_expr": clean_geo_expr,
    "fusions": clean_fusions,
    "cellosaurus": clean_cellosaurus,
    "metabolomics": clean_metabolomics,
    "mirna": clean_mirna,
    "mutations": clean_mutations,
    "signatures": clean_signature,
    "sample_info": clean_sample_info,
    "proteomics": clean_proteomics,
    "geo_info": clean_geo_info,
    "protein_matrix_averaged_20250211": clean_procan_raw_tsv,
    "model_list_20260709": clean_model_list,
}


def clean_all(tables: dict) -> dict:
    """
    Apply the appropriate cleaning function to every table in `tables`,
    looked up by dataset key in CLEANERS. Falls back to `base_clean`
    for any dataset without a specific cleaner.

    Note: 01_data_cleaning.py's clean_data() re-implements this same
    dispatch loop directly (with per-dataset try/except + logging), so
    this function is mainly for standalone/ad-hoc use rather than
    being called from the pipeline itself.

    Parameters
    ----------
    tables : dict
        Mapping of dataset key to raw :class:`pandas.DataFrame`.

    Returns
    -------
    dict
        Same keys, cleaned frames.

    Notes
    -----
    Iterates in dict order, so cellosaurus is not guaranteed to be
    cleaned first — which means ``NON_HUMAN_CVCLS`` may be empty when
    later tables are processed. The pipeline's own ``clean_data()``
    handles cellosaurus explicitly first for this reason, and also runs
    :func:`filter_non_human_rows`, which this function does not. So this
    is genuinely for ad-hoc use: its output is cleaned but not
    species-filtered.

    No error handling — one failing cleaner aborts the whole loop.
    """
    cleaned = {}
    for name, df in tables.items():
        cleaner = CLEANERS.get(name, base_clean)
        cleaned[name] = cleaner(df)
    return cleaned