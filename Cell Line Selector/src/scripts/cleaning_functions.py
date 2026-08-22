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
    """Strip the trailing .N version suffix from an Ensembl ID string."""
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
    """
    stripped = str(col).strip()
    return stripped == "" or stripped == "index" or bool(re.match(r"^unnamed:?_?\d*$", stripped))


def split_gene_uniprot(col: str):
    """
    Split a proteomics header into (gene_name, uniprot_id), regardless
    of which side of the parens the raw file puts the UniProt
    accession on. Returns (None, None) for headers that don't match
    the 'x_(y)' pattern (e.g. the model-id column).
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
    """
    # Promote the index to a column BEFORE base_clean, so its values
    # also get the standard lowercase/strip cleaning, and it isn't lost.
    if not isinstance(df.index, pd.RangeIndex):
        df = df.reset_index()

    df = base_clean(df)

    def rename_col(col: str) -> str:
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
    """
    df = base_clean(df)
    if "depmap_id" in df.columns:
        df = df.rename(columns={"depmap_id": "model_id"})
    return df


def clean_mirna(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the CCLE miRNA table.
    miRNAs-x-samples -> transpose to samples-x-miRNAs; id column -> 'ccle_id'.
    """
    df = base_clean(df)
    return transpose_with_id_col(df, "ccle_id")


def clean_mutations(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the somatic mutations table.
    'ensemblgeneid' -> 'gene_id'.
    """
    df = base_clean(df)
    if "ensemblgeneid" in df.columns:
        df = df.rename(columns={"ensemblgeneid": "gene_id"})
    return df

def clean_signature(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the signatures table.
    Normalise the model identifier to `model_id`.
    """
    df = base_clean(df)

    df = df.rename(columns={
        "modelid": "model_id",
        "depmap_id": "model_id"})

    return df


def clean_sample_info(df: pd.DataFrame) -> pd.DataFrame:

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
    """
    cleaned = {}
    for name, df in tables.items():
        cleaner = CLEANERS.get(name, base_clean)
        cleaned[name] = cleaner(df)
    return cleaned

