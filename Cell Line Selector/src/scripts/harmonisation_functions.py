"""
harmonisation_functions.py

Cross-table joining and identifier-reconciliation logic that runs
after cleaning: builds a unified cell-line roster and gene roster
from the cleaned tables, then uses the cell-line roster to attach a
consistent `model_id` to every table that doesn't already have one.

Layout of this file:
  1. Generic helpers (ACH/RRID/CVCL matching, aggregation)
  2. Cell-line roster (build_cell_line_roster)
  3. Gene roster + consistency check (build_gene_roster, check_gene_name_consistency)
  4. model_id attachment (attach_model_id, attach_model_id_everywhere, summarize_model_id_attachment)

Design notes
------------
* Multi-valued by default. A model_id legitimately has several names and
  accessions, and a gene_id several symbols, so roster fields hold
  semicolon-joined unique lists rather than a single chosen value. The
  ``(s)`` suffix in a column name marks this.
* Two normalisation regimes, deliberately not unified. Name matching for
  model_id attachment strips punctuation and case
  (:func:`_normalise_join_key`), so ``A-375`` and ``A375`` match; the
  Cellosaurus accession bridges match exact strings, which under-matches
  rather than risking a wrong accession. Gene symbols use a third form
  (:func:`_normalise_gene_symbol`, uppercase). Which applies is stated
  per function.
* Column names are resolved, not assumed. :func:`_resolve_column` allows
  case drift between sources, so a cleaned lowercase header and an
  original uppercase one both find their column.
* Bridges are inverted rather than re-matched. Once a roster column
  carries external IDs, the reverse mapping is derived from it
  (:func:`_invert_multivalue_column`) so the two directions cannot
  disagree.
* Idempotent where it matters. The COSMIC and model_id attach functions
  return early when their output column is already populated, so a
  re-run does not overwrite or duplicate work.

Order dependencies
------------------
:func:`add_cellosaurus_synonyms_to_cell_line_roster` before
:func:`expand_cvcl_ids_from_names`;
:func:`build_gene_roster` before
:func:`add_missing_proteomics_genes_to_gene_roster` before the COSMIC
bridges; :func:`add_cosmic_cnv_id_to_gene_roster` before
:func:`add_gene_id_to_cosmic_cna` (and likewise for the census pair).
Each is restated on the function itself.
"""

from typing import Optional
import re

import pandas as pd


# ============================================================
# Generic helpers
# ============================================================

# ACH-formatted model ID pattern, after base_clean's lowercasing.
ACH_ID_RE = re.compile(r"^ach-\d+$")

RRID_CVCL_RE = re.compile(r"(cvcl_[a-z0-9]+)", re.IGNORECASE)


def extract_cvcl_from_rrid(value) -> Optional[str]:
    """
    Pull the bare CVCL_xxxx accession out of an RRID-formatted value
    like 'rrid:cvcl_1234' or 'cvcl_1234'. Returns None if no CVCL
    pattern is found (e.g. an RRID pointing at a non-Cellosaurus
    resource).

    Parameters
    ----------
    value : any
        Candidate RRID string. Non-strings return None rather than
        raising, so this can be mapped over a mixed column.

    Returns
    -------
    str or None
        The CVCL token as it appeared, or None.

    Notes
    -----
    Searches rather than anchoring, so a CVCL embedded mid-string is
    still found. Case is preserved from the input — callers that need a
    stable key lowercase it themselves.
    """
    if not isinstance(value, str):
        return None
    match = RRID_CVCL_RE.search(value)
    return match.group(1) if match else None


def find_ach_column(df: pd.DataFrame, match_threshold: float = 0.9) -> Optional[str]:
    """
    Find the column in df whose non-null values are mostly ACH-XXXXXXX
    model IDs, regardless of what the column is actually named.
    match_threshold guards against false positives from a column that
    happens to contain a few ACH-looking strings by chance.

    This is a fallback for tables that carry an ACH-style ID under an
    unexpected column name. If the table already has a column literally
    named 'model_id', prefer checking that directly (see
    attach_model_id) rather than relying on this heuristic — a real
    model_id column with enough NaNs or format drift can fall below
    match_threshold and be missed here.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to inspect. Only object-dtype columns are considered.
    match_threshold : float, optional
        Minimum fraction of non-null values that must match the ACH
        pattern. Default 0.9.

    Returns
    -------
    str or None
        Name of the best-matching column, or None if none cleared the
        threshold.

    Notes
    -----
    The pattern is anchored and expects lowercase, since ``base_clean``
    lowercases values upstream. A table whose IDs were not passed through
    cleaning will not match.
    """
    best_col, best_ratio = None, 0.0
    for col in df.select_dtypes(include="object").columns:
        vals = df[col].dropna().astype(str)
        if vals.empty:
            continue
        ratio = vals.str.match(ACH_ID_RE).mean()
        if ratio > best_ratio:
            best_col, best_ratio = col, ratio
    return best_col if best_ratio >= match_threshold else None


def aggregate_by_key(df: pd.DataFrame, key_col: str, agg_cols: list) -> pd.DataFrame:
    """
    Group by `key_col` and collapse each of `agg_cols` into a sorted,
    unique, semicolon-joined string per key. Rows missing `key_col`
    are dropped before grouping.

    This is the shape every roster field takes: one row per key, with all
    values ever seen for it preserved rather than one arbitrarily chosen.

    Parameters
    ----------
    df : pandas.DataFrame
        Long frame with one row per (key, value) observation.
    key_col : str
        Column to group on.
    agg_cols : list of str
        Columns to collapse. Names absent from ``df`` are skipped.

    Returns
    -------
    pandas.DataFrame
        One row per key. If no requested column was present, the unique
        keys are returned alone.

    Notes
    -----
    Values are sorted for determinism, so the same input always produces
    the same string and reruns diff cleanly. Everything is cast to ``str``
    — a numeric ID becomes its string form.
    """
    df = df.dropna(subset=[key_col])

    def join_unique(s: pd.Series) -> str:
        vals = sorted(set(s.dropna().astype(str)))
        return "; ".join(vals)

    agg = {col: join_unique for col in agg_cols if col in df.columns}
    if not agg:
        return df[[key_col]].drop_duplicates()
    return df.groupby(key_col, as_index=False).agg(agg)


def _split_multivalue_string(value) -> list:
    """Split a semicolon/comma-delimited string or list-like value into items.

    The inverse of :func:`aggregate_by_key`'s join, tolerant of the
    variants that arrive from source tables: already-a-list, semicolon-
    joined, comma-joined, or a mix.

    Parameters
    ----------
    value : any
        String, list, tuple or set. NaN returns an empty list.

    Returns
    -------
    list of str
        Stripped tokens, with empties and the string forms of missing
        values (``nan``, ``none``, ``null``) discarded.
    """
    if pd.isna(value):
        return []
    if isinstance(value, (list, tuple, set)):
        items = value
    else:
        text = str(value).strip()
        if text == "":
            return []
        items = re.split(r"[;,]", text)
    cleaned = []
    for item in items:
        if item is None:
            continue
        token = str(item).strip()
        if token and token.lower() not in {"nan", "none", "null"}:
            cleaned.append(token)
    return cleaned


def _normalise_join_key(value) -> str:
    """Standardise name variations so A-375, A 375 and A375 all match.

    The permissive normalisation, used for cell-line name matching during
    model_id attachment. Contrast with the Cellosaurus accession bridges,
    which match exact strings on purpose.

    Parameters
    ----------
    value : any
        Name to normalise. NaN returns an empty string.

    Returns
    -------
    str
        Lowercase, alphanumeric-only key.
    """
    if pd.isna(value):
        return ""
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def _invert_multivalue_column(df: pd.DataFrame, key_col: str, multivalue_col: str) -> dict:
    """
    Given a DataFrame with `key_col` (e.g. gene_id) and `multivalue_col`
    holding a semicolon-joined list of external IDs (e.g. cosmic_cnv_id),
    build the INVERSE mapping: external_id -> semicolon-joined unique
    key_col values that reference it.

    Deriving the reverse direction from the forward bridge, rather than
    re-matching symbols backwards, guarantees the two directions agree.

    Parameters
    ----------
    df : pandas.DataFrame
        Table carrying the already-populated forward bridge.
    key_col : str
        Column holding the internal key, e.g. ``gene_id``.
    multivalue_col : str
        Column holding the semicolon-joined external IDs.

    Returns
    -------
    dict
        External ID to a semicolon-joined, sorted, unique list of keys.
        Many-to-many in both directions: one external ID can map back to
        several keys.
    """
    sub = df[[key_col, multivalue_col]].dropna(subset=[key_col])
    sub = sub.assign(**{multivalue_col: sub[multivalue_col].map(_split_multivalue_string)})
    sub = sub.explode(multivalue_col)
    sub = sub.dropna(subset=[multivalue_col])
    sub[multivalue_col] = sub[multivalue_col].astype(str).str.strip()
    sub = sub[sub[multivalue_col] != ""]

    return (
        sub.groupby(multivalue_col)[key_col]
        .agg(lambda s: "; ".join(sorted(set(str(v) for v in s.dropna()))))
        .to_dict()
    )


def canonical_ensg_id(value) -> Optional[str]:
    """Return a canonical ENSG ID without version suffixes, otherwise None.

    Parameters
    ----------
    value : any
        Candidate identifier. None and NaN are tolerated.

    Returns
    -------
    str or None
        Uppercase ENSG ID with any ``.N`` version stripped, or None if the
        value is not a well-formed ENSG ID.
    """
    if value is None or pd.isna(value):
        return None
    text = str(value).strip().upper()
    if not text:
        return None
    match = re.fullmatch(r"ENSG\d+(?:\.\d+)?", text)
    return text.split(".", 1)[0] if match else None


def filter_protein_coding_expression_tables(tables: dict, gene_roster: Optional[pd.DataFrame] = None) -> dict:
    """Drop non-protein-coding ENSG IDs from depmap_expr, hpa_rna and geo_expr.

    The roster's ``gene_id`` set defines what counts as protein-coding
    (see :func:`build_gene_roster`). Filtering here rather than downstream
    means the expression matrices carry only scoreable genes, which
    substantially narrows the widest tables in the project.

    Two shapes are handled: ``hpa_rna`` is long and filtered by row on its
    ``gene_id`` column; ``depmap_expr`` and ``geo_expr`` are wide and
    filtered by column, with surviving columns renamed to their canonical
    unversioned ENSG form.

    Parameters
    ----------
    tables : dict
        Cleaned tables. Not mutated — a shallow copy is returned.
    gene_roster : pandas.DataFrame, optional
        Roster to filter against. Falls back to ``tables["gene_roster"]``.

    Returns
    -------
    dict
        Copy of ``tables`` with the three expression tables replaced.
        Returned unchanged if no roster is available, the roster lacks
        ``gene_id``, or no valid ENSG IDs could be extracted from it.

    Notes
    -----
    Where a wide table carries the same gene under several versioned
    columns, only the first is kept — the alternative would be duplicate
    canonical column names. ``profile_id`` and ``gsm_id`` are preserved as
    the tables' key columns; every other non-ENSG column is dropped.
    """
    roster = gene_roster if gene_roster is not None else tables.get("gene_roster")
    if roster is None or roster.empty:
        return tables

    gene_id_col = "gene_id"
    if gene_id_col not in roster.columns:
        return tables

    valid = {
        canonical_ensg_id(v)
        for v in roster[gene_id_col].dropna().astype(str).tolist()
        if canonical_ensg_id(v) is not None
    }
    if not valid:
        return tables

    out = dict(tables)

    hpa = out.get("hpa_rna")
    if hpa is not None and isinstance(hpa, pd.DataFrame) and "gene_id" in hpa.columns:
        hpa = hpa.copy()
        hpa["gene_id"] = hpa["gene_id"].map(canonical_ensg_id)
        hpa = hpa[hpa["gene_id"].isin(valid)].dropna(subset=["gene_id"]).copy()
        out["hpa_rna"] = hpa

    for table_name in ("depmap_expr", "geo_expr"):
        df = out.get(table_name)
        if df is None or not isinstance(df, pd.DataFrame):
            continue
        df = df.copy()
        rename_map = {}
        keep_cols = []
        seen = set()
        for col in df.columns:
            if col in {"profile_id", "gsm_id"}:
                keep_cols.append(col)
                continue
            norm = canonical_ensg_id(col)
            if norm is None:
                continue
            if norm in valid:
                if norm not in seen:
                    keep_cols.append(col)
                    seen.add(norm)
                    rename_map[col] = norm
                else:
                    # keep the first valid versioned occurrence for each canonical gene ID
                    keep_cols = [c for c in keep_cols if c != col]
        filtered = df[keep_cols].copy()
        for old_name, new_name in rename_map.items():
            if old_name in filtered.columns:
                filtered = filtered.rename(columns={old_name: new_name})
        filtered = filtered.loc[:, ~filtered.columns.duplicated()]
        out[table_name] = filtered

    return out

# ============================================================
# Cell-line roster
# ============================================================

def add_cellosaurus_synonyms_to_cell_line_roster(
    cell_line_roster: pd.DataFrame,
    cellosaurus: pd.DataFrame,
    roster_cvcl_col: str = "cvcl_id(s)",
    roster_name_col: str = "cell_line_name(s)",
    cellosaurus_cvcl_col: str = "cvcl_id",
    cellosaurus_synonym_col: str = "synonyms",
) -> pd.DataFrame:
    """
    Add synonym names from Cellosaurus to the cell_line_roster using
    the CVCL bridge.

    The roster stores CVCL IDs as a semicolon-delimited string in
    `cvcl_id(s)`, while Cellosaurus stores one CVCL ID per row and a
    separate `synonyms` field that may itself be semicolon- or comma-
    delimited. We match on cvcl_id and append any synonym names to the
    corresponding roster `cell_line_name(s)` values while preserving
    the existing names.

    This is the accession-to-names direction. Run
    :func:`expand_cvcl_ids_from_names` afterwards for the reverse, so the
    names added here also get a chance to resolve further accessions.

    Parameters
    ----------
    cell_line_roster : pandas.DataFrame
        Roster from :func:`build_cell_line_roster`.
    cellosaurus : pandas.DataFrame
        Cleaned Cellosaurus table.
    roster_cvcl_col, roster_name_col : str, optional
        Roster columns, resolved case-insensitively.
    cellosaurus_cvcl_col, cellosaurus_synonym_col : str, optional
        Cellosaurus columns, resolved case-insensitively.

    Returns
    -------
    pandas.DataFrame
        Copy of the roster with widened ``cell_line_name(s)``. Row count
        is unchanged.

    Raises
    ------
    ValueError
        If either frame is None, or if a required column cannot be
        resolved on either side.

    Notes
    -----
    Accession matching is case-insensitive (keys are lowercased), but the
    synonym text is appended as written. Existing names are kept and come
    first; new ones are appended in Cellosaurus order and deduplicated on
    exact string.
    """
    if cell_line_roster is None:
        raise ValueError("cell_line_roster is required")
    if cellosaurus is None:
        raise ValueError("cellosaurus is required")

    roster_cvcl_col = _resolve_column(cell_line_roster, roster_cvcl_col)
    roster_name_col = _resolve_column(cell_line_roster, roster_name_col)
    cellosaurus_cvcl_col = _resolve_column(cellosaurus, cellosaurus_cvcl_col)
    cellosaurus_synonym_col = _resolve_column(cellosaurus, cellosaurus_synonym_col)

    if roster_cvcl_col is None or roster_name_col is None:
        raise ValueError("cell_line_roster must contain cvcl_id(s) and cell_line_name(s)")
    if cellosaurus_cvcl_col is None:
        raise ValueError("cellosaurus must contain a cvcl_id column")
    if cellosaurus_synonym_col is None:
        raise ValueError("cellosaurus must contain a synonyms column")

    roster = cell_line_roster.copy()

    # Build a lookup from CVCL ID -> list of synonyms from Cellosaurus.
    synonym_lookup = {}
    for _, row in cellosaurus[[cellosaurus_cvcl_col, cellosaurus_synonym_col]].dropna(subset=[cellosaurus_cvcl_col]).iterrows():
        cvcl_id = str(row[cellosaurus_cvcl_col]).strip()
        if not cvcl_id:
            continue
        synonyms = _split_multivalue_string(row[cellosaurus_synonym_col])
        if not synonyms:
            continue
        key = cvcl_id.lower()
        synonym_lookup.setdefault(key, [])
        for synonym in synonyms:
            if synonym not in synonym_lookup[key]:
                synonym_lookup[key].append(synonym)

    def enrich_names(row):
        existing = _split_multivalue_string(row[roster_name_col])
        cvcl_ids = _split_multivalue_string(row[roster_cvcl_col])
        seen = set()
        merged = []
        for value in existing + [v for cvcl_id in cvcl_ids for v in synonym_lookup.get(str(cvcl_id).strip().lower(), [])]:
            value = str(value).strip()
            if not value or value.lower() in {"nan", "none", "null"}:
                continue
            if value not in seen:
                seen.add(value)
                merged.append(value)
        return "; ".join(merged) if merged else row[roster_name_col]

    roster[roster_name_col] = roster.apply(enrich_names, axis=1)
    return roster

def expand_cvcl_ids_from_names(
    cell_line_roster: pd.DataFrame,
    cellosaurus: pd.DataFrame,
    roster_name_col: str = "cell_line_name(s)",
    roster_cvcl_col: str = "cvcl_id(s)",
    cellosaurus_name_col: str = "cell_line_name",
    cellosaurus_cvcl_col: str = "cvcl_id",
    cellosaurus_synonym_col: str = "synonyms",
) -> pd.DataFrame:
    """
    Extend each roster row's cvcl_id(s) with every CVCL accession that
    Cellosaurus associates with ANY of that row's cell_line_name(s).

    build_cell_line_roster() only picks up a CVCL when the source table
    already carried one, or when a Cellosaurus row joined on an exact
    name string. This pass closes the gap by also matching against
    Cellosaurus SYNONYMS, not just primary names.

    Matching is on the exact name string (whitespace-stripped only) --
    no case or punctuation normalisation. 'A-431' and 'A431' are treated
    as different names, so a line registered under one spelling in a
    source table only resolves if Cellosaurus lists that exact spelling.
    This is the conservative choice: it under-matches rather than
    risking a wrong accession.

    Existing cvcl_id(s) are preserved -- new accessions are appended,
    deduplicated, and semicolon-joined. Call AFTER
    add_cellosaurus_synonyms_to_cell_line_roster so the synonyms it
    added also get a chance to resolve accessions.

    Parameters
    ----------
    cell_line_roster : pandas.DataFrame
        Roster, ideally already widened with Cellosaurus synonyms.
    cellosaurus : pandas.DataFrame
        Cleaned Cellosaurus table.
    roster_name_col, roster_cvcl_col : str, optional
        Roster columns, resolved case-insensitively.
    cellosaurus_name_col, cellosaurus_cvcl_col, cellosaurus_synonym_col : str, optional
        Cellosaurus columns. The synonym column is optional; without it,
        only primary names are indexed.

    Returns
    -------
    pandas.DataFrame
        Copy of the roster with widened ``cvcl_id(s)``. Row count is
        unchanged.

    Raises
    ------
    ValueError
        If either frame is None, or if a required column cannot be
        resolved.

    Notes
    -----
    Prints three counts: the size of the name index, how many model_ids
    gained accessions and how many were added, and how many model_ids
    matched more than one CVCL. That last number is the one to watch — a
    model_id resolving to several accessions is expected for re-derived
    lines but is also what a bad name collision looks like.
    """
    if cell_line_roster is None:
        raise ValueError("cell_line_roster is required")
    if cellosaurus is None:
        raise ValueError("cellosaurus is required")

    roster_name_col = _resolve_column(cell_line_roster, roster_name_col)
    roster_cvcl_col = _resolve_column(cell_line_roster, roster_cvcl_col)
    cel_name_col = _resolve_column(cellosaurus, cellosaurus_name_col)
    cel_cvcl_col = _resolve_column(cellosaurus, cellosaurus_cvcl_col)
    cel_syn_col = _resolve_column(cellosaurus, cellosaurus_synonym_col)

    if roster_name_col is None or roster_cvcl_col is None:
        raise ValueError("cell_line_roster must contain cell_line_name(s) and cvcl_id(s)")
    if cel_name_col is None or cel_cvcl_col is None:
        raise ValueError("cellosaurus must contain cell_line_name and cvcl_id")

    # --- exact name -> set of CVCL accessions --------------------------
    name_to_cvcls = {}
    cols = [cel_name_col, cel_cvcl_col] + ([cel_syn_col] if cel_syn_col else [])
    for _, row in cellosaurus[cols].iterrows():
        cvcl = str(row[cel_cvcl_col]).strip()
        if not cvcl or cvcl.lower() == "nan":
            continue

        names = [row[cel_name_col]]
        if cel_syn_col:
            names.extend(_split_multivalue_string(row[cel_syn_col]))

        for nm in names:
            if pd.isna(nm):
                continue
            key = str(nm).strip()
            if key and key.lower() != "nan":
                name_to_cvcls.setdefault(key, set()).add(cvcl)

    print(f"  Cellosaurus name index: {len(name_to_cvcls):,} exact names -> CVCL accessions")

    # --- append matched CVCLs to each roster row -----------------------
    roster = cell_line_roster.copy()
    n_gained, n_ambiguous, total_added = 0, 0, 0

    def expand(row):
        nonlocal n_gained, n_ambiguous, total_added

        existing = _split_multivalue_string(row[roster_cvcl_col])
        seen = set(existing)
        merged = list(existing)

        matched = set()
        for nm in _split_multivalue_string(row[roster_name_col]):
            key = str(nm).strip()
            if key in name_to_cvcls:
                matched |= name_to_cvcls[key]

        added = 0
        for cvcl in sorted(matched):
            if cvcl not in seen:
                seen.add(cvcl)
                merged.append(cvcl)
                added += 1

        if added:
            n_gained += 1
            total_added += added
        if len(matched) > 1:
            n_ambiguous += 1

        return "; ".join(merged) if merged else row[roster_cvcl_col]

    roster[roster_cvcl_col] = roster.apply(expand, axis=1)

    print(f"  cvcl_id(s) expanded: {n_gained:,} / {len(roster):,} model_ids gained "
          f"{total_added:,} new accession(s) via exact name match")
    print(f"  {n_ambiguous:,} model_id(s) matched >1 CVCL")

    return roster

def build_cell_line_roster(tables: dict) -> pd.DataFrame:
    """
    One row per model_id, with every cell_line_name / cvcl_id / ccle_id
    / geo_accession seen for that model_id across all cleaned tables,
    collapsed into semicolon-joined unique lists.

    sample_info / metabolomics / depmap_profiles / fusions carry
    model_id directly. cellosaurus / hpa_desc / geo_info have no
    model_id -- bridged in via a match on cvcl_id first (more
    precise), falling back to cell_line_name.

    geo_info's 'cellosaurus_id' -> cvcl_id, 'cellline' -> cell_line_name
    are already renamed by clean_geo_info before this runs; its
    'geo_accession' column is carried through as its own field.

    stripped_cell_line_name (wherever present -- sample_info as a
    direct table, cellosaurus as a bridge table) is folded in as an
    additional cell_line_name value, so it lands in cell_line_name(s)
    alongside every other name variant.

    Note: a model_id only enters the roster if its table ALSO carries
    at least one of cell_line_name/cvcl_id/ccle_id/geo_accession on
    that row (a direct-table row with model_id and nothing else is
    skipped, since there's nothing to aggregate). If a table ever has
    model_id rows with no other identifier present, those model_ids
    won't appear in the roster at all -- worth checking against your
    data if roster row counts look lower than expected.

    Parameters
    ----------
    tables : dict
        Cleaned tables. Missing tables are skipped, so the roster is
        built from whichever sources are present.

    Returns
    -------
    pandas.DataFrame
        Columns ``model_id``, ``cell_line_name(s)``, ``cvcl_id(s)``,
        ``ccle_name(s)``, ``geo_accession(s)``, in that order. Columns
        with no contributing source are present but empty.

    Notes
    -----
    ``rrid`` columns on direct tables are mined for a CVCL accession via
    :func:`extract_cvcl_from_rrid` and used to fill gaps in ``cvcl_id``,
    never to overwrite an existing value.

    Bridge joins match exact strings, with no normalisation on either
    side — a name spelled differently between a bridge table and a direct
    table will not join. Both bridge directions are applied where
    possible and their results concatenated, so a bridge row can
    contribute through either route.
    """
    id_cols = ["cell_line_name", "cvcl_id", "ccle_id", "geo_accession"]

    direct_frames = []
    for name in ["sample_info", "metabolomics", "depmap_profiles", "fusions"]:
        df = tables.get(name)
        if df is None or "model_id" not in df.columns:
            continue
        df = df.copy()

        if "rrid" in df.columns:
            rrid_as_cvcl = df["rrid"].apply(extract_cvcl_from_rrid)
            if "cvcl_id" in df.columns:
                df["cvcl_id"] = df["cvcl_id"].fillna(rrid_as_cvcl)
            else:
                df["cvcl_id"] = rrid_as_cvcl

        cols = ["model_id"] + [c for c in id_cols if c in df.columns]
        if len(cols) > 1:
            direct_frames.append(df[cols].dropna(subset=["model_id"]))

        if "stripped_cell_line_name" in df.columns:
            extra_names = df[["model_id", "stripped_cell_line_name"]].rename(
                columns={"stripped_cell_line_name": "cell_line_name"}
            ).dropna(subset=["model_id"])
            direct_frames.append(extra_names)

    direct = (
        pd.concat(direct_frames, ignore_index=True)
        if direct_frames
        else pd.DataFrame(columns=["model_id"] + id_cols)
    )

    bridge_frames = []
    for name in ["cellosaurus", "hpa_desc", "geo_info"]:
        df = tables.get(name)
        if df is None:
            continue

        cols = [c for c in id_cols if c in df.columns]
        anchor_cols = [c for c in ["cell_line_name", "cvcl_id"] if c in cols]
        if cols and anchor_cols:
            bridge_frames.append(df[cols].dropna(how="all", subset=anchor_cols))

        if "stripped_cell_line_name" in df.columns:
            other_cols = [c for c in id_cols if c in df.columns and c != "cell_line_name"]
            extra = df[other_cols + ["stripped_cell_line_name"]].rename(
                columns={"stripped_cell_line_name": "cell_line_name"}
            ).dropna(subset=["cell_line_name"])
            bridge_frames.append(extra)

    if bridge_frames:
        bridge = pd.concat(bridge_frames, ignore_index=True)

        bridged_frames = []
        if "cvcl_id" in bridge.columns and "cvcl_id" in direct.columns and direct["cvcl_id"].notna().any():
            cvcl_lookup = direct.dropna(subset=["cvcl_id"])[["model_id", "cvcl_id"]].drop_duplicates()
            via_cvcl = bridge.dropna(subset=["cvcl_id"]).merge(cvcl_lookup, on="cvcl_id", how="inner")
            bridged_frames.append(via_cvcl)

        if "cell_line_name" in bridge.columns and "cell_line_name" in direct.columns and direct["cell_line_name"].notna().any():
            name_lookup = direct.dropna(subset=["cell_line_name"])[["model_id", "cell_line_name"]].drop_duplicates()
            via_name = bridge.dropna(subset=["cell_line_name"]).merge(name_lookup, on="cell_line_name", how="inner")
            bridged_frames.append(via_name)

        if bridged_frames:
            direct = pd.concat([direct] + bridged_frames, ignore_index=True)

    direct = direct.dropna(subset=["model_id"])
    roster = aggregate_by_key(direct, key_col="model_id", agg_cols=id_cols)
    roster = roster.rename(columns={
        "cell_line_name": "cell_line_name(s)",
        "cvcl_id": "cvcl_id(s)",
        "ccle_id": "ccle_name(s)",
        "geo_accession": "geo_accession(s)",
    })
    return roster.reindex(columns=[
        "model_id", "cell_line_name(s)", "cvcl_id(s)", "ccle_name(s)", "geo_accession(s)"
    ])


# ============================================================
# Gene roster + consistency check
# ============================================================

def _find_gene_name_col(df: pd.DataFrame) -> Optional[str]:
    """Look for a gene-symbol column under a few common cleaned names.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to inspect.

    Returns
    -------
    str or None
        First matching column from ``gene_name``, ``hugo_symbol``,
        ``hugosymbol``, ``gene_symbol``, ``symbol``, or None.

    Notes
    -----
    Exact, case-sensitive matching against cleaned lowercase names — use
    :func:`_resolve_column` where case drift is possible.
    """
    for candidate in ["gene_name", "hugo_symbol", "hugosymbol", "gene_symbol", "symbol"]:
        if candidate in df.columns:
            return candidate
    return None


def get_gene_name_pairs(df: pd.DataFrame, id_col: str = "gene_id", name_col: str = "gene_name") -> pd.DataFrame:
    """Extract deduped (id_col, name_col) pairs from a table, if both columns are present.

    Parameters
    ----------
    df : pandas.DataFrame or None
        Source table.
    id_col, name_col : str, optional
        Column names to pull.

    Returns
    -------
    pandas.DataFrame
        Deduplicated pairs, or an empty frame with those columns when the
        table is None or missing either column — so callers can
        concatenate the result unconditionally.
    """
    if df is None or id_col not in df.columns or name_col not in df.columns:
        return pd.DataFrame(columns=[id_col, name_col])
    return df[[id_col, name_col]].dropna(subset=[id_col, name_col]).drop_duplicates()


def _normalise_gene_symbol(value) -> str:
    """Normalize a gene symbol so roster aliases and CNA gene_symbol values match.

    The uppercase regime, used only for the COSMIC bridges. Distinct from
    :func:`_normalise_join_key` (lowercase, punctuation-stripped) — gene
    symbols carry meaningful hyphens, so only spaces are removed here.

    Parameters
    ----------
    value : any
        Symbol to normalise. NaN returns an empty string.

    Returns
    -------
    str
        Uppercase symbol with internal spaces removed.
    """
    if pd.isna(value):
        return ""
    return str(value).strip().upper().replace(" ", "")


def _resolve_column(df: pd.DataFrame, *candidates: str) -> Optional[str]:
    """Return a matching column name, ignoring case and allowing cleaned lowercase names.

    Lets a function accept the original header (``GENE_SYMBOL``) and find
    it after cleaning has lowercased it, without every call site knowing
    which stage the frame came from.

    Parameters
    ----------
    df : pandas.DataFrame or None
        Table to search. None returns None.
    *candidates : str
        Names to try, in preference order.

    Returns
    -------
    str or None
        The actual column name as it appears in ``df``, or None.

    Notes
    -----
    Exact matches are preferred over case-insensitive ones for each
    candidate in turn.
    """
    if df is None:
        return None
    lowered = {str(col).lower(): col for col in df.columns}
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
        if str(candidate).lower() in lowered:
            return lowered[str(candidate).lower()]
    return None


def _explode_gene_name_aliases(df: pd.DataFrame, id_col: str = "gene_id", name_col: str = "gene_name(s)") -> pd.DataFrame:
    """Explode a semicolon-delimited gene_name(s) field into one row per gene alias.

    Splits on both semicolons and commas and normalises each alias, giving
    the long form the COSMIC symbol bridges join on.

    Parameters
    ----------
    df : pandas.DataFrame or None
        Frame carrying the ID and the joined names.
    id_col, name_col : str, optional
        Column names.

    Returns
    -------
    pandas.DataFrame
        Columns ``[id_col, "gene_name"]``, one row per (id, alias),
        deduplicated. Empty with those columns if the input is unusable.
    """
    if df is None or id_col not in df.columns or name_col not in df.columns:
        return pd.DataFrame(columns=[id_col, "gene_name"])

    rows = []
    for _, row in df[[id_col, name_col]].dropna(subset=[id_col, name_col]).iterrows():
        gene_id = row[id_col]
        for token in str(row[name_col]).split(";"):
            for part in token.split(","):
                alias = _normalise_gene_symbol(part)
                if alias:
                    rows.append({id_col: gene_id, "gene_name": alias})

    return pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)


def add_cosmic_cnv_id_to_gene_roster(
    gene_roster: pd.DataFrame,
    cosmic_cna: pd.DataFrame,
    roster_name_col: str = "gene_name(s)",
    cosmic_symbol_col: str = "GENE_SYMBOL",
    cosmic_cnv_col: str = "COSMIC_CNV_ID",
    output_col: str = "cosmic_cnv_id",
) -> pd.DataFrame:
    """
    Add a `cosmic_cnv_id` column to the gene_roster by bridging each gene_id
    through the gene_symbol/gene_name(s) overlap with the COSMIC CNA table.

    The bridge is intentionally many-to-many: a gene symbol may map to several
    CNV IDs, and a roster gene_id may have several gene aliases. The final value
    is a semicolon-joined, unique list of matching COSMIC CNV IDs.

    Parameters
    ----------
    gene_roster : pandas.DataFrame
        Roster with ``gene_id`` and joined gene names.
    cosmic_cna : pandas.DataFrame
        COSMIC CNA table.
    roster_name_col : str, optional
        Roster name column, resolved case-insensitively.
    cosmic_symbol_col, cosmic_cnv_col : str, optional
        COSMIC columns, resolved case-insensitively.
    output_col : str, optional
        Column to write. Default ``"cosmic_cnv_id"``.

    Returns
    -------
    pandas.DataFrame
        Copy of the roster with the output column added. The column is
        all-null when the roster has no aliases or nothing matched.

    Raises
    ------
    ValueError
        If either frame is None, or a required column cannot be resolved.

    Notes
    -----
    Idempotent: if the output column already holds any non-empty value,
    the roster is returned untouched, so a re-run neither overwrites nor
    duplicates. Both sides are normalised with
    :func:`_normalise_gene_symbol` before joining.

    Run before :func:`add_gene_id_to_cosmic_cna`, which inverts this
    column rather than re-matching symbols.
    """
    if gene_roster is None:
        raise ValueError("gene_roster is required")
    if cosmic_cna is None:
        raise ValueError("cosmic_cna is required")

    gene_id_col = _resolve_column(gene_roster, "gene_id")
    roster_name_col = _resolve_column(gene_roster, roster_name_col)
    cosmic_symbol_col = _resolve_column(cosmic_cna, cosmic_symbol_col)
    cosmic_cnv_col = _resolve_column(cosmic_cna, cosmic_cnv_col)

    if gene_id_col is None:
        raise ValueError("gene_roster must contain a gene_id column")
    if cosmic_symbol_col is None or cosmic_cnv_col is None:
        raise ValueError(f"cosmic_cna must contain '{cosmic_symbol_col}' and '{cosmic_cnv_col}'")
    if roster_name_col is None:
        raise ValueError(f"gene_roster must contain '{roster_name_col}'")

    roster = gene_roster.copy()
    if output_col in roster.columns:
        filled = roster[output_col].fillna("").astype(str).str.strip()
        if not filled.eq("").all():
            return roster

    gene_aliases = _explode_gene_name_aliases(roster[[gene_id_col, roster_name_col]], id_col=gene_id_col, name_col=roster_name_col)
    if gene_aliases.empty:
        roster[output_col] = None
        return roster

    gene_aliases = gene_aliases.rename(columns={"gene_name": "gene_symbol"})
    gene_aliases["gene_symbol"] = gene_aliases["gene_symbol"].map(_normalise_gene_symbol)

    cosmic_lookup = cosmic_cna[[cosmic_symbol_col, cosmic_cnv_col]].copy()
    cosmic_lookup[cosmic_symbol_col] = cosmic_lookup[cosmic_symbol_col].map(_normalise_gene_symbol)
    cosmic_lookup = cosmic_lookup.dropna(subset=[cosmic_symbol_col, cosmic_cnv_col]).drop_duplicates()

    bridge = gene_aliases.merge(
        cosmic_lookup.rename(columns={cosmic_symbol_col: "gene_symbol"}),
        on="gene_symbol",
        how="inner",
    )
    if bridge.empty:
        roster[output_col] = None
        return roster

    cosmic_map = (
        bridge.groupby("gene_id")[cosmic_cnv_col]
        .agg(lambda s: "; ".join(sorted(set(str(v) for v in s.dropna()))))
        .rename(output_col)
    )
    roster = roster.merge(cosmic_map, on="gene_id", how="left")
    return roster


def add_cosmic_gene_id_to_gene_roster(
    gene_roster: pd.DataFrame,
    cosmic_gene_census: pd.DataFrame,
    roster_name_col: str = "gene_name(s)",
    cosmic_symbol_col: str = "GENE_SYMBOL",
    cosmic_gene_id_col: str = "COSMIC_GENE_ID",
    output_col: str = "cosmic_gene_id(s)",
) -> pd.DataFrame:
    """
    Add a `cosmic_gene_id(s)` column to the gene_roster using the shared
    gene_symbol <-> gene_name(s) bridge against the COSMIC Cancer Gene Census.

    Similar to the CNA bridge, this is intentionally many-to-many because a
    gene symbol can have multiple aliases and the roster may assemble several
    names per gene_id. We aggregate the matched gene IDs into a unique,
    semicolon-joined list per gene_id.

    Parameters
    ----------
    gene_roster : pandas.DataFrame
        Roster with ``gene_id`` and joined gene names.
    cosmic_gene_census : pandas.DataFrame
        COSMIC Cancer Gene Census table.
    roster_name_col : str, optional
        Roster name column, resolved case-insensitively.
    cosmic_symbol_col, cosmic_gene_id_col : str, optional
        COSMIC columns, resolved case-insensitively.
    output_col : str, optional
        Column to write. Default ``"cosmic_gene_id(s)"``.

    Returns
    -------
    pandas.DataFrame
        Copy of the roster with the output column added.

    Raises
    ------
    ValueError
        If either frame is None, or a required column cannot be resolved.

    Notes
    -----
    Idempotent on an already-populated output column, as with the CNA
    bridge. Run before :func:`add_gene_id_to_cosmic_gene_census`.
    """
    if gene_roster is None:
        raise ValueError("gene_roster is required")
    if cosmic_gene_census is None:
        raise ValueError("cosmic_gene_census is required")

    gene_id_col = _resolve_column(gene_roster, "gene_id")
    roster_name_col = _resolve_column(gene_roster, roster_name_col)
    cosmic_symbol_col = _resolve_column(cosmic_gene_census, cosmic_symbol_col)
    cosmic_gene_id_col = _resolve_column(cosmic_gene_census, cosmic_gene_id_col)

    if gene_id_col is None:
        raise ValueError("gene_roster must contain a gene_id column")
    if cosmic_symbol_col is None or cosmic_gene_id_col is None:
        raise ValueError(f"cosmic_gene_census must contain '{cosmic_symbol_col}' and '{cosmic_gene_id_col}'")
    if roster_name_col is None:
        raise ValueError(f"gene_roster must contain '{roster_name_col}'")

    roster = gene_roster.copy()
    if output_col in roster.columns:
        filled = roster[output_col].fillna("").astype(str).str.strip()
        if not filled.eq("").all():
            return roster

    gene_aliases = _explode_gene_name_aliases(roster[[gene_id_col, roster_name_col]], id_col=gene_id_col, name_col=roster_name_col)
    if gene_aliases.empty:
        roster[output_col] = None
        return roster

    gene_aliases = gene_aliases.rename(columns={"gene_name": "gene_symbol"})
    gene_aliases["gene_symbol"] = gene_aliases["gene_symbol"].map(_normalise_gene_symbol)

    gene_census_lookup = cosmic_gene_census[[cosmic_symbol_col, cosmic_gene_id_col]].copy()
    gene_census_lookup[cosmic_symbol_col] = gene_census_lookup[cosmic_symbol_col].map(_normalise_gene_symbol)
    gene_census_lookup = gene_census_lookup.dropna(subset=[cosmic_symbol_col, cosmic_gene_id_col]).drop_duplicates()

    bridge = gene_aliases.merge(
        gene_census_lookup.rename(columns={cosmic_symbol_col: "gene_symbol"}),
        on="gene_symbol",
        how="inner",
    )
    if bridge.empty:
        roster[output_col] = None
        return roster

    gene_id_map = (
        bridge.groupby("gene_id")[cosmic_gene_id_col]
        .agg(lambda s: "; ".join(sorted(set(str(v) for v in s.dropna()))))
        .rename(output_col)
    )
    roster = roster.merge(gene_id_map, on="gene_id", how="left")
    return roster


def add_gene_id_to_cosmic_cna(
    cosmic_cna: pd.DataFrame,
    gene_roster: pd.DataFrame,
    roster_gene_id_col: str = "gene_id",
    roster_cosmic_cnv_col: str = "cosmic_cnv_id",
    cosmic_cnv_col: str = "COSMIC_CNV_ID",
    output_col: str = "gene_id",
) -> pd.DataFrame:
    """
    Add a `gene_id` column to the COSMIC CNA table by INVERTING the
    already-populated gene_roster[roster_cosmic_cnv_col] bridge,
    rather than re-matching gene symbols in the reverse direction.

    This function keeps the project's canonical `gene_id` output, while
    also mirroring the value into `ensg_id` if that column is present
    on the table. That keeps older code paths expecting an ENSEMBL-style
    `ensg_id` column working without breaking the newer `gene_id` naming.

    Parameters
    ----------
    cosmic_cna : pandas.DataFrame
        COSMIC CNA table to annotate.
    gene_roster : pandas.DataFrame
        Roster whose CNV bridge column is already populated.
    roster_gene_id_col, roster_cosmic_cnv_col : str, optional
        Roster columns, resolved case-insensitively.
    cosmic_cnv_col : str, optional
        CNV ID column on the CNA table.
    output_col : str, optional
        Column to write. Default ``"gene_id"``.

    Returns
    -------
    pandas.DataFrame
        Copy of the CNA table with the output column added, and
        ``ensg_id`` gap-filled where that column exists.

    Raises
    ------
    ValueError
        If either frame is None, or a required column cannot be resolved.
        A missing roster bridge column names
        :func:`add_cosmic_cnv_id_to_gene_roster` in the message, since
        that is the fix.

    Notes
    -----
    Inverting the forward bridge guarantees the two directions cannot
    disagree — the alternative, re-matching symbols backwards, could
    produce an asymmetric mapping. Idempotent on an already-populated
    output column, though the ``ensg_id`` mirror is still applied in that
    branch.
    """
    if cosmic_cna is None:
        raise ValueError("cosmic_cna is required")
    if gene_roster is None:
        raise ValueError("gene_roster is required")

    roster_gene_id_col = _resolve_column(gene_roster, roster_gene_id_col)
    roster_cosmic_cnv_col = _resolve_column(gene_roster, roster_cosmic_cnv_col)
    cosmic_cnv_col = _resolve_column(cosmic_cna, cosmic_cnv_col)

    if roster_gene_id_col is None:
        raise ValueError("gene_roster must contain a gene_id column")
    if roster_cosmic_cnv_col is None:
        raise ValueError(
            f"gene_roster must contain '{roster_cosmic_cnv_col}' "
            f"(run add_cosmic_cnv_id_to_gene_roster first)"
        )
    if cosmic_cnv_col is None:
        raise ValueError(f"cosmic_cna must contain '{cosmic_cnv_col}'")

    df = cosmic_cna.copy()
    if output_col in df.columns:
        filled = df[output_col].fillna("").astype(str).str.strip()
        if not filled.eq("").all():
            # keep existing values but mirror into legacy ENSEMBL alias when present
            if "ensg_id" in df.columns:
                df["ensg_id"] = df["ensg_id"].fillna(df[output_col])
            return df

    cnv_to_gene_ids = _invert_multivalue_column(gene_roster, roster_gene_id_col, roster_cosmic_cnv_col)
    if not cnv_to_gene_ids:
        df[output_col] = None
        if "ensg_id" in df.columns:
            df["ensg_id"] = None
        return df

    lookup_keys = df[cosmic_cnv_col].astype(str).str.strip()
    df[output_col] = lookup_keys.map(cnv_to_gene_ids)
    if "ensg_id" in df.columns:
        df["ensg_id"] = df["ensg_id"].fillna(df[output_col])
    return df


def add_ensg_id_to_cosmic_cna(
    cosmic_cna: pd.DataFrame,
    gene_roster: pd.DataFrame,
    roster_gene_id_col: str = "gene_id",
    roster_cosmic_cnv_col: str = "cosmic_cnv_id",
    cosmic_cnv_col: str = "COSMIC_CNV_ID",
    output_col: str = "ensg_id",
) -> pd.DataFrame:
    """Backward-compatible wrapper for older code that expects `ensg_id`.

    Delegates to :func:`add_gene_id_to_cosmic_cna` with the output column
    renamed. See that function for the full parameter and return
    description.

    Returns
    -------
    pandas.DataFrame
        Copy of the CNA table with an ``ensg_id`` column.
    """
    return add_gene_id_to_cosmic_cna(
        cosmic_cna,
        gene_roster,
        roster_gene_id_col=roster_gene_id_col,
        roster_cosmic_cnv_col=roster_cosmic_cnv_col,
        cosmic_cnv_col=cosmic_cnv_col,
        output_col=output_col,
    )


def add_gene_id_to_cosmic_gene_census(
    cosmic_gene_census: pd.DataFrame,
    gene_roster: pd.DataFrame,
    roster_gene_id_col: str = "gene_id",
    roster_cosmic_gene_id_col: str = "cosmic_gene_id(s)",
    cosmic_gene_id_col: str = "COSMIC_GENE_ID",
    output_col: str = "gene_id",
) -> pd.DataFrame:
    """
    Add a `gene_id` column to the COSMIC Cancer Gene Census table by
    inverting the already-populated gene_roster[roster_cosmic_gene_id_col]
    bridge -- same approach as add_gene_id_to_cosmic_cna; see its
    docstring for why this avoids re-matching gene symbols.

    Mirror the value into `ensg_id` when that column exists so older
    projects expecting the ENSEMBL-style name still work.

    Parameters
    ----------
    cosmic_gene_census : pandas.DataFrame
        Census table to annotate.
    gene_roster : pandas.DataFrame
        Roster whose census bridge column is already populated.
    roster_gene_id_col, roster_cosmic_gene_id_col : str, optional
        Roster columns, resolved case-insensitively.
    cosmic_gene_id_col : str, optional
        COSMIC gene ID column on the census table.
    output_col : str, optional
        Column to write. Default ``"gene_id"``.

    Returns
    -------
    pandas.DataFrame
        Copy of the census table with the output column added.

    Raises
    ------
    ValueError
        If either frame is None, or a required column cannot be resolved.
        A missing roster bridge column names
        :func:`add_cosmic_gene_id_to_gene_roster` in the message.
    """
    if cosmic_gene_census is None:
        raise ValueError("cosmic_gene_census is required")
    if gene_roster is None:
        raise ValueError("gene_roster is required")

    roster_gene_id_col = _resolve_column(gene_roster, roster_gene_id_col)
    roster_cosmic_gene_id_col = _resolve_column(gene_roster, roster_cosmic_gene_id_col)
    cosmic_gene_id_col = _resolve_column(cosmic_gene_census, cosmic_gene_id_col)

    if roster_gene_id_col is None:
        raise ValueError("gene_roster must contain a gene_id column")
    if roster_cosmic_gene_id_col is None:
        raise ValueError(
            f"gene_roster must contain '{roster_cosmic_gene_id_col}' "
            f"(run add_cosmic_gene_id_to_gene_roster first)"
        )
    if cosmic_gene_id_col is None:
        raise ValueError(f"cosmic_gene_census must contain '{cosmic_gene_id_col}'")

    df = cosmic_gene_census.copy()
    if output_col in df.columns:
        filled = df[output_col].fillna("").astype(str).str.strip()
        if not filled.eq("").all():
            if "ensg_id" in df.columns:
                df["ensg_id"] = df["ensg_id"].fillna(df[output_col])
            return df

    gene_id_to_gene_ids = _invert_multivalue_column(gene_roster, roster_gene_id_col, roster_cosmic_gene_id_col)
    if not gene_id_to_gene_ids:
        df[output_col] = None
        if "ensg_id" in df.columns:
            df["ensg_id"] = None
        return df

    lookup_keys = df[cosmic_gene_id_col].astype(str).str.strip()
    df[output_col] = lookup_keys.map(gene_id_to_gene_ids)
    if "ensg_id" in df.columns:
        df["ensg_id"] = df["ensg_id"].fillna(df[output_col])
    return df


def add_ensg_id_to_cosmic_gene_census(
    cosmic_gene_census: pd.DataFrame,
    gene_roster: pd.DataFrame,
    roster_gene_id_col: str = "gene_id",
    roster_cosmic_gene_id_col: str = "cosmic_gene_id(s)",
    cosmic_gene_id_col: str = "COSMIC_GENE_ID",
    output_col: str = "ensg_id",
) -> pd.DataFrame:
    """Backward-compatible wrapper for older code that expects `ensg_id`.

    Delegates to :func:`add_gene_id_to_cosmic_gene_census` with the output
    column renamed. See that function for the full description.

    Returns
    -------
    pandas.DataFrame
        Copy of the census table with an ``ensg_id`` column.
    """
    return add_gene_id_to_cosmic_gene_census(
        cosmic_gene_census,
        gene_roster,
        roster_gene_id_col=roster_gene_id_col,
        roster_cosmic_gene_id_col=roster_cosmic_gene_id_col,
        cosmic_gene_id_col=cosmic_gene_id_col,
        output_col=output_col,
    )


def check_gene_name_consistency(tables: dict) -> pd.DataFrame:
    """
    Cross-check gene_id -> gene_name pairs between mutations+fusions
    and hpa_rna.

    One gene_id legitimately mapping to several gene_names (aliases,
    paralogs, annotation drift) is expected and NOT flagged here.
    Only gene_ids where the two sources share NO gene_name at all are
    returned -- that's the signature of an actual mismatch (e.g. a
    gene_id pointing at the wrong gene in one of the source files)
    rather than a normal naming alias.

    Read-only diagnostic -- doesn't modify or filter gene_roster.

    Parameters
    ----------
    tables : dict
        Cleaned tables. ``mutations``, ``fusions`` and ``hpa_rna`` are
        used; each may be absent.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_id``, ``names_mutations_fusions``,
        ``names_hpa_rna``, one row per disjoint gene_id. Empty when
        either side has no pairs at all, or when nothing is disjoint.

    Notes
    -----
    Only gene_ids present on *both* sides are compared — a gene known to
    one source alone is not a mismatch. Comparison is on exact strings,
    so a case or spacing difference between sources would read as
    disjoint; treat a large result as a normalisation problem before
    treating it as a data problem.

    Matters for :func:`build_gene_roster`, which aggregates names from
    these same sources without checking for conflict — a disjoint gene_id
    ends up with both names in ``gene_name(s)``.
    """
    mutations = tables.get("mutations")
    fusions = tables.get("fusions")
    hpa_rna = tables.get("hpa_rna")

    name_col = _find_gene_name_col(mutations) if mutations is not None else None
    primary = (
        get_gene_name_pairs(mutations, name_col=name_col)
        if name_col else pd.DataFrame(columns=["gene_id", "gene_name"])
    )
    if fusions is not None:
        for prefix in ["gene1", "gene2"]:
            n_col, i_col = f"{prefix}_name", f"{prefix}_ens_id"
            if n_col in fusions.columns and i_col in fusions.columns:
                extra = fusions[[i_col, n_col]].rename(columns={i_col: "gene_id", n_col: "gene_name"}).dropna()
                primary = pd.concat([primary, extra], ignore_index=True)

    secondary = get_gene_name_pairs(hpa_rna)

    if primary.empty or secondary.empty:
        return pd.DataFrame(columns=["gene_id", "names_mutations_fusions", "names_hpa_rna"])

    primary_sets = primary.groupby("gene_id")["gene_name"].apply(lambda s: set(s.dropna()))
    secondary_sets = secondary.groupby("gene_id")["gene_name"].apply(lambda s: set(s.dropna()))

    mismatches = []
    for gene_id in set(primary_sets.index) & set(secondary_sets.index):
        p_names, s_names = primary_sets[gene_id], secondary_sets[gene_id]
        if p_names.isdisjoint(s_names):
            mismatches.append({
                "gene_id": gene_id,
                "names_mutations_fusions": "; ".join(sorted(p_names)),
                "names_hpa_rna": "; ".join(sorted(s_names)),
            })

    return pd.DataFrame(mismatches)


def build_gene_roster(tables: dict) -> pd.DataFrame:
    """
    One row per protein-coding gene_id found in mutations.

    The mutations table's biotype column defines the gene universe: a
    gene_id is kept if it is flagged protein_coding on ANY of its
    mutation rows, and dropped otherwise. Restricting here is what makes
    the roster usable as the protein-coding filter for the expression
    tables (see filter_protein_coding_expression_tables).

    Fields assembled per gene_id:

    gene_name(s): from a gene-symbol column in mutations if one
    exists, topped up with fusions' gene1/gene2 name<->id pairs and
    hpa_rna's gene_id/gene_name pairs for any genes missing a name
    elsewhere. One gene_id can have multiple gene_names (aliases,
    annotation drift across sources) -- all are kept, not deduped
    down to one. Names that are themselves ENSG IDs are discarded, so a
    source that fell back to the ID as a label doesn't pollute the field.
    geo_gsm_id(s): GSM samples in geo_expr with a non-null value for
    that gene.
    profile_id(s): profiles in depmap_expr with a non-null value for
    that gene.
    uniprot_id(s): UniProt IDs from proteomics, bridged in via a
    gene_name text match (proteomics has no gene_id of its own) --
    see check_gene_name_consistency for how reliable that bridge is
    on your data.
    procan_uniprot_id: preferentially read straight from the raw ProCan
    header rows (row 0 accessions, row 1 symbols), falling back to
    proteomics_gene_map's procan_uniprot_ids column when the raw file
    isn't available.

    Parameters
    ----------
    tables : dict
        Cleaned tables. ``mutations`` is required; ``fusions``,
        ``hpa_rna``, ``geo_expr``, ``depmap_expr``,
        ``proteomics_gene_map`` and ``raw_procan`` are each optional and
        contribute a field when present.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_id``, ``gene_name(s)``, ``geo_gsm_id(s)``,
        ``profile_id(s)``, ``uniprot_id(s)``, ``procan_uniprot_id``.
        Fields with no contributing source are present but empty.

    Raises
    ------
    ValueError
        If ``mutations`` is absent or lacks ``gene_id``, or if the
        biotype column is missing — the message lists the available
        columns, since the expected name is a project assumption worth
        verifying.

    Warns
    -----
    Two module-level constants at the top of this function,
    ``PROTEIN_CODING_COL`` and ``PROTEIN_CODING_VALUE``, are marked in
    the source as needing verification against your cleaned mutations
    table. If either is wrong, the roster comes out empty or wrongly
    populated without any error being raised — check the row count.

    Notes
    -----
    Melting the full expression matrices is memory-heavy on large
    tables -- if depmap_expr/geo_expr are very wide, this may need
    chunking.

    depmap_expr's gene columns are lowercase while roster gene_ids are
    uppercase, so the melted IDs are upper-cased before joining. The
    UniProt bridge upper-cases both sides for the same reason.
    """
    PROTEIN_CODING_COL = "vepbiotype"        # <-- verify against your cleaned mutations columns
    PROTEIN_CODING_VALUE = "protein_coding"  # <-- verify the exact string used

    mutations = tables.get("mutations")
    if mutations is None or "gene_id" not in mutations.columns:
        raise ValueError("mutations table (with gene_id) is required to build gene_roster")
    if PROTEIN_CODING_COL not in mutations.columns:
        raise ValueError(
            f"Expected a '{PROTEIN_CODING_COL}' column in mutations to flag "
            f"protein-coding genes. Available columns: {list(mutations.columns)}"
        )

    # Keep every gene_id -- flag protein-coding status per gene_id
    # rather than dropping non-matches. A gene_id can appear on
    # multiple mutation rows with the same biotype, so this collapses
    # to one flag per gene_id (any() -- true if flagged protein_coding
    # on any row).
    biotype_flags = (
        mutations.dropna(subset=["gene_id"])
        .assign(is_protein_coding=lambda d: d[PROTEIN_CODING_COL] == PROTEIN_CODING_VALUE)
        .groupby("gene_id")["is_protein_coding"]
        .any()
        .reset_index()
    )
    protein_coding_ids = biotype_flags.loc[biotype_flags["is_protein_coding"], "gene_id"]
    roster = pd.DataFrame({"gene_id": protein_coding_ids})

    # --- gene_name(s) ---
    name_col = _find_gene_name_col(mutations)
    name_lookup = (
        mutations[["gene_id", name_col]].rename(columns={name_col: "gene_name"}).dropna(subset=["gene_id"])
        if name_col else pd.DataFrame(columns=["gene_id", "gene_name"])
    )

    fusions = tables.get("fusions")
    if fusions is not None:
        for prefix in ["gene1", "gene2"]:
            n_col, i_col = f"{prefix}_name", f"{prefix}_ens_id"
            if n_col in fusions.columns and i_col in fusions.columns:
                extra = fusions[[i_col, n_col]].rename(columns={i_col: "gene_id", n_col: "gene_name"}).dropna(subset=["gene_id"])
                name_lookup = pd.concat([name_lookup, extra], ignore_index=True)

    hpa_rna = tables.get("hpa_rna")
    if hpa_rna is not None:
        name_lookup = pd.concat([name_lookup, get_gene_name_pairs(hpa_rna)], ignore_index=True)

    name_lookup = name_lookup[
        ~name_lookup["gene_name"].astype(str).str.lower().str.startswith("ensg")
    ]

    gene_names = aggregate_by_key(name_lookup, key_col="gene_id", agg_cols=["gene_name"])
    roster = roster.merge(gene_names, on="gene_id", how="left")

    # --- geo_gsm_id(s) ---
    geo_expr = tables.get("geo_expr")
    if geo_expr is not None and "gsm_id" in geo_expr.columns:
        long = geo_expr.melt(id_vars="gsm_id", var_name="gene_id", value_name="value").dropna(subset=["value"])
        gsm_map = aggregate_by_key(long, key_col="gene_id", agg_cols=["gsm_id"])
        roster = roster.merge(gsm_map, on="gene_id", how="left")

    # --- profile_id(s) ---
    depmap_expr = tables.get("depmap_expr")
    if depmap_expr is not None and "profile_id" in depmap_expr.columns:
        long = depmap_expr.melt(id_vars="profile_id", var_name="gene_id", value_name="value").dropna(subset=["value"])
        long["gene_id"] = long["gene_id"].str.upper()  # depmap cols are lowercase; roster gene_ids are uppercase
        profile_map = aggregate_by_key(long, key_col="gene_id", agg_cols=["profile_id"])
        roster = roster.merge(profile_map, on="gene_id", how="left")

    # --- uniprot_id(s) (via gene_name bridge from proteomics) ---
    # Both sides are normalised to uppercase for a case-insensitive join.
    proteomics_map = tables.get("proteomics_gene_map")
    if proteomics_map is not None and not name_lookup.empty and "uniprot_id" in proteomics_map.columns:
        nl = name_lookup.copy()
        nl["gene_name"] = nl["gene_name"].astype(str).str.upper()
        pm = proteomics_map.copy()
        pm["gene_name"] = pm["gene_name"].astype(str).str.upper()
        bridge = nl.merge(pm[["gene_name", "uniprot_id"]], on="gene_name", how="inner")[["gene_id", "uniprot_id"]]
        uniprot_map = aggregate_by_key(bridge, key_col="gene_id", agg_cols=["uniprot_id"])
        roster = roster.merge(uniprot_map, on="gene_id", how="left")

    # --- procan_uniprot_id (direct from raw ProCan: row 0 = uniprot_id, row 1 = symbol) ---
    raw_procan = tables.get("raw_procan")
    if raw_procan is not None and not name_lookup.empty:
        _META = {"symbol", "model_name", "model_id", "uniprot_id", "nan", ""}
        uniprot_row = raw_procan.iloc[0].fillna("").astype(str).str.strip()
        symbol_row  = raw_procan.iloc[1].fillna("").astype(str).str.strip()
        pc_bridge = pd.DataFrame({"symbol": symbol_row, "procan_uniprot_id": uniprot_row})
        pc_bridge = pc_bridge[
            ~pc_bridge["symbol"].str.lower().isin(_META) &
            ~pc_bridge["procan_uniprot_id"].str.lower().isin(_META)
        ].drop_duplicates()
        # join via name_lookup (gene_id, gene_name) -- both sides already uppercase
        pc_merged = name_lookup.merge(
            pc_bridge, left_on="gene_name", right_on="symbol", how="inner"
        )[["gene_id", "procan_uniprot_id"]]
        procan_map = aggregate_by_key(pc_merged, key_col="gene_id", agg_cols=["procan_uniprot_id"])
        roster = roster.merge(procan_map, on="gene_id", how="left")
    elif proteomics_map is not None and not name_lookup.empty and "procan_uniprot_ids" in proteomics_map.columns:
        # fallback when raw ProCan is not available
        nl = name_lookup.copy()
        nl["gene_name"] = nl["gene_name"].astype(str).str.upper()
        pm_procan = proteomics_map.copy()
        pm_procan["gene_name"] = pm_procan["gene_name"].astype(str).str.upper()
        procan_bridge = nl.merge(
            pm_procan[["gene_name", "procan_uniprot_ids"]],
            on="gene_name",
            how="inner",
        )[["gene_id", "procan_uniprot_ids"]]
        procan_map = aggregate_by_key(procan_bridge, key_col="gene_id", agg_cols=["procan_uniprot_ids"])
        procan_map = procan_map.rename(columns={"procan_uniprot_ids": "procan_uniprot_id"})
        roster = roster.merge(procan_map, on="gene_id", how="left")

    roster = roster.rename(columns={
        "gene_name": "gene_name(s)",
        "gsm_id": "geo_gsm_id(s)",
        "profile_id": "profile_id(s)",
        "uniprot_id": "uniprot_id(s)",
    })

    return roster.reindex(columns=[
        "gene_id", "gene_name(s)", "geo_gsm_id(s)", "profile_id(s)", "uniprot_id(s)", "procan_uniprot_id"
    ])

def build_combined_proteomics_gene_map(
    tables: dict,
    proteomics_key: str = "proteomics_gene_map",
    procan_key: str = "procan_gene_map",
    procan_uniprot_output_col: str = "procan_uniprot_ids",
) -> pd.DataFrame:
    """
    Merge the two SEPARATELY-built protein maps -- proteomics_gene_map
    (gene_name, uniprot_id) from extract_proteomics_gene_map, and
    procan_gene_map (gene_name, uniprot_id) from extract_procan_gene_map
    -- into ONE table with gene_name, uniprot_id, procan_uniprot_ids.

    Neither source table has a procan_uniprot_ids column on its own --
    procan_gene_map's uniprot_id IS the ProCan-side accession, just
    under the wrong name for what build_gene_roster and
    add_missing_proteomics_genes_to_gene_roster expect. This is the
    ONE place that reconciles the naming into the single unified
    table the rest of the pipeline assumes already exists.

    Outer-joined on gene_name (lower/stripped) so a gene present in
    only one source still comes through, with the other side null.

    Raises only if NEITHER source table is available. If just one is
    missing, proceeds with the other alone and prints a warning.

    Parameters
    ----------
    tables : dict
        Cleaned tables containing one or both source maps.
    proteomics_key, procan_key : str, optional
        Keys of the two source maps.
    procan_uniprot_output_col : str, optional
        Name to give the ProCan accession column in the output.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_name``, ``uniprot_id``,
        ``procan_uniprot_ids`` — one row per gene name, accessions
        collapsed to semicolon-joined lists.

    Raises
    ------
    ValueError
        If neither source map is present.

    Notes
    -----
    Gene names are lowercased and stripped before joining, so the two
    maps align regardless of case. Prints the combined row count
    alongside each source's contribution — a combined count close to the
    sum of the two means the join found little overlap, which is worth
    investigating.
    """
    proteomics_map = tables.get(proteomics_key)
    procan_map = tables.get(procan_key)

    if proteomics_map is None and procan_map is None:
        raise ValueError(
            f"Neither '{proteomics_key}' nor '{procan_key}' found in tables -- "
            f"nothing to build a combined protein map from."
        )

    def prep(df: Optional[pd.DataFrame], uniprot_out_col: str) -> pd.DataFrame:
        if df is None or df.empty:
            return pd.DataFrame(columns=["gene_name", uniprot_out_col])
        name_col = _resolve_column(df, "gene_name")
        uni_col = _resolve_column(df, "uniprot_id")
        if name_col is None or uni_col is None:
            return pd.DataFrame(columns=["gene_name", uniprot_out_col])
        out = df[[name_col, uni_col]].rename(columns={name_col: "gene_name", uni_col: uniprot_out_col})
        out["gene_name"] = out["gene_name"].astype(str).str.strip().str.lower()
        out = out[out["gene_name"] != ""]
        return aggregate_by_key(out, key_col="gene_name", agg_cols=[uniprot_out_col])

    if proteomics_map is None:
        print(f"  [WARNING] '{proteomics_key}' not found -- combined map will have no uniprot_id data.")
    left = prep(proteomics_map, "uniprot_id")

    if procan_map is None:
        print(f"  [WARNING] '{procan_key}' not found -- combined map will have no {procan_uniprot_output_col} data.")
    right = prep(procan_map, procan_uniprot_output_col)

    combined = left.merge(right, on="gene_name", how="outer")
    print(f"  Combined protein map: {len(combined):,} unique gene_names "
          f"({left['gene_name'].nunique():,} from {proteomics_key}, "
          f"{right['gene_name'].nunique():,} from {procan_key})")
    return combined

def add_missing_proteomics_genes_to_gene_roster(
    gene_roster: pd.DataFrame,
    proteomics_gene_map: pd.DataFrame,
    roster_name_col: str = "gene_name(s)",
    proteomics_name_col: str = "gene_name",
    proteomics_uniprot_col: str = "uniprot_id",
    proteomics_procan_col: str = "procan_uniprot_ids",
    roster_uniprot_col: str = "uniprot_id(s)",
    roster_procan_col: str = "procan_uniprot_id",
) -> pd.DataFrame:
    """
    Append gene_names from proteomics_gene_map that AREN'T already
    covered anywhere in gene_roster's gene_name(s) column, as new
    rows -- carrying their uniprot_id(s) and procan_uniprot_id along.

    IMPORTANT: `proteomics_gene_map` passed here must be the COMBINED
    protein map (see build_combined_proteomics_gene_map) -- the
    two source tables 01_data_cleaning.py builds separately
    (proteomics_gene_map, procan_gene_map) each only have HALF the
    needed columns. Calling this against just one of them will still
    add rows, but they'll be missing whichever ID came from the
    other source.

    These new rows have NO gene_id -- there's no Ensembl gene_id
    known for a gene that only appears via a protein-level symbol
    match, so gene_id is left null. This is correct, not a gap to
    fill in.

    Matching against "already covered" is a direct lower/stripped
    string comparison on gene_name (both sides already went through
    base_clean's lowercasing) -- NOT the uppercase _normalise_gene_symbol
    used for the COSMIC bridge elsewhere.

    Call this AFTER build_gene_roster but BEFORE any COSMIC-bridge
    calls, so newly-added proteomics-only genes also get a chance to
    match a COSMIC symbol.

    Parameters
    ----------
    gene_roster : pandas.DataFrame
        Roster from :func:`build_gene_roster`.
    proteomics_gene_map : pandas.DataFrame
        The combined protein map. An empty map is a no-op.
    roster_name_col : str, optional
        Roster name column, resolved case-insensitively.
    proteomics_name_col, proteomics_uniprot_col, proteomics_procan_col : str, optional
        Columns on the protein map. The two ID columns are optional —
        whichever resolves is carried across.
    roster_uniprot_col, roster_procan_col : str, optional
        Destination columns on the roster, created if absent.

    Returns
    -------
    pandas.DataFrame
        Roster with the new rows appended. Row order is existing genes
        first, then the additions.

    Raises
    ------
    ValueError
        If either frame is None, or if the roster name column or the
        protein map's name column cannot be resolved.

    Notes
    -----
    Coverage is checked against every alias in ``gene_name(s)``, not just
    a primary name, so a protein whose symbol matches any existing alias
    is treated as covered. New rows are padded to the roster's full
    column set with nulls, so the concatenated frame keeps one schema.

    Prints the resolved column names and the number of rows added — if
    the added count is close to the size of the protein map, the coverage
    check probably failed to match anything and the normalisation is
    worth checking.
    """
    if gene_roster is None:
        raise ValueError("gene_roster is required")
    if proteomics_gene_map is None:
        raise ValueError("proteomics_gene_map is required")
    if proteomics_gene_map.empty:
        print("  proteomics_gene_map is empty -- nothing to backfill.")
        return gene_roster.copy()

    roster_name_col = _resolve_column(gene_roster, roster_name_col)
    proteomics_name_col = _resolve_column(proteomics_gene_map, proteomics_name_col)
    proteomics_uniprot_col = _resolve_column(proteomics_gene_map, proteomics_uniprot_col)
    proteomics_procan_col = _resolve_column(proteomics_gene_map, proteomics_procan_col)

    print(f"  Resolved columns -- roster name: {roster_name_col!r}, "
          f"proteomics name: {proteomics_name_col!r}, "
          f"uniprot: {proteomics_uniprot_col!r}, procan: {proteomics_procan_col!r}")

    if roster_name_col is None:
        raise ValueError(f"gene_roster must contain '{roster_name_col}'")
    if proteomics_name_col is None:
        raise ValueError(f"proteomics_gene_map must contain '{proteomics_name_col}'")

    roster = gene_roster.copy()

    covered = set()
    for value in roster[roster_name_col].dropna():
        for name in _split_multivalue_string(value):
            covered.add(name.strip().lower())

    proteomics = proteomics_gene_map.copy()
    proteomics["_norm_name"] = proteomics[proteomics_name_col].astype(str).str.strip().str.lower()
    proteomics = proteomics[proteomics["_norm_name"] != ""]

    missing = proteomics[~proteomics["_norm_name"].isin(covered)].copy()
    if missing.empty:
        print("  No proteomics-only gene_names to add -- everything already covered by gene_roster.")
        return roster

    missing[proteomics_name_col] = missing["_norm_name"]
    missing = missing.drop(columns=["_norm_name"])

    agg_cols = [c for c in [proteomics_uniprot_col, proteomics_procan_col] if c is not None]
    new_genes = aggregate_by_key(missing, key_col=proteomics_name_col, agg_cols=agg_cols)

    resolved_uniprot_col = _resolve_column(roster, roster_uniprot_col) or roster_uniprot_col
    resolved_procan_col = _resolve_column(roster, roster_procan_col) or roster_procan_col
    for col in [resolved_uniprot_col, resolved_procan_col]:
        if col not in roster.columns:
            roster[col] = None

    new_rows = pd.DataFrame({roster_name_col: new_genes[proteomics_name_col]})
    if proteomics_uniprot_col in new_genes.columns:
        new_rows[resolved_uniprot_col] = new_genes[proteomics_uniprot_col]
    if proteomics_procan_col in new_genes.columns:
        new_rows[resolved_procan_col] = new_genes[proteomics_procan_col]

    for col in roster.columns:
        if col not in new_rows.columns:
            new_rows[col] = None
    new_rows = new_rows[roster.columns]

    print(f"  Adding {len(new_rows):,} proteomics-only gene_name(s) not previously in gene_roster.")
    return pd.concat([roster, new_rows], ignore_index=True)

# ============================================================
# model_id attachment
# ============================================================

def attach_model_id_to_protein_matrix(
    protein_matrix: pd.DataFrame,
    cell_line_roster: pd.DataFrame,
    protein_name_col: str = "ccle_name",
    roster_name_col: str = "cell_line_name(s)",
    output_col: str = "model_id",
) -> pd.DataFrame:
    """
    Attach model_id to a protein matrix by matching the matrix's
    `ccle_name` column against the roster's name aliases. The roster stores
    several aliases for a model_id, including both cell_line_name(s) and
    ccle_name(s), so we match on both after normalising punctuation and spacing.

    One matrix row can match one or several roster model_ids; when there are
    several, we collapse them to a unique semicolon-joined model_id string so
    the output remains one row per protein-matrix row.

    Parameters
    ----------
    protein_matrix : pandas.DataFrame
        Wide protein matrix keyed by cell-line name.
    cell_line_roster : pandas.DataFrame
        Roster with ``model_id`` and its name columns.
    protein_name_col, roster_name_col : str, optional
        Columns to match on, resolved case-insensitively.
    output_col : str, optional
        Column to write. Default ``"model_id"``.

    Returns
    -------
    pandas.DataFrame
        Copy of the matrix with the output column added. Unmatched rows
        get ``pd.NA``. Row count is preserved.

    Raises
    ------
    ValueError
        If either frame is None, if a match column cannot be resolved, or
        if the roster has no usable name columns.

    Notes
    -----
    Row count is preserved deliberately — a wide protein matrix cannot be
    duplicated per ambiguous match without corrupting its shape, so
    ambiguity is collapsed into the value instead. This differs from
    :func:`attach_model_id`, which duplicates rows.

    Idempotent: returns unchanged if the output column exists with any
    non-null value.
    """
    if protein_matrix is None:
        raise ValueError("protein_matrix is required")
    if cell_line_roster is None:
        raise ValueError("cell_line_roster is required")

    protein_name_col = _resolve_column(protein_matrix, protein_name_col)
    roster_name_col = _resolve_column(cell_line_roster, roster_name_col)
    if protein_name_col is None:
        raise ValueError(f"protein_matrix must contain '{protein_name_col}'")
    if roster_name_col is None:
        raise ValueError(f"cell_line_roster must contain '{roster_name_col}'")

    df = protein_matrix.copy()
    if output_col in df.columns and not df[output_col].isna().all():
        return df

    roster_cols = [c for c in [roster_name_col, "ccle_name(s)"] if c in cell_line_roster.columns]
    if not roster_cols:
        raise ValueError("cell_line_roster has no usable name columns for protein-matrix matching")

    roster_lookup = cell_line_roster[["model_id"] + roster_cols].copy()
    roster_lookup = roster_lookup.melt(id_vars=["model_id"], value_vars=roster_cols, var_name="name_source", value_name="name")
    roster_lookup["name"] = roster_lookup["name"].map(_split_multivalue_string)
    roster_lookup = roster_lookup.explode("name").dropna(subset=["name"])
    roster_lookup["name"] = roster_lookup["name"].map(_normalise_join_key)
    roster_lookup = roster_lookup[roster_lookup["name"] != ""].drop_duplicates(subset=["name", "model_id"])

    name_to_model_ids = (
        roster_lookup.groupby("name", sort=False)["model_id"]
        .agg(lambda s: "; ".join(sorted({str(v) for v in s.dropna()})))
        .to_dict()
    )

    def lookup_model_ids(value):
        if pd.isna(value):
            return pd.NA
        key = _normalise_join_key(value)
        return name_to_model_ids.get(key, pd.NA)

    df[output_col] = (
        df[protein_name_col]
        .map(lookup_model_ids)
        .replace({"": pd.NA})
    )
    return df


def attach_model_id_to_model_list(
    model_list: pd.DataFrame,
    cell_line_roster: pd.DataFrame,
    rrid_col: str = "rrid",
    roster_cvcl_col: str = "cvcl_id(s)",
    output_col: str = "model_id",
) -> pd.DataFrame:
    """
    Attach model_id to a model-list table using the RRID column, which may
    contain values like 'RRID:CVCL_1234', 'cvcl_1234', or a string with a
    CVCL token embedded inside.

    The roster stores CVCL IDs in `cvcl_id(s)` as semicolon-delimited values,
    so we explode those and join one-to-many via CVCL ID.

    Parameters
    ----------
    model_list : pandas.DataFrame
        Table carrying RRIDs.
    cell_line_roster : pandas.DataFrame
        Roster with ``model_id`` and ``cvcl_id(s)``.
    rrid_col, roster_cvcl_col : str, optional
        Columns to match on, resolved case-insensitively.
    output_col : str, optional
        Column to write. Default ``"model_id"``.

    Returns
    -------
    pandas.DataFrame
        Copy of the table with the output column added. Unmatched rows
        get ``pd.NA``; row count is preserved.

    Raises
    ------
    ValueError
        If either frame is None, or a required column cannot be resolved.

    Notes
    -----
    Accession matching is on the lowercased CVCL token, so it is
    case-insensitive — unlike the Cellosaurus name bridges, which match
    exact strings. All CVCL tokens in a value are extracted, so an RRID
    field listing several accessions contributes all of their model_ids.
    """
    if model_list is None:
        raise ValueError("model_list is required")
    if cell_line_roster is None:
        raise ValueError("cell_line_roster is required")

    rrid_col = _resolve_column(model_list, rrid_col)
    roster_cvcl_col = _resolve_column(cell_line_roster, roster_cvcl_col)
    if rrid_col is None:
        raise ValueError(f"model_list must contain '{rrid_col}'")
    if roster_cvcl_col is None:
        raise ValueError(f"cell_line_roster must contain '{roster_cvcl_col}'")

    df = model_list.copy()
    if output_col in df.columns and not df[output_col].isna().all():
        return df

    roster_lookup = cell_line_roster[["model_id", roster_cvcl_col]].copy()
    roster_lookup[roster_cvcl_col] = roster_lookup[roster_cvcl_col].map(_split_multivalue_string)
    roster_lookup = roster_lookup.explode(roster_cvcl_col).dropna(subset=[roster_cvcl_col])
    roster_lookup[roster_cvcl_col] = (
        roster_lookup[roster_cvcl_col]
        .astype(str)
        .str.strip()
        .str.lower()
    )
    roster_lookup = roster_lookup.drop_duplicates(subset=[roster_cvcl_col, "model_id"])

    cvcl_to_model_ids = (
        roster_lookup.groupby(roster_cvcl_col, sort=False)["model_id"]
        .agg(lambda s: "; ".join(sorted({str(v) for v in s.dropna()})))
        .to_dict()
    )

    def lookup_model_ids(value):
        if pd.isna(value):
            return pd.NA

        found = []
        for cvcl in re.findall(r"cvcl_[a-z0-9]+", str(value).lower()):
            if cvcl in cvcl_to_model_ids:
                found.extend(str(cvcl_to_model_ids[cvcl]).split("; "))

        matches = sorted({m.strip() for m in found if str(m).strip()})
        return "; ".join(matches) if matches else pd.NA

    df[output_col] = df[rrid_col].map(lookup_model_ids).replace({"": pd.NA})
    return df


def attach_model_id(df: pd.DataFrame, roster: pd.DataFrame, table_name: str) -> tuple:
    """
    Add model_id to df by looking it up in cell_line_roster.

    Tables that already have a column literally named 'model_id' are
    left untouched. Tables without one are checked for an ACH-ID-shaped
    column under some other name (find_ach_column) as a fallback,
    before falling back further to a cell_line_name/cvcl_id/ccle_id
    lookup against the roster.

    When an identifier value matches MORE THAN ONE model_id, the row
    is duplicated -- once per matching model_id -- rather than left
    unassigned, since a genuinely ambiguous cell line legitimately
    belongs to more than one row here. Every such case is still
    logged in `flagged` so it's visible, even though it's kept.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to annotate.
    roster : pandas.DataFrame
        Cell-line roster to look up against.
    table_name : str
        Name recorded on any flagged rows, so a combined flag frame can
        be traced back to its source table.

    Returns
    -------
    df_with_id : pandas.DataFrame
        The table with ``model_id`` attached. Returned unchanged when the
        table already has an ID, has no usable lookup column, or the
        roster has no matching column for the one it has.
    flagged : pandas.DataFrame
        Columns ``table``, ``row_identifier``, ``matched_model_ids`` —
        one row per input row whose identifier was ambiguous. Empty when
        nothing was ambiguous or the table was skipped.

    Notes
    -----
    Row count is NOT preserved: an ambiguous match duplicates the row.
    This is deliberate here and is the opposite of
    :func:`attach_model_id_to_protein_matrix`'s behaviour. Any change in
    row count is reported by :func:`summarize_model_id_attachment`.

    Name matching goes through :func:`_normalise_join_key`, so
    punctuation and case differences do not block a match. The first
    usable lookup column wins, in the order ``cell_line_name``,
    ``cvcl_id``, ``ccle_id``, ``ccle_name``.
    """
    no_flags = pd.DataFrame(columns=["table", "row_identifier", "matched_model_ids"])

    if "model_id" in df.columns or find_ach_column(df) is not None:
        return df, no_flags

    lookup_col = next((c for c in ["cell_line_name", "cvcl_id", "ccle_id", "ccle_name"] if c in df.columns), None)
    if lookup_col is None:
        return df, no_flags

    roster_cols = []
    if lookup_col in {"cell_line_name", "ccle_name", "ccle_id"}:
        roster_cols += [c for c in ["cell_line_name(s)", "ccle_name(s)"] if c in roster.columns]
    if lookup_col == "cvcl_id":
        roster_cols += ["cvcl_id(s)"]
    if not roster_cols:
        return df, no_flags

    long = roster[["model_id"] + roster_cols].copy()
    long = long.melt(id_vars=["model_id"], value_vars=roster_cols, var_name="_source", value_name=lookup_col)
    long[lookup_col] = long[lookup_col].map(_split_multivalue_string)
    long = long.explode(lookup_col).dropna(subset=[lookup_col])
    long[lookup_col] = long[lookup_col].map(_normalise_join_key)
    long = long[long[lookup_col] != ""].drop_duplicates(subset=["model_id", lookup_col])

    df_lookup = df.copy()
    df_lookup["_join_key"] = df[lookup_col].fillna("").map(_normalise_join_key)
    long = long.rename(columns={lookup_col: "_join_key"})

    match_counts = long.groupby("_join_key")["model_id"].nunique()
    ambiguous_values = set(match_counts[match_counts > 1].index)

    flagged_mask = df_lookup["_join_key"].isin(ambiguous_values)
    flagged = pd.DataFrame({
        "table": table_name,
        "row_identifier": df_lookup.loc[flagged_mask, "_join_key"],
        "matched_model_ids": df_lookup.loc[flagged_mask, "_join_key"].map(
            lambda v: "; ".join(sorted(long.loc[long["_join_key"] == v, "model_id"].unique()))
        ),
    })

    df_with_id = df_lookup.merge(long, on="_join_key", how="left").drop(columns=["_join_key"])
    return df_with_id, flagged


def attach_model_id_everywhere(tables: dict, roster: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    """
    Run attach_model_id across every table in `tables`. Returns the
    updated tables dict and one combined DataFrame of all flagged
    ambiguous-identifier rows across all tables.

    Parameters
    ----------
    tables : dict
        Cleaned tables. Not mutated — a new dict is returned.
    roster : pandas.DataFrame
        Cell-line roster to look up against.

    Returns
    -------
    updated : dict
        Every table, annotated where possible and unchanged otherwise.
    combined_flags : pandas.DataFrame
        All ambiguous rows across all tables, with the source table
        named. Empty with the expected columns if nothing was flagged.

    Notes
    -----
    Tables that could not be annotated pass through unchanged rather than
    being dropped, so the returned dict always has the same keys as the
    input. Some tables will have gained rows — see
    :func:`attach_model_id`'s duplication behaviour.
    """
    updated = {}
    all_flags = []
    for name, df in tables.items():
        new_df, flagged = attach_model_id(df, roster, table_name=name)
        updated[name] = new_df
        if not flagged.empty:
            all_flags.append(flagged)

    combined_flags = pd.concat(all_flags, ignore_index=True) if all_flags else pd.DataFrame(
        columns=["table", "row_identifier", "matched_model_ids"]
    )
    return updated, combined_flags


def summarize_model_id_attachment(tables_before: dict, tables_after: dict) -> pd.DataFrame:
    """
    Per table: whether model_id was newly attached, and whether the
    attach merge duplicated any rows (a sign a lookup value matched
    more than one model_id and slipped through as extra rows rather
    than being caught by the ambiguous-value check).

    Parameters
    ----------
    tables_before : dict
        Tables as they were before attachment.
    tables_after : dict
        Tables after attachment. Keys absent here are skipped.

    Returns
    -------
    pandas.DataFrame
        One row per table with ``table``, ``had_ach_before``,
        ``model_id_attached``, ``rows_before``, ``rows_after``,
        ``duplicated_rows``.

    Notes
    -----
    ``duplicated_rows`` is the number to read first. A non-zero value is
    expected wherever ambiguous identifiers were found — cross-check it
    against the flag frame from :func:`attach_model_id_everywhere`. A
    non-zero count with no corresponding flags means rows were duplicated
    by a path the ambiguity check did not catch.

    This requires the caller to have deep-copied the tables before
    attachment; a shallow copy would make both dicts point at the same
    frames and the comparison meaningless.
    """
    rows = []
    for name, before in tables_before.items():
        after = tables_after.get(name)
        if after is None:
            continue
        had_ach_before = find_ach_column(before) is not None
        has_model_id_after = "model_id" in after.columns
        newly_attached = has_model_id_after and not had_ach_before

        row_count_before = len(before)
        row_count_after = len(after)
        duplicated_rows = row_count_after - row_count_before

        rows.append({
            "table": name,
            "had_ach_before": had_ach_before,
            "model_id_attached": newly_attached,
            "rows_before": row_count_before,
            "rows_after": row_count_after,
            "duplicated_rows": duplicated_rows,
        })

    return pd.DataFrame(rows)


def summarize_missing_gene_biotypes(tables: dict) -> dict:
    """Summarize genes missing from the roster by protein-coding status."""
    if not isinstance(tables, dict):
        return {}

    roster = tables.get("gene_roster")
    roster_ids = set()
    if roster is not None and isinstance(roster, pd.DataFrame) and "gene_id" in roster.columns:
        for value in roster["gene_id"].dropna().astype(str):
            norm = canonical_ensg_id(value)
            if norm is not None:
                roster_ids.add(norm)

    biotype_lookup = {}
    mutations = tables.get("mutations")
    if mutations is not None and isinstance(mutations, pd.DataFrame) and "gene_id" in mutations.columns:
        biotype_col = next((c for c in ["vepbiotype", "biotype"] if c in mutations.columns), None)
        if biotype_col is None:
            biotype_col = next((c for c in mutations.columns if "biotype" in str(c).lower()), None)
        if biotype_col is not None:
            for _, row in mutations[["gene_id", biotype_col]].dropna(subset=["gene_id"]).iterrows():
                norm = canonical_ensg_id(row["gene_id"])
                if norm is not None:
                    biotype_lookup[norm] = str(row[biotype_col]).strip().lower()

    def is_protein_coding(gene_id: str) -> bool:
        biotype = biotype_lookup.get(gene_id, "")
        if not biotype:
            return True
        if "protein_coding" in biotype:
            return True
        if any(token in biotype for token in ["lncrna", "non_coding", "pseudo", "processed", "misc_rna", "ncrna"]):
            return False
        return True

    summary = {}
    for table_name in ("depmap_expr", "geo_expr", "hpa_rna"):
        df = tables.get(table_name)
        if df is None or not isinstance(df, pd.DataFrame):
            continue

        missing = []
        seen = set()
        if "gene_id" in df.columns:
            for value in df["gene_id"].dropna().astype(str):
                norm = canonical_ensg_id(value)
                if norm is None or norm in roster_ids or norm in seen:
                    continue
                seen.add(norm)
                missing.append(norm)
        else:
            for col in df.columns:
                if col in {"profile_id", "gsm_id"}:
                    continue
                norm = canonical_ensg_id(col)
                if norm is None or norm in roster_ids or norm in seen:
                    continue
                seen.add(norm)
                missing.append(norm)

        protein_coding = [gid for gid in missing if is_protein_coding(gid)]
        non_protein_coding = [gid for gid in missing if gid not in set(protein_coding)]
        summary[table_name] = {
            "n_missing": len(missing),
            "n_protein_coding": len(protein_coding),
            "n_non_protein_coding": len(non_protein_coding),
            "protein_coding_examples": protein_coding,
        }

    return summary


def summarize_missing_model_species(tables: dict) -> dict:
    """Summarize whether missing model IDs are human, non-human, or unknown."""
    if not isinstance(tables, dict):
        return {}

    roster = tables.get("cell_line_roster")
    roster_ids = set()
    if roster is not None and isinstance(roster, pd.DataFrame) and "model_id" in roster.columns:
        roster_ids = {str(v).strip() for v in roster["model_id"].dropna().astype(str) if str(v).strip()}

    species_lookup = {}
    for table_name in ("cellosaurus", "model_list_20260709", "sample_info"):
        df = tables.get(table_name)
        if df is None or not isinstance(df, pd.DataFrame):
            continue

        model_col = next((c for c in ["model_id", "sample_id"] if c in df.columns), None)
        species_col = next((c for c in ["species_of_origin", "species"] if c in df.columns), None)
        if model_col is None or species_col is None:
            continue

        for _, row in df[[model_col, species_col]].dropna(subset=[model_col]).iterrows():
            mid = str(row[model_col]).strip()
            if mid:
                species_lookup[mid] = row[species_col]

    def classify(value):
        if value is None or pd.isna(value):
            return "unknown"
        text = str(value).strip()
        if not text:
            return "unknown"
        lowered = text.lower()
        if "human" in lowered or "homo sapiens" in lowered or "h. sapiens" in lowered or "9606" in lowered:
            return "human"
        return "non_human"

    summary = {}
    for name, df in tables.items():
        if df is None or not isinstance(df, pd.DataFrame):
            continue

        model_col = next((c for c in ["model_id", "sample_id"] if c in df.columns), None)
        if model_col is None:
            continue

        missing = []
        seen = set()
        for value in df[model_col].dropna().astype(str):
            mid = str(value).strip()
            if not mid or mid in roster_ids or mid in seen:
                continue
            seen.add(mid)
            missing.append(mid)

        stats = {"n_missing": len(missing), "n_human": 0, "n_non_human": 0, "n_unknown": 0}
        for mid in missing:
            status = classify(species_lookup.get(mid))
            stats[f"n_{status}"] += 1
        summary[name] = stats

    return summary