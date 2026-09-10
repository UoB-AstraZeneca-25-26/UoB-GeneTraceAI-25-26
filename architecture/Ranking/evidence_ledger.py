"""
Ranking/evidence_ledger.py
---------------------------
Shared derivation logic for evidence fields. Used by both build_evidence_ledger.py
(bulk) and explain_pair.py (single-pair live query) so the two never drift apart.
"""
import pandas as pd
import numpy as np


def derive_layers_present(row: pd.Series) -> list[str]:
    """
    List which omics layers have real evidence for one (gene, line) row.

    Parameters
    ----------
    row : pandas.Series
        One prediction row. RNA is assumed always present (this ledger
        only covers RNA-scored pairs); protein/mutation/fusion/CNA are
        each added only if that layer's own evidence field indicates a
        real measurement/alteration.

    Returns
    -------
    list of str
        Subset of ``["rna", "protein", "mutation", "fusion", "cna"]``,
        always starting with ``"rna"``.
    """
    layers = ["rna"]
    if pd.notna(row.get("prot_z_t")):
        layers.append("protein")
    if row.get("p_mutation", 0) > 0:
        layers.append("mutation")
    if row.get("p_fusion", 0) > 0:
        layers.append("fusion")
    if row.get("has_cna_alteration", False):
        layers.append("cna")
    return layers


def derive_driver_alteration(row: pd.Series) -> bool:
    """
    Read the has_driver_alteration flag off a row, defaulting to False.

    Parameters
    ----------
    row : pandas.Series
        One prediction row, expected to carry ``has_driver_alteration``
        from :mod:`driver_routing`.

    Returns
    -------
    bool
        ``row["has_driver_alteration"]`` if present, else ``False``.
    """
    return bool(row.get("has_driver_alteration", False))


def derive_cna_role_consistent(row: pd.Series) -> bool:
    """
    Check whether a CNA alteration matches its gene's expected direction.

    An oncogene is "consistent" if it has a copy-number amplification;
    a TSG is "consistent" if it does NOT (a TSG's expected alteration
    is loss, not gain). Any other combination — including genes with
    no CNA alteration at all when a TSG-consistent one would require
    its absence — falls through to False except the two matching
    cases above.

    Parameters
    ----------
    row : pandas.Series
        One prediction row, expected to carry ``gene_role`` and
        ``has_cna_alteration``.

    Returns
    -------
    bool
        True only for (oncogene, has amplification) or (TSG, no
        amplification); False otherwise, including missing/unknown
        ``gene_role``.
    """
    role = str(row.get("gene_role", ""))
    has_amp = row.get("has_cna_alteration", False)
    if "oncogene" in role and has_amp:
        return True
    if "tsg" in role and not has_amp:
        return True
    return False


def derive_essentiality_check(row: pd.Series) -> str:
    """
    Bucket a row's core_score into a coarse essentiality label.

    Parameters
    ----------
    row : pandas.Series
        One prediction row. Missing ``core_score`` defaults to 0.5
        ("possibly_essential") rather than being treated as an error —
        this function is a display-layer bucketing, not a scoring
        decision.

    Returns
    -------
    str
        ``"likely_essential"`` (score >= 0.8), ``"possibly_essential"``
        (0.5 <= score < 0.8), or ``"likely_non_essential"`` (score < 0.5).
    """
    score = row.get("core_score", 0.5)
    if score >= 0.8:
        return "likely_essential"
    if score >= 0.5:
        return "possibly_essential"
    return "likely_non_essential"


def derive_census_role(row: pd.Series, gene_lookup: pd.DataFrame) -> str:
    """
    Look up a row's gene_role directly from the gene lookup table.

    Parameters
    ----------
    row : pandas.Series
        One prediction row; only ``ensg_id`` is read.
    gene_lookup : pandas.DataFrame
        Must contain ``ensg_id`` and ``gene_role`` columns (e.g.
        ``reference/gene_lookup.parquet``).

    Returns
    -------
    str
        ``gene_lookup``'s ``gene_role`` for this row's ``ensg_id``, or
        ``"unknown"`` if the gene isn't found (or its ``gene_role`` is
        itself missing).
    """
    ensg = row.get("ensg_id", "")
    match = gene_lookup[gene_lookup["ensg_id"] == ensg]
    if match.empty:
        return "unknown"
    return str(match.iloc[0].get("gene_role", "unknown"))


def build_ledger_row(pred: pd.DataFrame, core: pd.DataFrame,
                     cna: pd.DataFrame, genes: pd.DataFrame,
                     lines: pd.DataFrame) -> pd.DataFrame:
    """
    Attach the shared evidence-ledger fields to a predictions table.

    Fills in ``has_cna_alteration``/``gene_role`` from ``cna``/``genes``
    if not already present on ``pred``, defaults ``p_mutation``/
    ``p_fusion`` to 0.0 where missing, then derives
    ``layers_present`` (stringified list, via
    :func:`derive_layers_present`), ``driver_alteration`` (via
    :func:`derive_driver_alteration`), and ``essentiality`` (via
    :func:`derive_essentiality_check`) as new columns.

    Parameters
    ----------
    pred : pandas.DataFrame
        Predictions table to annotate (not mutated in place — a copy
        is returned).
    core : pandas.DataFrame
        Unused directly in the current implementation; accepted for a
        consistent call signature with callers that pass it alongside
        the other tables.
    cna : pandas.DataFrame
        Must contain ``ensg_id``, ``model_id``, ``has_cna_alteration``
        — merged in only if ``pred`` doesn't already have that column.
    genes : pandas.DataFrame
        Must contain ``ensg_id``, ``hgnc_symbol`` (renamed to
        ``gene_symbol``), ``gene_role`` — merged in only if ``pred``
        doesn't already have ``gene_role``.
    lines : pandas.DataFrame
        Unused directly in the current implementation; accepted for
        the same call-signature reason as ``core``.

    Returns
    -------
    pandas.DataFrame
        Copy of ``pred`` with ``has_cna_alteration``, ``p_mutation``,
        ``p_fusion`` filled/defaulted, plus the three new derived
        columns described above.
    """
    df = pred.copy()
    if "has_cna_alteration" not in df.columns:
        df = df.merge(cna[["ensg_id","model_id","has_cna_alteration"]],
                      on=["ensg_id","model_id"], how="left")
    gene_meta = genes[["ensg_id","hgnc_symbol","gene_role"]].rename(
        columns={"hgnc_symbol":"gene_symbol"})
    if "gene_role" not in df.columns:
        df = df.merge(gene_meta, on="ensg_id", how="left")

    df["has_cna_alteration"] = df["has_cna_alteration"].fillna(False)
    df["p_mutation"]         = df.get("p_mutation", pd.Series(0.0, index=df.index)).fillna(0.0)
    df["p_fusion"]           = df.get("p_fusion", pd.Series(0.0, index=df.index)).fillna(0.0)

    df["layers_present"]   = df.apply(derive_layers_present, axis=1).apply(str)
    df["driver_alteration"]= df.apply(derive_driver_alteration, axis=1)
    df["essentiality"]     = df.apply(derive_essentiality_check, axis=1)

    return df
