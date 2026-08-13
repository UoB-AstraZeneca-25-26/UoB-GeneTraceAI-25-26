"""
metadata.py
-----------
Everything the warehouse knows about a gene or a cell line, in one place.

The scoring tables carry ids and numbers. The annotation lives in the
harmonisation warehouse (`v_gene`, `v_cell_line`, built by
build_warehouse_views.py). This module joins the two so any consumer -- the CLI,
the JSON export for the web app, a notebook -- describes an entity identically
instead of each re-deriving it.

Degrades gracefully: if the warehouse is absent the loaders return None and
callers fall back to ids alone, so the CLI never hard-fails on a missing
database.
"""

from __future__ import annotations

import os

import pandas as pd

DB_PATH = os.path.join("src", "pipeline", "outputs", "celllineselector.db")

# Field -> human label, in display order. Only fields worth a reader's attention
# are listed; the views carry more.
CELL_LINE_FIELDS = [
    ("canonical_name",        "name"),
    ("stripped_name",         "alt name"),
    ("lineage",               "lineage"),
    ("lineage_subtype",       "lineage subtype"),
    ("primary_disease",       "disease"),
    ("subtype",               "subtype"),
    ("tissue",                "tissue (Sanger)"),
    ("cancer_type",           "cancer type (Sanger)"),
    ("cancer_type_detail",    "cancer detail"),
    ("primary_or_metastasis", "primary/metastasis"),
    ("sample_collection_site", "collection site"),
    ("tissue_status",         "tissue status"),
    ("sex",                   "sex"),
    ("age",                   "age"),
    ("ethnicity",             "ethnicity"),
    ("model_type",            "model type"),
    ("growth_properties",     "growth"),
    ("msi_status",            "MSI status"),
    ("ploidy_wes",            "ploidy (WES)"),
    ("mutational_burden",     "mutational burden"),
    ("source",                "source"),
    ("gap_class",             "identity gap class"),
    ("in_roster",             "in DepMap roster"),
    ("source_wave",           "harmonisation wave"),
]

GENE_FIELDS = [
    ("hugo_symbol",        "symbol"),
    ("gene_full_name",     "name"),
    ("hgnc_id",            "HGNC id"),
    ("locus_type",         "locus type"),
    ("is_protein_coding",  "protein-coding"),
    ("gene_role",          "COSMIC role"),
    ("cgc_tier",           "CGC tier"),
    ("n_uniprot_ids",      "UniProt accessions"),
    ("n_lines_expressed",  "lines with expression"),
    ("median_log2tpm",     "median log2TPM"),
    ("top_vs_median_fold", "top vs median fold"),
    ("signal_spread",      "signal spread"),
]

MODALITY_COLS = ["depmap_expr", "geo_expr", "hpa_rna", "mutations", "fusions",
                 "proteomics", "procan_proteomics", "metabolomics", "mirna",
                 "signatures", "cosmic_cna"]

_CELL = None
_GENE = None


def _connect(db_path=DB_PATH):
    if not os.path.isfile(db_path):
        return None
    try:
        import duckdb
        return duckdb.connect(db_path, read_only=True)
    except Exception:
        return None


def load_cell_lines(db_path=DB_PATH):
    """v_cell_line indexed by model_id, or None if unavailable."""
    global _CELL
    if _CELL is not None:
        return _CELL
    con = _connect(db_path)
    if con is None:
        return None
    try:
        d = con.execute("SELECT * FROM v_cell_line").df()
    except Exception:
        return None
    finally:
        con.close()
    _CELL = d.set_index("model_id")
    return _CELL


def load_genes(db_path=DB_PATH):
    """v_gene indexed by gene_id, or None if unavailable."""
    global _GENE
    if _GENE is not None:
        return _GENE
    con = _connect(db_path)
    if con is None:
        return None
    try:
        d = con.execute("SELECT * FROM v_gene").df()
    except Exception:
        return None
    finally:
        con.close()
    _GENE = d.set_index("gene_id")
    return _GENE


def _clean(v):
    if v is None:
        return None
    if isinstance(v, float) and pd.isna(v):
        return None
    if isinstance(v, (list, tuple)) or hasattr(v, "tolist"):
        try:
            v = list(v)
        except Exception:
            return v
        return v or None
    if pd.isna(v):
        return None
    return v


def cell_line_meta(model_id, cells=None):
    """dict of every annotated field for one cell line. {} if unknown."""
    cells = cells if cells is not None else load_cell_lines()
    if cells is None or model_id not in cells.index:
        return {}
    row = cells.loc[model_id]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    out = {}
    for col, label in CELL_LINE_FIELDS:
        if col in row.index:
            v = _clean(row[col])
            if v is not None and v != "":
                out[label] = v
    out["_modalities"] = [c for c in MODALITY_COLS
                          if c in row.index and bool(row[c])]
    out["_n_modalities"] = _clean(row.get("n_modalities"))
    for idc in ("rrids", "sanger_ids", "cosmic_ids", "ccle_ids"):
        v = _clean(row.get(idc))
        if v:
            out[idc] = v
    return out


def gene_meta(gene_id, genes=None):
    """dict of every annotated field for one gene. {} if unknown."""
    genes = genes if genes is not None else load_genes()
    if genes is None or gene_id not in genes.index:
        return {}
    row = genes.loc[gene_id]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    out = {}
    for col, label in GENE_FIELDS:
        if col in row.index:
            v = _clean(row[col])
            if v is not None and v != "":
                out[label] = v
    aliases = _clean(row.get("gene_names"))
    if aliases and len(aliases) > 1:
        out["aliases"] = aliases
    return out


def format_block(meta, indent="  ", skip_private=True, width=22):
    """Render a metadata dict as aligned 'label: value' lines."""
    lines = []
    for k, v in meta.items():
        if skip_private and k.startswith("_"):
            continue
        if isinstance(v, list):
            v = ", ".join(str(x) for x in v[:6]) + (f" (+{len(v)-6} more)" if len(v) > 6 else "")
        if isinstance(v, float):
            v = f"{v:.4g}"
        lines.append(f"{indent}{k+':':<{width}} {v}")
    return "\n".join(lines)


def modality_line(meta, indent="  "):
    """One line naming which layers actually measured this cell line."""
    mods = meta.get("_modalities") or []
    n = meta.get("_n_modalities")
    if not mods:
        return f"{indent}layers measured:       (none recorded)"
    return (f"{indent}{'layers measured:':<22} {n}/11 — " + ", ".join(mods))
