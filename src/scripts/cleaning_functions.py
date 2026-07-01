"""
cleaning_functions.py
Reusable cleaning functions for each of the 14 multi-omics datasets.
"""

import re
import pandas as pd


def _clean_str_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Lowercase + strip/collapse whitespace on all object/string columns in-place."""
    str_cols = df.select_dtypes(include=["object", "string"]).columns
    for col in str_cols:
        df[col] = (
            df[col]
            .astype(str)
            .str.strip()
            .str.replace(r"\s+", " ", regex=True)
            .str.lower()
        )
        df[col] = df[col].replace("nan", pd.NA)
    return df


_ENS_VERSION = re.compile(r"(ens[gtp]\d+)\.\d+", re.IGNORECASE)
_ENS_CORE = re.compile(r"ens[gtp]\d+", re.IGNORECASE)


def clean_hpa_rna(df: pd.DataFrame) -> pd.DataFrame:
    """Clean HPA RNA cell line table."""
    df = df.copy()
    df.columns = [str(c).lower() for c in df.columns]
    df = _clean_str_cols(df)
    if "gene" in df.columns:
        df["gene"] = df["gene"].astype(str).str.replace(_ENS_VERSION, r"\1", regex=True)
        df["gene"] = df["gene"].replace("nan", pd.NA)
    return df


def clean_depmap_expr(df: pd.DataFrame) -> pd.DataFrame:
    """Clean DepMap expression matrix (rows = PR- IDs, columns = gene headers)."""
    df = df.copy()

    def _clean_header(col):
        c = str(col).lower().strip()
        c = re.sub(r"\s+", " ", c)
        m = _ENS_CORE.search(c)
        if m:
            ens = m.group(0).lower()
            ens = _ENS_VERSION.sub(r"\1", ens)
            return ens
        return c

    df.columns = [_clean_header(c) for c in df.columns]
    df.index = (
        df.index.astype(str).str.strip()
        .str.replace(r"\s+", " ", regex=True)
        .str.lower()
    )
    df = _clean_str_cols(df)
    return df


def clean_geo_expr(df: pd.DataFrame) -> pd.DataFrame:
    """Clean GEO expression table."""
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    df = _clean_str_cols(df)
    if "gene" in df.columns:
        df["gene"] = df["gene"].astype(str).str.replace(_ENS_VERSION, r"\1", regex=True)
        df["gene"] = df["gene"].replace("nan", pd.NA)
    return df


def clean_proteomics(df: pd.DataFrame):
    """
    Clean proteomics table.
    Returns (df_clean, protein_map) where protein_map maps uniprot_id -> gene_symbol.
    """
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    df = df.rename(columns={"unnamed: 0": "depmap_id"})
    df = _clean_str_cols(df)

    pattern = re.compile(r"^\s*([a-z0-9\-]+)\s*(?:\(([^)]*)\))?\s*$", re.IGNORECASE)
    id_cols = ["depmap_id"]
    map_rows, new_names = [], {}

    for col in df.columns:
        if col in id_cols:
            continue
        m = pattern.match(str(col))
        if m:
            uniprot_id = m.group(1).strip().lower()
            gene = m.group(2).strip().lower() if m.group(2) else pd.NA
        else:
            uniprot_id = str(col).strip().lower()
            gene = pd.NA
        new_names[col] = uniprot_id
        map_rows.append({"uniprot_id": uniprot_id, "gene_symbol": gene, "original_header": col})

    df = df.rename(columns=new_names)
    protein_map = pd.DataFrame(map_rows)
    return df, protein_map


def clean_fusions(df: pd.DataFrame) -> pd.DataFrame:
    """Clean fusions table, splitting gene1/gene2 Ensembl ID columns."""
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    df = df.rename(columns={"unnamed: 0": "fusion_index"})
    df = _clean_str_cols(df)

    split_pat = re.compile(r"^(.*?)\s*\(([^)]*)\)\s*$")

    def _split_gene_ens(val):
        if pd.isna(val):
            return (pd.NA, pd.NA)
        m = split_pat.match(str(val).strip())
        if m:
            gene = m.group(1).strip() or pd.NA
            ens = _ENS_VERSION.sub(r"\1", m.group(2).strip())
            ens = ens if ens else pd.NA
            return (gene, ens)
        return (str(val).strip() or pd.NA, pd.NA)

    for src, gene_col, ens_col in [
        ("gene1(ens id)", "gene1", "gene1_ens_id"),
        ("gene2(ens id)", "gene2", "gene2_ens_id"),
    ]:
        if src in df.columns:
            parsed = df[src].apply(_split_gene_ens)
            df[gene_col] = parsed.apply(lambda x: x[0])
            df[ens_col] = parsed.apply(lambda x: x[1])
            df = df.drop(columns=[src])
    return df


def clean_mutations(df: pd.DataFrame) -> pd.DataFrame:
    """Clean somatic mutations table."""
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    df = _clean_str_cols(df)
    for col in [c for c in ["ensemblgeneid", "ensemblfeatureid"] if c in df.columns]:
        df[col] = df[col].astype(str).str.replace(_ENS_VERSION, r"\1", regex=True)
        df[col] = df[col].replace("nan", pd.NA)
    return df


def clean_cellosaurus(df: pd.DataFrame) -> pd.DataFrame:
    """Clean Cellosaurus cell line dictionary."""
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    df = _clean_str_cols(df)
    df = df.rename(columns={
        "identifier (cell line name)": "cellosaurus_cell_line_name",
        "accession (cvcl_xxxx)": "cellosaurus_accession",
    })
    return df


def clean_depmap_profiles(df: pd.DataFrame) -> pd.DataFrame:
    """Clean DepMap omics profiles bridge table."""
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    return _clean_str_cols(df)


def clean_sample_info(df: pd.DataFrame) -> pd.DataFrame:
    """Clean DepMap sample info table."""
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    return _clean_str_cols(df)


def clean_geo_info(df: pd.DataFrame) -> pd.DataFrame:
    """Clean GEO series info table."""
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    return _clean_str_cols(df)


def clean_hpa_desc(df: pd.DataFrame) -> pd.DataFrame:
    """Clean HPA cell line descriptions table."""
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    return _clean_str_cols(df)


def clean_metabolomics(df: pd.DataFrame) -> pd.DataFrame:
    """Clean CCLE metabolomics table."""
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    return _clean_str_cols(df)


def clean_mirna(df: pd.DataFrame) -> pd.DataFrame:
    """Clean CCLE miRNA table."""
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    df = _clean_str_cols(df)
    if "name" in df.columns:
        df["name"] = df["name"].astype(str).str.replace(r"\.0+$", "", regex=True)
        df["name"] = df["name"].replace("nan", pd.NA)
    return df


def clean_signatures(df: pd.DataFrame) -> pd.DataFrame:
    """Clean global omics signatures table."""
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    df = df.rename(columns={"unnamed: 0": "signature_index"})
    return _clean_str_cols(df)


def clean_all(tables: dict):
    """
    Apply all cleaning functions to a dict of raw DataFrames.
    Expects keys matching load_raw_parquets() output.
    Returns a new dict of cleaned DataFrames (protein_map added separately).
    """
    cleaned = {}
    cleaned["hpa_rna"] = clean_hpa_rna(tables["hpa_rna"])
    cleaned["depmap_expr"] = clean_depmap_expr(tables["depmap_expr"])
    cleaned["geo_expr"] = clean_geo_expr(tables["geo_expr"])
    cleaned["proteomics"], cleaned["protein_map"] = clean_proteomics(tables["proteomics"])
    cleaned["fusions"] = clean_fusions(tables["fusions"])
    cleaned["mutations"] = clean_mutations(tables["mutations"])
    cleaned["cellosaurus"] = clean_cellosaurus(tables["cellosaurus"])
    cleaned["depmap_profiles"] = clean_depmap_profiles(tables["depmap_profiles"])
    cleaned["sample_info"] = clean_sample_info(tables["sample_info"])
    cleaned["geo_info"] = clean_geo_info(tables["geo_info"])
    cleaned["hpa_desc"] = clean_hpa_desc(tables["hpa_desc"])
    cleaned["metabolomics"] = clean_metabolomics(tables["metabolomics"])
    cleaned["mirna"] = clean_mirna(tables["mirna"])
    cleaned["signatures"] = clean_signatures(tables["signatures"])
    return cleaned
