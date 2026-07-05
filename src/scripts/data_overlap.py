"""
data_overlap.py
Unique contributions from teammate's 02_data_cleaning_pipeline.ipynb.

Covers three areas not already in cleaning_functions.py / data_utils.py:
  1. clean_transpose_depmap_expr  – alternative depmap_expr cleaner that
                                    transposes the matrix (genes as rows).
  2. Table reporting              – per-column summary + Excel export.
  3. Master spine ID utilities    – build master ID sets and stamp every table
                                    with six spine columns:
                                      master_ensembl_gene_id, master_ach_id,
                                      master_cvcl_id, master_pr_id,
                                      master_gsm_id, master_cellline_name.
"""

import re
import pandas as pd


# ---------------------------------------------------------------------------
# 1. ALTERNATIVE DEPMAP EXPR CLEANER (transpose form)
# ---------------------------------------------------------------------------
# Our existing clean_depmap_expr() cleans in-place (rows = PR- samples,
# cols = gene headers).  This version transposes so that genes become rows
# and PR- sample IDs become columns, then splits each gene header into
# separate 'gene' (symbol) and 'ensg_id' columns.
# ---------------------------------------------------------------------------

def clean_transpose_depmap_expr(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean + transpose the DepMap expression table.

    Input : rows = PR- profile IDs (index), cols = 'TSPAN6 (ENSG00000000003)'
    Output: rows = genes, cols = ['gene', 'ensg_id', <PR- sample columns ...>]
    """
    df = df.copy()

    ens_core    = re.compile(r"ens[gtp]\d+",          re.IGNORECASE)
    ens_version = re.compile(r"(ens[gtp]\d+)\.\d+",   re.IGNORECASE)

    df.index = (
        df.index.astype(str)
        .str.strip()
        .str.replace(r"\s+", " ", regex=True)
        .str.lower()
    )
    df.index.name = "sample_id"

    dft = df.T
    dft.index.name = "gene_header"
    dft = dft.reset_index()

    def split_header(h):
        h = str(h).strip().lower()
        h = re.sub(r"\s+", " ", h)

        m    = ens_core.search(h)
        ensg = ens_version.sub(r"\1", m.group(0)) if m else pd.NA

        if "(" in h:
            sym = h.split("(")[0].strip()
            sym = sym if sym else pd.NA
        else:
            sym = pd.NA if (m and h == m.group(0)) else h

        return pd.Series({"gene": sym, "ensg_id": ensg})

    split_cols  = dft["gene_header"].apply(split_header)
    sample_cols = [c for c in dft.columns if c != "gene_header"]
    return pd.concat([split_cols, dft[sample_cols]], axis=1)


# ---------------------------------------------------------------------------
# 2. TABLE REPORTING
# ---------------------------------------------------------------------------

def build_table_summary(df: pd.DataFrame, n_samples: int = 3) -> pd.DataFrame:
    """
    Per-column summary: column | dtype | sample_values | missing_% |
                        unique_values | primary_key | duplicates.
    """
    rows  = []
    n_rows = len(df)

    for col in df.columns:
        s           = df[col]
        n_missing   = s.isna().sum()
        n_present   = n_rows - n_missing
        n_unique    = s.nunique(dropna=True)
        n_dupes     = n_present - n_unique
        missing_pct = round((n_missing / n_rows) * 100, 2) if n_rows else 0
        is_pk       = (n_missing == 0) and (n_dupes == 0) and (n_rows > 0)
        sample      = s.dropna().astype(str).head(n_samples).tolist()

        rows.append({
            "column":        col,
            "dtype":         str(s.dtype),
            "sample_values": ", ".join(sample),
            "missing_%":     missing_pct,
            "unique_values": n_unique,
            "primary_key":   is_pk,
            "duplicates":    n_dupes,
        })

    return pd.DataFrame(rows)


def export_clean_report(tables: dict, out_path: str, n_samples: int = 3) -> None:
    """Export one Excel sheet per table with shape header + column summary."""
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for name, df in tables.items():
            summary    = build_table_summary(df, n_samples=n_samples)
            sheet_name = name[:31]

            header = pd.DataFrame({
                "info": [
                    f"TABLE: {name}",
                    f"shape: {df.shape[0]} rows x {df.shape[1]} cols",
                ]
            })
            header.to_excel(writer, sheet_name=sheet_name,
                            index=False, header=False, startrow=0)
            summary.to_excel(writer, sheet_name=sheet_name,
                             index=False, startrow=3)

    print(f"Written {len(tables)} sheets to: {out_path}")


# ---------------------------------------------------------------------------
# 3. MASTER SPINE ID UTILITIES
# ---------------------------------------------------------------------------
# Each ID type has two functions:
#   build_master_<x>_list(tables)  -> (id_set, breakdown_df)
#   add_master_<x>_id(tables)      -> updated tables dict
# ---------------------------------------------------------------------------

# --- 3a. ENSG (gene) IDs ---------------------------------------------------

ENSG_ID_COLUMNS = {
    "hpa_rna_clean":     ["gene"],
    "geo_expr_clean":    ["gene"],
    "mutations_clean":   ["ensemblgeneid"],
    "fusions_clean":     ["gene1_ens_id", "gene2_ens_id"],
    "depmap_expr_clean": ["ensg_id"],
}

GENE_SOURCE = {
    "hpa_rna_clean":         "gene",
    "geo_expr_clean":        "gene",
    "mutations_clean":       "ensemblgeneid",
    "fusions_clean":         "gene1_ens_id",
    "depmap_expr_clean":     "ensg_id",
    "proteomics_clean":      None,
    "cellosaurus_clean":     None,
    "depmap_profiles_clean": None,
    "sample_info_clean":     None,
    "geo_info_clean":        None,
    "hpa_desc_clean":        None,
    "metabolomics_clean":    None,
    "mirna_clean":           None,
    "signatures_clean":      None,
}

_FUSIONS_GENE2_COL = "gene2_ens_id"


def _extract_ensg_set(series: pd.Series) -> set:
    found = (
        series.dropna().astype(str)
        .str.extract(r"(ens[gtp]\d+)", flags=re.IGNORECASE)[0]
        .dropna().str.lower()
    )
    return set(found.unique())


def _extract_ensg(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.extract(r"(ens[gtp]\d+)", flags=re.IGNORECASE)[0]
        .str.lower()
    )


def build_master_gene_list(tables: dict, verbose: bool = True):
    all_ensg  = set()
    breakdown = []
    for name, df in tables.items():
        table_ids = set()
        for col in ENSG_ID_COLUMNS.get(name, []):
            if col in df.columns:
                table_ids |= _extract_ensg_set(df[col])
            elif verbose:
                gene_like = [c for c in df.columns
                             if "ens" in str(c).lower() or "gene" in str(c).lower()][:5]
                print(f"  {name}: expected column '{col}' not found. "
                      f"gene-like cols: {gene_like}")
        all_ensg |= table_ids
        breakdown.append({"table": name, "n_ensg_found": len(table_ids)})
        if verbose:
            print(f"{name:24s}: {len(table_ids):>6,} distinct ENSG")
    breakdown_df = (pd.DataFrame(breakdown)
                    .sort_values("n_ensg_found", ascending=False)
                    .reset_index(drop=True))
    return all_ensg, breakdown_df


def add_master_gene_id(tables: dict) -> dict:
    updated = {}
    for name, df in tables.items():
        df         = df.copy()
        source_col = GENE_SOURCE.get(name)

        if source_col and source_col in df.columns:
            df["master_ensembl_gene_id"] = _extract_ensg(df[source_col])
            if name == "fusions_clean" and _FUSIONS_GENE2_COL in df.columns:
                g2 = _extract_ensg(df[_FUSIONS_GENE2_COL])
                df["master_ensembl_gene_id"] = df["master_ensembl_gene_id"].fillna(g2)
        else:
            df["master_ensembl_gene_id"] = pd.NA

        updated[name] = df
        n = df["master_ensembl_gene_id"].notna().sum()
        print(f"{name:24s}: master_ensembl_gene_id filled in {n:>9,} of {len(df):>9,} rows")
    return updated


# --- 3b. ACH- (DepMap model) IDs ------------------------------------------

ACH_ID_COLUMNS = {
    "sample_info_clean":     ["depmap_id"],
    "depmap_profiles_clean": ["modelid"],
    "fusions_clean":         ["modelid"],
    "signatures_clean":      ["modelid"],
    "metabolomics_clean":    ["depmap_id"],
    "proteomics_clean":      ["depmap_id"],
}

ACH_SOURCE = {
    "sample_info_clean":     "depmap_id",
    "depmap_profiles_clean": "modelid",
    "fusions_clean":         "modelid",
    "signatures_clean":      "modelid",
    "metabolomics_clean":    "depmap_id",
    "proteomics_clean":      "depmap_id",
    "depmap_expr_clean":     None,
    "mutations_clean":       None,
    "geo_expr_clean":        None,
    "geo_info_clean":        None,
    "hpa_rna_clean":         None,
    "hpa_desc_clean":        None,
    "mirna_clean":           None,
    "cellosaurus_clean":     None,
}


def _extract_ach_set(series: pd.Series) -> set:
    found = (
        series.dropna().astype(str)
        .str.extract(r"(ach-\d+)", flags=re.IGNORECASE)[0]
        .dropna().str.lower()
    )
    return set(found.unique())


def _extract_ach(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.extract(r"(ach-\d+)", flags=re.IGNORECASE)[0]
        .str.lower()
    )


def build_master_cellline_list(tables: dict, verbose: bool = True):
    all_ach   = set()
    breakdown = []
    for name, df in tables.items():
        table_ids = set()
        for col in ACH_ID_COLUMNS.get(name, []):
            if col in df.columns:
                table_ids |= _extract_ach_set(df[col])
        all_ach |= table_ids
        breakdown.append({"table": name, "n_ach_found": len(table_ids)})
        if verbose:
            print(f"{name:24s}: {len(table_ids):>6,} distinct ACH-")
    breakdown_df = (pd.DataFrame(breakdown)
                    .sort_values("n_ach_found", ascending=False)
                    .reset_index(drop=True))
    return all_ach, breakdown_df


def add_master_ach_id(tables: dict) -> dict:
    updated = {}
    for name, df in tables.items():
        df         = df.copy()
        source_col = ACH_SOURCE.get(name)

        if source_col and source_col in df.columns:
            df["master_ach_id"] = _extract_ach(df[source_col])
        else:
            df["master_ach_id"] = pd.NA

        updated[name] = df
        n = df["master_ach_id"].notna().sum()
        print(f"{name:24s}: master_ach_id filled in {n:>9,} of {len(df):>9,} rows")
    return updated


# --- 3c. CVCL- (Cellosaurus) IDs ------------------------------------------

CVCL_ID_COLUMNS = {
    "cellosaurus_clean": ["cellosaurus_accession"],
    "sample_info_clean": ["rrid"],
    "hpa_desc_clean":    ["cellosaurus id"],
    "geo_info_clean":    ["cellosaurus_id"],
}

CVCL_SOURCE = {
    "cellosaurus_clean":     "cellosaurus_accession",
    "sample_info_clean":     "rrid",
    "hpa_desc_clean":        "cellosaurus id",
    "geo_info_clean":        "cellosaurus_id",
    "depmap_profiles_clean": None,
    "fusions_clean":         None,
    "signatures_clean":      None,
    "metabolomics_clean":    None,
    "proteomics_clean":      None,
    "depmap_expr_clean":     None,
    "mutations_clean":       None,
    "geo_expr_clean":        None,
    "hpa_rna_clean":         None,
    "mirna_clean":           None,
}


def _extract_cvcl_set(series: pd.Series) -> set:
    found = (
        series.dropna().astype(str)
        .str.extract(r"(cvcl_[a-z0-9]+)", flags=re.IGNORECASE)[0]
        .dropna().str.lower()
    )
    return set(found.unique())


def _extract_cvcl(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.extract(r"(cvcl_[a-z0-9]+)", flags=re.IGNORECASE)[0]
        .str.lower()
    )


def build_master_cvcl_list(tables: dict, verbose: bool = True):
    all_cvcl  = set()
    breakdown = []
    for name, df in tables.items():
        table_ids = set()
        for col in CVCL_ID_COLUMNS.get(name, []):
            if col in df.columns:
                table_ids |= _extract_cvcl_set(df[col])
            elif verbose:
                cands = [c for c in df.columns
                         if "cvcl" in str(c).lower() or "cellosaurus" in str(c).lower()]
                print(f"  {name}: '{col}' not found. candidates: {cands}")
        all_cvcl |= table_ids
        breakdown.append({"table": name, "n_cvcl_found": len(table_ids)})
        if verbose:
            print(f"{name:24s}: {len(table_ids):>6,} distinct CVCL-")
    breakdown_df = (pd.DataFrame(breakdown)
                    .sort_values("n_cvcl_found", ascending=False)
                    .reset_index(drop=True))
    return all_cvcl, breakdown_df


def add_master_cvcl_id(tables: dict) -> dict:
    updated = {}
    for name, df in tables.items():
        df         = df.copy()
        source_col = CVCL_SOURCE.get(name)

        if source_col and source_col in df.columns:
            df["master_cvcl_id"] = _extract_cvcl(df[source_col])
        else:
            df["master_cvcl_id"] = pd.NA

        updated[name] = df
        n = df["master_cvcl_id"].notna().sum()
        print(f"{name:24s}: master_cvcl_id filled in {n:>9,} of {len(df):>9,} rows")
    return updated


# --- 3d. PR- (sequencing profile) IDs ------------------------------------

PR_ID_COLUMNS = {
    "depmap_profiles_clean": ["profileid"],
    "mutations_clean":       ["profileid"],
    "signatures_clean":      ["sequencingid"],
    "depmap_expr_clean":     ["sample_id"],
}

PR_SOURCE = {
    "depmap_profiles_clean": "profileid",
    "mutations_clean":       "profileid",
    "signatures_clean":      "sequencingid",
    "depmap_expr_clean":     "sample_id",
    "sample_info_clean":     None,
    "fusions_clean":         None,
    "metabolomics_clean":    None,
    "proteomics_clean":      None,
    "cellosaurus_clean":     None,
    "geo_expr_clean":        None,
    "geo_info_clean":        None,
    "hpa_rna_clean":         None,
    "hpa_desc_clean":        None,
    "mirna_clean":           None,
}


def _extract_pr_set(series: pd.Series) -> set:
    found = (
        series.dropna().astype(str)
        .str.extract(r"(pr-[a-z0-9]+)", flags=re.IGNORECASE)[0]
        .dropna().str.lower()
    )
    return set(found.unique())


def _extract_pr(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.extract(r"(pr-[a-z0-9]+)", flags=re.IGNORECASE)[0]
        .str.lower()
    )


def build_master_profile_list(tables: dict, verbose: bool = True):
    all_pr    = set()
    breakdown = []
    for name, df in tables.items():
        table_ids = set()
        for col in PR_ID_COLUMNS.get(name, []):
            if col in df.columns:
                table_ids |= _extract_pr_set(df[col])
            elif verbose:
                cands = [c for c in df.columns
                         if "pr" in str(c).lower() or "profile" in str(c).lower()
                         or "sequencing" in str(c).lower() or "sample" in str(c).lower()]
                print(f"  {name}: '{col}' not found. candidates: {cands}")
        all_pr |= table_ids
        breakdown.append({"table": name, "n_pr_found": len(table_ids)})
        if verbose:
            print(f"{name:24s}: {len(table_ids):>6,} distinct PR-")
    breakdown_df = (pd.DataFrame(breakdown)
                    .sort_values("n_pr_found", ascending=False)
                    .reset_index(drop=True))
    return all_pr, breakdown_df


def add_master_pr_id(tables: dict) -> dict:
    updated = {}
    for name, df in tables.items():
        df         = df.copy()
        source_col = PR_SOURCE.get(name)

        if source_col and source_col in df.columns:
            df["master_pr_id"] = _extract_pr(df[source_col])
        else:
            df["master_pr_id"] = pd.NA

        updated[name] = df
        n = df["master_pr_id"].notna().sum()
        print(f"{name:24s}: master_pr_id filled in {n:>9,} of {len(df):>9,} rows")
    return updated


# --- 3e. GSM (GEO sample) IDs ---------------------------------------------

GSM_ID_SOURCE = {
    "geo_info_clean": ("value",  "geo_accession"),
    "geo_expr_clean": ("header", None),
}

GSM_SOURCE = {
    "geo_info_clean":        "geo_accession",
    "geo_expr_clean":        None,   # GSM in column headers, not rows
    "sample_info_clean":     None,
    "depmap_profiles_clean": None,
    "fusions_clean":         None,
    "signatures_clean":      None,
    "metabolomics_clean":    None,
    "proteomics_clean":      None,
    "cellosaurus_clean":     None,
    "depmap_expr_clean":     None,
    "mutations_clean":       None,
    "hpa_rna_clean":         None,
    "hpa_desc_clean":        None,
    "mirna_clean":           None,
}

_GSM_COL_PAT = re.compile(r"gsm\d+", re.IGNORECASE)


def _extract_gsm_set_from_values(series: pd.Series) -> set:
    found = (
        series.dropna().astype(str)
        .str.extract(r"(gsm\d+)", flags=re.IGNORECASE)[0]
        .dropna().str.lower()
    )
    return set(found.unique())


def _extract_gsm_set_from_headers(columns) -> set:
    ids = set()
    for c in columns:
        m = _GSM_COL_PAT.search(str(c))
        if m:
            ids.add(m.group(0).lower())
    return ids


def _extract_gsm(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.extract(r"(gsm\d+)", flags=re.IGNORECASE)[0]
        .str.lower()
    )


def build_master_gsm_list(tables: dict, verbose: bool = True):
    all_gsm   = set()
    breakdown = []
    for name, df in tables.items():
        table_ids = set()
        if name in GSM_ID_SOURCE:
            loc_type, col = GSM_ID_SOURCE[name]
            if loc_type == "header":
                table_ids |= _extract_gsm_set_from_headers(df.columns)
            elif col and col in df.columns:
                table_ids |= _extract_gsm_set_from_values(df[col])
            elif verbose and col:
                cands = [c for c in df.columns
                         if "gsm" in str(c).lower() or "geo" in str(c).lower()
                         or "accession" in str(c).lower()]
                print(f"  {name}: '{col}' not found. candidates: {cands}")
        all_gsm |= table_ids
        breakdown.append({"table": name, "n_gsm_found": len(table_ids)})
        if verbose:
            print(f"{name:24s}: {len(table_ids):>6,} distinct GSM")
    breakdown_df = (pd.DataFrame(breakdown)
                    .sort_values("n_gsm_found", ascending=False)
                    .reset_index(drop=True))
    return all_gsm, breakdown_df


def add_master_gsm_id(tables: dict) -> dict:
    updated = {}
    for name, df in tables.items():
        df         = df.copy()
        source_col = GSM_SOURCE.get(name)

        if source_col and source_col in df.columns:
            df["master_gsm_id"] = _extract_gsm(df[source_col])
        else:
            df["master_gsm_id"] = pd.NA

        updated[name] = df
        n = df["master_gsm_id"].notna().sum()
        print(f"{name:24s}: master_gsm_id filled in {n:>9,} of {len(df):>9,} rows")
    return updated


# --- 3f. Cell line NAME spine ---------------------------------------------

NAME_ID_COLUMNS = {
    "sample_info_clean":  ["cell_line_name"],
    "hpa_rna_clean":      ["cell line"],
    "hpa_desc_clean":     ["cell line"],
    "cellosaurus_clean":  ["cellosaurus_cell_line_name"],
    "metabolomics_clean": ["ccle_id"],
    "geo_info_clean":     ["cellline"],
}

NAME_SOURCE = {
    "sample_info_clean":     "cell_line_name",
    "hpa_rna_clean":         "cell line",
    "hpa_desc_clean":        "cell line",
    "cellosaurus_clean":     "cellosaurus_cell_line_name",
    "metabolomics_clean":    "cell_line_name",
    "geo_info_clean":        "cellline",
    "depmap_profiles_clean": None,
    "fusions_clean":         None,
    "signatures_clean":      None,
    "proteomics_clean":      None,
    "depmap_expr_clean":     None,
    "mutations_clean":       None,
    "geo_expr_clean":        None,
    "mirna_clean":           None,
}


def _clean_name_series(series: pd.Series) -> pd.Series:
    s = (series.astype(str)
         .str.strip()
         .str.replace(r"\s+", " ", regex=True)
         .str.lower())
    return s.replace("nan", pd.NA).replace("<na>", pd.NA).replace("", pd.NA)


def build_master_cellline_name_list(tables: dict, verbose: bool = True):
    all_names = set()
    breakdown = []
    for name, df in tables.items():
        table_names = set()
        for col in NAME_ID_COLUMNS.get(name, []):
            if col in df.columns:
                table_names |= set(_clean_name_series(df[col]).dropna().unique())
            elif verbose:
                cands = [c for c in df.columns
                         if "name" in str(c).lower() or "cell" in str(c).lower()
                         or "ccle" in str(c).lower()]
                print(f"  {name}: '{col}' not found. candidates: {cands}")

        # mirna: names live in column headers, not rows
        if name == "mirna_clean":
            header_names = {str(c).strip().lower() for c in df.columns} - {"name", "description"}
            table_names |= header_names

        all_names |= table_names
        breakdown.append({"table": name, "n_names_found": len(table_names)})
        if verbose:
            print(f"{name:24s}: {len(table_names):>6,} distinct names")

    breakdown_df = (pd.DataFrame(breakdown)
                    .sort_values("n_names_found", ascending=False)
                    .reset_index(drop=True))
    return all_names, breakdown_df


def add_master_cellline_name(tables: dict) -> dict:
    updated = {}
    for name, df in tables.items():
        df         = df.copy()
        source_col = NAME_SOURCE.get(name)

        if source_col and source_col in df.columns:
            df["master_cellline_name"] = _clean_name_series(df[source_col])
        else:
            df["master_cellline_name"] = pd.NA

        updated[name] = df
        n = df["master_cellline_name"].notna().sum()
        print(f"{name:24s}: master_cellline_name filled in {n:>9,} of {len(df):>9,} rows")
    return updated


# ---------------------------------------------------------------------------
# 4. CONVENIENCE: stamp all six spine columns in one call
# ---------------------------------------------------------------------------

def add_all_spine_ids(tables: dict) -> dict:
    """
    Chain all six add_master_* functions onto a tables dict.
    Returns updated dict with columns:
        master_ensembl_gene_id, master_ach_id, master_cvcl_id,
        master_pr_id, master_gsm_id, master_cellline_name
    """
    tables = add_master_gene_id(tables)
    tables = add_master_ach_id(tables)
    tables = add_master_cvcl_id(tables)
    tables = add_master_pr_id(tables)
    tables = add_master_gsm_id(tables)
    tables = add_master_cellline_name(tables)
    return tables


def print_spine_summary(tables: dict) -> None:
    """Print a compact summary of how many rows each spine column covers."""
    spine_cols = [
        ("gene",  "master_ensembl_gene_id"),
        ("ach",   "master_ach_id"),
        ("cvcl",  "master_cvcl_id"),
        ("pr",    "master_pr_id"),
        ("gsm",   "master_gsm_id"),
        ("name",  "master_cellline_name"),
    ]
    header = f"{'table':<24}" + "".join(f"{label:>10}" for label, _ in spine_cols)
    print("=" * (24 + 10 * len(spine_cols)))
    print(header)
    print("-" * (24 + 10 * len(spine_cols)))
    for tname, df in tables.items():
        vals = [
            df[col].notna().sum() if col in df.columns else 0
            for _, col in spine_cols
        ]
        print(f"{tname:<24}" + "".join(f"{v:>10,}" for v in vals))
    all_present = all(
        all(col in df.columns for _, col in spine_cols)
        for df in tables.values()
    )
    print(f"\nAll tables have all six spine columns: {all_present}")
