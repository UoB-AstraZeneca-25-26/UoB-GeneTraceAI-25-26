from __future__ import annotations
from typing import Any, Dict, List, Optional, Tuple
import pandas as pd


# ------------------------------------------------------------------ #
# Stage A: collapse aux table → coverage flags + optional carry cols  #
# ------------------------------------------------------------------ #
def _collapse_aux_table(
    df: pd.DataFrame,
    cfg: Dict[str, Any],
    dup_log_rows: List[dict],          
) -> pd.DataFrame:
    """
    Turn a multi-row-per-entity aux table into one row per entity.
    - Boolean coverage flags: has_<value>  (e.g. has_rna, has_wes, has_wgs)
    - Carry cols: pipe-joined per entity+collapse_on
      (e.g. modelcondition_rna, weskit_wes)
    """
    merge_key  = cfg["merge_key"]
    collapse_on = cfg.get("collapse_on")
    carry_cols  = cfg.get("carry_cols", [])          

    if collapse_on:
        # --- BUG 1 FIX: log duplicate merge_key + collapse_on pairs before pivot ---
        dup_mask = df.duplicated(subset=[merge_key, collapse_on], keep=False)
        n_dups = dup_mask.sum()
        if n_dups > 0:
            print(f"  [WARN] {n_dups} duplicate ({merge_key}, {collapse_on}) pairs "
                  f"found in aux table — pipe-joining string cols, any() for flags")
            for _, row in df.loc[dup_mask].iterrows():
                dup_log_rows.append({
                    "merge_key_value": row[merge_key],
                    "collapse_on_value": row[collapse_on],
                    "reason": f"duplicate {merge_key}+{collapse_on} pair in aux table",
                })

        # --- Boolean coverage flags ---
        flags = (
            df.assign(_present=True)
            .pivot_table(
                index=merge_key,
                columns=collapse_on,
                values="_present",
                aggfunc="any",
                fill_value=False,
            )
            .reset_index()
        )
        flags.columns = [
            f"has_{str(c).lower()}" if c != merge_key else c
            for c in flags.columns
        ]

        # --- BUG 2 FIX: carry cols — pivot with pipe-join for string values ---
        for col in carry_cols:
            if col not in df.columns:
                print(f"  [WARN] carry_col '{col}' not found in aux df — skipping")
                continue

            carried = (
                df.groupby([merge_key, collapse_on])[col]
                .apply(lambda x: "|".join(x.dropna().unique()))   
                .unstack(fill_value=pd.NA)
                .reset_index()
            )
            carried.columns = [
                f"{col}_{str(c).lower()}" if c != merge_key else c
                for c in carried.columns
            ]
            flags = flags.merge(carried, on=merge_key, how="left")

    else:
        flags = df[[merge_key]].drop_duplicates().reset_index(drop=True)
        flags["has_record"] = True

    assert flags[merge_key].duplicated().sum() == 0, (
        f"collapse step did not fully de-duplicate on '{merge_key}'"
    )
    return flags


# ------------------------------------------------------------------ #
# Stage B: cascade-match (unchanged — was correct)                    #
# ------------------------------------------------------------------ #
def _cascade_match(
    out: pd.DataFrame,
    reference_df: pd.DataFrame,
    ref_cfg: Dict[str, Any],
) -> pd.DataFrame:
    if "matched_id" not in out.columns:
        out["matched_id"] = pd.NA
        out["match_method"] = pd.NA

    id_col    = ref_cfg["id_col"]
    ref_index = reference_df.set_index(id_col)

    # 1. direct ID match
    direct_col = ref_cfg.get("direct_match_col")
    if direct_col and direct_col in out.columns:
        unmatched  = out["matched_id"].isna()
        is_direct  = unmatched & out[direct_col].isin(ref_index.index)
        out.loc[is_direct, "matched_id"]   = out.loc[is_direct, direct_col]
        out.loc[is_direct, "match_method"] = "direct_id"

    # 2. exact name match
    name_col     = ref_cfg.get("name_col")
    ref_name_col = ref_cfg.get("ref_name_col", name_col)
    if name_col and name_col in out.columns and ref_name_col:
        unmatched   = out["matched_id"].isna()
        name_lookup = reference_df.set_index(ref_name_col)[id_col]
        hit         = out.loc[unmatched, name_col].map(name_lookup)
        idx         = hit.dropna().index
        out.loc[idx, "matched_id"]   = hit.dropna()
        out.loc[idx, "match_method"] = "name_match"

    # 3. synonym fallback
    syn_col      = ref_cfg.get("synonym_col")
    alt_name_col = ref_cfg.get("alt_name_col", name_col)
    if syn_col and alt_name_col and alt_name_col in out.columns:
        unmatched = out["matched_id"].isna()
        delimiter = ref_cfg.get("synonym_delimiter", ";")
        syn_df    = reference_df[[id_col, syn_col]].dropna(subset=[syn_col])
        syn_df    = syn_df.assign(
            **{syn_col: syn_df[syn_col].str.split(delimiter)}
        ).explode(syn_col)
        syn_df[syn_col] = syn_df[syn_col].str.strip()
        syn_lookup = syn_df.set_index(syn_col)[id_col]
        hit  = out.loc[unmatched, alt_name_col].map(syn_lookup)
        idx  = hit.dropna().index
        out.loc[idx, "matched_id"]   = hit.dropna()
        out.loc[idx, "match_method"] = "synonym_fallback"

    return out


# ------------------------------------------------------------------ #
# Main entry point                                                     #
# ------------------------------------------------------------------ #
def build_lookup_table(
    dataframes: Dict[str, pd.DataFrame],
    identifiers: Dict[str, Dict[str, Any]],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Returns
    -------
    lookup_table   : one row per base entity, canonical IDs + flags
    unmapped_log   : rows with no reference match
    duplicate_log  : aux-table duplicate pairs 
    """
    base_name = next(n for n, c in identifiers.items() if c["role"] == "base")
    base_cfg  = identifiers[base_name]
    base_key  = base_cfg["key"]
    result    = dataframes[base_name].copy()

    dup_log_rows: List[dict] = []   
    flag_cols:    List[str]  = []

    # --- merge aux tables ---
    for name, cfg in identifiers.items():
        if cfg["role"] != "aux":
            continue
        # BUG 4 FIX: pass dup_log_rows into collapse
        collapsed    = _collapse_aux_table(dataframes[name], cfg, dup_log_rows)
        new_flag_cols = [c for c in collapsed.columns if c != cfg["merge_key"]]
        flag_cols.extend(new_flag_cols)

        result = result.merge(
            collapsed,
            left_on=cfg["base_key"],
            right_on=cfg["merge_key"],
            how="left",
        )
        for col in [c for c in new_flag_cols if result[col].dtype == bool]:
            result[col] = result[col].fillna(False)

    # --- cascade-match reference tables ---
    reference_names = [n for n, c in identifiers.items() if c["role"] == "reference"]
    for ref_name in reference_names:
        result = _cascade_match(result, dataframes[ref_name], identifiers[ref_name])

    # --- unmapped log ---
    if reference_names:
        unmapped_mask = result["matched_id"].isna()
        unmapped_log  = pd.DataFrame({
            "original_value":         result.loc[unmapped_mask, base_key],
            "dataset":                base_name,
            "resolution_step_reached":"synonym_fallback",
            "reason":                 "no direct_id / name / synonym match in any reference table",
        })
    else:
        unmapped_log = pd.DataFrame(
            columns=["original_value","dataset","resolution_step_reached","reason"]
        )

    duplicate_log = pd.DataFrame(dup_log_rows) if dup_log_rows else pd.DataFrame(
        columns=["merge_key_value","collapse_on_value","reason"]
    )

    if "matched_id" in result.columns:
        result = result.rename(columns={"matched_id": "cvcl_accession"})

    # --- assemble final lookup ---
    keep_cols = [base_key]
    for optional_col in ("name_col", "alt_name_col"):
        if optional_col in base_cfg:
            keep_cols.append(base_cfg[optional_col])
    keep_cols += flag_cols
    if reference_names:
        keep_cols += ["cvcl_accession", "match_method"]   # BUG 3 FIX: updated name
    keep_cols = [c for c in keep_cols if c in result.columns]

    lookup_table = result[keep_cols].copy()

    # --- summary ---
    print(f"Total entities       : {len(lookup_table)}")
    print(f"Aux duplicate pairs  : {len(duplicate_log)}")
    if reference_names:
        for method in ("direct_id", "name_match", "synonym_fallback"):
            print(f"Resolved via {method:20s}: {(lookup_table['match_method'] == method).sum()}")
        print(f"Unmapped             : {lookup_table['cvcl_accession'].isna().sum()}")

    return lookup_table, unmapped_log, duplicate_log

# ---------------------------------------------------------
# --- Sample Info ---

# track_a_tables = load_selected_parquets(
#     ["sample_info", "depmap_profiles", "cellosaurus"]
# )

# identifiers = {
#     "sample_info": {
#         "role":         "base",
#         "key":          "DepMap_ID",
#         "name_col":     "cell_line_name",
#         "alt_name_col": "stripped_cell_line_name",
#     },
#     "depmap_profiles": {
#         "role":        "aux",
#         "merge_key":   "modelid",
#         "base_key":    "DepMap_ID",
#         "collapse_on": "datatype",
#         "carry_cols":  ["modelcondition", "weskit"],  
#     },
#     "cellosaurus": {
#         "role":              "reference",
#         "id_col":            "Accession (CVCL_xxxx)",
#         "direct_match_col":  "RRID",
#         "name_col":          "cell_line_name",
#         "ref_name_col":      "Identifier (cell line name)",
#         "synonym_col":       "Synonyms",
#         "alt_name_col":      "stripped_cell_line_name",
#         "synonym_delimiter": ";",
#     },
# }

# lookup_table, unmapped_log, duplicate_log = build_lookup_table(
#     track_a_tables, identifiers
# )

# unmapped_log.to_csv("/data/logs/unmapped_celllines.csv", index=False)
# duplicate_log.to_csv("/data/logs/duplicate_pairs_depmap_profiles.csv", index=False)