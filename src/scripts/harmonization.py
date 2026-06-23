from __future__ import annotations
from typing import Any, Dict, List, Optional, Tuple
import pandas as pd


# ------------------------------------------------------------------ #
# Helper: parse delimiter-separated cross-reference columns           #
# ------------------------------------------------------------------ #
def _parse_xref_column(
    df: pd.DataFrame,
    id_col: str,
    xref_col: str,
    prefix: str,
    outer_delimiter: str = "||",
    inner_delimiter: str = ";",
) -> pd.Series:
    """
    Parse a delimited cross-reference column to extract IDs with a given prefix.

    Cellosaurus Cross-references format:
      "Wikidata; Q108819335 || DepMap; ACH-000001 || Sanger; SIDM00001"

    Args:
        df               : reference DataFrame containing id_col + xref_col
        id_col           : primary key column (e.g. "Accession (CVCL_xxxx)")
        xref_col         : cross-reference column (e.g. "Cross-references")
        prefix           : filter prefix  (e.g. "DepMap", "Sanger")
        outer_delimiter  : separates each db entry  (default " || ")
        inner_delimiter  : separates prefix from value  (default ";")

    Returns:
        pd.Series indexed by id_col → extracted xref value
    """
    xref_df = df[[id_col, xref_col]].dropna(subset=[xref_col]).copy()
    xref_df[xref_col] = xref_df[xref_col].str.split(outer_delimiter)
    xref_df = xref_df.explode(xref_col)
    xref_df[xref_col] = xref_df[xref_col].str.strip()
    xref_df = xref_df[xref_df[xref_col].str.startswith(prefix, na=False)].copy()
    xref_df["_extracted"] = (
        xref_df[xref_col].str.split(inner_delimiter).str[-1].str.strip()
    )
    # drop duplicates: keep first xref value per id
    xref_df = xref_df.drop_duplicates(subset=[id_col])
    return xref_df.set_index(id_col)["_extracted"]


# ------------------------------------------------------------------ #
# Stage A: collapse aux table → coverage flags + carry cols           #
# ------------------------------------------------------------------ #
def _collapse_aux_table(
    df: pd.DataFrame,
    cfg: Dict[str, Any],
    dup_log_rows: List[dict],
) -> pd.DataFrame:
    """
    Collapse a multi-row-per-entity aux table into one row per entity.

    Produces:
    - Boolean coverage flags:  has_<datatype>  (e.g. has_rna, has_wes, has_wgs)
    - Carry cols pipe-joined:  <col>_<datatype> (e.g. modelcondition_rna, weskit_wes)

    Config keys:
        merge_key   : column that uniquely identifies the entity in this df
        collapse_on : column whose values become flag suffixes  (e.g. "datatype")
        carry_cols  : list of string cols to pipe-join per entity+collapse_on
    """
    merge_key   = cfg["merge_key"]
    collapse_on = cfg.get("collapse_on")
    carry_cols  = cfg.get("carry_cols", [])

    if collapse_on:
        # -- Log duplicate (merge_key, collapse_on) pairs per integration contract --
        dup_mask = df.duplicated(subset=[merge_key, collapse_on], keep=False)
        n_dups   = dup_mask.sum()
        if n_dups > 0:
            print(f"  [WARN] {n_dups} duplicate ({merge_key}, {collapse_on}) pairs "
                  f"— pipe-joining string cols, any() for flags")
            for _, row in df.loc[dup_mask].iterrows():
                dup_log_rows.append({
                    "merge_key_value":    row[merge_key],
                    "collapse_on_value":  row[collapse_on],
                    "reason": f"duplicate {merge_key}+{collapse_on} pair in aux table",
                })

        # -- Boolean coverage flags --
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

        # -- Carry cols: pipe-join per (merge_key, collapse_on) --
        for col in carry_cols:
            if col not in df.columns:
                print(f"  [WARN] carry_col '{col}' not in aux df — skipping")
                continue
            carried = (
                df.groupby([merge_key, collapse_on])[col]
                .apply(lambda x: "|".join(x.dropna().astype(str).unique()))
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
        f"_collapse_aux_table: '{merge_key}' still has duplicates after pivot"
    )
    return flags


# ------------------------------------------------------------------ #
# Stage A2: recover unmatched aux rows via cellosaurus xref bridge    #
# ------------------------------------------------------------------ #
def _recover_aux_via_xref_bridge(
    result: pd.DataFrame,
    aux_collapsed: pd.DataFrame,
    cfg: Dict[str, Any],
    reference_dfs: Dict[str, pd.DataFrame],
    flag_cols: List[str],
) -> pd.DataFrame:
    """
    Recovery step for aux table rows that failed the direct key join.

    Concrete case (from insights):
      sample_info.DepMap_ID ↔ depmap_profiles.modelid failed for 91 rows.
      Recovery: sample_info.RRID → cellosaurus.Accession
                → cellosaurus["Cross-references"] (DepMap; ACH-XXXXXX)
                → match recovered ACH to depmap_profiles.modelid  (+19)

    Config key "xref_bridge" (optional):
        reference      : name of the reference df to use (e.g. "cellosaurus")
        src_col        : column in base table  (e.g. "RRID")
        ref_id_col     : canonical key in reference df  (e.g. "Accession (CVCL_xxxx)")
        ref_xref_col   : cross-ref column in reference df  (e.g. "Cross-references")
        xref_prefix    : prefix to filter cross-refs  (e.g. "DepMap")
    """
    bridge_cfg = cfg.get("xref_bridge")
    if not bridge_cfg:
        return result

    ref_name    = bridge_cfg["reference"]
    src_col     = bridge_cfg["src_col"]
    ref_id_col  = bridge_cfg["ref_id_col"]
    ref_xref_col = bridge_cfg["ref_xref_col"]
    xref_prefix = bridge_cfg["xref_prefix"]
    base_key    = cfg["base_key"]
    merge_key   = cfg["merge_key"]

    if ref_name not in reference_dfs:
        print(f"  [WARN] xref_bridge reference '{ref_name}' not in dataframes — skipping")
        return result

    ref_df = reference_dfs[ref_name]

    # Build: CVCL → ACH (from cellosaurus cross-references)
    cvcl_to_ach = _parse_xref_column(
        ref_df, ref_id_col, ref_xref_col, prefix=xref_prefix
    )

    # Identify unmatched rows (flag cols are all False/NaN from the left join)
    bool_flags    = [c for c in flag_cols if str(result[c].dtype) == "bool"]
    unmatched_mask = result[bool_flags].eq(False).all(axis=1) if bool_flags else result[flag_cols[0]].isna()

    unmatched = result.loc[unmatched_mask].copy()
    if unmatched.empty or src_col not in unmatched.columns:
        return result

    # Map: src_col (RRID) → CVCL → ACH
    recovered_ach = unmatched[src_col].map(cvcl_to_ach)

    # Merge recovered ACH IDs against the collapsed aux table
    recover_df = pd.DataFrame({
        base_key:   unmatched[base_key].values,
        merge_key:  recovered_ach.values,
    }).dropna(subset=[merge_key])

    if recover_df.empty:
        return result

    recovered_flags = recover_df.merge(
        aux_collapsed, on=merge_key, how="inner"
    ).drop(columns=[merge_key])

    n_recovered = len(recovered_flags)
    print(f"  [BRIDGE] xref_bridge ({xref_prefix}) recovered {n_recovered} rows "
          f"that failed direct key join")

    # Write recovered values back into result
    result = result.set_index(base_key)
    recovered_flags = recovered_flags.set_index(base_key)

    for col in flag_cols:
        if col in recovered_flags.columns:
            result.loc[recovered_flags.index, col] = recovered_flags[col]

    result["bridge_method"] = result.get("bridge_method", pd.NA)
    result.loc[recovered_flags.index, "bridge_method"] = f"xref_{xref_prefix.lower()}"

    return result.reset_index()


# ------------------------------------------------------------------ #
# Stage B: cascade-match base against reference tables                #
# ------------------------------------------------------------------ #
def _cascade_match(
    out: pd.DataFrame,
    reference_df: pd.DataFrame,
    ref_cfg: Dict[str, Any],
) -> pd.DataFrame:
    """
    Fill matched_id / match_method for unmatched rows, in priority order:

    1. direct_id           — base.direct_match_col ↔ reference.id_col
    2. secondary_accession — base.direct_match_col ↔ exploded secondary accessions
    3. name_match          — base.name_col ↔ reference.ref_name_col
    4. synonym_fallback    — base.alt_name_col ↔ exploded reference.synonym_col
    5. sanger_xref         — base.src_xref_col ↔ reference.xref_col (prefix filter)

    From insights:
      Chain 1+2 → 1,812 + 1 = 1,813  (RRID direct + secondary accession)
      Chain 3   → +8                   (name match)
      Chain 5   → +17                  (Sanger xref via cross-references)
      Total     → 1,838 / 1,840 = 99.89%

    NOTE: All .map() calls use plain dicts (not pd.Series) to avoid
    InvalidIndexError on PyArrow-backed string columns (ArrowExtensionArray
    routes through get_indexer which enforces unique-index strictly).
    Ambiguous keys (duplicated after explode) are always dropped with
    keep=False — a silent wrong assignment is worse than a missed match.
    """
    if "matched_id" not in out.columns:
        out = out.copy()
        out["matched_id"]   = pd.NA
        out["match_method"] = pd.NA

    id_col    = ref_cfg["id_col"]
    ref_index = reference_df.set_index(id_col)

    # ---- Step 1: direct ID match ---------------------------------- #
    direct_col = ref_cfg.get("direct_match_col")
    if direct_col and direct_col in out.columns:
        unmatched = out["matched_id"].isna()
        is_direct = unmatched & out[direct_col].isin(ref_index.index)
        out.loc[is_direct, "matched_id"]   = out.loc[is_direct, direct_col]
        out.loc[is_direct, "match_method"] = "direct_id"
        print(f"  Step 1 direct_id          : +{is_direct.sum()}")

    # ---- Step 2: secondary accession match ------------------------ #
    sec_acc_col = ref_cfg.get("secondary_accession_col")
    if sec_acc_col and direct_col and sec_acc_col in reference_df.columns:
        unmatched = out["matched_id"].isna()
        if unmatched.any():
            sec_delimiter = ref_cfg.get("secondary_accession_delimiter", ";")
            sec_df = reference_df[[id_col, sec_acc_col]].dropna(subset=[sec_acc_col])
            sec_df = sec_df.assign(
                **{sec_acc_col: sec_df[sec_acc_col].str.split(sec_delimiter)}
            ).explode(sec_acc_col)
            sec_df[sec_acc_col] = sec_df[sec_acc_col].str.strip()

            # Drop ambiguous: same secondary accession → >1 primary CVCL
            ambiguous = sec_df.duplicated(subset=[sec_acc_col], keep=False)
            if ambiguous.sum():
                print(f"  [WARN] Step 2: dropping "
                      f"{sec_df.loc[ambiguous, sec_acc_col].nunique()} ambiguous "
                      f"secondary accession(s) ({ambiguous.sum()} rows)")
            sec_df = sec_df[~ambiguous]

            # Use dict — bypasses ArrowExtensionArray.map → get_indexer path
            sec_lookup = dict(zip(sec_df[sec_acc_col], sec_df[id_col]))
            hit = out.loc[unmatched, direct_col].map(sec_lookup)
            idx = hit.dropna().index
            out.loc[idx, "matched_id"]   = hit.dropna()
            out.loc[idx, "match_method"] = "secondary_accession"
            print(f"  Step 2 secondary_accession: +{len(idx)}")

    # ---- Step 3: exact name match --------------------------------- #
    name_col     = ref_cfg.get("name_col")
    ref_name_col = ref_cfg.get("ref_name_col", name_col)
    if name_col and name_col in out.columns and ref_name_col:
        unmatched = out["matched_id"].isna()
        if unmatched.any():
            name_counts = reference_df[ref_name_col].value_counts()
            safe_names  = name_counts[name_counts == 1].index
            safe_ref    = reference_df.loc[
                reference_df[ref_name_col].isin(safe_names)
            ]
            # Use dict — bypasses ArrowExtensionArray.map → get_indexer path
            safe_lookup = dict(zip(safe_ref[ref_name_col], safe_ref[id_col]))
            hit = out.loc[unmatched, name_col].map(safe_lookup)
            idx = hit.dropna().index
            out.loc[idx, "matched_id"]   = hit.dropna()
            out.loc[idx, "match_method"] = "name_match"
            print(f"  Step 3 name_match         : +{len(idx)}")

    # ---- Step 4: synonym fallback --------------------------------- #
    syn_col      = ref_cfg.get("synonym_col")
    alt_name_col = ref_cfg.get("alt_name_col", name_col)
    if syn_col and alt_name_col and alt_name_col in out.columns:
        unmatched = out["matched_id"].isna()
        if unmatched.any():
            delimiter = ref_cfg.get("synonym_delimiter", ";")
            syn_df    = reference_df[[id_col, syn_col]].dropna(subset=[syn_col])
            syn_df    = syn_df.assign(
                **{syn_col: syn_df[syn_col].str.split(delimiter)}
            ).explode(syn_col)
            syn_df[syn_col] = syn_df[syn_col].str.strip()

            # Drop ambiguous: same synonym token → >1 accession
            ambiguous = syn_df.duplicated(subset=[syn_col], keep=False)
            if ambiguous.sum():
                print(f"  [WARN] Step 4: dropping "
                      f"{syn_df.loc[ambiguous, syn_col].nunique()} ambiguous "
                      f"synonym token(s) ({ambiguous.sum()} rows)")
            syn_df = syn_df[~ambiguous]

            # Use dict — bypasses ArrowExtensionArray.map → get_indexer path
            syn_lookup = dict(zip(syn_df[syn_col], syn_df[id_col]))
            hit = out.loc[unmatched, alt_name_col].map(syn_lookup)
            idx = hit.dropna().index
            out.loc[idx, "matched_id"]   = hit.dropna()
            out.loc[idx, "match_method"] = "synonym_fallback"
            print(f"  Step 4 synonym_fallback   : +{len(idx)}")

    # ---- Step 5: cross-reference / Sanger xref ------------------- #
    xref_col     = ref_cfg.get("xref_col")
    src_xref_col = ref_cfg.get("src_xref_col")
    xref_prefix  = ref_cfg.get("xref_prefix")
    if xref_col and src_xref_col and xref_prefix and src_xref_col in out.columns:
        unmatched = out["matched_id"].isna()
        if unmatched.any():
            xref_to_cvcl = _parse_xref_column(
                reference_df, id_col, xref_col, prefix=xref_prefix
            )
            # Invert CVCL → ext_id to ext_id → CVCL
            ext_series = pd.Series(
                xref_to_cvcl.index.values,
                index=xref_to_cvcl.values,
            )
            # Drop ambiguous: same ext_id → >1 CVCL after inversion
            ext_dup = ext_series.index.duplicated(keep=False)
            if ext_dup.sum():
                print(f"  [WARN] Step 5: dropping "
                      f"{ext_series[ext_dup].index.nunique()} ambiguous "
                      f"xref value(s) ({ext_dup.sum()} rows)")
            ext_series = ext_series[~ext_dup]

            # Use dict — bypasses ArrowExtensionArray.map → get_indexer path
            ext_to_cvcl = ext_series.to_dict()
            hit = out.loc[unmatched, src_xref_col].map(ext_to_cvcl)
            idx = hit.dropna().index
            out.loc[idx, "matched_id"]   = hit.dropna()
            out.loc[idx, "match_method"] = f"xref_{xref_prefix.lower()}"
            print(f"  Step 5 xref_{xref_prefix.lower():<15}: +{len(idx)}")

    return out


# ------------------------------------------------------------------ #
# Main entry point                                                     #
# ------------------------------------------------------------------ #
def build_lookup_table(
    dataframes: Dict[str, pd.DataFrame],
    identifiers: Dict[str, Dict[str, Any]],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build a unified cell-line lookup table from three Track A datasets.

    Resolution chains (from EDA insights):

    Cellosaurus ↔ Sample_info (99.89% — 1,838/1,840):
      direct_id           : 1,812 via RRID ↔ Accession
      secondary_accession : +1    via exploded secondary accessions
      name_match          : +8    via safe exact name (excl. 69 ambiguous)
      xref_sanger         : +17   via Cross-references Sanger → Sanger_Model_ID
      irrecoverable       : 2     (both are NaN rows — no identifier)

    Sample_info ↔ DepMap Profiles (96.09% — 1,768/1,840):
      direct_id           : 1,749 via DepMap_ID ↔ modelid
      xref_bridge         : +19   via RRID → CVCL → cellosaurus["DepMap"] → ACH
      irrecoverable       : 72    (release version gap)

    Returns
    -------
    lookup_table   : one row per model_id, all canonical IDs + coverage flags
    unmapped_log   : rows unresolvable after all chains  (contract spec format)
    duplicate_log  : aux-table duplicate pairs logged before pivot
    """
    # ---- Identify base, aux, reference tables ---- #
    base_name = next(n for n, c in identifiers.items() if c["role"] == "base")
    base_cfg  = identifiers[base_name]
    base_key  = base_cfg["key"]
    result    = dataframes[base_name].copy()

    dup_log_rows: List[dict] = []
    flag_cols:    List[str]  = []
    aux_collapsed_store: Dict[str, pd.DataFrame] = {}

    reference_dfs = {
        n: dataframes[n]
        for n, c in identifiers.items()
        if c["role"] == "reference"
    }

    # ---- Stage A: merge aux tables ---- #
    for name, cfg in identifiers.items():
        if cfg["role"] != "aux":
            continue
        print(f"\n[AUX] Collapsing '{name}'...")
        collapsed = _collapse_aux_table(dataframes[name], cfg, dup_log_rows)
        aux_collapsed_store[name] = collapsed

        new_flag_cols = [c for c in collapsed.columns if c != cfg["merge_key"]]
        flag_cols.extend(new_flag_cols)

        result = result.merge(
            collapsed,
            left_on=cfg["base_key"],
            right_on=cfg["merge_key"],
            how="left",
        )
        # Fill boolean flags; leave string carry cols as NaN (meaningful: not profiled)
        for col in new_flag_cols:
            if result[col].dtype == bool or str(result[col].dtype) == "bool":
                result[col] = result[col].fillna(False)

        direct_match = result[cfg["base_key"]].isin(dataframes[name][cfg["merge_key"]])
        print(f"  Direct key match: {direct_match.sum()} / {len(result)} "
              f"({100*direct_match.sum()/len(result):.2f}%)")

        # Stage A2: xref bridge recovery for unmatched aux rows
        result = _recover_aux_via_xref_bridge(
            result, collapsed, cfg, reference_dfs, new_flag_cols
        )

    # ---- Stage B: cascade-match reference tables ---- #
    reference_names = [n for n, c in identifiers.items() if c["role"] == "reference"]
    for ref_name in reference_names:
        print(f"\n[REF] Cascade-matching against '{ref_name}'...")
        result = _cascade_match(result, dataframes[ref_name], identifiers[ref_name])

    # ---- Build unmapped log (per contract spec) ---- #
    if reference_names:
        unmapped_mask = result["matched_id"].isna()
        unmapped_log  = pd.DataFrame({
            "original_value":          result.loc[unmapped_mask, base_key],
            "dataset":                 base_name,
            "resolution_step_reached": "xref_sanger",
            "reason":                  "no match via direct_id / secondary_accession / "
                                       "name_match / synonym_fallback / xref_sanger",
        })
    else:
        unmapped_log = pd.DataFrame(
            columns=["original_value","dataset","resolution_step_reached","reason"]
        )

    duplicate_log = (
        pd.DataFrame(dup_log_rows) if dup_log_rows
        else pd.DataFrame(columns=["merge_key_value","collapse_on_value","reason"])
    )

    # ---- Rename matched_id → cvcl_accession (per contract schema) ---- #
    if "matched_id" in result.columns:
        result = result.rename(columns={"matched_id": "cvcl_accession"})

    # ---- Assemble final lookup table ---- #
    keep_cols = [base_key]
    for optional_col in ("name_col", "alt_name_col"):
        col_name = base_cfg.get(optional_col)
        if col_name and col_name in result.columns:
            keep_cols.append(col_name)
    keep_cols += flag_cols
    if reference_names:
        keep_cols += ["cvcl_accession", "match_method"]
    if "bridge_method" in result.columns:
        keep_cols.append("bridge_method")
    keep_cols = [c for c in keep_cols if c in result.columns]

    lookup_table = result[keep_cols].copy()

    # ---- Summary ---- #
    print(f"\n{'='*55}")
    print(f"Track A Lookup Table — Build Summary")
    print(f"{'='*55}")
    print(f"Total cell lines            : {len(lookup_table)}")
    print(f"Aux duplicate pairs logged  : {len(duplicate_log)}")
    if reference_names:
        for method in ("direct_id","secondary_accession","name_match",
                       "synonym_fallback","xref_sanger"):
            n = (lookup_table["match_method"] == method).sum()
            if n > 0:
                print(f"Resolved via {method:<22}: {n}")
    if "bridge_method" in lookup_table.columns:
        n = lookup_table["bridge_method"].notna().sum()
        print(f"Recovered via xref bridge   : {n}")
    print(f"Unmapped (irrecoverable)    : "
          f"{lookup_table['cvcl_accession'].isna().sum() if reference_names else 'N/A'}")
    print(f"Match rate                  : "
          f"{100*(1 - lookup_table['cvcl_accession'].isna().mean()):.2f}%"
          if reference_names else "")
    print(f"{'='*55}")

    return lookup_table, unmapped_log, duplicate_log
