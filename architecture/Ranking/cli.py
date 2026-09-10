"""
Ranking/cli.py  —  GeneTraceAI query interface
------------------------------------------------
"""
import sys
from pathlib import Path

import duckdb
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import PREDICTIONS, GENE_LKP, CELL_LKP, DB

TOP_N = 30
DIV   = "-" * 64
BOLD  = "=" * 64

_NEIGHBOURS_PATH = (Path(__file__).resolve().parent.parent
                    / "cell_similarity" / "outputs" / "rna_neighbours.parquet")
_SIM_TOP_K = 5
_SIM_NOTE  = "+0.186 drug-response concordance vs same-tissue baseline"


# ── helpers ────────────────────────────────────────────────────────────────

def _parse_args(raw: list[str]) -> tuple[list[str], list[str]]:
    """Split a raw arg list into positional args and --lineage values.

    --lineage lung --lineage skin  ->  lineages = ["lung", "skin"]
    Lineage terms match against both `lineage` and `lineage_subtype` columns
    (case-insensitive, substring match), so "nsclc" or "lung" both work.
    """
    positional, lineages = [], []
    i = 0
    while i < len(raw):
        if raw[i] == "--lineage" and i + 1 < len(raw):
            lineages.append(raw[i + 1].strip().lower())
            i += 2
        else:
            positional.append(raw[i])
            i += 1
    return positional, lineages


def _lineage_model_ids(terms: list[str]) -> set[str] | None:
    """Return model_ids matching any of the lineage terms, or None (no filter).

    Matches against lineage OR lineage_subtype (case-insensitive substring).
    Returns None when terms is empty so callers can skip the filter entirely.
    """
    if not terms:
        return None
    if not DB.exists():
        print("[WARNING] DuckDB not found — lineage filter skipped.")
        return None
    try:
        con = duckdb.connect(str(DB), read_only=True)
        # Build a LIKE clause for each term against both columns
        clauses = " OR ".join(
            f"(lower(lineage) LIKE '%{t}%' OR lower(lineage_subtype) LIKE '%{t}%')"
            for t in terms
        )
        rows = con.execute(
            f"SELECT lower(model_id) AS model_id FROM sample_info WHERE {clauses}"
        ).fetchall()
        con.close()
        matched = {r[0] for r in rows}
        if not matched:
            print(f"[WARNING] No lines matched lineage filter: {terms!r}")
        else:
            print(f"[LINEAGE] {', '.join(terms)!r} -> {len(matched):,} lines in scope")
        return matched
    except Exception as exc:
        print(f"[WARNING] Lineage filter failed ({exc}) — no filter applied.")
        return None


def _apply_lineage(df: pd.DataFrame, model_ids: set[str] | None,
                   col: str = "model_id") -> pd.DataFrame:
    """Filter a DataFrame to the lineage model_id set. No-op when model_ids is None."""
    if model_ids is None:
        return df
    return df[df[col].str.lower().isin(model_ids)]


def _load_genes() -> pd.DataFrame:
    return pd.read_parquet(GENE_LKP, columns=["ensg_id", "hgnc_symbol",
                                               "chromosomal_location"])


def _y_linked_set(gene_lkp: pd.DataFrame) -> set:
    """ENSG IDs on the Y chromosome, derived from HGNC chromosomal_location."""
    loc = gene_lkp["chromosomal_location"].fillna("")
    return set(gene_lkp.loc[loc.str.startswith("Y"), "ensg_id"])


def _line_sex() -> pd.DataFrame:
    """Return model_id → sex from sample_info (lowercase, null→'unknown')."""
    if not DB.exists():
        return pd.DataFrame(columns=["model_id", "sex"])
    try:
        con = duckdb.connect(str(DB), read_only=True)
        df = con.execute(
            "SELECT lower(model_id) AS model_id, "
            "lower(coalesce(sex, 'unknown')) AS sex FROM sample_info"
        ).df()
        con.close()
        return df
    except Exception:
        return pd.DataFrame(columns=["model_id", "sex"])


def _apply_sex_guard(df: pd.DataFrame, y_set: set,
                     score_col: str, ensg_ids: list[str]) -> tuple[pd.DataFrame, int]:
    """
    For any gene in ensg_ids that is Y-linked, set score_col to NaN for
    female and unknown-sex lines — structural absence, not biological signal.
    Returns (filtered_df, n_nulled).
    Consistent with the orphan-pair NaN pattern in driver_routing.py.
    """
    affected_genes = [e for e in ensg_ids if e in y_set]
    if not affected_genes:
        return df, 0

    sex_df = _line_sex()
    df = df.merge(sex_df, on="model_id", how="left")
    df["sex"] = df["sex"].fillna("unknown")

    non_male = df["sex"].isin(["female", "unknown"])
    n_nulled = int(non_male.sum())
    df.loc[non_male, score_col] = float("nan")
    return df, n_nulled


def _load_lines(cols=("model_id", "cell_line_name")) -> pd.DataFrame:
    df = pd.read_parquet(CELL_LKP, columns=list(cols))
    df["model_id"] = df["model_id"].str.lower()
    return df


def _resolve(query: str, gene_lkp: pd.DataFrame) -> tuple[str, str]:
    """Returns (ensg_id, hgnc_symbol). Accepts symbol or ENSG id."""
    q = query.strip().upper()
    m = gene_lkp[gene_lkp["hgnc_symbol"].str.upper() == q]
    if not m.empty:
        return m.iloc[0]["ensg_id"], m.iloc[0]["hgnc_symbol"]
    m2 = gene_lkp[gene_lkp["ensg_id"].str.upper() == q]
    if not m2.empty:
        return m2.iloc[0]["ensg_id"], m2.iloc[0]["hgnc_symbol"]
    raise SystemExit(f"Gene not found: {query!r}")


def _similar_lines(model_id: str, lines_df: pd.DataFrame) -> list[dict]:
    """Top-_SIM_TOP_K RNA neighbours for model_id. Returns [] if file absent."""
    if not _NEIGHBOURS_PATH.exists():
        return []
    try:
        nbrs = pd.read_parquet(_NEIGHBOURS_PATH)
        sub = (nbrs[nbrs["model_id"].str.lower() == model_id.lower()]
               .sort_values("rank")
               .head(_SIM_TOP_K)
               .merge(
                   lines_df[["model_id", "cell_line_name"]]
                   .rename(columns={"model_id": "neighbour"}),
                   on="neighbour", how="left"
               ))
        return sub.to_dict("records")
    except Exception:
        return []


def _print_similar(neighbours: list[dict], label: str = "") -> None:
    if not neighbours:
        return
    tag = f" for {label}" if label else ""
    print(f"\n  RNA alternatives{tag}  [{_SIM_NOTE}]:")
    for r in neighbours:
        name = str(r.get("cell_line_name", "?"))[:27]
        nbr  = str(r.get("neighbour", "?")).upper()
        sim  = float(r.get("similarity", float("nan")))
        print(f"    {nbr:<13}  {name:<28}  sim={sim:.3f}")


def _check_outputs():
    for p in [PREDICTIONS, GENE_LKP, CELL_LKP]:
        if not p.exists():
            raise SystemExit(f"Missing: {p}\nRun scoring stages first (python run_all.py).")


def _lineage_meta(model_id: str) -> dict:
    """Fetch lineage + clinical fields from DuckDB."""
    if not DB.exists():
        return {}
    try:
        con = duckdb.connect(str(DB), read_only=True)
        row = con.execute(
            "SELECT lineage, lineage_subtype, primary_disease, sex, age, "
            "default_growth_pattern, primary_or_metastasis "
            "FROM main.sample_info WHERE model_id = ?",
            [model_id.lower()]
        ).fetchone()
        con.close()
        if row:
            keys = ["lineage", "lineage_subtype", "primary_disease", "sex",
                    "age", "growth_pattern", "primary_or_metastasis"]
            return {k: v for k, v in zip(keys, row) if v is not None}
    except Exception:
        pass
    return {}


# ── commands ───────────────────────────────────────────────────────────────

def cmd_gene(args: list) -> None:
    args, lineages = _parse_args(args)
    if not args:
        raise SystemExit("Usage: cli.py gene <GENE> [ACH-ID] [--lineage <tissue>]")
    _check_outputs()

    genes = _load_genes()
    y_set = _y_linked_set(genes)
    ensg, sym = _resolve(args[0], genes)
    lineage_ids = _lineage_model_ids(lineages)
    pred = pd.read_parquet(PREDICTIONS)
    pred = _apply_lineage(pred, lineage_ids)
    sub  = pred[pred.ensg_id == ensg].copy()
    if sub.empty:
        raise SystemExit(f"No predictions for {sym} ({ensg})")

    sub = sub.sort_values("core_score", ascending=False).reset_index(drop=True)
    lines_full = pd.read_parquet(CELL_LKP)
    lines_full["model_id"] = lines_full["model_id"].str.lower()
    sub = sub.merge(lines_full[["model_id", "cell_line_name"]], on="model_id", how="left")

    # ── Detail view ─────────────────────────────────────────────────────────
    if len(args) >= 2:
        ach = args[1].strip()
        row = sub[sub.model_id.str.lower() == ach.lower()]
        if row.empty:
            raise SystemExit(f"No prediction found for {sym} + {ach}")
        row = row.iloc[0]
        rank_pos = sub[sub.model_id.str.lower() == ach.lower()].index[0] + 1

        meta = lines_full[lines_full.model_id.str.lower() == ach.lower()]
        meta = meta.iloc[0].to_dict() if not meta.empty else {}
        lin  = _lineage_meta(ach)

        print(f"\n{BOLD}")
        print(f"  Gene:         {sym}  ({ensg})")
        print(f"  Cell Line:    {meta.get('cell_line_name', '?')}  ({row['model_id'].upper()})")
        print(BOLD)
        print(f"  Score:        {row['core_score']:.4f}")
        print(f"  Tier:         {row.get('confidence_tier', '?')}")
        print(f"  Rank:         {rank_pos} / {len(sub):,}")
        n_lay = int(row.get("n_layers", 0))
        layer_label = "RNA + protein" if n_lay == 2 else "RNA only"
        print(f"  Data layers:  {n_lay} ({layer_label})")
        print(f"  Driver alt:   {'Yes' if row.get('has_driver_alteration', False) else 'No'}")
        if float(row.get("p_mutation", 0)) > 0:
            print(f"  p_mutation:   {float(row['p_mutation']):.3f}")
        if float(row.get("p_fusion", 0)) > 0:
            print(f"  p_fusion:     {float(row['p_fusion']):.3f}")
        if row.get("has_cna_alteration", False):
            print(f"  CNA:          Yes")
        print(DIV)
        print(f"  Lineage:      {lin.get('lineage', '?')}")
        if lin.get("lineage_subtype"):
            print(f"  Subtype:      {lin['lineage_subtype']}")
        if lin.get("primary_disease"):
            print(f"  Disease:      {lin['primary_disease']}")
        if lin.get("primary_or_metastasis"):
            print(f"  Site:         {lin['primary_or_metastasis']}")
        if lin.get("sex"):
            print(f"  Sex:          {lin['sex']}")
        if lin.get("age"):
            print(f"  Age:          {lin['age']}")
        if lin.get("growth_pattern"):
            print(f"  Growth:       {lin['growth_pattern']}")
        cellosaurus_dis = meta.get("diseases", "")
        if cellosaurus_dis:
            print(f"  Cellosaurus:  {str(cellosaurus_dis)[:80]}")
        nbrs = _similar_lines(row["model_id"], lines_full[["model_id", "cell_line_name"]])
        _print_similar(nbrs)
        print(BOLD + "\n")
        return

    # ── List view ────────────────────────────────────────────────────────────
    lineage_tag = f"  |  Lineage: {', '.join(lineages)}" if lineages else ""
    print(f"\nGene: {sym}  ({ensg}){lineage_tag}")
    print(DIV)
    print(f"  {'ACH ID':<13} {'Cell Line Name':<30} {'Score':>6}  Tier")
    print(DIV)
    for _, row in sub.head(TOP_N).iterrows():
        name = str(row.get("cell_line_name", "?"))[:29]
        print(f"  {row['model_id']:<13} {name:<30} {row['core_score']:>6.3f}  {row.get('confidence_tier', '?')}")
    print(f"\n  Showing top {min(TOP_N, len(sub))} of {len(sub):,} lines")
    if not sub.empty:
        top = sub.iloc[0]
        nbrs = _similar_lines(top["model_id"], lines_full[["model_id", "cell_line_name"]])
        _print_similar(nbrs, label=f"{top['model_id'].upper()} (rank 1)")


def cmd_genes(gene_queries: list) -> None:
    gene_queries, lineages = _parse_args(gene_queries)
    if len(gene_queries) < 2:
        raise SystemExit("Usage: cli.py genes <GENE_A> <GENE_B> [GENE_C ...] [--lineage <tissue>]")
    _check_outputs()

    genes = _load_genes()
    y_set = _y_linked_set(genes)
    resolved = [_resolve(q, genes) for q in gene_queries]
    ensgs = [e for e, _ in resolved]
    syms  = [s for _, s in resolved]

    lineage_ids = _lineage_model_ids(lineages)
    pred  = pd.read_parquet(PREDICTIONS, columns=["model_id", "ensg_id", "core_score"])
    pred  = _apply_lineage(pred, lineage_ids)
    lines = _load_lines()

    wide = None
    for ensg, sym in zip(ensgs, syms):
        s = pred[pred.ensg_id == ensg][["model_id", "core_score"]].rename(columns={"core_score": sym})
        wide = s if wide is None else wide.merge(s, on="model_id", how="inner")

    if wide is None or wide.empty:
        raise SystemExit("No shared cell lines found across all requested genes.")

    # Sex guard: Y-linked genes are structurally absent in female/unknown-sex lines.
    # Null joint_score for affected pairs rather than letting structural absence
    # depress co-selection scores (false negatives).
    y_in_query = [e for e in ensgs if e in y_set]
    if y_in_query:
        y_syms = [syms[ensgs.index(e)] for e in y_in_query]
        sex_df = _line_sex()
        wide = wide.merge(sex_df, on="model_id", how="left")
        wide["sex"] = wide["sex"].fillna("unknown")
        non_male = wide["sex"].isin(["female", "unknown"])
        n_nulled = int(non_male.sum())
        for sym in y_syms:
            wide.loc[non_male, sym] = float("nan")
        y_names = ", ".join(syms[ensgs.index(e)] for e in y_in_query)
        print(f"\n[SEX GUARD] {y_names} is Y-linked. "
              f"{n_nulled} female/unknown-sex pairs set to NaN "
              f"(structural absence, not biological co-selection).")
        wide = wide.drop(columns=["sex"])

    JOINT_FLOOR = 0.50
    wide["joint_score"] = wide[syms].min(axis=1, skipna=False)
    wide = wide[wide["joint_score"].notna() & (wide["joint_score"] >= JOINT_FLOOR)]
    wide = wide.sort_values("joint_score", ascending=False)
    wide = wide.merge(lines, on="model_id", how="left").reset_index(drop=True)

    sym_hdrs = "  ".join(f"{s[:7]:>7}" for s in syms)
    lineage_tag = f"  |  Lineage: {', '.join(lineages)}" if lineages else ""
    print(f"\nCo-selection: {' AND '.join(syms)}{lineage_tag}")
    print(f"Score = min({', '.join(syms)})  |  Floor = {JOINT_FLOOR}")
    print(DIV)
    print(f"  {'ACH ID':<13} {'Cell Line Name':<28} {'Joint':>6}  {sym_hdrs}")
    print(DIV)
    for _, row in wide.head(TOP_N).iterrows():
        name = str(row.get("cell_line_name", "?"))[:27]
        scores = "  ".join(f"{row[s]:>7.3f}" for s in syms)
        print(f"  {row['model_id']:<13} {name:<28} {row['joint_score']:>6.3f}  {scores}")
    print(f"\n  Lines clearing floor: {len(wide):,}")


def cmd_exclude(args: list) -> None:
    args, lineages = _parse_args(args)
    if len(args) < 2:
        raise SystemExit("Usage: cli.py exclude <GENE_A> <GENE_B> [--lineage <tissue>]")
    _check_outputs()

    genes = _load_genes()
    y_set = _y_linked_set(genes)
    ensg_a, sym_a = _resolve(args[0], genes)
    ensg_b, sym_b = _resolve(args[1], genes)

    lineage_ids = _lineage_model_ids(lineages)
    pred  = pd.read_parquet(PREDICTIONS, columns=["model_id", "ensg_id", "core_score"])
    pred  = _apply_lineage(pred, lineage_ids)
    lines = _load_lines()

    a = pred[pred.ensg_id == ensg_a][["model_id", "core_score"]].rename(columns={"core_score": "score_a"})
    b = pred[pred.ensg_id == ensg_b][["model_id", "core_score"]].rename(columns={"core_score": "score_b"})

    merged = a.merge(b, on="model_id", how="inner")

    # Sex guard: if gene_B is Y-linked, female and unknown-sex lines are structurally
    # absent for gene_B — not biologically non-essential. Set selectivity to NaN and
    # exclude from ranking (consistent with orphan-pair handling in driver_routing.py).
    n_nulled = 0
    if ensg_b in y_set:
        sex_df = _line_sex()
        merged = merged.merge(sex_df, on="model_id", how="left")
        merged["sex"] = merged["sex"].fillna("unknown")
        non_male = merged["sex"].isin(["female", "unknown"])
        n_nulled = int(non_male.sum())
        merged.loc[non_male, "score_b"] = float("nan")
        merged = merged.drop(columns=["sex"])
        print(f"\n[SEX GUARD] {sym_b} is Y-linked (chromosomal absence in XX lines). "
              f"{n_nulled} female/unknown-sex pairs excluded from ranking "
              f"(selectivity=NaN — structural absence, not biological signal).")

    merged["selectivity"] = merged["score_a"] * (1.0 - merged["score_b"])
    # NaN selectivity (from guard above) sorts to bottom naturally with na_position='last'
    merged = merged.sort_values("selectivity", ascending=False, na_position="last")
    scored = merged[merged["selectivity"].notna()]
    merged = merged.merge(lines, on="model_id", how="left").reset_index(drop=True)

    print(f"\nExclusion: {sym_a} HIGH  ×  {sym_b} LOW")
    print(f"Selectivity = score_{sym_a} × (1 − score_{sym_b})")
    if n_nulled:
        print(f"[Note: {n_nulled} female/unknown-sex pairs hidden — Y-linked gene_B]")
    print(DIV)
    print(f"  {'ACH ID':<13} {'Cell Line Name':<28} {sym_a[:7]:>7}  {sym_b[:7]:>7}  {'Select.':>7}")
    print(DIV)
    display = merged[merged["selectivity"].notna()]
    for _, row in display.head(TOP_N).iterrows():
        name = str(row.get("cell_line_name", "?"))[:27]
        print(f"  {row['model_id']:<13} {name:<28} {row['score_a']:>7.3f}  {row['score_b']:>7.3f}  {row['selectivity']:>7.3f}")
    print(f"\n  Ranked pairs: {len(scored):,}  |  Excluded (sex guard): {n_nulled:,}")
    if not display.empty:
        top = display.iloc[0]
        nbrs = _similar_lines(top["model_id"], lines)
        _print_similar(nbrs, label=f"{top['model_id'].upper()} (rank 1)")


# ── entry point ────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)

    cmd  = sys.argv[1].lower()
    rest = sys.argv[2:]

    dispatch = {"gene": cmd_gene, "genes": cmd_genes, "exclude": cmd_exclude}
    if cmd not in dispatch:
        print(f"Unknown command: {cmd!r}")
        print("Valid commands: gene | genes | exclude")
        sys.exit(1)
    dispatch[cmd](rest)


if __name__ == "__main__":
    main()
