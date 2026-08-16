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


# ── helpers ────────────────────────────────────────────────────────────────

def _load_genes() -> pd.DataFrame:
    return pd.read_parquet(GENE_LKP, columns=["ensg_id", "hgnc_symbol"])


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
    if not args:
        raise SystemExit("Usage: cli.py gene <GENE> [ACH-ID]")
    _check_outputs()

    genes = _load_genes()
    ensg, sym = _resolve(args[0], genes)
    pred = pd.read_parquet(PREDICTIONS)
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
        print(BOLD + "\n")
        return

    # ── List view ────────────────────────────────────────────────────────────
    print(f"\nGene: {sym}  ({ensg})")
    print(DIV)
    print(f"  {'ACH ID':<13} {'Cell Line Name':<30} {'Score':>6}  Tier")
    print(DIV)
    for _, row in sub.head(TOP_N).iterrows():
        name = str(row.get("cell_line_name", "?"))[:29]
        print(f"  {row['model_id']:<13} {name:<30} {row['core_score']:>6.3f}  {row.get('confidence_tier', '?')}")
    print(f"\n  Showing top {min(TOP_N, len(sub))} of {len(sub):,} lines")


def cmd_genes(gene_queries: list) -> None:
    if len(gene_queries) < 2:
        raise SystemExit("Usage: cli.py genes <GENE_A> <GENE_B> [GENE_C ...]")
    _check_outputs()

    genes = _load_genes()
    resolved = [_resolve(q, genes) for q in gene_queries]
    ensgs = [e for e, _ in resolved]
    syms  = [s for _, s in resolved]

    pred  = pd.read_parquet(PREDICTIONS, columns=["model_id", "ensg_id", "core_score"])
    lines = _load_lines()

    wide = None
    for ensg, sym in zip(ensgs, syms):
        s = pred[pred.ensg_id == ensg][["model_id", "core_score"]].rename(columns={"core_score": sym})
        wide = s if wide is None else wide.merge(s, on="model_id", how="inner")

    if wide is None or wide.empty:
        raise SystemExit("No shared cell lines found across all requested genes.")

    JOINT_FLOOR = 0.50
    wide["joint_score"] = wide[syms].min(axis=1)
    wide = wide[wide["joint_score"] >= JOINT_FLOOR].sort_values("joint_score", ascending=False)
    wide = wide.merge(lines, on="model_id", how="left").reset_index(drop=True)

    sym_hdrs = "  ".join(f"{s[:7]:>7}" for s in syms)
    print(f"\nCo-selection: {' AND '.join(syms)}")
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
    if len(args) < 2:
        raise SystemExit("Usage: cli.py exclude <GENE_A> <GENE_B>")
    _check_outputs()

    genes = _load_genes()
    ensg_a, sym_a = _resolve(args[0], genes)
    ensg_b, sym_b = _resolve(args[1], genes)

    pred  = pd.read_parquet(PREDICTIONS, columns=["model_id", "ensg_id", "core_score"])
    lines = _load_lines()

    a = pred[pred.ensg_id == ensg_a][["model_id", "core_score"]].rename(columns={"core_score": "score_a"})
    b = pred[pred.ensg_id == ensg_b][["model_id", "core_score"]].rename(columns={"core_score": "score_b"})

    merged = a.merge(b, on="model_id", how="inner")
    merged["selectivity"] = merged["score_a"] * (1.0 - merged["score_b"])
    merged = merged.sort_values("selectivity", ascending=False)
    merged = merged.merge(lines, on="model_id", how="left").reset_index(drop=True)

    print(f"\nExclusion: {sym_a} HIGH  ×  {sym_b} LOW")
    print(f"Selectivity = score_{sym_a} × (1 − score_{sym_b})")
    print(DIV)
    print(f"  {'ACH ID':<13} {'Cell Line Name':<28} {sym_a[:7]:>7}  {sym_b[:7]:>7}  {'Select.':>7}")
    print(DIV)
    for _, row in merged.head(TOP_N).iterrows():
        name = str(row.get("cell_line_name", "?"))[:27]
        print(f"  {row['model_id']:<13} {name:<28} {row['score_a']:>7.3f}  {row['score_b']:>7.3f}  {row['selectivity']:>7.3f}")
    print(f"\n  Total pairs evaluated: {len(merged):,}")


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
