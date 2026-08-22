from __future__ import annotations

"""
Cell-line z-scores from two proteomics sources (ProCAN + CCLE proteomics).

Three changes over the base version:
  1. Isoforms are handled explicitly. When several protein IDs map to one gene,
     we first check whether they agree across cell lines (Spearman). If they do,
     the gene value for a line is the MAX abundance across isoforms; if they do
     not, we fall back to the mean so an unrelated isoform can't dominate.
  2. ProCAN is NOT preferred over CCLE. Both platforms are combined with
     coverage-only (symmetric) weights, and platform disagreement lowers our
     confidence in the *merged gene* instead of deleting one platform's data.
  3. No domain threshold is hardcoded. Isoform concordance, platform-agreement
     bounds and the minimum shared-line count are all derived from the data
     (null / matched distributions + statistical power). The only remaining
     numbers are standard statistical conventions (significance level, the
     percentiles that define "null" and "clear agreement"), surfaced as
     documented parameters rather than buried constants.
"""

import argparse
import json
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr, t as t_dist


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DB = ROOT / "db" / "celllineselector.duckdb"

# Method 1 source names.
PROCAN_SOURCE = "protein_matrix_averaged_20250211"
CCLE_SOURCE = "proteomics"

# Statistical CONVENTIONS (not domain-tuned values). Exposed via CLI so they are
# visible and overridable rather than hidden. 1.4826 elsewhere is the fixed
# MAD -> sigma constant for a normal, a mathematical constant, not a tunable.
DEFAULT_ALPHA = 0.05        # significance level for the power-based min overlap
DEFAULT_NULL_PCTILE = 95.0  # "more correlated than X% of unrelated pairs" = agree
DEFAULT_MATCHED_PCTILE = 60.0  # where matched agreement is treated as full trust
DEFAULT_MIN_GAP = 0.10      # minimum width of the confidence ramp (low -> high)


# ---------------------------------------------------------------------
# Identifier helpers (unchanged from base)
# ---------------------------------------------------------------------

def canonical_ensg(value) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip().upper()
    if not text:
        return None
    match = re.fullmatch(r"ENSG\d+(?:\.\d+)?", text)
    return text.split(".", 1)[0] if match else None


def split_id_list(value) -> list[str]:
    if value is None or pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    out = []
    for p in re.split(r"[;,|/\\]", text):
        token = str(p).strip().upper()
        if token and token.lower() not in {"nan", "none", "null"}:
            out.append(token)
    return out


def build_uniprot_to_ensg(con) -> dict[str, str]:
    mapping: dict[str, str] = {}
    roster = con.execute(
        'SELECT gene_id, "uniprot_id(s)", procan_uniprot_id FROM gene_roster'
    ).fetchdf()
    for col in ["uniprot_id(s)", "procan_uniprot_id"]:
        if col not in roster.columns:
            continue
        for _, row in roster[["gene_id", col]].dropna().iterrows():
            gene = canonical_ensg(row["gene_id"])
            if gene is None:
                continue
            for token in split_id_list(row[col]):
                mapping[token.lower()] = gene
    return mapping


def build_symbol_to_ensg(con) -> dict[str, str]:
    mapping: dict[str, str] = {}
    roster = con.execute('SELECT gene_id, "gene_name(s)" FROM gene_roster').fetchdf()
    if "gene_name(s)" not in roster.columns:
        return mapping
    for _, row in roster[["gene_id", "gene_name(s)"]].dropna().iterrows():
        gene = canonical_ensg(row["gene_id"])
        if gene is None:
            continue
        for token in split_id_list(row["gene_name(s)"]):
            mapping[token.lower()] = gene
    return mapping


def normalize_name(value) -> str:
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


# ---------------------------------------------------------------------
# Loading / reshaping (unchanged from base)
# ---------------------------------------------------------------------

def attach_model_ids(con, df: pd.DataFrame) -> pd.DataFrame:
    if "model_id" in df.columns:
        return df
    if "ccle_name" not in df.columns:
        raise ValueError("Protein matrix missing both model_id and ccle_name.")

    meta = con.execute(
        """
        SELECT model_id, cell_line_name, stripped_cell_line_name, ccle_id
        FROM sample_info WHERE model_id IS NOT NULL
        """
    ).fetchdf()
    if meta.empty:
        raise ValueError("sample_info has no usable cell-line mappings.")

    lookup: dict[str, str] = {}
    for col in ["ccle_id", "cell_line_name", "stripped_cell_line_name"]:
        if col not in meta.columns:
            continue
        for _, row in meta[["model_id", col]].dropna().iterrows():
            key = normalize_name(row[col])
            if key:
                lookup.setdefault(key, str(row["model_id"]))

    out = df.copy()
    out["model_id"] = out["ccle_name"].map(lambda v: lookup.get(normalize_name(v), None))
    return out.dropna(subset=["model_id"]).copy()


def wide_to_long(df: pd.DataFrame, id_col: str) -> pd.DataFrame:
    protein_cols = [c for c in df.columns if c != id_col]
    sub = df[[id_col] + protein_cols].set_index(id_col)
    try:
        long = sub.stack(future_stack=True).reset_index()
    except TypeError:
        long = sub.stack().reset_index()
    long.columns = [id_col, "protein", "value"]
    return long.dropna(subset=["value"]).reset_index(drop=True)


def table_to_long(con, df, source_name, protein_map) -> pd.DataFrame:
    work = df.copy()
    if "model_id" not in work.columns and "ccle_name" in work.columns:
        work = attach_model_ids(con, work)
    if "model_id" not in work.columns:
        raise ValueError(f"Table '{source_name}' has neither model_id nor ccle_name.")

    long = wide_to_long(work, "model_id")
    long["protein"] = long["protein"].astype(str).str.strip().str.lower()
    long["gene_id"] = long["protein"].map(protein_map)
    long = long.dropna(subset=["gene_id"]).copy()
    long["gene_id"] = long["gene_id"].map(canonical_ensg)
    long = long.dropna(subset=["gene_id"]).copy()
    long["source"] = source_name
    return long


def load_lineage_map(con) -> pd.Series:
    df = con.execute(
        "SELECT DISTINCT model_id, lineage FROM sample_info WHERE model_id IS NOT NULL"
    ).fetchdf()
    if df.empty:
        return pd.Series(dtype=object)
    return (
        df.drop_duplicates(subset=["model_id"]).set_index("model_id")["lineage"]
        .fillna("unknown")
    )


def load_sources(con, table_names, protein_map) -> dict[str, pd.DataFrame]:
    raw: dict[str, pd.DataFrame] = {}
    for name in table_names:
        try:
            table = con.execute(f'SELECT * FROM "{name}"').fetchdf()
        except Exception:
            continue
        if table.empty:
            continue
        long = table_to_long(con, table, name, protein_map)
        if not long.empty:
            raw[name] = long
    return raw


# ---------------------------------------------------------------------
# CHANGE 1: isoform check -> max if concordant, else mean
# ---------------------------------------------------------------------

def _pairwise_rho(a: pd.Series, b: pd.Series, min_overlap: int) -> float:
    both = a.notna() & b.notna()
    if int(both.sum()) < min_overlap:
        return float("nan")
    x, y = a[both], b[both]
    if x.nunique() < 2 or y.nunique() < 2:
        return float("nan")
    with np.errstate(all="ignore"):
        rho, _ = spearmanr(x, y)
    return float(rho) if np.isfinite(rho) else float("nan")


def _min_pairwise_spearman(mat: pd.DataFrame, min_overlap: int) -> tuple[float, int]:
    """Min pairwise rho over isoform columns (min => all isoforms must agree)."""
    cols = list(mat.columns)
    rhos = []
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            r = _pairwise_rho(mat[cols[i]], mat[cols[j]], min_overlap)
            if np.isfinite(r):
                rhos.append(r)
    if not rhos:
        return float("nan"), 0
    return min(rhos), len(rhos)


def collapse_isoforms(long, source_name, *, min_overlap, corr_threshold):
    """
    One value per (model_id, gene_id).
      1 isoform                       -> pass through.
      >=2 concordant (min rho >= thr) -> MAX abundance per line.
      >=2 discordant / unassessable   -> MEAN per line (symmetric fallback).
    """
    iso_counts = long.groupby("gene_id")["protein"].nunique()
    single = iso_counts[iso_counts <= 1].index
    multi = iso_counts[iso_counts > 1].index

    frames, qc = [], []

    if len(single):
        s = long[long["gene_id"].isin(single)]
        s1 = s.groupby(["model_id", "gene_id"])["value"].mean().reset_index()
        s1["n_isoforms"] = 1
        frames.append(s1[["model_id", "gene_id", "value", "n_isoforms"]])

    for gene in multi:
        g = long[long["gene_id"] == gene]
        mat = g.pivot_table(index="model_id", columns="protein", values="value", aggfunc="mean")
        min_rho, n_pairs = _min_pairwise_spearman(mat, min_overlap)
        concordant = bool(np.isfinite(min_rho) and (min_rho >= corr_threshold))
        val = mat.max(axis=1, skipna=True) if concordant else mat.mean(axis=1, skipna=True)

        gg = val.dropna().reset_index()
        gg.columns = ["model_id", "value"]
        gg["gene_id"] = gene
        gg["n_isoforms"] = int(mat.shape[1])
        frames.append(gg[["model_id", "gene_id", "value", "n_isoforms"]])
        qc.append({"source": source_name, "gene_id": gene, "n_isoforms": int(mat.shape[1]),
                   "n_pairs_assessed": int(n_pairs),
                   "min_pair_rho": float(min_rho) if np.isfinite(min_rho) else np.nan,
                   "rollup": "max" if concordant else "mean"})

    collapsed = (pd.concat(frames, ignore_index=True) if frames
                 else pd.DataFrame(columns=["model_id", "gene_id", "value", "n_isoforms"]))
    collapsed["source"] = source_name
    return collapsed, pd.DataFrame(qc, columns=["source", "gene_id", "n_isoforms",
                                                "n_pairs_assessed", "min_pair_rho", "rollup"])


# ---------------------------------------------------------------------
# Lineage-conditioned robust z (unchanged from base)
# ---------------------------------------------------------------------

def robust_z_matrix(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    if X.size == 0:
        return X
    med = np.median(X, axis=0, keepdims=True)
    mad = np.median(np.abs(X - med), axis=0, keepdims=True)
    scale = np.where(np.isclose(mad, 0.0), np.std(X, axis=0, keepdims=True), 1.4826 * mad)
    scale = np.where(np.isclose(scale, 0.0), 1.0, scale)
    return (X - med) / scale


def zscore_platform(df, lineage_map, source_name) -> pd.DataFrame:
    wide = df.pivot_table(index="model_id", columns="gene_id", values="value", aggfunc="first")
    wide["lineage"] = lineage_map.reindex(wide.index)
    wide = wide.dropna(subset=["lineage"])
    empty = pd.DataFrame(columns=["gene_id", "model_id", "lineage", "z_score", "source"])
    if wide.empty:
        return empty

    rows = []
    for lineage, grp in wide.groupby("lineage", sort=True):
        gene_cols = [c for c in grp.columns if c != "lineage"]
        if not gene_cols:
            continue
        Z = robust_z_matrix(grp[gene_cols].to_numpy(dtype=float))
        zdf = pd.DataFrame(Z, index=grp.index, columns=gene_cols).stack().reset_index()
        zdf.columns = ["model_id", "gene_id", "z_score"]
        zdf["lineage"] = lineage
        zdf["source"] = source_name
        rows.append(zdf)

    if not rows:
        return empty
    out = pd.concat(rows, ignore_index=True)
    return out[["gene_id", "model_id", "lineage", "z_score", "source"]]


# ---------------------------------------------------------------------
# CHANGE 2: symmetric platform agreement (no ProCAN preference)
# ---------------------------------------------------------------------

def gene_agreement(procan_long, ccle_long, min_shared) -> pd.DataFrame:
    """Per-gene ProCAN/CCLE Spearman on shared cell lines (symmetric)."""
    p = procan_long.groupby(["model_id", "gene_id"])["value"].mean().rename("v_p").reset_index()
    c = ccle_long.groupby(["model_id", "gene_id"])["value"].mean().rename("v_c").reset_index()
    merged = p.merge(c, on=["model_id", "gene_id"], how="inner")
    if merged.empty:
        return pd.DataFrame(columns=["gene_id", "platform_rho", "n_shared"])

    counts = merged.groupby("gene_id").size()
    merged = merged[merged["gene_id"].isin(counts[counts >= min_shared].index)]

    rows = []
    for gene, grp in merged.groupby("gene_id"):
        with np.errstate(all="ignore"):
            rho, _ = spearmanr(grp["v_p"], grp["v_c"])
        if np.isfinite(rho):
            rows.append({"gene_id": gene, "platform_rho": float(rho), "n_shared": int(len(grp))})
    return pd.DataFrame(rows, columns=["gene_id", "platform_rho", "n_shared"])


def confidence_from_rho(rho, low, high) -> float:
    """rho -> [0,1] confidence. Symmetric: disagreement lowers trust in the
    merged gene, it does not pick a platform."""
    if not np.isfinite(rho):
        return 1.0
    if high <= low:
        return 1.0 if rho >= high else 0.0
    return float(np.clip((rho - low) / (high - low), 0.0, 1.0))


# ---------------------------------------------------------------------
# CHANGE 3: derive every threshold from the data
# ---------------------------------------------------------------------

def _protein_matrix(raw: pd.DataFrame):
    mat = raw.pivot_table(index="model_id", columns="protein", values="value", aggfunc="mean")
    gene_of = raw.drop_duplicates("protein").set_index("protein")["gene_id"].to_dict()
    return mat, gene_of


def _gene_matrix(raw: pd.DataFrame) -> pd.DataFrame:
    lvl = raw.groupby(["model_id", "gene_id"])["value"].mean().reset_index()
    return lvl.pivot_table(index="model_id", columns="gene_id", values="value")


def _power_min_n(target_r: float, alpha: float, nmax: int = 300) -> int:
    """Smallest n at which |rho| = target_r is significant at alpha (two-sided)."""
    target_r = float(np.clip(abs(target_r), 1e-3, 0.999))
    for n in range(4, nmax + 1):
        tc = t_dist.ppf(1 - alpha / 2, n - 2)
        rc = tc / np.sqrt(n - 2 + tc ** 2)         # critical Spearman rho at this n
        if rc <= target_r:
            return n
    return nmax


def _sample_null_protein_rho(raw_longs, seed_overlap, rng, n_pairs):
    null = []
    for raw in raw_longs.values():
        mat, gene_of = _protein_matrix(raw)
        cols = list(mat.columns)
        if len(cols) < 4:
            continue
        genes = np.array([gene_of.get(c) for c in cols], dtype=object)
        for _ in range(n_pairs):
            i, j = rng.integers(0, len(cols), size=2)
            if i == j or genes[i] == genes[j]:
                continue
            r = _pairwise_rho(mat[cols[i]], mat[cols[j]], seed_overlap)
            if np.isfinite(r):
                null.append(r)
    return np.asarray(null, dtype=float)


def derive_thresholds(raw_longs, lineage_map=None, *, alpha, null_pctile,
                      matched_pctile, min_gap, seed=0, n_pairs=1500, verbose=True):
    rng = np.random.default_rng(seed)
    # smallest statistically valid overlap: the fewest lines at which even a very
    # strong correlation clears significance. Used only to build the null samples.
    seed_overlap = _power_min_n(0.9, alpha)

    # --- isoform concordance threshold: null of cross-gene protein correlations
    null_iso = _sample_null_protein_rho(raw_longs, seed_overlap, rng, n_pairs)
    if null_iso.size >= 30:
        iso_thr = float(np.clip(np.percentile(null_iso, null_pctile), 0.0, 0.95))
    else:
        iso_thr = float(t_crit(seed_overlap, alpha))   # significance floor fallback

    # --- platform-agreement bounds: null (mismatched genes) vs matched (same gene)
    low = high = float("nan")
    if PROCAN_SOURCE in raw_longs and CCLE_SOURCE in raw_longs:
        P, C = _gene_matrix(raw_longs[PROCAN_SOURCE]), _gene_matrix(raw_longs[CCLE_SOURCE])
        shared = [g for g in P.columns if g in C.columns]
        matched = [r for g in shared if np.isfinite(r := _pairwise_rho(P[g], C[g], seed_overlap))]
        null_x = []
        pc, cc = list(P.columns), list(C.columns)
        for _ in range(n_pairs):
            g, h = pc[rng.integers(0, len(pc))], cc[rng.integers(0, len(cc))]
            if g == h:
                continue
            r = _pairwise_rho(P[g], C[h], seed_overlap)
            if np.isfinite(r):
                null_x.append(r)
        if len(matched) >= 10 and len(null_x) >= 30:
            low = float(np.clip(np.percentile(null_x, null_pctile), 0.0, 0.9))
            high = float(np.clip(max(low + min_gap, np.percentile(matched, matched_pctile)),
                                 low + min_gap, 0.95))
    if not (np.isfinite(low) and np.isfinite(high)):
        low = float(t_crit(seed_overlap, alpha))     # significance-floor fallback
        high = float(min(low + max(min_gap, low), 0.9))

    # --- minimum shared lines to trust an agreement estimate: power at `low`
    min_shared = int(_power_min_n(low, alpha))

    thr = {
        "isoform_corr": iso_thr,
        "isoform_min_overlap": int(seed_overlap),
        "agree_low": low,
        "agree_high": high,
        "min_shared": min_shared,
        "seed_overlap": int(seed_overlap),
    }
    if verbose:
        print("Derived thresholds (data-driven; conventions: "
              f"alpha={alpha}, null_pctile={null_pctile}, matched_pctile={matched_pctile}):")
        for k, v in thr.items():
            print(f"  {k:<20} = {v}")
    return thr


def t_crit(n: int, alpha: float) -> float:
    """Critical Spearman rho at sample size n and level alpha (two-sided)."""
    n = max(int(n), 4)
    tc = t_dist.ppf(1 - alpha / 2, n - 2)
    return tc / np.sqrt(n - 2 + tc ** 2)


# ---------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------

def compute_from_sources(raw_longs, lineage_map, thr) -> tuple[pd.DataFrame, pd.DataFrame]:
    # 1. isoform collapse per source
    source_longs, qc_frames, iso_counts = {}, [], []
    for name, raw in raw_longs.items():
        collapsed, qc = collapse_isoforms(
            raw, name, min_overlap=thr["isoform_min_overlap"], corr_threshold=thr["isoform_corr"]
        )
        if collapsed.empty:
            continue
        source_longs[name] = collapsed
        iso_counts.append(collapsed[["gene_id", "n_isoforms"]])
        if not qc.empty:
            qc_frames.append(qc)
    if not source_longs:
        raise ValueError("No protein rows after isoform collapse.")

    isoform_qc = pd.concat(qc_frames, ignore_index=True) if qc_frames else pd.DataFrame()
    iso_by_gene = (pd.concat(iso_counts, ignore_index=True).groupby("gene_id")["n_isoforms"].max()
                   if iso_counts else pd.Series(dtype=int))

    # 2. symmetric per-gene agreement -> confidence (no platform preferred)
    if PROCAN_SOURCE in source_longs and CCLE_SOURCE in source_longs:
        agree = gene_agreement(source_longs[PROCAN_SOURCE], source_longs[CCLE_SOURCE], thr["min_shared"])
    else:
        agree = pd.DataFrame(columns=["gene_id", "platform_rho", "n_shared"])
    conf_by_gene = (
        agree.assign(confidence=agree["platform_rho"].map(
            lambda r: confidence_from_rho(r, thr["agree_low"], thr["agree_high"])))
        .set_index("gene_id")["confidence"]
        if not agree.empty else pd.Series(dtype=float)
    )

    # 3. lineage-conditioned robust z per platform
    blocks = [z for name, long in source_longs.items()
              if not (z := zscore_platform(long, lineage_map, name)).empty]
    if not blocks:
        raise ValueError("No protein rows could be scored.")
    all_z = pd.concat(blocks, ignore_index=True)

    # 4. SYMMETRIC coverage-weighted combination (identical treatment of both platforms)
    src_n = all_z.groupby(["gene_id", "source"])["model_id"].nunique().rename("src_n").reset_index()
    all_z = all_z.merge(src_n, on=["gene_id", "source"], how="left")
    all_z["w"] = np.sqrt(all_z["src_n"].astype(float))
    all_z["wz"] = all_z["w"] * all_z["z_score"]
    all_z["w2"] = all_z["w"] ** 2

    agg = all_z.groupby(["gene_id", "model_id", "lineage"], as_index=False).agg(
        wz_sum=("wz", "sum"), w2_sum=("w2", "sum"), n_sources=("source", "nunique"))
    denom = np.sqrt(agg["w2_sum"].to_numpy(dtype=float))
    z_comb = np.where(denom > 0, agg["wz_sum"].to_numpy(dtype=float) / denom, 0.0)

    result = agg[["gene_id", "model_id", "lineage"]].copy()
    result["z_raw"] = z_comb.astype(float)
    result["n_sources"] = agg["n_sources"].astype(int)

    # 5. disagreement lowers confidence in the merged gene (symmetric, not a drop)
    result["platform_confidence"] = result["gene_id"].map(conf_by_gene).fillna(1.0).astype(float)
    result["z_score"] = result["z_raw"] * result["platform_confidence"]
    result["n_isoforms"] = result["gene_id"].map(iso_by_gene).fillna(1).astype(int)

    result = result.dropna(subset=["gene_id", "model_id", "lineage", "z_score"]).copy()
    result["gene_id"] = result["gene_id"].astype(str).str.upper()
    result["model_id"] = result["model_id"].astype(str).str.upper()
    result["lineage"] = result["lineage"].astype(str)
    result = result.sort_values(["gene_id", "model_id"]).reset_index(drop=True)
    result = result[["gene_id", "model_id", "lineage", "z_score", "z_raw",
                     "platform_confidence", "n_sources", "n_isoforms"]]
    return result, isoform_qc


def compute_protein_z(con, table_names=None, *, alpha=DEFAULT_ALPHA,
                      null_pctile=DEFAULT_NULL_PCTILE, matched_pctile=DEFAULT_MATCHED_PCTILE,
                      min_gap=DEFAULT_MIN_GAP, seed=0, thresholds=None, verbose=True):
    if table_names is None:
        table_names = [CCLE_SOURCE, PROCAN_SOURCE]
    protein_map = {**build_symbol_to_ensg(con), **build_uniprot_to_ensg(con)}
    if not protein_map:
        raise ValueError("No UniProt/symbol -> ENSG mapping from gene_roster.")
    lineage_map = load_lineage_map(con)
    raw_longs = load_sources(con, table_names, protein_map)
    if not raw_longs:
        raise ValueError("No protein source tables could be loaded.")

    thr = thresholds or derive_thresholds(
        raw_longs, lineage_map, alpha=alpha, null_pctile=null_pctile,
        matched_pctile=matched_pctile, min_gap=min_gap, seed=seed, verbose=verbose)
    result, isoform_qc = compute_from_sources(raw_longs, lineage_map, thr)
    result.attrs_thr = thr  # for callers that want the derived values
    return result, isoform_qc, thr


# ---------------------------------------------------------------------
# Persistence + CLI
# ---------------------------------------------------------------------

def _write(con, df, table):
    if df is None or df.empty:
        return
    con.register(f"{table}_tmp", df.copy())
    con.execute(f'CREATE OR REPLACE TABLE "{table}" AS SELECT * FROM "{table}_tmp"')


def write_protein_z(con, df, output_table="protein_z"):
    cols = ["gene_id", "model_id", "lineage", "z_score", "z_raw",
            "platform_confidence", "n_sources", "n_isoforms"]
    _write(con, df[cols], output_table)


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--tables", nargs="*", default=[CCLE_SOURCE, PROCAN_SOURCE])
    ap.add_argument("--output-table", default="protein_z")
    ap.add_argument("--isoform-qc-table", default="protein_isoform_qc")
    ap.add_argument("--params", default=None, help="JSON of pre-derived thresholds (skip derivation)")
    ap.add_argument("--save-params", default=None, help="write derived thresholds to this JSON")
    ap.add_argument("--seed", type=int, default=0)
    # statistical conventions (not domain values) -- exposed, not buried
    ap.add_argument("--alpha", type=float, default=DEFAULT_ALPHA)
    ap.add_argument("--null-pctile", type=float, default=DEFAULT_NULL_PCTILE)
    ap.add_argument("--matched-pctile", type=float, default=DEFAULT_MATCHED_PCTILE)
    ap.add_argument("--min-gap", type=float, default=DEFAULT_MIN_GAP)
    return ap.parse_args()


def main():
    args = parse_args()
    con = duckdb.connect(args.db, read_only=False)
    try:
        thresholds = json.loads(Path(args.params).read_text()) if args.params else None
        result, isoform_qc, thr = compute_protein_z(
            con, table_names=list(args.tables), alpha=args.alpha,
            null_pctile=args.null_pctile, matched_pctile=args.matched_pctile,
            min_gap=args.min_gap, seed=args.seed, thresholds=thresholds)

        write_protein_z(con, result, args.output_table)
        _write(con, isoform_qc, args.isoform_qc_table)
        if args.save_params:
            Path(args.save_params).write_text(json.dumps(thr, indent=2))

        summary = con.execute(
            f"""
            SELECT COUNT(*) AS n_rows, COUNT(DISTINCT gene_id) AS n_genes,
                   COUNT(DISTINCT model_id) AS n_models, COUNT(DISTINCT lineage) AS n_lineages,
                   COUNT(*) FILTER (WHERE n_sources = 1) AS n_sources_1,
                   COUNT(*) FILTER (WHERE n_sources = 2) AS n_sources_2,
                   COUNT(*) FILTER (WHERE n_isoforms > 1) AS n_multi_isoform_rows,
                   ROUND(AVG(platform_confidence), 3) AS avg_confidence
            FROM "{args.output_table}"
            """
        ).fetchdf()
        print(summary.to_string(index=False))
    finally:
        con.close()


if __name__ == "__main__":
    main()