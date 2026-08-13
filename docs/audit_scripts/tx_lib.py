"""Harness: the notebook's own functions, extracted verbatim, pointed at the real DB.

Source: C:\\Users\\vigne\\Downloads\\02_transcriptonomics.ipynb
Nothing here is production code -- it exists so the notebook's behaviour can be
measured and so corrected variants can be compared against it like for like.
"""
import numpy as np, pandas as pd, duckdb
from scipy import stats

DB_PATH = r"C:\Disertation\UoB-GeneTraceAI-25-26\src\pipeline\outputs\celllineselector.db"

# ---- constants, cell 1 verbatim ----
MIN_PEERS_FOR_LINEAGE = 15
EXPRESSED_MIN  = 1.0
PRIOR_STRENGTH = 1.0
I2_CUT         = 50
MAD_FLOOR       = 0.1
MIN_PEERS_FOR_Z = 3
SILENT_FRAC     = 0.20
HPA_FOLD   = 5.0
DELTA_CUT  = 0.5
TAU_BROAD, TAU_SPECIFIC = 0.30, 0.60
NO_LINEAGE      = "(no lineage recorded)"
SILENT_SOURCES  = ["depmap_expr", "hpa_rna"]
FLOOR_FRAC_LOW  = 0.02
TAU_SOURCES = ["depmap_expr", "hpa_rna"]


def connect():
    return duckdb.connect(DB_PATH, read_only=True)


def cols_of(con, t):
    return [r[0] for r in con.execute(f'DESCRIBE "{t}"').fetchall()]


def table_exists(con, t):
    return con.execute("SELECT count(*) FROM information_schema.tables "
                       "WHERE lower(table_name)=?", [t.lower()]).fetchone()[0] > 0


def has_col(con, t, c):
    return table_exists(con, t) and c in cols_of(con, t)


def resolve_gene(con, gene):
    g = str(gene).strip().lower()
    if g.startswith("ensg"):
        q = "SELECT gene_id, gene_names, uniprot_ids, hugo_symbol FROM gene WHERE gene_id = ?"
    else:
        q = ("SELECT gene_id, gene_names, uniprot_ids, hugo_symbol FROM gene "
             "WHERE list_contains(list_transform(gene_names, x -> lower(x)), ?)")
    row = con.execute(q, [g]).fetchone()
    if row is None:
        raise ValueError(f"'{gene}' not found in gene table (as id or name)")
    return row[0], list(row[1] or []), list(row[2] or []), row[3]


def lineage_map(con, table="sample_info", col="lineage"):
    d = con.execute(f'SELECT model_id, "{col}" AS lineage FROM "{table}" '
                    f'WHERE "{col}" IS NOT NULL').df()
    return d.drop_duplicates("model_id").set_index("model_id")["lineage"]


_LIN_CACHE = {}


def fetch_gene_expression(con, gene, log_geo=True, log_hpa=True, verbose=True,
                          geo_by_series=False):
    """Notebook cell 3, verbatim, plus an OPTIONAL geo_by_series switch used only
    by A3 to expose the per-series structure the notebook collapses away."""
    gid, gnames, _, hugo = resolve_gene(con, gene)
    frames = []

    if has_col(con, "depmap_expr", gid):
        d = con.execute(f'SELECT model_id, "{gid}" AS value FROM depmap_expr '
                        f'WHERE "{gid}" IS NOT NULL').df()
        d["n_samples"] = 1; d["source"] = "depmap_expr"; frames.append(d)

    if has_col(con, "geo_expr", gid):
        if geo_by_series:
            d = con.execute(
                f'SELECT gi.gse_id AS gse_id, ge.model_id AS model_id, '
                f'median(ge."{gid}") AS value, count(*) AS n_samples '
                f'FROM geo_expr ge LEFT JOIN geo_info gi ON lower(ge.sample)=lower(gi.geo_accession) '
                f'WHERE ge."{gid}" IS NOT NULL GROUP BY 1,2').df()
            if log_geo:
                d["value"] = np.log2(d["value"].clip(lower=0) + 1)
            d["gse_id"] = d["gse_id"].fillna("(no gse_id)")
            d["source"] = "geo_expr:" + d["gse_id"]
            frames.append(d.drop(columns=["gse_id"]))
        else:
            d = con.execute(f'SELECT model_id, median("{gid}") AS value, count(*) AS n_samples '
                            f'FROM geo_expr WHERE "{gid}" IS NOT NULL GROUP BY model_id').df()
            if log_geo:
                d["value"] = np.log2(d["value"].clip(lower=0) + 1)
            d["source"] = "geo_expr"; frames.append(d)

    vcol = next((c for c in ("ntpm", "ptpm", "tpm") if has_col(con, "hpa_rna", c)), None)
    if vcol:
        d = con.execute(f'SELECT model_id, avg("{vcol}") AS value, count(*) AS n_samples '
                        f'FROM hpa_rna WHERE lower(gene)=? GROUP BY model_id', [gid]).df()
        if log_hpa:
            d["value"] = np.log2(d["value"].clip(lower=0) + 1)
        d["source"] = "hpa_rna"; frames.append(d)

    df = pd.concat(frames, ignore_index=True)
    if "lm" not in _LIN_CACHE:
        _LIN_CACHE["lm"] = lineage_map(con)
    df["lineage"] = df["model_id"].map(_LIN_CACHE["lm"])
    if verbose:
        print(f"{gene} -> {gid} ({hugo or ', '.join(gnames)})")
        print(df.groupby("source").agg(n_lines=("model_id", "nunique"),
              median=("value", "median"), min=("value", "min"),
              max=("value", "max")).round(2).to_string())
    return df, gid


# ---------------- cell 8 verbatim ----------------
def robust_z(x, v, mad_floor=MAD_FLOOR, min_n=MIN_PEERS_FOR_Z,
             expressed_min=EXPRESSED_MIN, silent_frac=SILENT_FRAC):
    v = np.asarray(v, float)
    if len(v) < min_n:
        return np.nan, "too few peers"
    if (v > expressed_min).mean() < silent_frac:
        return np.nan, "gene silent in lineage"
    mad = stats.median_abs_deviation(v, scale="normal")
    if mad <= 0:
        return np.nan, "no spread"
    return (x - np.median(v)) / max(mad, mad_floor), "ok"


def empirical_p(x, v, tail="two"):
    v = np.asarray(v); n = len(v)
    hi = ((v >= x).sum() + 1) / (n + 1); lo = ((v <= x).sum() + 1) / (n + 1)
    return min(1.0, 2 * min(hi, lo)) if tail == "two" else (hi if tail == "hi" else lo)


def stouffer(zs, weights=None):
    z = np.asarray([v for v in zs if np.isfinite(v)], float)
    if not len(z):
        return np.nan, np.nan, 0
    w = np.ones(len(z)) if weights is None else np.asarray(
        [w for w, v in zip(weights, zs) if np.isfinite(v)], float)
    Z = (w * z).sum() / np.sqrt((w ** 2).sum())
    return Z, 2 * (1 - stats.norm.cdf(abs(Z))), len(z)


def stouffer_corr(zs, weights=None, R=None):
    """CORRECTED combination -- Strube (1985) / Hartung (1999) generalised Stouffer.
    Denominator becomes sqrt(w' R w) instead of sqrt(w'w)."""
    keep = [i for i, v in enumerate(zs) if np.isfinite(v)]
    z = np.asarray([zs[i] for i in keep], float)
    if not len(z):
        return np.nan, np.nan, 0
    w = (np.ones(len(z)) if weights is None
         else np.asarray([weights[i] for i in keep], float))
    if R is None:
        Rk = np.eye(len(z))
    else:
        Rk = np.asarray(R)[np.ix_(keep, keep)]
    denom = np.sqrt(max(w @ Rk @ w, 1e-12))
    Z = (w * z).sum() / denom
    return Z, 2 * (1 - stats.norm.cdf(abs(Z))), len(z)


def heterogeneity(zs):
    z = np.asarray([v for v in zs if np.isfinite(v)], float); k = len(z)
    if k < 2:
        return np.nan, np.nan, np.nan
    Q = ((z - z.mean()) ** 2).sum()
    return Q, 1 - stats.chi2.cdf(Q, k - 1), (max(0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0)


def shrink(z, k, prior_strength=1.0):
    return z * k / (k + prior_strength) if np.isfinite(z) else np.nan


def score_gene_lineage(df, min_peers=MIN_PEERS_FOR_LINEAGE, expressed_min=EXPRESSED_MIN,
                       weight_by_samples=True, prior_strength=1.0,
                       replicate_warn=1.0, silent_sources=None, verbose=True,
                       R=None, use_empirical_z=False):
    """Notebook cell 8, verbatim, plus two OPTIONAL switches used only to measure
    the corrected variants:
        R                -- source correlation matrix for stouffer_corr (C1/A4)
        use_empirical_z  -- combine z=Phi^-1(1-p_emp/2)*sign instead of robust z (C1)
    With R=None and use_empirical_z=False the behaviour is bit-identical to the
    notebook."""
    sil_src = silent_sources if silent_sources is not None else SILENT_SOURCES
    d = df.dropna(subset=["value"]).copy()
    d["lineage"] = d["lineage"].fillna(NO_LINEAGE)
    if "n_samples" not in d.columns:
        d["n_samples"] = 1

    d = (d.groupby(["model_id", "source", "lineage"], as_index=False)
           .agg(value=("value", "median"), n_samples=("n_samples", "sum")))

    recs = []
    for (src, lin), g in d.groupby(["source", "lineage"]):
        v = g["value"].to_numpy()
        expressed = v[v > expressed_min]
        small = len(v) < min_peers
        for _, r in g.iterrows():
            z, why = robust_z(r.value, v, expressed_min=expressed_min)
            recs.append({"model_id": r.model_id, "source": src, "lineage": lin,
                         "value": r.value, "n_peers": len(v), "small_lineage": small,
                         "lineage_median": np.median(v), "z": z, "z_reason": why,
                         "pct_expressed": (np.mean(expressed < r.value) * 100
                                           if len(expressed) else np.nan),
                         "p_emp": empirical_p(r.value, v),
                         "n_samples": r.n_samples,
                         "frac_expressed": float((v > expressed_min).mean())})
    per = pd.DataFrame(recs)

    if use_empirical_z:
        # C1: convert the rank-based empirical p into a standardised statistic,
        # signed by the direction of the deviation.
        sgn = np.sign(per["value"] - per["lineage_median"]).replace(0, 1.0)
        pe = per["p_emp"].clip(lower=1e-12, upper=1 - 1e-12)
        per["z_use"] = sgn * stats.norm.isf(pe / 2.0)
        per.loc[~np.isfinite(per["z"]), "z_use"] = np.nan   # keep the same guards
    else:
        per["z_use"] = per["z"]

    src_order = sorted(per.source.unique())
    idx = {s: i for i, s in enumerate(src_order)}

    out = []
    for mid, g in per.groupby("model_id"):
        if R is None:
            zs = g.z_use.tolist()
            w = np.sqrt(g.n_samples.astype(float)).tolist() if weight_by_samples else None
            Z, p, k = stouffer(zs, w)
        else:
            full = [np.nan] * len(src_order)
            fw = [1.0] * len(src_order)
            for _, r in g.iterrows():
                full[idx[r.source]] = r.z_use
                fw[idx[r.source]] = np.sqrt(float(r.n_samples)) if weight_by_samples else 1.0
            Z, p, k = stouffer_corr(full, fw, R)
        Q, pq, I2 = heterogeneity(g.z_use.tolist())
        ps = g.p_emp.dropna().to_numpy()
        fisher_p = (stats.combine_pvalues(ps, method="fisher")[1] if len(ps) > 1
                    else (ps[0] if len(ps) else np.nan))
        gs = g[g.source.isin(sil_src)]
        out.append({"model_id": mid, "lineage": g.lineage.iloc[0],
                    "z_reasons": "; ".join(sorted(set(g.z_reason))),
                    "lineage_silent": bool(len(gs) and (gs.frac_expressed < SILENT_FRAC).all()),
                    "n_sources": k, "n_peers": int(g.n_peers.max()),
                    "n_profiles": int(g.n_samples.sum()),
                    "small_lineage": bool(g.small_lineage.any()),
                    "mean_value": round(g.value.mean(), 2),
                    "z_combined": Z, "z_shrunk": shrink(Z, k, prior_strength),
                    "p_combined": p, "fisher_p": fisher_p, "I2": I2,
                    "sources_agree": (not np.isfinite(I2)) or I2 < I2_CUT,
                    "frac_lineage_expressed": (round(gs.frac_expressed.min(), 3)
                                               if len(gs) else np.nan)})
    res = pd.DataFrame(out)

    ok = res.p_combined.notna()
    res.loc[ok, "q_value"] = stats.false_discovery_control(res.loc[ok, "p_combined"])

    def verdict(r):
        if not np.isfinite(r.z_shrunk):   return f"Insufficient Data ({r.z_reasons})"
        if not r.sources_agree:           return "Conflicting Sources"
        strong = (r.q_value < .05) and abs(r.z_shrunk) > 1.5
        d_ = "HIGH" if r.z_shrunk > 0 else "LOW"
        if strong and r.n_sources >= 2 and not r.small_lineage: return d_
        if strong and r.n_sources < 2:    return f"{d_} (Single Source)"
        if strong and r.small_lineage:    return f"{d_} (Small Lineage)"
        if r.z_shrunk >  1.0: return "High (Weak)"
        if r.z_shrunk < -1.0: return "Low (weak)"
        return "Typical"

    res["verdict"] = res.apply(verdict, axis=1)
    return res.sort_values("z_shrunk", ascending=False).reset_index(drop=True), per
