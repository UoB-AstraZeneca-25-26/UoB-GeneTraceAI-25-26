"""
Measure every number the report figures need, straight from the production
parquet outputs, and cache to outputs/report_figures/_stats.json.

Nothing here is hard-coded from a doc. If a figure quotes a number, that
number was computed by this script on this run, and `_stats.json` records
the provenance so it can be re-checked.

Run:  python src/figures/compute_stats.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import spearmanr, mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (ROOT, PIPE, REF, DATA_CLEAN, TRACK, VALID, OUT, STATS,
                    load_universe, load_gdsc_truth, hit_at_k)

RNG = np.random.default_rng(42)
S = {}


def section(name):
    print(f"\n{'=' * 62}\n{name}\n{'=' * 62}")


# ─────────────────────────────────────────────────────────────────────────
# 1. Core entity counts and the n_layers distribution
# ─────────────────────────────────────────────────────────────────────────
def counts():
    section("1. Entity counts / n_layers distribution")
    pf = pq.ParquetFile(PIPE / "core_score.parquet")
    n_rows = pf.metadata.num_rows

    nl = pd.read_parquet(PIPE / "core_score.parquet", columns=["n_layers"])["n_layers"]
    dist = nl.value_counts().sort_index()

    cov_gene = pd.read_parquet(ROOT / "outputs/coverage_tables/coverage_per_gene.parquet")
    cov_line = pd.read_parquet(ROOT / "outputs/coverage_tables/coverage_per_cell_line.parquet")

    lookup_lines = pq.ParquetFile(REF / "cell_line_lookup.parquet").metadata.num_rows
    lookup_genes = pq.ParquetFile(REF / "gene_lookup.parquet").metadata.num_rows
    universe_genes = len(load_universe())
    harm_rows = pq.ParquetFile(PIPE / "harmonised_enriched.parquet").metadata.num_rows

    S["counts"] = {
        "core_score_rows": int(n_rows),
        "scored_genes": int(len(cov_gene)),
        "scored_cell_lines": int(len(cov_line)),
        "lookup_cell_lines": int(lookup_lines),
        "lookup_genes": int(lookup_genes),
        "universe_genes": int(universe_genes),
        "harmonised_rows": int(harm_rows),
        "n_layers_dist": {str(k): int(v) for k, v in dist.items()},
        "pct_2layer": float(100 * dist.get(2, 0) / n_rows),
        "flags_rows": int(pq.ParquetFile(PIPE / "flags_with_driver.parquet").metadata.num_rows),
        "cna_rows": int(pq.ParquetFile(PIPE / "cna_flags.parquet").metadata.num_rows),
        "curated_genes": int(pq.ParquetFile(PIPE / "gene_regime.parquet").metadata.num_rows),
    }
    print(json.dumps(S["counts"], indent=2))
    del nl


# ─────────────────────────────────────────────────────────────────────────
# 2. Which modalities each cell line actually has (drives the UpSet figure)
# ─────────────────────────────────────────────────────────────────────────
def modality_presence():
    section("2. Modality presence per cell line")
    lookup = pd.read_parquet(REF / "cell_line_lookup.parquet")[["model_id"]]
    lookup["model_id"] = lookup["model_id"].str.lower()
    universe = set(lookup["model_id"])

    def ids(series):
        return set(series.dropna().astype(str).str.lower()) & universe

    prof = pd.read_parquet(REF / "depmap_profiles.parquet")
    sets = {}
    sets["RNA"] = ids(prof[prof["datatype"] == "rna"]["modelid"])
    sets["Protein"] = ids(pd.read_parquet(TRACK / "proteomics.parquet")["depmap_id"])
    sets["Mutation"] = ids(pd.read_parquet(
        TRACK / "mutations_collapsed.parquet", columns=["model_id"])["model_id"])
    sets["Fusion"] = ids(pd.read_parquet(
        TRACK / "fusions_gene_level.parquet", columns=["model_id"])["model_id"])
    sets["CNA"] = ids(pd.read_parquet(
        PIPE / "cna_flags.parquet", columns=["model_id"])["model_id"])
    sets["Signatures"] = ids(pd.read_parquet(
        TRACK / "signatures_model_level.parquet", columns=["model_id"])["model_id"])
    sets["Metabolomics"] = ids(pd.read_parquet(TRACK / "metabolomics.parquet",
                                               columns=["model_id"])["model_id"])
    sets["miRNA"] = ids(pd.read_parquet(TRACK / "mirna_model_level.parquet",
                                        columns=["model_id"])["model_id"])
    bridge = pd.read_parquet(VALID / "id_bridge.parquet")
    sets["Chronos"] = ids(bridge["model_id"])
    sets["GDSC drug"] = ids(pd.read_parquet(VALID / "gdsc_ach.parquet",
                                            columns=["model_id"])["model_id"])

    order = list(sets)
    mat = pd.DataFrame({k: lookup["model_id"].isin(v).values for k, v in sets.items()},
                       index=lookup["model_id"].values)

    # Combination counts, biggest first -- the UpSet bar series.
    combo = (mat.apply(lambda r: tuple(sorted(c for c in order if r[c])), axis=1)
                .value_counts())
    combos = [{"layers": list(k), "n": int(v)} for k, v in combo.head(22).items()]

    S["modality"] = {
        "order": order,
        "set_sizes": {k: int(len(v)) for k, v in sets.items()},
        "universe": int(len(universe)),
        "combinations": combos,
        "n_modalities_hist": {str(k): int(v) for k, v in
                              mat.sum(axis=1).value_counts().sort_index().items()},
    }
    print("set sizes:", S["modality"]["set_sizes"])
    print("top combos:", combos[:6])


# ─────────────────────────────────────────────────────────────────────────
# 3. Gene-axis reach: which sources can answer a query keyed on ensg_id
# ─────────────────────────────────────────────────────────────────────────
def gene_axis_reach():
    section("3. Gene-axis reach per source")
    universe = set(load_universe()["ensg_id"])

    # (label, path, how the gene axis is obtained, native key)
    specs = [
        ("Expression (DepMap)", DATA_CLEAN / "depmap_expr_clean.parquet", "columns", "ENSG (native)"),
        ("Proteomics (CCLE/Gygi)", TRACK / "proteomics.parquet", "uniprot_cols", "UniProt -> ENSG"),
        ("Mutations", TRACK / "mutations_collapsed.parquet", "ensg_col", "ENSG (native)"),
        ("Fusions", TRACK / "fusions_gene_level.parquet", "ensg_col", "ENSG (native)"),
        ("CNA flags", PIPE / "cna_flags.parquet", "ensg_col", "ENSG (native)"),
        ("Chronos essentiality", VALID / "chronos_long.parquet", "ensg_col", "ENSG (native)"),
        ("GDSC drug targets", VALID / "gdsc_pairs.parquet", "target_ensg", "symbol -> ENSG"),
        ("HPA RNA", DATA_CLEAN / "hpa_rna_clean.parquet", "probe", "symbol/ENSG"),
        ("GEO expression", DATA_CLEAN / "geo_expr_clean.parquet", "probe", "probe/symbol"),
        ("Metabolomics", TRACK / "metabolomics.parquet", "none", "model_id only"),
        ("miRNA", TRACK / "mirna_model_level.parquet", "none", "model_id only"),
        ("Signatures (MSI/CIN/ploidy)", TRACK / "signatures_model_level.parquet", "none", "model_id only"),
        ("Cellosaurus", TRACK / "Cellosaurus.parquet", "none", "cell-line identity"),
        ("DepMap sample_info", TRACK / "sample_info.parquet", "none", "cell-line identity"),
    ]

    rows = []
    for label, path, how, native in specs:
        genes = set()
        note = ""
        try:
            if how == "columns":
                cols = pq.ParquetFile(path).schema_arrow.names
                genes = {c.split(".")[0].lower() for c in cols
                         if c.lower().startswith("ensg")}
            elif how == "uniprot_cols":
                u = load_universe()
                up = (u[["ensg_id", "uniprot_ids"]].dropna(subset=["uniprot_ids"])
                      .assign(k=lambda x: x["uniprot_ids"].str.strip().str.lower()))
                m = dict(zip(up["k"], up["ensg_id"]))
                cols = pq.ParquetFile(path).schema_arrow.names
                genes = {m[c.lower()] for c in cols if c.lower() in m}
                note = "reachable only after the UniProt->ENSG bridge"
            elif how == "ensg_col":
                col = pd.read_parquet(path, columns=["ensg_id"])["ensg_id"]
                genes = set(col.dropna().astype(str).str.split(".").str[0].str.lower())
            elif how == "target_ensg":
                col = pd.read_parquet(path, columns=["target_ensg"])["target_ensg"]
                genes = set(col.dropna().astype(str).str.lower())
                note = "curated drug->target map only"
            elif how == "probe":
                cols = pq.ParquetFile(path).schema_arrow.names
                low = [c.lower() for c in cols]
                if "ensg_id" in low:
                    col = pd.read_parquet(path, columns=[cols[low.index("ensg_id")]]).iloc[:, 0]
                    genes = set(col.dropna().astype(str).str.split(".").str[0].str.lower())
                else:
                    genes = {c.split(".")[0].lower() for c in cols if c.lower().startswith("ensg")}
                    note = "no ensg_id column -- symbol/probe axis only"
            elif how == "none":
                genes = set()
                note = "no gene axis: cannot answer a per-gene query at all"
        except Exception as exc:                      # noqa: BLE001
            note = f"unreadable: {exc}"

        inuni = genes & universe
        rows.append({"source": label, "native_key": native,
                     "genes_any": int(len(genes)),
                     "genes_in_universe": int(len(inuni)),
                     "note": note})
        print(f"  {label:<30} {len(inuni):>7,} / {len(universe):,}  {note}")

    S["gene_reach"] = {"universe": int(len(universe)), "sources": rows}
    S["gene_reach"]["n_with_reach"] = int(sum(r["genes_in_universe"] > 0 for r in rows))
    S["gene_reach"]["n_total"] = len(rows)


# ─────────────────────────────────────────────────────────────────────────
# 4. Inter-layer correlation -- the evidence for the correlation penalty
# ─────────────────────────────────────────────────────────────────────────
def interlayer_correlation():
    section("4. Inter-layer correlation matrix")
    epc = pd.read_parquet(VALID / "expr_protein_correlation_per_gene.parquet")
    # pre-v2 validation artefact, still keyed on UPPERCASE ENSG. Its index is
    # intersected with chron["ensg_id"] below, which IS lowercased -- leaving it
    # uppercase makes that intersection empty, and rho_EC / rho_PC come out null
    # over 0 genes rather than raising.
    epc.index = epc.index.astype(str).str.split(".").str[0].str.lower()
    rho_ep = float(epc["rho"].median())

    # Expression vs Chronos, and Protein vs Chronos, per gene on shared pairs.
    chron = pd.read_parquet(VALID / "chronos_long.parquet")
    bridge = pd.read_parquet(VALID / "id_bridge.parquet")
    chron = chron.merge(bridge, on="sanger_model_id", how="inner")
    chron["model_id"] = chron["model_id"].str.lower()
    chron["ensg_id"] = chron["ensg_id"].astype(str).str.split(".").str[0].str.lower()

    # Sample genes to keep this tractable; report n.
    genes = sorted(set(epc.index) & set(chron["ensg_id"].unique()))
    sample = list(RNG.choice(genes, size=min(600, len(genes)), replace=False))
    chron_s = chron[chron["ensg_id"].isin(set(sample))]
    chron_p = chron_s.pivot_table(index="model_id", columns="ensg_id",
                                  values="essentiality", aggfunc="mean")
    del chron, chron_s

    from common import load_expr_wide, load_prot_wide
    expr = load_expr_wide(genes=sample)
    prot = load_prot_wide(genes=sample)

    def per_gene_rho(A, B):
        shared_g = [g for g in A.columns if g in B.columns]
        shared_m = A.index.intersection(B.index)
        out = []
        for g in shared_g:
            a = A.loc[shared_m, g]
            b = B.loc[shared_m, g]
            ok = a.notna() & b.notna()
            if ok.sum() >= 30:
                out.append(spearmanr(a[ok], b[ok]).statistic)
        return np.array([x for x in out if np.isfinite(x)])

    r_ec = per_gene_rho(expr, chron_p)
    r_pc = per_gene_rho(prot, chron_p)

    S["interlayer"] = {
        "rho_EP_median": rho_ep,
        "rho_EP_mean": float(epc["rho"].mean()),
        "rho_EP_n_genes": int(len(epc)),
        "rho_EP_q25": float(epc["rho"].quantile(0.25)),
        "rho_EP_q75": float(epc["rho"].quantile(0.75)),
        "rho_EP_pct_below_0": float(100 * (epc["rho"] < 0).mean()),
        "rho_EC_median": float(np.median(r_ec)) if len(r_ec) else None,
        "rho_EC_n_genes": int(len(r_ec)),
        "rho_PC_median": float(np.median(r_pc)) if len(r_pc) else None,
        "rho_PC_n_genes": int(len(r_pc)),
        "rho_EP_hist": np.histogram(epc["rho"].values, bins=40, range=(-1, 1))[0].tolist(),
    }
    print(json.dumps({k: v for k, v in S["interlayer"].items() if k != "rho_EP_hist"},
                     indent=2))

    # Correlation-penalised weights that follow from these numbers.
    r = {"EP": abs(rho_ep),
         "EC": abs(S["interlayer"]["rho_EC_median"] or 0.0),
         "PC": abs(S["interlayer"]["rho_PC_median"] or 0.0)}
    mean_rho = {"Expression": (r["EP"] + r["EC"]) / 2,
                "Protein": (r["EP"] + r["PC"]) / 2,
                "Chronos": (r["EC"] + r["PC"]) / 2}
    raw = {k: 1 / (1 + v) for k, v in mean_rho.items()}
    tot = sum(raw.values())
    S["interlayer"]["weights_three_layer"] = {k: v / tot for k, v in raw.items()}
    print("implied 3-layer weights:", S["interlayer"]["weights_three_layer"])


# ─────────────────────────────────────────────────────────────────────────
# 5. Normalisation comparison: raw vs percentile vs min-max
# ─────────────────────────────────────────────────────────────────────────
def normalisation():
    section("5. Normalisation comparison")
    from common import load_expr_wide, load_prot_wide
    universe = sorted(load_universe()["ensg_id"])
    sample = list(RNG.choice(universe, size=300, replace=False))

    expr = load_expr_wide(genes=sample)
    prot = load_prot_wide(genes=sample)

    def summarise(df, label):
        raw = df.values.ravel()
        raw = raw[np.isfinite(raw)]
        pct = df.rank(axis=0, pct=True, na_option="keep").values.ravel()
        pct = pct[np.isfinite(pct)]
        mn, mx = df.min(axis=0), df.max(axis=0)
        span = (mx - mn).replace(0, np.nan)
        mm = ((df - mn) / span).values.ravel()
        mm = mm[np.isfinite(mm)]
        return {
            "label": label,
            "n": int(len(raw)),
            "raw_hist": np.histogram(raw, bins=60)[0].tolist(),
            "raw_edges": np.histogram(raw, bins=60)[1].tolist(),
            "raw_min": float(raw.min()), "raw_max": float(raw.max()),
            "pct_hist": np.histogram(pct, bins=60, range=(0, 1))[0].tolist(),
            "mm_hist": np.histogram(mm, bins=60, range=(0, 1))[0].tolist(),
            "mm_median": float(np.median(mm)),
            "mm_iqr": float(np.percentile(mm, 75) - np.percentile(mm, 25)),
            "pct_iqr": float(np.percentile(pct, 75) - np.percentile(pct, 25)),
        }

    S["normalisation"] = {"expr": summarise(expr, "Expression (log2 TPM+1)"),
                          "prot": summarise(prot, "Proteomics (log2 ratio)"),
                          "n_genes_sampled": len(sample)}
    for k in ("expr", "prot"):
        d = S["normalisation"][k]
        print(f"  {d['label']}: raw [{d['raw_min']:.2f}, {d['raw_max']:.2f}]  "
              f"minmax IQR={d['mm_iqr']:.3f}  pct IQR={d['pct_iqr']:.3f}")


# ─────────────────────────────────────────────────────────────────────────
# 6. core_score vs n_layers -- the honesty figure
# ─────────────────────────────────────────────────────────────────────────
def score_vs_layers():
    section("6. core_score / stratum_rank vs n_layers")
    df = pd.read_parquet(PIPE / "core_score.parquet",
                         columns=["core_score", "n_layers", "stratum_rank"])
    out = {}
    for col in ("core_score", "stratum_rank"):
        out[col] = {}
        for nl, grp in df.groupby("n_layers"):
            v = grp[col].dropna().values
            if len(v) > 400_000:
                v = RNG.choice(v, 400_000, replace=False)
            out[col][str(int(nl))] = {
                "n": int(len(grp)),
                "mean": float(np.mean(v)),
                "quantiles": [float(np.percentile(v, q))
                              for q in (5, 25, 50, 75, 95)],
                "sample": v[:8000].round(4).tolist(),
            }
    S["score_vs_layers"] = out
    for col in out:
        for nl, d in out[col].items():
            print(f"  {col} n_layers={nl}: n={d['n']:,} mean={d['mean']:.4f} "
                  f"median={d['quantiles'][2]:.4f}")
    del df


# ─────────────────────────────────────────────────────────────────────────
# 7. TSG mechanism check: PTEN homozygous-deleted vs diploid
# ─────────────────────────────────────────────────────────────────────────
def tsg_deletion():
    section("7. TSG deletion vs diploid score distributions")
    gl = pd.read_parquet(REF / "gene_lookup.parquet")[["ensg_id", "hgnc_symbol", "gene_role"]]
    gl["ensg_id"] = gl["ensg_id"].astype(str).str.split(".").str[0].str.lower()
    sym2ensg = dict(zip(gl["hgnc_symbol"], gl["ensg_id"]))

    cna = pd.read_parquet(PIPE / "cna_flags.parquet")
    cna["model_id"] = cna["model_id"].str.lower()

    results = {}
    for symbol in ["PTEN", "TP53", "CDKN2A", "RB1", "SMAD4", "STK11"]:
        ensg = sym2ensg.get(symbol)
        if ensg is None:
            continue
        sc = pd.read_parquet(PIPE / "core_score.parquet",
                             filters=[("ensg_id", "==", ensg)],
                             columns=["model_id", "ensg_id", "core_score", "n_layers"])
        if sc.empty:
            continue
        sc["model_id"] = sc["model_id"].str.lower()

        # cna_flags only holds pairs that have CNA context at all, so the
        # comparison arm is "every other scored line for this gene" --
        # lines with a retained copy, plus lines with no CNA call. Calling
        # that arm "diploid" would overclaim; it is "not homozygous-deleted".
        g = cna[cna["ensg_id"] == ensg][["model_id", "min_cn", "max_cn", "cna_type"]]
        m = sc.merge(g, on="model_id", how="left")

        is_del = m["min_cn"].notna() & (m["min_cn"] <= 0)
        deleted = m.loc[is_del, "core_score"].dropna().values
        retained = m.loc[~is_del, "core_score"].dropna().values
        no_call = m.loc[m["min_cn"].isna(), "core_score"].dropna().values

        entry = {
            "ensg_id": ensg,
            "gene_role": str(gl.loc[gl["ensg_id"] == ensg, "gene_role"].iloc[0]),
            "n_scored": int(len(m)),
            "n_with_cna_call": int(m["min_cn"].notna().sum()),
            "deleted_n": int(len(deleted)), "deleted_mean": _m(deleted),
            "retained_n": int(len(retained)), "retained_mean": _m(retained),
            "no_call_n": int(len(no_call)), "no_call_mean": _m(no_call),
            # The production table serves TSGs inverted (1 - core_score), so
            # report both scales -- the report quotes the served number.
            "deleted_mean_inverted": None if not len(deleted) else float(1 - np.mean(deleted)),
            "retained_mean_inverted": None if not len(retained) else float(1 - np.mean(retained)),
            "deleted_sample": deleted.round(4).tolist()[:2000],
            "retained_sample": retained.round(4).tolist()[:3000],
        }
        if len(deleted) >= 5 and len(retained) >= 5:
            u = mannwhitneyu(deleted, retained, alternative="two-sided")
            entry["mannwhitney_p"] = float(u.pvalue)
            entry["delta"] = entry["deleted_mean"] - entry["retained_mean"]
            entry["delta_inverted"] = (entry["deleted_mean_inverted"]
                                       - entry["retained_mean_inverted"])
        results[symbol] = entry
        print(f"  {symbol} ({ensg}): deleted n={len(deleted)} mean={_m(deleted):.4f} "
              f"(inverted {entry['deleted_mean_inverted']:.4f})  "
              f"retained n={len(retained)} mean={_m(retained):.4f} "
              f"(inverted {entry['retained_mean_inverted']:.4f})  "
              f"p={entry.get('mannwhitney_p')}")

    S["tsg_deletion"] = results


def _m(a):
    return float(np.mean(a)) if len(a) else None


# ─────────────────────────────────────────────────────────────────────────
# 8. Routing statistic distribution (gene-class routing)
# ─────────────────────────────────────────────────────────────────────────
def routing():
    section("8. Routing statistic distribution")
    reg = pd.read_parquet(PIPE / "gene_regime.parquet")
    n_lines = pd.read_parquet(ROOT / "outputs/coverage_tables/coverage_per_cell_line.parquet")
    threshold = 2 * (20 / len(n_lines))

    # How the full 28.4M-row table is actually routed (class + regime_source).
    pr = pd.read_parquet(PIPE / "predictions_with_confidence.parquet",
                         columns=["class", "regime_source"])
    class_rows = pr["class"].value_counts()
    source_rows = pr["regime_source"].value_counts()
    del pr

    gl = pd.read_parquet(REF / "gene_lookup.parquet")[["ensg_id", "hgnc_symbol"]]
    gl["ensg_id"] = gl["ensg_id"].astype(str).str.split(".").str[0].str.lower()
    reg = reg.merge(gl, on="ensg_id", how="left")

    S["routing"] = {
        "threshold": float(threshold),
        "n_scored_lines": int(len(n_lines)),
        "curated": [
            {"ensg_id": r.ensg_id, "symbol": r.hgnc_symbol,
             "stat": float(r.hit_at_20_train), "class": r["class"],
             "n_sensitive": int(r.n_sensitive_gdsc)}
            for _, r in reg.iterrows()
        ],
        "class_counts_curated": {k: int(v) for k, v in reg["class"].value_counts().items()},
        "class_rows_full": {str(k): int(v) for k, v in class_rows.items()},
        "regime_source_rows_full": {str(k): int(v) for k, v in source_rows.items()},
    }
    print("curated class counts:", S["routing"]["class_counts_curated"])
    print("full-table class rows:", S["routing"]["class_rows_full"])
    print("regime_source rows:", S["routing"]["regime_source_rows_full"])
    print(f"threshold (2x random baseline) = {threshold:.4f}")


# ─────────────────────────────────────────────────────────────────────────
# 9. Recipe variants: mean rho with bootstrap CI (paired over drug-targets)
# ─────────────────────────────────────────────────────────────────────────
def recipes():
    section("9. Recipe variants -- rho with bootstrap CI")
    d = pd.read_parquet(VALID / "phase3_fair4way_results.parquet")
    g = d[(d["arm"] == "gdsc") & d["spearman_rho"].notna()]

    wide = g.pivot_table(index=["drug_id", "target"], columns="variant",
                         values="spearman_rho")
    wide = wide.dropna()
    variants = list(wide.columns)
    n = len(wide)
    print(f"  paired over {n} drug-target pairs, {len(variants)} variants")

    B = 5000
    idx = RNG.integers(0, n, size=(B, n))
    arr = wide.values
    boot = arr[idx].mean(axis=1)          # (B, n_variants)

    res = []
    for i, v in enumerate(variants):
        pt = float(arr[:, i].mean())
        lo, hi = np.percentile(boot[:, i], [2.5, 97.5])
        res.append({"variant": v, "rho": pt, "lo": float(lo), "hi": float(hi),
                    "median_rho": float(np.median(arr[:, i])),
                    "excludes_zero": bool(lo > 0 or hi < 0)})
        print(f"  {v:<32} rho={pt:+.4f}  95% CI [{lo:+.4f}, {hi:+.4f}]"
              f"{'  EXCLUDES 0' if lo > 0 or hi < 0 else ''}")

    # Paired delta of each variant against the production recipe.
    prod = "corr_penalized_weighted_sum"
    pi = variants.index(prod)
    deltas = []
    for i, v in enumerate(variants):
        if v == prod:
            continue
        d_arr = arr[:, pi] - arr[:, i]
        bl, bh = np.percentile(boot[:, pi] - boot[:, i], [2.5, 97.5])
        deltas.append({"vs": v, "delta": float(d_arr.mean()),
                       "lo": float(bl), "hi": float(bh),
                       "significant": bool(bl > 0 or bh < 0)})

    S["recipes"] = {"n_pairs": int(n), "variants": res,
                    "paired_vs_production": deltas,
                    "production": prod, "n_boot": B}

    # The full 12-way normalisation x combination grid (median rho only --
    # per-pair values were not retained for these).
    dm = pd.read_parquet(VALID / "phase3_decision_matrix_full.parquet")
    S["recipes"]["grid12"] = dm[dm["arm"] == "gdsc"].to_dict("records")


# ─────────────────────────────────────────────────────────────────────────
# 10. Per-gene hit@20 on the held-out split (consistent denominator)
# ─────────────────────────────────────────────────────────────────────────
def per_gene_hits():
    section("10. Per-gene hit@20, held-out genes")
    ev = pd.read_parquet(PIPE / "stage4_eval_consistent_denom.parquet")
    gl = pd.read_parquet(REF / "gene_lookup.parquet")[["ensg_id", "hgnc_symbol"]]
    gl["ensg_id"] = gl["ensg_id"].astype(str).str.split(".").str[0].str.lower()
    ev = ev.merge(gl, on="ensg_id", how="left")

    with open(PIPE / "stage4_eval_summary.json", encoding="utf-8") as fh:
        summary = json.load(fh)

    S["per_gene_hits"] = {
        "genes": ev[["ensg_id", "hgnc_symbol", "class", "known",
                     "hits_flat", "hits_driver", "hits_any",
                     "hr_flat", "hr_driver", "hr_any"]].to_dict("records"),
        "summary": summary["overall"],
        "per_class": summary["per_class"],
        "variance_reduction": summary["variance_reduction"],
        "denominator": summary["denominator"],
        "random_baseline": float(20 / len(pd.read_parquet(
            ROOT / "outputs/coverage_tables/coverage_per_cell_line.parquet"))),
    }
    print(f"  {len(ev)} held-out genes; denominators {ev['known'].min()}-{ev['known'].max()}")
    print(f"  random baseline = {S['per_gene_hits']['random_baseline']:.4f}")


# ─────────────────────────────────────────────────────────────────────────
# 11. Does the confidence tier predict a better outcome?
# ─────────────────────────────────────────────────────────────────────────
def confidence_vs_outcome():
    section("11. Confidence tier vs observed GDSC hit rate")
    truth = load_gdsc_truth()
    genes = sorted(truth)
    print(f"  {len(genes)} genes with GDSC ground truth")

    pred = pd.read_parquet(
        PIPE / "predictions_with_confidence.parquet",
        filters=[("ensg_id", "in", genes)],
        columns=["model_id", "ensg_id", "confidence", "stratum_rank",
                 "class", "rank_basis"])
    pred["model_id"] = pred["model_id"].str.lower()
    pred["sensitive"] = [m in truth.get(g, ()) for m, g in
                         zip(pred["model_id"], pred["ensg_id"])]

    tiers = pred.groupby("confidence").agg(
        n_pairs=("sensitive", "size"),
        n_sensitive=("sensitive", "sum"),
        n_genes=("ensg_id", "nunique"))
    tiers["rate"] = tiers["n_sensitive"] / tiers["n_pairs"]

    # Wilson interval -- the pairs are not independent, so treat as indicative.
    z = 1.96
    p, nn = tiers["rate"], tiers["n_pairs"]
    denom = 1 + z ** 2 / nn
    centre = (p + z ** 2 / (2 * nn)) / denom
    half = z * np.sqrt(p * (1 - p) / nn + z ** 2 / (4 * nn ** 2)) / denom
    tiers["lo"], tiers["hi"] = centre - half, centre + half

    # Top-20 precision per tier: does a high-confidence gene rank better?
    top20 = []
    for g, grp in pred.groupby("ensg_id"):
        known = truth.get(g, set())
        ranked = grp.sort_values("stratum_rank", ascending=False)["model_id"].tolist()
        hr = hit_at_k(ranked, known, 20)
        if np.isfinite(hr):
            top20.append({"ensg_id": g,
                          "confidence": grp["confidence"].mode().iloc[0],
                          "hit_at_20": float(hr), "known": int(len(known & set(ranked)))})

    t20 = pd.DataFrame(top20)
    by_tier = (t20.groupby("confidence")["hit_at_20"]
               .agg(["count", "mean", "median"]).to_dict("index"))

    S["confidence"] = {
        "pair_level": tiers.reset_index().to_dict("records"),
        "gene_level_hit20": {k: {kk: float(vv) for kk, vv in v.items()}
                             for k, v in by_tier.items()},
        "gene_level_points": top20,
        "overall_rate": float(pred["sensitive"].mean()),
        "n_pairs": int(len(pred)),
        "caveat": ("pairs within a gene are not independent; intervals are "
                   "indicative, and tiers differ in gene composition"),
    }
    print(tiers.to_string())
    print("gene-level hit@20 by tier:", by_tier)


# ─────────────────────────────────────────────────────────────────────────
# 12. Ablation: what each component adds to held-out hit@20
# ─────────────────────────────────────────────────────────────────────────
def ablation():
    section("12. Ablation ladder")
    with open(PIPE / "stage4_eval_summary.json", encoding="utf-8") as fh:
        summ = json.load(fh)
    with open(PIPE / "test_run_tsg_two_hit_results.json", encoding="utf-8") as fh:
        tsg = json.load(fh)
    with open(PIPE / "test_run_signature_discount_results.json", encoding="utf-8") as fh:
        sig = json.load(fh)

    ev = pd.read_parquet(PIPE / "stage4_eval_consistent_denom.parquet")

    # Paired bootstrap on the per-gene hit rates for the two live steps.
    def paired_ci(a, b, B=5000):
        n = len(a)
        idx = RNG.integers(0, n, size=(B, n))
        diffs = (a[idx] - b[idx]).mean(axis=1)
        return float((a - b).mean()), float(np.percentile(diffs, 2.5)), \
            float(np.percentile(diffs, 97.5))

    flat = ev["hr_flat"].values
    driver = ev["hr_driver"].values
    anyalt = ev["hr_any"].values
    d_driver = paired_ci(driver, flat)
    d_any = paired_ci(anyalt, flat)

    steps = [
        {"label": "Flat core_score\n(abundance only)",
         "value": float(flat.mean()), "delta": None, "lo": None, "hi": None,
         "kind": "base"},
        {"label": "+ driver gating\n(Stage 4, production)",
         "value": float(driver.mean()), "delta": d_driver[0],
         "lo": d_driver[1], "hi": d_driver[2], "kind": "live"},
        {"label": "+ CNA modifier\n(Stage 5)",
         "value": float(driver.mean()), "delta": 0.0, "lo": 0.0, "hi": 0.0,
         "kind": "inert",
         "note": "CNA never enters the sort key -- confidence label only"},
        {"label": "+ Chronos check\n(Stage 5)",
         "value": float(driver.mean()), "delta": 0.0, "lo": 0.0, "hi": 0.0,
         "kind": "inert",
         "note": "read-only per-gene validation; rejected as a ranking layer"},
        {"label": "+ TSG role prior\n(score inversion)",
         "value": float(tsg["hit_at_20_mean"]["flat_low_abundance_TSG"]),
         "delta": float(tsg["hit_at_20_mean"]["flat_low_abundance_TSG"]
                        - tsg["hit_at_20_mean"]["flat_high_abundance"]),
         "lo": None, "hi": None, "kind": "harmful",
         "note": "measured on 14 TSG genes, not the held-out split"},
    ]
    S["ablation"] = {
        "steps": steps,
        "any_alteration_reference": {"value": float(anyalt.mean()),
                                     "delta": d_any[0], "lo": d_any[1], "hi": d_any[2]},
        "variance_reduction_pct": summ["variance_reduction"]["reduction_pct"],
        "delta_std_any": summ["variance_reduction"]["delta_std_any_alt"],
        "delta_std_driver": summ["variance_reduction"]["delta_std_driver"],
        "n_genes": int(len(ev)),
        "tsg_detail": tsg["hit_at_20_mean"],
        "signature_discount": sig["paired_bootstrap_any_adj_minus_any"],
    }
    for s in steps:
        print(f"  {s['label'][:34]:<36} {s['value']:.5f}  "
              f"delta={s['delta'] if s['delta'] is None else round(s['delta'], 5)}")


# ─────────────────────────────────────────────────────────────────────────
# 13. Negative results panel
# ─────────────────────────────────────────────────────────────────────────
def negatives():
    section("13. Negative results")
    def rd(name):
        with open(PIPE / name, encoding="utf-8") as fh:
            return json.load(fh)

    sig = rd("test_run_signature_discount_results.json")
    met_mir = rd("test_run_metabolomics_mirna_results.json")
    met_conf = rd("test_run_metabolomics_lineage_confound_results.json")
    mir_conf = rd("test_run_mirna_lineage_confound_results.json")
    ploidy = rd("test_run_ploidy_normalised_cna_results.json")
    proj = rd("test_run_project_score_replication_results.json")

    met_rho = [v["rho"] for v in met_mir["metabolomics_test"].values()]

    S["negatives"] = {
        "signature_discount": {
            "title": "Mutation-signature discount",
            "hypothesis": "Discounting hypermutated lines shrinks routing variance",
            "point": sig["paired_bootstrap_any_adj_minus_any"]["point_estimate"],
            "lo": sig["paired_bootstrap_any_adj_minus_any"]["ci_95_low"],
            "hi": sig["paired_bootstrap_any_adj_minus_any"]["ci_95_high"],
            "before": sig["hit_at_20"]["any_alteration_baseline"],
            "after": sig["hit_at_20"]["any_alteration_DISCOUNT"],
            "std_before": sig["activation_driven_delta_std"]["any_alteration"],
            "std_after": sig["activation_driven_delta_std"]["any_alteration_discount"],
            "verdict": "inert -- CI includes 0, variance did not shrink",
        },
        "metabolomics": {
            "title": "Metabolomics as a global per-line scalar",
            "anova_r2": met_conf["anova_r2"],
            "anova_p": met_conf["anova_p"],
            "n_lineages": met_conf["n_lineage_groups"],
            "per_gene_rho": met_rho,
            "verdict": "rejected -- gene-agnostic scalar, weak lineage proxy",
        },
        "mirna": {
            "title": "miRNA-as-inhibitor dampening",
            "anova_r2": mir_conf["anova_r2"],
            "anova_p": mir_conf["anova_p"],
            "n_lineages": mir_conf["n_lineage_groups"],
            "cases": met_mir["mirna_inhibitor_test"],
            "verdict": "rejected -- rho moved the wrong way in both test genes",
        },
        "ploidy": {
            "title": "Ploidy-normalised CNA thresholds",
            "rho_before": ploidy["PRIMARY_ploidy_correlation"]["absolute (current)"]["rho"],
            "p_before": ploidy["PRIMARY_ploidy_correlation"]["absolute (current)"]["p"],
            "rho_after": ploidy["PRIMARY_ploidy_correlation"]["ploidy-relative"]["rho"],
            "p_after": ploidy["PRIMARY_ploidy_correlation"]["ploidy-relative"]["p"],
            "wgd_before": ploidy["PRIMARY_wgd_bias"]["absolute (current)"]["ratio"],
            "wgd_after": ploidy["PRIMARY_wgd_bias"]["ploidy-relative"]["ratio"],
            "blast_radius": ploidy["SECONDARY_confidence_blast_radius"],
            "verdict": "artifact removed, WGD excess survives -- precision, not performance",
        },
        "project_score": {
            "title": "Abundance vs dependency, replicated on Project Score",
            "mean_rho": proj["project_score"]["mean_rho"],
            "chronos_mean_rho": proj["chronos"]["mean_rho"],
            "rho_of_rhos": proj["cross_dataset_agreement"]["rho_of_rhos"],
            "n_genes": proj["coverage"]["genes_with_min_n"],
            "verdict": "replicates -- abundance does not track dependency",
        },
    }
    print(json.dumps({k: v.get("verdict") for k, v in S["negatives"].items()}, indent=2))


# ─────────────────────────────────────────────────────────────────────────
def main():
    steps = [counts, modality_presence, gene_axis_reach, interlayer_correlation,
             normalisation, score_vs_layers, tsg_deletion, routing, recipes,
             per_gene_hits, confidence_vs_outcome, ablation, negatives]
    only = sys.argv[1:] or None
    for fn in steps:
        if only and fn.__name__ not in only:
            continue
        try:
            fn()
        except Exception as exc:                       # noqa: BLE001
            import traceback
            print(f"\n!! {fn.__name__} FAILED: {exc}")
            traceback.print_exc()
            S.setdefault("_errors", {})[fn.__name__] = str(exc)

    if STATS.exists() and only:
        with open(STATS, encoding="utf-8") as fh:
            prev = json.load(fh)
        prev.update(S)
        payload = prev
    else:
        payload = S
    with open(STATS, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    print(f"\nWrote {STATS}  ({STATS.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
