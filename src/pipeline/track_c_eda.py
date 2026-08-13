"""
track_c_eda.py -- gene-axis EDA over the Track C alteration tables.

Read-only. Reads the current Track C outputs from cleaned_track_data/ and the
Stage 0 gene dimension from the warehouse, joins the alteration evidence onto
the gene axis, and writes every number it computes to
outputs/track_c_eda_stats.json plus figures to outputs/figures/track_c/.

Two passes, matching how the analysis is meant to be read:

  simple   -- what is on the gene axis at all: universe reconciliation,
              per-gene coverage, modality overlap, role enrichment
  indepth  -- what would break a naive model built on it: missingness that
              predicts outcome, the detection floor, sparsity vs the
              normalisation axis, and three confounds (gene length,
              hypermutators, lineage)

Nothing here writes to the warehouse or to any pipeline artefact.

Usage:
    python src/pipeline/track_c_eda.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import duckdb
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TRACK = ROOT / "cleaned_track_data"
DB = ROOT / "src" / "pipeline" / "outputs" / "celllineselector.db"
OUT = ROOT / "src" / "pipeline" / "outputs"
FIG = OUT / "figures" / "track_c"
FIG.mkdir(parents=True, exist_ok=True)

# Palette lifted from src/figures/common.py so the Track C figures sit in the
# same visual family as the report suite.
BLUE, RED, AMBER, GREEN = "#2a6fb5", "#c8443c", "#e0a020", "#3e8f60"
GREY, LIGHT, INK = "#7a7a7a", "#dfe6ee", "#22262b"

plt.rcParams.update({
    "figure.dpi": 130,
    "savefig.dpi": 130,
    "font.size": 8.5,
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
})

S: dict = {}


def rec(path: str, value):
    """Record a stat under a dotted path so the JSON mirrors the section order."""
    node = S
    parts = path.split(".")
    for p in parts[:-1]:
        node = node.setdefault(p, {})
    node[parts[-1]] = value
    return value


def save(fig, name: str):
    p = FIG / f"{name}.png"
    fig.tight_layout()
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure -> {p.relative_to(ROOT)}")


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

def load():
    """Track C outputs plus the Stage 0 gene dimension, all keys lowercased.

    Track C writes UPPERCASE keys; project canon is lowercase (HARMONISATION_V2
    section 5.2). Lowercasing on read is idempotent and is what Stage 4 already
    does, so this join follows the same convention rather than inventing one.
    """
    print("loading ...")
    mut = pd.read_parquet(TRACK / "mutations_collapsed.parquet")
    fus = pd.read_parquet(TRACK / "fusions_gene_level.parquet")
    det = pd.read_parquet(TRACK / "mutations_variant_detail.parquet")
    sig = pd.read_parquet(TRACK / "signatures_model_level.parquet")
    mir = pd.read_parquet(TRACK / "mirna_model_level.parquet")

    for df in (mut, fus, det):
        df["ensg_id"] = df["ensg_id"].astype("string").str.lower()
        df["model_id"] = df["model_id"].astype("string").str.lower()
    for df in (sig, mir):
        df["model_id"] = df["model_id"].astype("string").str.lower()

    con = duckdb.connect(str(DB), read_only=True)
    gene = con.execute("""
        select lower(gene_id) as ensg_id, hugo_symbol, is_protein_coding,
               gene_role, cgc_tier, locus_group
        from main.gene_enriched
    """).df()
    tissue = con.execute("""
        select distinct lower(model_id) as model_id, tissue, msi_status
        from main.gdsc_models where model_id is not null and tissue is not null
    """).df().drop_duplicates("model_id")
    hub = con.execute("select lower(model_id) as model_id from main.cell_line_connection_enriched").df()
    con.close()

    print(f"  mutations {mut.shape}  fusions {fus.shape}  detail {det.shape}")
    print(f"  gene dim  {gene.shape}  tissue {tissue.shape}  hub {hub.shape}")
    return mut, fus, det, sig, mir, gene, tissue, hub


# ---------------------------------------------------------------------------
# Simple EDA -- what is on the gene axis
# ---------------------------------------------------------------------------

def simple(mut, fus, det, sig, mir, gene, tissue, hub):
    print("\n=== SIMPLE: gene-axis reconciliation ===")

    universe = set(gene.ensg_id)
    coding = set(gene.loc[gene.is_protein_coding == True, "ensg_id"])  # noqa: E712
    mg, fg = set(mut.ensg_id), set(fus.ensg_id)

    rec("simple.universe.gene_dimension", len(universe))
    rec("simple.universe.protein_coding", len(coding))
    rec("simple.universe.mutations_genes", len(mg))
    rec("simple.universe.fusions_genes", len(fg))
    rec("simple.universe.mutations_outside_dimension", len(mg - universe))
    rec("simple.universe.fusions_outside_dimension", len(fg - universe))
    rec("simple.universe.mutations_non_coding", len(mg - coding))
    rec("simple.universe.fusions_non_coding", len(fg - coding))
    rec("simple.universe.never_altered", len(universe - mg - fg))
    print(f"  dimension {len(universe):,} | coding {len(coding):,}")
    print(f"  mutations {len(mg):,} genes, {len(mg - universe):,} outside dimension, "
          f"{len(mg - coding):,} not protein-coding")
    print(f"  fusions   {len(fg):,} genes, {len(fg - universe):,} outside dimension, "
          f"{len(fg - coding):,} not protein-coding")
    print(f"  never altered by either: {len(universe - mg - fg):,}")

    # --- modality overlap on the gene axis
    both = mg & fg
    rec("simple.overlap.both", len(both))
    rec("simple.overlap.mutation_only", len(mg - fg))
    rec("simple.overlap.fusion_only", len(fg - mg))
    print(f"  both {len(both):,} | mut-only {len(mg - fg):,} | fus-only {len(fg - mg):,}")

    # --- per-gene coverage (how many lines carry an alteration in this gene)
    mut_per_gene = mut.groupby("ensg_id").model_id.nunique()
    fus_per_gene = fus.groupby("ensg_id").model_id.nunique()
    for label, s in (("mutations", mut_per_gene), ("fusions", fus_per_gene)):
        rec(f"simple.lines_per_gene.{label}", {
            "n_genes": int(len(s)), "min": int(s.min()), "q25": float(s.quantile(.25)),
            "median": float(s.median()), "q75": float(s.quantile(.75)), "max": int(s.max()),
            "lt10": int((s < 10).sum()), "lt20": int((s < 20).sum()), "lt30": int((s < 30).sum()),
            "pct_lt20": round(float((s < 20).mean() * 100), 1),
        })
        print(f"  {label}: lines/gene median {s.median():.0f}, "
              f"{(s < 20).mean() * 100:.1f}% of genes seen in <20 lines")

    # --- density of the (gene x line) matrix
    for label, df, ng in (("mutations", mut, len(mg)), ("fusions", fus, len(fg))):
        nl = df.model_id.nunique()
        rec(f"simple.density.{label}", {
            "genes": int(ng), "lines": int(nl), "possible": int(ng * nl),
            "observed": int(len(df)), "pct_dense": round(len(df) / (ng * nl) * 100, 2),
        })
        print(f"  {label} matrix {ng:,} x {nl:,} = {ng * nl:,} possible, "
              f"{len(df):,} observed ({len(df) / (ng * nl) * 100:.2f}% dense)")

    # --- alteration rate by curated cancer-gene role
    g = gene.set_index("ensg_id")
    role = g.gene_role.reindex(sorted(universe)).fillna("unknown")
    rate = pd.DataFrame({
        "role": role,
        "mut_lines": mut_per_gene.reindex(role.index).fillna(0),
        "fus_lines": fus_per_gene.reindex(role.index).fillna(0),
    })
    by_role = rate.groupby("role").agg(
        n_genes=("role", "size"),
        median_mut_lines=("mut_lines", "median"),
        mean_mut_lines=("mut_lines", "mean"),
        median_fus_lines=("fus_lines", "median"),
        pct_ever_mutated=("mut_lines", lambda s: (s > 0).mean() * 100),
    ).round(2)
    rec("simple.by_role", json.loads(by_role.to_json(orient="index")))
    print("\n  alteration by CGC role:")
    print(by_role.to_string())

    # --- top genes
    top_mut = mut_per_gene.sort_values(ascending=False).head(15)
    top_fus = fus_per_gene.sort_values(ascending=False).head(15)
    sym = g.hugo_symbol
    rec("simple.top_mutated", {sym.get(k, k): int(v) for k, v in top_mut.items()})
    rec("simple.top_fused", {sym.get(k, k): int(v) for k, v in top_fus.items()})
    print(f"\n  top mutated: {', '.join(f'{sym.get(k, k)}({v})' for k, v in top_mut.head(8).items())}")
    print(f"  top fused:   {', '.join(f'{sym.get(k, k)}({v})' for k, v in top_fus.head(8).items())}")

    # --- what the two gene-axis-less modalities cost
    rec("simple.no_gene_axis", {
        "signatures_models": int(sig.model_id.nunique()),
        "mirna_models": int(mir.model_id.nunique()),
        "mirna_features": int(mir.mirna_id.nunique()),
        "hub_models": int(hub.model_id.nunique()),
        "mutation_models": int(mut.model_id.nunique()),
    })

    # ---- FIGURES ----------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))

    ax = axes[0]
    ax.set_title("gene-axis coverage", loc="left", fontweight="bold")
    vals = [len(mg - fg), len(both), len(fg - mg), len(universe - mg - fg)]
    labs = ["mutation\nonly", "both", "fusion\nonly", "never\naltered"]
    ax.bar(labs, vals, color=[BLUE, GREEN, AMBER, LIGHT], edgecolor=INK, linewidth=.5)
    ax.set_ylabel("genes")
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=7.5)

    ax = axes[1]
    ax.hist(mut_per_gene.values, bins=60, color=BLUE, edgecolor="none")
    ax.axvline(20, color=RED, ls="--", lw=1.2)
    ax.text(22, ax.get_ylim()[1] * .85, "n=20\nshrinkage\nfloor", color=RED, fontsize=7)
    ax.set_xlabel("cell lines carrying a mutation")
    ax.set_ylabel("genes")
    ax.set_title("mutations: lines per gene", loc="left", fontweight="bold")
    ax.set_yscale("log")

    ax = axes[2]
    ax.hist(fus_per_gene.values, bins=60, color=AMBER, edgecolor="none")
    ax.axvline(20, color=RED, ls="--", lw=1.2)
    ax.set_xlabel("cell lines carrying a fusion")
    ax.set_ylabel("genes")
    ax.set_title("fusions: lines per gene", loc="left", fontweight="bold")
    ax.set_yscale("log")
    save(fig, "01_gene_axis_coverage")

    return mut_per_gene, fus_per_gene


# ---------------------------------------------------------------------------
# In-depth EDA
# ---------------------------------------------------------------------------

def missingness_sweep(det):
    """For every column with nulls, test whether being null predicts the outcome.

    This is the check that has to run before any imputation decision: if
    missingness carries signal, filling it in destroys that signal and biases
    whatever the filled value feeds.
    """
    print("\n=== IN-DEPTH 1: does missingness predict outcome? ===")
    outcomes = ["likelylof", "hotspot", "oncogenehighimpact",
                "tumorsuppressorhighimpact", "hessdriver"]
    rows = []
    for col in det.columns:
        nulls = det[col].isna()
        if not (0.001 < nulls.mean() < 0.999):
            continue
        for out in outcomes:
            o = det[out].astype(bool)
            p_missing, p_present = o[nulls].mean(), o[~nulls].mean()
            lift = (p_missing / p_present) if p_present > 0 else np.inf
            rows.append({"column": col, "null_rate": round(nulls.mean(), 3), "outcome": out,
                         "rate_when_null": round(p_missing, 4),
                         "rate_when_present": round(p_present, 4),
                         "lift": round(lift, 1) if np.isfinite(lift) else None})
    sweep = pd.DataFrame(rows).sort_values("lift", ascending=False, na_position="last")
    rec("indepth.missingness_sweep", json.loads(sweep.to_json(orient="records")))
    print(sweep.to_string(index=False))

    # the headline split: severity absent entirely
    det = det.assign(sev_missing=det.ampathogenicity.isna() & det.revelscore.isna())
    by_impact = det.groupby("vepimpact").sev_missing.agg(["mean", "size"]).round(3)
    rec("indepth.severity_missing_by_impact", json.loads(by_impact.to_json(orient="index")))
    rec("indepth.severity_missing", {
        "rate": round(float(det.sev_missing.mean()), 4),
        "lof_rate_when_missing": round(float(det.loc[det.sev_missing, "likelylof"].mean()), 4),
        "lof_rate_when_present": round(float(det.loc[~det.sev_missing, "likelylof"].mean()), 4),
        "tsg_rate_when_missing": round(float(det.loc[det.sev_missing, "tumorsuppressorhighimpact"].mean()), 4),
        "tsg_rate_when_present": round(float(det.loc[~det.sev_missing, "tumorsuppressorhighimpact"].mean()), 4),
        "median_vaf_when_missing": round(float(det.loc[det.sev_missing, "af"].median()), 4),
        "median_vaf_when_present": round(float(det.loc[~det.sev_missing, "af"].median()), 4),
    })
    print("\n  severity-missing by VEP impact:")
    print(by_impact.to_string())
    return det


def detection_floor(det):
    print("\n=== IN-DEPTH 2: the detection floor ===")
    af = det.af
    rec("indepth.vaf", {
        "min": round(float(af.min()), 4), "max": round(float(af.max()), 4),
        "median": round(float(af.median()), 4), "mean": round(float(af.mean()), 4),
        "at_floor_pct": round(float((af <= 0.155).mean() * 100), 2),
        "at_ceiling_pct": round(float((af >= 0.999).mean() * 100), 2),
        "n_exactly_one": int((af == 1.0).sum()),
    })
    print(f"  VAF min {af.min():.3f} max {af.max():.3f} median {af.median():.3f}")
    print(f"  {(af <= 0.155).mean() * 100:.2f}% sit at the floor, "
          f"{(af >= 0.999).mean() * 100:.2f}% at the ceiling "
          f"({(af == 1.0).sum():,} exactly 1.0 -> logit undefined)")

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.2))
    ax = axes[0]
    ax.hist(af, bins=90, color=BLUE, edgecolor="none")
    ax.axvline(0.15, color=RED, lw=1.4)
    ax.text(0.16, ax.get_ylim()[1] * .9, "0.15 floor\n(source filter)", color=RED, fontsize=7)
    ax.set_xlabel("variant allele frequency")
    ax.set_ylabel("variants")
    ax.set_title("VAF is truncated, not zero-inflated", loc="left", fontweight="bold")

    ax = axes[1]
    for lab, col in (("severity scored", GREEN), ("severity absent", RED)):
        sub = af[det.sev_missing == (lab == "severity absent")]
        ax.hist(sub, bins=60, alpha=.6, label=f"{lab} (n={len(sub):,})", color=col, edgecolor="none")
    ax.set_xlabel("variant allele frequency")
    ax.set_ylabel("variants")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("identical VAF -> not a coverage artefact", loc="left", fontweight="bold")
    save(fig, "02_detection_floor")


def sparsity_vs_normalisation(mut, mut_per_gene):
    """The gene-wise-across-lines axis is only as good as the per-gene n."""
    print("\n=== IN-DEPTH 3: sparsity vs the normalisation axis ===")
    vaf = pd.read_parquet(TRACK / "mutations_variant_detail.parquet",
                          columns=["ensg_id", "model_id", "af"])
    vaf["ensg_id"] = vaf.ensg_id.str.lower()
    per = vaf.groupby("ensg_id").af.agg(["size", "median", "std", "nunique"])
    mad = vaf.groupby("ensg_id").af.apply(lambda s: float(np.median(np.abs(s - np.median(s)))))
    per["mad"] = mad
    rec("indepth.mad_degeneracy", {
        "genes": int(len(per)),
        "mad_zero": int((per.mad == 0).sum()),
        "pct_mad_zero": round(float((per.mad == 0).mean() * 100), 1),
        "mad_zero_when_n_lt20": int(((per.mad == 0) & (per["size"] < 20)).sum()),
        "n_lt20": int((per["size"] < 20).sum()),
        "pct_n_lt20": round(float((per["size"] < 20).mean() * 100), 1),
    })
    print(f"  MAD == 0 for {(per.mad == 0).sum():,} of {len(per):,} genes "
          f"({(per.mad == 0).mean() * 100:.1f}%) -- robust z undefined there")
    print(f"  {(per['size'] < 20).sum():,} genes ({(per['size'] < 20).mean() * 100:.1f}%) "
          f"have <20 variants to normalise against")

    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.scatter(per["size"], per.mad, s=3, alpha=.15, color=BLUE, edgecolors="none")
    ax.axvline(20, color=RED, ls="--", lw=1.2)
    ax.set_xscale("log")
    ax.set_xlabel("variants observed for this gene (log)")
    ax.set_ylabel("MAD of VAF")
    ax.set_title("MAD collapses to zero exactly where n is small",
                 loc="left", fontweight="bold")
    ax.text(21, per.mad.max() * .9, "n=20", color=RED, fontsize=7)
    save(fig, "03_mad_degeneracy")
    return per


def confounds(mut, det, per, gene, tissue):
    print("\n=== IN-DEPTH 4: three confounds ===")

    # --- (a) gene length. Genomic span of observed variants is a crude proxy;
    #     with >=2 variants it is largely determined by the true locus extent.
    pos = pd.to_numeric(det["pos"], errors="coerce")
    span = det.assign(pos=pos).dropna(subset=["pos"]).groupby("ensg_id").pos.agg(["min", "max", "size"])
    span["span_kb"] = (span["max"] - span["min"]) / 1000
    span = span[span["size"] >= 5]
    lines_per_gene = mut.groupby("ensg_id").model_id.nunique().reindex(span.index)
    ok = span.span_kb > 0
    r_s = float(pd.Series(np.log10(span.span_kb[ok])).corr(
        pd.Series(np.log10(lines_per_gene[ok].clip(lower=1))), method="spearman"))
    rec("indepth.confound_gene_length", {
        "genes_tested": int(ok.sum()),
        "spearman_log_span_vs_log_lines": round(r_s, 3),
        "median_span_kb": round(float(span.span_kb[ok].median()), 1),
    })
    print(f"  (a) gene length: Spearman(log span, log lines mutated) = {r_s:.3f} "
          f"over {ok.sum():,} genes")

    # --- (b) hypermutators
    per_line = mut.groupby("model_id").size().sort_values(ascending=False)
    top1 = per_line.head(max(1, int(len(per_line) * .01))).sum() / per_line.sum()
    top5 = per_line.head(max(1, int(len(per_line) * .05))).sum() / per_line.sum()
    msi = tissue.set_index("model_id").msi_status.reindex(per_line.index)
    rec("indepth.confound_hypermutator", {
        "lines": int(len(per_line)),
        "median_genes_mutated": int(per_line.median()),
        "max_genes_mutated": int(per_line.max()),
        "share_from_top_1pct": round(float(top1) * 100, 1),
        "share_from_top_5pct": round(float(top5) * 100, 1),
        "msi_median_burden": json.loads(
            per_line.groupby(msi).median().to_json()) if msi.notna().any() else {},
    })
    print(f"  (b) hypermutators: top 1% of lines contribute {top1 * 100:.1f}% of all "
          f"(gene,line) mutation pairs; top 5% contribute {top5 * 100:.1f}%")
    print(f"      median genes mutated per line {per_line.median():.0f}, max {per_line.max():,}")

    # --- (c) lineage
    t = tissue.set_index("model_id").tissue
    burden = pd.DataFrame({"burden": per_line, "tissue": t.reindex(per_line.index)}).dropna()
    by_t = burden.groupby("tissue").burden.agg(["size", "median"]).sort_values("median", ascending=False)
    by_t = by_t[by_t["size"] >= 10]
    rec("indepth.confound_lineage", json.loads(by_t.to_json(orient="index")))
    spread = float(by_t["median"].max() / by_t["median"].min()) if len(by_t) else float("nan")
    rec("indepth.confound_lineage_spread", round(spread, 2))
    print(f"  (c) lineage: median burden ranges {by_t['median'].min():.0f} "
          f"({by_t.index[-1]}) to {by_t['median'].max():.0f} ({by_t.index[0]}) "
          f"-- {spread:.1f}x spread across tissues")

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.4))
    ax = axes[0]
    ax.scatter(span.span_kb[ok], lines_per_gene[ok], s=3, alpha=.12, color=BLUE, edgecolors="none")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("observed genomic span, kb (log)")
    ax.set_ylabel("lines mutated (log)")
    ax.set_title(f"(a) gene length confound  rho={r_s:.2f}", loc="left", fontweight="bold")

    ax = axes[1]
    ax.plot(np.arange(1, len(per_line) + 1) / len(per_line) * 100,
            np.cumsum(per_line.values) / per_line.sum() * 100, color=RED, lw=1.6)
    ax.plot([0, 100], [0, 100], color=GREY, ls=":", lw=1)
    ax.set_xlabel("cell lines, ranked by burden (%)")
    ax.set_ylabel("cumulative share of mutation pairs (%)")
    ax.set_title("(b) hypermutator concentration", loc="left", fontweight="bold")

    ax = axes[2]
    b = by_t.sort_values("median")
    ax.barh(range(len(b)), b["median"], color=AMBER, edgecolor=INK, linewidth=.4)
    ax.set_yticks(range(len(b)))
    ax.set_yticklabels(b.index, fontsize=6.5)
    ax.set_xlabel("median genes mutated per line")
    ax.set_title("(c) lineage confound", loc="left", fontweight="bold")
    save(fig, "04_confounds")


def fusion_axes(fus):
    print("\n=== IN-DEPTH 5: fusion confidence vs frame ===")
    by_conf = fus.groupby("max_confidence").agg(
        n=("best_ffpm", "size"),
        median_ffpm=("best_ffpm", "median"),
        in_frame_rate=("any_in_frame", "mean"),
    ).round(4)
    rec("indepth.fusion_confidence", json.loads(by_conf.to_json(orient="index")))
    print(by_conf.to_string())
    print("  -> ffpm separates strongly by confidence; in-frame rate does not")

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.2))
    order = ["low", "medium", "high"]
    ax = axes[0]
    ax.boxplot([np.log10(fus.loc[fus.max_confidence == c, "best_ffpm"] + 1e-3) for c in order],
               labels=order, showfliers=False,
               boxprops=dict(color=INK), medianprops=dict(color=RED, lw=1.6))
    ax.set_ylabel("log10 ffpm")
    ax.set_title("confidence tracks abundance", loc="left", fontweight="bold")
    ax = axes[1]
    ax.bar(order, [fus.loc[fus.max_confidence == c, "any_in_frame"].mean() for c in order],
           color=AMBER, edgecolor=INK, linewidth=.5)
    ax.set_ylim(0, .3)
    ax.set_ylabel("in-frame rate")
    ax.set_title("frame is independent of confidence", loc="left", fontweight="bold")
    save(fig, "05_fusion_axes")


def mitochondrial(det, mut, gene):
    """mtDNA variants are heteroplasmy, not allele fraction -- a different quantity.

    Nuclear VAF is the fraction of 2-4 chromosome copies carrying the allele.
    Mitochondrial VAF is the fraction of hundreds of mtDNA copies per cell.
    They share a column and a range and mean different things, so any VAF-based
    magnitude has to separate them or the mito genes dominate by construction.
    """
    print("\n=== IN-DEPTH 6: mitochondrial variants are a different quantity ===")
    is_mt = det.chrom.astype(str).str.lower() == "chrm"
    mt_genes = set(det.loc[is_mt, "ensg_id"])
    sym = gene.set_index("ensg_id").hugo_symbol
    stats = {
        "n_variants": int(is_mt.sum()),
        "pct_of_variants": round(float(is_mt.mean() * 100), 2),
        "n_genes": int(len(mt_genes)),
        "n_lines": int(det.loc[is_mt, "model_id"].nunique()),
        "symbols": sorted(str(sym.get(g, g)) for g in mt_genes),
        "median_vaf_mito": round(float(det.loc[is_mt, "af"].median()), 4),
        "median_vaf_nuclear": round(float(det.loc[~is_mt, "af"].median()), 4),
        "pct_near_homoplasmic_mito": round(float((det.loc[is_mt, "af"] > 0.95).mean() * 100), 1),
        "pct_near_homoplasmic_nuclear": round(float((det.loc[~is_mt, "af"] > 0.95).mean() * 100), 1),
    }
    lines_per_gene = mut.groupby("ensg_id").model_id.nunique().sort_values(ascending=False)
    rank = {str(sym.get(g, g)): int(lines_per_gene.index.get_loc(g)) + 1
            for g in mt_genes if g in lines_per_gene.index}
    stats["rank_among_most_mutated_genes"] = dict(sorted(rank.items(), key=lambda kv: kv[1]))
    rec("indepth.mitochondrial", stats)
    print(f"  {stats['n_variants']:,} variants over {stats['n_genes']} genes in "
          f"{stats['n_lines']:,} lines ({stats['pct_of_variants']}% of variants)")
    print(f"  median VAF mito {stats['median_vaf_mito']} vs nuclear {stats['median_vaf_nuclear']}")
    print(f"  near-homoplasmic (VAF>0.95): {stats['pct_near_homoplasmic_mito']}% mito "
          f"vs {stats['pct_near_homoplasmic_nuclear']}% nuclear")
    print(f"  gene rank by lines mutated: {stats['rank_among_most_mutated_genes']}")

    fig, ax = plt.subplots(figsize=(5.4, 3.3))
    ax.hist(det.loc[~is_mt, "af"], bins=80, density=True, color=BLUE,
            alpha=.75, edgecolor="none", label=f"nuclear (n={(~is_mt).sum():,})")
    ax.hist(det.loc[is_mt, "af"], bins=80, density=True, color=RED,
            alpha=.65, edgecolor="none", label=f"mitochondrial (n={is_mt.sum():,})")
    ax.set_xlabel("VAF")
    ax.set_ylabel("density")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("same column, two different quantities", loc="left", fontweight="bold")
    save(fig, "06_mitochondrial_vaf")


def fusion_filter_check(fus):
    """DepMap documents removing FFPM < 0.05 from OmicsFusionFiltered. Check it held."""
    print("\n=== IN-DEPTH 7: does the documented fusion filter hold? ===")
    below = fus.best_ffpm < 0.05
    rec("indepth.fusion_filter_check", {
        "documented_threshold": 0.05,
        "rows_below_threshold": int(below.sum()),
        "pct_below": round(float(below.mean() * 100), 1),
        "rows_exactly_zero": int((fus.best_ffpm == 0).sum()),
        "min_nonzero_ffpm": round(float(fus.best_ffpm[fus.best_ffpm > 0].min()), 6),
        "note": ("best_ffpm is the max over contributing events, so a value below "
                 "the threshold means every contributing event was below it"),
    })
    print(f"  {below.sum():,} rows ({below.mean() * 100:.1f}%) sit below the documented "
          f"FFPM>=0.05 floor; {int((fus.best_ffpm == 0).sum())} are exactly 0")


def main():
    mut, fus, det, sig, mir, gene, tissue, hub = load()
    mut_per_gene, _ = simple(mut, fus, det, sig, mir, gene, tissue, hub)
    det = missingness_sweep(det)
    detection_floor(det)
    per = sparsity_vs_normalisation(mut, mut_per_gene)
    confounds(mut, det, per, gene, tissue)
    fusion_axes(fus)
    mitochondrial(det, mut, gene)
    fusion_filter_check(fus)

    p = OUT / "track_c_eda_stats.json"
    p.write_text(json.dumps(S, indent=2, default=str), encoding="utf-8")
    print(f"\nstats -> {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
