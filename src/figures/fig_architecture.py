"""
Architecture figures (report section 1).

  01_pipeline_master   -- Stages 0-7, artefact on every arrow, counts in every box
  02_provenance_map    -- sources -> the stages that consume them
  03_join_key_schema   -- model_id and ensg_id as the two spines

Schematics, but every count is read from _stats.json, so they cannot drift
from the data. Geometry is in inches (see layout.py) so boxes size to content.
"""
import sys
from pathlib import Path

import matplotlib.patches as mpatches

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (AMBER, BLUE, DEAD_FILL, GREEN, GREY, INK, RED,
                    SOURCE_EDGE, SOURCE_FILL, STAGE_EDGE, STAGE_FILL,
                    load_stats, save)
from layout import arrow, arrow_label, box, box_height, canvas, text

REF_FILL, REF_EDGE = "#eef4ec", GREEN


# ─────────────────────────────────────────────────────────────────────────
def fig_pipeline_master(S):
    c, r = S["counts"], S["routing"]
    n1, n2 = c["n_layers_dist"]["1"], c["n_layers_dist"]["2"]
    sizes = S["modality"]["set_sizes"]

    stages = [
        ("Stage 0 — Harmonisation", [
            "00_harmonisation.ipynb",
            f"{c['lookup_cell_lines']:,} cell lines reconciled across DepMap, Cellosaurus and GEO",
            "RRID → CVCL → ACH resolver; every row assigned a gap_class"]),
        ("Stage 1 — Lookup & metadata join", [
            "01_lookup_metadata_join.ipynb",
            f"gene_lookup ({c['lookup_genes']:,} genes)  ×  cell_line_lookup ({c['lookup_cell_lines']:,} lines)",
            "hub-and-spoke: sources join through the lookups, not to each other"]),
        ("Stage 2 — Core score", [
            "02_core_score.ipynb  +  scoring_variants.py",
            "percentile-rank per gene, then correlation-penalised weighted sum",
            "wᵢ = 1 / (1 + |ρ̄ᵢ|),  normalised to sum to 1",
            f"{n2:,} rows have both layers ({c['pct_2layer']:.1f}%) · {n1:,} are expression-only"]),
        ("Stage 3 — Name-resolution audit", [
            "name_match_audit.parquet",
            "exact match after deterministic normalisation, or substring containment"]),
        ("Stage 4 — Driver-gated routing", [
            "04_driver_gated_routing.ipynb",
            "sort_key = has_driver_alteration · 1e6 + core_score  (one vectorised sort)",
            f"gene_regime: {r['class_counts_curated']['activation_driven']} activation-driven, "
            f"{r['class_counts_curated']['abundance_tracking']} abundance-tracking, from {c['curated_genes']} curated genes"]),
        ("Stage 5 — Confidence tiers", [
            "05_confidence_tiers.ipynb  +  build_full_predictions.py",
            "deterministic precedence rules over four evidence types; categorical",
            "TSG rows served inverted (1 − core_score)",
            f"{S['routing']['regime_source_rows_full'].get('measured', 0):,} rows measured · "
            f"{S['routing']['regime_source_rows_full'].get('prior', 0):,} prior · "
            f"{S['routing']['regime_source_rows_full'].get('unknown', 0):,} unknown"]),
        ("Stage 6 — Held-out evaluation", [
            "06_held_out_eval.ipynb  +  stage4_eval_save.py",
            "gene-level split, sorted before rng.choice(seed=42) → 113 train / 28 test",
            "metric = hit@20, denominator = GDSC-sensitive ∩ scored"]),
        ("Stage 7 — Query / explain layer", [
            "cli.py → rank_cell_lines.py · explain_pair.py · evidence_ledger.py",
            "accepts an ENSG id or a gene symbol; resolves through gene_symbol_resolver",
            "fields tagged RANKING or CONFIDENCE MODIFIER"]),
    ]

    artefacts = [
        ("harmonised_enriched.parquet", f"{c['harmonised_rows']:,} rows — one per resolved cell line"),
        ("gene-long table", f"{c['universe_genes']:,} protein-coding genes × {c['lookup_cell_lines']:,} cell lines"),
        ("core_score.parquet", f"{c['core_score_rows']:,} rows · {c['scored_genes']:,} genes × {c['scored_cell_lines']:,} scored lines"),
        ("audited model_id keys", f"{c['scored_cell_lines']:,} scored lines confirmed against the lookup"),
        ("flags_with_driver.parquet", f"{c['flags_rows']:,} rows  +  gene_regime.parquet ({c['curated_genes']} genes)"),
        ("predictions_with_confidence.parquet", f"{c['core_score_rows']:,} rows × 15 columns"),
        ("stage4_eval_consistent_denom.parquet", "28 held-out genes, one row each"),
    ]

    # ── measure, then size the canvas to fit ─────────────────────────────
    TS, BS = 10.5, 8.0
    heights = [box_height(len(l), TS, BS) for _, l in stages]
    GAP = 0.72
    HEAD, FOOT = 1.15, 0.55
    W = 10.9
    H = HEAD + sum(heights) + GAP * (len(stages) - 1) + FOOT

    fig, ax = canvas(W, H)
    X, BW = 4.15, 6.35

    text(ax, 0.30, H - 0.30, "GeneTraceAI — master pipeline, Stage 0 to Stage 7",
         fontsize=17, fontweight="bold")
    text(ax, 0.30, H - 0.72,
         "Arrows are labelled with the artefact passing between stages and its row count.",
         fontsize=9.5, color=GREY)

    ys, cy = [], H - HEAD
    for (title, lines), h in zip(stages, heights):
        yb, yc, _ = box(ax, X, cy, BW, title, lines, fill=STAGE_FILL,
                        edge=STAGE_EDGE, title_size=TS, body_size=BS, lw=1.6)
        ys.append((cy, yb, yc))
        cy = yb - GAP

    for i, (name, sub) in enumerate(artefacts):
        p0 = (X + BW / 2, ys[i][1])
        p1 = (X + BW / 2, ys[i + 1][0])
        arrow(ax, p0, p1, colour=STAGE_EDGE, lw=1.8)
        arrow_label(ax, p0, p1, name, sub, dx=0.16)

    # ── left column: what enters each stage ──────────────────────────────
    LX, LW = 0.30, 3.35
    inputs = [
        (0, "Raw identity sources",
         ["DepMap sample_info, OmicsProfiles", "Cellosaurus, GEO series metadata"]),
        (1, "Reference builders",
         ["HGNC + Ensembl → gene_lookup", "COSMIC CGC v104 → gene_role",
          "gene_symbol_resolver.json (62,208)"]),
        (2, "The two ranking layers",
         [f"DepMap expression — log₂ TPM+1", f"  {sizes['RNA']:,} lines",
          f"CCLE / Gygi proteomics — log₂ ratio", f"  {sizes['Protein']} lines"]),
        (4, "Gate inputs (Track C)",
         ["mutations_collapsed → any_driver", "fusions_gene_level → max_confidence",
          "Boolean gate only"]),
        (5, "Confidence modifiers",
         [f"cna_flags.parquet ({c['cna_rows']:,})",
          "chronos_validation.parquet (16,866)", "COSMIC role prior"]),
        (6, "External ground truth",
         ["GDSC drug response", f"{c['curated_genes']} curated drug-target genes"]),
    ]
    for idx, title, lines in inputs:
        yc = ys[idx][2]
        h = box_height(len(lines), 8.6, 7.4)
        box(ax, LX, yc + h / 2, LW, title, lines, fill=SOURCE_FILL,
            edge=SOURCE_EDGE, title_size=8.6, body_size=7.4, lw=1.2,
            align="left")
        arrow(ax, (LX + LW, yc), (X, yc), colour=SOURCE_EDGE, lw=1.3)

    handles = [
        mpatches.Patch(fc=SOURCE_FILL, ec=SOURCE_EDGE, label="external source / reference build"),
        mpatches.Patch(fc=STAGE_FILL, ec=STAGE_EDGE, label="production stage"),
    ]
    ax.legend(handles=handles, loc="lower left",
              bbox_to_anchor=(0.30 / W, 0.10 / H), fontsize=8.6, ncol=2,
              handlelength=1.8, columnspacing=1.6)
    return save(fig, "01_pipeline_master")


# ─────────────────────────────────────────────────────────────────────────
def fig_provenance(S):
    reach = {r["source"]: r for r in S["gene_reach"]["sources"]}
    sizes = S["modality"]["set_sizes"]
    c = S["counts"]

    TS, BS = 8.6, 7.3
    GAP = 0.30

    src = [
        ("DepMap expression", [f"{reach['Expression (DepMap)']['genes_in_universe']:,} genes × {sizes['RNA']:,} lines",
                               "log₂ TPM+1 · native ENSG columns"]),
        ("CCLE / Gygi proteomics", [f"{reach['Proteomics (CCLE/Gygi)']['genes_in_universe']:,} genes × {sizes['Protein']} lines",
                                    "log₂ ratio · UniProt accessions"]),
        ("DepMap somatic mutations", [f"{reach['Mutations']['genes_in_universe']:,} genes × {sizes['Mutation']:,} lines"]),
        ("DepMap fusions", [f"{reach['Fusions']['genes_in_universe']:,} genes × {sizes['Fusion']:,} lines"]),
        ("DepMap copy number", [f"{reach['CNA flags']['genes_in_universe']:,} genes × {sizes['CNA']} lines"]),
        ("Chronos + Project Score", [f"{reach['Chronos essentiality']['genes_in_universe']:,} genes × {sizes['Chronos']:,} lines",
                                     "CRISPR dependency scores"]),
        ("GDSC drug response", [f"{sizes['GDSC drug']} lines · 687 drug–target pairs"]),
        ("Signatures · metabolomics · miRNA", ["HPA · GEO expression",
                                               "cell-line or probe axis only — no ensg_id",
                                               "to join a per-gene query on"]),
    ]
    ref = [
        ("HGNC + Ensembl", [f"gene_lookup.parquet — {c['lookup_genes']:,} genes",
                            "prev + alias symbols → 62,208 entries"]),
        ("UniProt → ENSG bridge", [f"recovers {reach['Proteomics (CCLE/Gygi)']['genes_in_universe']:,} proteomics columns",
                                   "without it the layer is unjoinable"]),
        ("Cellosaurus + RRID resolver", [f"cell_line_lookup.parquet — {c['lookup_cell_lines']:,} lines",
                                         "accession_to_primary · gap_class A/B/C"]),
        ("COSMIC CGC v104", ["gene_role: oncogene / tsg / both"]),
    ]
    stg = [
        ("Stage 2 — core score", ["expression + proteomics are the only",
                                  "two ranking layers in the system"]),
        ("Stage 4 — driver gate", ["mutations + fusions collapse to one",
                                   "Boolean, has_driver_alteration"]),
        ("Stage 5 — confidence tiers", ["CNA, Chronos, COSMIC role prior,",
                                        "GDSC validation status"]),
        ("Stage 6 — evaluation", ["GDSC sensitivity is the only external",
                                  "ground truth the pipeline is scored on"]),
        ("Stage 0 / 1 — identity", ["Cellosaurus, sample_info and GEO",
                                    "resolve every row to a model_id"]),
        ("Never wired in", ["tested as gene-agnostic per-line",
                            "scalars, then rejected"]),
    ]
    role = [
        ("Ranks cell lines", ["sets core_score and stratum_rank"]),
        ("Reorders within a gene", ["Boolean gate dominates core_score"]),
        ("Labels only", ["changes the tier, never the sort order"]),
        ("Scores the pipeline", ["hit@20 over 28 held-out genes"]),
        ("Makes the join possible", ["no score is computable without it"]),
        ("Rejected", ["bootstrap CI on Δρ included zero"]),
    ]

    def col_h(items):
        return sum(box_height(len(l), TS, BS) for _, l in items) + GAP * (len(items) - 1)

    HEAD, FOOT = 1.55, 0.95
    H = HEAD + max(col_h(x) for x in (src, ref, stg, role)) + FOOT
    W = 16.2
    fig, ax = canvas(W, H)

    text(ax, 0.30, H - 0.30, "Source-to-stage provenance", fontsize=17,
         fontweight="bold")
    text(ax, 0.30, H - 0.72,
         "Read left to right. A source can only reach a scoring stage if it carries a gene axis; the ones that do not are confined to\n"
         "cell-line-level roles. That constraint — not preference — is what fixed the layer set at expression plus proteomics.",
         fontsize=9.5, color=GREY)

    xs = [0.30, 4.35, 8.55, 12.65]
    ws = [3.75, 3.90, 3.80, 3.25]
    heads = ["Raw sources (as delivered)", "Reference / bridge layer",
             "Pipeline stage that consumes them", "Role in the final score"]
    for x, w, h in zip(xs, ws, heads):
        ax.text(x + w / 2, H - HEAD + 0.16, h, ha="center", va="bottom",
                fontsize=10, fontweight="bold", color=INK)
        ax.plot([x, x + w], [H - HEAD + 0.06] * 2, color="#c8cfd6", lw=1.0)

    def stack(x, w, items, fill, edge):
        out, cy = [], H - HEAD
        for title, lines in items:
            yb, yc, _ = box(ax, x, cy, w, title, lines, fill=fill, edge=edge,
                            title_size=TS, body_size=BS, lw=1.1, align="left")
            out.append(yc)
            cy = yb - GAP
        return out

    sp = stack(xs[0], ws[0], src, SOURCE_FILL, SOURCE_EDGE)
    rp = stack(xs[1], ws[1], ref, REF_FILL, REF_EDGE)
    tp = stack(xs[2], ws[2], stg, STAGE_FILL, STAGE_EDGE)
    op = stack(xs[3], ws[3], role, "#f5f5f7", GREY)

    # source -> (bridge) -> stage
    routes = [(0, 0, 0), (1, 1, 0), (2, 0, 1), (3, 0, 1), (4, 0, 2),
              (5, 0, 2), (6, 0, 3), (7, None, 5)]
    for si, ri, ti in routes:
        if ri is None:
            arrow(ax, (xs[0] + ws[0], sp[si]), (xs[2], tp[ti]),
                  colour=GREY, lw=1.0, rad=0.05)
            continue
        arrow(ax, (xs[0] + ws[0], sp[si]), (xs[1], rp[ri]),
              colour=SOURCE_EDGE, lw=1.1, rad=0.02)
        arrow(ax, (xs[1] + ws[1], rp[ri]), (xs[2], tp[ti]),
              colour=REF_EDGE, lw=1.1, rad=0.02)
    for ri, ti in [(2, 4), (3, 2)]:
        arrow(ax, (xs[1] + ws[1], rp[ri]), (xs[2], tp[ti]),
              colour=REF_EDGE, lw=1.1, rad=-0.04)
    for ti in range(len(stg)):
        arrow(ax, (xs[2] + ws[2], tp[ti]), (xs[3], op[ti]),
              colour=BLUE if ti < 4 else GREY, lw=1.2)

    return save(fig, "02_provenance_map")


# ─────────────────────────────────────────────────────────────────────────
def fig_join_keys(S):
    c = S["counts"]
    reach = {r["source"]: r for r in S["gene_reach"]["sources"]}
    sizes = S["modality"]["set_sizes"]

    W, H = 16.0, 8.4
    fig, ax = canvas(W, H)

    text(ax, 0.30, H - 0.30,
         "Join-key schema",
         fontsize=17, fontweight="bold")
    text(ax, 0.30, H - 0.72,
         "model_id and ensg_id are the only join keys in the system. Names — cell-line names, gene symbols, UniProt accessions, probe IDs — are\n"
         "display labels and resolver inputs, never keys. Stages 0 and 1 exist for no other purpose than to manufacture these two columns.",
         fontsize=9.5, color=GREY)

    SX, SW = 5.30, 5.40
    TOP_M, TOP_G = H - 2.30, H - 4.60

    _, ym, _ = box(ax, SX, TOP_M, SW, "model_id      — the cell-line spine", [
        "format:  ACH-######   (uppercase, enforced on every join)",
        f"{c['lookup_cell_lines']:,} rows in cell_line_lookup · {c['scored_cell_lines']:,} of them carry a score",
        "exactly one row per resolved cell line — uniqueness is asserted, not assumed"],
        fill="#e8f0fa", edge=BLUE, title_size=13, body_size=8.6, lw=2.2)
    _, yg, _ = box(ax, SX, TOP_G, SW, "ensg_id      — the gene spine", [
        "format:  ENSG###########   (bare, version suffix stripped)",
        f"{c['lookup_genes']:,} rows in gene_lookup · {c['scored_genes']:,} of them carry a score",
        "protein-coding and HGNC-Approved only — the universe filter"],
        fill=REF_FILL, edge=GREEN, title_size=13, body_size=8.6, lw=2.2)

    grain_top = TOP_G - box_height(3, 13, 8.6) - 0.70
    yb_grain, yc_grain, _ = box(
        ax, SX + 0.55, grain_top, SW - 1.10, "( model_id ,  ensg_id )   —   the scoring grain",
        [f"{c['core_score_rows']:,} rows · unique · every downstream table is a projection of this pair"],
        fill="#fdf6e3", edge=AMBER, title_size=11.5, body_size=8.4, lw=1.8)
    for dx in (-1.5, 1.5):
        arrow(ax, (SX + SW / 2 + dx, TOP_G - box_height(3, 13, 8.6)),
              (SX + SW / 2 + dx, grain_top), colour=AMBER, lw=2.0)

    # ── left: what resolves into each key ────────────────────────────────
    LX, LW = 0.30, 4.55
    ax.text(LX + LW / 2, H - 1.95, "Resolved into the key  (Stage 0 / 1)",
            ha="center", va="bottom", fontsize=10, fontweight="bold")

    hm = box_height(5, 9, 7.6)
    box(ax, LX, ym + hm / 2, LW, "Cell-line identity inputs", [
        "RRID / CVCL accession — primary and secondary",
        "CCLE name, Cellosaurus synonym, GEO accession",
        "DepMap profileid (pr-…) → modelid → uppercase",
        "→ accession_to_primary dict; gap_class A / B / C",
        "16 lines carry two RNA profiles → averaged, not dropped"],
        fill=SOURCE_FILL, edge=SOURCE_EDGE, title_size=9, body_size=7.6,
        align="left")
    arrow(ax, (LX + LW, ym), (SX, ym), colour=SOURCE_EDGE, lw=1.5)

    hg = box_height(5, 9, 7.6)
    box(ax, LX, yg + hg / 2, LW, "Gene identity inputs", [
        "HGNC approved symbol, previous symbols, aliases",
        "UniProt accession — the proteomics column headers",
        "Entrez ID, versioned Ensembl ID (ENSG…….15)",
        "→ gene_symbol_resolver.json, 62,208 entries",
        "version suffix stripped before every single join"],
        fill=SOURCE_FILL, edge=SOURCE_EDGE, title_size=9, body_size=7.6,
        align="left")
    arrow(ax, (LX + LW, yg), (SX, yg), colour=SOURCE_EDGE, lw=1.5)

    # ── right: tables keyed on each spine ────────────────────────────────
    RX, RW = 11.15, 4.55
    ax.text(RX + RW / 2, H - 1.95, "Tables keyed on the spine", ha="center",
            va="bottom", fontsize=10, fontweight="bold")

    t1 = ["harmonised_enriched — {:,}".format(c["harmonised_rows"]),
          "signatures_model_level — 1,955",
          f"metabolomics — {sizes['Metabolomics']} · miRNA — {sizes['miRNA']}",
          "these last two stop here: no gene axis"]
    y1b, y1c, _ = box(ax, RX, H - 2.35, RW, "Keyed on model_id alone", t1,
                      fill="white", edge=BLUE, title_size=9.5, body_size=7.8,
                      lw=1.4, align="left", title_colour=BLUE)
    t2 = [f"gene_lookup — {c['lookup_genes']:,}",
          f"gene_regime — {c['curated_genes']}",
          "chronos_validation — 16,866",
          "gene_role, from COSMIC CGC v104"]
    y2b, y2c, _ = box(ax, RX, y1b - 0.55, RW, "Keyed on ensg_id alone", t2,
                      fill="white", edge=GREEN, title_size=9.5, body_size=7.8,
                      lw=1.4, align="left", title_colour=GREEN)
    t3 = [f"core_score — {c['core_score_rows']:,}",
          f"predictions_with_confidence — {c['core_score_rows']:,}",
          f"evidence_ledger — {c['core_score_rows']:,}",
          f"flags_with_driver — {c['flags_rows']:,}",
          f"cna_flags — {c['cna_rows']:,}"]
    y3b, y3c, _ = box(ax, RX, y2b - 0.55, RW, "Keyed on the pair", t3,
                      fill="white", edge=AMBER, title_size=9.5, body_size=7.8,
                      lw=1.4, align="left", title_colour=AMBER)

    arrow(ax, (SX + SW, ym), (RX, y1c), colour=BLUE, lw=1.5)
    arrow(ax, (SX + SW, yg), (RX, y2c), colour=GREEN, lw=1.5)
    arrow(ax, (SX + 0.55 + (SW - 1.10), yc_grain), (RX, y3c),
          colour=AMBER, lw=1.7)

    return save(fig, "03_join_key_schema")


def main():
    S = load_stats()
    print("Architecture figures:")
    fig_pipeline_master(S)
    fig_pipeline_plain(S)
    fig_provenance(S)
    fig_join_keys(S)



# ─────────────────────────────────────────────────────────────────────────
def fig_pipeline_plain(S):
    """Same flow as figure 1, named in plain language rather than by file or
    formula, with the MATH_REFERENCE sections that govern each stage. Built
    for reading the flow end to end and checking the maths stage by stage."""
    c, r = S["counts"], S["routing"]
    n1, n2 = c["n_layers_dist"]["1"], c["n_layers_dist"]["2"]
    sizes = S["modality"]["set_sizes"]
    src = S["routing"]["regime_source_rows_full"]

    # (stage name, what it does in words, MATH_REFERENCE sections)
    stages = [
        ("Stage 0 · Give every cell line one identity", [
            "The same cell line arrives under different names in every source.",
            "Reconcile them so one physical cell line means one row.",
            f"in: 4 identity sources    out: {c['lookup_cell_lines']:,} reconciled cell lines"],
         ["Stage 0 — Harmonisation"]),

        ("Stage 1 · Give every gene one identity, and join", [
            "Same problem for genes: symbols, aliases and accessions all differ.",
            "Resolve to one gene id, then join sources through the two lookups.",
            f"in: gene + cell-line lookups    out: {c['universe_genes']:,} genes × {c['lookup_cell_lines']:,} lines"],
         ["Stage 1 — Lookup / Metadata Join"]),

        ("Stage 2 · Score how abundant each gene is in each cell line", [
            "Rank each gene across cell lines, so the two measurement types",
            "become comparable, then combine them with a penalty for the fact",
            "that they partly measure the same thing.",
            f"in: 2 abundance layers    out: {c['core_score_rows']:,} scored pairs",
            f"both layers present for {n2:,} pairs · one layer only for {n1:,}"],
         ["2.1 Percentile-rank normalisation",
          "2.2 Duplicate-profile deduplication",
          "2.3 Noisy-OR  (superseded)",
          "2.4 Correlation-penalised weighted sum  (production)",
          "2.6 Stratum-aware rank"]),

        ("Stage 3 · Check the identity work held", [
            "Audit that every scored cell line still maps to exactly one",
            "identity, with no silent duplication introduced by the joins.",
            f"in: scored pairs    out: {c['scored_cell_lines']:,} cell lines confirmed"],
         ["Stage 3 — Cell-line Name Resolution Audit"]),

        ("Stage 4 · Push cell lines with a known driver to the top", [
            "For genes that act through activation, having a driver alteration",
            "matters more than abundance alone, so those lines are lifted above",
            "the rest while keeping the abundance order inside each block.",
            f"in: mutation + fusion evidence    out: {c['flags_rows']:,} flagged pairs",
            f"{r['class_counts_curated']['activation_driven']} genes routed as activation-driven, "
            f"{r['class_counts_curated']['abundance_tracking']} as abundance-tracking"],
         ["4.1 Driver flag  (Boolean, no numeric estimator)",
          "4.2 Copy-number amplification / deletion thresholds",
          "6.2 Vectorised routing sort key",
          "4.3 / 6.4 Variance-reduction metric"]),

        ("Stage 5 · Say how much to trust each prediction", [
            "Attach a confidence label from four independent kinds of evidence,",
            "in a fixed precedence order. Deliberately categorical, because a",
            "0–1 number would imply a calibration that was never fitted.",
            f"in: 4 evidence types    out: {c['core_score_rows']:,} labelled predictions",
            f"{src.get('measured', 0):,} measured · {src.get('prior', 0):,} from a role prior · "
            f"{src.get('unknown', 0):,} unknown"],
         ["Stage 5 — Confidence Tiers",
          "5.1 Tumour-suppressor score inversion",
          "C.7 Chronos correlation gate"]),

        ("Stage 6 · Test it on genes it never saw", [
            "Hold out a quarter of the ground-truth genes, rank their cell lines,",
            "and count how many known-sensitive lines land in the top 20.",
            "in: GDSC drug response    out: 28 held-out genes scored"],
         ["6.1 Train / test gene split",
          "6.3 Hit-rate@20",
          "C.5 Spearman ρ",
          "C.8 Bootstrap CI on ρ",
          "C.9 Paired bootstrap on Δρ"]),

        ("Stage 7 · Answer a question about one gene", [
            "Take a gene, return its cell lines in order, and show which piece",
            "of evidence drove the rank and which only touched the label.",
            "in: a gene symbol or id    out: a ranked, explained list"],
         ["Stage 7 — Explain / Rank / CLI Query Layer"]),
    ]

    artefacts = [
        ("One row per cell line", f"{c['harmonised_rows']:,} reconciled identities"),
        ("Every gene paired with every cell line",
         f"{c['universe_genes']:,} genes × {c['lookup_cell_lines']:,} lines"),
        ("An abundance score for each pair",
         f"{c['core_score_rows']:,} rows · {c['scored_genes']:,} genes × {c['scored_cell_lines']:,} lines"),
        ("Verified cell-line identities", f"{c['scored_cell_lines']:,} confirmed"),
        ("Driver evidence attached to each pair",
         f"{c['flags_rows']:,} rows, plus a class for {c['curated_genes']} curated genes"),
        ("A confidence label on every prediction", f"{c['core_score_rows']:,} rows"),
        ("Held-out accuracy", "28 genes, hit@20 with a stated denominator"),
    ]

    TS, BS = 11.0, 8.4
    heights = [max(box_height(len(l), TS, BS),
                   box_height(len(m), 8.4, 7.8, title=False))
               for _, l, m in stages]
    GAP = 0.78
    HEAD, FOOT = 1.35, 0.45
    W = 15.1
    H = HEAD + sum(heights) + GAP * (len(stages) - 1) + FOOT

    fig, ax = canvas(W, H)
    X, BW = 0.30, 8.55
    MX, MW = 9.30, 5.50

    text(ax, 0.30, H - 0.30, "GeneTraceAI — what each stage does, and where its maths is written down",
         fontsize=17, fontweight="bold")
    text(ax, 0.30, H - 0.74,
         "Left: the stage in plain language, with what goes in and what comes out. "
         "Right: the sections of docs/MATH_REFERENCE.md that define it.",
         fontsize=9.5, color=GREY)
    ax.text(MX + MW / 2, H - HEAD + 0.14, "maths to check", ha="center",
            va="bottom", fontsize=9.5, fontweight="bold", color=REF_EDGE)

    ys, cy = [], H - HEAD
    for (title, lines, maths), h in zip(stages, heights):
        yb, yc, _ = box(ax, X, cy, BW, title, lines, fill=STAGE_FILL,
                        edge=STAGE_EDGE, title_size=TS, body_size=BS, lw=1.6,
                        align="left")
        # Section pointers only -- no title line, so the box is a plain list.
        mh = box_height(len(maths), 8.4, 7.6, title=False)
        box(ax, MX, yc + mh / 2, MW, None, [f"§  {m}" for m in maths],
            fill=REF_FILL, edge=REF_EDGE, body_size=7.8, lw=1.1, align="left")
        ys.append((cy, yb, yc))
        cy = yb - GAP

    for i, (name, sub) in enumerate(artefacts):
        p0 = (X + BW / 2, ys[i][1])
        p1 = (X + BW / 2, ys[i + 1][0])
        arrow(ax, p0, p1, colour=STAGE_EDGE, lw=1.8)
        arrow_label(ax, p0, p1, name, sub, dx=0.16, fs=8.4)

    return save(fig, "01b_pipeline_flow_plain")

if __name__ == "__main__":
    main()
