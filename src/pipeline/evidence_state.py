"""
evidence_state.py
-----------------
Decides whether a ranking should be shown at all, and says so outright when it
should not.

WHY THIS EXISTS
    core_score is a WITHIN-GENE percentile. Every gene is rescaled onto the same
    0-1 spread, so a housekeeping gene whose expression barely moves across cell
    lines produces a top-10 that reads exactly as confidently as a sharply
    tissue-restricted gene's. GAPDH's best line scores 0.9989 with its top line
    only 2.8x above a typical one; ERBB2's best scores 0.9975 at 76x. The numbers
    are indistinguishable; the biology is not.

    The pipeline already measures the difference (gene_dispersion.signal_spread)
    and the CLI already prints a warning. A warning under a confident-looking
    table is not enough: the reader still sees a ranked list and a 0.9989. This
    module turns that measurement into a STATE that suppresses or qualifies the
    ranking, so "there is nothing to rank here" is the headline rather than a
    footnote.

DESIGN PROVENANCE
    Follows CellLineSelector - Confidence Engine: Design Record (v1):

      C1  ordinal bands (High/Moderate/Low) PLUS distinct non-ordinal states,
          rather than collapsing everything onto one axis. Design Record uses
          `Conflicting` for source disagreement; this module adds NO_EVIDENCE and
          UNINFORMATIVE for the two "cannot answer" cases, which are categorically
          different from "answer is Low".

      C7  role separation -- the verdict is reported ALONGSIDE the ranking and is
          never summed into core_score. Nothing here changes a score.

      §5  "Missing modalities widen uncertainty, never lower the band
          (missing != against)." A gene with no protein layer is NOT low
          confidence for that reason; it is less corroborated. Absence of
          evidence is reported as absence, never as evidence against.

      L4  the Design Record explicitly flags that the reliability proxy diverges
          "in a patterned way at ubiquitously-high and flat-profile genes".
          UNINFORMATIVE is that population, named and handled.

    Evidence is TALLIED, never averaged (C3): each state below is triggered by a
    named, independently checkable fact, not by a weighted blend.

STATES
    NO_EVIDENCE    nothing measured -- no layer covers this gene (or this pair).
                   Ranking is not shown. This is the literal "no info" case.
    UNINFORMATIVE  measured, but the gene does not discriminate between cell
                   lines. A ranking exists arithmetically and means nothing.
                   Ranking is shown only behind an explicit refusal banner.
    CONFLICTING    sources disagree (abundance ranking runs opposite to measured
                   CRISPR dependency). Routed to its own state, not averaged away.
    RANKED         the ordering is worth acting on; the existing ordinal
                   confidence tier (high/moderate/low) applies per row.

THRESHOLDS
    Data-derived, and documented as design choices rather than literature values
    -- the same stance the Design Record takes in §6. `signal_spread` is computed
    in 02_core_score.ipynb Cell 7c as p99-p50 of raw log2(TPM+1), banded against
    genes that are actually expressed. This module consumes that banding; it does
    not invent a second one.

    `frac_expressed` / SILENT_FRAC (Cell 7c-silence, added per
    docs/TRANSCRIPTOMICS_DECISIONS.md Q4) is a second, independent flat-signal
    check and is NOT redundant with signal_spread. Measured
    (docs/audit_scripts/d_f2.py PART F2): of 5,847 genes failing the silence
    guard, signal_spread calls 2,420 of them "wide" and 2,091 exceed
    FLAT_FOLD_CEILING -- a dispersion check does not merely miss this
    population, it actively certifies it as well-dispersed. Both checks OR
    together into UNINFORMATIVE.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

# A gene whose top line is under this fold-change above a typical line cannot
# meaningfully separate cell lines, whatever its percentiles look like. Set to
# agree with gene_dispersion's own "narrow" band rather than as a free parameter.
FLAT_FOLD_CEILING = 3.0

# A gene expressed (log2TPM > EXPRESSED_MIN) in fewer than this fraction of
# scored lines cannot support a meaningful ranking, however wide its
# dynamic_range -- signal_spread alone misses most of this population (see
# THRESHOLDS above). Values match the transcriptomics notebook's own constants
# (02_transcriptonomics.ipynb cell 1 / docs/audit_scripts/tx_lib.py) rather
# than being re-derived here.
EXPRESSED_MIN = 1.0
SILENT_FRAC = 0.20

# Below this many scored lines, a within-gene percentile is not a stable ranking.
MIN_LINES_FOR_RANKING = 30

NO_EVIDENCE = "NO_EVIDENCE"
UNINFORMATIVE = "UNINFORMATIVE"
CONFLICTING = "CONFLICTING"
RANKED = "RANKED"


@dataclass
class Verdict:
    """What the architecture is willing to claim about one gene."""
    state: str
    headline: str                       # one line, said outright
    detail: list[str] = field(default_factory=list)
    show_ranking: bool = True           # False => suppress the table entirely
    qualify_ranking: bool = False       # True  => show, but behind a refusal banner
    facts: dict = field(default_factory=dict)

    @property
    def is_answerable(self) -> bool:
        return self.state == RANKED


def gene_verdict(ensg_id, symbol, n_candidates, dispersion_row=None,
                 contradicted=False, gene_class=None):
    """
    Classify a gene BEFORE any ranking is presented.

    dispersion_row: the gene's row of gene_dispersion.parquet (or None if the
    gene has no expression layer at all -- which is itself informative).
    """
    facts = {"ensg_id": ensg_id, "symbol": symbol, "n_candidates": int(n_candidates or 0)}

    # ---- 1. no evidence at all ------------------------------------------
    if not n_candidates:
        return Verdict(
            NO_EVIDENCE,
            f"NO INFORMATION: no omics layer covers {symbol} in any cell line.",
            ["Nothing is scored for this gene, so there is no ranking to give.",
             "This is an absence of measurement, not a low result -- the pipeline "
             "is not saying this gene is unimportant, it is saying it has no data."],
            show_ranking=False, facts=facts)

    if n_candidates < MIN_LINES_FOR_RANKING:
        return Verdict(
            NO_EVIDENCE,
            f"TOO FEW LINES: {symbol} is scored in only {n_candidates} cell "
            f"line(s) -- below the {MIN_LINES_FOR_RANKING} needed for a percentile.",
            ["core_score is a rank across cell lines. With this few, the rank is "
             "an artefact of who happened to be measured.",
             "Reported as unanswerable rather than as a short ranking."],
            show_ranking=False, facts=facts)

    # ---- 2. measured, but flat ------------------------------------------
    spread = fold = median_tpm = frac_expressed = None
    if dispersion_row is not None and len(dispersion_row):
        r = dispersion_row.iloc[0] if hasattr(dispersion_row, "iloc") else dispersion_row
        spread = r.get("signal_spread")
        fold = r.get("top_vs_median_fold")
        median_tpm = r.get("median_log2tpm")
        frac_expressed = r.get("frac_expressed")
        facts.update({"signal_spread": spread, "top_vs_median_fold": fold,
                      "median_log2tpm": median_tpm, "frac_expressed": frac_expressed})

    # Defer to gene_dispersion's OWN band. It is derived in 02_core_score Cell 7c
    # against the distribution of genes that are actually expressed, which is a
    # far better reference than any absolute fold-change this module could pick:
    # BRAF sits at 2.9x and is banded 'moderate', GAPDH at 2.8x and 'narrow'.
    # A flat cutoff near 3 would refuse to rank BRAF, the project's own anchor
    # gene. The absolute ceiling is therefore a FALLBACK for when the band is
    # missing, not a second gate.
    is_silent = frac_expressed is not None and frac_expressed < SILENT_FRAC
    if spread is not None:
        is_flat = (spread == "narrow") or is_silent
    else:
        is_flat = (fold is not None and fold < FLAT_FOLD_CEILING) or is_silent

    if is_flat:
        fold_s = f"{fold:.1f}x" if fold is not None else "n/a"
        frac_s = f"{100*frac_expressed:.1f}%" if frac_expressed is not None else "n/a"

        # An activation_driven gene is not ranked on abundance in the first place
        # -- Stage 4 sorts it on the driver flag, with core_score only as support.
        # A flat expression profile therefore does not invalidate its ranking the
        # way it does for an abundance_tracking gene. Say that, rather than
        # refusing a ranking that never rested on the flat quantity.
        if gene_class == "activation_driven":
            return Verdict(
                RANKED,
                f"{symbol}: flat abundance ({fold_s}), but this gene is ranked on "
                "driver alterations, not abundance.",
                ["Expression barely varies across cell lines, so core_score alone "
                 "would not separate them.",
                 "The ordering below is driven by the mutation/fusion driver flag; "
                 "core_score only breaks ties within it."],
                show_ranking=True, facts=facts)

        detail = [
            f"The best cell line expresses {symbol} only {fold_s} above a typical one. "
            "Across cell lines this gene is effectively constant.",
        ]
        if is_silent:
            detail.append(
                f"Expressed (log2TPM > {EXPRESSED_MIN}) in only {frac_s} of scored "
                f"lines -- below the {100*SILENT_FRAC:.0f}% silence-guard floor. A wide "
                "dynamic_range does not rule this out: a gene can look well-dispersed by "
                "p99-p50 and still be silent almost everywhere, with signal confined to "
                "a handful of lines that a spread statistic alone would certify as fine.")
        detail += [
            "core_score is a within-gene percentile, so it still spans 0-1 and the top "
            "row still reads ~0.99. That number reflects the rescaling, not a real gap.",
            "Any ordering below is dominated by measurement noise. Do not select cell "
            "lines on it.",
        ]
        if median_tpm is not None and median_tpm > 8:
            detail.append(
                f"Median expression is high (log2TPM {median_tpm:.1f}) and near-uniform -- "
                "the signature of a housekeeping gene.")
        return Verdict(
            UNINFORMATIVE,
            f"NO DISCRIMINATING SIGNAL: {symbol} does not separate cell lines.",
            detail, show_ranking=True, qualify_ranking=True, facts=facts)

    # ---- 3. sources disagree --------------------------------------------
    if contradicted:
        return Verdict(
            CONFLICTING,
            f"CONFLICTING EVIDENCE: for {symbol}, abundance ranking runs OPPOSITE "
            "to measured CRISPR dependency.",
            ["The lines this ranking puts on top are, if anything, the ones least "
             "dependent on the gene.",
             "Reported as its own state rather than folded into a lower tier -- a "
             "contradiction is not the same as weak support."],
            show_ranking=True, qualify_ranking=True, facts=facts)

    # ---- 4. answerable ---------------------------------------------------
    fold_s = f"{fold:.1f}x" if fold is not None else "n/a"
    return Verdict(
        RANKED,
        f"{symbol}: ranking is interpretable ({n_candidates:,} cell lines, "
        f"top line {fold_s} above typical).",
        [], show_ranking=True, facts=facts)


def pair_verdict(symbol, model_id, row=None, gene_v=None):
    """Classify one (gene, cell line) pair. Absence is reported as absence."""
    if row is None or (hasattr(row, "__len__") and len(row) == 0):
        return Verdict(
            NO_EVIDENCE,
            f"NO INFORMATION: no layer measures {symbol} in {model_id}.",
            ["This pair was never scored. The pipeline has nothing to say about it.",
             "Not a low score -- no score."],
            show_ranking=False)
    if gene_v is not None and gene_v.state in (UNINFORMATIVE, CONFLICTING):
        return Verdict(
            gene_v.state,
            f"{gene_v.headline} Any position given for {model_id} inherits that.",
            gene_v.detail, show_ranking=True, qualify_ranking=True)
    return Verdict(RANKED, "", [], show_ranking=True)


def load_dispersion(base="src/pipeline/outputs"):
    """gene_dispersion, keyed on lowercase ensg_id. Returns None if absent."""
    try:
        d = pd.read_parquet(f"{base}/gene_dispersion.parquet")
    except (FileNotFoundError, OSError):
        return None
    d["ensg_id"] = d["ensg_id"].astype("string").str.lower()
    return d.set_index("ensg_id")


def render(v: Verdict, indent="  "):
    """One concise line for the CLI: [STATE] headline. Full detail/facts
    remain on the Verdict object for callers that want them (export_web.py
    serialises v.detail into JSON as-is; this is terminal-only)."""
    if v.state == RANKED and not v.detail:
        return ""
    tag = f"[{v.state}] " if v.state != RANKED else ""
    return f"{tag}{v.headline}"
