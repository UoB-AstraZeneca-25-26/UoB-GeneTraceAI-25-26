"""
gate_audit/common.py
--------------------
Single source of truth for the gate audit. Every script in this folder imports
its loaders, its scoring and its metrics from here, so the five audits cannot
silently diverge from each other or from the original
`src/pipeline/test_run_abundance_gate.py`.

Nothing in this folder writes to `src/pipeline/outputs/`. All output goes to
`gate_audit/outputs/`.

WHAT IS REPRODUCED FROM THE ORIGINAL, EXACTLY
---------------------------------------------
  * van der Waerden layer transform          -- `rank(method="min")`, norm.ppf
  * Stouffer combination with rho=0.46, w_E=0.987, w_P=0.381
  * McClish-standardised partial AUROC
  * decile assignment                        -- `rankdata(..., "ordinal")`
  * DEP_THRESHOLD = -0.5, MIN_LINES = 100, MIN_POS/MIN_NEG = 10

WHAT IS NEW, AND WHY IT IS PARAMETERISED
----------------------------------------
The original hard-codes two choices that the audit exists to interrogate:

  1. TIE HANDLING. `vdw()` uses `rank(method="min")`, so every line in a tie
     block receives the SAME z. `decile_index()` then uses
     `rankdata(method="ordinal")`, which breaks those ties by array order and
     therefore SPREADS one tie block across several deciles. That is the
     "ordinal / arbitrary tiebreak" branch of the concern this audit tests.
     `decile_index(tie_mode=...)` exposes all three behaviours:

        "ordinal"  reproduce the headline exactly (arbitrary spread)
        "average"  a whole tie block lands in the decile of its mean rank
        "drop"     tied lines are excluded from the decile profile entirely

  2. GATE REGION. `partial_auc(lo, hi)` was fixed at [0.8, 1.0]. Script 04
     sweeps it.

GENE VALIDITY
-------------
`docs/TRANSCRIPTOMICS_SPEC.md` F2 defines it, and this module implements that
definition unchanged:

    frac_expressed(g) = |{lines with log2(TPM+1) > EXPRESSED_MIN}| / |lines|
    gene is INVALID (UNINFORMATIVE) when frac_expressed(g) < SILENT_FRAC

with EXPRESSED_MIN = 1.0 and SILENT_FRAC = 0.20, which is the cut that produced
5,847 of 19,173 genes flagged (30.5%) and therefore 13,326 valid genes.
`frac_expressed` is computed over the FULL DepMap expression panel, not over the
label-restricted subset, so the numbers reconcile with that document.

SIGN CONVENTIONS
----------------
  score        higher = higher abundance = stronger candidate
  essentiality more NEGATIVE = more essential (checked, not assumed)
  positive     essentiality <= DEP_THRESHOLD
  ratio < 1 in decile 1 = DEPLETION = what the gate hypothesis predicts
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import norm, rankdata

warnings.filterwarnings("ignore", category=RuntimeWarning)

# --------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parent.parent.parent
AUDIT = Path(__file__).resolve().parent
OUT = AUDIT / "outputs"
OUT.mkdir(exist_ok=True)

REF = ROOT / "reference"
DATA_CLEAN = ROOT / "data" / "parquet" / "data_clean"
CLEANED_TRACK = ROOT / "cleaned_track_data"
PIPE_OUT = ROOT / "architecture" / "outputs"
CHRONOS_DIR = ROOT / "data" / "DepMap_Chronos"

CHRONOS_FILES = {
    "achilles": CHRONOS_DIR / "GeneFitnessEffect_Chronos_Achilles.hdf5",   # Broad
    "score": CHRONOS_DIR / "GeneFitnessEffect_Chronos_Score.hdf5",         # Sanger
}

PANEL = OUT / "panel.parquet"
PANEL_META = OUT / "panel_meta.json"

# --------------------------------------------------------------- constants
RHO_GLOBAL = 0.46
W_E, W_P = 0.987, 0.381
DEP_THRESHOLD = -0.5
MIN_LINES = 100
MIN_POS = 10
MIN_NEG = 10
SEED = 42

EXPRESSED_MIN = 1.0          # log2(TPM+1); TPM >= 1 detection convention
SILENT_FRAC = 0.20           # the cut that yields 13,326 valid genes

GATE_REGION = (0.8, 1.0)     # the headline region, swept in 04
RANKER_REGION = (0.0, 0.2)

# Negative-control families: not expressed in cancer lines, not essential.
# KRTAP / VN1R are excluded from FPR estimates (multi-targeting paralogy), as in
# src/pipeline/test_run_gate_mechanism.py.
CONTROL_FAMILIES_CLEAN = (
    ("OR", lambda s: bool(re.fullmatch(r"OR\d+[A-Z]\d*", s))),
    ("TAS2R", lambda s: s.startswith("TAS2R")),
)


# =============================================================== math
def vdw(s: pd.Series, tie_method: str = "min") -> pd.Series:
    """van der Waerden / inverse-normal score. `tie_method="min"` reproduces the
    original exactly; tied values then share one z."""
    n = int(s.notna().sum())
    if n == 0:
        return pd.Series(dtype=float, index=s.index)
    return pd.Series(norm.ppf(s.rank(method=tie_method) / (n + 1)), index=s.index)


def stouffer(zE: pd.Series, zP: pd.Series,
             rho: float = RHO_GLOBAL, wE: float = W_E, wP: float = W_P) -> pd.Series:
    out = np.full(len(zE), np.nan)
    e_ok, p_ok = zE.notna().values, zP.notna().values
    both = e_ok & p_ok
    ze, zp = zE.values, zP.values
    out[e_ok & ~p_ok] = ze[e_ok & ~p_ok]
    out[p_ok & ~e_ok] = zp[p_ok & ~e_ok]
    out[both] = (wE * ze[both] + wP * zp[both]) / np.sqrt(
        wE ** 2 + wP ** 2 + 2 * wE * wP * rho)
    return pd.Series(out, index=zE.index)


def decile_index(score: np.ndarray, n_bins: int = 10, tie_mode: str = "ordinal"):
    """Bin `score` into `n_bins` equal-count strata. Returns (bin_index, keep_mask).

    tie_mode
      "ordinal"  rankdata(method="ordinal") -- the ORIGINAL. Ties are broken by
                 array order, so a tie block of 800 lines is smeared across
                 however many bins 800 lines happen to span. Deterministic, but
                 the assignment carries no information.
      "average"  rankdata(method="average") -- every member of a tie block gets
                 the block's mean rank and so lands in ONE bin together. A block
                 larger than a bin therefore overflows that bin; bin sizes stop
                 being equal, which is the honest consequence of the tie mass.
      "drop"     tied values are removed before binning. Bins are then computed
                 on distinct values only. Loses lines, but every remaining
                 assignment is meaningful.
    """
    score = np.asarray(score, dtype=float)
    keep = ~np.isnan(score)
    if tie_mode == "drop":
        vals, counts = np.unique(score[keep], return_counts=True)
        tied_vals = set(vals[counts > 1].tolist())
        keep = keep & ~np.isin(score, list(tied_vals) or [np.inf])
    s = score[keep]
    n = len(s)
    if n == 0:
        return np.full(len(score), -1), keep
    method = "average" if tie_mode == "average" else "ordinal"
    r = rankdata(s, method=method)
    q = np.minimum(((np.ceil(r) - 1) * n_bins // n).astype(int), n_bins - 1)
    out = np.full(len(score), -1)
    out[keep] = q
    return out, keep


def roc_curve(score: np.ndarray, pos: np.ndarray):
    """FPR, TPR. Equal scores collapse into a single ROC step, so tie blocks are
    handled correctly here regardless of `decile_index`'s tie_mode."""
    order = np.argsort(-score, kind="mergesort")
    s, p = score[order], pos[order]
    P, N = p.sum(), (~p).sum()
    distinct = np.r_[np.diff(s) != 0, True]
    tp = np.cumsum(p)[distinct]
    fp = np.cumsum(~p)[distinct]
    return np.r_[0, fp / N], np.r_[0, tp / P]


def partial_auc(score: np.ndarray, pos: np.ndarray, lo: float, hi: float) -> float:
    """McClish (1989) standardised partial AUROC over FPR in [lo, hi]. 0.5 = chance."""
    ok = ~np.isnan(score)
    score, pos = score[ok], pos[ok]
    if pos.sum() < 2 or (~pos).sum() < 2:
        return np.nan
    fpr, tpr = roc_curve(score, pos)
    grid = np.linspace(lo, hi, 501)
    t = np.interp(grid, fpr, tpr)
    area = np.trapezoid(t, grid)
    a_min = (hi ** 2 - lo ** 2) / 2.0
    a_max = hi - lo
    return 0.5 * (1.0 + (area - a_min) / (a_max - a_min))


def auroc(score: np.ndarray, pos: np.ndarray) -> float:
    ok = ~np.isnan(score)
    score, pos = score[ok], pos[ok]
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 < 1 or n0 < 1:
        return np.nan
    r = rankdata(score)
    return (r[pos].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def tie_stats(v: np.ndarray) -> dict:
    """Tie structure of one gene's value vector."""
    v = np.asarray(v, dtype=float)
    v = v[~np.isnan(v)]
    n = len(v)
    if n == 0:
        return {"n": 0}
    vals, counts = np.unique(v, return_counts=True)
    tied = counts[counts > 1]
    biggest = int(counts.max())
    return {
        "n": n,
        "n_distinct": int(len(vals)),
        "tie_mass": float(tied.sum() / n),              # share of lines in ANY tie
        "largest_block": biggest,
        "largest_block_frac": float(biggest / n),
        "floor_block": int((v == v.min()).sum()),
        "floor_block_frac": float((v == v.min()).mean()),
        "n_tie_blocks": int((counts > 1).sum()),
    }


# =============================================================== loaders
def load_gene_lookup():
    """protein-coding + HGNC-Approved universe, plus symbol and UniProt maps."""
    gl = pd.read_parquet(REF / "gene_lookup.parquet",
                         columns=["ensg_id", "hgnc_symbol", "biotype",
                                  "hgnc_status", "uniprot_ids"])
    gl["ensg_id"] = gl.ensg_id.astype("string").str.split(".").str[0].str.lower()
    uni = gl[(gl.biotype == "protein_coding") & (gl.hgnc_status == "Approved")].copy()
    symbol_of = dict(zip(uni.ensg_id, uni.hgnc_symbol.astype(str).str.upper()))
    sym2ensg = {s: e for e, s in symbol_of.items()}
    u2e = {}
    for e, u in zip(uni.ensg_id, uni.uniprot_ids):
        if isinstance(u, str):
            for acc in u.split("|"):
                acc = acc.strip().lower()
                if acc:
                    u2e.setdefault(acc, e)
    return uni, set(uni.ensg_id), symbol_of, sym2ensg, u2e


def load_chronos(which: str, sym2ensg: dict, verbose: bool = True):
    """Real DepMap Chronos gene effect. Returns (wide DataFrame lines x ensg,
    scale-check dict). Aborts if the scale check fails -- the check whose absence
    let a negated Project Score matrix pass as Chronos project-wide.

    which: "achilles" (Broad) or "score" (Sanger).
    """
    import h5py
    cp = CHRONOS_FILES[which]
    if not cp.exists():
        raise SystemExit(f"{cp} not found.")
    with h5py.File(cp, "r") as fh:
        D = fh["data"][:]
        lines_h = [x.decode() for x in fh["dim_0"][:]]
        genes_h = [x.decode() for x in fh["dim_1"][:]]
    keep_i, keep_e = [], []
    for i, c in enumerate(genes_h):
        e = sym2ensg.get(c.split(" (")[0].strip().upper())
        if e:
            keep_i.append(i)
            keep_e.append(e)
    W = pd.DataFrame(D[:, keep_i], index=[m.lower() for m in lines_h], columns=keep_e)
    W = W.T.groupby(level=0).mean().T

    ess = set(pd.read_csv(CHRONOS_DIR / "ReferenceEssentials.csv").iloc[:, 0]
              .astype(str).str.split(" (", regex=False).str[0].str.upper())
    non = set(pd.read_csv(CHRONOS_DIR / "ReferenceNonEssentials.csv").iloc[:, 0]
              .astype(str).str.split(" (", regex=False).str[0].str.upper())
    e_ess = sorted({sym2ensg[g] for g in ess if g in sym2ensg} & set(W.columns))
    e_non = sorted({sym2ensg[g] for g in non if g in sym2ensg} & set(W.columns))
    chk = {
        "file": cp.name,
        "n_lines": int(W.shape[0]),
        "n_genes": int(W.shape[1]),
        "median_all": float(np.nanmedian(W.values)),
        "median_essential": float(np.nanmedian(W[e_ess].values)),
        "median_nonessential": float(np.nanmedian(W[e_non].values)),
        "label_fpr": float(np.nanmean(W[e_non].values <= DEP_THRESHOLD)),
        "n_ref_nonessential": len(e_non),
    }
    if verbose:
        print(f"  chronos[{which}] {cp.name}: {W.shape[0]} lines x {W.shape[1]:,} genes")
        print(f"    scale check: median {chk['median_all']:+.4f} (expect ~0); "
              f"essentials {chk['median_essential']:+.4f} (~-1); "
              f"non-essentials {chk['median_nonessential']:+.4f} (~0)")
        print(f"    LABEL FPR at {DEP_THRESHOLD}: {chk['label_fpr']:.4f} "
              f"({len(e_non)} curated non-essentials)")
    if not (chk["median_essential"] < -0.5 < chk["median_nonessential"] + 0.5):
        raise SystemExit(f"Chronos scale check FAILED for {which}. Aborting.")
    return W, chk


def _wide_from_symbol_matrix(df, sym2ensg, index_name="model"):
    """DepMap 24Q4 and DEMETER2 matrices label columns 'SYMBOL (ENTREZ)'. Map to
    ensg, drop unmapped, and average any ensg reached by two symbols."""
    cols, keep = [], []
    for c in df.columns:
        e = sym2ensg.get(str(c).split(" (")[0].strip().upper())
        if e:
            keep.append(c)
            cols.append(e)
    W = df[keep]
    W.columns = cols
    W = W.T.groupby(level=0).mean().T
    W.index = [str(m).lower() for m in W.index]
    return W


def _scale_check(W, sym2ensg, tag, verbose=True, mode="absolute",
                 target_fpr=None):
    """Guard against a mis-signed or mis-scaled label entering the audit.

    mode="absolute"     require essentials < -0.5 < non-essentials + 0.5. Correct
                        for Chronos-scaled CRISPR, where a gene effect of -1 is
                        the common-essential anchor.
    mode="directional"  require only that essentials sit clearly below
                        non-essentials, and CALIBRATE a threshold instead of
                        assuming -0.5 transfers.

    Why the second mode exists: RNAi knockdown is PARTIAL where CRISPR knockout
    is complete, so DEMETER2 effect sizes are compressed -- curated essentials
    land near -0.35, not -1. Applying DEP_THRESHOLD = -0.5 to DEMETER2 would not
    be a strict threshold, it would be a DIFFERENT threshold, and every
    prevalence would silently shift. The scale is bridged instead by matching the
    label's false-positive rate on curated NON-essentials to `target_fpr`, which
    is the one quantity the two assays should agree on: how often a gene known
    not to be essential is nonetheless called dependent.
    """
    ess = set(pd.read_csv(CHRONOS_DIR / "ReferenceEssentials.csv").iloc[:, 0]
              .astype(str).str.split(" (", regex=False).str[0].str.upper())
    non = set(pd.read_csv(CHRONOS_DIR / "ReferenceNonEssentials.csv").iloc[:, 0]
              .astype(str).str.split(" (", regex=False).str[0].str.upper())
    e_ess = sorted({sym2ensg[g] for g in ess if g in sym2ensg} & set(W.columns))
    e_non = sorted({sym2ensg[g] for g in non if g in sym2ensg} & set(W.columns))
    chk = {"tag": tag, "n_lines": int(W.shape[0]), "n_genes": int(W.shape[1]),
           "median_all": float(np.nanmedian(W.values)),
           "median_essential": float(np.nanmedian(W[e_ess].values)),
           "median_nonessential": float(np.nanmedian(W[e_non].values)),
           "label_fpr": float(np.nanmean(W[e_non].values <= DEP_THRESHOLD)),
           "n_ref_nonessential": len(e_non)}
    if verbose:
        print(f"  {tag}: {W.shape[0]} lines x {W.shape[1]:,} genes")
        print(f"    scale check: median {chk['median_all']:+.4f} (expect ~0); "
              f"essentials {chk['median_essential']:+.4f}; "
              f"non-essentials {chk['median_nonessential']:+.4f}")

    if mode == "absolute":
        if not (chk["median_essential"] < -0.5 < chk["median_nonessential"] + 0.5):
            raise SystemExit(f"scale check FAILED for {tag}. Aborting.")
        chk["threshold"] = DEP_THRESHOLD
        if verbose:
            print(f"    LABEL FPR at {DEP_THRESHOLD}: {chk['label_fpr']:.4f}")
        return chk

    # directional: separation must be real, but the cut is calibrated not assumed
    sep = chk["median_nonessential"] - chk["median_essential"]
    chk["separation"] = float(sep)
    if not (chk["median_essential"] < chk["median_nonessential"] and sep > 0.15):
        raise SystemExit(
            f"scale check FAILED for {tag}: reference sets do not separate "
            f"(essentials {chk['median_essential']:+.4f} vs non-essentials "
            f"{chk['median_nonessential']:+.4f}). Aborting.")
    v = W[e_non].values
    v = v[~np.isnan(v)]
    thr = float(np.quantile(v, target_fpr)) if target_fpr else DEP_THRESHOLD
    chk["threshold"] = thr
    chk["threshold_target_fpr"] = target_fpr
    chk["label_fpr_at_threshold"] = float(np.nanmean(W[e_non].values <= thr))
    if verbose:
        print(f"    separation (non-ess - ess): {sep:+.4f}")
        print(f"    DEP_THRESHOLD {DEP_THRESHOLD} does NOT transfer to this scale")
        print(f"    calibrated threshold {thr:+.4f} "
              f"(matches non-essential FPR {chk['label_fpr_at_threshold']:.4f} "
              f"to target {target_fpr})")
    return chk


DEPMAP_24Q4 = ROOT / "data" / "DepMap_24Q4"
DEMETER2_DIR = ROOT / "data" / "DEMETER2"


def load_crispr_24q4(sym2ensg, verbose=True):
    """DepMap 24Q4 CRISPRGeneEffect (Chronos). Rows are ModelID directly, so no
    profile collapse is needed. Replaces the 20Q2-vintage Zenodo bundle."""
    f = DEPMAP_24Q4 / "CRISPRGeneEffect.csv"
    if not f.exists():
        raise SystemExit(f"{f} not found.")
    df = pd.read_csv(f, index_col=0, low_memory=False)
    W = _wide_from_symbol_matrix(df, sym2ensg)
    return W, _scale_check(W, sym2ensg, "CRISPR 24Q4 (Chronos)", verbose)


def load_expression_24q4(sym2ensg, verbose=True):
    """DepMap 24Q4 OmicsExpressionProteinCodingGenesTPMLogp1, log2(TPM+1), keyed
    by ModelID. Same units as the 20Q2 matrix the pipeline uses, so
    EXPRESSED_MIN = 1.0 carries over unchanged."""
    f = DEPMAP_24Q4 / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    if not f.exists():
        raise SystemExit(f"{f} not found.")
    df = pd.read_csv(f, index_col=0, low_memory=False)
    E = _wide_from_symbol_matrix(df, sym2ensg)
    if verbose:
        print(f"  expression 24Q4: {E.shape[0]} lines x {E.shape[1]:,} genes")
    return E


def load_demeter2(sym2ensg, verbose=True, target_fpr=0.0063):
    """DEMETER2 combined RNAi gene dependency (Achilles shRNA + DRIVE +
    Marcotte). A THIRD perturbation type: knockdown, not knockout, not drug.

    The matrix is genes x CCLE-name; it is transposed and bridged to ModelID via
    DepMap 24Q4 Model.csv `CCLEName`. Lines that do not bridge are DROPPED, never
    defaulted -- the audit-3 rule.
    """
    f = DEMETER2_DIR / "D2_combined_gene_dep_scores.csv"
    if not f.exists():
        raise SystemExit(f"{f} not found.")
    df = pd.read_csv(f, index_col=0, low_memory=False).T      # -> lines x genes
    W = _wide_from_symbol_matrix(df, sym2ensg)

    mf = DEPMAP_24Q4 / "Model.csv"
    if not mf.exists():
        raise SystemExit(f"{mf} not found (needed to bridge CCLE names).")
    md = pd.read_csv(mf, low_memory=False)
    ccle_col = "CCLEName" if "CCLEName" in md.columns else "CCLE_Name"
    bridge = {str(c).strip().lower(): str(m).lower()
              for c, m in zip(md[ccle_col], md["ModelID"]) if isinstance(c, str)}
    mapped = [bridge.get(x) for x in W.index]
    ok = [i for i, v in enumerate(mapped) if v]
    if verbose:
        print(f"  DEMETER2: {W.shape[0]} lines, bridged to ModelID: {len(ok)} "
              f"({W.shape[0]-len(ok)} dropped, never defaulted)")
    W = W.iloc[ok]
    W.index = [mapped[i] for i in ok]
    W = W.groupby(level=0).mean()
    return W, _scale_check(W, sym2ensg, "DEMETER2 (RNAi)", verbose,
                           mode="directional", target_fpr=target_fpr)


def expression_columns():
    names = pq.ParquetFile(DATA_CLEAN / "depmap_expr_clean.parquet").schema_arrow.names
    return {c.split(".")[0].lower(): c for c in names if c != "index"}


def load_expression(genes, expr_col=None):
    """model_id x ensg matrix of log2(TPM+1), profiles collapsed to model by mean.
    This is the layer the gate runs on and the layer whose ties are in question."""
    expr_col = expr_col or expression_columns()
    genes = [g for g in genes if g in expr_col]
    prof = pd.read_parquet(REF / "depmap_profiles.parquet")
    rp = prof[prof.datatype == "rna"][["profileid", "modelid"]].copy()
    rp["model_id"] = rp.modelid.str.lower()
    rp["profileid"] = rp.profileid.astype(str)
    E = pq.read_table(DATA_CLEAN / "depmap_expr_clean.parquet",
                      columns=["index"] + [expr_col[g] for g in genes]).to_pandas()
    if "index" not in E.columns:
        E = E.reset_index()
    E["profileid"] = E["index"].astype(str)
    E = E.drop(columns=["index"]).merge(rp[["profileid", "model_id"]],
                                        on="profileid", how="inner")
    E = E.drop(columns=["profileid"]).groupby("model_id").mean()
    E.columns = [c.split(".")[0].lower() for c in E.columns]
    return E


def load_proteomics(genes, valid_ensg, u2e):
    """CCLE (Track C) and ProCan protein layers, keyed model_id x accession,
    plus the ensg -> accession maps. Mirrors the original loader."""
    import duckdb
    prot_names = pq.ParquetFile(CLEANED_TRACK / "proteomics.parquet").schema_arrow.names
    key_col = "model_id" if "model_id" in prot_names else "depmap_id"
    ccle_by_ensg = {}
    for c in prot_names:
        if c == key_col:
            continue
        a = c.strip().lower()
        e = u2e.get(a) or u2e.get(a.split("-")[0])
        if e in valid_ensg:
            ccle_by_ensg.setdefault(e, []).append(c)

    con = duckdb.connect(str(PIPE_OUT / "celllineselector.db"), read_only=True)
    pc_names = [r[0] for r in con.execute("DESCRIBE procan_proteomics").fetchall()]
    PC_META = {"gdsc_model_name", "sanger_model_id", "model_id",
               "matched_via", "n_model_id", "is_ambiguous"}
    procan_by_ensg = {}
    for c in pc_names:
        if c in PC_META:
            continue
        e = u2e.get(c.strip().lower())
        if e in valid_ensg:
            procan_by_ensg.setdefault(e, []).append(c)

    need_ccle = sorted({a for g in genes for a in ccle_by_ensg.get(g, [])})
    C = pq.read_table(CLEANED_TRACK / "proteomics.parquet",
                      columns=[key_col] + need_ccle).to_pandas()
    C.index = C[key_col].str.lower().values
    C = C.drop(columns=[key_col]).astype(float).groupby(level=0).mean()

    need_pc = sorted({a for g in genes for a in procan_by_ensg.get(g, [])})
    if need_pc:
        q = ", ".join(f'"{a}"' for a in need_pc)
        P = con.execute(f"SELECT model_id, {q} FROM procan_proteomics "
                        "WHERE model_id IS NOT NULL").df()
        P.index = P.model_id.str.lower().values
        P = P.drop(columns=["model_id"]).astype(float).groupby(level=0).mean()
    else:
        P = pd.DataFrame()
    con.close()
    return C, P, ccle_by_ensg, procan_by_ensg


def frac_expressed(E: pd.DataFrame) -> pd.Series:
    """docs/TRANSCRIPTOMICS_SPEC.md F2, verbatim. Computed over the whole panel."""
    return (E > EXPRESSED_MIN).sum(axis=0) / E.notna().sum(axis=0)


def load_lineage() -> pd.Series:
    """model_id -> lineage label, for the cluster bootstrap. First NCIT term of
    `diseases` in the cell-line lookup, matching
    src/pipeline/cell_line_grouping.py:disease_short_label. Lines with no
    disease term are pooled into one 'unknown' cluster, which is the
    conservative choice: it lowers the effective sample size rather than
    inflating it with singletons."""
    cl = pd.read_parquet(REF / "cell_line_lookup.parquet",
                         columns=["model_id", "diseases", "category"])

    def short(d):
        if pd.isna(d) or not d:
            return "unknown"
        first = str(d).split("||")[0].strip()
        parts = [p.strip() for p in first.split(";")]
        return parts[-1] if parts and parts[-1] else "unknown"

    cl["lineage"] = cl.diseases.map(short)
    cl["model_id"] = cl.model_id.astype(str).str.lower()
    return cl.drop_duplicates("model_id").set_index("model_id")["lineage"]


def control_fpr(Wide: pd.DataFrame, symbol_of: dict) -> dict:
    """Label false-positive rate from non-expressed, non-essential gene families.
    Paralogy-contaminated families are excluded (see test_run_gate_mechanism.py)."""
    out = {}
    keep = set()
    for name, fn in CONTROL_FAMILIES_CLEAN:
        genes = [e for e in Wide.columns
                 if symbol_of.get(e) and fn(symbol_of[e])]
        if genes:
            v = Wide[genes].values
            out[name] = {"n_genes": len(genes),
                         "dep_rate": float(np.nanmean(v <= DEP_THRESHOLD))}
            keep |= set(genes)
    if keep:
        v = Wide[sorted(keep)].values
        out["combined"] = {"n_genes": len(keep),
                           "dep_rate": float(np.nanmean(v <= DEP_THRESHOLD))}
    return out


# =============================================================== panel
def read_panel() -> pd.DataFrame:
    if not PANEL.exists():
        raise SystemExit(
            "gate_audit/outputs/panel.parquet not found.\n"
            "Run  python gate_audit/00_build_panel.py  first.")
    return pd.read_parquet(PANEL)


def gene_frame(panel: pd.DataFrame, g: str) -> pd.DataFrame:
    return panel[panel.ensg_id == g]


def eligible(d: pd.DataFrame) -> bool:
    return (len(d) >= MIN_LINES
            and int(d.dep.sum()) >= MIN_POS
            and int((~d.dep).sum()) >= MIN_NEG)


def banner(title: str):
    print("=" * 78)
    print(title)
    print("=" * 78)


# =============================================================== gate readout
def gate_readout(panel: pd.DataFrame, score_col: str = "Z",
                 n_bins: int = 10, tie_mode: str = "ordinal",
                 region=GATE_REGION, npv_cuts=(5, 10, 20)) -> pd.DataFrame:
    """Per-gene gate statistics, one row per gene. This is the ONE definition of
    "the gate result" used by 02, 04 and 05, so those three cannot disagree
    about what they are measuring.

    Columns: base_rate, d1_ratio, d_last_ratio, upper_sd, pauc_gate,
             pauc_ranker, auc, npv{c}, npv_base, floor_rate, floor_ratio,
             n_floor, d1_hyper_p.
    """
    from scipy.stats import hypergeom

    rows = []
    for g, d in panel.groupby("ensg_id", sort=False):
        s = d[score_col].values.astype(float)
        p = d.dep.values.astype(bool)
        e = d.expr.values.astype(float)
        ok = ~np.isnan(s)
        if ok.sum() < MIN_LINES or p[ok].sum() < MIN_POS or (~p[ok]).sum() < MIN_NEG:
            continue
        s, p, e = s[ok], p[ok], e[ok]
        base = float(p.mean())
        rec = {"ensg_id": g, "n_lines": int(len(s)), "n_pos": int(p.sum()),
               "base_rate": base}

        q, keep = decile_index(s, n_bins, tie_mode)
        qv, pv = q[keep], p[keep]
        rates = np.array([pv[qv == i].mean() if (qv == i).sum() else np.nan
                          for i in range(n_bins)])
        ratios = rates / base
        rec["d1_ratio"] = float(ratios[0])
        rec["d_last_ratio"] = float(ratios[-1])
        rec["upper_sd"] = float(np.nanstd(ratios[1:]))
        rec["profile"] = [float(x) for x in ratios]   # list, so parquet accepts it
        n1 = int((qv == 0).sum())
        obs = int(pv[qv == 0].sum())
        rec["d1_hyper_p"] = float(hypergeom.cdf(obs, len(pv), int(pv.sum()), n1))

        rec["pauc_gate"] = partial_auc(s, p, *region)
        rec["pauc_ranker"] = partial_auc(s, p, *RANKER_REGION)
        rec["auc"] = auroc(s, p)

        order = np.argsort(s, kind="mergesort")
        for c in npv_cuts:
            k = max(1, int(round(len(s) * c / 100)))
            rec[f"npv{c}"] = float(1.0 - p[order[:k]].mean())
        rec["npv_base"] = float(1.0 - base)

        fl = e == np.nanmin(e)
        rec["n_floor"] = int(fl.sum())
        rec["floor_share"] = float(fl.mean())
        if fl.sum() >= 20:
            rec["floor_rate"] = float(p[fl].mean())
            rec["floor_ratio"] = rec["floor_rate"] / base
        rows.append(rec)
    return pd.DataFrame(rows)


def summarise_gate(R: pd.DataFrame, n_bins: int = 10) -> dict:
    """Collapse a gate_readout() frame to the headline numbers, with the
    significance tests the panel will ask about."""
    from scipy.stats import wilcoxon
    P = np.vstack(R.profile.values) if len(R) else np.zeros((0, n_bins))
    med_prof = [float(np.nanmedian(P[:, i])) for i in range(n_bins)]
    out = {"n_genes": int(len(R)), "median_profile": med_prof}
    if not len(R):
        return out
    d1 = R.d1_ratio.dropna()
    out["decile1_ratio"] = float(d1.median())
    out["decile1_wilcoxon_p"] = float(wilcoxon(d1 - 1.0).pvalue) if len(d1) > 5 else None
    out["frac_genes_depleted"] = float((d1 < 1).mean())
    out["upper_profile_sd"] = float(np.nanstd(med_prof[1:]))
    out["upper_profile_slope"] = float(
        np.polyfit(np.arange(2, n_bins + 1), med_prof[1:], 1)[0])
    pg = R.pauc_gate.dropna()
    pr = R.pauc_ranker.dropna()
    out["pauc_gate"] = float(pg.median())
    out["pauc_gate_p"] = float(wilcoxon(pg - 0.5).pvalue) if len(pg) > 5 else None
    out["pauc_gate_frac_above_chance"] = float((pg > 0.5).mean())
    out["pauc_ranker"] = float(pr.median())
    out["full_auc"] = float(R.auc.median())
    for c in (5, 10, 20):
        if f"npv{c}" in R:
            out[f"npv{c}_uplift"] = float((R[f"npv{c}"] - R.npv_base).median())
    F = R.dropna(subset=["floor_ratio"]) if "floor_ratio" in R else R.iloc[:0]
    if len(F) > 5:
        out["floor"] = {
            "n_genes": int(len(F)),
            "median_ratio": float(F.floor_ratio.median()),
            "median_rate": float(F.floor_rate.median()),
            "median_base": float(F.base_rate.median()),
            "frac_below_1": float((F.floor_ratio < 1).mean()),
            "wilcoxon_p": float(wilcoxon(F.floor_ratio - 1.0).pvalue),
        }
    return out
