"""
Shared style, paths and loaders for the report figure suite.

Every figure in outputs/report_figures/ is generated from this package.
Numbers are never hard-coded in the plotting code -- they come from
_stats.json, which compute_stats.py writes straight out of the production
parquet outputs. That way a figure can never silently drift from the data.
"""
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(r"C:\Disertation\UoB-GeneTraceAI-25-26")
PIPE = ROOT / "src" / "pipeline" / "outputs"
REF = ROOT / "reference"
DATA_CLEAN = ROOT / "data" / "parquet" / "data_clean"
TRACK = ROOT / "cleaned_track_data"
VALID = ROOT / "validation" / "prepared"
OUT = ROOT / "outputs" / "report_figures"
STATS = OUT / "_stats.json"

OUT.mkdir(parents=True, exist_ok=True)

# ── Palette ───────────────────────────────────────────────────────────────
# Sequential-safe, colour-blind-tolerant, prints legibly in greyscale.
BLUE = "#2a6fb5"      # primary / RNA / production arm
RED = "#c8443c"       # rejected / negative / protein
AMBER = "#e0a020"     # caution / modifier
GREEN = "#3e8f60"     # validated / positive result
GREY = "#7a7a7a"
LIGHT = "#dfe6ee"
INK = "#22262b"

STAGE_FILL = "#e8f0fa"
STAGE_EDGE = BLUE
SOURCE_FILL = "#f3ede2"
SOURCE_EDGE = "#b08a4a"
DEAD_FILL = "#f6e3e2"
DEAD_EDGE = RED

plt.rcParams.update({
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "#9aa2ab",
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#e6e9ed",
    "grid.linewidth": 0.7,
    "font.size": 9.5,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "legend.frameon": False,
    "figure.dpi": 110,
})

DPI = 200


def save(fig, name):
    """Write a figure to outputs/report_figures/<name>.png at print DPI."""
    path = OUT / f"{name}.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    kb = path.stat().st_size // 1024
    print(f"  saved {path.name}  ({kb} KB)")
    return path


def load_stats():
    with open(STATS, encoding="utf-8") as fh:
        return json.load(fh)


def caption(fig, text, y=-0.02):
    """One-line provenance note under a figure -- what it was measured from."""
    fig.text(0.5, y, text, ha="center", va="top", fontsize=7.5, color=GREY,
             style="italic", wrap=True)


# ── Loaders (the same joins the pipeline itself uses) ─────────────────────
def load_universe():
    """Protein-coding + HGNC-Approved gene universe, matching 02_core_score."""
    gl = pd.read_parquet(REF / "gene_lookup.parquet")
    u = gl[(gl["biotype"] == "protein_coding") &
           (gl["hgnc_status"] == "Approved")].copy()
    u["ensg_id"] = u["ensg_id"].astype("string").str.split(".").str[0].str.lower()
    return u


def load_expr_wide(genes=None):
    """(model_id x ensg_id) log2 TPM+1 matrix, profile IDs resolved to ACH-,
    duplicate profiles averaged -- reproduces 02_core_score cells 4 and 7."""
    prof = pd.read_parquet(REF / "depmap_profiles.parquet")
    rna = prof[prof["datatype"] == "rna"][["profileid", "modelid"]].copy()
    rna["model_id"] = rna["modelid"].str.lower()

    expr = pd.read_parquet(DATA_CLEAN / "depmap_expr_clean.parquet")
    expr.columns = [c.split(".")[0].lower() for c in expr.columns]
    if genes is not None:
        keep = [c for c in expr.columns if c in set(genes)]
        expr = expr[keep]
    expr.index.name = "profileid"
    expr = expr.reset_index().merge(
        rna[["profileid", "model_id"]], on="profileid", how="inner")
    expr = expr.drop(columns=["profileid"]).groupby("model_id").mean()
    return expr


def load_prot_wide(genes=None):
    """(model_id x ensg_id) proteomics log-ratio matrix, UniProt -> ENSG
    mapped -- reproduces 02_core_score cell 5."""
    u = load_universe()
    up2ensg = (u[["ensg_id", "uniprot_ids"]].dropna(subset=["uniprot_ids"])
               .assign(uniprot_id=lambda x: x["uniprot_ids"].str.strip().str.lower())
               .drop_duplicates(subset=["uniprot_id"])[["uniprot_id", "ensg_id"]])

    raw = pd.read_parquet(TRACK / "proteomics.parquet")
    key = "model_id" if "model_id" in raw.columns else "depmap_id"
    long = raw.melt(id_vars=[key], var_name="uniprot_id", value_name="log_ratio")
    long = long.dropna(subset=["log_ratio"])
    long["model_id"] = long[key].str.lower()
    long = long.merge(up2ensg, on="uniprot_id", how="inner")
    if genes is not None:
        long = long[long["ensg_id"].isin(set(genes))]
    wide = (long.groupby(["model_id", "ensg_id"])["log_ratio"].mean()
                .unstack("ensg_id"))
    return wide


def load_gdsc_truth():
    """gene -> set of GDSC-sensitive model_ids, keyed on ensg_id."""
    pairs = pd.read_parquet(VALID / "gdsc_pairs.parquet")
    gdsc = pd.read_parquet(VALID / "gdsc_ach.parquet")
    gt = (pairs[["drug_id", "target_ensg"]]
          .merge(gdsc[gdsc["sensitive"] == True][["drug_id", "model_id"]],
                 on="drug_id")
          .rename(columns={"target_ensg": "ensg_id"}))
    # BOTH keys, not just model_id. The returned dict is keyed on ensg_id and is
    # looked up with ids that come from core_score / predictions, which are
    # lowercase -- an uppercase key here matches nothing and the failure is
    # silent: empty frames, then a KeyError several functions downstream.
    gt["model_id"] = gt["model_id"].str.lower()
    gt["ensg_id"] = gt["ensg_id"].astype("string").str.split(".").str[0].str.lower()
    return gt.groupby("ensg_id")["model_id"].apply(lambda s: set(s)).to_dict()


def hit_at_k(ranked_model_ids, known, k=20):
    """hit@k with the consistent denominator: |known n scored|, not |known|.
    ranked_model_ids must already be sorted best-first."""
    scored = known & set(ranked_model_ids)
    if not scored:
        return np.nan
    return len(set(ranked_model_ids[:k]) & known) / len(scored)
