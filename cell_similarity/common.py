"""
cell_similarity/common.py
-------------------------
Single source of truth for the cell-line similarity graph. Loaders, the
similarity metric, and the neighbour machinery live here so the three build /
validate / query scripts cannot drift apart.

Read-only with respect to the pipeline. Writes only into
`cell_similarity/outputs/`.

THREE AXES, KEPT SEPARATE ON PURPOSE
------------------------------------
    rna     DepMap 24Q4 OmicsExpressionProteinCodingGenesTPMLogp1  (1,673 lines)
    metab   CCLE metabolomics, 225 metabolites                     (  928 lines)
    mirna   CCLE miRNA, model-level                                (  ~950 lines)

They are NOT fused. Each produces its own line x line similarity matrix and its
own neighbour list. Fusing before knowing what each contributes would make it
impossible to say afterwards which axis carried a result — and the axes have
very different coverage, so a fused score would silently be "RNA wherever the
other two are missing".

THE METRIC, AND WHY
-------------------
For each axis, in order:

  1. Drop features missing in more than `MAX_FEATURE_MISSING` of lines.
  2. Median-impute what remains (never mean — metabolite and miRNA
     distributions are skewed).
  3. Keep the `TOP_VAR` most variable features. Similarity between cell lines is
     carried by the features that differ between them; invariant features add a
     constant to every pair and compress the dynamic range of the whole matrix.
  4. Standardise each feature to z (so a high-variance feature does not dominate
     purely through scale).
  5. Pearson correlation between LINES over those features.

Pearson on standardised features is cosine similarity on centred data, which is
the conventional choice for expression-profile similarity and is what the
lineage-clustering literature uses. Spearman is available via `method="spearman"`
and is reported alongside in 01 as a robustness check, because a rank metric is
less sensitive to the imputation in step 2.

A NOTE THAT DECIDES HOW EVERYTHING HERE IS READ
------------------------------------------------
Cell-line similarity is dominated by lineage. That is not a flaw, it is the
strongest real signal in the data — but it means "close pairs behave similarly"
is nearly guaranteed to pass as stated, because close pairs are mostly
same-tissue pairs and same-tissue lines share dependencies for reasons that have
nothing to do with this graph. Every validation in `02` is therefore run twice:
pooled, and WITHIN LINEAGE. Only the within-lineage version tells you whether the
graph adds anything a tissue label does not already give you for free.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)

# Deliberate reuse: the DepMap 24Q4 / DEMETER2 loaders and the scale checks live
# in gate_audit/common.py. Importing them is better than copying them, because a
# copy would drift and the scale checks are the guard against a mis-scaled label.
# Loaded by explicit path under a distinct module name -- both folders have a
# `common.py`, and a plain `import common` resolves to whichever directory is
# first on sys.path (i.e. this one).
import importlib.util as _ilu  # noqa: E402

_ga_spec = _ilu.spec_from_file_location("gate_audit_common",
                                        ROOT / "gate_audit" / "common.py")
GA = _ilu.module_from_spec(_ga_spec)
sys.modules["gate_audit_common"] = GA
_ga_spec.loader.exec_module(GA)

CLEANED = ROOT / "cleaned_track_data"
DATA_CLEAN = ROOT / "data" / "parquet" / "data_clean"

MAX_FEATURE_MISSING = 0.20
TOP_VAR = 2000
K_DEFAULT = 25
SEED = 42

AXES = ("rna", "metab", "mirna")
AXIS_LABEL = {"rna": "RNA expression", "metab": "metabolomics", "mirna": "miRNA"}


def banner(title: str):
    print("=" * 78)
    print(title)
    print("=" * 78)


# =============================================================== loaders
# DepMap 24Q4 OncotreeLineage -> the 20Q2-era sample_info vocabulary. Explicit
# because the two sources name the same tissues differently; anything not listed
# is snake-cased as-is. Curation step, recorded rather than hidden.
_ONCOTREE_MAP = {
    "Bowel": "colorectal", "CNS/Brain": "central_nervous_system",
    "Lymphoid": "lymphocyte", "Myeloid": "blood",
    "Head and Neck": "upper_aerodigestive", "Esophagus/Stomach": "gastric",
    "Biliary Tract": "bile_duct", "Ampulla of Vater": "bile_duct",
    "Soft Tissue": "soft_tissue", "Peripheral Nervous System":
        "peripheral_nervous_system", "Bladder/Urinary Tract": "urinary_tract",
    "Uterus": "uterus", "Ovary/Fallopian Tube": "ovary",
    "Pleura": "pleura", "Eye": "eye", "Skin": "skin", "Lung": "lung",
    "Breast": "breast", "Kidney": "kidney", "Pancreas": "pancreas",
    "Liver": "liver", "Bone": "bone", "Thyroid": "thyroid",
    "Prostate": "prostate", "Cervix": "cervix", "Testis": "testis",
    "Adrenal Gland": "adrenal_cortex", "Vulva/Vagina": "vulva",
    "Fibroblast": "fibroblast", "Muscle": "muscle", "Embryonal": "embryo",
    "Other": "unknown",
}


def load_lineage() -> pd.DataFrame:
    """model_id -> lineage / lineage_subtype / lineage_source.

    Primary source is the pipeline's own cleaned `sample_info.parquet` (30-value
    controlled vocabulary). That table is 20Q2-vintage and covers 1,840 lines;
    DepMap 24Q4 expression adds **246 lines it has never heard of**, which
    previously fell through as UNLABELLED — 14.7% of the panel, and the worst
    cluster in every query.

    Those are filled from `data/DepMap_24Q4/Model.csv` (`OncotreeLineage`),
    mapped onto the sample_info vocabulary by `_ONCOTREE_MAP`. `lineage_source`
    records which table each label came from, because the two vocabularies are
    harmonised, not identical.
    """
    si = pd.read_parquet(CLEANED / "sample_info.parquet",
                         columns=["model_id", "lineage", "lineage_subtype"])
    si["model_id"] = si.model_id.astype(str).str.lower()
    si = si.dropna(subset=["model_id"]).drop_duplicates("model_id")
    si = si.set_index("model_id")
    si["lineage_source"] = "sample_info_20Q2"

    mf = GA.DEPMAP_24Q4 / "Model.csv"
    if mf.exists():
        md = pd.read_csv(mf, low_memory=False,
                         usecols=["ModelID", "OncotreeLineage",
                                  "OncotreeSubtype"])
        md["model_id"] = md.ModelID.astype(str).str.lower()
        md = md.drop_duplicates("model_id").set_index("model_id")
        md["lineage"] = md.OncotreeLineage.map(
            lambda v: _ONCOTREE_MAP.get(v, str(v).lower().replace(" ", "_")
                                        .replace("/", "_"))
            if pd.notna(v) else None)
        add = md.index.difference(si.index)
        extra = pd.DataFrame({
            "lineage": md.loc[add, "lineage"],
            "lineage_subtype": md.loc[add, "OncotreeSubtype"],
            "lineage_source": "model_csv_24Q4",
        })
        si = pd.concat([si, extra])
    return si


def load_rna() -> pd.DataFrame:
    """DepMap 24Q4 expression, lines x ensg, log2(TPM+1). Keyed by ModelID."""
    _, _, _, sym2ensg, _ = GA.load_gene_lookup()
    return GA.load_expression_24q4(sym2ensg, verbose=False)


def load_metab() -> pd.DataFrame:
    """CCLE metabolomics, lines x metabolite. 928 lines, 225 metabolites."""
    d = pd.read_parquet(CLEANED / "metabolomics.parquet")
    d = d[d.model_id.notna()].copy()
    d["model_id"] = d.model_id.astype(str).str.lower()
    feats = [c for c in d.columns if c not in ("ccle_id", "model_id")]
    X = d.set_index("model_id")[feats].astype(float)
    return X.groupby(level=0).mean()


def load_mirna() -> pd.DataFrame:
    """CCLE miRNA, long -> lines x miRNA."""
    d = pd.read_parquet(CLEANED / "mirna_model_level.parquet")
    d = d[d.model_id.notna()].copy()
    d["model_id"] = d.model_id.astype(str).str.lower()
    X = d.pivot_table(index="model_id", columns="mirna_id",
                      values="expression_value", aggfunc="mean")
    return X.astype(float)


LOADERS = {"rna": load_rna, "metab": load_metab, "mirna": load_mirna}


# =============================================================== similarity
def prepare(X: pd.DataFrame, top_var: int = TOP_VAR,
            max_missing: float = MAX_FEATURE_MISSING, verbose: bool = True):
    """Drop sparse features, median-impute, keep the most variable, z-score."""
    n0 = X.shape[1]
    keep = X.isna().mean(axis=0) <= max_missing
    X = X.loc[:, keep]
    X = X.fillna(X.median(axis=0))
    X = X.loc[:, X.std(axis=0) > 0]
    v = X.var(axis=0).sort_values(ascending=False)
    X = X[v.index[:min(top_var, len(v))]]
    Z = (X - X.mean(axis=0)) / X.std(axis=0)
    if verbose:
        print(f"    features {n0:,} -> {X.shape[1]:,} kept "
              f"(dropped {int((~keep).sum()):,} with >{max_missing:.0%} missing)")
    return Z


def similarity(X: pd.DataFrame, method: str = "pearson", **kw) -> pd.DataFrame:
    """line x line similarity. Pearson on z-scored features = cosine on centred
    data. Spearman ranks features within each line first, which is robust to the
    median imputation but slower."""
    Z = prepare(X, **kw)
    if method == "spearman":
        Z = Z.rank(axis=1)
        Z = (Z.sub(Z.mean(axis=1), axis=0)).div(Z.std(axis=1), axis=0)
        M = Z.values
        S = (M @ M.T) / M.shape[1]
    else:
        M = Z.values
        M = M - M.mean(axis=1, keepdims=True)
        norm = np.linalg.norm(M, axis=1, keepdims=True)
        norm[norm == 0] = 1.0
        M = M / norm
        S = M @ M.T
    np.fill_diagonal(S, np.nan)          # a line is not its own neighbour
    return pd.DataFrame(S, index=X.index, columns=X.index)


def neighbours(S: pd.DataFrame, k: int = K_DEFAULT) -> pd.DataFrame:
    """Long-form k nearest neighbours: model_id, neighbour, rank, similarity."""
    rows = []
    idx = S.index.to_numpy()
    V = S.values
    kk = min(k, V.shape[1] - 1)
    for i in range(V.shape[0]):
        v = V[i]
        ok = ~np.isnan(v)
        if ok.sum() < 1:
            continue
        order = np.argsort(-np.where(ok, v, -np.inf))[:kk]
        for r, j in enumerate(order, start=1):
            rows.append((idx[i], idx[j], r, float(v[j])))
    return pd.DataFrame(rows, columns=["model_id", "neighbour", "rank",
                                       "similarity"])


def upper_pairs(S: pd.DataFrame, max_pairs: int | None = None, seed: int = SEED):
    """Unique unordered pairs as (i_idx, j_idx, sim), optionally subsampled.
    Pair counts are quadratic — 1,673 lines is 1.4M pairs — so most tests
    subsample rather than materialise everything."""
    n = S.shape[0]
    iu, ju = np.triu_indices(n, k=1)
    v = S.values[iu, ju]
    ok = ~np.isnan(v)
    iu, ju, v = iu[ok], ju[ok], v[ok]
    if max_pairs and len(v) > max_pairs:
        rng = np.random.default_rng(seed)
        sel = rng.choice(len(v), max_pairs, replace=False)
        iu, ju, v = iu[sel], ju[sel], v[sel]
    return iu, ju, v


# =============================================================== artefacts
def axis_path(axis: str, what: str) -> Path:
    return OUT / f"{axis}_{what}.parquet"


def save_axis(axis: str, S: pd.DataFrame, NN: pd.DataFrame):
    S.to_parquet(axis_path(axis, "similarity"))
    NN.to_parquet(axis_path(axis, "neighbours"), index=False)


def load_axis_similarity(axis: str) -> pd.DataFrame:
    p = axis_path(axis, "similarity")
    if not p.exists():
        raise SystemExit(f"{p} not found. Run cell_similarity/01_build_axes.py first.")
    return pd.read_parquet(p)


def load_axis_neighbours(axis: str) -> pd.DataFrame:
    p = axis_path(axis, "neighbours")
    if not p.exists():
        raise SystemExit(f"{p} not found. Run cell_similarity/01_build_axes.py first.")
    return pd.read_parquet(p)


def available_axes() -> list:
    return [a for a in AXES if axis_path(a, "similarity").exists()]


# =============================================================== dependency
def common_essentials(sym2ensg: dict, which: str = "inferred") -> set:
    """DepMap's own common-essential list, as ensg ids.

    `inferred`  CRISPRInferredCommonEssentials.csv (1,523 genes) -- DepMap's
                data-driven call, the right list for "remove the housekeeping
                floor".
    `achilles`  AchillesCommonEssentialControls.csv (1,247 genes) -- the curated
                control set.
    """
    f = ("CRISPRInferredCommonEssentials.csv" if which == "inferred"
         else "AchillesCommonEssentialControls.csv")
    p = GA.DEPMAP_24Q4 / f
    if not p.exists():
        raise SystemExit(f"{p} not found.")
    syms = (pd.read_csv(p).iloc[:, 0].astype(str)
            .str.split(" (", regex=False).str[0].str.upper())
    return {sym2ensg[s] for s in syms if s in sym2ensg}


def dependency_correlation(W: pd.DataFrame) -> pd.DataFrame:
    """line x line correlation of dependency profiles, diagonal NaN."""
    Wz = W.sub(W.mean(axis=1), axis=0)
    M = np.nan_to_num(Wz.values, nan=0.0)
    norm = np.linalg.norm(M, axis=1, keepdims=True)
    norm[norm == 0] = 1.0
    M = M / norm
    D = pd.DataFrame(M @ M.T, index=W.index, columns=W.index)
    np.fill_diagonal(D.values, np.nan)
    return D


def load_rnai(verbose: bool = True) -> pd.DataFrame:
    """DEMETER2 RNAi gene dependency, lines x ensg. A different perturbation
    from CRISPR, and — critically for validating an RNA-derived graph — one that
    does not share CRISPR's expression confound."""
    _, _, _, sym2ensg, _ = GA.load_gene_lookup()
    W, chk = GA.load_demeter2(sym2ensg, verbose=verbose)
    return W


def load_dependency(min_lines: int = 100, selective_only: bool = True,
                    verbose: bool = True) -> pd.DataFrame:
    """CRISPR 24Q4 gene effect, lines x ensg, for the week-2 concordance test.

    `selective_only` keeps genes whose effect actually VARIES across lines.
    Pan-essential genes (ribosome, proteasome, spliceosome) are near -1
    everywhere; including them makes every pair of cell lines correlate highly
    for a reason that has nothing to do with similarity, and would let the test
    pass no matter what the graph looked like.
    """
    _, _, _, sym2ensg, _ = GA.load_gene_lookup()
    W, chk = GA.load_crispr_24q4(sym2ensg, verbose=verbose)
    if selective_only:
        sd = W.std(axis=0)
        rng_ = W.max(axis=0) - W.min(axis=0)
        keep = (sd > 0.20) & (W.mean(axis=0) > -0.75) & rng_.notna()
        if verbose:
            print(f"    selective genes: {int(keep.sum()):,} of {W.shape[1]:,} "
                  f"(sd>0.20 and not pan-essential)")
        W = W.loc[:, keep]
    return W
