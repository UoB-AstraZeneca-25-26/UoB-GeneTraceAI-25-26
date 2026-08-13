"""
expression_fusion/patch_core_score_orthogonal.py
-------------------------------------------------
Patches 02_core_score.ipynb to use the orthogonal RNA-protein combination
instead of Noisy-OR + correlation-penalised mean (Cells 8 and 8b).

WHAT THIS ADDS / CHANGES
-------------------------

Cell 7g (existing, UPDATED):
  - Removes rank(pct=True) — keeps rna_z_t as raw z-scores in rna_z_wide.
  - Still writes expr_pct (for backward-compat reference) but marks it unused.

Cell 7h (NEW, inserted after 7g):
  - Loads expression_fusion/outputs/bulk_prot_z.parquet -> prot_z_wide.
  - Same shape convention as rna_z_wide: model_id (UPPER) x ensg_id (lower).

Cell 7i (NEW, inserted after 7h):
  - Per-gene rho: Pearson correlation between rna_z and prot_z across the
    cell lines that have BOTH. Shrunk toward RHO_PRIOR = 0.373 when overlap
    < N_MIN_RHO = 30 (mirrors MAD_FLOOR in the RNA scorer).
  - Orthogonal combination:
      prot_resid_z = (prot_z - rho_g * rna_z) / sqrt(1 - rho_g^2)
      core_z = (W_RNA * rna_z + W_PROT * prot_resid_z) / norm
    where W_RNA = sqrt(2.0) [RNA Kish n_eff], W_PROT = sqrt(1.37) [prot Kish n_eff].
  - Single-source cells: core_z = rna_z_t (RNA only) or prot_z_t (protein only).
  - Final: core_score = scipy.stats.norm.cdf(core_z) maps z -> (0,1),
    preserving magnitude unlike rank(pct=True).
  - Sets n_layers (0/1/2) for stratum_rank in Cell 9.

Cell 8 + 8b: superseded. Left intact in the notebook (read-only history), but
their core_score output is overwritten by Cell 7i before Cell 9 runs.

WHY core_score = Phi(core_z) RATHER THAN rank(pct=True)
---------------------------------------------------------
Phi(z) = scipy.stats.norm.cdf(z) is monotone -> preserves ranking.
At the same time it preserves magnitudes: a cell at z=4 scores 0.99997,
a cell at z=1 scores 0.841, a cell at z=0 scores 0.500. These distances
are information -- percentile rank collapses them to 0.99 vs 0.50.

PRECONDITIONS
-------------
  Run 05_bulk_rna_scorer.py first  -> expression_fusion/outputs/bulk_rna_z.parquet
  Run 06_bulk_protein_scorer.py    -> expression_fusion/outputs/bulk_prot_z.parquet
  Then execute 02_core_score.ipynb from the top.
"""

import json
from pathlib import Path

NB  = Path(r"C:\Disertation\UoB-GeneTraceAI-25-26\src\pipeline\02_core_score.ipynb")
OUT_RNA  = r"C:\Disertation\UoB-GeneTraceAI-25-26\expression_fusion\outputs\bulk_rna_z.parquet"
OUT_PROT = r"C:\Disertation\UoB-GeneTraceAI-25-26\expression_fusion\outputs\bulk_prot_z.parquet"

nb = json.loads(NB.read_text(encoding="utf-8"))


# ──────────────────────────────────────────────────────────── helpers
def _make_code(source_lines: list[str]) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": source_lines[0].strip().lstrip("# -").strip()[:40].replace(" ", "-").lower(),
        "metadata": {},
        "outputs": [],
        "source": [l if l.endswith("\n") else l + "\n" for l in source_lines[:-1]]
               + [source_lines[-1].rstrip("\n")],
    }


def _make_md(source_lines: list[str]) -> dict:
    return {
        "cell_type": "markdown",
        "id": source_lines[0].strip().lstrip("# ").strip()[:40].replace(" ", "-").lower(),
        "metadata": {},
        "source": [l if l.endswith("\n") else l + "\n" for l in source_lines[:-1]]
               + [source_lines[-1].rstrip("\n")],
    }


# ──────────────────────────────── 1. Update Cell 7g: keep z-scores, drop rank()
RANK_LINE   = "    expr_pct = _rna_wide.rank(pct=True, method=\"min\")"
RANK_REPLACE = """\
    # Keep z-scores as rna_z_wide -- the orthogonal combination in Cell 7i
    # uses these directly. expr_pct kept for reference only (unused in scoring).
    rna_z_wide = _rna_wide.copy()
    expr_pct   = _rna_wide.rank(pct=True, method=\"min\")  # reference; not fed to Cell 8"""

cell7g_idx = None
for i, cell in enumerate(nb["cells"]):
    src = "".join(cell.get("source", []))
    if "# -- Cell 7g --" in src and cell["cell_type"] == "code":
        cell7g_idx = i
        if RANK_LINE in src:
            new_src = src.replace(RANK_LINE, RANK_REPLACE)
            cell["source"] = [l + "\n" for l in new_src.splitlines()]
            if cell["source"]:
                cell["source"][-1] = cell["source"][-1].rstrip("\n")
            cell["outputs"] = []
            cell["execution_count"] = None
            print(f"Cell {i}: updated Cell 7g to expose rna_z_wide (z-scores).")
        else:
            print(f"Cell {i}: Cell 7g found but rank() line not in expected form — skipping 7g edit.")
        break

if cell7g_idx is None:
    raise SystemExit("Could not find Cell 7g. Run patch_core_score.py first to add it.")


# ──────────────────────────────────── 2. Build Cell 7h — load prot z-scores
CELL_7H_MD = _make_md([
    "-- Cell 7h -- Load bulk protein z-scores (orthogonal combination input)\n",
    "\n",
    "Loads `bulk_prot_z.parquet` produced by\n",
    "`expression_fusion/06_bulk_protein_scorer.py`. Columns are lineage-\n",
    "conditioned Stouffer-combined z-scores from ProCAN (DIA-MS) + CCLE (TMT),\n",
    "normalised within each source and lineage before combining — the same\n",
    "approach as the RNA scorer.\n",
])

CELL_7H_CODE = _make_code([
    "# -- Cell 7h -- Protein z-scores from bulk scorer\n",
    "import pathlib as _pl\n",
    f'_prot_z_path = _pl.Path(r"{OUT_PROT}")\n',
    "\n",
    "if not _prot_z_path.exists():\n",
    "    print(\"bulk_prot_z.parquet not found -- protein layer will be skipped.\")\n",
    "    print(\"Run: cd expression_fusion && python 06_bulk_protein_scorer.py\")\n",
    "    prot_z_wide = None\n",
    "else:\n",
    "    _prot_z = pd.read_parquet(_prot_z_path)\n",
    "    _prot_z[\"ensg_id\"] = _prot_z[\"gene_id\"].str.lower()  # already lower; belt+braces\n",
    "    _prot_z_pivot = _prot_z.pivot(\n",
    "        index=\"model_id\", columns=\"ensg_id\", values=\"prot_z_t\"\n",
    "    )\n",
    "    _prot_z_pivot.index = _prot_z_pivot.index.str.upper()  # UPPER ACH- to match core score\n",
    "    prot_z_wide = _prot_z_pivot\n",
    "\n",
    "    # n_sources audit (1=ProCAN-only or CCLE-only, 2=both agreed)\n",
    "    _prot_nsrc = _prot_z.pivot(\n",
    "        index=\"model_id\", columns=\"ensg_id\", values=\"n_sources\"\n",
    "    )\n",
    "    _prot_nsrc.index = _prot_nsrc.index.str.upper()\n",
    "    prot_n_sources = _prot_nsrc\n",
    "\n",
    "    print(f\"prot_z_wide loaded: {prot_z_wide.shape}\")\n",
    "    print(f\"  model_ids: {prot_z_wide.shape[0]:,}\")\n",
    "    print(f\"  genes:     {prot_z_wide.shape[1]:,}\")\n",
    "    print(f\"  n_sources 1={int((_prot_z.n_sources==1).sum()):,}  \"\n",
    "          f\"2={int((_prot_z.n_sources==2).sum()):,}\")\n",
])


# ──────────────────────────────────── 3. Build Cell 7i — orthogonal combination
CELL_7I_MD = _make_md([
    "-- Cell 7i -- Orthogonal RNA-protein combination\n",
    "\n",
    "Replaces Noisy-OR (Cell 8) and correlation-penalised mean (Cell 8b).\n",
    "\n",
    "**Why orthogonalise instead of plain weighted sum?**\n",
    "RNA and protein are correlated (rho ~ 0.37-0.46). A weighted sum treats\n",
    "them as independent evidence, inflating the combined score. Partial\n",
    "residualisation removes the RNA-explained variance from the protein signal\n",
    "before adding it -- the remaining protein contribution is genuinely new\n",
    "information (post-transcriptional regulation, protein stability, etc.).\n",
    "\n",
    "**Combination formula** for cells with both RNA and protein:\n",
    "```\n",
    "prot_resid_z = (prot_z - rho_g * rna_z) / sqrt(1 - rho_g^2)\n",
    "core_z = (W_RNA * rna_z + W_PROT * prot_resid_z) / sqrt(W_RNA^2 + W_PROT^2)\n",
    "core_score = Phi(core_z)   # standard normal CDF -> (0, 1)\n",
    "```\n",
    "where rho_g is estimated per gene from overlapping lines (shrunk toward\n",
    "the global prior of 0.373 when fewer than 30 lines overlap).\n",
    "\n",
    "**Single-source cells:** `core_z = rna_z_t` (RNA only) or `prot_z_t`\n",
    "(protein only). `n_layers` (0/1/2) is set accordingly for `stratum_rank`.\n",
])

CELL_7I_CODE = _make_code([
    "# -- Cell 7i -- Orthogonal RNA-protein combination -> core_score\n",
    "from scipy import stats as _scipy_stats\n",
    "\n",
    "# Kish-corrected layer weights\n",
    "# RNA: W_RNA is PROVISIONAL. n_eff=2.0 requires rho_bar=0.25 across all 3\n",
    "# source pairs, but DepMap-HPA alone is 0.812, so rho_bar >= 0.25 is impossible.\n",
    "# True n_eff is estimated < 1.60 (see EXPR_FUSION_PREREGISTRATION.md T3).\n",
    "# Replace sqrt(2.00) with sqrt(n_eff_actual) once T3 has been run.\n",
    "# Protein: 2 sources, ProCAN-CCLE rho=0.373 -> n_eff = 2/(1+0.373) = 1.45\n",
    "_W_RNA  = np.sqrt(2.00)   # PROVISIONAL — awaiting T3 correction\n",
    "_W_PROT = np.sqrt(1.45)   # 1.204\n",
    "_W_NORM = np.sqrt(_W_RNA**2 + _W_PROT**2)\n",
    "\n",
    "_RHO_PRIOR = 0.373   # median ProCAN-CCLE Spearman rho from overlap test\n",
    "_N_MIN_RHO = 30      # shrink toward prior below this overlap count\n",
    "\n",
    "# ── align on common axes ──────────────────────────────────────────────────\n",
    "_rna = rna_z_wide  # model_id (UPPER) x ensg_id (lower)\n",
    "\n",
    "if prot_z_wide is None:\n",
    "    # Protein scorer not run: RNA-only score\n",
    "    print(\"WARNING: prot_z_wide is None. Scoring RNA layer only.\")\n",
    "    _all_models = sorted(_rna.index)\n",
    "    _all_genes  = sorted(_rna.columns)\n",
    "    core_z      = _rna.reindex(index=_all_models, columns=_all_genes)\n",
    "    n_layers    = (~core_z.isna()).astype(int)\n",
    "    rho_audit   = pd.Series(_RHO_PRIOR, index=_all_genes, name=\"rho_g\")\n",
    "else:\n",
    "    _prot = prot_z_wide\n",
    "\n",
    "    _all_models = sorted(set(_rna.index) | set(_prot.index))\n",
    "    _all_genes  = sorted(set(_rna.columns) | set(_prot.columns))\n",
    "\n",
    "    _R = _rna.reindex(index=_all_models, columns=_all_genes)   # NaN where absent\n",
    "    _P = _prot.reindex(index=_all_models, columns=_all_genes)\n",
    "\n",
    "    # ── per-gene rho estimation ───────────────────────────────────────────\n",
    "    # Only computable for genes in both layers; use prior for RNA-only genes.\n",
    "    _common_genes   = _rna.columns.intersection(_prot.columns)\n",
    "    _overlap_models = _rna.index.intersection(_prot.index)\n",
    "\n",
    "    _r_ov = _R.loc[_overlap_models, _common_genes].values.astype(float)  # (n_ov, n_cg)\n",
    "    _p_ov = _P.loc[_overlap_models, _common_genes].values.astype(float)\n",
    "\n",
    "    _both_valid = np.isfinite(_r_ov) & np.isfinite(_p_ov)\n",
    "    _n_valid    = _both_valid.sum(axis=0)  # per gene\n",
    "\n",
    "    # Pearson r per gene, vectorised\n",
    "    _r_m = np.where(_both_valid, _r_ov, np.nan)\n",
    "    _p_m = np.where(_both_valid, _p_ov, np.nan)\n",
    "    _r_c = _r_m - np.nanmean(_r_m, axis=0)\n",
    "    _p_c = _p_m - np.nanmean(_p_m, axis=0)\n",
    "    _num  = np.nansum(_r_c * _p_c,     axis=0)\n",
    "    _den  = np.sqrt(np.nansum(_r_c**2, axis=0) * np.nansum(_p_c**2, axis=0))\n",
    "    _rho_raw = np.where(_den > 0, _num / np.where(_den > 0, _den, 1.0), np.nan)\n",
    "\n",
    "    # James-Stein shrinkage toward prior\n",
    "    _lam     = np.where(_n_valid >= _N_MIN_RHO,\n",
    "                        _N_MIN_RHO / (_n_valid + _N_MIN_RHO), 1.0)\n",
    "    _rho_cg  = _lam * _RHO_PRIOR + (1 - _lam) * np.where(np.isfinite(_rho_raw),\n",
    "                                                           _rho_raw, _RHO_PRIOR)\n",
    "    _rho_cg  = np.clip(_rho_cg, -0.99, 0.99)\n",
    "\n",
    "    # Broadcast rho to all_genes (RNA-only genes get the prior)\n",
    "    rho_audit = pd.Series(_RHO_PRIOR, index=_all_genes, name=\"rho_g\")\n",
    "    rho_audit[_common_genes] = _rho_cg\n",
    "\n",
    "    print(f\"Per-gene rho: mean={rho_audit.mean():.3f}  \"\n",
    "          f\"median={rho_audit.median():.3f}  \"\n",
    "          f\"std={rho_audit.std():.3f}\")\n",
    "    print(f\"Genes with n_overlap >= {_N_MIN_RHO}: \"\n",
    "          f\"{int((_n_valid >= _N_MIN_RHO).sum()):,} / {len(_common_genes):,}\")\n",
    "\n",
    "    # ── orthogonal combination ────────────────────────────────────────────\n",
    "    # Align rho to gene axis of _R, _P\n",
    "    _rho_arr = rho_audit.reindex(_all_genes).values  # (n_genes,)\n",
    "\n",
    "    _R_vals = _R.values  # (n_models, n_genes)\n",
    "    _P_vals = _P.values\n",
    "\n",
    "    # Residual: protein signal orthogonal to RNA\n",
    "    _prot_resid    = _P_vals - _rho_arr[np.newaxis, :] * _R_vals\n",
    "    _resid_scale   = np.sqrt(1.0 - _rho_arr**2)[np.newaxis, :]  # (1, n_genes)\n",
    "    _resid_scale   = np.where(_resid_scale > 0, _resid_scale, 1.0)\n",
    "    _prot_resid_z  = _prot_resid / _resid_scale\n",
    "\n",
    "    # Weighted combination where BOTH layers present\n",
    "    _rna_present  = np.isfinite(_R_vals)\n",
    "    _prot_present = np.isfinite(_P_vals)\n",
    "    _both_present = _rna_present & _prot_present\n",
    "\n",
    "    _core_z = np.full_like(_R_vals, np.nan)\n",
    "\n",
    "    # Both: orthogonal weighted combination\n",
    "    _core_z[_both_present] = (\n",
    "        _W_RNA  * _R_vals[_both_present]\n",
    "        + _W_PROT * _prot_resid_z[_both_present]\n",
    "    ) / _W_NORM\n",
    "\n",
    "    # RNA only\n",
    "    _rna_only = _rna_present & ~_prot_present\n",
    "    _core_z[_rna_only] = _R_vals[_rna_only]\n",
    "\n",
    "    # Protein only (rare: a model_id in protein but not RNA)\n",
    "    _prot_only = ~_rna_present & _prot_present\n",
    "    _core_z[_prot_only] = _P_vals[_prot_only]\n",
    "\n",
    "    core_z   = pd.DataFrame(_core_z, index=_all_models, columns=_all_genes)\n",
    "    n_layers = _rna_present.astype(int) + _prot_present.astype(int)\n",
    "    n_layers = pd.DataFrame(n_layers, index=_all_models, columns=_all_genes)\n",
    "\n",
    "# ── map z-scores to (0, 1) via standard normal CDF ───────────────────────\n",
    "# Phi(z) preserves all magnitude information: z=4 -> 0.99997, z=1 -> 0.841.\n",
    "# Unlike rank(pct=True), two lines at z=4 vs z=1 are NOT equidistant.\n",
    "core_score = core_z.apply(\n",
    "    lambda col: _scipy_stats.norm.cdf(col.values), axis=0\n",
    ")\n",
    "core_score = pd.DataFrame(\n",
    "    core_score, index=core_z.index, columns=core_z.columns\n",
    ")\n",
    "\n",
    "# Propagate NaN (no evidence) through Phi\n",
    "core_score[core_z.isna()] = np.nan\n",
    "\n",
    "_n_scored = core_score.notna().sum().sum()\n",
    "_n_total  = core_score.size\n",
    "print(f\"\\ncore_score shape  : {core_score.shape}\")\n",
    "print(f\"non-NaN cells     : {_n_scored:,} / {_n_total:,} \"\n",
    "      f\"({_n_scored/_n_total*100:.1f}%)\")\n",
    "print(f\"score range       : [{core_score.stack().min():.4f}, \"\n",
    "      f\"{core_score.stack().max():.4f}]\")\n",
    "print(f\"n_layers distribution:\")\n",
    "print(n_layers.stack().value_counts().sort_index().to_string())\n",
])


# ──────────────────────────────────────── 4. Find insertion point and patch
# Insert 7h and 7i after Cell 7g (cell7g_idx).
# Cells 8 and 8b remain in the notebook but their core_score output is
# overwritten by Cell 7i before Cell 9 runs.

insert_at = cell7g_idx + 1
nb["cells"].insert(insert_at,     CELL_7H_CODE)
nb["cells"].insert(insert_at,     CELL_7H_MD)
nb["cells"].insert(insert_at + 2, CELL_7I_CODE)
nb["cells"].insert(insert_at + 2, CELL_7I_MD)

print(f"Inserted Cells 7h (md+code) and 7i (md+code) after Cell 7g.")
print(f"Notebook now has {len(nb['cells'])} cells.")
print()
print("WHAT TO DO NEXT:")
print("  1. python expression_fusion/05_bulk_rna_scorer.py   (if not already done)")
print("  2. python expression_fusion/06_bulk_protein_scorer.py")
print("  3. Run 02_core_score.ipynb top-to-bottom in Jupyter")
print()
print("CELLS 8 AND 8b remain in the notebook as historical record.")
print("Their core_score output is overwritten by Cell 7i before Cell 9 runs.")

NB.write_text(json.dumps(nb, indent=1), encoding="utf-8")
print(f"\nPatched notebook written to {NB}")
