"""
Scoring/core_score_variant_A_current.py  (TEST HARNESS, documentation-only --
not a computation. Not wired into run_all.py.)
-----------------------------------------------------------------------------
Variant A of the three-way residualization comparison IS the live, already-
promoted core_score.py formula:

    prot_resid   = prot_z - beta_g * rna_z
    prot_resid_z = standardize_with_shrinkage_and_clip(prot_resid)   [per gene,
                   w=n/(n+25) toward panel-wide SD, clip +/-5]
    rna_std_z    = standardize_by_stratum(rna_z)                     [per
                   (gene, n_sources), w=n/(n+25) toward gene-global SD when
                   n<5, clip +/-5]
    core_z       = (W_rna*rna_std_z + W_prot*prot_resid_z) / norm

No separate computation is needed or run here -- its output is simply the
live final_pipeline/outputs/core_score.parquet, already verified structurally
(row/gene/line counts, n_layers=1 bit-identity, spike figures) in the prior
promotion task. This file exists purely so the three-way comparison has a
named, documented counterpart to _variant_B_flipped.py and
_variant_C_no_residualization.py, and so "Variant A" unambiguously maps to a
specific file/commit state (core_score.py as promoted, PROT_RESID_N0=25,
RNA_STRATUM_N0=25, both WINSOR bounds=5.0) rather than a moving target.

Reads from (for reference, not re-run): final_pipeline/outputs/core_score.parquet
"""
from pathlib import Path

CORE_SCORE_LIVE = Path(__file__).resolve().parent.parent / "outputs" / "core_score.parquet"

if __name__ == "__main__":
    print(__doc__)
    print(f"Variant A == {CORE_SCORE_LIVE}  (exists: {CORE_SCORE_LIVE.exists()})")
