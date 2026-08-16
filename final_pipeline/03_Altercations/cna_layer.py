"""
03_Altercations/cna_layer.py
------------------------------
CNA confidence modifier flags from the harmonised COSMIC CNA warehouse table.

CNA is a CONFIDENCE MODIFIER, not a ranking layer:
  oncogene / activation_driven  -> amplification raises confidence
  tsg / loss_of_function        -> deletion raises confidence
  unknown_dual                  -> CNA surfaced as context only

Ploidy normalisation:
  The warehouse table main.cosmic_cna contains COSMIC's own cna_call column
  ('amplification' / 'deletion'), which is produced by COSMIC's internal
  ploidy-aware calling pipeline. Using cna_call rather than re-thresholding
  TOTAL_CN avoids the absolute-threshold ploidy confound: a tetraploid line's
  baseline CN=4 is correctly treated as normal by COSMIC, whereas a fixed
  threshold of 2.5 would call it amplified.

  Diagnostic (2026-08): median per-line TOTAL_CN = 3.0, with 33% of lines
  having mean CN > 4.5 (near-tetraploid). Using COSMIC's pre-called cna_call
  makes amplification/deletion assignments ploidy-relative by construction.

Reads from:  celllineselector.db  (main.cosmic_cna — already harmonised)
             reference/gene_lookup.parquet
Writes to:   final_pipeline/outputs/cna_flags.parquet
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import DB, GENE_LKP, CNA_FLAGS
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
import common as C


def run():
    print("=" * 70)
    print("03_Altercations — CNA layer (COSMIC ploidy-aware calls)")
    print("=" * 70)

    con = C.connect()
    cna = con.execute(
        """
        SELECT model_id, lower(gene_id) AS ensg_id, total_cn, cna_call
        FROM main.cosmic_cna
        WHERE model_id IS NOT NULL
          AND gene_id  IS NOT NULL
          AND cna_call IS NOT NULL
          AND is_ambiguous = FALSE
        """
    ).df()
    con.close()

    print(f"cosmic_cna rows (unambiguous): {len(cna):,}")
    print(f"  lines: {cna.model_id.nunique():,}  |  genes: {cna.ensg_id.nunique():,}")
    print(f"  call distribution:\n{cna.cna_call.value_counts().to_string()}")

    # Gene role for gating — lowercase to match warehouse key convention
    gl = pd.read_parquet(GENE_LKP, columns=["ensg_id", "gene_role"])
    gl["ensg_id"] = gl["ensg_id"].str.lower()
    cna = cna.merge(gl, on="ensg_id", how="left")
    cna["gene_role"] = cna["gene_role"].fillna("unknown")

    # Use COSMIC's calls directly — already ploidy-aware
    cna["is_amplification"] = cna["cna_call"] == "amplification"
    cna["is_deletion"]      = cna["cna_call"] == "deletion"

    # Role-gated flags
    cna["has_cna_alteration"] = (
        (cna["is_amplification"] & cna["gene_role"].isin(["oncogene", "both"])) |
        (cna["is_deletion"]      & cna["gene_role"].isin(["tsg", "both"]))
    )
    cna["has_cna_context"] = (
        (cna["is_amplification"] | cna["is_deletion"]) & ~cna["has_cna_alteration"]
    )

    flags = (
        cna.groupby(["model_id", "ensg_id"])
        .agg(
            has_cna_alteration=("has_cna_alteration", "any"),
            has_cna_context   =("has_cna_context",    "any"),
            total_cn_mean     =("total_cn",            "mean"),
        )
        .reset_index()
    )
    flags.to_parquet(CNA_FLAGS, index=False)

    print(f"\nFlags: {len(flags):,} (gene, line) pairs  ->  {CNA_FLAGS}")
    print(f"  has_cna_alteration: {flags.has_cna_alteration.sum():,}")
    print(f"  has_cna_context:    {flags.has_cna_context.sum():,}")


if __name__ == "__main__":
    run()
