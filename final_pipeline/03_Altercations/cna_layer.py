"""
03_Altercations/cna_layer.py
------------------------------
CNA confidence modifier flags from COSMIC CellLines CNA data.

CNA is a CONFIDENCE MODIFIER, not a ranking layer:
  oncogene / activation_driven  -> amplification (CN > 2.5) raises confidence
  tsg / loss_of_function        -> deletion (CN < 1.5) raises confidence
  unknown_dual                  -> CNA surfaced as context only

Reads from:  data/COSMIC/CellLinesProject_CompleteCNA_v104_GRCh37.tsv
             data/GDSC/model_list_20260709.csv
             reference/gene_lookup.parquet
Writes to:   final_pipeline/outputs/cna_flags.parquet
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import REPO, GENE_LKP, CNA_FLAGS, COSMIC_CNA, GDSC_MODEL

AMP_THRESHOLD = 2.5
DEL_THRESHOLD = 1.5


def run():
    print("=" * 70)
    print("03_Altercations — CNA layer")
    print("=" * 70)

    print("Loading COSMIC CNA...")
    cna = pd.read_csv(COSMIC_CNA, sep="\t",
                      usecols=["GENE_SYMBOL", "SAMPLE_NAME", "TOTAL_CN", "MUT_TYPE"])
    print(f"  Raw rows: {len(cna):,}  |  genes: {cna.GENE_SYMBOL.nunique():,}")

    model = pd.read_csv(GDSC_MODEL, usecols=["model_name", "BROAD_ID"]).dropna(subset=["BROAD_ID"])
    model["model_id"]    = model["BROAD_ID"].str.lower()
    model["sample_name"] = model["model_name"].str.strip()
    model = model.drop_duplicates(subset="sample_name", keep="first")

    cna = cna.merge(model[["sample_name", "model_id"]],
                    left_on="SAMPLE_NAME", right_on="sample_name", how="left")
    cna = cna[cna["model_id"].notna()].copy()
    print(f"  After bridge: {len(cna):,} rows  |  lines: {cna.model_id.nunique():,}")

    gl = pd.read_parquet(GENE_LKP, columns=["ensg_id", "hgnc_symbol", "gene_role"])
    gl["hgnc_symbol"] = gl["hgnc_symbol"].astype("string").str.lower()  # for COSMIC symbol join only
    cna["GENE_SYMBOL"] = cna["GENE_SYMBOL"].astype("string").str.strip().str.lower()

    cna = cna.merge(gl, left_on="GENE_SYMBOL", right_on="hgnc_symbol", how="left")
    cna = cna[cna["ensg_id"].notna()].copy()
    cna["gene_role"] = cna["gene_role"].fillna("unknown")
    print(f"  After ENSG join: {len(cna):,} rows  |  genes: {cna.ensg_id.nunique():,}")

    cna["is_amplification"] = cna["TOTAL_CN"] > AMP_THRESHOLD
    cna["is_deletion"]      = cna["TOTAL_CN"] < DEL_THRESHOLD

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
            total_cn_mean     =("TOTAL_CN",            "mean"),
        )
        .reset_index()
    )
    flags.to_parquet(CNA_FLAGS, index=False)
    print(f"\nFlags: {len(flags):,} (gene, line) pairs  ->  {CNA_FLAGS}")
    print(f"  has_cna_alteration: {flags.has_cna_alteration.sum():,}")
    print(f"  has_cna_context:    {flags.has_cna_context.sum():,}")


if __name__ == "__main__":
    run()
