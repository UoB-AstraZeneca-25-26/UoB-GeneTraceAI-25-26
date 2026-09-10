"""
Ranking/evidence_ledger.py
---------------------------
Shared derivation logic for evidence fields. Used by both build_evidence_ledger.py
(bulk) and explain_pair.py (single-pair live query) so the two never drift apart.
"""
import pandas as pd
import numpy as np


def derive_layers_present(row: pd.Series) -> list[str]:
    layers = ["rna"]
    if pd.notna(row.get("prot_z_t")):
        layers.append("protein")
    if row.get("p_mutation", 0) > 0:
        layers.append("mutation")
    if row.get("p_fusion", 0) > 0:
        layers.append("fusion")
    if row.get("has_cna_alteration", False):
        layers.append("cna")
    return layers


def derive_driver_alteration(row: pd.Series) -> bool:
    return bool(row.get("has_driver_alteration", False))


def derive_cna_role_consistent(row: pd.Series) -> bool:
    role = str(row.get("gene_role", ""))
    has_amp = row.get("has_cna_alteration", False)
    if "oncogene" in role and has_amp:
        return True
    if "tsg" in role and not has_amp:
        return True
    return False


def derive_essentiality_check(row: pd.Series) -> str:
    score = row.get("core_score", 0.5)
    if score >= 0.8:
        return "likely_essential"
    if score >= 0.5:
        return "possibly_essential"
    return "likely_non_essential"


def derive_census_role(row: pd.Series, gene_lookup: pd.DataFrame) -> str:
    ensg = row.get("ensg_id", "")
    match = gene_lookup[gene_lookup["ensg_id"] == ensg]
    if match.empty:
        return "unknown"
    return str(match.iloc[0].get("gene_role", "unknown"))


def build_ledger_row(pred: pd.DataFrame, core: pd.DataFrame,
                     cna: pd.DataFrame, genes: pd.DataFrame,
                     lines: pd.DataFrame) -> pd.DataFrame:
    df = pred.copy()
    if "has_cna_alteration" not in df.columns:
        df = df.merge(cna[["ensg_id","model_id","has_cna_alteration"]],
                      on=["ensg_id","model_id"], how="left")
    gene_meta = genes[["ensg_id","hgnc_symbol","gene_role"]].rename(
        columns={"hgnc_symbol":"gene_symbol"})
    if "gene_role" not in df.columns:
        df = df.merge(gene_meta, on="ensg_id", how="left")

    df["has_cna_alteration"] = df["has_cna_alteration"].fillna(False)
    df["p_mutation"]         = df.get("p_mutation", pd.Series(0.0, index=df.index)).fillna(0.0)
    df["p_fusion"]           = df.get("p_fusion", pd.Series(0.0, index=df.index)).fillna(0.0)

    df["layers_present"]   = df.apply(derive_layers_present, axis=1).apply(str)
    df["driver_alteration"]= df.apply(derive_driver_alteration, axis=1)
    df["essentiality"]     = df.apply(derive_essentiality_check, axis=1)

    return df
