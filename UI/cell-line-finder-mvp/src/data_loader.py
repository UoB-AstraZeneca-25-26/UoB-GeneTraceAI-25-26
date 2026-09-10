import streamlit as st
import pandas as pd
import numpy as np

@st.cache_data
def load_mock_data():
    """Generates and caches mock multi-omics data for the MVP."""
    np.random.seed(42)
    cell_lines = ["HCC827", "A549", "MCF7", "HeLa", "PC3", "K562", "Jurkat", "HepG2", "U87", "HT29"]
    tissues = ["Lung", "Lung", "Breast", "Cervix", "Prostate", "Blood", "Blood", "Liver", "Brain", "Colon"]
    genes = ["EGFR", "BRCA1", "TP53", "KRAS", "PTEN", "MYC", "BCR::ABL1"]
    
    data = []
    for i, cl in enumerate(cell_lines):
        for gene in genes:
            tpm = np.random.uniform(0, 150) if gene != "BCR::ABL1" else (200 if cl == "K562" else 0)
            protein = tpm * np.random.uniform(0.5, 1.2) 
            has_mutation = np.random.choice([True, False], p=[0.2, 0.8])
            
            data.append({
                "Cell Line": cl,
                "Tissue": tissues[i],
                "Gene": gene,
                "TPM_Expression": round(tpm, 2),
                "Protein_Expression": round(protein, 2),
                "Has_Somatic_Variant": has_mutation
            })
            
    return pd.DataFrame(data)