import pandas as pd

def run_scoring_pipeline(df, targets, exclusions, lineage="Any", top_k=5):
    """Filters data and calculates confidence scores for ranking."""
    
    # 1. Filter by lineage (Tissue) before doing any math
    if lineage != "Any":
        df = df[df["Tissue"] == lineage]
        
    # 2. Filter for targets
    target_df = df[df["Gene"].isin(targets)].copy()
    
    # 3. Aggregate target data per cell line
    grouped = target_df.groupby(["Cell Line", "Tissue"]).agg(
        Mean_TPM=('TPM_Expression', 'mean'),
        Mutations=('Has_Somatic_Variant', 'sum')
    ).reset_index()
    
    # 4. Calculate exclusions and final score
    exclusion_df = df[df["Gene"].isin(exclusions)]
    for idx, row in grouped.iterrows():
        cl = row["Cell Line"]
        ex_exp = exclusion_df[exclusion_df["Cell Line"] == cl]["TPM_Expression"].mean()
        if pd.isna(ex_exp): 
            ex_exp = 0
            
        score = min(100, max(0, (row["Mean_TPM"] / 1.5) - (ex_exp / 2)))
        grouped.at[idx, "Confidence Score (%)"] = round(score, 1)
        grouped.at[idx, "Exclusion Penalty"] = round(ex_exp, 2)
        
    # 5. Sort by score and enforce the Top K limit
    ranked_df = grouped.sort_values(by="Confidence Score (%)", ascending=False)
    return ranked_df.head(top_k)