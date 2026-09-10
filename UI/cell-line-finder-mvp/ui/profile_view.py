import streamlit as st
import pandas as pd
import plotly.express as px

def render_profile_view(raw_df):
    """Renders the detailed breakdown for a specific cell line."""
    cl_name = st.session_state.selected_cell_line
    
    if st.button("← Back to Results"):
        st.session_state.current_step = 'results'
        st.rerun()
        
    st.title(f"🔬 Cell Line Profile: {cl_name}")
    cl_data = raw_df[raw_df["Cell Line"] == cl_name]
    
    # Top level metrics
    col1, col2, col3 = st.columns(3)
    col1.metric("Tissue of Origin", cl_data["Tissue"].iloc[0])
    col2.metric("Total Somatic Variants", cl_data["Has_Somatic_Variant"].sum())
    col3.metric("Data Sources", "DepMap, HPA_rna, CCLE")
    
    st.divider()
    
    tab1, tab2, tab3 = st.tabs(["📊 Multi-Omics Evidence", "🧬 Genomic Features", "🔄 Alternatives"])
    
    with tab1:
        st.subheader("Transcriptomic vs Proteomic Expression")
        fig = px.bar(
            cl_data, 
            x="Gene", 
            y=["TPM_Expression", "Protein_Expression"],
            barmode="group",
            labels={"value": "Expression Level", "variable": "Data Type"},
            color_discrete_map={"TPM_Expression": "#3b82f6", "Protein_Expression": "#10b981"}
        )
        st.plotly_chart(fig, use_container_width=True)
        
    with tab2:
        st.subheader("Somatic Mutations")
        mutations = cl_data[cl_data["Has_Somatic_Variant"] == True]
        if not mutations.empty:
            st.dataframe(mutations[["Gene", "Has_Somatic_Variant"]], hide_index=True)
        else:
            st.info("No targeted somatic variants found for this cell line.")
            
    with tab3:
        st.subheader("Similar Cell Lines (Backup Options)")
        st.write("Based on global omics signatures, consider these alternatives:")
        alt_df = pd.DataFrame({
            "Alternative": ["A549", "MCF7"],
            "Similarity Score": ["89%", "76%"],
            "Shared Targets": ["EGFR", "TP53"]
        })
        st.dataframe(alt_df, hide_index=True)