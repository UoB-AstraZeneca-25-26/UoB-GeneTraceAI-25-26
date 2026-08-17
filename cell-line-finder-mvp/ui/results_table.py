import streamlit as st
from src.scoring import run_scoring_pipeline

def render_results_view(raw_df):
    """Renders the ranked recommendations and selection panel."""
    st.title("📊 Ranked Recommendations")
    
    # Extract all parameters from the query
    targets = st.session_state.query_params.get('targets', [])
    exclusions = st.session_state.query_params.get('exclusions', [])
    lineage = st.session_state.query_params.get('lineage', 'Any')
    top_k = st.session_state.query_params.get('top_k', 5)
    
    # Update the caption to show the current lineage filter
    st.caption(f"**Targets:** {', '.join(targets)} | **Exclusions:** {', '.join(exclusions) if exclusions else 'None'} | **Lineage:** {lineage}")
    
    # Pass the new arguments into the scoring pipeline
    results_df = run_scoring_pipeline(raw_df, targets, exclusions, lineage, top_k)
    
    # Safety check: if no cell lines match the filters
    if results_df.empty:
        st.warning("No cell lines matched your specific criteria. Try adjusting your lineage or target genes.")
        return
    
    col1, col2 = st.columns([3, 1])
    
    with col1:
        st.dataframe(
            results_df[['Cell Line', 'Tissue', 'Confidence Score (%)', 'Mean_TPM', 'Mutations']],
            use_container_width=True,
            hide_index=True
        )
        
    with col2:
        st.markdown("### Deep Dive")
        st.markdown("Select a cell line to view integrated evidence.")
        selected = st.selectbox("Choose Cell Line", results_df["Cell Line"].tolist())
        
        if st.button("Inspect Profile", type="primary", use_container_width=True):
            st.session_state.selected_cell_line = selected
            st.session_state.current_step = 'profile'
            st.rerun()