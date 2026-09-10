import streamlit as st

def render_about_view():
    st.title("About Cell Line Finder")

    st.write(
        """
        **Cell Line Finder** is a prototype decision-support
        platform for selecting human cell lines using integrated
        multi-omics evidence.
        """
    )

    st.subheader("Purpose")
    st.write(
        """
        Cell line choice can influence assay relevance,
        experimental robustness, and reproducibility.

        The platform aims to integrate molecular evidence from
        public datasets and provide ranked cell-line
        recommendations with interpretable confidence metrics.
        """
    )

    st.subheader("Planned Data Sources")
    st.markdown(
        """
        - DepMap
        - Human Protein Atlas (HPA)
        - GEO
        - RNA expression
        - Protein expression
        - Mutations
        - Gene fusions
        - Molecular signatures
        - miRNA
        - Metabolomics
        """
    )

    st.subheader("Workflow")
    st.code(
        """
User Input
↓
Evidence Retrieval
↓
Evidence Aggregation
↓
Ranking Engine
↓
Confidence Engine
↓
Recommended Cell Lines
↓
Comparison / Alternatives
↓
Export
        """
    )

    st.info(
        "Current version uses mock data. "
        "The next step is connecting the UI to the real "
        "ranking and evidence pipeline."
    )