import streamlit as st

def render_query_view(genes):
    st.title("🧬 Cell Line Finder")
    
    st.write(
        """
        Identify cell lines that best match your biological
        targets and experimental requirements.
        """
    )

    st.divider()

    with st.form("query_form"):
        st.subheader("1. Define your biological target")
        targets = st.multiselect(
            "Target genes",
            [
                "EGFR", "KRAS", "BRAF", "TP53", 
                "PIK3CA", "ALK", "MET", "ERBB2"
            ],
            default=["EGFR"],
        )

        st.subheader("2. Define exclusion criteria")
        exclusions = st.multiselect(
            "Exclude genes",
            [
                "BRAF", "TP53", "KRAS", 
                "PIK3CA", "ALK", "MET"
            ],
        )

        st.subheader("3. Biological context")
        lineage = st.selectbox(
            "Tissue / lineage",
            [
                "Any", "Lung", "Breast", 
                "Colon", "Liver", "Blood"
            ],
        )

        st.divider()

        st.subheader("Search settings")
        top_k = st.slider(
            "Number of recommendations",
            min_value=3,
            max_value=10,
            value=5,
        )

        st.write("")
        
        # The submit button takes the place of your old button
        submitted = st.form_submit_button(
            "🔍 Find Cell Lines", 
            type="primary", 
            use_container_width=True
        )

        if submitted:
            if not targets:
                st.warning("Please select at least one target gene.")
            else:
                # Save all new parameters to the session state
                st.session_state.query_params = {
                    'targets': targets, 
                    'exclusions': exclusions,
                    'lineage': lineage,
                    'top_k': top_k
                }
                st.session_state.current_step = 'results'
                st.rerun()