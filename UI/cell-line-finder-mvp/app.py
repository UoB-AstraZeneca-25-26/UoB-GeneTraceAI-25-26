import streamlit as st
from src.data_loader import load_mock_data
from ui.query_builder import render_query_view
from ui.results_table import render_results_view
from ui.profile_view import render_profile_view
from ui.about_page import render_about_view

# Page configuration must be the first Streamlit command
st.set_page_config(page_title="Cell Line Finder", layout="wide", page_icon="🧬")

def init_session_state():
    """Initializes global state variables for routing and data sharing."""
    if 'current_step' not in st.session_state: 
        st.session_state.current_step = 'query'
    if 'query_params' not in st.session_state: 
        st.session_state.query_params = {}
    if 'selected_cell_line' not in st.session_state: 
        st.session_state.selected_cell_line = None

def main():
    init_session_state()
    
    # Load data globally so it's available to all views
    raw_df = load_mock_data()
    all_genes = sorted(raw_df["Gene"].unique().tolist())
    
    # Global Sidebar
    with st.sidebar:
        st.header("Control Panel")
        if st.button("🔄 Start New Query", use_container_width=True):
            st.session_state.current_step = 'query'
            st.session_state.query_params = {}
            st.session_state.selected_cell_line = None
            st.rerun()

        if st.button("ℹ️ About", use_container_width=True):
            st.session_state.current_step = 'about'
            st.rerun()
            
    # View Router
    if st.session_state.current_step == 'query':
        render_query_view(all_genes)
    elif st.session_state.current_step == 'results':
        render_results_view(raw_df)
    elif st.session_state.current_step == 'profile':
        render_profile_view(raw_df)
    elif st.session_state.current_step == 'about':
        render_about_view()

if __name__ == "__main__":
    main()