# Cell Line Finder (MVP)

An integrated multi-omics cell line selection platform built with Streamlit. This MVP allows users to specify molecular targets and exclusion criteria to identify optimal human cell lines, providing a ranked recommendation with confidence metrics.

## Folder Structure Overview

```
cell-line-finder-mvp/
├── README.md
├── requirements.txt
├── app.py
├── src/
│   ├── __init__.py
│   ├── data_loader.py
│   └── scoring.py
└── ui/
    ├── __init__.py
    ├── query_builder.py
    ├── results_table.py
    └── profile_view.py
```

(Note: Create the ```__init__.py``` files as empty files. They tell Python to treat the `src` and `ui` directories as importable modules.)

## Architecture
The application uses a layered architecture to keep the Streamlit codebase clean and maintainable:
- **`app.py`**: The main router and state manager.
- **`src/`**: Pure Python business logic and data handling.
- **`ui/`**: Modular Streamlit UI components for rendering specific views.

## Setup & Installation

1. Create a virtual environment (recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

2. Install the required dependencies
    ```bash
    pip install -r requirements.txt
    ```

3. Run the application
    ```bash
    streamlit run app.py
    ```