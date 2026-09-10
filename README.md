# GeneTraceAI

## Folder structure

```
UoB-GeneTraceAI-25-26/
│
├── architecture/               — the production pipeline (renamed from final_pipeline)
│   ├── 00_harmonisation/          Stage 0: builds the DuckDB warehouse (celllineselector.db)
│   ├── 01_Transcriptomics/        Stage 1: RNA scoring
│   ├── 02_Proteinomics/           Stage 2: protein scoring
│   ├── 03_Altercations/           Stage 3: CNA / fusion / mutation scoring
│   ├── Scoring/                   core_score, confidence tiers, driver routing
│   ├── Ranking/                   ranking CLI and evidence-ledger logic
│   ├── api/                       FastAPI app served in production (app.py, routes.py, ranking.py, bedrock.py)
│   ├── utils/                     shared helpers (common.py, robust_z_matrix, etc.)
│   ├── reference/                 lookup exports used by the live API
│   ├── results/                   pipeline run outputs
│   ├── cell_similarity/           RNA-similarity clustering (k-NN neighbours used by /gene/detail)
│   ├── expression_fusion/         expression + fusion combination work
│   ├── docs/                      method specs, scoring/ranking declarations, math reference, diagrams
│   ├── config.py, run_all.py, Dockerfile, requirements.txt
│   └── harvest_geo_meta.py, audit_warehouse.py, apply_transcriptomics_renames.py
│
├── Testing and validation/     — everything that checks the pipeline's output, not part of it
│   ├── Testing/                   run_tests.py
│   ├── Validation/                eval.py (stage-4 held-out evaluation)
│   ├── gate_audit/                six-audit re-examination of the depletion gate
│   ├── diagnostics/                ~50 standalone diagnostic/audit scripts (T*, V*, W*, X* series)
│   └── Confidence intervals/      Monte Carlo confidence-interval + rank-stability scoring
│
├── AIAgent/                    — the LLM agent service
│   ├── api/                       FastAPI app, routes, settings, security
│   ├── AgentDevelopment/          LLM provider wiring, function calling, gene index
│   ├── Tests/                     unit + lambda tests
│   ├── Schema/                    prompt schemas
│   ├── Knowledge/                 math_reference.md (source-of-truth citations for agent answers)
│   ├── lambda_handler.py          AWS Lambda entrypoint (Mangum wrapping the FastAPI app)
│   ├── Dockerfile, requirement_agents.txt, requirements.lock.txt
│   └── CLEANUP_AUDIT.md, agent_api.md, agent_details.md, plan.md
│
├── UI/                          — all three front-end implementations
│   ├── cell-line-finder-react/    the primary React/Vite web app
│   ├── cell-line-finder-mvp/      Streamlit prototype UI
│   └── genetraceai/                third, separate React app
│
├── EDA/                        — exploratory data analysis
│   ├── BasicEDA.ipynb
│   └── Initial Data Audit and Mapping/   data audit/mapping notebooks + resolution logs
│
├── Cell Line Selector/          — Triveni's original harmonisation/transcriptomics/proteomics work
│   ├── pipelines/                  numbered pipeline stages (00-06)
│   ├── src/scripts/                cleaning, harmonisation, transcriptomics, proteomics modules
│   ├── figs/                       pipeline diagrams
│   └── test_transcriptomics_gene_filter.py, transcriptomics_stale.py, transcriptomics_z_audit.py
│
├── reference/                   — canonical lookup tables (gene_lookup, cell_line_lookup, resolution logs)
├── research-evidence/           — reference.bib
│
└── README.md, pyrightconfig.json, missing_protein_coding_585_buckets.csv, merge_probe.txt
```

### Where each contributor's work landed

| Contributor | Work | Location |
|---|---|---|
| Chai | Production pipeline, AI agent, React frontend | `architecture/`, `AIAgent/`, `UI/cell-line-finder-react/` |
| Triveni | Original harmonisation, transcriptomics, proteomics pipeline | `Cell Line Selector/`, proteomics files in `architecture/02_Proteinomics/` |
| Musa | Streamlit MVP, Monte Carlo confidence intervals | `UI/cell-line-finder-mvp/`, `Testing and validation/Confidence intervals/` |
| Thiruvel | AWS Lambda deployment, cleanup audit | `AIAgent/lambda_handler.py`, `AIAgent/CLEANUP_AUDIT.md` |
| Team (gate/diagnostics) | Gate mechanism audit, ~50 diagnostic scripts | `Testing and validation/gate_audit/`, `Testing and validation/diagnostics/` |
