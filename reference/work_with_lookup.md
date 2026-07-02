## What are the lookup tables?

There are two shared tables sitting in the `reference/` folder:

- `gene_lookup.parquet` — a dictionary of ~19,446 human protein-coding genes
- `cell_line_lookup.parquet` — a dictionary of 1,840 cancer cell lines

Think of them as the **single source of truth** for how genes and cell lines are identified across this project. Every dataset in this project uses different identifiers — some use gene symbols, some use Ensembl IDs, some use cell line names, some use accession numbers. The lookup tables resolve all of that into two canonical keys that everyone agrees on:

- **`ensg_id`** — bare Ensembl gene ID (e.g. `ENSG00000141510`) — the canonical gene key
- **`model_id`** — DepMap ACH- ID (e.g. `ACH-000016`) — the canonical cell line key

Your job when cleaning your datasets is to **add these two columns** to your data. Once every dataset has `ensg_id` and `model_id`, everything can be joined together.

---

## Loading the tables

Always load from the `reference/` folder at the start of your cleaning notebook:

```python
import pandas as pd
import os

REF = "reference/"  # adjust path to match your setup

gene_lookup      = pd.read_parquet(os.path.join(REF, "gene_lookup.parquet"))
cell_line_lookup = pd.read_parquet(os.path.join(REF, "cell_line_lookup.parquet"))
```

---

## What's inside each table

### gene_lookup

| Column | What it contains |
|---|---|
| `ensg_id` | Bare Ensembl ID — **primary key** |
| `hgnc_symbol` | Current approved gene symbol (e.g. `TP53`) |
| `prev_symbols` | Old symbols, pipe-delimited (e.g. `LFS1\|P53`) |
| `alias_symbols` | Alternative names, pipe-delimited |
| `hgnc_id` | HGNC accession (e.g. `HGNC:11998`) |
| `biotype` | Always `protein_coding` |

### cell_line_lookup

| Column | What it contains |
|---|---|
| `model_id` | ACH- DepMap ID — **primary key** |
| `cell_line_name` | Full cell line name (e.g. `MCF7`) |
| `stripped_cell_line_name` | Simplified name, no hyphens/spaces (e.g. `MCF7`) |
| `cvcl_accession` | Cellosaurus accession (e.g. `CVCL_0031`) |
| `synonyms` | Alternative names, pipe-delimited |
| `rrid` | RRID identifier (e.g. `CVCL_0031`) |

---

## Using gene_lookup in your dataset

### Situation A — Your dataset already has bare ENSG IDs

This is the simplest case. Just validate that your IDs exist in the lookup and rename the column:

```python
# Keep only genes present in the lookup (filters out non-protein-coding)
df = df[df["your_gene_column"].isin(gene_lookup["ensg_id"])].copy()

# Rename to canonical column name
df = df.rename(columns={"your_gene_column": "ensg_id"})
```

### Situation B — Your dataset has ENSG IDs with version suffixes (e.g. ENSG00000141510.15)

Strip the suffix first, then validate:

```python
df["ensg_id"] = df["your_gene_column"].str.replace(r"\.\d+$", "", regex=True)
df = df[df["ensg_id"].isin(gene_lookup["ensg_id"])].copy()
```

### Situation C — Your dataset has gene symbols instead of ENSG IDs

Build a reverse lookup from symbol → ensg_id and map across:

```python
symbol_to_ensg = dict(zip(gene_lookup["hgnc_symbol"], gene_lookup["ensg_id"]))

df["ensg_id"] = df["your_symbol_column"].map(symbol_to_ensg)

# Log anything that didn't resolve
unmapped_genes = df[df["ensg_id"].isna()]["your_symbol_column"]
unmapped_genes.to_csv("logs/unmapped_genes_yourdataset.csv", index=False)
```

---

## Using cell_line_lookup in your dataset

Cell line identifiers vary more than gene identifiers across datasets. There are five resolution situations — find which one applies to your dataset and follow that pattern.

---

### Chain 1 — Your dataset already has ACH- model IDs

The simplest case. Just validate:

```python
valid_ids = set(cell_line_lookup["model_id"])

# Flag any that don't match
orphaned = df[~df["your_id_column"].isin(valid_ids)]
orphaned[["your_id_column"]].to_csv("logs/unmapped_cells_yourdataset.csv", index=False)

# Keep matched rows and rename
df = df[df["your_id_column"].isin(valid_ids)].copy()
df = df.rename(columns={"your_id_column": "model_id"})
```

**Applies to:** proteomics, metabolomics, fusions, signatures

---

### Chain 2 — Your dataset has PR- profile IDs

PR- IDs are profile-level, not cell-line-level. You need `depmap_profiles` to bridge them to ACH- IDs first:

```python
depmap_profiles = pd.read_parquet("path/to/8_DepMap_OmicsProfiles.parquet")

# For RNA datasets: filter to rna datatype
# For mutation datasets: filter to wes or wgs
rna_lookup = (depmap_profiles[depmap_profiles["Datatype"] == "rna"]
              .set_index("ProfileID")["ModelID"]
              .to_dict())

df["model_id"] = df["your_profile_id_column"].map(rna_lookup)

# Log unresolved
unmapped = df[df["model_id"].isna()]
unmapped[["your_profile_id_column"]].to_csv("logs/unmapped_cells_yourdataset.csv", index=False)
```

**Applies to:** depmap_expr (use `rna`), mutations (use `wes`/`wgs`)

---

### Chain 3 — Your dataset has cell line names

Try matching in order — exact name first, then stripped name, then CVCL:

```python
# Build dictionaries for each fallback step
name_to_id     = dict(zip(cell_line_lookup["cell_line_name"],
                          cell_line_lookup["model_id"]))
stripped_to_id = dict(zip(cell_line_lookup["stripped_cell_line_name"],
                          cell_line_lookup["model_id"]))
rrid_to_id     = dict(zip(cell_line_lookup["rrid"],
                          cell_line_lookup["model_id"]))

# Step 1: exact name match
df["model_id"] = df["your_name_column"].map(name_to_id)

# Step 2: stripped name match for anything still null
still_null = df["model_id"].isna()
df.loc[still_null, "model_id"] = df.loc[still_null, "your_name_column"].map(stripped_to_id)

# Step 3: CVCL match via hpa_desc (for HPA datasets only)
# see hpa_desc["Cellosaurus ID"] → map to rrid_to_id

# Log permanently unresolved
unmapped = df[df["model_id"].isna()]
unmapped[["your_name_column"]].to_csv("logs/unmapped_cells_yourdataset.csv", index=False)
```

**Applies to:** hpa_rna, hpa_desc  
**Expected recovery rate:** ~91.5% (102 permanently unresolvable from EDA)

---

### Chain 4 — Your dataset has GSM accession IDs as columns

GSM IDs need to be mapped through geo_info to get cell line names, then resolved:

```python
geo_info = pd.read_parquet("path/to/10_GEOInfo.parquet")

# Build GSM → model_id via cellosaurus_id
gsm_to_cvcl = dict(zip(geo_info["Geo_accession"], geo_info["Cellosaurus_ID"]))
rrid_to_id  = dict(zip(cell_line_lookup["rrid"], cell_line_lookup["model_id"]))

gsm_cols = [c for c in df.columns if c.startswith("GSM")]

gsm_to_model = {}
unresolved_gsm = []

for gsm in gsm_cols:
    cvcl = gsm_to_cvcl.get(gsm)
    mid  = rrid_to_id.get(cvcl) if cvcl else None
    if mid:
        gsm_to_model[gsm] = mid
    else:
        unresolved_gsm.append(gsm)

# Log unresolved
pd.Series(unresolved_gsm).to_csv("logs/unmapped_cells_geo_expr.csv", index=False)
```

**Applies to:** geo_expr  
**Expected:** all 3,267 GSMs found in geo_info; 0 unmapped at GSM step

---

### Chain 5 — Your dataset has CCLE-style column headers

CCLE names follow the pattern `CELLNAME_TISSUE` (e.g. `MCF7_BREAST`). Strip the tissue suffix and match:

```python
stripped_to_id = dict(zip(cell_line_lookup["stripped_cell_line_name"],
                          cell_line_lookup["model_id"]))

id_cols = ["Name", "Description"]  # columns that are NOT cell lines
cell_cols = [c for c in df.columns if c not in id_cols]

col_to_model = {}
unresolved = []

for col in cell_cols:
    base = col.split("_")[0].lower()
    mid  = stripped_to_id.get(base)
    if mid:
        col_to_model[col] = mid
    else:
        unresolved.append(col)

# Log unresolved column names
pd.Series(unresolved).to_csv("logs/unmapped_cells_mirna.csv", index=False)
```

**Applies to:** mirna  
**Note:** Match rate not yet quantified — run and record before Week 3 handoff

---

## The unmapped log — non-negotiable

Every cleaning function must write unresolved rows to a log. Never drop silently.

```python
os.makedirs("logs", exist_ok=True)

unmapped.to_csv("logs/unmapped_{your_dataset_name}.csv", index=False)
```

Log columns should include: `original_value`, `dataset`, `resolution_step_reached`

---

## What your cleaned dataset should look like at the end

Every dataset, regardless of source, should come out of cleaning with at minimum:

```
ensg_id          | model_id    | [your dataset columns]
ENSG00000141510  | ACH-000016  | value=5.2
ENSG00000141510  | ACH-000032  | value=3.1
```

Datasets without a gene dimension (e.g. signatures, metabolomics) only need `model_id`.  
Datasets without a cell line dimension (none in this project) would only need `ensg_id`.

Once your dataset has these canonical keys, it can be joined to any other dataset in the project without any further ID resolution.

---

## Quick reference — which chain applies to your dataset

| Dataset | Gene resolution | Cell line resolution |
|---|---|---|
| depmap_expr | Extract ENSG from `SYMBOL (ENSG...)` headers | Chain 2 — PR- via depmap_profiles |
| hpa_rna | Bare ENSG — Situation A | Chain 3 — cell line names |
| hpa_desc | No gene column | Chain 3 — cell line names |
| geo_expr | Bare ENSG — Situation A | Chain 4 — GSM accessions |
| mutations | Bare ENSG in `EnsemblGeneID` — Situation A | Chain 2 — PR- via depmap_profiles |
| fusions | ENSG with version suffix — Situation B | Chain 1 — direct ACH- |
| signatures | No gene column | Chain 1 — direct ACH- |
| mirna | No ENSG (miRNA IDs) | Chain 5 — CCLE column headers |
| proteomics | UniProt IDs (not ENSG) | Chain 1 — direct ACH- |
| metabolomics | No gene column | Chain 1 — direct ACH- |
| geo_info | No gene column | Chain 4 — GSM accessions |