# GeneTraceAI — The DuckDB Query Layer

How every dataset in this repository is exposed as a SQL-queryable catalog, why
that layer exists, and what it deliberately does *not* change.

Companion documents: [GLOSSARY.md](GLOSSARY.md) (what the quantities mean),
[MATH_REFERENCE.md](MATH_REFERENCE.md) (the estimators), and
[../integration_contract.md](../integration_contract.md) (the key and schema
contract the catalog inherits).

**All timings and sizes below were measured on this machine against this
repository's own tables**, not quoted from documentation. They were taken on
2026-08-01; the command that produces each is given so they can be re-measured.

---

## 1. What DuckDB is

DuckDB is an **embedded analytical database**. Two words carry the meaning:

- **Embedded** — it is a library, not a server. `pip install duckdb`, `import
  duckdb`, done. There is no daemon to start, no port, no user accounts, no
  connection string. In this respect it is like SQLite; the whole database is
  one file on disk, or no file at all if you work purely in memory.
- **Analytical** — it is built for scanning and aggregating large columnar
  tables (the kind of query that reads millions of rows and returns a summary),
  rather than for the small single-row reads and writes a transactional database
  like PostgreSQL optimises for. It stores and processes data column by column
  and executes queries in parallel across CPU cores.

The property that matters most for this project is that **DuckDB reads Parquet
files directly as if they were tables**. It does not need the data imported,
copied, or converted first. `SELECT * FROM read_parquet('core_score.parquet')`
works on a file this pipeline already produces, unchanged.

That is the whole basis of the design below: this repository is already a
collection of Parquet files, so it is already a DuckDB database. The catalog
just gives those files names.

---

## 2. Why the pipeline needs it

Before this layer, every consumer of a table did the same thing: load the entire
Parquet file into a pandas DataFrame with `pd.read_parquet`, then filter it in
Python. That is fine for a 141-row table and increasingly wasteful for the large
ones.

`src/pipeline/outputs/core_score.parquet` holds **28,403,110 rows × 5 columns** —
one row per (gene, cell line) pair. Answering *"what are the top cell lines for
BRAF?"* touches 1,485 of those rows. Measured cost of the two routes:

| Route | Time | Peak memory |
|---|---|---|
| `pd.read_parquet` (full load, then filter) | 3.01 s | 4.18 GB |
| DuckDB filtered scan | 0.06 s | negligible |

```bash
python -c "import time,pandas as pd; t=time.time(); df=pd.read_parquet('src/pipeline/outputs/core_score.parquet'); print(time.time()-t, df.shape, df.memory_usage(deep=True).sum()/1e9)"
```

The 50× gap is not clever engineering — it is that DuckDB reads only the column
chunks whose Parquet metadata says they can contain the requested gene, and
never materialises the other 28.4 million rows. The 4.18 GB figure is the more
important one: it is what forces the notebooks to work one table at a time, and
it is why joining three large tables in pandas is awkward on a laptop.

The second motivation is **joins**. This project's whole structure is a join
problem — `ensg_id` for genes, `model_id` for cell lines, across five zones of
tables (see [integration_contract.md](../integration_contract.md) §3–4). Those
joins are currently written as chains of `pd.merge` calls whose intermediate
results all live in memory. In SQL they are one statement, and the database
chooses the execution order.

---

## 3. The catalog

`src/scripts/duckdb_catalog.py` registers **87 views** over every Parquet and CSV
file in the repository, grouped into five schemas — one per zone of the pipeline:

| Schema | Directory | Views | Holds |
|---|---|---:|---|
| `raw` | `data/parquet/raw_data` | 14 | the source extracts, uncleaned |
| `clean` | `data/parquet/data_clean` | 15 | per-dataset cleaned tables |
| `track` | `cleaned_track_data` | 15 | Track A–D outputs |
| `ref` | `reference` | 9 | lookup tables — the join anchors |
| `out` | `src/pipeline/outputs` | 34 | Stage 0–7 scoring outputs, incl. `variants/` |

So `ref.gene_lookup` is the gene lookup table, `out.core_score` is the Stage 2
score, `raw.depmap_expr` is the untouched DepMap expression extract. The zone
prefix makes the provenance of a table visible at the point of use, which
matters here because several names recur across zones — `sample_info`,
`proteomics` and `metabolomics` each exist in three different zones at three
different stages of cleaning, and `clean.proteomics` is emphatically not
`track.proteomics`.

### Views, not copies

Every registered object is a **view**, defined as:

```sql
CREATE OR REPLACE VIEW ref.gene_lookup AS
SELECT * FROM read_parquet('.../reference/gene_lookup.parquet');
```

Nothing is imported. The catalog file `genetrace.duckdb` is **7.2 MB** of view
definitions against 10.6 GB of Parquet and CSV; those files remain the single
source of truth, and the catalog is a disposable index over them. It is
`.gitignore`d for exactly that reason — it is rebuilt, never authored.

The consequence worth internalising: **a view is re-bound by DuckDB on every
query**, so when a pipeline stage rewrites its Parquet output, the next query
sees the new data with no re-import and no rebuild step. This was verified
rather than assumed — a file rewritten with different rows *and an added column*
was reflected immediately in the existing view. The catalog only falls out of
step when a file is **added or deleted**, which is the one case
`build_catalog()` handles by re-registering and dropping orphans.

### Naming

The 14 raw extracts carry ordinal filenames (`2_DepMap_OmicsExpressionAllGenes
TPMLogp1Profile.parquet`), which are not usable as SQL identifiers and not
readable either. Those are mapped to the same short names the existing
`data_utils.load_raw_parquets` already uses — `depmap_expr`, `sample_info`,
`signatures` — so a name means the same thing in both access paths. Everything
else is auto-named from its filename (lowercased, non-alphanumerics collapsed to
`_`), with the redundant `_clean` suffix dropped inside the `clean` schema, and
files in `outputs/variants/` prefixed `variants_`. A new file dropped into any
zone is picked up automatically; no registration list needs editing.

---

## 4. Using it

### From a notebook or script

```python
import sys; sys.path.insert(0, 'src/scripts')
from duckdb_catalog import q, tables, schema_of

tables()                      # every view with its row and column count
schema_of('out.core_score')   # columns and types

q("""
  SELECT g.hgnc_symbol, c.cell_line_name, c.model_id, s.core_score
  FROM out.core_score s
  JOIN ref.gene_lookup      g USING (ensg_id)
  JOIN ref.cell_line_lookup c USING (model_id)
  WHERE g.hgnc_symbol = $sym
  ORDER BY s.core_score DESC
  LIMIT 10
""", sym='BRAF')
```

`q()` returns an ordinary **pandas DataFrame**, so everything downstream —
matplotlib, existing helper functions, `.to_parquet()` — is unchanged. DuckDB
replaces the loading and joining, not the analysis.

Pass values as named parameters (`$sym`) rather than formatting them into the
SQL string. It is the habit that keeps a gene symbol from ever being able to
alter the statement.

That three-table join over 28.4M rows runs in **1.21 s** cold, **0.20 s** warm.

### From the command line

```bash
python src/scripts/duckdb_catalog.py --list
```

```bash
python src/scripts/duckdb_catalog.py --schema out.evidence_ledger
```

```bash
python src/scripts/duckdb_catalog.py --sql "SELECT class, count(*) FROM out.gene_regime GROUP BY 1 ORDER BY 2 DESC"
```

```bash
python src/scripts/duckdb_catalog.py --build
```

`--build` is the only one that writes; run it after adding or deleting a data
file. Everything else reads.

### Connection behaviour

`connect()` defaults to `refresh="auto"`: it compares the files on disk against
the registered views and rebuilds **only** on a mismatch. Measured — an
unchanged repository connects in **0.097 s**; a genuine rebuild costs **11.7 s**,
because registering 87 views means reading 87 Parquet schemas, and two of those
are very wide (`clean.depmap_expr` has 53,962 columns, `clean.proteomics`
12,559). Hence the check rather than an unconditional rebuild.

`connect(read_only=True)` opens without writing — necessary because **DuckDB
permits one writer at a time**, so a second process (or a second notebook
kernel) must open read-only while the first holds the file.

---

## 5. What this does not change

This layer is **purely additive**. Nothing in the existing pipeline was
rewritten to use it:

- `data_utils.load_raw_parquets` / `load_clean_parquets` still work exactly as
  before, and the notebooks that call them are untouched.
- No Stage 0–7 script or notebook has a new dependency.
- No Parquet file was moved, renamed, or rewritten.
- The canonical keys and schema rules in
  [integration_contract.md](../integration_contract.md) are unaffected — the
  catalog exposes the tables, it does not define them. `ensg_id` and `model_id`
  remain the join keys; names remain display labels.

The catalog is therefore reversible: deleting `genetrace.duckdb` and
`duckdb_catalog.py` returns the repository to its previous state exactly.

### Known limits

- **One writer.** Two notebook kernels cannot both hold the catalog for writing.
  Use `read_only=True` in the second, or let each build its own catalog file.
- **Wide tables are awkward in SQL.** `clean.depmap_expr` (53,962 columns) and
  `clean.proteomics` (12,559) are stored gene-per-column. SQL is built for
  many rows, not many columns, so `SELECT *` on these is unpleasant and
  selecting a gene means quoting a column name. These tables are better reached
  through the long-form harmonised outputs in `out.`, which are already one row
  per (gene, cell line).
- **The catalog is not a schema.** It exposes what exists, including dead
  branches — `out.layer3_*` are the rejected Chronos-as-third-layer tables and
  are registered simply because the files are still on disk. Registration is not
  endorsement; see [MATH_REFERENCE.md](MATH_REFERENCE.md) and `DECISIONS.md` for
  which outputs are live.
- **`reference/gene_symbol_resolver.json`** is not registered. It is a 62,209-key
  nested mapping, not a tabular file; it stays a JSON load.

---

## 6. A second DuckDB surface: the harmonisation warehouse

Everything above describes the **catalog** — a set of names registered over
parquet files that already exist, holding no data of its own. As of 2026-08-03
there is a second, different DuckDB surface in the repo, and the two should not
be confused.

`src/pipeline/outputs/celllineselector.db` is a **materialised warehouse**: 28
tables physically stored inside a 1.5 GB database file, written by
`00_harmonisation.ipynb` and `00b_enriched_harmonisation.ipynb`. See
[HARMONISATION_V2.md](HARMONISATION_V2.md).

| | catalog (`src/scripts/duckdb_catalog.py`) | warehouse (`celllineselector.db`) |
|---|---|---|
| holds data | no — views over parquet | yes — tables on disk |
| built by | one script, seconds | two notebooks, ~5 min |
| covers | the whole repo's parquet files | Stage 0 inputs and identity/gene dimensions |
| source of truth | the parquet files | the parquet exports it writes |

The warehouse exists because harmonisation needs SQL *while it runs*, not just to
query afterwards. Its gene dimension is assembled from five sources whose gene
keys sit in different shapes — column names in the wide expression matrices,
column values in the long ones, two columns per row in `fusions` — and
`list_sort(list_distinct(list(...)))` over a `UNION ALL` expresses that in a way
a pandas `groupby` chain does not.

Both surfaces are gitignored and both are regenerated from source. Neither is a
source of truth: the parquet files are.

---

## 7. Where this could go next

Recorded as options, not commitments — none of these have been done:

1. **Migrate one stage's joins to SQL.** `01_lookup_metadata_join` is the natural
   candidate: it is almost entirely `pd.merge` chains against `ref.` tables, so
   the translation is mechanical and the result is checkable against the existing
   output row-for-row.
2. **Back the Stage 7 query layer with SQL.** `rank_cell_lines.py` and
   `explain_pair.py` currently load `predictions_with_confidence.parquet`
   (28.4M rows, 15 columns) to answer single-gene questions. This is the place
   where the 0.06 s vs 3.01 s gap is felt by a user waiting on the CLI.
3. **Express the integration contract as constraints.** The non-negotiable rules
   in §2 of the contract — primary keys never null, never duplicated — are three
   lines of SQL each and could run as an assertion over the catalog rather than
   living only as prose.
