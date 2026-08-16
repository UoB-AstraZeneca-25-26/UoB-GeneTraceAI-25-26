# DEFENCE_CARD — Sex-Linked Gene Guard

**Date:** 2026-08-16  
**Scope:** `final_pipeline/Ranking/cli.py` — `cmd_exclude` and `cmd_genes`  
**Fault class:** Absence-as-negative, instance 8 (see `MULTIGENE_SPEC.md` §A1)

---

## Four-move structure

### Move 1 — What it does

Y-linked genes (ENSG IDs with `chromosomal_location` starting "Y" in
`gene_lookup.parquet`) are structurally absent from XX (female / unknown-sex)
cell lines. Before this fix, the scoring pipeline assigned these lines a low
`core_score` for Y-linked genes — correctly, in the sense that the essentiality
screen found no dependency. But the ranking CLI consumed that low score as if
it were a biological signal ("this gene is non-essential in this line") rather
than a structural absence signal ("this gene does not exist in this line").

In `cmd_exclude(GENE_A, GENE_B)`, a low `score_b` for gene_B inflates
selectivity (`score_a × (1 − score_b)`), so female lines falsely appeared as
ideal cell lines for a GENE_A × GENE_B exclusion experiment where gene_B is
Y-linked. Example: 17 of the top-30 selectivity pairs for DDX3Y as GENE_B were
female cell lines. After the fix: 0.

In `cmd_genes`, a low Y-linked gene score depresses the `joint_score` (min of
all genes), producing false negatives — male lines with genuine co-dependency
were ranked below female lines where the score was structurally zero.

**The fix:** `score_b = NaN` (exclude) and per-gene `score = NaN` (co-select)
for female / unknown-sex lines where the query gene is Y-linked.
`selectivity = NaN` sorts to last position — excluded from display, not set to
zero or to a penalised rank. Consistent with the existing orphan-pair NaN
pattern in `driver_routing.py:100`.

---

### Move 2 — What was compared

Three alternative guard designs were considered:

| Option | Description | Rejected because |
|---|---|---|
| **NaN exclusion** ✓ | `score_b = NaN` → `selectivity = NaN` → sorted last | Consistent with existing orphan-pair handling; explicit in output |
| Penalised rank | Push to rank N+1 explicitly | Invents a rank that doesn't exist; user sees a number that implies measurement |
| Reason-code column | Add `exclusion_reason='y_linked_sex_guard'` alongside score | Richer but would require CLI output schema change; printout already shows `[SEX GUARD]` banner |
| 0.0 fill | Set `score_b = 0.0` | Conflates "structural absence" with "measured essential in every line" — worse than NaN |

---

### Move 3 — What decided it

The NaN approach was chosen because:

1. It is the existing project convention for "not evaluable" pairs (orphan
   pairs in `driver_routing.py`). Using the same representation keeps the
   downstream guarantee: `NaN selectivity → excluded`, no special-casing needed.

2. The CLI already prints a visible `[SEX GUARD]` banner before the ranked
   table, reporting the count of excluded pairs and the reason. The researcher
   cannot miss the exclusion.

3. Sex annotation in `sample_info`: male 991, female 747, unknown 102 (5.5%).
   Unknown-sex lines are conservatively included in the guard — they may be
   female, and the cost of a false exclusion (slightly smaller ranked set) is
   lower than the cost of a false positive (a female line appearing as top
   selectivity for a Y-linked exclusion gene).

4. Y-linked gene identity comes from `gene_lookup.parquet`'s
   `chromosomal_location` column (HGNC-sourced, format "Yq11.221" / "Yp11.2"),
   via `str.startswith("Y")`. Not a hardcoded list — so it won't silently go
   stale when new Y-linked genes enter the scored universe.

---

### Move 4 — What it costs

- **Runtime:** one additional DuckDB query per CLI invocation (`sample_info`
  sex lookup) and one parquet read that already happened (`gene_lookup`). Both
  are sub-second. No impact on throughput.

- **Coverage reduction:** 102 unknown-sex lines are conservatively excluded when
  a Y-linked gene_B is queried. These lines remain visible in `cmd_gene` (single
  gene view); the guard applies only to the pairwise commands where structural
  absence would produce false positives/negatives.

- **Stage 6 validation:** 0 Y-linked genes appear in `gdsc_scored_ready.parquet`
  (confirmed). Hit@20 figures in `ALTERATION_DEFENCE.md` are unchanged.

---

## Limitations

### X-linked dosage-sensitive genes — NOT auto-guarded

X-linked genes with dosage-sensitive functions (e.g. KDM6A, KDM5C, ATRX,
MED12, BRWD3, HUWE1 — those on Xp/Xq arms that escape X-inactivation) are
structurally present in female cells but expressed at different dosage due to
incomplete X-inactivation (XCI) or XCI escape.

**Why not auto-guarded:**

- X-inactivation is gene-specific and variable. Blanket sex-guarding of
  X-linked genes would suppress real signal from lines where XCI escape allows
  normal dosage. Unlike Y-linkage (binary: present / absent), XCI is a
  continuous and gene-specific quantity.

- Cancer cell lines show widespread skewed XCI that is not annotated in DepMap's
  `sample_info`. Guarding without per-gene XCI measurements would introduce noise
  proportional to how wrong the escape-status assumption is.

- The correct treatment requires per-gene XCI escape classifications (Cotton
  et al. 2013, PMID 23352436; Tukiainen et al. 2017, PMID 29022598), not
  available in the current pipeline inputs.

**What to do instead:** when querying an X-linked dosage-sensitive gene, note
the effect is present in both sexes at different magnitude. The researcher
should stratify by sex in the output and interpret female/male sub-rankings
separately. The `cmd_gene` output shows sex annotation in the detail view
(`--detail ACH-ID`); the `cmd_exclude` and `cmd_genes` outputs show sex only
when a guard fires.

**Identified X-linked dosage-sensitive candidates** (escape confirmed in Cotton
et al. / Tukiainen et al. and present in scored universe — verify before
guarding): KDM6A, KDM5C, ATRX, MED12, BRWD3, HUWE1, KDM5D-homologue region.
This list is not exhaustive. Decision deferred to researcher.
