# MULTIGENE_COVERAGE — the measurements behind A4 and A5

Every figure here is computed by `src/pipeline/measure_multigene_coverage.py`
and stored in `src/pipeline/outputs/multigene_coverage_results.json`. Spec:
[MULTIGENE_SPEC.md](MULTIGENE_SPEC.md).

## The panel

| quantity | value | source |
|---|---|---|
| cell lines | **1,673** | DepMap 24Q4 `OmicsExpressionProteinCodingGenesTPMLogp1.csv` |
| genes | **18,623** | same |
| genes passing the validity guard (`frac_expressed ≥ 0.20`) | **13,047 (70.1%)** | computed |
| lines ever sequenced (alteration layer) | **1,424 (85.1%)** | `mutations_collapsed.parquet` |

The 70.1% reconciles with the 30.5% failure rate in
`docs/TRANSCRIPTOMICS_SPEC.md` F2 — that figure is on a 19,173-gene universe,
this on 18,623.

**249 of 1,673 lines have expression but were never sequenced.** Those lines are
UNKNOWN for every `altered()` term. Under boolean logic they would pass every
`NOT altered(...)` filter ever written.

---

## A4 — joint coverage, measured rather than assumed

### Independence was tested

The brief expected coverage to compound worse than multiplicatively, because
genes co-vary in which lines were profiled. **It does not.** 300 valid genes, 200
random draws per arity, line-level resolution:

| terms | observed | independence predicts | ratio |
|---|---|---|---|
| 1 | 85.2% | 85.2% | 1.00× |
| 2 | **73.7%** | 73.5% | 1.00× |
| 3 | **62.4%** | 62.2% | 1.00× |
| 4 | **59.7%** | 59.3% | 1.01× |

Mixed-type, `expressed(A) AND NOT altered(B)`, 150 random pairs:

| | value |
|---|---|
| `expressed(A)` alone | 87.0% |
| `NOT altered(B)` alone | 85.1% |
| **joint, observed** | **74.2%** |
| independence predicts | 74.1% |

The independence approximation is safe here. That is a measurement, not an
assumption, and it is the reason it can now be used.

### Two denominators, not one

The brief's "~70% for one term, ~48% for two" and the table above are measuring
different things, and both are needed:

| | what it is | value |
|---|---|---|
| **gene-level** | chance a randomly chosen gene is usable at all | **70.1%** |
| **line-level** | share of lines resolved, given the gene passes | **85.2%** |

A query on a random gene faces both. Combining them, a two-term query on two
randomly chosen genes resolves roughly `0.701² × 0.737 ≈ 36%` of the panel —
close to the brief's ~48% intuition but lower, because the gene-level guard bites
twice before the line-level loss applies.

### Worked example: `expressed(EGFR) AND NOT altered(KRAS)`

```
COVERAGE  1,355 of 1,673 lines resolved (81.0%)
  expressed(EGFR)      1,512 resolved (90.4%)   161 unknown -- within measurement uncertainty
  NOT altered(KRAS)    1,424 resolved (85.1%)   249 unknown -- line was never sequenced

MATCHES 792    POSSIBLE 318    EXCLUDED 563
```

The 318 Possible lines split as 249 never-sequenced (KRAS) and the remainder
sitting inside the expression uncertainty band. Both would have passed silently
under boolean `AND NOT`.

---

## A5 — abstention trigger rates

Sampled across the gene universe, 3,000 random genes per condition.

### Single-gene queries

| condition | meaning | rate |
|---|---|---|
| `SILENT_GENE` | `expressed(G)`, `frac_expressed < 0.20` — UNKNOWN everywhere, constrains nothing | **29.6%** |
| `NEAR_UNIVERSAL_GENE` | `NOT expressed(G)`, G expressed in > 95% of lines — excludes almost everything | **43.8%** |
| `TRIVIALLY_SATISFIED` | `NOT altered(G)`, fewer than 10 lines carry a driver alteration | **78.2%** |

**The 78.2% is the most important number in this document for query design.**
Only **4,173 of 18,623 genes** have ≥ 10 lines with a recorded driver alteration.
For the other 78%, `NOT altered(G)` is satisfied by essentially every sequenced
line — the clause looks like a constraint, costs 15% of coverage by pulling in
the never-sequenced lines as UNKNOWN, and narrows nothing. The layer says so
rather than returning a large meaningless match set.

### Two-gene queries

`expressed(A) AND expressed(B)` on two random genes: **at least one abstention
fires in 53.2%** of 1,000 sampled pairs.

More than half of arbitrary two-gene queries are degenerate in at least one term.
This is a property of the data, not of the logic, and it is the reason abstention
is a first-class output rather than a footnote.

---

## Per-lineage resolution

`expressed(EGFR) AND NOT altered(KRAS)`, lineages with ≥ 10 lines, sorted by
resolution. A heavily-unresolved lineage is itself information (A6).

| lineage | n | matches | possible | excluded | resolved |
|---|---|---|---|---|---|
| unknown | 249 | 29 | **156** | 64 | **37.3%** |
| fibroblast | 38 | 25 | 13 | 0 | 65.8% |
| skin | 87 | 32 | 27 | 28 | 69.0% |
| bone | 41 | 15 | 11 | 15 | 73.2% |
| eye | 17 | 6 | 4 | 7 | 76.5% |
| peripheral_nervous_system | 36 | 25 | 7 | 4 | 80.6% |
| … | | | | | |
| liver | 24 | 24 | 0 | 0 | 100.0% |
| pancreas | 53 | 4 | 0 | 49 | 100.0% |
| thyroid | 18 | 17 | 0 | 1 | 100.0% |
| urinary_tract | 37 | 36 | 0 | 1 | 100.0% |

Two readings a user needs:

- **`unknown` at 37.3%** is the worst cluster and it is an artefact of metadata,
  not biology: those 249 lines lack a lineage label *and* are disproportionately
  the never-sequenced ones. Presented as a cluster so it cannot be mistaken for a
  biological group.
- **`pancreas` at 100% resolved with 4 matches and 49 excluded** is a
  fully-informative negative — KRAS is altered in nearly every pancreatic line,
  so the exclusion is real and complete. Contrast with a lineage that has few
  matches *because* it is unresolved. The three-column output distinguishes them;
  a two-bucket output could not.

---

## Timing (A9)

| | value |
|---|---|
| cold load, 507 MB CSV | 13,835 ms |
| cold load, cached parquet | **3,144 ms** |
| warm query, median of 5 | **27 ms** |
| matrix | 1,673 × 18,623 |

Boolean operations over the matrix were never the cost; the source CSV parse was.
Cached to `cell_similarity/outputs/expression_24q4.parquet` on first use.
Real-time confirmed.

---

## Reproduce

```bash
python src/pipeline/measure_multigene_coverage.py
```
