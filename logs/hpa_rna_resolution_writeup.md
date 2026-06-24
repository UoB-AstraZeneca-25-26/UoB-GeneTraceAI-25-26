# hpa_rna — Cell Line Key Resolution

**Owner:** Musa (Track B)
**File:** `hpa_rna` (Human Protein Atlas RNA expression)
**Status:** Cell line key resolution complete. Value transform (Step 5) not yet started.

---

## 1. Why this work was needed

`hpa_rna` carries cell line identity as a free-text, hyphenated name (e.g. `A-431`, `SK-BR-3`) rather than a structured ID like DepMap's ACH- accessions. The existing `clean_hpa_rna()` function only lowercased and whitespace-collapsed this column — it did not normalise punctuation, brackets, or other formatting variants, which meant a meaningful share of names would silently fail to join against `sample_info` without further work. This task built and validated a four-tier resolution chain to get every HPA cell line name to a canonical `model_id` (DepMap ACH-) wherever one genuinely exists, and to clearly identify and document the cases where it does not.

## 2. Files used

| File | Role |
|---|---|
| `hpa_rna` | The file being resolved |
| `hpa_desc` | Tier 3 fallback — provides CVCL accessions for HPA cell lines not resolvable by name alone |
| `sample_info` | Canonical DepMap registry every tier ultimately confirms against (`cell_line_name`, `stripped_cell_line_name`, `rrid`) |
| `cellosaurus` | Tier 4 fallback — name authority used when `sample_info` and `hpa_desc` both fail |

## 3. Resolution chain

```
Tier 1: normalised name → exact match vs sample_info.cell_line_name
Tier 2: normalised name → match vs sample_info.stripped_cell_line_name
Tier 3: normalised name → CVCL via hpa_desc → match vs sample_info.rrid
Tier 4: normalised name → CVCL via cellosaurus (disambiguated, see below)
        → match vs sample_info.rrid
        → if CVCL found but no sample_info match: identified, no DepMap equivalent
Unresolved → logged with name preserved
```

## 4. Problems found and fixed

### 4.1 — Punctuation normalisation gap
`clean_hpa_rna()` only lowercased cell line names; hyphens, slashes, and underscores were left untouched. Names like `23132/87` survived cleaning with punctuation intact and would not exact-match `sample_info`. Fixed by writing a dedicated `normalise_base()` function that strips hyphens, slashes, underscores, and whitespace before any comparison. This single fix moved a substantial number of cell lines from needing the expensive Tier 3/4 fallback to resolving immediately at Tier 1.

### 4.2 — Bracketed annotations not handled
A small number of names carry a descriptive annotation in brackets, e.g. `BJ [human fibroblast]`. The base normalisation left brackets and their contents in place, so these names matched nothing. Fixed by stripping bracket content for the base-name comparison.

### 4.3 — Bracket-stripping introduced a silent collision risk (significant finding)
Stripping bracket content blindly created a new problem: Cellosaurus contains multiple, biologically distinct cell lines that share the same base name and are only distinguished by their bracket content — for example `BJ [human fibroblast]`, `BJ [human B-cell IHW]`, and `BJ [human pancreatic adenocarcinoma]` are three unrelated cell lines that all collapse to the same key `"bj"` once brackets are stripped. An audit of the full Cellosaurus table found **730 base names** with this same multi-entry collision risk.

Naively taking the first match (`.iloc[0]`) in these cases silently returns *a* result without confirming it is the *correct* result. This was caught when resolving `BJ [human fibroblast]` returned the wrong CVCL — a different cell line entirely — leading to an apparently successful but incorrect match.

**Fix:** added bracket-content-aware disambiguation. When a base name has multiple Cellosaurus candidates, the original bracket content from the HPA name is used to select the correct one. If disambiguation still cannot isolate a single match, the row is logged rather than guessed.

This fix changed the result set meaningfully: several names that had previously appeared to resolve successfully via the naive Cellosaurus match were revealed, once disambiguated correctly, to have no DepMap equivalent at all. The corrected result is smaller in raw resolved count than the naive first pass, but is verified and trustworthy rather than silently wrong.

### 4.4 — `+` character not stripped by normalisation
Two remaining unresolved names (`BJ hTERT+ SV40 Large T+`, `BJ hTERT+ SV40 Large T+ RasG12V`) contained a `+` character that the normalisation function did not strip. Added `+` to the stripped-character set and re-checked both names against Cellosaurus and `sample_info` directly — confirmed no match exists under any reasonable normalisation, ruling this out as the cause of their non-resolution.

## 5. Final result

| Category | Count | % | Resolution method |
|---|---|---|---|
| Resolved — usable for confidence scoring | 1,104 | 91.5% | `sample_info_exact` (1,024), `sample_info_stripped` (37), `hpa_desc_cvcl_fallback` (43) |
| Identified, no DepMap equivalent | 97 | 8.0% | `cellosaurus_name_fallback` — real cell line, valid CVCL, but no DepMap ACH- counterpart |
| Genuinely unresolved | 5 | 0.4% | No match found via any tier after exhausting normalisation options |
| **Total** | **1,206** | **100%** | |

### The 5 unresolved cases, individually investigated

| Name | Likely explanation |
|---|---|
| `asc2telo differentiated` | Engineered/differentiated derivative of a base line, not separately catalogued |
| `bj htert+ sv40 large t+` | Engineered subline of BJ-hTERT (already resolved separately); checked with `+` stripped, no match in Cellosaurus or `sample_info` |
| `bj htert+ sv40 large t+ rasg12v` | Further-engineered subline (adds RAS-G12V); same check, no match |
| `hhstec` | Likely a primary cell type abbreviation, not an established/catalogued line |
| `hskmc` | Likely a primary cell type abbreviation, not an established/catalogued line |

Each was checked against both Cellosaurus and `sample_info` with full normalisation (bracket-content disambiguation, `+` stripping) before being classified as genuinely unresolved — not assumed.

## 6. Outputs produced

```
hpa_rna_cell_line_crosswalk.csv       — 1,104 rows, has model_id, ready for integration
hpa_rna_no_depmap_equivalent.csv      — 97 rows, identified with CVCL, no model_id
unmapped_hpa_rna.csv                  — 5 rows, genuinely unresolved, with name preserved
```

## 7. Methodological note worth sharing with the team

The bracket-collision bug (Section 4.3) is a general risk for anyone doing name-based fuzzy matching against Cellosaurus elsewhere in the project, not just for `hpa_rna`. A naive name-stripping match can return a confident, plausible-looking, but incorrect result without raising any error — this is more dangerous than a clean failure, because it can silently introduce a wrong cell line identity into the dataset. Recommend any other track using Cellosaurus name matching as a fallback check their own matches for the same collision pattern before trusting results.

## 8. Next step

Step 5 — apply the value transform to `hpa_rna`: drop `tpm` and `ptpm`, keep only `ntpm`, apply `log2(ntpm+1)`. This is independent of the identity work above and is being deferred until the cell line crosswalk was fully locked in, per the project's rule of resolving identity before transforming values.
