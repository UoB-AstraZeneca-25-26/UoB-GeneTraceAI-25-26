# Flag for Track A — sample_info version mismatch affecting depmap_expr resolution

**From:** Musa (Track B)
**Re:** `sample_info` appears to be a stale/older DepMap snapshot than `depmap_profiles` and `signatures`
**Status update:** investigated further since first draft — found a partial fix on my end, narrowed down what's genuinely blocked

---

## Summary

While resolving `depmap_expr` cell line IDs (PR- → ACH- via `depmap_profiles`, confirmed against `sample_info`), 67 of my cell lines failed to match `sample_info.depmap_id`. Investigation across `signatures` and `cellosaurus` confirmed this is a real data-freshness gap in `sample_info`, not a bug in my resolution logic — and let me recover 26 of the 67 myself using an additional fallback path. **41 remain genuinely blocked and need an updated `sample_info` from Track A.**

This affects every track's cell line resolution, since `sample_info` is the canonical registry everyone confirms against — flagging in case others are seeing the same gap without realising it yet.

---

## Evidence

**Max ACH- ID number per file (proxy for how recent each export is):**

| File | Max ACH- ID |
|---|---|
| `sample_info` | ACH-002926 |
| `depmap_profiles` | ACH-003161 |
| `signatures` | ACH-003480 |

`sample_info` tops out well below both other files — it simply does not extend far enough to contain a meaningful chunk of currently-circulating ACH- IDs.

**The 67 originally unresolved IDs split into two patterns:**

### Group A — 25 IDs, tight consecutive block (ACH-003132–ACH-003161)
```
ach-003132, ach-003133, ach-003134, ach-003135, ach-003136,
ach-003138, ach-003139, ach-003141, ach-003142, ach-003143,
ach-003145, ach-003147, ach-003148, ach-003149, ach-003150,
ach-003152, ach-003153, ach-003154, ach-003155, ach-003156,
ach-003157, ach-003158, ach-003159, ach-003160, ach-003161
```
- Exceed `sample_info`'s max ID (002926) entirely — cannot possibly be present in this snapshot.
- Absent from `signatures` too, despite `signatures` extending to 003480 (well past this range).
- **Checked against `cellosaurus` cross-references — 0 of 25 recovered.**
- All three independent sources (`sample_info`, `signatures`, `cellosaurus`) agree these simply don't exist anywhere yet except `depmap_profiles`. This is the newest batch of models — too new to have propagated anywhere else.
- **Reason code:** `unresolvable_too_new_confirmed_3way`

### Group B — 42 IDs, scattered widely (ACH-001437–ACH-002925)
```
ach-001437, ach-001438, ach-001481, ach-001508, ach-001537,
ach-001672, ach-001679, ach-001691, ach-001693, ach-001705,
ach-001847, ach-001854, ach-001855, ach-001971, ach-001975,
ach-001986, ach-001990, ach-002035, ach-002070, ach-002485,
ach-002486, ach-002497, ach-002522, ach-002523, ach-002524,
ach-002526, ach-002531, ach-002533, ach-002535, ach-002539,
ach-002650, ach-002654, ach-002780, ach-002781, ach-002784,
ach-002787, ach-002799, ach-002801, ach-002806, ach-002921,
ach-002922, ach-002925
```
- All confirmed present in `signatures` — genuinely real, valid DepMap models, not a join error on my end.
- Spread across ~1,500 ID numbers in small clusters, rather than one block — suggests `sample_info` has missed individual additions across **multiple** DepMap update cycles, not one release boundary.
- **Checked against `cellosaurus` cross-references (`cross-references` column, parsed `depmap; ach-XXXXXX` pattern) — 26 of 42 recovered**, each with a full cell line name and CVCL accession obtained directly from Cellosaurus (e.g. `MAVER-1 → CVCL_1831`, `SNU-1327 → CVCL_5022`). These 26 are now fully resolved in my crosswalk output, tagged `resolution_method = cellosaurus_xref_fallback`.
- **Checked `cellosaurus.comments` for DepMap "discontinued" flags — 0 of 67 confirmed discontinued.** Ruled out as an explanation; these are not deliberately retired models.
- 16 of the 42 remain unresolved even after the Cellosaurus fallback.
- **Reason code (remaining 16):** `verified_real_missing_from_sample_info`

---

## Net result

| Status | Count |
|---|---|
| Originally orphaned against `sample_info` | 67 |
| Resolved via Cellosaurus cross-reference fallback (name + CVCL obtained) | 26 |
| Still unresolved — too new, confirmed absent from 3 independent sources (Group A) | 25 |
| Still unresolved — real but missing from `sample_info` (Group B remainder) | 16 |
| **Total still blocked, pending updated `sample_info`** | **41** |

---

## What I'm doing about it now

- The 26 recovered via Cellosaurus are added to my crosswalk output directly, with `resolution_method = cellosaurus_xref_fallback` recorded so the data lineage is transparent about how they differ from the primary `sample_info` chain.
- The remaining 41 are logged to `unmapped_depmap_expr.csv` with the two reason codes above, not retried further — per the project's "log, don't force" rule.

## What I need from Track A

1. **Confirm `sample_info`'s export/release date** against `depmap_profiles` and `signatures` — the evidence above strongly suggests it's the oldest of the three.
2. **If possible, source a current `sample_info` export** (DepMap's public portal — depmap.org/portal/download — note the file may now be called `Model.csv` in recent releases rather than `sample_info.csv`). This should resolve most or all of the remaining 41:
   - The 16 from Group B are already independently confirmed real (via `signatures`) — a refreshed file should pick them up immediately.
   - The 25 from Group A are likely the newest models DepMap has added — worth checking if the newest available release includes them, but they may need a more recent release than whatever Track A currently has access to.
3. Worth checking whether this same staleness affects **other tracks' resolution chains** — anything cross-referencing `sample_info.depmap_id` (mutations, fusions, proteomics, metabolomics, instability) likely has the same blind spot to some degree, just not yet quantified by their owners.
4. **Sharing the technique:** the Cellosaurus `cross-references` field (parsed for `depmap; ach-XXXXXX` pattern) recovered 26 of my 42 scattered-gap orphans directly, with full name + CVCL. Worth other tracks checking their own `sample_info` orphans against this same field before assuming a hard blocker — it may recover some of their gaps too without needing a file update at all.

This is no longer fully blocking my `depmap_expr` work — 26/67 are now resolved via the Cellosaurus fallback. The remaining 41 are logged and documented; happy to re-run against an updated `sample_info` whenever it's available, but not waiting idle in the meantime.
