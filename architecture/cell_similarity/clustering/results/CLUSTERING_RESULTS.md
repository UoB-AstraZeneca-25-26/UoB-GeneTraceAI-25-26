# Clustering Analysis of the RNA Similarity Graph
**Date:** 2026-08-17 (v8.0); updated 2026-08-17 (v8.1 follow-up: F1 bootstrap, F2 duplicate check, F3 cluster #2 vet)
**Status:** Exploratory result — does not feed into `core_score` or the ranking pipeline

---

## What was tested

Unsupervised clustering of the validated RNA similarity matrix (1,673 cell lines ×
1,673 cell lines, Pearson correlation over the 2,000 most variable DepMap 24Q4
expression features, z-scored per gene). The question: do clusters correspond to
tissue of origin, or do they cut across tissue boundaries in a way that is
biologically coherent?

**Methods:**
- Hierarchical clustering (average linkage) on distance = 1 − similarity, cutting
  at six cluster counts: 24, 29, 34, 39, 44, 51 (tissue count n=34 ± ~30%)
- Louvain community detection on a k-nearest-neighbour graph at k = 10, 15, 20, 30
  (10 clusterings total)
- Concordance metrics: Adjusted Rand Index (ARI) and Normalised Mutual Information
  (NMI) against DepMap lineage labels (same source as Stage 1 of the pipeline)
- Null baseline: 1,000 label-shuffle permutations per clustering; observed ARI/NMI
  compared against this distribution

**Tissue distribution:** 34 distinct labels; 8 tissues with < 10 lines
(unknown=9, testis=7, muscle=5, embryo=4, hair=2, adrenal\_cortex=2, pleura=1,
epidermoid\_carcinoma=1) — noted, not dropped; included in metrics but too small
for per-tissue interpretation.

---

## What was found

### Concordance with tissue labels

All 10 clusterings are above the null at 100% (observed ARI and NMI > every
permuted baseline across 1,000 shuffles). The signal is real. But the magnitude
tells the story:

| Method | k | Clusters | ARI | NMI |
|---|---|---|---|---|
| HC (average) | 24 | 24 | 0.219 | 0.497 |
| HC (average) | 29 | 29 | 0.219 | 0.498 |
| HC (average) | 34 | 34 | 0.219 | 0.504 |
| HC (average) | 39 | 39 | 0.223 | 0.509 |
| HC (average) | 44 | 44 | 0.271 | 0.524 |
| HC (average) | 51 | 51 | 0.282 | 0.535 |
| Louvain | 10 | 18 | 0.265 | 0.551 |
| Louvain | 15 | 15 | 0.240 | 0.524 |
| Louvain | 20 | 14 | 0.255 | 0.522 |
| Louvain | 30 | 11 | 0.247 | 0.504 |
| **Null mean** | — | — | **~0.000** | **~0.04–0.10** |

**ARI = 0.22–0.28** across both methods. ARI = 1.0 is perfect tissue recovery;
ARI = 0.0 is chance. The graph is tissue-informed — tissue-of-origin is the
dominant single factor — but clusters are not tissue-defined: the RNA similarity
graph organises lines in ways that substantially diverge from the tissue partition.

**NMI ≈ 0.50–0.55** means clusters share roughly half the information content with
tissue labels. The other half is structure the tissue label does not capture.

Both methods agree on the range. This agreement across independent algorithms is
the primary claim.

### Clusters that are nearly pure (≥90% one tissue)

At k=34 (hierarchical), 9 of 34 clusters are ≥90% single-tissue. At k=51, 12 of
51 are nearly pure. These confirm that the graph strongly recovers the most
molecularly distinctive lineages (skin/melanoma, fibroblast, blood, colorectal,
kidney at finer cuts).

---

## Bootstrap stability (F1 follow-up)

**Criterion (pre-stated):** A reference cluster is considered "reproduced" in a
bootstrap replicate if there exists a replicate cluster C′ such that:

> |R\_obs ∩ C′| / |R\_obs| ≥ 0.50

where R\_obs = reference cluster lines that appear in this bootstrap sample. Applied
across 200 replicates of line-level bootstrapping (sample with replacement) ×
2 methods (Louvain k=20, HC average k=34). Replicates with < 5 observable members
are excluded from the count.

| Cluster | Louvain k=20 (200 replicates) | HC k=34 (200 replicates) |
|---|---|---|
| Hematological (n=65, blood) | **100.0%** (200/200) | **100.0%** (200/200) |
| Gynecological epithelial (n=120) | **98.5%** (197/200) | **33.0%** (66/200) |
| Squamous carcinoma (n=188) | **100.0%** (200/200) | **100.0%** (200/200) |
| Mesenchymal/neuroectodermal (n=385) | **100.0%** (200/200) | **100.0%** (200/200) |

**Hematological note:** In Louvain\_k20, haematopoietic lines are split across three
sub-communities by exact subtype (clusters 3, 5, 7: pure blood / blood+lymphocyte /
lymphocyte+plasma\_cell). The F1 test used the largest (n=65, pure blood). The
blood+lymphocyte+plasma\_cell supercluster visible in HC at coarse cuts is a
hierarchical artefact of merging these three; it is not a finding to report as novel.

**Gynecological cluster — ESCALATE condition and resolution:**
The HC reproduction rate is 33%, triggering the escalation condition in v8.1 §2.
Diagnosis: at k=34, HC splits the 120-line community into two adjacent HC branches:
- HC cluster 3 (n=59 in full panel): kidney-heavy (40 kidney, 9 ovary, 3 uterus, 3
  urinary\_tract) — captures 44% of gynaecological lines
- HC cluster 9 (n=56 in full panel): ovary+uterus-heavy (32 ovary, 22 uterus) —
  captures 43% of gynaecological lines

Neither branch alone reaches the 50% threshold. In bootstrap replicates, the
split is not reproducible at a fixed k=34 — some replicates one branch absorbs
more than 50%, others don't. The community structure is real (Louvain 98.5%
bootstrap), but the HC split diagnoses a genuine internal structure: kidney lines
and ovary+uterus lines are more strongly connected within their respective
sub-groups than between them, while still being drawn into the same Louvain
community by their cross-tissue edges. This is correctly described as a **graph
community spanning two adjacent HC branches**, not as a single undifferentiated
cluster. The characterisation has been updated accordingly (see below).

---

## Duplicate / technical replicate check (F2 follow-up)

DepMap Model.csv was checked for duplicate ModelIDs, quality/contamination flags,
and any problematic-line signals across all four reference clusters. No
contamination or duplication columns are present in 24Q4 Model.csv (Cellosaurus
cross-check is not automatable without web access and is noted as a gap). No
duplicate ModelIDs were found within any cluster. All four clusters are **clean**
by the checks that can be performed from on-disk metadata.

---

## Recurring cross-tissue clusters

### 1. Hematological sub-communities

As noted above, this is a taxonomy artefact in the HC output: HC at coarse cuts
merges blood, lymphocyte, and plasma\_cell (separately labelled but RNA-similar)
into one supercluster. In Louvain, they are three distinct sub-communities. Not
reported as a novel cross-tissue biological finding.

---

### 2. Gynecological / renal epithelial community *(moderate–high confidence)*

**Exact composition (Louvain\_k20, cluster 9):** n = 120 lines

| Tissue | Count | % |
|---|---|---|
| ovary | 39 | 32.5% |
| kidney | 36 | 30.0% |
| uterus | 23 | 19.2% |
| liver | 7 | 5.8% |
| bile\_duct | 4 | 3.3% |
| urinary\_tract | 3 | 2.5% |
| lung | 3 | 2.5% |
| thyroid | 2 | 1.7% |
| breast + pancreas + soft\_tissue | 3 | 2.5% |

Core three tissues (ovary+kidney+uterus): **98 lines (81.7%)**.
Within the core three, no single tissue dominates: ovary 39.8%, kidney 36.7%,
uterus 23.5%. Three-way mix confirmed.

**Bootstrap stability:** Louvain 98.5% (197/200); HC 33% (66/200). The HC rate
reflects the internal structure of the community (adjacent HC branches for kidney
vs ovary+uterus) rather than absence of signal. Characterisation: **graph community
at moderate–high confidence**.

**Cluster-defining genes (F3.3 — top genes by Cohen's d, cluster vs rest of panel):**

| Gene | Cohen's d | Known biology relevant to this cluster |
|---|---|---|
| PAX8 | +3.38 | Marks carcinomas of **nephric duct** (kidney, ureter) and **Müllerian duct** (uterus, cervix, fallopian tube) origin; used as a pathological marker for renal cell carcinoma and ovarian/endometrial carcinoma |
| HAVCR1 (KIM-1) | +3.21 | Kidney Injury Molecule 1; marker of proximal tubular epithelium and clear cell RCC |
| KCNJ16 | +3.27 | Renal potassium channel; kidney tubular expression |
| CDH16 | +2.83 | KSP-cadherin; expressed almost exclusively in kidney tubular epithelium |
| CLDN16 | +2.81 | Claudin-16; tight junction protein specific to the thick ascending limb of Henle's loop |
| HNF1B | +2.54 | Transcription factor required for kidney and Müllerian duct development; **specifically overexpressed in clear cell carcinoma of both ovarian and renal types** |
| FXYD2 | +2.43 | Renal tubular Na/K-ATPase regulator |
| SLC17A3 / SLC17A1 | +2.33 / +2.27 | Sodium-phosphate transporters; kidney proximal tubule |

The top-upregulated genes in this cluster form a coherent **renal tubular /
Müllerian epithelial transcriptional program**. PAX8 is the molecular unifier: it
marks carcinomas derived from both the nephric duct (kidney) and the Müllerian duct
(ovary, uterus, fallopian tube). HNF1B specifically bridges the kidney and ovarian
clear cell connection.

**Literature context (F3.4):**
- PAX8 is well established as a pathological marker for renal cell carcinoma (RCC)
  and ovarian/uterine carcinoma. Its expression in kidney and Müllerian structures
  reflects a shared developmental ancestry: both the nephric and Müllerian ducts
  arise from the **urogenital ridge**, providing a developmental origin for the
  transcriptional program this cluster appears to reflect.
- HNF1B is specifically overexpressed in clear cell carcinoma of the ovary and is
  used as an IHC marker for this subtype; clear cell RCC also has elevated HNF1B
  relative to other RCC subtypes.
- This grouping is consistent with independently-established cell-of-origin
  structure in TCGA's pan-cancer analyses, which separately identify pan-kidney and
  pan-gynecological cancers as coherent molecular groupings driven by shared
  developmental lineage (Hoadley et al., 2018). That analysis does not itself report
  a single merged kidney/ovary/uterus cluster — it treats pan-kidney and
  pan-gynecological as distinct groupings — so this should be read as convergent
  support for the PAX8/HNF1B-driven grouping observed here, not as direct
  replication of it. The two analyses also differ methodologically: Hoadley et al.
  used integrative multi-platform clustering (aneuploidy, methylation, mRNA, miRNA,
  and protein arrays) on ~10,000 tumour samples, whereas this analysis uses RNA-only
  similarity clustering on ~1,673 cell lines.
- **The puzzle:** kidney does NOT share Müllerian developmental origin with ovary and
  uterus. The connection appears to be transcriptional convergence driven by shared
  PAX8 / HNF1B expression programs, not shared developmental lineage. Whether the
  full kidney lineage is driving this or specifically clear-cell-histology lines is
  unresolvable from lineage labels alone — histology metadata is required for that
  question.

**What is claimed:** these 120 lines are RNA-similar across tissue boundaries, and
the similarity is dominated by a PAX8+ / HNF1B+ renal/Müllerian epithelial
transcriptional signature. No further mechanistic claim is made.

**What is NOT claimed:** that this represents a shared molecular subtype in the
oncological sense; that all kidney lines or all ovarian lines behave this way; that
the HC fragmentation is irrelevant (it is informative — the community has internal
kidney vs ovary+uterus structure). This cluster should be discussed as a candidate
finding, not a definitive one, given the HC stability caveat.

---

### 3. Squamous carcinoma community *(high confidence)*

**Bootstrap stability:** Louvain 100% (200/200); HC 100% (200/200).

~188 lines across all Louvain settings (upper\_aerodigestive 34%, esophagus 12%,
lung 11%, urinary\_tract 11%, cervix 6%). Head and neck squamous, esophageal
squamous, squamous-subtype lung, urothelial, and cervical carcinoma lines group
together. These share a squamous differentiation program and are known to form a
Pan-Squamous molecular cluster in TCGA (TCGA Research Network, Cancer Cell 2018).
The RNA similarity graph recovers this known grouping without being told it. This
is a well-established molecular grouping; reporting it as a finding here provides
independent confirmation from a different dataset (DepMap expression) and method.

---

### 4. Mesenchymal / neuroectodermal supercluster *(high confidence, heterogeneous)*

**Bootstrap stability:** Louvain 100% (200/200); HC 100% (200/200).

~289–385 lines depending on method (central\_nervous\_system 22–27%, fibroblast
10–13%, soft\_tissue 11%, lung 9–14%, bone 6–8%). CNS, soft tissue sarcoma, bone
sarcoma, and fibroblast lines cluster together — sharing neuroectodermal or
mesenchymal developmental origin. The lung lines within this cluster likely
represent a histological subset (neuroendocrine or sarcomatoid), but this cannot
be confirmed from lineage labels alone. Noted without mechanistic interpretation.
The lung representation (9–14%) and the heterogeneity of this cluster mean it is
the least diagnostically specific of the four, despite its high bootstrap
stability.

---

## Sanity checks summary

| Check | Hematological | Gynecological | Squamous | Mesenchymal |
|---|---|---|---|---|
| Y-linked gene drive | None | None | None | None |
| Duplicate ModelID | None | None | None | None |
| Quality flags in Model.csv | None | None | None | None |
| Bootstrap (Louvain) | 100% | 98.5% | 100% | 100% |
| Bootstrap (HC) | 100% | **33%** — see note | 100% | 100% |
| HC fragmentation diagnosis | n/a | Adjacent branches, methodological | n/a | n/a |

---

## Context: what this means for the pipeline

The clustering result is consistent with the similarity graph's validated drug-
response concordance (+0.186 over same-tissue baseline). The same-tissue baseline
corresponds to the tissue-recapitulating structure (ARI ~0.25); the benefit over
that baseline corresponds to the cross-tissue structure the clustering characterises.

This does not change any pipeline output. Core scores, driver-gating, ranking, and
the sex guard are unaffected. The clustering is a descriptive characterisation of
the validated similarity graph.

---

## What this is NOT

- This analysis does not claim shared molecular subtypes without independent
  evidence. Gene lists are reported as observed; mechanistic interpretation is
  reserved for claims that are either well-established (squamous cluster) or flagged
  as hypotheses (gynecological epithelial).
- This is not defended as an architectural decision. It is exploratory output.
- The HC instability of the gynecological cluster (33%) is not explained away — it
  is reported, diagnosed, and incorporated into the characterisation. Presenting
  only the Louvain result would misrepresent the evidence.
- The duplicate-line check (v8.1 F2) covered exact ModelID duplication and known
  quality-flag columns in DepMap 24Q4 Model.csv, but did not cross-reference
  Cellosaurus for documented derivative or related lines (different passage, same
  patient, or documented cross-contamination history). This remains an open
  verification step. If any of the 120 gynecological cluster lines are
  Cellosaurus-flagged, the composition may shift.

---

## Artefact locations

| File | Contents |
|---|---|
| `results/concordance_scores.csv` | ARI/NMI per clustering with null baseline stats |
| `results/cluster_assignments.csv` | Per-line cluster label, all methods |
| `results/cluster_compositions.csv` | Per-cluster tissue breakdown, all methods |
| `results/full_results.json` | Complete C1–C4 structured output |
| `results/followup_results.json` | F1 bootstrap rates, F2 duplicate check, F3 cluster #2 detail |
| `results/cluster2_defining_genes.csv` | Top 20 up/down genes by Cohen's d (cluster #2 vs panel) |
| `cluster_analysis.py` | Reproducible v8.0 main analysis |
| `followup_analysis.py` | Reproducible v8.1 follow-up analysis |

---

## References

Hoadley, K.A., Yau, C., Hinoue, T., Wolf, D.M., Lazar, A.J., Drill, E., Shen, R.,
Taylor, A.M., Cherniack, A.D., Thorsson, V., Akbani, R., Bowlby, R., Wong, C.K.,
Wiznerowicz, M., Sanchez-Vega, F., Robertson, A.G., Schneider, B.G., Lawrence,
M.S., Noushmehr, H., Malta, T.M., Stuart, J.M., Benz, C.C. and Laird, P.W. (2018)
'Cell-of-origin patterns dominate the molecular classification of 10,000 tumors from
33 types of cancer', *Cell*, 173(2), pp. 291–304.e6.
doi: 10.1016/j.cell.2018.03.022.
