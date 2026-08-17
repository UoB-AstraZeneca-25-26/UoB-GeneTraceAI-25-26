# GeneTraceAI Pipeline Architecture v2

---

## Exploratory extension: RNA similarity graph clustering (2026-08-17)

The validated RNA similarity graph (1,673 lines, Pearson correlation over the 2,000
most variable DepMap 24Q4 expression features, validated at +0.186 drug-response
concordance over same-tissue baseline; `cell_similarity/08_drug_arm.py`) was
subjected to unsupervised clustering as a descriptive analysis of the graph's
internal structure. Hierarchical clustering (average linkage) and Louvain community
detection on a k-nearest-neighbour graph were both run across a range of parameter
settings; concordance with tissue labels was measured by Adjusted Rand Index and
Normalised Mutual Information against a 1,000-permutation null baseline. The result
is mixed: tissue-of-origin is the dominant organising signal (ARI ≈ 0.22–0.28,
NMI ≈ 0.50–0.55, both methods, all above null at 100%), but the graph also resolves
several cross-tissue communities — a squamous carcinoma community (head-and-neck,
esophagus, urothelial, squamous lung; ~188 lines; 100%/100% bootstrap stability by
both methods) and a renal/Müllerian epithelial community (ovary, kidney, uterus;
120 lines; defined by a PAX8+/HNF1B+ transcriptional program; 98.5% Louvain
bootstrap but 33% HC bootstrap — the HC fragmentation reflects internal
kidney-vs-ovary+uterus structure rather than absence of signal). Bootstrap stability
was verified with 200 line-resampled replicates per method; no duplicate-line artefacts
were found. This exploratory analysis does not feed back into `core_score`,
driver-gating, or the ranking layer; it is purely additive, does not modify any
ranking output, and is documented separately from the defended architectural
decisions in `docs/DEFENCE_CARD.md`. Full results and follow-up sanity checks are in
`cell_similarity/clustering/results/CLUSTERING_RESULTS.md`.
