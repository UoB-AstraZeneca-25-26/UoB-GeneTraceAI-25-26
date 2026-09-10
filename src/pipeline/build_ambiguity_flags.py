"""
Per-mechanism gene ambiguity / reliability flags.  (Q8: per-mechanism, not one
column -- the mechanisms have different remedies.)

Currently builds:
  polya_depleted_flag      C2, replication-dependent histones, poly-A depletion
  multi_accession_flag     T3, genes whose protein value is an average over
                           several UniProt accessions before ranking
  ambiguous_accession_flag T2/Q1, genes whose protein accession maps to >1 gene

NOT built here (needs network): paralog_ambiguous_flag -- see C1 /
build_paralog_flags.py. HGNC gene_group is deliberately NOT used as a proxy:
it is functional curation, not sequence identity.

Read-only against the warehouse; writes outputs/gene_ambiguity_flags.parquet.
"""
from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ambiguity_detectors import assert_detectors_agree, composite_sql  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "outputs"
con = duckdb.connect(str(OUT / "celllineselector.db"), read_only=True)
assert_detectors_agree(con)   # C6: python and SQL detectors must match
hg = pd.read_csv(ROOT / "data" / "gene_with_protein_product.txt", sep="\t",
                 dtype=str, low_memory=False)
ens = {s: str(e).lower() for s, e in zip(hg["symbol"], hg["ensembl_gene_id"])
       if isinstance(e, str)}

meta = {}

# ------------------------------------------------------------------ C2
# Replication-dependent (clustered) histones. These are stem-loop terminated,
# not polyadenylated (Marzluff et al. 2008, Nat Rev Genet 9:843-54), and CCLE
# RNA-seq selects polyadenylated mRNA with oligo-dT beads (Ghandi et al. 2019,
# Nature 569:503-508). The depletion is therefore protocol-determined: it is a
# constant offset applied to every cell line, so NO within-gene percentile
# check can detect it -- the percentile is invariant to any transform that
# preserves the within-gene ordering, and this one does.
fam = defaultdict(set)
for sym, grp in zip(hg["symbol"], hg["gene_group"]):
    if isinstance(grp, str):
        for g in grp.split("|"):
            if g.strip():
                fam[g.strip()].add(sym)
hist = set().union(*[s for g, s in fam.items() if "histone" in g.lower()])
rep_dep = {s for s in hist
           if re.fullmatch(r"H[1-4]C\d+[A-Z]?", s)
           or re.fullmatch(r"H2[AB][CU]\d+[A-Z]?", s)}
polya = {ens[s] for s in rep_dep if s in ens}
meta["polya_depleted_flag"] = {
    "n_symbols": len(rep_dep), "n_with_ensg": len(polya),
    "mechanism": "oligo-dT poly-A selection cannot capture stem-loop "
                 "terminated replication-dependent histone transcripts",
    "protocol_citation": "Ghandi et al. 2019, Nature 569:503-508 (CCLE "
                         "TruSeq, non-strand-specific, oligo-dT selection)",
    "biology_citation": "Marzluff et al. 2008, Nat Rev Genet 9:843-54",
    "measured_effect": {"median_log2tpm_histones": 0.6957,
                        "median_log2tpm_coding_background": 2.3305,
                        "delta": -1.6347, "mannwhitney_p": 5.435e-03,
                        "frac_below_background_median": 0.823,
                        "source": "test_run_identifier_ambiguity_audit.py T5"},
    "undetectable_by": "any within-gene percentile check -- the offset is "
                       "identical across all 1,479 cell lines, so it cannot "
                       "change a within-gene rank",
    "quantifier_note": "DepMap TPM comes from RSEM in unstranded mode "
                       "(GENCODE v38, GRCh38); unstranded loses the antisense "
                       "information that would help disambiguate overlaps",
}

# ------------------------------------------------------------------ T3
multi_acc = con.execute("""
    SELECT gene_id, n_uniprot_ids FROM gene_enriched WHERE n_uniprot_ids > 1
""").df()
meta["multi_accession_flag"] = {
    "n_genes": len(multi_acc), "max_accessions": int(multi_acc.n_uniprot_ids.max()),
    "issue": "the protein layer averages the accessions BEFORE ranking "
             "(average-then-rank). rank-then-average would give a different "
             "answer. The choice is undocumented.",
    "aggravating": "the number of accessions contributing varies per cell "
                   "line (PLEC: 3/4/5/6 across 375 lines), so the averaged "
                   "quantity is not homogeneous down the column",
}

# ------------------------------------------------------------------ Q1
# Accessions that HGNC maps to more than one approved symbol. HGNC is a PROXY
# for UniProt here (Q7 remains open) -- so this flag may UNDERCOUNT.
inv = defaultdict(set)
for sym, ids in zip(hg["symbol"], hg["uniprot_ids"]):
    if isinstance(ids, str):
        for a in re.split(r"[|,;]", ids):
            if a.strip():
                inv[a.strip().lower()].add(sym)
ambiguous_acc = {a for a, s in inv.items() if len(s) > 1}

# The HGNC inversion alone is NOT sufficient. q9bte6 (AARSD1 /
# PTGES3L-AARSD1) is single-symbol in HGNC but composite in the ProCan header,
# so it is invisible to the nomenclature route and visible only at source.
# Union both detectors, or the flag silently undercounts by one gene / 707 rows.
composite_acc = set(con.execute(f"""
    SELECT lower(uniprot_id) AS u FROM procan_protein_map
    WHERE {composite_sql("gene_symbol")}
    UNION
    SELECT lower(uniprot_id) FROM protein_map
    WHERE {composite_sql("gene_symbol")}
""").df()["u"])
ambiguous_acc |= composite_acc

live = con.execute("""
    SELECT gene_id, hugo_symbol, lower(u) AS u
    FROM (SELECT gene_id, hugo_symbol, unnest(uniprot_ids) AS u
          FROM gene_enriched WHERE uniprot_ids IS NOT NULL)
""").df()
hit = live[live["u"].isin(ambiguous_acc)]
meta["ambiguous_accession_flag"] = {
    "n_genes": int(hit.gene_id.nunique()),
    "accessions": sorted(hit.u.unique()),
    "genes": sorted(hit.hugo_symbol.dropna().unique()),
    "issue": "the protein measurement is not gene-specific: one accession is "
             "the product of several genes with identical or near-identical "
             "sequence, and the gene it is attributed to was chosen "
             "arbitrarily upstream",
    "detectors": ["HGNC uniprot_ids inverted (accession under >1 approved "
                  "symbol)",
                  "composite symbol field in the raw ProCan / CCLE headers "
                  "(four-delimiter scan)"],
    "n_accessions_hgnc_route": len({a for a, s in inv.items() if len(s) > 1}),
    "n_accessions_composite_route": len(composite_acc),
    "proxy_warning": "detected against the local HGNC release, NOT UniProt. "
                     "Q7 (direct UniProt verification) is open, so this is a "
                     "lower bound. Neither detector subsumes the other: HGNC "
                     "alone misses q9bte6, the header scan alone misses "
                     "p68431/p84243/p04908.",
}

# ------------------------------------------------------------------ assemble
genes = con.execute("SELECT gene_id, hugo_symbol FROM gene_enriched").df()
genes["polya_depleted_flag"] = genes.gene_id.isin(polya)
genes["multi_accession_flag"] = genes.gene_id.isin(set(multi_acc.gene_id))
genes["ambiguous_accession_flag"] = genes.gene_id.isin(set(hit.gene_id))
genes["paralog_ambiguous_flag"] = pd.NA          # C1, network required
genes["any_flag"] = genes[["polya_depleted_flag", "multi_accession_flag",
                           "ambiguous_accession_flag"]].any(axis=1)

p = OUT / "gene_ambiguity_flags.parquet"
genes.to_parquet(p, index=False)
(OUT / "gene_ambiguity_flags_meta.json").write_text(
    json.dumps(meta, indent=1, default=str), encoding="utf-8")

print(f"genes: {len(genes):,}")
for c in ("polya_depleted_flag", "multi_accession_flag",
          "ambiguous_accession_flag", "any_flag"):
    print(f"  {c:<28s} {int(genes[c].sum()):>6,}")
print("  paralog_ambiguous_flag        <NA>  (C1, needs Ensembl paralogy)")
print(f"\n-> {p}")
print(f"-> {OUT / 'gene_ambiguity_flags_meta.json'}")
con.close()
