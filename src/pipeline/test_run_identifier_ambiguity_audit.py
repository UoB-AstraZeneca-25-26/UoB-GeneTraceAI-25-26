"""
Identifier-ambiguity audit (T1-T8), read-only.

Motivation
----------
The ProCan protein matrix carries 18 columns whose `symbol` header names 2-14
genes for a single UniProt accession (P62805 -> 14 histone H4 genes). Stage 0
stores that composite verbatim in `procan_protein_map.gene_symbol`, so every
uniqueness check on the MAP tables passes: 'h4c1;h4c2;...' is one distinct
string. The ambiguity is only visible in the raw header. This harness
generalises that miss into eight tests over every raw source and its ingested
counterpart.

Test taxonomy (the failure classes the tests are looking for)
  Class A  composite identifier -- one row/column names several entities
  Class B  arbitrary resolution -- upstream silently picked one winner
  Class C  systematic measurement ambiguity -- paralogy, protocol depletion

Tests
  T1  composite identifier scan over raw identifier columns
  T2  raw-vs-warehouse cardinality reconciliation (+ Q1: accessions that map
      to >1 gene, measured against the local HGNC release as a UniProt proxy)
  T3  reverse direction: one gene, many features
  T4  paralog-family annotation from HGNC gene_group
  T5  replication-dependent histone / poly-A depletion, measured
  T6  identifier namespace collisions across tables
  T7  multi-target reagent provenance (Chronos release, GDSC drug ids)
  T8  frozen regression assertions

Nothing here writes to the warehouse. Results land in
outputs/test_run_identifier_ambiguity_audit_results.json.
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT = Path(__file__).resolve().parent / "outputs"
DB = OUT / "celllineselector.db"

CCLE_PROT = DATA / "parquet" / "4_Harmonized_MS_CCLE_Gygi_subsetted.parquet"
PROCAN = DATA / "proteomics_procan" / "Protein_matrix_averaged_20250211.tsv"
HGNC_RAW = DATA / "gene_with_protein_product.txt"
CGC = DATA / "COSMIC" / "Cosmic_CancerGeneCensus_v104_GRCh38.tsv"
COSMIC_GENES = DATA / "COSMIC" / "Cosmic_Genes_v104_GRCh38.tsv"
GDSC_CMPD = DATA / "GDSC" / "screened_compounds_rel_8.5.csv"
GDSC_MODELS = DATA / "GDSC" / "model_list_20260709.csv"
CHRONOS = DATA / "DepMap_Chronos"

# C6: the delimiter set is defined ONCE, in ambiguity_detectors, and the SQL
# form is derived from the same constant. Restating it here is what produced
# the original false pass (LIKE '%;%' missing the CCLE slash composite).
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ambiguity_detectors import (DELIMS, assert_detectors_agree,   # noqa: E402
                                 composite_sql, n_entities)

R = {}          # results
con = duckdb.connect(str(DB), read_only=True)


def head(t, n=72):
    print()
    print("=" * n)
    print(t)
    print("=" * n)


def sub(t):
    print()
    print(f"--- {t}")


# ---------------------------------------------------------------- helpers
def procan_header():
    """(accessions, symbols) from the three-row ProCan header block."""
    with open(PROCAN, encoding="utf-8") as f:
        l0 = next(f).rstrip("\n").split("\t")
        l1 = next(f).rstrip("\n").split("\t")
    return l0[2:], l1[2:]


def ccle_header():
    """[(accession, symbol_or_None)] from the Gygi parquet column names."""
    cols = [c[0] for c in con.execute(
        f"DESCRIBE SELECT * FROM read_parquet('{CCLE_PROT.as_posix()}')").fetchall()]
    out = []
    for c in cols[1:]:
        m = re.fullmatch(r"(\S+) \((.*)\)", c)
        out.append((m.group(1), m.group(2)) if m else (c, None))
    return out


def parquet_cols(p):
    return [c[0] for c in con.execute(
        f"DESCRIBE SELECT * FROM read_parquet('{Path(p).as_posix()}')").fetchall()]


def tsv_header(p, sep="\t"):
    with open(p, encoding="utf-8", errors="replace") as f:
        return next(csv.reader(f, delimiter=sep))


def scan_composites(pairs, label):
    """pairs: iterable of (identifier, value_that_might_be_composite)."""
    hits = []
    for ident, val in pairs:
        if val is None:
            continue
        s = str(val)
        if DELIMS.search(s):
            hits.append({"identifier": ident, "value": s,
                         "n_entities": n_entities(s)})
    hits.sort(key=lambda h: -h["n_entities"])
    print(f"  {label:<46s} composites = {len(hits):>5,}")
    return hits


# ================================================================ T1
head("T1 -- Composite identifier scan (Class A), raw sources pre-ingest")

t1 = {}

sub("proteomics: symbol fields in the raw headers")
p_acc, p_sym = procan_header()
t1["procan_symbol_row"] = scan_composites(zip(p_acc, p_sym), "ProCan header row 1 (symbol)")
c_hdr = ccle_header()
t1["ccle_symbol_header"] = scan_composites(
    ((a, s) for a, s in c_hdr), "CCLE/Gygi header (symbol in parens)")

sub("HGNC raw: the nomenclature fields that are legitimately multi-valued")
hg = pd.read_csv(HGNC_RAW, sep="\t", dtype=str, low_memory=False)
for col in ("uniprot_ids", "alias_symbol", "prev_symbol", "ensembl_gene_id",
            "entrez_id", "refseq_accession", "gene_group"):
    if col in hg.columns:
        hits = scan_composites(zip(hg["symbol"], hg[col]), f"hgnc.{col}")
        t1[f"hgnc_{col}"] = hits[:40]

sub("COSMIC / GDSC / Cellosaurus identifier columns")
if CGC.exists():
    cgc = pd.read_csv(CGC, sep="\t", dtype=str, low_memory=False)
    idcol = "GENE_SYMBOL" if "GENE_SYMBOL" in cgc.columns else cgc.columns[0]
    for col in [c for c in ("GENE_SYMBOL", "SYNONYMS", "TRANSCRIPT_ACCESSION",
                            "TIER") if c in cgc.columns]:
        t1[f"cgc_{col}"] = scan_composites(zip(cgc[idcol], cgc[col]), f"cgc.{col}")[:40]

cmpd = pd.read_csv(GDSC_CMPD, dtype=str)
for col in ("DRUG_NAME", "SYNONYMS", "TARGET"):
    t1[f"gdsc_{col}"] = scan_composites(zip(cmpd["DRUG_ID"], cmpd[col]),
                                        f"gdsc_compounds.{col}")[:40]

cl = con.execute("SELECT cellosaurus_accession a, cellosaurus_cell_line_name n, "
                 "synonyms s FROM cellosaurus").df()
t1["cellosaurus_name"] = scan_composites(zip(cl.a, cl.n), "cellosaurus.name")[:40]

sub("warehouse map tables (where the original miss was invisible)")
pm = con.execute("SELECT uniprot_id, gene_symbol FROM protein_map").df()
t1["protein_map"] = scan_composites(zip(pm.uniprot_id, pm.gene_symbol), "protein_map.gene_symbol")
ppm = con.execute("SELECT uniprot_id, gene_symbol FROM procan_protein_map").df()
t1["procan_protein_map"] = scan_composites(zip(ppm.uniprot_id, ppm.gene_symbol),
                                           "procan_protein_map.gene_symbol")
R["T1"] = {k: {"n": len(v), "examples": v[:20]} for k, v in t1.items()}


# ================================================================ T2
head("T2 -- Raw-vs-warehouse cardinality reconciliation (Class B)")

# gene_enriched.uniprot_ids is the ONLY route from a gene to a protein column.
linked = con.execute("""
    SELECT lower(u) AS u, gene_id, hugo_symbol
    FROM (SELECT gene_id, hugo_symbol, unnest(uniprot_ids) AS u
          FROM gene_enriched WHERE uniprot_ids IS NOT NULL)
""").df()
linked_acc = set(linked["u"])

recon = []


def reconcile(src, col, raw_ids, map_rows, linked_ids):
    row = {"source": src, "id_col": col,
           "n_raw_distinct": len(set(raw_ids)),
           "n_warehouse_rows": map_rows,
           "n_warehouse_linked": len(set(raw_ids) & linked_ids)}
    row["gap"] = row["n_raw_distinct"] - row["n_warehouse_linked"]
    row["FLAG"] = row["gap"] > 0
    recon.append(row)
    return row


ccle_acc = [a.lower() for a, _ in c_hdr]
procan_acc = [a.lower() for a in p_acc]
reconcile("CCLE/Gygi proteomics", "uniprot accession", ccle_acc,
          con.execute("SELECT count(*) FROM protein_map").fetchone()[0], linked_acc)
reconcile("ProCan proteomics", "uniprot accession", procan_acc,
          con.execute("SELECT count(*) FROM procan_protein_map").fetchone()[0], linked_acc)

# RNA layers reconcile on ENSG against `gene`, not on accession.
# The raw DepMap matrix mixes two header conventions in ONE file:
#   'TSPAN6 (ENSG00000000003)'  for annotated genes
#   'ENSG00000288723'           bare, for the unannotated tail
# Reading only the bare form silently sees 422 of 53,961 columns.
gene_ids = set(con.execute("SELECT gene_id FROM gene").df()["gene_id"])
ENSG_RE = re.compile(r"(ensg\d+)")


def ensg_from(cols):
    out = []
    for c in cols:
        m = ENSG_RE.search(str(c).lower())
        if m:
            out.append(m.group(1))
    return out


dm_cols = parquet_cols(DATA / "parquet" /
                       "2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.parquet")
dm_ensg = ensg_from(dm_cols)
dm_named = sum(1 for c in dm_cols if " (" in str(c))
dm_bare = sum(1 for c in dm_cols if ENSG_RE.fullmatch(str(c).lower()) is not None)
print(f"  [raw DepMap header] {len(dm_cols):,} columns = "
      f"{dm_named:,} 'SYMBOL (ENSG)' + {dm_bare:,} bare ENSG + "
      f"{len(dm_cols)-dm_named-dm_bare:,} other")

# GEO is TRANSPOSED: gene-rows x GSM-columns. Its identifiers live in a column.
geo_ids = con.execute(
    f"""SELECT DISTINCT lower("Gene") AS g
        FROM read_parquet('{(DATA / "parquet" / "3_GEOexpression.parquet").as_posix()}')
        WHERE "Gene" IS NOT NULL""").df()["g"].tolist()
geo_ensg = ensg_from(geo_ids)
print(f"  [raw GEO] {len(geo_ids):,} gene rows, {len(geo_ensg):,} parse as ENSG")

reconcile("DepMap RNA-seq", "ensg_id", dm_ensg,
          con.execute("SELECT count(*) FROM gene").fetchone()[0], gene_ids)
reconcile("GEO expression", "ensg_id", geo_ensg,
          con.execute("SELECT count(*) FROM gene").fetchone()[0], gene_ids)

# What the warehouse actually kept, vs what the raw file offered
wh_dm = len([c[0] for c in con.execute("DESCRIBE depmap_expr").fetchall()]) - 1
wh_geo = len([c[0] for c in con.execute("DESCRIBE geo_expr").fetchall()]) - 1
print(f"  [ingest attrition] depmap_expr keeps {wh_dm:,} of {len(dm_ensg):,} raw "
      f"features ({100*wh_dm/max(1,len(dm_ensg)):.1f}%)")
print(f"  [ingest attrition] geo_expr    keeps {wh_geo:,} of {len(geo_ensg):,} raw "
      f"features ({100*wh_geo/max(1,len(geo_ensg)):.1f}%)")
R["T2_rna_attrition"] = {
    "depmap_raw_columns": len(dm_cols), "depmap_raw_ensg": len(dm_ensg),
    "depmap_header_named": dm_named, "depmap_header_bare": dm_bare,
    "depmap_warehouse_features": wh_dm,
    "geo_raw_gene_rows": len(geo_ids), "geo_raw_ensg": len(geo_ensg),
    "geo_warehouse_features": wh_geo,
}
hpa_g = set(con.execute("SELECT DISTINCT gene FROM hpa_rna").df()["gene"])
reconcile("HPA RNA", "ensg_id", hpa_g,
          con.execute("SELECT count(*) FROM gene").fetchone()[0], gene_ids)
cna_g = set(con.execute("SELECT DISTINCT gene_id FROM cosmic_cna "
                        "WHERE gene_id IS NOT NULL").df()["gene_id"])
cna_sym = set(con.execute("SELECT DISTINCT cosmic_gene_symbol FROM cosmic_cna").df()
              ["cosmic_gene_symbol"])
recon.append({"source": "COSMIC CNA", "id_col": "cosmic_gene_symbol",
              "n_raw_distinct": len(cna_sym),
              "n_warehouse_rows": len(cna_sym),
              "n_warehouse_linked": len(cna_g),
              "gap": len(cna_sym) - len(cna_g),
              "FLAG": len(cna_sym) - len(cna_g) > 0})

rec = pd.DataFrame(recon)
print(rec.to_string(index=False))
R["T2_reconciliation"] = recon

sub("classify the CCLE + ProCan gaps")


def classify_gap(acc_list, sym_lookup, src):
    """unmatched (safe) | composite (Class A) | resolved-to-one-winner (Class B)"""
    lost = [a for a in acc_list if a not in linked_acc]
    cls = Counter()
    detail = defaultdict(list)
    for a in lost:
        s = sym_lookup.get(a)
        if s is None or str(s).strip() == "" or str(s).lower() == "nan":
            cls["no_symbol_in_source"] += 1
        elif DELIMS.search(str(s)):
            cls["composite_ClassA"] += 1
            detail["composite_ClassA"].append(f"{a} -> {s}")
        else:
            cls["symbol_present_but_unmatched"] += 1
            detail["symbol_present_but_unmatched"].append(f"{a} -> {s}")
    print(f"  {src}: {len(lost):,} accessions unreachable from any gene")
    for k, v in cls.most_common():
        print(f"      {k:<34s} {v:>6,}")
    return {"n_lost": len(lost), "classes": dict(cls),
            "examples": {k: v[:12] for k, v in detail.items()}}


R["T2_ccle_gap"] = classify_gap(ccle_acc, {a.lower(): s for a, s in c_hdr},
                                "CCLE/Gygi")
R["T2_procan_gap"] = classify_gap(procan_acc,
                                  {a.lower(): s for a, s in zip(p_acc, p_sym)},
                                  "ProCan")

sub("Q1 -- accessions mapping to >1 gene (HGNC release as UniProt proxy)")
# HGNC's uniprot_ids is a per-symbol list. Invert it: an accession listed under
# two approved symbols is an accession that cannot identify a single gene.
inv = defaultdict(set)
for sym, ids in zip(hg["symbol"], hg["uniprot_ids"]):
    if isinstance(ids, str):
        for a in re.split(r"[|,;]", ids):
            a = a.strip().lower()
            if a:
                inv[a].add(sym)
multi = {a: sorted(s) for a, s in inv.items() if len(s) > 1}
print(f"  HGNC accessions total                : {len(inv):,}")
print(f"  accessions under >1 approved symbol  : {len(multi):,}")
ccle_multi = {a: multi[a] for a in set(ccle_acc) & set(multi)}
procan_multi = {a: multi[a] for a in set(procan_acc) & set(multi)}
ccle_multi_linked = {a: v for a, v in ccle_multi.items() if a in linked_acc}
procan_multi_linked = {a: v for a, v in procan_multi.items() if a in linked_acc}
print(f"  of CCLE's {len(set(ccle_acc)):,} accessions   : {len(ccle_multi):,} ambiguous "
      f"({len(ccle_multi_linked):,} live in gene_enriched)")
print(f"  of ProCan's {len(set(procan_acc)):,} accessions : {len(procan_multi):,} ambiguous "
      f"({len(procan_multi_linked):,} live in gene_enriched)")
for a, v in sorted(ccle_multi_linked.items())[:15]:
    print(f"      {a:<12s} -> {', '.join(v)}")
R["T2_Q1"] = {
    "proxy": "HGNC gene_with_protein_product.txt uniprot_ids column, inverted",
    "hgnc_accessions": len(inv),
    "hgnc_accessions_multi_gene": len(multi),
    "ccle_n_accessions": len(set(ccle_acc)),
    "ccle_ambiguous": len(ccle_multi),
    "ccle_ambiguous_linked": len(ccle_multi_linked),
    "procan_n_accessions": len(set(procan_acc)),
    "procan_ambiguous": len(procan_multi),
    "procan_ambiguous_linked": len(procan_multi_linked),
    "ccle_ambiguous_linked_examples": {a: v for a, v in
                                       sorted(ccle_multi_linked.items())[:40]},
}


# ================================================================ T3
head("T3 -- Reverse direction: one gene, many features")

t3 = {}
for label, q in [
    ("protein_map (CCLE accessions per gene)",
     """SELECT g.gene_id, count(*) n, list(m.uniprot_id) f
        FROM protein_map m JOIN gene g ON g.hugo_symbol = m.gene_symbol
        GROUP BY 1 HAVING count(*) > 1 ORDER BY 2 DESC"""),
    ("procan_protein_map (ProCan accessions per gene)",
     """SELECT g.gene_id, count(*) n, list(m.uniprot_id) f
        FROM procan_protein_map m JOIN gene g ON g.hugo_symbol = m.gene_symbol
        GROUP BY 1 HAVING count(*) > 1 ORDER BY 2 DESC"""),
    ("gene_enriched.uniprot_ids (what production actually averages)",
     """SELECT gene_id, n_uniprot_ids n, uniprot_ids f FROM gene_enriched
        WHERE n_uniprot_ids > 1 ORDER BY 2 DESC"""),
]:
    d = con.execute(q).df()
    print(f"  {label:<52s} genes with >1 feature = {len(d):>6,}"
          f"   max = {int(d['n'].max()) if len(d) else 0}")
    t3[label] = {"n_genes": len(d),
                 "max_features": int(d["n"].max()) if len(d) else 0,
                 "top": [{"gene_id": r.gene_id, "n": int(r.n),
                          "features": [str(x) for x in list(r.f)]}
                         for r in d.head(8).itertuples()]}

sub("modalities where a per-gene aggregation rule is required but absent")
# miRNA columns are mirbase-style ids with no gene mapping ingested at all;
# methylation is not ingested. Both are recorded as gaps, not measured.
mir_cols = [c[0] for c in con.execute("DESCRIBE mirna").fetchall()][1:]
print(f"  mirna: {len(mir_cols):,} feature columns, no gene mapping table in the "
      f"warehouse -> aggregation rule undefined")
meth = (DATA / "COSMIC" / "Cosmic_CompleteDifferentialMethylation_v104_GRCh38.tsv")
print(f"  methylation: source present ({meth.exists()}), not ingested -> "
      f"CpG-per-gene rule undefined")
t3["undefined_aggregation"] = {"mirna_features": len(mir_cols),
                               "methylation_source_present": meth.exists(),
                               "methylation_ingested": False}
R["T3"] = t3


# ================================================================ T4
head("T4 -- Paralog-family annotation from HGNC gene_group (Class C)")

fam = defaultdict(set)
for sym, grp in zip(hg["symbol"], hg["gene_group"]):
    if isinstance(grp, str):
        for g in grp.split("|"):
            g = g.strip()
            if g:
                fam[g].add(sym)
multi_fam = {g: s for g, s in fam.items() if len(s) > 1}
in_family = set().union(*multi_fam.values()) if multi_fam else set()
print(f"  HGNC gene groups                    : {len(fam):,}")
print(f"  groups with >1 member               : {len(multi_fam):,}")
print(f"  symbols in >=1 multi-member group   : {len(in_family):,} of {len(hg):,} "
      f"({100*len(in_family)/len(hg):.1f}%)")

scored = con.execute("SELECT DISTINCT hugo_symbol FROM gene_enriched "
                     "WHERE hugo_symbol IS NOT NULL").df()["hugo_symbol"]
scored = set(scored)
hit = {s.lower() for s in in_family} & scored
print(f"  warehouse genes in a paralog family : {len(hit):,} of {len(scored):,} "
      f"({100*len(hit)/max(1,len(scored)):.1f}%)")

sub("families the audit predicted would be large")
pred = ["Histones", "Hemoglobin", "Cancer/testis antigen", "MAGE", "Tubulins",
        "Ribosomal proteins", "Ubiquitin", "Survival motor neuron"]
fam_report = []
for p in pred:
    match = {g: s for g, s in multi_fam.items() if p.lower() in g.lower()}
    for g, s in sorted(match.items(), key=lambda kv: -len(kv[1]))[:3]:
        n_scored = len({x.lower() for x in s} & scored)
        print(f"  {g[:52]:<52s} members={len(s):>4}  in warehouse={n_scored:>4}")
        fam_report.append({"group": g, "members": len(s), "in_warehouse": n_scored})
R["T4"] = {"hgnc_groups": len(fam), "multi_member_groups": len(multi_fam),
           "symbols_in_family": len(in_family),
           "warehouse_genes_in_family": len(hit),
           "warehouse_genes_total": len(scored),
           "predicted_families": fam_report,
           "note": "HGNC gene_group is a curated family grouping, NOT a CDS "
                   "identity threshold. >95% CDS identity needs Ensembl "
                   "paralogy (not available offline). This is an upper-bound "
                   "screen, not the T4 spec."}


# ================================================================ T5
head("T5 -- Replication-dependent histones / poly-A depletion (Class C)")

hist_groups = {g: s for g, s in fam.items() if "histone" in g.lower()}
hist_syms = set().union(*hist_groups.values()) if hist_groups else set()
# Replication-dependent = the clustered H1/H2A/H2B/H3/H4 genes (H*C*, H2A*, etc.)
# as opposed to the replication-INdependent variants (H2AZ1, H3-3A, H1-10...).
rep_dep = {s for s in hist_syms if re.fullmatch(r"H[1-4]C\d+[A-Z]?", s)
           or re.fullmatch(r"H2[AB][CU]\d+[A-Z]?", s)}
print(f"  HGNC histone groups                  : {len(hist_groups)}")
print(f"  histone symbols                      : {len(hist_syms)}")
print(f"  replication-dependent (clustered)    : {len(rep_dep)}")

ens = dict(zip(hg["symbol"], hg["ensembl_gene_id"]))
rep_ensg = {str(ens[s]).lower() for s in rep_dep if isinstance(ens.get(s), str)}
raw_set = set(dm_ensg)
# The measurement has to run on the INGESTED matrix, whose columns are bare ENSG.
wh_set = {c[0] for c in con.execute("DESCRIBE depmap_expr").fetchall()} - {"profileid"}
rep_in_raw = sorted(rep_ensg & raw_set)
rep_in_dm = sorted(rep_ensg & wh_set)
print(f"  in the raw DepMap matrix             : {len(rep_in_raw)} of {len(rep_ensg)}")
print(f"  surviving into depmap_expr           : {len(rep_in_dm)} of {len(rep_ensg)}")

t5 = {"histone_groups": len(hist_groups), "histone_symbols": len(hist_syms),
      "replication_dependent": len(rep_dep),
      "rep_dep_with_ensg": len(rep_ensg),
      "rep_dep_in_raw_matrix": len(rep_in_raw),
      "rep_dep_in_warehouse": len(rep_in_dm)}

if rep_in_dm:
    # Median log2(TPM+1) of replication-dependent histones vs a background of
    # 600 random protein-coding genes. Poly-A selection predicts depletion.
    # Background must be protein-coding too, or it fills with unexpressed
    # non-coding genes and the comparison flatters the histones.
    coding = {str(ens[s]).lower() for s, lg in zip(hg["symbol"], hg["locus_group"])
              if isinstance(ens.get(s), str) and lg == "protein-coding gene"}
    bg = sorted((wh_set & coding) - rep_ensg)
    rng = np.random.default_rng(0)
    bg = [bg[i] for i in rng.choice(len(bg), size=600, replace=False)]

    def med(cols):
        sel = ", ".join(f'median("{c}")' for c in cols)
        v = con.execute(f"SELECT {sel} FROM depmap_expr").fetchone()
        return np.array([x for x in v if x is not None], dtype=float)

    h, b = med(rep_in_dm), med(bg)
    from scipy.stats import mannwhitneyu
    u = mannwhitneyu(h, b, alternative="less")
    frac_below_med = float((h < np.median(b)).mean())
    print(f"  median log2(TPM+1), rep-dep histones : {np.median(h):.4f}  "
          f"(n={len(h)} genes)")
    print(f"  median log2(TPM+1), coding background: {np.median(b):.4f}  "
          f"(n={len(b)} genes)")
    print(f"  delta                                : {np.median(h)-np.median(b):+.4f}")
    print(f"  Mann-Whitney (histones < background) : U={u.statistic:.0f}  "
          f"p={u.pvalue:.3e}")
    print(f"  histones below the background median : {100*frac_below_med:.1f}%")
    print(f"  histones below background quartile 1 : "
          f"{100*float((h < np.percentile(b, 25)).mean()):.1f}%")
    t5.update({"median_histone": float(np.median(h)),
               "median_background_coding": float(np.median(b)),
               "delta": float(np.median(h) - np.median(b)),
               "mannwhitney_U": float(u.statistic), "mannwhitney_p": float(u.pvalue),
               "frac_histones_below_bg_median": frac_below_med,
               "frac_histones_below_bg_q1":
                   float((h < np.percentile(b, 25)).mean()),
               "note": "Protocol-determined and identical across every cell "
                       "line, so it is invisible to any within-gene percentile "
                       "check by construction. Depends on Q2 (was the DepMap "
                       "library poly-A selected?) -- unverified here."})
R["T5"] = t5


# ================================================================ T6
head("T6 -- Identifier namespace collisions")

PATTERNS = [("ACH", r"^ach-\d+"), ("SIDM", r"^sidm\d+"), ("SIDG", r"^sidg\d+"),
            ("ENSG", r"^ensg\d+"), ("COSG", r"^cosg\d+"), ("COSM", r"^cosm\d+"),
            ("CVCL", r"^cvcl_"), ("PR", r"^pr-"), ("GSM", r"^gsm\d+"),
            ("HGNC", r"^hgnc:"), ("UNIPROT", r"^[opq][0-9][a-z0-9]{3}[0-9]$|"
                                             r"^[a-n,r-z][0-9]{5}$")]


def id_kind(vals):
    kinds = Counter()
    for v in vals:
        # list-valued columns (gene.uniprot_ids et al.) must be flattened, or
        # every one of them classifies as 'other' and fakes a collision
        items = list(v) if isinstance(v, (list, np.ndarray)) else [v]
        for it in items:
            s = str(it).strip().lower()
            for name, pat in PATTERNS:
                if re.match(pat, s):
                    kinds[name] += 1
                    break
            else:
                kinds["other"] += 1
    return kinds.most_common(1)[0][0] if kinds else "empty"


tables = [r[0] for r in con.execute("SHOW TABLES").fetchall()]
# only narrow tables -- the wide omics matrices have 8k-20k feature columns
narrow = [t for t in tables
          if len(con.execute(f"DESCRIBE {t}").fetchall()) < 200]
where = defaultdict(list)
for t in narrow:
    cols = [c[0] for c in con.execute(f"DESCRIBE {t}").fetchall()]
    for c in cols:
        if re.search(r"(^|_)(id|ids|name|names|symbol|accession)s?$", c) or \
           c in ("gene", "sample", "model_id", "profileid"):
            try:
                vals = con.execute(f'SELECT "{c}" FROM {t} '
                                   f'WHERE "{c}" IS NOT NULL LIMIT 200').df()[c]
            except Exception:
                continue
            if len(vals) == 0:
                continue
            where[c].append({"table": t, "kind": id_kind(vals),
                             "sample": str(vals.iloc[0])[:28]})

collisions = {c: v for c, v in where.items()
              if len({x["kind"] for x in v}) > 1}
print(f"  identifier-ish columns examined : {len(where)}")
print(f"  column names holding >1 id system: {len(collisions)}")
for c, v in sorted(collisions.items(), key=lambda kv: -len(kv[1])):
    kinds = sorted({x["kind"] for x in v})
    print(f"\n  {c}   -> {kinds}")
    for x in v:
        print(f"      {x['table']:<34s} {x['kind']:<9s} e.g. {x['sample']}")
R["T6"] = {"columns_examined": len(where),
           "colliding_columns": {c: v for c, v in collisions.items()}}

sub("the known ProCan case, checked at source")
proc_hdr_row2 = None
with open(PROCAN, encoding="utf-8") as f:
    next(f); next(f)
    proc_hdr_row2 = next(f).rstrip("\n").split("\t")[:2]
print(f"  raw ProCan header row 2 names its key columns {proc_hdr_row2}")
print(f"  but that 'model_id' holds SIDM ids, while model_id everywhere in the "
      f"warehouse holds ACH")
R["T6"]["procan_raw_key_columns"] = proc_hdr_row2


# ================================================================ T7
head("T7 -- Multi-target reagent provenance")

t7 = {}
sub("Chronos release identification")
try:
    import h5py
    for f in sorted(CHRONOS.glob("*.hdf5")):
        with h5py.File(f, "r") as h:
            keys = list(h.keys())
            attrs = {k: str(v) for k, v in h.attrs.items()}
        print(f"  {f.name}")
        print(f"      datasets: {keys}")
        print(f"      root attrs: {attrs if attrs else '(none)'}")
        t7[f.name] = {"datasets": keys, "attrs": attrs}
except Exception as e:                                   # pragma: no cover
    print(f"  h5py unavailable: {e}")
print("  NOTE: no release/version string is embedded in these files. Q3 cannot "
      "be answered from the data -- it needs the DepMap release notes.")

sub("GDSC drug identifier duplication")
dup_name = cmpd.groupby(cmpd["DRUG_NAME"].str.lower())["DRUG_ID"].nunique()
dup_name = dup_name[dup_name > 1]
print(f"  compounds: {len(cmpd):,} rows, {cmpd.DRUG_ID.nunique():,} DRUG_IDs, "
      f"{cmpd.DRUG_NAME.str.lower().nunique():,} distinct names")
print(f"  drug names carrying >1 DRUG_ID  : {len(dup_name)}")
for n, k in dup_name.sort_values(ascending=False).head(12).items():
    ids = sorted(cmpd[cmpd.DRUG_NAME.str.lower() == n].DRUG_ID.unique())
    sites = sorted(cmpd[cmpd.DRUG_NAME.str.lower() == n].SCREENING_SITE.unique())
    print(f"      {n:<24s} ids={ids} sites={sites}")
t7["gdsc"] = {"rows": len(cmpd), "n_drug_ids": int(cmpd.DRUG_ID.nunique()),
              "n_names": int(cmpd.DRUG_NAME.str.lower().nunique()),
              "names_with_multiple_ids": len(dup_name),
              "examples": {str(n): sorted(cmpd[cmpd.DRUG_NAME.str.lower() == n]
                                          .DRUG_ID.unique())
                           for n in dup_name.sort_values(ascending=False)
                           .head(20).index}}

sub("miRNA feature namespace")
print(f"  mirna columns look like {mir_cols[:3]} -- these are not mirbase "
      f"accessions and carry no family annotation, so sequence-identical "
      f"family members cannot be detected from the warehouse")
t7["mirna"] = {"n_features": len(mir_cols), "example_ids": mir_cols[:5],
               "family_annotation_present": False}
R["T7"] = t7


# ================================================================ T8
head("T8 -- Frozen regression assertions")

# C6: the SQL predicate is DERIVED from the same constant the T1 scan used, so
# the two cannot drift. Guarded by a fixture test that fails the build if they
# ever do -- 'cnk3/ipcef1' is one of the fixtures.
assert_detectors_agree(con)
print("  detector equivalence (python vs SQL) : PASS")
DELIM_SQL = composite_sql("gene_symbol")
asserts = [
    ("procan composite symbols",
     f"SELECT count(*) FROM procan_protein_map WHERE {DELIM_SQL}",
     len(t1["procan_protein_map"])),
    ("protein_map composite symbols",
     f"SELECT count(*) FROM protein_map WHERE {DELIM_SQL}",
     len(t1["protein_map"])),
    ("protein_map symbol-less accessions",
     "SELECT count(*) FROM protein_map WHERE gene_symbol IS NULL", None),
    ("accessions claimed by >1 gene in gene_enriched",
     """SELECT count(*) FROM (SELECT u FROM (SELECT gene_id, unnest(uniprot_ids) u
        FROM gene_enriched) GROUP BY u HAVING count(DISTINCT gene_id) > 1)""", 0),
    ("genes with >1 accession averaged before ranking",
     "SELECT count(*) FROM gene_enriched WHERE n_uniprot_ids > 1", None),
]
frozen = []
for label, q, expect in asserts:
    got = con.execute(q).fetchone()[0]
    ok = "" if expect is None else ("PASS" if got == expect else "FAIL")
    print(f"  {label:<52s} = {got:>6}   {ok}")
    frozen.append({"assertion": label, "sql": " ".join(q.split()),
                   "observed": int(got), "expected": expect})
R["T8"] = frozen


# ================================================================ save
outp = OUT / "test_run_identifier_ambiguity_audit_results.json"
outp.write_text(json.dumps(R, indent=1, default=str), encoding="utf-8")
head("done")
print(f"results -> {outp}")
con.close()
