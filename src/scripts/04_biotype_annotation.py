"""
04_biotype_annotation.py
CellLineSelector - Feature Engineering: gene biotype annotation (closes L4, activates E6).

WHY
---
Decision E6 (Expression Decision Note) stores all ~53,961 DepMap genes but makes the
recommender's DEFAULT gene universe PROTEIN-CODING, with non-coding genes queryable but
flagged. That default needs a biotype per ENSG. Biotype cannot be derived from expression
data - it is external annotation - so this is the single static lookup E6 committed to:
a pinned Ensembl release GTF, joined on ENSG. Pinning the release makes it reproducible and
citable ("Ensembl release N, GRCh38").

WHAT IT DOES (measures/annotates only - no method choice)
---------------------------------------------------------
1. Parses the `gene` lines of an Ensembl (or GENCODE) GTF -> gene_id, gene_biotype, gene_name.
2. Normalises ids on both sides (upper-case, strip version suffix) and left-joins to your ENSG set.
3. Emits per-gene: gene_biotype (raw), biotype_coarse, gene_name, is_protein_coding,
   in_default_universe. Unmatched ids -> 'unmatched' (conservatively NOT in the default universe).
4. Reports match rate and the biotype distribution, and saves features/gene_biotype.parquet.

GET THE GTF (one download; pin the release for provenance)
----------------------------------------------------------
Ensembl (matches DepMap's ENSG namespace directly):
    wget https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
GENCODE (equivalent; uses `gene_type` + versioned ENSG - the parser handles both):
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz

release-111 is a concrete example that works; for the tightest match rate, pick the release
closest to your DepMap expression vintage (see DepMap release notes). The reported match rate
below tells you if the release is a poor fit - if it is low (<95% for real genes), bump/switch release.

USAGE (in your Jupyter session)
-------------------------------
    import importlib
    bt = importlib.import_module("04_biotype_annotation")

    ensg_ids = [c for c in depmap_expr.columns if str(c).lower().startswith("ensg")]
    ann = bt.build_biotype_annotation(ensg_ids,
                                      gtf_path="Homo_sapiens.GRCh38.111.gtf.gz")
    # ann is one row per gene; also saved to features/gene_biotype.parquet

ALTERNATIVE (no manual download): bt.biotype_via_biomart(ensg_ids)  # needs `pip install pybiomart`

Author: (CellLineSelector team) | Depends on: pandas; pybiomart optional (fallback path)
"""

from __future__ import annotations

import gzip
import os
import re

import pandas as pd

_GID = re.compile(r'gene_id "([^"]+)"')
_BT = re.compile(r'gene_(?:biotype|type) "([^"]+)"')   # Ensembl=gene_biotype, GENCODE=gene_type
_NM = re.compile(r'gene_name "([^"]+)"')

# small-ncRNA and other Ensembl biotypes grouped for a clean coarse category
_SMALL_NC = {"mirna", "snrna", "snorna", "rrna", "misc_rna", "scarna", "vault_rna",
             "ribozyme", "srp_rna", "trna", "scrna", "srna", "y_rna", "rrna_pseudogene",
             "mt_rrna", "mt_trna"}
_LNC = {"lncrna", "lincrna", "antisense", "processed_transcript", "sense_intronic",
        "sense_overlapping", "macro_lncrna", "bidirectional_promoter_lncrna",
        "3prime_overlapping_ncrna", "non_coding"}


def _norm_ensg(s: pd.Series) -> pd.Series:
    """Upper-case and strip Ensembl version suffix (ENSG000...12 -> ENSG000...)."""
    return s.astype(str).str.strip().str.upper().str.split(".").str[0]


def coarse_biotype(bt: str) -> str:
    """Collapse Ensembl's many biotypes into a small, interpretable vocabulary."""
    b = (bt or "").lower()
    if b == "unmatched":
        return "unmatched"
    if b == "protein_coding":
        return "protein_coding"
    if "pseudogene" in b:
        return "pseudogene"
    if b in _LNC or "lncrna" in b:
        return "lncRNA"
    if b in _SMALL_NC:
        return "small_ncRNA"
    if b.startswith(("ig_", "tr_")):
        return "IG_TR_gene"
    return "other"


def parse_gtf_genes(gtf_path: str) -> pd.DataFrame:
    """Return DataFrame[gene_id, gene_biotype, gene_name] from the `gene` lines of a GTF.

    Streams the file and keeps only feature=='gene' rows, so memory stays small even for a
    full ~50 MB gzip. Works for Ensembl (gene_biotype) and GENCODE (gene_type) GTFs.
    """
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF not found: {gtf_path} - see the download commands in the "
                                f"module docstring")
    op = gzip.open if str(gtf_path).endswith(".gz") else open
    rows = []
    with op(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            # cheap pre-filter before splitting: the 3rd tab-field must be 'gene'
            p = line.split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            attr = p[8]
            gid = _GID.search(attr)
            if not gid:
                continue
            bt = _BT.search(attr)
            nm = _NM.search(attr)
            rows.append((gid.group(1), bt.group(1) if bt else "unknown",
                         nm.group(1) if nm else None))
    df = pd.DataFrame(rows, columns=["gene_id", "gene_biotype", "gene_name"])
    if df.empty:
        raise ValueError(f"No `gene` feature rows parsed from {gtf_path} - is this a gene GTF?")
    return df


def build_biotype_annotation(ensg_ids, gtf_path: str,
                             out: str = "features/gene_biotype.parquet") -> pd.DataFrame:
    """Join biotype onto your ENSG set and emit the annotation table used by E6.

    Parameters
    ----------
    ensg_ids : iterable of str   - your gene ids (e.g. depmap_expr ENSG columns)
    gtf_path : str               - path to a pinned Ensembl/GENCODE GTF (.gtf or .gtf.gz)
    out      : str               - parquet output path

    Returns one row per input gene with: gene (original id), gene_biotype, biotype_coarse,
    gene_name, is_protein_coding, in_default_universe.
    """
    genes = parse_gtf_genes(gtf_path)
    genes["_key"] = _norm_ensg(genes["gene_id"])
    genes = genes.drop_duplicates("_key", keep="first")

    ann = pd.DataFrame({"gene": pd.Index(list(ensg_ids)).astype(str)})
    ann["_key"] = _norm_ensg(ann["gene"])
    merged = ann.merge(genes[["_key", "gene_biotype", "gene_name"]], on="_key", how="left")
    merged["gene_biotype"] = merged["gene_biotype"].fillna("unmatched")
    merged["biotype_coarse"] = merged["gene_biotype"].map(coarse_biotype)
    merged["is_protein_coding"] = merged["gene_biotype"].eq("protein_coding")
    merged["in_default_universe"] = merged["is_protein_coding"]     # E6 default view
    merged = merged.drop(columns="_key")

    _report(merged, gtf_path)

    if out:
        os.makedirs(os.path.dirname(out), exist_ok=True)
        merged.to_parquet(out, index=False)
        print(f"[biotype] wrote {out}")
    return merged


def _report(merged: pd.DataFrame, gtf_path: str):
    n = len(merged)
    matched = int((merged["gene_biotype"] != "unmatched").sum())
    pc = int(merged["is_protein_coding"].sum())
    print(f"[biotype] source: {os.path.basename(gtf_path)}")
    print(f"[biotype] genes: {n:,} | matched {matched:,} ({matched/n*100:.1f}%) | "
          f"unmatched {n-matched:,}")
    print(f"[biotype] PROTEIN-CODING (default universe): {pc:,} ({pc/n*100:.1f}%)")
    print("[biotype] coarse biotype distribution:")
    dist = merged["biotype_coarse"].value_counts()
    for k, v in dist.items():
        print(f"           {k:<16} {v:>7,} ({v/n*100:5.1f}%)")
    if matched / n < 0.95:
        print("[biotype] NOTE: match rate <95% - the GTF release may not match your DepMap "
              "vintage; try the release closest to it (see DepMap release notes).")


def biotype_via_biomart(ensg_ids, out: str = "features/gene_biotype.parquet") -> pd.DataFrame:
    """Fallback: query Ensembl BioMart live (needs `pip install pybiomart`, internet).

    Less reproducible than a pinned GTF (BioMart tracks current Ensembl), so prefer
    build_biotype_annotation for the dissertation record. Handy for a quick first pass.
    """
    from pybiomart import Server                                    # noqa: local import
    server = Server(host="http://www.ensembl.org")
    ds = server.marts["ENSEMBL_MART_ENSEMBL"].datasets["hsapiens_gene_ensembl"]
    bm = ds.query(attributes=["ensembl_gene_id", "gene_biotype", "external_gene_name"])
    bm.columns = ["gene_id", "gene_biotype", "gene_name"]
    bm["_key"] = _norm_ensg(bm["gene_id"])
    bm = bm.drop_duplicates("_key")

    ann = pd.DataFrame({"gene": pd.Index(list(ensg_ids)).astype(str)})
    ann["_key"] = _norm_ensg(ann["gene"])
    merged = ann.merge(bm[["_key", "gene_biotype", "gene_name"]], on="_key", how="left")
    merged["gene_biotype"] = merged["gene_biotype"].fillna("unmatched")
    merged["biotype_coarse"] = merged["gene_biotype"].map(coarse_biotype)
    merged["is_protein_coding"] = merged["gene_biotype"].eq("protein_coding")
    merged["in_default_universe"] = merged["is_protein_coding"]
    merged = merged.drop(columns="_key")
    _report(merged, "BioMart (live, current Ensembl)")
    if out:
        os.makedirs(os.path.dirname(out), exist_ok=True)
        merged.to_parquet(out, index=False)
        print(f"[biotype] wrote {out}")
    return merged


if __name__ == "__main__":
    print(__doc__)
