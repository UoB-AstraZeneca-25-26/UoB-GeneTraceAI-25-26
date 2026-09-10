"""
build_altered_uncertainty_list.py
---------------------------------
Generates the pre-committed uncertainty list for Fix 1.

THE RULE, fixed before any call is inspected:

    CGC Tier 1
    AND ROLE_IN_CANCER contains "oncogene"
    AND MUTATION_TYPES contains "O"

Source: data/COSMIC/Cosmic_CancerGeneCensus_v104_GRCh37.tsv (768 genes, Tier 1 =
592). NOTE: the GRCh38 copy of the same file in this repo carries only 99 genes
and is missing EGFR, TP53, KRAS, PTEN and CTNNB1 -- it is a fragment and is not
used. `build_gene_roles.py` also reads the GRCh37 file, so this is consistent
with `reference/gene_lookup.parquet`.

WHY THESE THREE CONDITIONS
--------------------------
  Tier 1        the census's own high-confidence set
  oncogene      an in-frame indel that matters is an ACTIVATING event; for a
                pure tumour suppressor, loss-of-function is already caught by
                the HIGH-impact classes (nonsense, frameshift, splice)
  "O" in
  MUTATION_TYPES
                CGC codes mutation classes as A (amplification), D (large
                deletion), F (frameshift), Mis (missense), N (nonsense),
                O (other), S (splice), T (translocation). There is no dedicated
                in-frame-indel code, so in-frame indels and complex coding
                events fall under "O". INTERPRETATION, not a documented CGC
                statement -- flagged as such in the spec.

The rule excludes KRAS (MUTATION_TYPES = "Mis" only), which is correct: a
MODERATE-impact KRAS row is not a missed activating indel, because KRAS is not
altered by that mechanism. It includes EGFR ("A, O, Mis") and CTNNB1
("Mis, O, T").

WHAT THE LIST IS NOT
--------------------
Gene-level, not variant-level. It flags genes where MODERATE-impact calls are
untrustworthy. It does not identify activating variants and makes no claim about
any specific mutation. A weaker instrument, labelled as one.

    python src/pipeline/build_altered_uncertainty_list.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CGC = ROOT / "data" / "COSMIC" / "Cosmic_CancerGeneCensus_v104_GRCh37.tsv"
OUT = ROOT / "src" / "pipeline" / "outputs" / "altered_uncertainty_list.json"

d = pd.read_csv(CGC, sep="\t")
print(f"CGC GRCh37: {len(d)} genes, Tier 1 = {int((d.TIER == 1).sum())}")

role = d.ROLE_IN_CANCER.fillna("").astype(str).str.lower()
mtyp = d.MUTATION_TYPES.fillna("").astype(str)


def has_O(s):
    """'O' as a standalone comma-separated token, not the O inside 'Mis'/'Other'."""
    return "O" in [t.strip() for t in s.split(",")]


sel = (d.TIER == 1) & role.str.contains("oncogene") & mtyp.map(has_O)
genes = sorted(set(d.loc[sel, "GENE_SYMBOL"].astype(str).str.upper()))
print(f"Tier 1 AND oncogene AND 'O' in MUTATION_TYPES: {len(genes)} genes")

gl = pd.read_parquet(ROOT / "reference" / "gene_lookup.parquet",
                     columns=["ensg_id", "hgnc_symbol"])
gl["hgnc_symbol"] = gl.hgnc_symbol.astype(str).str.upper()
gl["ensg_id"] = gl.ensg_id.astype("string").str.split(".").str[0].str.lower()
m = gl[gl.hgnc_symbol.isin(genes)].drop_duplicates("hgnc_symbol")
ensg = sorted(set(m.ensg_id))
missing = sorted(set(genes) - set(m.hgnc_symbol))
print(f"mapped to ensg: {len(ensg)}   unmapped symbols: {len(missing)}")
if missing:
    print("  unmapped:", ", ".join(missing[:20]))

payload = {
    "rule": "CGC Tier 1 AND role contains 'oncogene' AND 'O' in MUTATION_TYPES",
    "source": str(CGC.relative_to(ROOT)),
    "cgc_rows": int(len(d)),
    "n_symbols": len(genes),
    "n_ensg": len(ensg),
    "unmapped_symbols": missing,
    "symbols": genes,
    "ensg": ensg,
    "caveat": ("gene-level, not variant-level; flags genes where MODERATE-impact "
               "calls are untrustworthy, does not identify activating variants"),
}
OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(f"\nwrote {OUT}")
print("\nfirst 30 symbols:")
print("  " + ", ".join(genes[:30]))
