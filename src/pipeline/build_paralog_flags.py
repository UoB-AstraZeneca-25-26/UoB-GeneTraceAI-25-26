"""
C1 -- paralog ambiguity flags from Ensembl paralogy (real CDS identity).

Why not HGNC gene_group: that is functional curation, not sequence identity. It
put 80.2% of the warehouse in a "family", which is useless as a flag. Ensembl
`perc_id` is the actual query/target CDS identity, so the threshold means
something.

Why the 95 / 100 split (Q4). At 100% identity no information exists that could
assign a read to one member rather than another -- RSEM's EM allocation
(Li & Dewey 2011) emits a confident-looking per-gene TPM that is a probabilistic
split of literally unassignable reads. That is not a degraded measurement, it is
not a measurement. Below 100% some reads span a distinguishing position, so the
estimate carries real signal and degrades smoothly. One threshold flattens a
qualitative boundary into a quantitative one.

    paralog_identical  (max perc_id == 100)  -> no information exists
    paralog_ambiguous  (95 <= max < 100)     -> degraded but informative

Both directions are considered: a read is unassignable if EITHER direction is
identical, so the flag uses max(perc_id_target, perc_id_query).

Reference input: data/ensembl_paralogy/ (see PROVENANCE.json). Read-only against
the warehouse; rewrites outputs/gene_ambiguity_flags.parquet with the paralog
columns filled in.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import duckdb
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ambiguity_detectors import assert_detectors_agree  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "outputs"
PARALOGY = ROOT / "data" / "ensembl_paralogy"
SRC = PARALOGY / "hsapiens_paralogy_biomart_20260806.parquet"

IDENTICAL = 100.0
AMBIGUOUS = 95.0

con = duckdb.connect(str(OUT / "celllineselector.db"), read_only=True)
assert_detectors_agree(con)

p = pd.read_parquet(SRC)
p["gene"] = p["gene"].str.lower()
p["paralog"] = p["paralog"].str.lower()
# a read is unassignable if EITHER direction is identical
p["identity"] = p[["perc_id_target", "perc_id_query"]].max(axis=1)
print(f"paralog pairs: {len(p):,}   genes with >=1 paralog: {p.gene.nunique():,}")

hi = p[p.identity >= AMBIGUOUS]
agg = hi.groupby("gene").agg(
    max_paralog_identity=("identity", "max"),
    n_high_identity_paralogs=("paralog", "nunique"),
).reset_index().rename(columns={"gene": "gene_id"})
ident = (p[p.identity >= IDENTICAL].groupby("gene")["paralog"].nunique()
         .rename("n_identical_paralogs").reset_index()
         .rename(columns={"gene": "gene_id"}))
agg = agg.merge(ident, on="gene_id", how="left")
agg["n_identical_paralogs"] = agg["n_identical_paralogs"].fillna(0).astype(int)
print(f"genes with a paralog >= {AMBIGUOUS}% identity: {len(agg):,}")
print(f"genes with a paralog == 100% identity     : "
      f"{int((agg.n_identical_paralogs > 0).sum()):,}")

# ------------------------------------------------------------- join + verify
flags = pd.read_parquet(OUT / "gene_ambiguity_flags.parquet")
flags = flags.drop(columns=[c for c in flags.columns
                            if c.startswith("paralog_")
                            or c in ("max_paralog_identity",
                                     "n_high_identity_paralogs",
                                     "n_identical_paralogs")], errors="ignore")
flags = flags.merge(agg, on="gene_id", how="left")
flags["n_high_identity_paralogs"] = flags["n_high_identity_paralogs"].fillna(0).astype(int)
flags["n_identical_paralogs"] = flags["n_identical_paralogs"].fillna(0).astype(int)
flags["paralog_identical_flag"] = flags["n_identical_paralogs"] > 0
flags["paralog_ambiguous_flag"] = (flags["max_paralog_identity"] >= AMBIGUOUS) & \
                                  (~flags["paralog_identical_flag"])

# disposition (Q4). polya is deliberately NOT a discount: it is a constant
# offset across every cell line, so within-gene ranks are untouched.
def disposition(r):
    if r.paralog_identical_flag:
        return "exclude"
    if r.paralog_ambiguous_flag:
        return "discount"
    if r.polya_depleted_flag or r.multi_accession_flag or r.ambiguous_accession_flag:
        return "annotate"
    return "ok"


flags["validity_disposition"] = flags.apply(disposition, axis=1)
flags["any_flag"] = flags["validity_disposition"] != "ok"

print("\ndisposition:")
print(flags.validity_disposition.value_counts().to_string())

# ---------------------------------------------------------- sanity checks
print("\nsanity checks (must all be flagged):")
checks = {"H4C1": "h4c1", "H3C1": "h3c1", "H2AC4": "h2ac4", "H2BC5": "h2bc5",
          "HBA1": "hba1", "HBA2": "hba2", "SMN1": "smn1", "SMN2": "smn2"}
neg = {"BRAF": "braf", "PLEC": "plec", "TP53": "tp53", "ERBB2": "erbb2"}
ok = True
for lab, sym in checks.items():
    r = flags[flags.hugo_symbol == sym]
    if r.empty:
        print(f"  {lab:<7s} NOT IN WAREHOUSE")
        continue
    r = r.iloc[0]
    good = bool(r.paralog_identical_flag or r.paralog_ambiguous_flag)
    ok &= good
    print(f"  {lab:<7s} identical={r.paralog_identical_flag!s:<5s} "
          f"max_id={r.max_paralog_identity!s:<9s} n_hi={r.n_high_identity_paralogs:<3} "
          f"-> {r.validity_disposition:<9s} {'PASS' if good else 'FAIL'}")
print("negative controls (must NOT be paralog-flagged):")
for lab, sym in neg.items():
    r = flags[flags.hugo_symbol == sym]
    if r.empty:
        continue
    r = r.iloc[0]
    good = not (r.paralog_identical_flag or r.paralog_ambiguous_flag)
    ok &= good
    print(f"  {lab:<7s} max_id={r.max_paralog_identity!s:<9s} "
          f"-> {r.validity_disposition:<9s} {'PASS' if good else 'FAIL'}")

flags.to_parquet(OUT / "gene_ambiguity_flags.parquet", index=False)
meta = json.loads((OUT / "gene_ambiguity_flags_meta.json").read_text(encoding="utf-8"))
meta["paralog_flags"] = {
    "source": json.loads((PARALOGY / "PROVENANCE.json").read_text(encoding="utf-8")),
    "thresholds": {"paralog_identical": IDENTICAL, "paralog_ambiguous": AMBIGUOUS},
    "direction_rule": "max(perc_id_target, perc_id_query) -- a read is "
                      "unassignable if either direction is identical",
    "n_identical": int(flags.paralog_identical_flag.sum()),
    "n_ambiguous": int(flags.paralog_ambiguous_flag.sum()),
    "dispositions": flags.validity_disposition.value_counts().to_dict(),
    "rationale_citations": ["Li & Dewey 2011, BMC Bioinformatics 12:323 (RSEM "
                            "EM allocation of multi-mapping reads)",
                           "Robert & Watson 2015, Genome Biology 16:177"],
    "formalisation_note": "the 95/100 disposition split and the validity "
                          "channel are this project's own formalisation; no "
                          "published standard exists for the remedy.",
    "sanity_checks_passed": bool(ok),
}
(OUT / "gene_ambiguity_flags_meta.json").write_text(
    json.dumps(meta, indent=1, default=str), encoding="utf-8")
print(f"\nsanity checks overall: {'PASS' if ok else 'FAIL'}")
print(f"-> {OUT / 'gene_ambiguity_flags.parquet'}")
con.close()
