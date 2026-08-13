"""
Builds reference/gene_lookup.parquet with a gene_role column from COSMIC CGC v104.
Run once whenever COSMIC data or gene_lookup changes:

    python src/pipeline/build_gene_roles.py

Then rebuild predictions:

    python src/pipeline/build_full_predictions.py
"""
import pandas as pd
from pathlib import Path

REF = Path("reference")
COSMIC_CGC = Path("data/COSMIC/Cosmic_CancerGeneCensus_v104_GRCh37.tsv")

gl = pd.read_parquet(REF / "gene_lookup.parquet")

# ── Load COSMIC CGC ──────────────────────────────────────────────────────────
cgc = pd.read_csv(COSMIC_CGC, sep="\t", usecols=["GENE_SYMBOL", "ROLE_IN_CANCER", "SYNONYMS"])
cgc.columns = ["hgnc_symbol", "role_raw", "synonyms"]
cgc["hgnc_symbol"] = cgc["hgnc_symbol"].str.lower()

def normalise_role(r):
    if pd.isna(r):
        return "unknown"
    r = r.lower()
    has_tsg = "tsg" in r
    has_onc = "oncogene" in r
    if has_tsg and has_onc:
        return "both"
    if has_tsg:
        return "tsg"
    if has_onc:
        return "oncogene"
    return "unknown"

cgc["gene_role"] = cgc["role_raw"].apply(normalise_role)
cgc = cgc[cgc["gene_role"] != "unknown"]

# Primary join: COSMIC GENE_SYMBOL → gene_lookup hgnc_symbol (case-normalised)
gl_upper = gl[["ensg_id", "hgnc_symbol"]].copy()
gl_upper["_sym_upper"] = gl_upper["hgnc_symbol"].str.lower()

role_dedup = cgc[["hgnc_symbol", "gene_role"]].drop_duplicates("hgnc_symbol")

# Pass 1 — direct symbol match (join on upper-cased symbol)
gl_sym = gl_upper[["ensg_id", "_sym_upper"]].rename(columns={"_sym_upper": "hgnc_symbol"})
role_df = role_dedup.merge(gl_sym, on="hgnc_symbol", how="inner")
mapped_ensg = set(role_df["ensg_id"])

# Pass 2 — alias fallback: strip version suffix from ENSG IDs in SYNONYMS column
unmatched = role_dedup[~role_dedup["hgnc_symbol"].isin(role_df["hgnc_symbol"])]
if len(unmatched):
    # explode SYNONYMS and try bare ENSG ID match
    syn = cgc[["hgnc_symbol","gene_role","synonyms"]].copy()
    syn = syn[syn["hgnc_symbol"].isin(unmatched["hgnc_symbol"])]
    syn = syn.dropna(subset=["synonyms"])
    syn = syn.assign(syn_token=syn["synonyms"].str.split(",")).explode("syn_token")
    syn["syn_token"] = syn["syn_token"].str.strip()
    # strip version (e.g. ensg00000148584.14 -> ensg00000148584)
    syn["ensg_bare"] = syn["syn_token"].str.extract(r"^(ENSG\d+)")
    ensg_hits = syn.dropna(subset=["ensg_bare"]).merge(
        gl_upper[["ensg_id"]], left_on="ensg_bare", right_on="ensg_id", how="inner"
    )[["hgnc_symbol","gene_role","ensg_id"]].drop_duplicates("ensg_id")
    ensg_hits = ensg_hits[~ensg_hits["ensg_id"].isin(mapped_ensg)]
    role_df = pd.concat([role_df, ensg_hits], ignore_index=True)

print(f"Pass 1 (symbol): {len(set(role_df['ensg_id']) - (set(role_df['ensg_id']) - mapped_ensg))} genes")
print(f"Pass 2 (ENSG alias): {len(role_df) - len(mapped_ensg)} additional genes")

print(f"COSMIC CGC genes mapped to ensg_id: {len(role_df)}")
print(role_df.gene_role.value_counts().to_dict())

# Join onto gene_lookup
gl = gl.drop(columns=["gene_role"], errors="ignore")
gl_out = gl.merge(role_df[["ensg_id", "gene_role"]], on="ensg_id", how="left")
gl_out["gene_role"] = gl_out["gene_role"].fillna("unknown")

print("\nFull lookup role distribution:")
print(gl_out.gene_role.value_counts().to_dict())

for sym in ["PTEN", "TP53", "RB1", "BRCA1", "APC", "BRAF", "KRAS", "BCL2", "NOTCH1"]:
    row = gl_out[gl_out.hgnc_symbol == sym]
    role = row["gene_role"].iloc[0] if len(row) else "NOT FOUND"
    print(f"  {sym}: {role}")

gl_out.to_parquet(REF / "gene_lookup.parquet", index=False)
print("\ngene_lookup.parquet updated with COSMIC gene_role column.")
print("Next: python src/pipeline/build_full_predictions.py")
