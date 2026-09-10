"""
build_evidence_ledger.py
--------------------------
Bulk evidence ledger for all 28.4M (model_id, ensg_id) rows in
predictions_with_confidence.parquet. A projection of existing columns
(Stage 1-6), not new modelling -- see evidence_ledger.py for the field
derivations shared with explain_pair.py's single-pair path.

Output: src/pipeline/outputs/evidence_ledger.parquet
Columns: model_id, ensg_id, confidence (sort key, unchanged), contradicted,
         layers_present, n_layers, driver_alteration, cna_role_consistent,
         essentiality_check, census_role, validation_source

Batched by gene (same BATCH=200 pattern as build_full_predictions.py) to
keep memory bounded at this row count.
"""
import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
import sys
sys.path.insert(0, "src/pipeline")
from evidence_ledger import (
    rna_covered_model_ids, derive_layers_present, derive_driver_alteration,
    derive_cna_role_consistent, derive_essentiality_check, derive_census_role,
    derive_validation_source, derive_signal_spread, LEDGER_FIELD_ORDER,
)

OUTPUTS = Path("src/pipeline/outputs")

print("Loading RNA-coverage set...")
rna_covered = rna_covered_model_ids()
print(f"  {len(rna_covered):,} cell lines with an RNA profile")

print("Loading chronos_validation (per-gene rho/pval)...")
cv = pd.read_parquet(OUTPUTS / "chronos_validation.parquet")
cv_map = cv.set_index("ensg_id")[["rho", "pval"]].to_dict("index")

gene_list = pd.read_parquet(OUTPUTS / "core_score.parquet", columns=["ensg_id"])["ensg_id"].unique()
print(f"{len(gene_list):,} genes to process")

BATCH = 200
writer = None
n_contradicted_rows = 0

for batch_start in range(0, len(gene_list), BATCH):
    batch_genes = gene_list[batch_start:batch_start + BATCH]

    batch_genes_list = list(batch_genes)
    pred_b = pd.read_parquet(
        OUTPUTS / "predictions_with_confidence.parquet",
        columns=["model_id", "ensg_id", "confidence", "n_layers", "class",
                 "regime_source", "has_cna_alteration", "gene_role", "chronos_check",
                 "signal_spread", "gene_top_vs_median_fold"],
        filters=[("ensg_id", "in", batch_genes_list)],
    )

    flags_b = pd.read_parquet(
        OUTPUTS / "flags_with_driver.parquet",
        columns=["model_id", "ensg_id", "mut_driver", "fusion_driver", "has_cna_context"],
        filters=[("ensg_id", "in", batch_genes_list)],
    )

    c = pred_b.merge(flags_b, on=["model_id", "ensg_id"], how="left")
    c["mut_driver"] = c["mut_driver"].fillna(False)
    c["fusion_driver"] = c["fusion_driver"].fillna(False)
    c["has_cna_context"] = c["has_cna_context"].fillna(False)

    c["layers_present"] = [
        derive_layers_present(nl, mid, rna_covered)
        for nl, mid in zip(c["n_layers"], c["model_id"])
    ]
    # Overwrites the raw band with its explained form (same field name, ledger text)
    c["signal_spread"] = [
        derive_signal_spread(sp, fold)
        for sp, fold in zip(c["signal_spread"], c["gene_top_vs_median_fold"])
    ]
    c["driver_alteration"] = [
        derive_driver_alteration(m, f) for m, f in zip(c["mut_driver"], c["fusion_driver"])
    ]
    c["cna_role_consistent"] = [
        derive_cna_role_consistent(gr, cls, cna_alt, cna_ctx)
        for gr, cls, cna_alt, cna_ctx in zip(
            c["gene_role"], c["class"], c["has_cna_alteration"], c["has_cna_context"])
    ]
    c["essentiality_check"] = [
        derive_essentiality_check(
            chk, cv_map.get(g, {}).get("rho"), cv_map.get(g, {}).get("pval")
        )
        for chk, g in zip(c["chronos_check"], c["ensg_id"])
    ]
    c["census_role"] = [derive_census_role(gr) for gr in c["gene_role"]]
    c["validation_source"] = [derive_validation_source(rs) for rs in c["regime_source"]]
    c["contradicted"] = c["chronos_check"] == "inverted"
    n_contradicted_rows += int(c["contradicted"].sum())

    out_cols = ["model_id", "ensg_id", "confidence", "contradicted"] + LEDGER_FIELD_ORDER
    table = pa.Table.from_pandas(c[out_cols], preserve_index=False)
    if writer is None:
        writer = pq.ParquetWriter(OUTPUTS / "evidence_ledger.parquet", table.schema)
    writer.write_table(table)

    done = min(batch_start + BATCH, len(gene_list))
    if (batch_start // BATCH) % 10 == 0:
        print(f"  {done:>6}/{len(gene_list):,} genes written", flush=True)

if writer:
    writer.close()

print(f"\nDone. contradicted rows written: {n_contradicted_rows:,}")

print("\nVerifying...")
led = pd.read_parquet(OUTPUTS / "evidence_ledger.parquet")
print(f"Shape: {led.shape}")
print(f"contradicted=True: {led['contradicted'].sum():,} rows "
      f"({led.loc[led['contradicted'], 'ensg_id'].nunique():,} distinct genes)")
print("\nlayers_present distribution:")
print(led["layers_present"].value_counts().to_string())
print("\nconfidence x contradicted cross-tab:")
print(pd.crosstab(led["confidence"], led["contradicted"]).to_string())
braf = led[led.ensg_id == "ensg00000157764"].head(3)
print("\nBRAF sample rows:")
print(braf.to_string(index=False))
print("Done.")
