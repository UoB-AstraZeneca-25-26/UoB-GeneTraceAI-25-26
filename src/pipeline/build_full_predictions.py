"""
Extends predictions_with_confidence.parquet from 141 curated genes
to all 19,176 scored genes. Reads core_score and flags in gene-batches
to keep memory low, writes parquet via pyarrow in streaming fashion.
"""
import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path

OUTPUTS = Path("src/pipeline/outputs")

# Load lightweight tables in full
regime = pd.read_parquet(OUTPUTS / "gene_regime.parquet")
regime_map = regime.set_index("ensg_id")[["class","regime_source"]].to_dict("index")

# Gene list only from core (no data yet)
gene_list = pd.read_parquet(
    OUTPUTS / "core_score.parquet", columns=["ensg_id"]
)["ensg_id"].unique()
print(f"{len(gene_list):,} genes to process")

OUT_COLS = ["model_id","ensg_id","core_score","n_layers","stratum_rank",
            "class","regime_source","rank_basis","has_driver_alteration",
            "confidence","confidence_reason"]

writer = None
BATCH = 200  # genes per batch

for batch_start in range(0, len(gene_list), BATCH):
    batch_genes = gene_list[batch_start:batch_start+BATCH]

    # Read only this batch's rows from both parquets
    core_b = pd.read_parquet(OUTPUTS / "core_score.parquet").query("ensg_id in @batch_genes")
    core_b["model_id"] = core_b["model_id"].str.upper()

    flags_b = pd.read_parquet(
        OUTPUTS / "flags_with_driver.parquet",
        columns=["model_id","ensg_id","has_driver_alteration"]
    ).query("ensg_id in @batch_genes")
    flags_b["model_id"] = flags_b["model_id"].str.upper()

    c = core_b.merge(flags_b, on=["model_id","ensg_id"], how="left")
    c["has_driver_alteration"] = c["has_driver_alteration"].fillna(False)

    c["class"]         = c["ensg_id"].map(lambda g: regime_map.get(g, {}).get("class", "unknown"))
    c["regime_source"] = c["ensg_id"].map(lambda g: regime_map.get(g, {}).get("regime_source", "unknown"))

    c["rank_basis"] = np.where(
        c["class"] == "abundance_tracking", "core_score",
        np.where(c["has_driver_alteration"], "driver_flag+score", "score_only"))

    tier = np.full(len(c), "moderate", dtype=object)
    tier = np.where(
        (c.regime_source == "measured") & (c["class"] == "abundance_tracking") & (c.n_layers == 2),
        "high", tier)
    tier = np.where(
        (c["class"] == "activation_driven") & (c.rank_basis == "score_only"),
        "low", tier)
    tier = np.where(c.regime_source != "measured", "unknown", tier)
    c["confidence"] = tier

    c["confidence_reason"] = np.where(
        c.confidence == "high",
        "High: abundance_tracking gene (empirically validated), scored on "
            + c.n_layers.astype(str) + " layers.",
        np.where(
            c.confidence == "low",
            "Low: activation-driven gene dependency; ranked on score alone "
            "with no supporting driver alteration in this line.",
            np.where(
                c.confidence == "unknown",
                "Unknown: no validation data for this gene; scoring regime "
                "unclassified, ranking shown but not regime-weighted.",
                np.where(
                    c.rank_basis == "driver_flag+score",
                    "Moderate: activation-driven gene; ranked on driver alteration "
                    "present in this line, with score as support.",
                    "Moderate: " + c["class"] + " gene, ranked on "
                        + c.rank_basis + ", " + c.n_layers.astype(str) + "-layer coverage."
                )
            )
        )
    )

    table = pa.Table.from_pandas(c[OUT_COLS], preserve_index=False)
    if writer is None:
        writer = pq.ParquetWriter(OUTPUTS / "predictions_with_confidence.parquet", table.schema)
    writer.write_table(table)

    done = min(batch_start + BATCH, len(gene_list))
    if (batch_start // BATCH) % 5 == 0:
        print(f"  {done:>6}/{len(gene_list):,} genes written", flush=True)

if writer:
    writer.close()

print("Verifying...")
pred = pd.read_parquet(OUTPUTS / "predictions_with_confidence.parquet")
print(f"Shape: {pred.shape}  |  genes: {pred.ensg_id.nunique():,}")
print("Confidence dist:", pred.confidence.value_counts().to_dict())
shank3 = pred[pred.ensg_id == "ENSG00000251322"]
print(f"SHANK3 rows: {len(shank3)}")
if not shank3.empty:
    print(shank3[["model_id","core_score","class","confidence","rank_basis"]].head(3).to_string(index=False))
print("Done.")
