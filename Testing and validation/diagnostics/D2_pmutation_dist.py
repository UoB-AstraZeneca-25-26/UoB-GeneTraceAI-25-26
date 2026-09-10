"""p_mutation distribution and ceiling fraction (Decision 5)."""
import sys
import numpy as np
import pandas as pd
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
MUT_SRC = REPO / "cleaned_track_data/mutations_collapsed.parquet"

P0_VEP=2.5; K_VEP=2.0; P0_PATH=0.5; K_PATH=2.0; P0_BURDEN=3.0; K_BURDEN=1.5
DRIVER_BOOST=0.15; MUT_GATE=0.5

def hill(x, p0, k):
    xk = np.power(np.clip(x, 0, None), k)
    return xk / (xk + p0**k)

mut = pd.read_parquet(MUT_SRC)
mut = mut.copy()
mut["variant_burden"] = mut["variant_count"].clip(upper=5)
mut["p_vep"]    = hill(mut["max_vep_rank"].astype(float), P0_VEP, K_VEP)
mut["p_path"]   = hill(mut["max_pathogenicity"].astype(float), P0_PATH, K_PATH)
mut["p_burden"] = hill(mut["variant_burden"].astype(float), P0_BURDEN, K_BURDEN)
p_base = (1.0 - (1.0-mut["p_vep"].fillna(0)) *
                (1.0-mut["p_path"].fillna(0)) *
                (1.0-mut["p_burden"].fillna(0)))
driver_flag = mut["any_driver"].fillna(False).astype(int)
mut["p_mutation"] = (p_base + DRIVER_BOOST * driver_flag).clip(upper=1.0)
mut["mut_driver"] = mut["p_mutation"] >= MUT_GATE

print("=== ALL RECORDS ===")
print(mut["p_mutation"].describe().round(4))
print(f"At ceiling (p_mutation=1.0): {(mut.p_mutation==1.0).mean()*100:.2f}%")
print(f"p_mutation >= 0.95:          {(mut.p_mutation>=0.95).mean()*100:.2f}%")
print(f"p_mutation >= 0.85:          {(mut.p_mutation>=0.85).mean()*100:.2f}%")
print(f"p_mutation >= 0.5 (driver):  {mut.mut_driver.mean()*100:.2f}%")
print(f"p_base >= 0.85 (boost saturates): {(p_base>=0.85).mean()*100:.2f}%")
print()

print("=== p_base DISTRIBUTION (before driver boost) ===")
print(p_base.describe().round(4))
print()

print("=== ROLE-ANNOTATED GENES ONLY (oncogene/tsg/both) ===")
gl = pd.read_parquet(REPO / "reference/gene_lookup.parquet", columns=["ensg_id","gene_role"])
mut2 = mut.merge(gl, on="ensg_id", how="left")
role_mut = mut2[mut2.gene_role.isin(["oncogene","tsg","both"])]
print(f"n records: {len(role_mut):,}  ({role_mut.ensg_id.nunique()} genes)")
print(f"mut_driver rate: {role_mut.mut_driver.mean()*100:.2f}%")
print(f"p_mutation at ceiling: {(role_mut.p_mutation==1.0).mean()*100:.2f}%")
print()

print("=== CONTRIBUTION PER CHANNEL ===")
mut3 = mut.copy()
mut3["only_vep"]     = (mut3["p_vep"] >= 0.5) & (mut3["p_path"] < 0.5) & (mut3["p_burden"] < 0.5)
mut3["only_path"]    = (mut3["p_path"] >= 0.5) & (mut3["p_vep"] < 0.5) & (mut3["p_burden"] < 0.5)
mut3["only_burden"]  = (mut3["p_burden"] >= 0.5) & (mut3["p_vep"] < 0.5) & (mut3["p_path"] < 0.5)
mut3["multi_chan"]   = ((mut3["p_vep"] >= 0.5).astype(int) +
                        (mut3["p_path"] >= 0.5).astype(int) +
                        (mut3["p_burden"] >= 0.5).astype(int)) >= 2
passing = mut3[mut3["mut_driver"]]
print(f"Among mut_driver=TRUE records ({len(passing):,}):")
print(f"  Only p_vep >= 0.5:    {passing.only_vep.mean()*100:.1f}%")
print(f"  Only p_path >= 0.5:   {passing.only_path.mean()*100:.1f}%")
print(f"  Only p_burden >= 0.5: {passing.only_burden.mean()*100:.1f}%")
print(f"  Multi-channel:        {passing.multi_chan.mean()*100:.1f}%")
