import pandas as pd, numpy as np, json

core  = pd.read_parquet("src/pipeline/outputs/core_score.parquet",
                        columns=["model_id","ensg_id","core_score","n_layers"])
regime = pd.read_parquet("src/pipeline/outputs/gene_regime.parquet")
flags  = pd.read_parquet("src/pipeline/outputs/flags_with_driver.parquet",
                         columns=["model_id","ensg_id","has_driver_alteration","has_alteration"])
gdsc   = pd.read_parquet("validation/prepared/gdsc_scored_ready.parquet",
                         columns=["model_id","target_ensg","sensitive"])

core["model_id"]  = core["model_id"].str.lower()
flags["model_id"] = flags["model_id"].str.lower()
gdsc["model_id"]  = gdsc["model_id"].str.lower()
gdsc["target_ensg"] = gdsc["target_ensg"].astype("string").str.lower()

# Sorted before sampling — deterministic across Python runs
curated_genes = sorted(set(regime.ensg_id) & set(gdsc.target_ensg))
rng = np.random.default_rng(42)
test_genes  = set(rng.choice(curated_genes, size=max(1, len(curated_genes)//5), replace=False))
train_genes = [g for g in curated_genes if g not in test_genes]

print(f"curated: {len(curated_genes)}, test: {len(test_genes)}, train: {len(train_genes)}")
print("BCL2 in test:", "ensg00000171791" in test_genes)
print("BRAF in test:", "ensg00000157764" in test_genes)

core_t  = core[core.ensg_id.isin(test_genes)].copy()
flags_t = flags[flags.ensg_id.isin(test_genes)].copy()
gdsc_t  = (gdsc[gdsc.target_ensg.isin(test_genes) & gdsc.sensitive]
           .rename(columns={"target_ensg": "ensg_id"}))

core_t = core_t.merge(regime[["ensg_id","class"]], on="ensg_id", how="left")
core_t = core_t.merge(
    flags_t[["model_id","ensg_id","has_driver_alteration","has_alteration"]],
    on=["model_id","ensg_id"], how="left")
core_t["has_driver_alteration"] = core_t["has_driver_alteration"].fillna(False)
core_t["has_alteration"]        = core_t["has_alteration"].fillna(False)

sens_pairs = gdsc_t[["ensg_id","model_id"]].drop_duplicates().copy()
sens_pairs["is_sensitive"] = True
core_t = core_t.merge(sens_pairs, on=["ensg_id","model_id"], how="left")
core_t["is_sensitive"] = core_t["is_sensitive"].fillna(False)

core_t["rank_flat"] = core_t.groupby("ensg_id")["core_score"].rank(ascending=False, method="first")

core_t["sort_driver"] = np.where(
    core_t["class"] == "abundance_tracking", core_t["core_score"],
    core_t["has_driver_alteration"].astype(float) * 1e6 + core_t["core_score"])
core_t["rank_driver"] = core_t.groupby("ensg_id")["sort_driver"].rank(ascending=False, method="first")

core_t["sort_any"] = np.where(
    core_t["class"] == "abundance_tracking", core_t["core_score"],
    core_t["has_alteration"].astype(float) * 1e6 + core_t["core_score"])
core_t["rank_any"] = core_t.groupby("ensg_id")["sort_any"].rank(ascending=False, method="first")

known = core_t.groupby("ensg_id")["is_sensitive"].sum()
per_gene = pd.DataFrame({
    "hits_flat":   core_t[core_t.rank_flat   <= 20].groupby("ensg_id")["is_sensitive"].sum(),
    "hits_driver": core_t[core_t.rank_driver <= 20].groupby("ensg_id")["is_sensitive"].sum(),
    "hits_any":    core_t[core_t.rank_any    <= 20].groupby("ensg_id")["is_sensitive"].sum(),
    "known": known,
})
per_gene = per_gene[per_gene.known > 0].merge(regime[["ensg_id","class"]], on="ensg_id")
per_gene["hr_flat"]      = per_gene.hits_flat   / per_gene.known
per_gene["hr_driver"]    = per_gene.hits_driver / per_gene.known
per_gene["hr_any"]       = per_gene.hits_any    / per_gene.known
per_gene["delta_driver"] = per_gene.hr_driver - per_gene.hr_flat
per_gene["delta_any"]    = per_gene.hr_any    - per_gene.hr_flat

activ = per_gene[per_gene["class"] == "activation_driven"]
abund = per_gene[per_gene["class"] == "abundance_tracking"]

print()
print("=== Per-class means (denominator = known) ===")
print(per_gene.groupby("class")[["hr_flat","hr_any","hr_driver","delta_driver","delta_any"]]
      .mean().round(4).to_string())
print()
print("Variance reduction (activation_driven delta std):")
print(f"  any-alt:  {activ.delta_any.std():.4f}")
print(f"  driver:   {activ.delta_driver.std():.4f}")
print(f"  pct:      {(1-activ.delta_driver.std()/activ.delta_any.std())*100:.1f}%")
print()

# Spotlight: BCL2 and BRAF — evaluate from full data regardless of split
for name, ensg in [("BCL2","ensg00000171791"), ("BRAF","ensg00000157764")]:
    in_test = ensg in test_genes
    gdsc_g  = gdsc[gdsc.target_ensg == ensg]
    core_g  = core[core.ensg_id == ensg]
    flags_g = flags[flags.ensg_id == ensg]
    sens    = set(gdsc_g[gdsc_g.sensitive].model_id)
    kn      = sens & set(core_g.model_id)
    rc      = regime.set_index("ensg_id").loc[ensg, "class"]
    m = core_g.merge(flags_g[["model_id","has_driver_alteration"]], on="model_id", how="left")
    m["has_driver_alteration"] = m["has_driver_alteration"].fillna(False)
    if rc == "abundance_tracking":
        m = m.sort_values("core_score", ascending=False)
    else:
        m = m.sort_values(["has_driver_alteration","core_score"], ascending=[False,False])
    flat_hr = len(set(core_g.sort_values("core_score",ascending=False).head(20).model_id) & sens) / len(kn)
    dr_hr   = len(set(m.head(20).model_id) & sens) / len(kn)
    print(f"{name} ({rc}, {'test' if in_test else 'train'}): "
          f"flat={flat_hr:.4f}  driver={dr_hr:.4f}  delta={dr_hr-flat_hr:+.4f}  known={len(kn)}")

per_gene.to_parquet("src/pipeline/outputs/stage4_eval_consistent_denom.parquet", index=True)

summary = {
    "denominator": "known = sensitive_in_GDSC intersect scored_in_core_score",
    "split_note": "curated_genes sorted() before rng.choice(seed=42) — reproducible",
    "seed": 42,
    "n_curated": len(curated_genes),
    "n_test": len(test_genes),
    "n_train": len(train_genes),
    "test_genes": sorted(list(test_genes)),
    "per_class": {
        "abundance_tracking": {
            "hr_flat":         round(abund.hr_flat.mean(), 6),
            "hr_any_alt":      round(abund.hr_any.mean(), 6),
            "hr_driver":       round(abund.hr_driver.mean(), 6),
            "delta_driver":    round(abund.delta_driver.mean(), 6),
            "delta_driver_std":round(float(abund.delta_driver.std()), 6),
            "n_genes": int(len(abund)),
        },
        "activation_driven": {
            "hr_flat":         round(activ.hr_flat.mean(), 6),
            "hr_any_alt":      round(activ.hr_any.mean(), 6),
            "hr_driver":       round(activ.hr_driver.mean(), 6),
            "delta_driver":    round(activ.delta_driver.mean(), 6),
            "delta_driver_std":round(activ.delta_driver.std(), 6),
            "delta_any_std":   round(activ.delta_any.std(), 6),
            "n_genes": int(len(activ)),
        },
    },
    "overall": {
        "hr_flat":         round(per_gene.hr_flat.mean(), 6),
        "hr_driver":       round(per_gene.hr_driver.mean(), 6),
        "hr_any_alt":      round(per_gene.hr_any.mean(), 6),
        "driver_lift_pct": round(
            (per_gene.hr_driver.mean()-per_gene.hr_flat.mean())/per_gene.hr_flat.mean()*100, 2),
    },
    "variance_reduction": {
        "delta_std_any_alt": round(activ.delta_any.std(), 6),
        "delta_std_driver":  round(activ.delta_driver.std(), 6),
        "reduction_pct":     round((1-activ.delta_driver.std()/activ.delta_any.std())*100, 1),
    },
    "spotlight": {
<<<<<<< HEAD
        "BCL2_ENSG00000171791": {"class": "abundance_tracking", "in_test_set": "ensg00000171791" in test_genes},
        "BRAF_ENSG00000157764": {"class": "activation_driven",  "in_test_set": "ensg00000157764" in test_genes},
=======
        "BCL2_ENSG00000171791": {"class": "abundance_tracking", "in_test_set": "ENSG00000171791" in test_genes},
        "BRAF_ENSG00000157764": {"class": "activation_driven",  "in_test_set": "ENSG00000157764" in test_genes},
>>>>>>> 20fd4e72bc1fd7b405f50772ccc06f338c5a0bb4
    },
}
with open("src/pipeline/outputs/stage4_eval_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print()
print("Saved: stage4_eval_consistent_denom.parquet")
print("Saved: stage4_eval_summary.json")
