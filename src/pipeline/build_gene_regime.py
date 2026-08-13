"""
build_gene_regime.py
--------------------
Writes src/pipeline/outputs/gene_regime.parquet — the missing producer.

    python src/pipeline/build_gene_regime.py            # build
    python src/pipeline/build_gene_regime.py --compare  # build + diff vs the archived file

WHY THIS EXISTS
    gene_regime.parquet had no producer anywhere in the repository. It was a
    static file of untraceable provenance, read by notebooks 05 and 06,
    build_full_predictions.py, build_cna_layer.py, stage4_eval_save.py, two
    test_run_* diagnostics and both figure scripts. A clean rebuild could not
    regenerate it, so the pipeline could not be re-run past Stage 4, and
    src/figures/compute_stats.py could not run at all.

    Worse, the archived copy became actively misleading once Stage 2 was re-run:
    its regime labels were measured against a core_score that no longer exists.
    Restoring it would have paired new scores with old labels, silently.

WHAT A REGIME IS
    regime_class records how a gene's biology relates to abundance:

      abundance_tracking  more protein means more drug target — abundance alone
                          should predict response, so rank on core_score
      activation_driven   what matters is whether the gene is ACTIVATED (mutation
                          or fusion), not how much of it there is, so rank on the
                          driver flag first and core_score second

    The class is what selects the routing strategy in Stage 4 and gates the
    confidence tier in Stage 5, which is why a stale label is not cosmetic.

HOW THE CLASS IS DERIVED
    For each curated gene, rank every scored cell line by core_score alone
    (the "flat" ranking — no driver flag), and ask how many of the lines GDSC
    calls drug-sensitive land in the top k:

        hit_at_20_train = |{sensitive lines in top-k by core_score}|
                          / |{sensitive lines present in core_score}|

    i.e. recall@k over the SCORED universe. A gene whose abundance ranking
    recovers sensitive lines better than chance is abundance-tracking; a gene
    where abundance alone tells you little must be driven by activation instead.

    The cut is expressed as a multiple of the random baseline rather than as an
    absolute number. Picking k lines at random out of the n scored for a gene
    recovers a fraction k/n of the sensitive ones in expectation, so:

        baseline_g   = k / n_scored_g
        class        = abundance_tracking  if hit_at_20 >= LIFT * baseline_g
                       activation_driven   otherwise

  ┌─ RECONSTRUCTED, NOT RECOVERED ────────────────────────────────────────────┐
  │ The original rule was never written down — not in MATH_REFERENCE.md, not  │
  │ GLOSSARY.md, not DECISIONS.md, and not in any script. What IS recoverable │
  │ from discarded/outputs/gene_regime.parquet is that the two classes are    │
  │ PERFECTLY separated by hit_at_20_train:                                   │
  │                                                                           │
  │     activation_driven   max = 0.026667                                    │
  │     abundance_tracking  min = 0.027027                                    │
  │                                                                           │
  │ so the class was certainly a threshold on that column, and the cut lay in │
  │ (0.026667, 0.027027). Two candidates fall inside that gap and cannot be   │
  │ told apart from the data alone:                                          │
  │                                                                           │
  │     2 x 20/1485 = 0.026936   (twice random, on the old 1,485-line universe)│
  │     0.027                    (a round hand-picked number)                 │
  │                                                                           │
  │ LIFT = 2.0 over a per-gene baseline is implemented because it reproduces  │
  │ the old boundary AND stays meaningful when coverage changes — the scored  │
  │ universe has already moved from 1,485 to 1,479-1,746 lines per gene, and  │
  │ a hard-coded 0.027 would quietly mean something different on it.          │
  │                                                                           │
  │ Both are settable: --lift, or --threshold for an absolute cut.            │
  └───────────────────────────────────────────────────────────────────────────┘

  ┌─ KNOWN METHODOLOGICAL WEAKNESS (inherited, not introduced) ───────────────┐
  │ The class is fitted on ALL curated genes using their own GDSC sensitive   │
  │ lines. stage4_eval_save.py and 06_held_out_eval.ipynb then hold out 20%   │
  │ of those same genes and evaluate on them — but those genes' labels were   │
  │ derived from the data being tested. That is label leakage, and it will    │
  │ flatter the held-out numbers.                                            │
  │                                                                           │
  │ The column name `hit_at_20_train` suggests the original author knew the   │
  │ measurement was in-sample. This script reproduces that behaviour faithfully│
  │ rather than silently changing the methodology; fixing it means assigning   │
  │ regimes inside the training fold only, which changes what Stage 5/6 mean  │
  │ and is a decision for the author, not this rebuild.                       │
  └───────────────────────────────────────────────────────────────────────────┘

OUTPUT SCHEMA
    The five original columns, in their original order and dtypes, so every
    existing reader works unchanged:

        ensg_id, hit_at_20_train, class, regime_source, n_sensitive_gdsc

    plus diagnostic columns the original lacked. Every downstream reader selects
    columns explicitly, so these are additive and safe:

        n_scored, n_known, hits_at_k, baseline_at_k, lift_over_random
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / "src/pipeline/outputs/core_score.parquet"
GDSC = ROOT / "validation/prepared/gdsc_scored_ready.parquet"
OUT = ROOT / "src/pipeline/outputs/gene_regime.parquet"
ARCHIVED = ROOT / "discarded/outputs/gene_regime.parquet"

K_DEFAULT = 20
LIFT_DEFAULT = 2.0


def load_inputs():
    if not CORE.is_file():
        raise SystemExit(f"{CORE} not found — run 02_core_score.ipynb first.")
    if not GDSC.is_file():
        raise SystemExit(f"{GDSC} not found — this is the GDSC ground truth.")

    core = pd.read_parquet(CORE, columns=["model_id", "ensg_id", "core_score"])
    gdsc = pd.read_parquet(GDSC, columns=["model_id", "target_ensg", "sensitive"])

    # canonical key case is lowercase, matching the harmonisation warehouse
    for df, cols in [(core, ("model_id", "ensg_id")),
                     (gdsc, ("model_id", "target_ensg"))]:
        for c in cols:
            df[c] = df[c].astype("string").str.lower()
    return core, gdsc


def build(core, gdsc, k=K_DEFAULT, lift=LIFT_DEFAULT, threshold=None):
    # The curated set is every gene GDSC has a drug target for. Verified against
    # the archived file: exactly the same 141 genes, no minimum-evidence cut.
    sens = (gdsc[gdsc["sensitive"]]
            .drop_duplicates(["target_ensg", "model_id"])
            .rename(columns={"target_ensg": "ensg_id"})[["ensg_id", "model_id"]])
    n_sens = sens.groupby("ensg_id").size().rename("n_sensitive_gdsc")
    curated = sorted(n_sens.index)
    print(f"curated genes (GDSC targets with >=1 sensitive line): {len(curated):,}")

    c = core[core["ensg_id"].isin(set(curated))].copy()
    c = c.merge(sens.assign(is_sensitive=True), on=["ensg_id", "model_id"], how="left")
    c["is_sensitive"] = c["is_sensitive"].astype("boolean").fillna(False).astype(bool)

    # "first" matches stage4_eval_save.py, so the two agree on tie handling
    c["rank_flat"] = c.groupby("ensg_id")["core_score"].rank(ascending=False, method="first")

    n_scored = c.groupby("ensg_id").size().rename("n_scored")
    # denominator is sensitive AND scored -- a sensitive line the pipeline never
    # scored cannot be recovered by any ranking, so counting it would penalise
    # the gene for a coverage gap rather than for a bad ranking
    n_known = c.groupby("ensg_id")["is_sensitive"].sum().rename("n_known")
    hits = (c[c["rank_flat"] <= k].groupby("ensg_id")["is_sensitive"].sum()
            .rename("hits_at_k"))

    g = pd.concat([n_scored, n_known, hits], axis=1).reindex(curated)
    g[["n_known", "hits_at_k"]] = g[["n_known", "hits_at_k"]].fillna(0)
    g["n_scored"] = g["n_scored"].fillna(0)

    g["hit_at_20_train"] = (g["hits_at_k"] / g["n_known"].replace(0, np.nan)).fillna(0.0)
    g["baseline_at_k"] = k / g["n_scored"].replace(0, np.nan)
    g["lift_over_random"] = g["hit_at_20_train"] / g["baseline_at_k"]

    if threshold is not None:
        cut = pd.Series(float(threshold), index=g.index)
        rule = f"absolute hit@{k} >= {threshold}"
    else:
        cut = lift * g["baseline_at_k"]
        rule = f"hit@{k} >= {lift} x random baseline (k / n_scored, per gene)"

    g["class"] = np.where(g["hit_at_20_train"] >= cut,
                          "abundance_tracking", "activation_driven")
    g["regime_source"] = "measured"
    g = g.join(n_sens)

    out = (g.reset_index().rename(columns={"index": "ensg_id"})
           [["ensg_id", "hit_at_20_train", "class", "regime_source", "n_sensitive_gdsc",
             "n_scored", "n_known", "hits_at_k", "baseline_at_k", "lift_over_random"]]
           .sort_values("ensg_id").reset_index(drop=True))
    out["n_sensitive_gdsc"] = out["n_sensitive_gdsc"].astype(int)
    for col in ("n_scored", "n_known", "hits_at_k"):
        out[col] = out[col].astype(int)

    print(f"rule: {rule}")
    print(f"\nclass counts:\n{out['class'].value_counts().to_string()}")
    print(f"\nhit@{k} by class:")
    print(out.groupby("class")["hit_at_20_train"]
          .describe()[["count", "min", "mean", "max"]].to_string())
    print(f"\ngenes scoring exactly 0 at k={k}: {int((out.hit_at_20_train == 0).sum())}")
    print(f"lines scored per gene: {out.n_scored.min():,}-{out.n_scored.max():,}")
    return out


def compare_to_archived(new):
    if not ARCHIVED.is_file():
        print(f"\n[compare] {ARCHIVED} not present — skipping")
        return
    old = pd.read_parquet(ARCHIVED)
    old["ensg_id"] = old["ensg_id"].astype("string").str.lower()
    m = (old[["ensg_id", "hit_at_20_train", "class", "n_sensitive_gdsc"]]
         .rename(columns={"hit_at_20_train": "hit_old", "class": "class_old",
                          "n_sensitive_gdsc": "n_sens_old"})
         .merge(new[["ensg_id", "hit_at_20_train", "class", "n_sensitive_gdsc"]]
                .rename(columns={"hit_at_20_train": "hit_new", "class": "class_new",
                                 "n_sensitive_gdsc": "n_sens_new"}),
                on="ensg_id", how="outer", indicator=True))

    print("\n" + "=" * 68)
    print("COMPARISON WITH THE ARCHIVED FILE")
    print("=" * 68)
    print(f"genes: archived {len(old)}, rebuilt {len(new)}, "
          f"in both {(m._merge == 'both').sum()}")
    both = m[m._merge == "both"]
    print(f"n_sensitive_gdsc identical : {int((both.n_sens_old == both.n_sens_new).sum())}"
          f"/{len(both)}  <- confirms the gene set and the sensitivity count")
    agree = int((both.class_old == both.class_new).sum())
    print(f"class agreement            : {agree}/{len(both)} ({agree / len(both) * 100:.1f}%)")
    print(f"hit@20 Spearman            : "
          f"{both[['hit_old', 'hit_new']].corr(method='spearman').iloc[0, 1]:.3f}")
    print("\nconfusion (rows = archived, cols = rebuilt):")
    print(pd.crosstab(both.class_old, both.class_new,
                      margins=True, margins_name="TOTAL").to_string())
    print("\nDisagreement is EXPECTED, not a defect: the archived labels were measured\n"
          "against the pre-v2 core_score, which had 28.4M rows over 1,485 lines and no\n"
          "ProCan proteomics. These are measured against 29.8M rows over 1,746 lines.\n"
          "That re-measurement is the entire reason for rebuilding rather than restoring.")


def main():
    ap = argparse.ArgumentParser(description="Build gene_regime.parquet")
    ap.add_argument("--k", type=int, default=K_DEFAULT, help="top-k for the hit rate")
    ap.add_argument("--lift", type=float, default=LIFT_DEFAULT,
                    help="multiple of the random baseline that marks abundance_tracking")
    ap.add_argument("--threshold", type=float, default=None,
                    help="absolute hit-rate cut, overriding --lift (e.g. 0.027)")
    ap.add_argument("--compare", action="store_true",
                    help="diff the result against discarded/outputs/gene_regime.parquet")
    ap.add_argument("--dry-run", action="store_true", help="do not write the parquet")
    args = ap.parse_args()

    core, gdsc = load_inputs()
    out = build(core, gdsc, k=args.k, lift=args.lift, threshold=args.threshold)

    if args.compare:
        compare_to_archived(out)

    if args.dry_run:
        print("\n[dry run] not written")
        return
    out.to_parquet(OUT, index=False)
    print(f"\nWritten : {OUT}")
    print(f"Shape   : {out.shape}")
    print(f"Columns : {out.columns.tolist()}")


if __name__ == "__main__":
    main()
