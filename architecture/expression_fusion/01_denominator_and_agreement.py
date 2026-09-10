"""
expression_fusion/01_denominator_and_agreement.py
----------------------------------------------------
D1 (denominator) and D2 (detection floor), run for real against this repo's
warehouse, on a random sample of genes -- before any combination method is
chosen. Per the standing rule elsewhere in this repo: measure before you build.

A CORRECTION TO THE "BEFORE / AFTER" D1 AS FIRST PROPOSED
------------------------------------------------------------
The original framing was: "recompute percentiles restricted to the
intersection panel, re-measure cross-source Spearman rho before vs after; if
it jumps, the disagreement was never biological."

That comparison is a no-op under Spearman correlation specifically, and it is
worth being precise about why, because the error is easy to repeat elsewhere
in this project. Percentile-rank is a strictly monotonic function of the raw
value (ties aside). Spearman rho between two variables depends only on the
RELATIVE ORDER of the paired observations being compared -- and restricting a
full-panel percentile to a subset I preserves exactly the same relative order
among members of I as computing the percentile fresh on I would. So:

    rho( percentile_full(X)|_I , percentile_full(Y)|_I )
        ==  rho( percentile_I(X) , percentile_I(Y) )
        ==  rho( X|_I , Y|_I )                              (raw values)

"Before" and "after" are mathematically identical under Spearman. Measuring
them separately would either produce two identical numbers (if implemented
correctly) or two numbers that differ only due to a tie-handling or
implementation bug -- neither result answers the intended question.

What DOES change with the denominator, and is the operationally relevant
quantity, is the NUMERIC PERCENTILE VALUE for a given line -- because
percentile = rank / n, and n shrinks when the panel shrinks, even though
relative order is unchanged. That value is exactly what a tolerance-based
consumer reads (e.g. BUILD_SPEC.md's proposed GEO-corroboration rule,
"agrees within tau percentile units"). So this script measures two different,
non-redundant things:

  AGREEMENT   Spearman rho between sources on the intersection panel, per
              gene. One number (no before/after distinction needed) -- this
              is "how much do the sources actually agree, biologically."
  DENOMINATOR How much a line's percentile VALUE shifts, for the same line
              and the same source, between "percentile against that source's
              own full panel" (== today's E) and "percentile against the
              3-way intersection panel." This is the real before/after, and
              it is what a switch to a common-denominator design (M1) would
              change even though it would NOT change any Spearman-based
              agreement number.

Also included: a lightweight per-source floor check (D2), reusing
EXPRESSED_MIN / FLOOR_FRAC_LOW from common.py, since it falls out of the same
fetch almost for free.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

import common as C

N_GENES = 250

C.banner(f"D1/D2 -- denominator + agreement + floor, {N_GENES} genes, seed={C.SEED}")

con = C.connect()
genes = C.sample_genes(con, N_GENES)
print(f"sampled {len(genes)} genes from a universe of "
      f"{len(C.gene_universe(con)):,} genes present in all three sources\n")

long = C.fetch_all(con, genes.gene_id.tolist())
lineage = C.load_lineage(con)
print()

# =============================================================== D2: floor
C.banner("D2 -- per-source detection floor, across the sampled genes")
floor_rows = []
for src, g in long.groupby("source"):
    v = g["value"]
    frac_low = (v < C.EXPRESSED_MIN).mean()
    floor_rows.append({
        "source": src, "n_values": len(v),
        "min": round(float(v.min()), 3), "p01": round(float(v.quantile(.01)), 3),
        "frac_below_expressed_min": round(float(frac_low), 4),
        "floored": bool(v.min() > C.EXPRESSED_MIN and frac_low < C.FLOOR_FRAC_LOW),
    })
floor_df = pd.DataFrame(floor_rows)
print(floor_df.to_string(index=False))
floor_df.to_parquet(C.OUT / "d2_source_floor.parquet", index=False)

# =============================================================== D1: agreement + denominator
C.banner("D1 -- cross-source agreement and denominator effect, per gene")

wide_by_gene = {g: df.pivot(index="model_id", columns="source", values="value")
                for g, df in long.groupby("gene_id")}

agree_rows, denom_rows = [], []
rng = np.random.default_rng(C.SEED)

for gid, w in wide_by_gene.items():
    n_full = w.notna().sum()          # per-source panel size for this gene
    inter = w.dropna()                # 3-way intersection panel
    n_inter = len(inter)

    row = {"gene_id": gid,
           "n_depmap": int(n_full.get("depmap_expr", 0)),
           "n_hpa": int(n_full.get("hpa_rna", 0)),
           "n_geo": int(n_full.get("geo_expr", 0)),
           "n_intersection": n_inter}

    if n_inter >= 15:  # MIN_PEERS_FOR_LINEAGE-scale floor, below this rho is unstable
        for a, b in [("depmap_expr", "hpa_rna"), ("depmap_expr", "geo_expr"),
                     ("hpa_rna", "geo_expr")]:
            rho, p = stats.spearmanr(inter[a], inter[b])
            row[f"rho_{a}_{b}"] = round(float(rho), 4)
    agree_rows.append(row)

    # denominator effect: for each source, percentile against its OWN full
    # panel (today's E) vs percentile against the 3-way intersection only,
    # evaluated on the same lines (the intersection), same gene.
    if n_inter >= 15:
        for src in C.SOURCES:
            full_vals = w[src].dropna()
            pct_full_on_full = full_vals.rank(pct=True)
            pct_full_restricted = pct_full_on_full.loc[inter.index]
            pct_inter_only = inter[src].rank(pct=True)
            shift = (pct_full_restricted - pct_inter_only).abs()
            denom_rows.append({
                "gene_id": gid, "source": src,
                "n_full_panel": len(full_vals), "n_intersection": n_inter,
                "median_abs_pct_shift": float(shift.median()),
                "p90_abs_pct_shift": float(shift.quantile(0.90)),
                "max_abs_pct_shift": float(shift.max()),
            })

agree_df = pd.DataFrame(agree_rows)
denom_df = pd.DataFrame(denom_rows)
agree_df.to_parquet(C.OUT / "d1_agreement.parquet", index=False)
denom_df.to_parquet(C.OUT / "d1_denominator_shift.parquet", index=False)

usable = agree_df[agree_df.n_intersection >= 15]
print(f"genes with >=15 lines in the 3-way intersection: {len(usable)} of {len(agree_df)}\n")

print("panel sizes (median [IQR]):")
for col in ["n_depmap", "n_hpa", "n_geo", "n_intersection"]:
    s = agree_df[col]
    print(f"  {col:16s} {s.median():6.1f}  [{s.quantile(.25):.0f} - {s.quantile(.75):.0f}]")

print("\ncross-source agreement, Spearman rho on the intersection panel (median [IQR]):")
for pair in ["rho_depmap_expr_hpa_rna", "rho_depmap_expr_geo_expr", "rho_hpa_rna_geo_expr"]:
    s = usable[pair].dropna()
    if len(s):
        print(f"  {pair:28s} n_genes={len(s):4d}  median={s.median():.3f}  "
              f"IQR=[{s.quantile(.25):.3f} - {s.quantile(.75):.3f}]")

print("\ndenominator effect -- |percentile against own full panel - percentile against "
      "3-way intersection|, same line, same gene, same source (median [IQR] of the per-gene medians):")
for src, g in denom_df.groupby("source"):
    s = g["median_abs_pct_shift"]
    p90 = g["p90_abs_pct_shift"]
    print(f"  {src:12s} n_genes={len(g):4d}  median shift={s.median():.4f}  "
          f"IQR=[{s.quantile(.25):.4f} - {s.quantile(.75):.4f}]  "
          f"p90-of-p90={p90.quantile(.9):.4f}")

print(f"\nwrote {C.OUT / 'd1_agreement.parquet'}")
print(f"wrote {C.OUT / 'd1_denominator_shift.parquet'}")
print(f"wrote {C.OUT / 'd2_source_floor.parquet'}")
con.close()
