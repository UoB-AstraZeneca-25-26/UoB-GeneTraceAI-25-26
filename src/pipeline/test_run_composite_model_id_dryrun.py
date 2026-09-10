"""
Composite model_id: diagnosis, candidate fix, and DRY-RUN impact. Read-only.

Root cause
----------
`gdsc_models.broad_ach` is a delimiter-joined field IN THE GDSC SOURCE. The
resolver copied it verbatim into `model_id` and recorded `matched_via = 'ach'`,
so the value looks resolved and is not. `explode_model_id` never ran on this
table. `gdsc_models.rrid` carries the same defect ('cvcl_2717;cvcl_1888').

    sidm00096 sr      broad_ach = ach-000338;ach-000338            <- literal dup
    sidm01095 sjrh30  broad_ach = ach-000833;ach-001741;ach-001189
    sidm00400 sc-1    broad_ach = ach-002392;ach-002303   rrid = cvcl_2717;cvcl_1888

From there the composite propagated: gdsc_models -> procan_proteomics (resolved
via SIDM) -> core_score, where it appears as three protein-only phantom lines
carrying 15,736 rows across 6,612 genes.

Disposition is NOT uniform -- this is why Q1/Q2 had to be answered first:

  ach-000338;ach-000338            DELETE. ach-000338 already carries its own
                                   ProCan row and 5,270 two-layer rows. The
                                   phantom duplicates protein that is already
                                   correctly attributed.
  ach-000833;ach-001741;ach-001189 DELETE. ach-001189 (CVCL_0041, same line as
                                   ach-000833) already carries its own ProCan
                                   row. ach-001741 does not exist in DepMap.
  ach-002392;ach-002303            RENAME -> ach-002392. This phantom is the
                                   ONLY carrier: ach-002392 has 0 core_score
                                   rows and 0 ProCan rows today. ach-002303 and
                                   cvcl_2717 are absent from the DepMap roster,
                                   so ach-002392 (CVCL_1888, SC-1) is the sole
                                   resolvable target. Deleting would LOSE a
                                   real cell line; renaming recovers it.

So: 2 deletes, 1 rename. A uniform rule would be wrong either way -- delete-all
loses SC-1, rename-all duplicates two lines.
"""
from __future__ import annotations

import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ambiguity_detectors import assert_detectors_agree, composite_sql  # noqa: E402

OUT = Path(__file__).resolve().parent / "outputs"
con = duckdb.connect(str(OUT / "celllineselector.db"), read_only=True)
assert_detectors_agree(con)

DELETE = ["ach-000338;ach-000338", "ach-000833;ach-001741;ach-001189"]
RENAME = {"ach-002392;ach-002303": "ach-002392"}
ALL = DELETE + list(RENAME)


def head(t):
    print()
    print("=" * 74)
    print(t)
    print("=" * 74)


# ------------------------------------------------------------------ 1. spread
head("1 -- where the composite reaches")
rows = []
for t in [r[0] for r in con.execute("SHOW TABLES").fetchall()]:
    cols = {c[0]: str(c[1]) for c in con.execute(f"DESCRIBE {t}").fetchall()}
    for c, ty in cols.items():
        if "VARCHAR" not in ty:
            continue
        if c not in ("model_id", "broad_ach", "rrid"):
            continue
        n = con.execute(f"SELECT count(*) FROM {t} "
                        f"WHERE {composite_sql(chr(34)+c+chr(34))}").fetchone()[0]
        if n:
            rows.append({"table": t, "column": c, "rows": n})
spread = pd.DataFrame(rows).sort_values("rows", ascending=False)
print(spread.to_string(index=False))

# ------------------------------------------------- 2. what each phantom holds
head("2 -- what each phantom carries, and what its real counterpart has")
core = pd.read_parquet(OUT / "core_score.parquet")
for p in ALL:
    d = core[core.model_id == p]
    tgt = RENAME.get(p)
    print(f"\n  {p}")
    print(f"    core_score rows {len(d):>6,}   n_layers {d.n_layers.value_counts().to_dict()}"
          f"   genes {d.ensg_id.nunique():,}")
    for part in sorted(set(p.split(";"))):
        e = core[core.model_id == part]
        pc = con.execute("SELECT count(*) FROM procan_proteomics WHERE model_id=?",
                         [part]).fetchone()[0]
        exists = con.execute("SELECT count(*) FROM sample_info WHERE model_id=?",
                             [part]).fetchone()[0]
        print(f"      {part:<14s} in_roster={bool(exists)!s:<5s} core rows={len(e):>6,} "
              f"procan={pc}  n_layers={e.n_layers.value_counts().to_dict()}")
    print(f"    disposition: {'RENAME -> ' + tgt if tgt else 'DELETE'}")

# --------------------------------------------- 3. dry-run percentile impact
head("3 -- DRY RUN: ProCan protein percentile, before vs after")
# The phantoms are protein-only rows, so they inflate the ProCan percentile
# denominator. Rebuild the affected columns with the corrected line set.
affected = sorted(core[core.model_id.isin(ALL)].ensg_id.unique())
print(f"  genes touched by a phantom: {len(affected):,}")

g2u = con.execute("""
    SELECT gene_id, lower(u) AS u
    FROM (SELECT gene_id, unnest(uniprot_ids) AS u FROM gene_enriched
          WHERE uniprot_ids IS NOT NULL)
""").df()
procan_cols = {c[0] for c in con.execute("DESCRIBE procan_proteomics").fetchall()}
g2u = g2u[g2u.u.isin(procan_cols)]
usable = g2u[g2u.gene_id.isin(affected)].groupby("gene_id")["u"].apply(list)
print(f"  of those, measurable in ProCan: {len(usable):,}")

rng = np.random.default_rng(0)
sample = list(usable.index)
if len(sample) > 400:
    sample = [sample[i] for i in rng.choice(len(sample), 400, replace=False)]
print(f"  sample drawn: {len(sample)} genes")


def pct(s):
    s = pd.to_numeric(s, errors="coerce").dropna()
    return s.rank(pct=True, method="min")


recs = []
for gid in sample:
    accs = usable.loc[gid]
    expr = " + ".join(f'coalesce("{u}", 0)' for u in accs)
    cnt = " + ".join(f'CASE WHEN "{u}" IS NULL THEN 0 ELSE 1 END' for u in accs)
    d = con.execute(f'SELECT model_id, ({expr}) s, ({cnt}) n FROM procan_proteomics '
                    f'WHERE model_id IS NOT NULL').df()
    d = d[d["n"] > 0]
    if d.empty:
        continue
    before = (d["s"] / d["n"]).groupby(d["model_id"]).mean()

    # corrected line set: drop the deletes, relabel the rename
    d2 = d[~d.model_id.isin(DELETE)].copy()
    d2["model_id"] = d2["model_id"].replace(RENAME)
    after = (d2["s"] / d2["n"]).groupby(d2["model_id"]).mean()

    pb, pa = pct(before), pct(after)
    common = pb.index.intersection(pa.index)
    if len(common) == 0:
        continue
    delta = (pa[common] - pb[common]).abs()
    recs.append({"gene_id": gid, "n_before": len(pb), "n_after": len(pa),
                 "n_changed": int((delta > 1e-12).sum()),
                 "max_abs_dpct": float(delta.max()),
                 "median_abs_dpct": float(delta.median()),
                 "sc1_recovered": "ach-002392" in pa.index})

r = pd.DataFrame(recs)
print(f"\n  genes evaluated: {len(r)}")
print(f"  denominator n: before median {r.n_before.median():.0f}, "
      f"after median {r.n_after.median():.0f}  "
      f"(delta {r.n_after.median()-r.n_before.median():+.0f})")
print(f"  lines whose percentile moves, per gene: median {r.n_changed.median():.0f}, "
      f"max {r.n_changed.max()}")
print(f"  |delta pct| per gene: median-of-medians {r.median_abs_dpct.median():.5f}, "
      f"median-of-max {r.max_abs_dpct.median():.5f}, worst {r.max_abs_dpct.max():.5f}")
print(f"  genes where the DKW band (+/-0.044) is exceeded: "
      f"{int((r.max_abs_dpct > 0.044).sum())} of {len(r)}")
print(f"  genes where SC-1 (ach-002392) becomes scoreable: "
      f"{int(r.sc1_recovered.sum())} of {len(r)}")

# ------------------------------------------------------ 4. net row accounting
head("4 -- net effect on core_score rows")
del_rows = int(core.model_id.isin(DELETE).sum())
ren_rows = int(core.model_id.isin(RENAME).sum())
print(f"  duplicate protein-only rows removed : {del_rows:,}")
print(f"  rows relabelled to ach-002392       : {ren_rows:,}")
print(f"  net change                          : {-del_rows:+,}")
print(f"  cell lines recovered                : 1 (SC-1 / ach-002392, "
      f"currently 0 core_score rows)")
print(f"  phantom lines eliminated            : 3")

r.to_parquet(OUT / "test_run_composite_model_id_dryrun.parquet", index=False)
print(f"\nper-gene detail -> "
      f"{OUT / 'test_run_composite_model_id_dryrun.parquet'}")
print("\nNOTHING WAS WRITTEN to the warehouse or to core_score.parquet.")
con.close()
