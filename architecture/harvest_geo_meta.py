"""
harvest_geo_meta.py
--------------------
Builds geo_meta/*.txt (per-GSM SOFT metadata) from NCBI GEO, needed by
the new transcriptomics pipeline's --build-geo-platform step to split
GEO samples into platform x processing-family methods.

`harvest_geo_methods()` below is copied verbatim from Triveni's
05_transcriptomics_stale.py (origin/triveni, Cell Line Selector/pipelines/),
minus the two src.scripts imports (get_logger, log_error) that aren't
needed for this function specifically -- kept here as a standalone
script rather than pulling in her whole package, per instruction to
leave the folder structure alone for now.

Usage:
    python harvest_geo_meta.py --db outputs/celllineselector.db --destdir geo_meta
    python harvest_geo_meta.py --db outputs/celllineselector.db --destdir geo_meta --limit 5   # smoke test
"""
from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

try:
    import GEOparse  # type: ignore[import-not-found]
except ImportError:
    GEOparse = None


def harvest_geo_methods(gsms: list, destdir: str = "./geo_meta", verbose: bool = True) -> pd.DataFrame:
    """
    Fetch per-GSM metadata from NCBI and classify the normalisation
    method. GEO archives; it does not compute -- each GSM carries
    whatever its submitter uploaded, so the method must be read per
    sample rather than assumed for a whole series.

    Checks destdir for an already-downloaded file per GSM first, so
    reruns only hit the network for genuinely new accessions.
    """
    if GEOparse is None:
        raise ImportError("GEOparse is required. Install it with: pip install GEOparse")

    os.makedirs(destdir, exist_ok=True)

    rows = []
    n_cached = 0
    for i, g in enumerate(gsms, 1):
        try:
            existing = glob.glob(os.path.join(destdir, f"{g.upper()}*"))
            if existing:
                s = GEOparse.get_GEO(filepath=existing[0], silent=True)
                n_cached += 1
            else:
                s = GEOparse.get_GEO(geo=g.upper(), destdir=destdir, silent=True)
            m = s.metadata
            rows.append({
                "gsm": g.lower(),
                "platform": m.get("platform_id", [None])[0],
                "processing": (m.get("data_processing", [""])[0] or "")[:200],
                "value_def": next((v for k, v in s.columns.description.items()
                                    if k.upper() == "VALUE"), None),
                "has_abs_call": "ABS_CALL" in s.table.columns,
                "series": "; ".join(m.get("series_id", [])),
                "title": (m.get("title", [""])[0] or "")[:120],
            })
        except Exception as e:
            rows.append({"gsm": g.lower(), "platform": None,
                         "processing": f"FETCH FAILED: {e}"})
        if verbose and i % 200 == 0:
            print(f"  {i:,} / {len(gsms):,}  ({n_cached:,} from local cache)")

    meta = pd.DataFrame(rows)

    # gcRMA MUST be tested before RMA -- the string "gcRMA" contains "RMA"
    meta["method"] = np.select(
        [meta.processing.str.contains("gcrma", case=False, na=False),
         meta.processing.str.contains("rma|robust multi", case=False, na=False),
         meta.processing.str.contains("mas5|mas 5|gcos", case=False, na=False),
         meta.processing.str.contains("quantile", case=False, na=False)],
        ["gcRMA", "RMA", "MAS5", "quantile"], default="other")

    if verbose:
        print(f"\n{len(meta):,} GSMs ({n_cached:,} loaded locally)\n")
        print(meta.platform.value_counts().to_string())
        print()
        print(meta.method.value_counts().to_string())
        failed = meta.processing.str.startswith("FETCH FAILED").sum()
        if failed:
            print(f"FAILED to fetch: {failed:,}")
    return meta


def get_gsm_list(db_path: str) -> list[str]:
    con = duckdb.connect(db_path, read_only=True)
    try:
        gsms = con.execute(
            'SELECT DISTINCT upper(trim("sample")) AS gsm FROM geo_expr WHERE "sample" IS NOT NULL'
        ).fetchdf()["gsm"].tolist()
    finally:
        con.close()
    return sorted(gsms)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, help="Path to celllineselector.db")
    ap.add_argument("--destdir", default="geo_meta", help="Output dir for per-GSM SOFT files")
    ap.add_argument("--limit", type=int, default=None,
                    help="Harvest only the first N GSMs (smoke test)")
    args = ap.parse_args()

    gsms = get_gsm_list(args.db)
    print(f"{len(gsms):,} distinct GSM accessions in geo_expr")
    if args.limit:
        gsms = gsms[:args.limit]
        print(f"--limit set: harvesting only {len(gsms)}")

    meta = harvest_geo_methods(gsms, destdir=args.destdir)
    out_csv = Path(args.destdir).parent / "geo_meta_harvest_summary.csv"
    meta.to_csv(out_csv, index=False)
    print(f"\nharvest summary -> {out_csv}")


if __name__ == "__main__":
    main()
