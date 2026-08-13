"""Profile every parquet under data/ and src/**/outputs/ without loading full frames."""
import sys, os, json, glob, datetime
import pyarrow.parquet as pq

ROOT = r"C:\Disertation\UoB-GeneTraceAI-25-26"

targets = []
for pat in [
    "data/**/*.parquet",
    "src/**/outputs/**/*.parquet",
    "outputs/**/*.parquet",
    "cleaned_track_data/**/*.parquet",
    "reference/**/*.parquet",
]:
    targets += glob.glob(os.path.join(ROOT, pat), recursive=True)
targets = sorted(set(targets))

out = []
for p in targets:
    rel = os.path.relpath(p, ROOT).replace("\\", "/")
    try:
        f = pq.ParquetFile(p)
        md = f.metadata
        sch = f.schema_arrow
        rec = {
            "path": rel,
            "rows": md.num_rows,
            "cols": md.num_columns,
            "size_mb": round(os.path.getsize(p) / 1048576, 2),
            "mtime": datetime.datetime.fromtimestamp(os.path.getmtime(p)).strftime("%Y-%m-%d %H:%M"),
            "columns": [(n, str(t)) for n, t in zip(sch.names, sch.types)],
        }
    except Exception as e:
        rec = {"path": rel, "error": repr(e)}
    out.append(rec)

with open(os.path.join(os.path.dirname(__file__), "parquet_profile.json"), "w") as fh:
    json.dump(out, fh, indent=1)

for r in out:
    if "error" in r:
        print(f"ERROR {r['path']}: {r['error']}")
        continue
    print(f"\n=== {r['path']}  rows={r['rows']:,} cols={r['cols']} {r['size_mb']}MB  mtime={r['mtime']}")
    print("    " + "; ".join(f"{n}:{t}" for n, t in r["columns"][:60]))
    if len(r["columns"]) > 60:
        print(f"    ... +{len(r['columns'])-60} more cols")
