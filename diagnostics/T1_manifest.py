"""
T1 — Environment and provenance manifest.
Writes diagnostics/out/T1_manifest.txt
Run from repo root: python diagnostics/T1_manifest.py
"""
import hashlib, os, sys, time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT  = REPO / "diagnostics" / "out"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO / "final_pipeline"))
from config import DB, CORE_SCORE, PROT_TIER, FUSIONS_SCR, MUTATIONS_SCR

lines = []

def log(s=""):
    print(s)
    lines.append(s)

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def file_entry(label, path):
    p = Path(path)
    if not p.exists():
        log(f"[PROVENANCE] {label}: MISSING — {path}")
        return False
    sz   = p.stat().st_size
    mt   = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(p.stat().st_mtime))
    h    = sha256(p)
    log(f"[PROVENANCE] {label}")
    log(f"    path   = {p}")
    log(f"    sha256 = {h}")
    log(f"    size   = {sz:,} bytes")
    log(f"    mtime  = {mt}")
    return True

# ── 1. Versions ────────────────────────────────────────────────────────────────
log("=" * 72)
log("T1 — ENVIRONMENT AND PROVENANCE MANIFEST")
log("=" * 72)
log()
log("[ENV] Python: " + sys.version.replace("\n", " "))
import pandas as pd, numpy as np, scipy, duckdb
log(f"[ENV] pandas  {pd.__version__}")
log(f"[ENV] numpy   {np.__version__}")
log(f"[ENV] scipy   {scipy.__version__}")
log(f"[ENV] duckdb  {duckdb.__version__}")
log()

# ── 2. File manifest ────────────────────────────────────────────────────────────
log("=" * 72)
log("FILE MANIFEST")
log("=" * 72)

all_present = True
files_to_check = [
    ("Chronos Achilles HDF5",  REPO / "data/DepMap_Chronos/GeneFitnessEffect_Chronos_Achilles.hdf5"),
    ("Chronos Score HDF5",     REPO / "data/DepMap_Chronos/GeneFitnessEffect_Chronos_Score.hdf5"),
    ("chronos_long.parquet",   REPO / "validation/prepared/chronos_long.parquet"),
    ("DEMETER2 Achilles CSV",  REPO / "data/DEMETER2/D2_Achilles_gene_dep_scores.csv"),
    ("DEMETER2 DRIVE CSV",     REPO / "data/DEMETER2/D2_DRIVE_gene_dep_scores.csv"),
    ("DuckDB warehouse",       DB),
    ("core_score.parquet",     CORE_SCORE),
    ("protein_platform_tier",  PROT_TIER),
    ("fusions_scores",         FUSIONS_SCR),
    ("mutations_scores",       MUTATIONS_SCR),
]

# GDSC files
gdsc_dir = REPO / "data/GDSC"
if gdsc_dir.exists():
    for gf in sorted(gdsc_dir.iterdir()):
        files_to_check.append((f"GDSC/{gf.name}", gf))
else:
    log(f"[PROVENANCE] data/GDSC/: MISSING — directory does not exist")
    all_present = False

for label, path in files_to_check:
    ok = file_entry(label, path)
    if not ok:
        all_present = False
    log()

# ── 3. DuckDB warehouse tables ─────────────────────────────────────────────────
log("=" * 72)
log("DUCKDB WAREHOUSE — TABLES")
log("=" * 72)
if DB.exists():
    con = duckdb.connect(str(DB), read_only=True)
    tables = con.execute("SHOW TABLES").df()
    log(f"Tables ({len(tables)}):")
    for tbl in tables["name"]:
        try:
            n = con.execute(f"SELECT count(*) FROM main.{tbl}").fetchone()[0]
            cols = con.execute(f"DESCRIBE main.{tbl}").df()
            col_str = ", ".join(f"{r['column_name']}:{r['column_type']}" for _, r in cols.iterrows())
            log(f"  {tbl}: {n:,} rows | {col_str[:120]}")
        except Exception as e:
            log(f"  {tbl}: ERROR — {e}")
    con.close()
else:
    log("DuckDB warehouse MISSING — cannot inspect tables")
    all_present = False

log()
log("=" * 72)
if all_present:
    log("PASS — all expected files present")
else:
    log("FAIL — one or more files MISSING (see above)")
log("=" * 72)

# ── Write manifest ─────────────────────────────────────────────────────────────
out_path = OUT / "T1_manifest.txt"
out_path.write_text("\n".join(lines), encoding="utf-8")
print(f"\nManifest written: {out_path}")
