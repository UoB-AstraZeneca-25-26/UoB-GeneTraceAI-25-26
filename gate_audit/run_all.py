"""
gate_audit/run_all.py
---------------------
Run the whole audit in order, from the repo root:

    python gate_audit/run_all.py

Options:
    --skip-panel        reuse an existing outputs/panel.parquet (00 is the slow step)
    --quick             smaller gene sample and fewer bootstrap replicates
    --only 03 05        run just those steps (00 is implied if the panel is absent)

Order matters. 00 builds the panel every later step reads; 01-06 are independent
of each other given the panel, but they are numbered in the order the findings
should be read, and each one's docstring assumes the earlier ones.

Every script is read-only with respect to the pipeline. Nothing here writes to
src/pipeline/outputs/ or to any production table.
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

STEPS = [
    ("00", "00_build_panel.py", "build the scored panel", []),
    ("01", "01_tie_mass.py", "tie mass and tie handling", []),
    ("02", "02_gate_valid_genes.py", "gate on valid genes only", []),
    ("03", "03_denominator.py", "the denominator / profiling confound", []),
    ("04", "04_cluster_bootstrap.py", "cluster bootstrap for an honest p", []),
    ("05", "05_region_sweep.py", "sweep the post-hoc free parameters", []),
    ("06", "06_screen_replication.py", "Broad vs Sanger on real Chronos", []),
    ("07", "07_independent_labels.py", "GDSC2 as an independent perturbation", []),
    ("08", "08_release_refresh.py", "rebuild on DepMap 24Q4", []),
    ("09", "09_perturbation_types.py", "CRISPR vs RNAi vs drug", []),
]

QUICK = {
    "00": ["--n-genes", "3000"],
    "03": ["--n-genes", "800"],
    "04": ["--n-boot", "200", "--n-genes", "600"],
    "08": ["--n-genes", "1500"],
    "09": ["--n-genes", "1500"],
}
FULL = {
    "00": ["--n-genes", "20000"],
    "03": ["--n-genes", "2000"],
    "04": ["--n-boot", "1000", "--n-genes", "3900"],
    "08": ["--n-genes", "4000"],
    "09": ["--n-genes", "4000"],
}

ap = argparse.ArgumentParser()
ap.add_argument("--skip-panel", action="store_true")
ap.add_argument("--quick", action="store_true")
ap.add_argument("--only", nargs="*", default=None)
args = ap.parse_args()

extra = QUICK if args.quick else FULL
panel = HERE / "outputs" / "panel.parquet"

todo = []
for num, script, desc, _ in STEPS:
    if num == "00":
        if args.skip_panel and panel.exists():
            continue
        if args.only and "00" not in args.only and panel.exists():
            continue
    elif args.only and num not in args.only:
        continue
    todo.append((num, script, desc))

print("=" * 78)
print("GATE AUDIT -- running", ", ".join(n for n, _, _ in todo))
print("=" * 78)

failed = []
for num, script, desc in todo:
    cmd = [sys.executable, str(HERE / script)] + extra.get(num, [])
    print(f"\n\n{'#' * 78}\n# {num}  {desc}\n# {' '.join(cmd[1:])}\n{'#' * 78}\n")
    t0 = time.time()
    rc = subprocess.call(cmd, cwd=str(HERE.parent))
    dt = time.time() - t0
    if rc != 0:
        failed.append(num)
        print(f"\n!! step {num} FAILED (exit {rc}) after {dt:.0f}s")
        if num == "00":
            print("!! the panel is required by every later step -- stopping")
            break
    else:
        print(f"\n-- step {num} done in {dt:.0f}s")

print("\n" + "=" * 78)
if failed:
    print(f"FINISHED WITH FAILURES: {', '.join(failed)}")
    sys.exit(1)
print("ALL STEPS COMPLETED")
print(f"results in {HERE / 'outputs'}")
print("=" * 78)
