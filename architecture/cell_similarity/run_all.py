"""
cell_similarity/run_all.py
--------------------------
Run the build and both validations in order, from the repo root:

    python cell_similarity/run_all.py

01 builds the three axes and every later step reads them, so `--skip-build`
reuses `outputs/*_similarity.parquet`. 04 is a query interface, not a batch step,
so it is not run here — see the README for how to call it.
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

STEPS = [
    ("01", "01_build_axes.py", "week 1 -- build the three axes, kept separate"),
    ("02", "02_validate_dependency.py", "week 2 -- do similar lines behave similarly?"),
    ("03", "03_validate_tissue.py", "week 2 -- do same-tissue lines come out close?"),
    ("05", "05_robustness.py", "cluster bootstrap, RNAi, common essentials"),
    ("06", "06_mirna_stability.py", "miRNA metric instability -- fix or demote"),
    ("07", "07_reliability.py", "reliability, disattenuation, power control"),
    ("08", "08_drug_arm.py", "third arm -- GDSC2 drug response"),
]

ap = argparse.ArgumentParser()
ap.add_argument("--skip-build", action="store_true")
ap.add_argument("--only", nargs="*", default=None)
args = ap.parse_args()

built = (HERE / "outputs" / "rna_similarity.parquet").exists()
todo = []
for num, script, desc in STEPS:
    if num == "01" and args.skip_build and built:
        continue
    if args.only and num not in args.only:
        continue
    todo.append((num, script, desc))

print("=" * 78)
print("CELL SIMILARITY -- running", ", ".join(n for n, _, _ in todo))
print("=" * 78)

failed = []
for num, script, desc in todo:
    print(f"\n\n{'#' * 78}\n# {num}  {desc}\n{'#' * 78}\n")
    t0 = time.time()
    rc = subprocess.call([sys.executable, str(HERE / script)], cwd=str(HERE.parent))
    dt = time.time() - t0
    if rc != 0:
        failed.append(num)
        print(f"\n!! step {num} FAILED (exit {rc}) after {dt:.0f}s")
        if num == "01":
            print("!! the axes are required by every later step -- stopping")
            break
    else:
        print(f"\n-- step {num} done in {dt:.0f}s")

print("\n" + "=" * 78)
if failed:
    print(f"FINISHED WITH FAILURES: {', '.join(failed)}")
    sys.exit(1)
print("ALL STEPS COMPLETED")
print(f"results in {HERE / 'outputs'}")
print("\nQuery the graph with:")
print("  python cell_similarity/04_query.py --gene TP53 --line ACH-000552")
print("=" * 78)
