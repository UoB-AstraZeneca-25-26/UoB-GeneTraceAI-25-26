"""
final_pipeline/run_all.py
--------------------------
End-to-end pipeline runner. Run from the repo root:

  python final_pipeline/run_all.py

Stages run in order. Each stage imports its own script and calls run().
Skip completed stages with --from <stage_number>:

  python final_pipeline/run_all.py --from 4   # restart from 03_Altercations

Outputs land in final_pipeline/outputs/.
"""
import sys, time, argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

STAGES = [
    (0, "00_harmonisation", "Check warehouse exists"),
    (1, "01_Transcriptomics", "RNA z-scores (transcriptomics.py, full genome)  ~2h"),
    (2, "02_Proteinomics",    "Protein tier + z-scores (M2)              ~15-20 min"),
    (3, "03_Altercations",    "Mutations + Fusions + CNA scoring"),
    (4, "Scoring",            "core_score + driver routing + tiers"),
    (5, "Ranking",            "Evidence ledger (all pairs)"),
    (6, "Validation",         "Held-out evaluation (hit@20)"),
    (7, "Testing",            "Diagnostic tests on outputs"),
]


def run_stage_0():
    from importlib import import_module
    mod = import_module("00_harmonisation.run")
    mod.run()


def run_stage_1():
    from importlib import import_module
    mod = import_module("01_Transcriptomics.run_transcriptomics_stage")
    mod.run()


def run_stage_2():
    from importlib import import_module
    import_module("02_Proteinomics.platform_tier").run()
    import_module("02_Proteinomics.run_protein_m2_stage").run()


def run_stage_3():
    from importlib import import_module
    import_module("03_Altercations.mutations_scoring").run()
    import_module("03_Altercations.fusions_scoring").run()
    import_module("03_Altercations.cna_layer").run()


def run_stage_4():
    from importlib import import_module
    import_module("Scoring.core_score").run()
    import_module("Scoring.driver_routing").run()
    import_module("Scoring.confidence_tiers").run()


def run_stage_5():
    from importlib import import_module
    import_module("Ranking.build_evidence_ledger").run()


def run_stage_6():
    from importlib import import_module
    import_module("Validation.eval").run()


def run_stage_7():
    from importlib import import_module
    import_module("Testing.run_tests")


RUNNERS = [run_stage_0, run_stage_1, run_stage_2, run_stage_3,
           run_stage_4, run_stage_5, run_stage_6, run_stage_7]


def main():
    parser = argparse.ArgumentParser(description="GeneTraceAI end-to-end pipeline")
    parser.add_argument("--from", dest="from_stage", type=int, default=0,
                        help="Start from this stage number (0-7)")
    parser.add_argument("--only", dest="only_stage", type=int, default=None,
                        help="Run only this stage number")
    args = parser.parse_args()

    print("\n" + "=" * 70)
    print("GeneTraceAI — Full Pipeline")
    print("=" * 70)
    for n, name, desc in STAGES:
        status = "  " if n >= args.from_stage else "SKIP"
        print(f"  [{status}] Stage {n}: {name:<22s} {desc}")
    print()

    t_total = time.time()
    for n, name, desc in STAGES:
        if args.only_stage is not None and n != args.only_stage:
            continue
        if n < args.from_stage:
            continue

        print(f"\n{'='*70}")
        print(f"Stage {n}: {name}")
        print(f"{'='*70}")
        t = time.time()
        try:
            RUNNERS[n]()
            print(f"\nStage {n} complete ({time.time()-t:.0f}s)")
        except SystemExit as e:
            print(f"\nStage {n} STOPPED: {e}")
            print("Fix the issue above and rerun with --from", n)
            sys.exit(1)

    print(f"\n{'='*70}")
    print(f"Pipeline complete. Total time: {time.time()-t_total:.0f}s")
    print(f"Outputs in: final_pipeline/outputs/")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
