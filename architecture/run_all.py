"""
architecture/run_all.py
-------------------------
End-to-end pipeline runner. Run from the repo root:

  python architecture/run_all.py

Stages run in order. Each stage imports its own script and calls run().
Skip completed stages with --from <stage_number>:

  python architecture/run_all.py --from 4   # restart from 03_Altercations

Outputs land in architecture/outputs/.
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
    """
    Stage 0: verify the Stage-0 DuckDB warehouse exists.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Delegates to ``00_harmonisation.run.run()``, which raises if the
    warehouse (``celllineselector.db``) is missing — this stage does
    not build the warehouse itself, only checks for it.
    """
    from importlib import import_module
    mod = import_module("00_harmonisation.run")
    mod.run()


def run_stage_1():
    """
    Stage 1: compute RNA z-scores for the full genome (~2h).

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Delegates to ``01_Transcriptomics.run_transcriptomics_stage.run()``.
    Writes ``bulk_rna_z.parquet``.
    """
    from importlib import import_module
    mod = import_module("01_Transcriptomics.run_transcriptomics_stage")
    mod.run()


def run_stage_2():
    """
    Stage 2: compute protein platform tiers and z-scores (~15-20 min).

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Runs two sub-stages in order: ``02_Proteinomics.platform_tier``
    (writes ``protein_platform_tier.parquet``), then
    ``02_Proteinomics.run_protein_m2_stage`` (writes
    ``bulk_prot_z.parquet``) — the tier table must exist first since
    the M2 scorer reads it.
    """
    from importlib import import_module
    import_module("02_Proteinomics.platform_tier").run()
    import_module("02_Proteinomics.run_protein_m2_stage").run()


def run_stage_3():
    """
    Stage 3: score mutations, fusions, and copy-number alterations.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Runs ``03_Altercations.mutations_scoring``,
    ``fusions_scoring``, then ``cna_layer`` in that order (each
    independent of the others' outputs, but run sequentially for
    consistent logging). Writes ``mutations_scores.parquet``,
    ``fusions_scores.parquet``, ``cna_flags.parquet``.
    """
    from importlib import import_module
    import_module("03_Altercations.mutations_scoring").run()
    import_module("03_Altercations.fusions_scoring").run()
    import_module("03_Altercations.cna_layer").run()


def run_stage_4():
    """
    Stage 4: combine RNA/protein into core_score, then tier by driver evidence.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Runs, in order (each depends on the previous one's output):
    ``Scoring.core_score`` (writes ``core_score.parquet``), then
    ``Scoring.driver_routing`` (writes ``flags_with_driver.parquet``),
    then ``Scoring.confidence_tiers`` (writes
    ``predictions_with_confidence.parquet``).
    """
    from importlib import import_module
    import_module("Scoring.core_score").run()
    import_module("Scoring.driver_routing").run()
    import_module("Scoring.confidence_tiers").run()


def run_stage_5():
    """
    Stage 5: build the full evidence ledger for every scored pair.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Delegates to ``Ranking.build_evidence_ledger.run()``. Writes
    ``evidence_ledger.parquet``.
    """
    from importlib import import_module
    import_module("Ranking.build_evidence_ledger").run()


def run_stage_6():
    """
    Stage 6: held-out evaluation (hit@20).

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Delegates to ``Validation.eval.run()``. Writes
    ``stage4_eval_summary.json``.
    """
    from importlib import import_module
    import_module("Validation.eval").run()


def run_stage_7():
    """
    Stage 7: run diagnostic tests against the pipeline's outputs.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Imports ``Testing.run_tests``, which runs its checks as
    import-time side effects (unlike every other stage, this module
    is not called via a ``.run()`` method — importing it is running
    it).
    """
    from importlib import import_module
    import_module("Testing.run_tests")


RUNNERS = [run_stage_0, run_stage_1, run_stage_2, run_stage_3,
           run_stage_4, run_stage_5, run_stage_6, run_stage_7]


def main():
    """
    Run the pipeline stages selected by ``--from``/``--only``, in order.

    Parameters
    ----------
    None
        Reads ``--from <int>`` (default 0, run this stage onward) and
        ``--only <int>`` (default None, run every stage in range) from
        ``sys.argv``.

    Returns
    -------
    None
        Prints a stage plan, then runs each selected stage's runner
        from ``RUNNERS``, printing per-stage and total elapsed time.

    Notes
    -----
    If a stage's runner raises ``SystemExit`` (every stage's ``run()``
    does this when a required input file is missing), this function
    prints the error, tells the user to rerun with ``--from <n>`` from
    that stage, and calls ``sys.exit(1)`` — it does not continue to
    later stages after a failure.
    """
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
    print(f"Outputs in: architecture/outputs/")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
