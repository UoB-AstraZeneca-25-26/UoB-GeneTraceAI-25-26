"""
logging_utils.py

Shared logging setup for the pipeline. Errors are written to a dated
log file with full context (file, folder, extension, read settings,
step, error type/message, traceback) so an admin can diagnose a bad
file without needing to reproduce the failure themselves.

Import get_logger() to set up a named logger per script, then use
log_file_error() for failures tied to a specific raw/output file, or
log_error() for failures not tied to a file (e.g. a dataset-level
cleaning error). Every pipeline script sharing this module writes to
the same day's log file (logs/pipeline_<date>.log), so all stages'
errors end up in one place.

Design notes
------------
* Two levels, two destinations. The console carries INFO and above, so
  a run reads as progress; the file carries WARNING and above, so the
  log is signal rather than a transcript.
* Errors are logged as ``key=value`` pairs rather than prose, so a log
  line can be grepped or split on ``|`` without parsing English.
* Structured context is what makes a failure reproducible: the read
  settings recorded alongside a parse error are usually enough to
  identify the cause without re-running anything.
* The dated filename means one file per day across every stage, so the
  ordering between stages of a single run is preserved.
* No side effects at import beyond defining the paths — the log
  directory is created on first :func:`get_logger` call, not here.

Log location
------------
``<project root>/logs/pipeline_<YYYY-MM-DD>.log``, where the project
root is resolved two levels up from this file.
"""

import logging
import sys
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LOG_DIR = PROJECT_ROOT / "logs"


def get_logger(name: str) -> logging.Logger:
    """
    Return a logger that writes:
      - INFO and above to the console
      - WARNING and above to logs/pipeline_<date>.log

    Safe to call multiple times with the same name — handlers are
    only attached once.

    The split is deliberate: the console gets a bare message format so
    progress output reads cleanly, while the file gets timestamp, level
    and logger name so entries from different stages of the same day's
    run can be told apart afterwards.

    Parameters
    ----------
    name : str
        Logger name, conventionally the calling script's stem (e.g.
        ``"01_data_cleaning"``). It appears in every file log line, so it
        is what identifies which stage raised an entry.

    Returns
    -------
    logging.Logger
        Configured logger. The same instance is returned on repeat calls
        with the same name.

    Notes
    -----
    The logger's own level is DEBUG, with filtering left to the two
    handlers — so a future handler can be added at any level without
    reconfiguring the logger. The log directory is created here if it
    does not exist. The date is resolved at call time, so a run spanning
    midnight writes its later entries to the next day's file.
    """
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger

    logger.setLevel(logging.DEBUG)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_file = LOG_DIR / f"pipeline_{datetime.now():%Y-%m-%d}.log"

    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setLevel(logging.WARNING)
    file_handler.setFormatter(logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter("%(message)s"))

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    return logger


def _log_structured(logger: logging.Logger, details: dict) -> None:
    """Shared formatter/emitter for log_file_error and log_error.

    Renders the details dict as pipe-separated ``key=value`` pairs and
    emits it at ERROR with the traceback attached.

    Parameters
    ----------
    logger : logging.Logger
        Logger from :func:`get_logger`.
    details : dict
        Context to record. Insertion order is preserved in the output, so
        the callers put the identifying fields first.

    Returns
    -------
    None

    Notes
    -----
    ``exc_info=True`` reads the exception currently being handled, so this
    must be called from inside an ``except`` block for the traceback to
    appear. Called outside one, the entry is still written but carries no
    traceback — which is why the ``error`` argument is also recorded
    explicitly by both callers rather than relied on from ``exc_info``.
    """
    detail_str = " | ".join(f"{k}={v}" for k, v in details.items())
    logger.error(detail_str, exc_info=True)


def log_file_error(logger: logging.Logger, *, file_path, step: str, error: Exception, **context) -> None:
    """
    Log a structured error about a specific raw data file: path,
    extension, folder, the pipeline step that failed, the exception
    type/message, and any extra context (sep, skiprows, index_col,
    key, etc.) passed as kwargs. Full traceback is included.

    The path is split into name, folder and extension as separate fields
    as well as being kept whole, so entries can be filtered by any of
    them — all failures in one folder, or all ``.gct`` failures — without
    string surgery on the full path.

    Parameters
    ----------
    logger : logging.Logger
        Logger from :func:`get_logger`.
    file_path : str or Path
        The file being read or written when the failure occurred.
        Coerced to :class:`~pathlib.Path`.
    step : str
        Function or stage name that failed, e.g.
        ``"convert_text_files_to_csv"``. Recorded first, since it is the
        usual entry point when reading the log.
    error : Exception
        The caught exception; its type and message are recorded.
    **context
        Any extra fields worth recording — typically the read settings
        (``sep``, ``skiprows``, ``index_col``, ``key``) for a parse
        failure, or the dataset name and shape for a write failure.

    Returns
    -------
    None
        The entry is written to the log as a side effect.

    Notes
    -----
    Keyword-only after ``logger``, so a call site cannot accidentally
    transpose ``step`` and ``file_path``. Call from inside the ``except``
    block that caught ``error`` so the traceback is captured.
    """
    file_path = Path(file_path)
    _log_structured(logger, {
        "step": step,
        "file": file_path.name,
        "folder": str(file_path.parent),
        "extension": file_path.suffix,
        "full_path": str(file_path),
        "error_type": type(error).__name__,
        "error_message": str(error),
        **context,
    })


def log_error(logger: logging.Logger, *, step: str, error: Exception, **context) -> None:
    """
    Log a structured error not tied to a specific raw file (e.g. a
    dataset-level failure during cleaning or spot-checking). Includes
    the step, error type/message, any extra context kwargs, and the
    full traceback.

    Use this where no single file is at fault — a roster build, a
    cross-table join, a diagnostic over several datasets. Where one file
    can be named, :func:`log_file_error` records more and is preferred.

    Parameters
    ----------
    logger : logging.Logger
        Logger from :func:`get_logger`.
    step : str
        Function or stage name that failed, e.g. ``"build_gene_roster"``.
    error : Exception
        The caught exception; its type and message are recorded.
    **context
        Any extra fields worth recording — typically ``dataset``, or
        ``source_dir`` for a load failure.

    Returns
    -------
    None
        The entry is written to the log as a side effect.

    Notes
    -----
    Keyword-only after ``logger``. Call from inside the ``except`` block
    that caught ``error`` so the traceback is captured.
    """
    _log_structured(logger, {
        "step": step,
        "error_type": type(error).__name__,
        "error_message": str(error),
        **context,
    })