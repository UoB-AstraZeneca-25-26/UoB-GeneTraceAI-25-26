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
    """Shared formatter/emitter for log_file_error and log_error."""
    detail_str = " | ".join(f"{k}={v}" for k, v in details.items())
    logger.error(detail_str, exc_info=True)


def log_file_error(logger: logging.Logger, *, file_path, step: str, error: Exception, **context) -> None:
    """
    Log a structured error about a specific raw data file: path,
    extension, folder, the pipeline step that failed, the exception
    type/message, and any extra context (sep, skiprows, index_col,
    key, etc.) passed as kwargs. Full traceback is included.
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
    """
    _log_structured(logger, {
        "step": step,
        "error_type": type(error).__name__,
        "error_message": str(error),
        **context,
    })