import sys
from pathlib import Path

import pytest

# Windows consoles default stdout to the cp1252 codepage, which can't encode
# Unicode punctuation (e.g. narrow no-break spaces) that Groq models return.
# AgentExecutor(verbose=True) prints raw LLM output straight to stdout, so
# without this it crashes mid-test instead of just mis-rendering a character.
for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8", errors="replace")

_AIAGENT_DIR = Path(__file__).resolve().parent.parent
for subdir in ("AgentDevelopment", "Schema"):
    path = str(_AIAGENT_DIR / subdir)
    if path not in sys.path:
        sys.path.insert(0, path)


@pytest.fixture(autouse=True)
def _fresh_http_client_per_test():
    """http_utils.get_http_client() is a process-wide singleton — correct
    for production (one event loop for the app's whole life) but unsafe
    across pytest-asyncio's per-test event loops: a client created in one
    test's loop still holds keep-alive connections when the next test's
    loop tries to reuse it, and closing a connection bound to an
    already-closed loop raises `RuntimeError: Event loop is closed`
    (reproducible: Tests/test_gene_alias.py's live tests fail in this
    exact way when run together, but pass individually).

    There's no safe way to gracefully aclose() a client whose event loop is
    already gone, so just drop the reference — the next test creates a
    fresh client bound to its own loop on first use.
    """
    yield
    import http_utils

    http_utils._client = None
