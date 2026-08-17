import sys
from pathlib import Path

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
