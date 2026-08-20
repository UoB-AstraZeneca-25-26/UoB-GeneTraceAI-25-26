"""Puts AgentDevelopment/ and Schema/ on sys.path.

They're plain script folders, not installable packages, so anything that
wants to `import FunctionCalling` / `import BusinessFlow` needs this run
first. Mirrors Tests/conftest.py's setup for the non-test runtime.
Importing this module is the side effect — nothing else to call.
"""

import sys
from pathlib import Path

_AIAGENT_DIR = Path(__file__).resolve().parent.parent
for _subdir in ("AgentDevelopment", "Schema"):
    _path = str(_AIAGENT_DIR / _subdir)
    if _path not in sys.path:
        sys.path.insert(0, _path)
