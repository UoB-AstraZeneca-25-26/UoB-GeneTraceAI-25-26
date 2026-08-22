#!/usr/bin/env bash
set -euo pipefail

if [ -d ".venv" ]; then
  echo ".venv already exists. To (re)install, remove it or activate it with: source .venv/bin/activate"
  exit 0
fi

python3 -m venv .venv
# shellcheck source=/dev/null
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt

echo "Virtual environment created and dependencies installed.\nActivate it with: source .venv/bin/activate"
