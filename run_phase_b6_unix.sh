#!/usr/bin/env bash
set -euo pipefail
if [[ -x .venv/bin/python ]]; then
  PY=.venv/bin/python
else
  PY=python3
fi
"$PY" scripts/phase_b6_periodicity_atlas.py
