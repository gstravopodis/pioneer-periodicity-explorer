#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [ -x .venv/bin/python ]; then
  PY=.venv/bin/python
else
  PY=python3
fi
"$PY" scripts/phase_b5_periodicity.py
printf '\nPhase B5 complete. See results/phase_b5_periodicity/periodicity_summary.json\n'
