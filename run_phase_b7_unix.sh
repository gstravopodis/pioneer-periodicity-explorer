#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PY=".venv/bin/python"
if [[ ! -x "$PY" ]]; then
  echo "Missing .venv. Create the environment and install requirements first." >&2
  exit 1
fi
"$PY" scripts/phase_b7_sampling_tracking.py
