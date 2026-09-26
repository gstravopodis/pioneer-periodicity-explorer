#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PY=".venv/bin/python"
if [[ ! -x "$PY" ]]; then
  echo "Virtual environment not found. Run the setup used for earlier phases first."
  exit 1
fi
"$PY" scripts/phase_b9_spectral_significance.py --simulations 1999
