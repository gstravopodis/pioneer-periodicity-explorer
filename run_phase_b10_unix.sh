#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PYTHON=".venv/bin/python"
if [[ ! -x "$PYTHON" ]]; then
  echo "Virtual environment not found. Run the setup used for earlier phases first."
  exit 1
fi
"$PYTHON" scripts/phase_b10_intermittent_significance.py --simulations 1999
