#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [ ! -x .venv/bin/python ]; then
  echo "Missing .venv. Run run_phase_a_unix.sh first." >&2
  exit 1
fi
.venv/bin/python scripts/phase_b2_screening.py --all-channels
