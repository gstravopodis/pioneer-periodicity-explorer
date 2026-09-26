#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [ ! -x .venv/bin/python ]; then
  echo "Missing .venv. Run run_phase_a_unix.sh first."
  exit 1
fi
echo "Phase B3 first-pass significance: 199 surrogates per channel/variant."
echo "This may take several minutes."
.venv/bin/python scripts/phase_b3_significance.py --simulations 199
