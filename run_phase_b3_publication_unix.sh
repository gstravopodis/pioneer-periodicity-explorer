#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [ ! -x .venv/bin/python ]; then
  python3 -m venv .venv
fi
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/phase_b3_significance.py --simulations 1999 --out results/phase_b3_publication
python scripts/review_b3.py --input results/phase_b3_publication/significance_summary.json --out results/phase_b3_publication
