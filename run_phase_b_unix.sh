#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ ! -x .venv/bin/python ]]; then
  echo "Existing .venv not found. Run Phase A once first."
  exit 1
fi
.venv/bin/python scripts/calendar_reanalysis.py --channel p_11_20_mev
echo "Phase B1 complete. See results/phase_b_calendar/calendar_reanalysis_summary.json"
