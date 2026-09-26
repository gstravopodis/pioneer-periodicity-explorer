@echo off
setlocal
cd /d "%~dp0"
if not exist .venv (
  py -m venv .venv
)
call .venv\Scripts\activate
python -m pip install -r requirements.txt
if exist data\nasa\processed\p10_cpi_daily.csv if exist data\nasa\processed\p11_cpi_daily.csv (
  echo Existing processed NASA CPI files found; skipping re-download.
  python scripts\run_phase_a.py --channel p_11_20_mev --skip-download
) else (
  python scripts\run_phase_a.py --channel p_11_20_mev
)
pause
