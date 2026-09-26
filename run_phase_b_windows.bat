@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo Existing .venv not found. Run run_phase_a_windows.bat once first.
  pause
  exit /b 1
)
if not exist "data\nasa\processed\p10_cpi_daily.csv" (
  echo P10 processed CPI data missing. Run Phase A first.
  pause
  exit /b 1
)
if not exist "data\nasa\processed\p11_cpi_daily.csv" (
  echo P11 processed CPI data missing. Run Phase A first.
  pause
  exit /b 1
)
".venv\Scripts\python.exe" scripts\calendar_reanalysis.py --channel p_11_20_mev
if errorlevel 1 (
  echo Phase B1 calendar-aware analysis failed.
  pause
  exit /b 1
)
echo.
echo Phase B1 complete. See results\phase_b_calendar\calendar_reanalysis_summary.json
pause
