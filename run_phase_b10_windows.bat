@echo off
setlocal
cd /d "%~dp0"

if not exist ".venv\Scripts\python.exe" (
  echo Virtual environment not found. Run the setup used for earlier phases first.
  pause
  exit /b 1
)

echo Running Phase B10 intermittent same-window/same-period significance...
echo This uses 1999 AR(1) surrogates and may take a while.
.venv\Scripts\python.exe scripts\phase_b10_intermittent_significance.py --simulations 1999
if errorlevel 1 (
  echo Phase B10 failed.
  pause
  exit /b 1
)

echo.
echo Phase B10 complete. See results\phase_b10_intermittent_significance\intermittent_significance_summary.json
pause
