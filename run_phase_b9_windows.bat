@echo off
setlocal
cd /d "%~dp0"

if not exist ".venv\Scripts\python.exe" (
  echo Virtual environment not found. Run the setup used for earlier phases first.
  pause
  exit /b 1
)

echo Running Phase B9 publication-resolution candidate spectral significance...
.venv\Scripts\python.exe scripts\phase_b9_spectral_significance.py --simulations 1999
if errorlevel 1 (
  echo Phase B9 failed.
  pause
  exit /b 1
)

echo.
echo Phase B9 complete. See results\phase_b9_spectral_significance\spectral_significance_summary.json
pause
