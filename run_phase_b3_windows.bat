@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo Missing .venv. Run run_phase_a_windows.bat first.
  pause
  exit /b 1
)
echo Phase B3 first-pass significance: 199 surrogates per channel/variant.
echo This may take several minutes.
".venv\Scripts\python.exe" scripts\phase_b3_significance.py --simulations 199
if errorlevel 1 (
  echo Phase B3 failed with exit code %errorlevel%.
  pause
  exit /b %errorlevel%
)
pause
