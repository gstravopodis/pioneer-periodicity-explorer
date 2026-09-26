@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo Missing .venv. Run run_phase_a_windows.bat first.
  pause
  exit /b 1
)
".venv\Scripts\python.exe" scripts\phase_b2_screening.py --all-channels
if errorlevel 1 (
  echo Phase B2 failed with exit code %errorlevel%.
  pause
  exit /b %errorlevel%
)
pause
