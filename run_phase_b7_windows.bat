@echo off
setlocal
cd /d "%~dp0"
if not exist .venv\Scripts\python.exe (
  echo Missing .venv. Create the environment and install requirements first.
  pause
  exit /b 1
)
.venv\Scripts\python.exe scripts\phase_b7_sampling_tracking.py
if errorlevel 1 (
  echo.
  echo Phase B7 failed.
  pause
  exit /b 1
)
echo.
echo Phase B7 complete. See results\phase_b7_sampling_tracking\sampling_tracking_summary.json
pause
