@echo off
setlocal
cd /d "%~dp0"
if not exist .venv\Scripts\python.exe (
  echo Missing .venv. Create the environment and install requirements first.
  pause
  exit /b 1
)
.venv\Scripts\python.exe scripts\phase_b8_resolution_stress.py
if errorlevel 1 (
  echo.
  echo Phase B8 failed.
  pause
  exit /b 1
)
echo.
echo Phase B8 complete. See results\phase_b8_resolution_stress\resolution_stress_summary.json
pause
