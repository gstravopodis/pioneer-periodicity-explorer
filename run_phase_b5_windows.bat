@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo Virtual environment not found. Run a previous phase launcher first or create .venv.
  exit /b 1
)
.venv\Scripts\python.exe scripts\phase_b5_periodicity.py
if errorlevel 1 exit /b %errorlevel%
echo.
echo Phase B5 complete. See results\phase_b5_periodicity\periodicity_summary.json
pause
