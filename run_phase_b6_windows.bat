@echo off
setlocal
if not exist .venv\Scripts\python.exe (
  echo Virtual environment not found. Run the earlier setup/Phase A launcher first.
  exit /b 1
)
.venv\Scripts\python.exe scripts\phase_b6_periodicity_atlas.py
if errorlevel 1 exit /b %errorlevel%
echo.
echo Phase B6 complete. See results\phase_b6_periodicity_atlas\periodicity_atlas_summary.json
pause
