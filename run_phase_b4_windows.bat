@echo off
setlocal
cd /d "%~dp0"
if not exist .venv\Scripts\python.exe (
  py -m venv .venv 2>nul || python -m venv .venv
)
call .venv\Scripts\activate.bat
python -m pip install -r requirements.txt
if errorlevel 1 goto :fail

echo.
echo Running Phase B4 full-lag morphology and time localization...
python scripts\phase_b4_localization.py
if errorlevel 1 goto :fail

echo.
echo Phase B4 complete.
echo See results\phase_b4_localization\localization_summary.json
echo See results\phase_b4_localization\sliding_window_summary.csv
pause
exit /b 0

:fail
echo.
echo Phase B4 failed with exit code %errorlevel%.
pause
exit /b %errorlevel%
