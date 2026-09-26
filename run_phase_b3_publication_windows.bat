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
echo Running publication-resolution Phase B3 with 1999 surrogates per test...
echo This is substantially heavier than the 199-simulation screening run.
python scripts\phase_b3_significance.py --simulations 1999 --out results\phase_b3_publication
if errorlevel 1 goto :fail
python scripts\review_b3.py --input results\phase_b3_publication\significance_summary.json --out results\phase_b3_publication
if errorlevel 1 goto :fail

echo.
echo Publication-resolution B3 complete.
echo See results\phase_b3_publication\significance_summary.json
echo See results\phase_b3_publication\candidate_review.json
pause
exit /b 0

:fail
echo.
echo Phase B3 publication run failed with exit code %errorlevel%.
pause
exit /b %errorlevel%
