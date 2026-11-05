@echo off
setlocal
set "PYTHONUTF8=1"
set "OPENBLAS_NUM_THREADS=1"
set "OMP_NUM_THREADS=1"
if defined VPSC_PYTHON goto launch
if exist "%~dp0..\.venv\Scripts\python.exe" (
  set "VPSC_PYTHON=%~dp0..\.venv\Scripts\python.exe"
) else if exist "E:\ProgramData\anaconda3\python.exe" (
  set "VPSC_PYTHON=E:\ProgramData\anaconda3\python.exe"
) else (
  set "VPSC_PYTHON=python"
)
:launch
"%VPSC_PYTHON%" -X utf8 %*
exit /b %ERRORLEVEL%
