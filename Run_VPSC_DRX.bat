@echo off
setlocal
cd /d "%~dp0"
call "%~dp0app\python.bat" simulate.py %*
if errorlevel 1 pause
