@echo off
setlocal
cd /d "%~dp0"
call "%~dp0app\python.bat" control_panel.py
if errorlevel 1 pause
