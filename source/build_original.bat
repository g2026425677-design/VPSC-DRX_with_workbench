@echo off
setlocal
cd /d "%~dp0"
if not defined VPSC_GFORTRAN set "VPSC_GFORTRAN=gfortran"
if not exist build_stock mkdir build_stock
cd build_stock
"%VPSC_GFORTRAN%" -fdefault-real-8 -std=legacy -static ..\stock_vpsc\vpsc8dim.for ..\stock_vpsc\vpsc8main.for ..\stock_vpsc\vpsc8sub.for -o ..\..\bin\vpsc_original.new.exe
if errorlevel 1 exit /b 1
"%VPSC_GFORTRAN%" -fdefault-real-8 -std=legacy -static ..\stock_vpsc\vpsc8dim.for ..\stock_vpsc\vpsc8main.for ..\stock_vpsc\vpsc8sub.for ..\stock_observer.f90 -Wl,--wrap=write_shear_activity_ -o ..\..\bin\vpsc_original_report.new.exe
if errorlevel 1 exit /b 1
move /y ..\..\bin\vpsc_original.new.exe ..\..\bin\vpsc_original.exe
if errorlevel 1 exit /b 1
move /y ..\..\bin\vpsc_original_report.new.exe ..\..\bin\vpsc_original_report.exe
exit /b %ERRORLEVEL%

