@echo off
setlocal
cd /d "%~dp0"
if not defined VPSC_GFORTRAN set "VPSC_GFORTRAN=gfortran"
"%VPSC_GFORTRAN%" -fdefault-real-8 -std=legacy -O2 -static -fcheck=bounds -fbacktrace -ffpe-trap=invalid,zero,overflow vpsc8dim.for vpsc8alloc.f90 vpsc8drx.f90 vpsc8main.for vpsc8sub.for -o ..\bin\vpsc_drx.new.exe
if errorlevel 1 exit /b 1
move /y ..\bin\vpsc_drx.new.exe ..\bin\vpsc_drx.exe
if errorlevel 1 exit /b 1
call "%~dp0build_original.bat"
exit /b %ERRORLEVEL%
