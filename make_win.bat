@echo off
setlocal EnableExtensions
REM Windows counterpart of the Makefile.
REM
REM Usage:  make_win.bat [daily] [clean] [wheel] [app-install] [uninstall] [test]
REM   daily   use C:\Program Files\ChimeraX-Daily instead of the stable install
REM The actions run in the order listed above, whatever order you give them in.
REM e.g.    make_win.bat clean app-install test

set "BUNDLE_BASE_NAME=ISOLDE-DNA_RESTRAINTS"
set "CX_APP=C:\Program Files\ChimeraX"
set DO_CLEAN=
set DO_WHEEL=
set DO_INSTALL=
set DO_UNINSTALL=
set DO_TEST=
for %%A in (%*) do (
    if /i "%%A"=="daily"       set "CX_APP=C:\Program Files\ChimeraX-Daily"
    if /i "%%A"=="clean"       set DO_CLEAN=1
    if /i "%%A"=="wheel"       set DO_WHEEL=1
    if /i "%%A"=="app-install" set DO_INSTALL=1
    if /i "%%A"=="uninstall"   set DO_UNINSTALL=1
    if /i "%%A"=="test"        set DO_TEST=1
)

set "CX=%CX_APP%\bin\ChimeraX-console.exe"
if not exist "%CX%" (
    echo ERROR: ChimeraX not found at "%CX%" 1>&2
    exit /b 1
)

if defined DO_CLEAN     "%CX%" --nogui --safemode --exit --cmd "devel clean ."
if defined DO_WHEEL     "%CX%" --nogui --safemode --exit --cmd "devel build ."
if defined DO_INSTALL   "%CX%" --nogui --safemode --exit --cmd "devel install ."
if defined DO_UNINSTALL "%CX%" --nogui --safemode --exit --cmd "toolshed uninstall %BUNDLE_BASE_NAME%"
if defined DO_TEST      "%CX%" --nogui --exit --script tests\test_commands.py

endlocal
