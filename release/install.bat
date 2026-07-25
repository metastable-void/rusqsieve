@echo off
setlocal

where powershell.exe >nul 2>&1
if errorlevel 1 (
    echo install.bat: Windows PowerShell is required. 1>&2
    exit /b 1
)

powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0install.ps1" %*
exit /b %errorlevel%

