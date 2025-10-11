@echo off
echo Generating QR code...
echo.

REM Render the main static presentation (no Shiny)
call C:\ProgramData\anaconda3\Library\bin\conda.bat run -n quarto python create_qr.py

if errorlevel 1 (
    echo Generation failed
    pause
    exit /b 1
)

echo.
echo QR code generated successfully!
echo.
pause
