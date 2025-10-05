@echo off
echo Publishing static Quarto presentation...
echo.

REM Render the main static presentation (no Shiny)
call C:\ProgramData\anaconda3\Library\bin\conda.bat run -n quarto quarto render 01-einleitung.qmd

if errorlevel 1 (
    echo Rendering failed
    pause
    exit /b 1
)

echo.
echo Static presentation published successfully!
echo Open 01-einleitung.html in your browser.
echo.
pause
