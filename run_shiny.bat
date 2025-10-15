@echo off
echo Starting Quarto Shiny app (NOT REQUIRED!)...
echo.
echo This will open a browser with the interactive simulation.
echo Press Ctrl+C to stop the server.
echo.

C:\ProgramData\anaconda3\Library\bin\conda.bat run -n quarto quarto preview vortrag.qmd

pause
