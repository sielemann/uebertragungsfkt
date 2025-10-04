@echo off
echo Starting Quarto Shiny app...
echo.
echo This will open a browser with the interactive simulation.
echo Press Ctrl+C to stop the server.
echo.

C:\ProgramData\anaconda3\Library\bin\conda.bat run -n quarto quarto preview 02-simulation.qmd

pause
