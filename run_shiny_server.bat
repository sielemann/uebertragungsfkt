@echo off
echo Starting standalone Shiny server...
echo.
echo Server will run on http://localhost:8000
echo Press Ctrl+C to stop.
echo.

C:\ProgramData\anaconda3\Library\bin\conda.bat run -n quarto quarto serve 02-simulation.qmd --port 8000

pause
