@echo off
echo Starting standalone Shiny server (NOT REQUIRED!)...
echo.
echo Server will run on http://localhost:8000
echo Press Ctrl+C to stop.
echo.

C:\ProgramData\anaconda3\Library\bin\conda.bat run -n quarto quarto serve vortrag.qmd --port 8000

pause
