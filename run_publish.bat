@echo off
echo Running Python command...
C:\ProgramData\anaconda3\Library\bin\conda.bat run -n quarto python -m publish
rem C:\ProgramData\anaconda3\Library\bin\conda.bat run -n quarto python -m pytest tests -vv
rem C:\ProgramData\anaconda3\Library\bin\conda.bat run -n quarto python -m pytest -v -m tf_test tests/
if errorlevel 1 (
    echo Test execution failed
    pause
    exit /b 1
)
echo All tests completed
pause
