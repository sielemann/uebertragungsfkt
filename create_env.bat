@echo off

rem Create the environment
conda create --name quarto python

rem Activate the environment
conda activate quarto

rem Run installer from https://quarto.org/docs/get-started/ (do NOT get quarto-cli package via pip!)

rem Install packages per requirements.txt
rem pip freeze > requirements.txt
rem pip install -r requirement.txt

rem In VS Code, to select the quarto env and preview, press Ctrl + Shift + P to bring up "Command Palette", and select "Select Python interpreter". Then, select the Conda env "quarto"