@echo off

where conda >nul 2>nul 
if %errorlevel% == 0 (
    echo Conda is already installed
) else (
    :: Specify Miniconda version
    set MINICONDA_VERSION=Miniconda3-latest
    set MINICONDA_SH=%MINICONDA_VERSION%-Windows-x86_64.exe

    :: Download Miniconda
    powershell -Command "(New-Object System.Net.WebClient).DownloadFile('https://repo.anaconda.com/miniconda/%MINICONDA_SH%', '%MINICONDA_SH%')"

    :: Install Miniconda
    start /wait "" %MINICONDA_SH% /InstallationType=JustMe /RegisterPython=0 /S /D=%UserProfile%\Miniconda3

    :: Initialize Conda
    call %UserProfile%\Miniconda3\Scripts\activate.bat
    call conda init

    echo Miniconda installed.
)

:: ----------------------------------------------------------------------------
:: Activate Miniconda
call %UserProfile%\Miniconda3\Scripts\activate.bat

:: Install Mamba from Conda-Forge
call conda install mamba -n base -c conda-forge

:: ----------------------------------------------------------------------------
:: Activate conda
call %UserProfile%\Miniconda3\Scripts\activate.bat

call conda install -c conda-forge m2w64-gcc
call conda install -c conda-forge imagemagick
call conda install -c msys2 m2w64-gcc-fortran

:: Create conda environment from YAML
call conda env create -f environment.yml

echo Done! Type conda activate seidart to access the seidart environment.
