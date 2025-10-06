@echo off
REM TB Resistance Predictor - Demo Launcher (Windows)
REM Usage: run_demo.bat

echo.
echo ====================================
echo TB Resistance Predictor - Demo
echo ====================================
echo.

REM Check Python
where python >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] Python not found. Please install Python 3.8+
    pause
    exit /b 1
)

for /f "tokens=2" %%i in ('python --version') do set PYTHON_VERSION=%%i
echo [OK] Python %PYTHON_VERSION% found
echo.

REM Check venv
if not exist ".venv\" (
    echo [INFO] Creating virtual environment...
    python -m venv .venv
    echo [OK] Virtual environment created
    echo.
)

REM Activate venv
echo Activating virtual environment...
call .venv\Scripts\activate.bat

REM Check dependencies
python -c "import streamlit" 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo [INFO] Installing dependencies...
    pip install -q streamlit pandas numpy matplotlib
    echo [OK] Dependencies installed
    echo.
) else (
    echo [OK] Dependencies already installed
    echo.
)

REM Check files
if not exist "tb_resistance_mvp\drug_filtering_demo.py" (
    echo [ERROR] drug_filtering_demo.py not found
    pause
    exit /b 1
)

if not exist "ui_streamlit_demo.py" (
    echo [ERROR] ui_streamlit_demo.py not found
    pause
    exit /b 1
)

REM Launch
echo.
echo ====================================
echo Launching demo...
echo ====================================
echo.
echo The demo will open in your browser at:
echo   http://localhost:8501
echo.
echo Press Ctrl+C to stop the server
echo.

streamlit run ui_streamlit_demo.py --server.port 8501 --server.headless true

REM Cleanup
call .venv\Scripts\deactivate.bat