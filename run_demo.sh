#!/bin/bash

# TB Resistance Predictor - Demo Launcher
# Usage: ./run_demo.sh

set -e  # Exit on error

echo "🧬 TB Resistance Predictor - Demo Version"
echo "=========================================="
echo ""

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check Python
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}❌ Python 3 not found. Please install Python 3.8+${NC}"
    exit 1
fi

PYTHON_VERSION=$(python3 --version | cut -d " " -f 2)
echo -e "${GREEN}✅ Python ${PYTHON_VERSION} found${NC}"

# Check if venv exists
if [ ! -d ".venv" ]; then
    echo ""
    echo -e "${YELLOW}⚠️  Virtual environment not found. Creating...${NC}"
    python3 -m venv .venv
    echo -e "${GREEN}✅ Virtual environment created${NC}"
fi

# Activate venv
echo ""
echo "Activating virtual environment..."
source .venv/bin/activate

# Check if dependencies installed
if ! python -c "import streamlit" 2>/dev/null; then
    echo ""
    echo -e "${YELLOW}⚠️  Dependencies not found. Installing...${NC}"
    pip install -q streamlit pandas numpy matplotlib
    echo -e "${GREEN}✅ Dependencies installed${NC}"
else
    echo -e "${GREEN}✅ Dependencies already installed${NC}"
fi

# Run tests (optional)
if [ "$1" == "--test" ]; then
    echo ""
    echo "Running tests..."
    python test_demo.py
    echo ""
fi

# Check if demo files exist
if [ ! -f "tb_resistance_mvp/drug_filtering_demo.py" ]; then
    echo -e "${RED}❌ drug_filtering_demo.py not found${NC}"
    echo "Please ensure you have all demo files in tb_resistance_mvp/"
    exit 1
fi

if [ ! -f "ui_streamlit_demo.py" ]; then
    echo -e "${RED}❌ ui_streamlit_demo.py not found${NC}"
    exit 1
fi

# Launch Streamlit
echo ""
echo "=========================================="
echo -e "${GREEN}🚀 Launching demo...${NC}"
echo "=========================================="
echo ""
echo "The demo will open in your browser at:"
echo "  👉 http://localhost:8501"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

streamlit run ui_streamlit_demo.py --server.port 8501 --server.headless true

# Cleanup on exit
deactivate