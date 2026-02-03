#!/bin/bash
#
# Launch Jupyter Notebook for Interactive Parser Tutorial
#

echo "========================================="
echo "Launching Jupyter Notebook"
echo "========================================="
echo ""

# Activate conda environment
source /shared/software/miniconda3/latest/etc/profile.d/conda.sh
conda activate opfi

# Check if jupyter is installed
if ! command -v jupyter &> /dev/null; then
    echo "Installing jupyter..."
    pip install jupyter notebook
fi

# Launch jupyter in the project directory
echo "Starting Jupyter Notebook..."
echo "The notebook will open in your browser"
echo ""
echo "To stop the server, press Ctrl+C"
echo ""

jupyter notebook Interactive_Parser_Tutorial.ipynb

# Alternative: Use JupyterLab
# jupyter lab Interactive_Parser_Tutorial.ipynb
