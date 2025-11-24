#!/bin/bash
#SBATCH --job-name=plot_cell_type_L2_new_markers_dotplot
#SBATCH --output=./logs/plot_cell_type_L2_new_markers_dotplot.out
#SBATCH --error=./logs/plot_cell_type_L2_new_markers_dotplot.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=82G
#SBATCH --time=120:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Set strict error handling
set -euo pipefail

# Base directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"

# Change to working directory
cd "$BASE_DIR"

# Load conda environment
source ~/.bashrc

conda activate seurat_full2

python plot_cell_type_check_layers.py
python plot_cell_type_L2_new_markers_dotplot.py