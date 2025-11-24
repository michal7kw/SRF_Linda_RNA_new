#!/bin/bash
#SBATCH --job-name=plot_cell_type_L1_new_markers_all_methods
#SBATCH --output=./logs/plot_cell_type_L1_new_markers_all_methods_%a.out
#SBATCH --error=./logs/plot_cell_type_L1_new_markers_all_methods_%a.err
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=82G
#SBATCH --time=120:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Master script to generate all visualization types for cell type L1 markers
# This script runs all available visualization methods based on scanpy best practices
#
# Visualization types:
# 1. Dotplot - Shows fraction expressing (size) and mean expression (color)
# 2. Matrixplot - Heatmap of mean expression per cell type
# 3. Stacked Violin - Distribution of expression across cell types
# 4. Heatmap - Individual cell expression patterns (subsampled)
#

set -e  # Exit on error

# Load conda environment
source ~/.bashrc

conda activate seurat_full2

# Set working directory to script location
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
cd "$SCRIPT_DIR"

echo "Working directory: $(pwd)"
echo ""

# Define arrays of scripts and their names (full paths)
SCRIPTS=(
    "${SCRIPT_DIR}/plot_cell_type_L1_new_markers_dotplot.py"
    "${SCRIPT_DIR}/plot_cell_type_L1_new_markers_matrixplot.py"
    "${SCRIPT_DIR}/plot_cell_type_L1_new_markers_stacked_violin.py"
    "${SCRIPT_DIR}/plot_cell_type_L1_new_markers_heatmap.py"
)

NAMES=(
    "DOTPLOTS"
    "MATRIXPLOTS"
    "STACKED VIOLINS"
    "HEATMAPS"
)

OUTPUT_DIRS=(
    "dotplots/L1/"
    "matrixplots/L1/"
    "stacked_violins/L1/"
    "heatmaps/L1/"
)

# Get the script and name for this array task
SCRIPT=${SCRIPTS[$SLURM_ARRAY_TASK_ID]}
NAME=${NAMES[$SLURM_ARRAY_TASK_ID]}
OUTPUT_DIR=${OUTPUT_DIRS[$SLURM_ARRAY_TASK_ID]}

echo "========================================================================"
echo "CELL TYPE L1 MARKER VISUALIZATION - ${NAME}"
echo "Job array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "========================================================================"
echo ""
echo "Running: ${SCRIPT}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "========================================================================"
echo ""

# Run the visualization script
python "$SCRIPT"

echo ""
echo "========================================================================"
echo "${NAME} COMPLETED"
echo "========================================================================"
echo ""
echo "Output saved to: ${OUTPUT_DIR}"
echo "Files include PNG (300 DPI) and PDF (vector) formats"
echo ""
echo "========================================================================"
