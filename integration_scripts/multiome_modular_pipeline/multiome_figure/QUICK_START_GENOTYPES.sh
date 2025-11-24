#!/bin/bash

#SBATCH --job-name=QUICK_START_GENOTYPES
#SBATCH --output=logs/QUICK_START_GENOTYPES.out
#SBATCH --error=logs/QUICK_START_GENOTYPES.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=workq
#SBATCH --account=kubacki.michal

# Load conda environment
source ~/.bashrc
conda activate r-complex-vis 

################################################################################
# QUICK START: Multiome Figure Generation (Genotype-Stratified)
################################################################################
#
# This script runs the complete multiome figure generation pipeline for
# BOTH Nestin and Emx1 genotypes automatically.
#
# Prerequisites:
#   1. Signac pipeline Steps 1-3 completed
#   2. Results in signac_results/celltype_results/
#   3. Conda environment with R, Seurat, Signac, ComplexHeatmap
#
# Usage:
#   bash QUICK_START_GENOTYPES.sh
#
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

echo "================================================================================"
echo "  MULTIOME FIGURE GENERATION - GENOTYPE-STRATIFIED ANALYSIS"
echo "================================================================================"
echo ""

# Change to multiome_figure directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/multiome_figure"
cd "$SCRIPT_DIR"

echo "Working directory: $SCRIPT_DIR"
echo ""

# Check prerequisites
echo "=== Checking prerequisites ==="
echo ""

# Check if Seurat object exists
SEURAT_OBJ="../signac_results/integrated_seurat_processed.rds"
if [ ! -f "$SEURAT_OBJ" ]; then
    echo "✗ Error: Seurat object not found: $SEURAT_OBJ"
    echo "  Please run signac_02_process_and_reduce.R first"
    exit 1
fi
echo "✓ Found Seurat object"

# Check if DEG results exist
DEG_DIR="../signac_results/celltype_results/DEG"
if [ ! -d "$DEG_DIR" ] || [ -z "$(ls -A $DEG_DIR 2>/dev/null)" ]; then
    echo "✗ Error: No DEG results found in $DEG_DIR"
    echo "  Please run signac_03_celltype_analysis.R first"
    exit 1
fi
echo "✓ Found DEG results"

# Count result files
N_DEG=$(ls -1 $DEG_DIR/*.csv 2>/dev/null | wc -l)
N_DA=$(ls -1 ../signac_results/celltype_results/DA/*.csv 2>/dev/null | wc -l)
N_LINKS=$(ls -1 ../signac_results/celltype_results/peak_gene_links/*.csv 2>/dev/null | wc -l)

echo "  DEG files: $N_DEG"
echo "  DA files: $N_DA"
echo "  Link files: $N_LINKS"
echo ""

if [ "$N_DEG" -eq 0 ] || [ "$N_DA" -eq 0 ] || [ "$N_LINKS" -eq 0 ]; then
    echo "✗ Error: Missing result files"
    echo "  Please run the complete Signac pipeline (Steps 1-3)"
    exit 1
fi

# Activate conda environment (if needed)
if [ -n "${CONDA_DEFAULT_ENV:-}" ]; then
    echo "✓ Conda environment active: $CONDA_DEFAULT_ENV"
else
    echo "Note: No conda environment detected. Assuming R packages are available."
fi
echo ""

################################################################################
# Step 1: Prepare Data
################################################################################

echo "================================================================================"
echo "STEP 1: Prepare Multiome Heatmap Data"
echo "================================================================================"
echo ""
echo "This step:"
echo "  - Identifies concordant chromatin-driven DEGs for each genotype"
echo "  - Selects top cell types by abundance"
echo "  - Extracts RNA expression and ATAC gene activity signals"
echo "  - Calculates per-cell-type, per-condition averages"
echo ""
echo "Expected runtime: ~10-15 minutes"
echo ""

Rscript prepare_multiome_heatmap_data_genotypes.R

if [ $? -ne 0 ]; then
    echo ""
    echo "✗ Step 1 failed. Check error messages above."
    exit 1
fi

echo ""
echo "✓ Step 1 complete"
echo ""

################################################################################
# Step 2: Create Visualizations
################################################################################

echo "================================================================================"
echo "STEP 2: Create Heatmap Visualizations"
echo "================================================================================"
echo ""
echo "This step:"
echo "  - Creates separate heatmaps for each genotype"
echo "  - Shows RNA expression and ATAC gene activity side-by-side"
echo "  - Saves publication-quality PDF and PNG files"
echo ""
echo "Expected runtime: ~2-5 minutes"
echo ""

Rscript plot_multiome_heatmap_genotypes.R

if [ $? -ne 0 ]; then
    echo ""
    echo "✗ Step 2 failed. Check error messages above."
    exit 1
fi

echo ""
echo "✓ Step 2 complete"
echo ""

################################################################################
# Summary
################################################################################

echo "================================================================================"
echo "✓ MULTIOME FIGURE GENERATION COMPLETE!"
echo "================================================================================"
echo ""
echo "Output location: $SCRIPT_DIR/output_genotypes/"
echo ""
echo "Files created:"
echo "  For each genotype (Nestin, Emx1):"
echo "    - multiome_heatmap_<genotype>.pdf (publication-quality)"
echo "    - multiome_heatmap_<genotype>.png (presentation-ready)"
echo "    - chromatin_driven_degs_<genotype>.csv (gene list with statistics)"
echo ""
echo "Data files (for custom analysis):"
echo "    - multiome_heatmap_data_<genotype>.rds (R object with matrices)"
echo ""

# List actual output files
if [ -d "output_genotypes" ]; then
    echo "Generated files:"
    ls -lh output_genotypes/*.pdf output_genotypes/*.csv 2>/dev/null | awk '{print "    " $9 " (" $5 ")"}'
fi

echo ""
echo "Next steps:"
echo "  1. Review heatmaps in output_genotypes/"
echo "  2. Check chromatin_driven_degs_*.csv for gene lists"
echo "  3. Use these figures in your publication/presentation"
echo ""
echo "For custom analysis, load the .rds files in R:"
echo "  data <- readRDS('multiome_heatmap_data_Nestin.rds')"
echo ""
