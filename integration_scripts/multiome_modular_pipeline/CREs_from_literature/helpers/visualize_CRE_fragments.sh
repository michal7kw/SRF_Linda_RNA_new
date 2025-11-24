#!/bin/bash
#SBATCH --job-name=visualize_CREs
#SBATCH --output=logs/run_visualize_CRE_fragments.out
#SBATCH --error=logs/run_visualize_CRE_fragments.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Visualize ATAC Coverage at Literature CRE Locations
#
# This script visualizes ATAC-seq coverage at CRE locations from literature
# to analyze chromatin accessibility patterns.
#
# DEPENDENCIES:
# - R script called: visualize_CRE_fragments.R (main visualization script)
# - Conda environment: sc-chromatin2 (contains R packages and reticulate)
# - Reticulate Python configuration for R-Python integration
# - ATAC-seq fragment data from previous pipeline steps
# - Literature CRE coordinates from external sources
# - Gene annotation data from reference genomes
#
# INPUT FILES:
# - visualize_CRE_fragments.R: R script to visualize CRE coverage
# - ATAC-seq fragment data
# - Literature CRE coordinates
# - Gene annotation data
#
# OUTPUT FILES:
# - Coverage visualization plots (PNG format)
# - logs/run_visualize_CRE_fragments.out: Execution log file
# - logs/run_visualize_CRE_fragments.err: Error log file
#
################################################################################

# Create logs directory
mkdir -p logs

echo "========================================"
echo "Visualize ATAC Coverage at Literature CRE Locations"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate sc-chromatin2

# Set reticulate Python
export RETICULATE_PYTHON="/beegfs/scratch/ric.broccoli/kubacki.michal/conda/envs/sc-chromatin2/bin/python"
export RETICULATE_CONDA="/beegfs/scratch/ric.broccoli/kubacki.michal/conda/bin/conda"

echo "Environment:"
echo "  Conda env: sc-chromatin2"
echo "  R version: $(Rscript --version 2>&1)"
echo "  Working directory: $(pwd)"
echo ""

# Run visualization script
echo "Running fragment coverage visualization..."
Rscript visualize_CRE_fragments.R

EXIT_CODE=$?

echo ""
echo "========================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Visualization completed successfully"
else
    echo "✗ Visualization failed with exit code: $EXIT_CODE"
fi
echo "Finished: $(date)"
echo "========================================"

exit $EXIT_CODE
