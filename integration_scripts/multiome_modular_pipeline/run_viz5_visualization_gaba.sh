#!/bin/bash
#SBATCH --job-name=viz5_gaba
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/run_viz5_visualization_gaba.out
#SBATCH --error=logs/run_viz5_visualization_gaba.err

################################################################################
# SLURM Job Script: Peak-Gene Linkage Visualization for GABA Neurons
#
# This script creates visualization plots from existing linkage analysis results:
# - Coverage plots with genomic arcs showing peak-gene links
# - Fragment abundance QC plots
# - Correlation scatter plots
#
# Prerequisites:
#   - Must run analyze_peak_gene_linkage_gaba_UPDATED.R first --> run_peak_gene_linkage_gaba.sh
#   - Requires existing linkage CSV files in signac_results_L1/peak_gene_linkage_analysis/GABA/
#
# Resource requirements:
#   - 4 CPUs for parallel processing
#   - 32GB RAM (lighter than linkage analysis)
#   - 2 hours runtime
#
# Usage:
#   sbatch run_visualization_gaba.sh
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit if any command in pipe fails

# Create logs directory if it doesn't exist
mkdir -p logs

echo "================================================================================"
echo "PEAK-GENE LINKAGE VISUALIZATION: GABA NEURONS"
echo "================================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "Running on: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 32GB"
echo "================================================================================"
echo ""

# Activate conda environment
echo "Activating conda environment: scatac_analysis"
source ~/.bashrc
conda activate scatac_analysis

# Verify environment
echo ""
echo "Environment check:"
echo "  R version: $(R --version | head -1)"
echo "  Conda env: $CONDA_DEFAULT_ENV"
echo ""

# Set reticulate Python (for Seurat)
export RETICULATE_PYTHON="/beegfs/scratch/ric.broccoli/kubacki.michal/conda/envs/scatac_analysis/bin/python"
export RETICULATE_CONDA="/beegfs/scratch/ric.broccoli/kubacki.michal/conda/condabin/conda"

# Change to script directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline

################################################################################
# Visualization
################################################################################

echo "================================================================================"
echo "CREATING VISUALIZATION PLOTS"
echo "================================================================================"
echo ""

Rscript scripts/visualize_peak_gene_linkage_gaba.R

if [ $? -ne 0 ]; then
    echo ""
    echo "✗ Error: Visualization failed!"
    exit 1
fi

echo ""
echo "✓ Visualization completed successfully"
echo ""

################################################################################
# Summary
################################################################################

echo "================================================================================"
echo "VISUALIZATION COMPLETE"
echo "================================================================================"
echo ""
echo "Finished: $(date)"
echo ""

# Print output file locations
RESULTS_DIR="signac_results_L1/peak_gene_linkage_analysis/GABA"

echo "Output files:"
echo "  Visualizations: ${RESULTS_DIR}/visualizations/"
echo ""

# Count output files
if [ -d "${RESULTS_DIR}/visualizations" ]; then
    n_viz_files=$(ls ${RESULTS_DIR}/visualizations/*.pdf 2>/dev/null | wc -l)
    echo "  Generated visualization files: $n_viz_files"
fi

echo ""
echo "✓ All visualization tasks completed successfully!"
echo "================================================================================"
