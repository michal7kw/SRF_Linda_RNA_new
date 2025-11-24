#!/bin/bash
#SBATCH --job-name=linkage_gc_sep
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/linkage_gc_separate.out
#SBATCH --error=logs/linkage_gc_separate.err

################################################################################
# SLURM Job Script: Peak-Gene Linkage Analysis for Granule Cells (Separate Genotypes)
#
# This script runs the complete peak-gene linkage workflow for granule cells
# with separate analysis for each genotype (no combined controls):
# 1. Computes peak-gene correlations for genes in genes_gc.txt
# 2. Analyzes Nestin and Emx1 separately (Ctrl vs Mut within each genotype)
# 3. Creates visualization plots (coverage, fragment abundance, scatter)
#
# Resource requirements:
#   - 8 CPUs for parallel processing
#   - 64GB RAM for Seurat object and correlation computations
#   - 4 hours runtime (should complete in ~1-2 hours)
#
# Usage:
#   sbatch run_peak_gene_linkage_gc_separate.sh
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit if any command in pipe fails

# Create logs directory if it doesn't exist
mkdir -p logs

echo "================================================================================"
echo "PEAK-GENE LINKAGE ANALYSIS: GRANULE CELLS (SEPARATE GENOTYPES)"
echo "================================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "Running on: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 64GB"
echo "================================================================================"
echo ""

# Activate conda environment
echo "Activating conda environment: scatac_analysis"
source ~/.bashrc
conda activate sc-chromatin2

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
# Step 1: Linkage Analysis (Separate Genotypes)
################################################################################

echo "================================================================================"
echo "STEP 1: PEAK-GENE LINKAGE ANALYSIS (SEPARATE GENOTYPES)"
echo "================================================================================"
echo ""

Rscript ./scripts/analyze_peak_gene_linkage_gc_separate_UPDATED.R

if [ $? -ne 0 ]; then
    echo ""
    echo "✗ Error: Linkage analysis failed!"
    exit 1
fi

echo ""
echo "✓ Linkage analysis completed successfully"
echo ""

################################################################################
# Step 2: Visualization
################################################################################

echo "================================================================================"
echo "STEP 2: VISUALIZATION"
echo "================================================================================"
echo ""

Rscript visualize_peak_gene_linkage_gc_separate.R

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
echo "WORKFLOW COMPLETE"
echo "================================================================================"
echo ""
echo "Finished: $(date)"
echo ""

# Print output file locations
RESULTS_DIR="signac_results_L1/peak_gene_linkage_analysis_separate/GC"

echo "Output files:"
echo "  Linkage results: ${RESULTS_DIR}/"
echo "  Summary statistics: ${RESULTS_DIR}/peak_gene_linkage_summary_separate.csv"
echo "  By genotype: ${RESULTS_DIR}/by_genotype/"
echo "  Visualizations: ${RESULTS_DIR}/visualizations/"
echo ""

# Count output files
if [ -d "$RESULTS_DIR" ]; then
    if [ -d "${RESULTS_DIR}/by_genotype" ]; then
        n_link_files=$(ls ${RESULTS_DIR}/by_genotype/*_peak_gene_links.csv 2>/dev/null | wc -l)
        echo "  Generated linkage files (by genotype): $n_link_files"
    fi

    if [ -d "${RESULTS_DIR}/visualizations" ]; then
        n_viz_files=$(ls ${RESULTS_DIR}/visualizations/*.pdf 2>/dev/null | wc -l)
        echo "  Generated visualization files: $n_viz_files"
    fi
fi

echo ""
echo "✓ All tasks completed successfully!"
echo "================================================================================"
