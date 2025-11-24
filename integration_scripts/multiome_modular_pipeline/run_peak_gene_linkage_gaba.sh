#!/bin/bash
#SBATCH --job-name=gaba_linkage
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/run_peak_gene_linkage_gaba.out
#SBATCH --error=logs/run_peak_gene_linkage_gaba.err

################################################################################
# SLURM Job Script: Peak-Gene Linkage Analysis for GABA Neurons
#
# This script computes peak-gene correlations for genes in genes_inter.txt.
# Visualization step has been separated into run_visualization_gaba.sh.
#
# Resource requirements:
#   - 8 CPUs for parallel processing
#   - 64GB RAM for Seurat object and correlation computations
#   - 4 hours runtime (should complete in ~1-2 hours)
#
# Usage:
#   sbatch run_peak_gene_linkage_gaba.sh      # Linkage analysis only
#   sbatch run_visualization_gaba.sh          # Visualization (after linkage)
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit if any command in pipe fails

# Create logs directory if it doesn't exist
mkdir -p logs

echo "================================================================================"
echo "PEAK-GENE LINKAGE ANALYSIS: GABA NEURONS"
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
# Step 1: Linkage Analysis
################################################################################

echo "================================================================================"
echo "STEP 1: PEAK-GENE LINKAGE ANALYSIS"
echo "================================================================================"
echo ""

Rscript ./scripts/analyze_peak_gene_linkage_gaba_UPDATED.R

if [ $? -ne 0 ]; then
    echo ""
    echo "✗ Error: Linkage analysis failed!"
    exit 1
fi

echo ""
echo "✓ Linkage analysis completed successfully"
echo ""

################################################################################
# Summary
################################################################################

echo "================================================================================"
echo "LINKAGE ANALYSIS COMPLETE"
echo "================================================================================"
echo ""
echo "Finished: $(date)"
echo ""

# Print output file locations
RESULTS_DIR="signac_results_L1/peak_gene_linkage_analysis/GABA"

echo "Output files:"
echo "  Linkage results: ${RESULTS_DIR}/"
echo "  Summary statistics: ${RESULTS_DIR}/peak_gene_linkage_summary.csv"
echo ""

# Count output files
if [ -d "$RESULTS_DIR" ]; then
    n_link_files=$(ls ${RESULTS_DIR}/*_peak_gene_links.csv 2>/dev/null | wc -l)
    echo "  Generated linkage files: $n_link_files"
fi

echo ""
echo "✓ Linkage analysis completed successfully!"
echo ""
echo "Next steps:"
echo "  - To create visualizations, run: sbatch run_visualization_gaba.sh"
echo "================================================================================"
