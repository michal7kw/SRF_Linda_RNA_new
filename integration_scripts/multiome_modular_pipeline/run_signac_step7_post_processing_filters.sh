#!/bin/bash

################################################################################
# SLURM Submission Script for Post-Processing Filters
################################################################################

#SBATCH --job-name=run_signac_step7_post_processing_filters
#SBATCH --output=logs/run_signac_step7_post_processing_filters.out
#SBATCH --error=logs/run_signac_step7_post_processing_filters.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --partition=workq
#SBATCH --account=kubacki.michal

################################################################################
# Configuration
################################################################################

# Base directory
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"

# Create logs directory
mkdir -p ${BASE_DIR}/logs

# Load conda environment
source ~/.bashrc
conda activate sc-chromatin2

# Set number of threads
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "=================================================================================="
echo "  POST-PROCESSING FILTERS"
echo "=================================================================================="
echo ""
echo "Applies recommended filters to L1 results:"
echo "  1. Filters extreme log2FC values in DA peaks"
echo "  2. Increases correlation threshold for peak-gene links (r > 0.2)"
echo "  3. Applies FDR correction to correlations"
echo ""
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo "Start time: $(date)"

echo ""
echo "This script applies publication-quality filters to:"
echo "  1. Differential Accessibility (DA) results"
echo "  2. Peak-Gene Linkage results"
echo ""
echo "Filters applied:"
echo "  - Remove extreme log2FC values (|log2FC| > 8)"
echo "  - Remove binary ON/OFF peaks (min.pct < 20%)"
echo "  - Apply stricter correlation threshold (r > 0.2)"
echo "  - Apply FDR correction (q < 0.05)"
echo ""
echo "================================================================================"
echo ""

################################################################################
# Run Post-Processing
################################################################################

START_TIME=$(date +%s)

Rscript ${BASE_DIR}/signac_07_POST_PROCESSING_FILTERS.R 2>&1

EXIT_CODE=$?
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "✗ Post-processing failed with exit code ${EXIT_CODE}"
    exit 1
fi

echo ""
echo "✓ Post-processing filters complete (elapsed time: ${ELAPSED}s)"
echo ""

################################################################################
# Summary
################################################################################

echo "================================================================================"
echo "POST-PROCESSING COMPLETE!"
echo "================================================================================"
echo ""
echo "End time: $(date)"
echo ""
echo "High-confidence filtered results saved to:"
echo "  ${BASE_DIR}/signac_results/celltype_results/filtered/"
echo ""
echo "Files created:"
echo "  - signac_results/celltype_results/filtered/*_DA_peaks.csv"
echo "  - signac_results/celltype_results/filtered/*_peak_gene_links.csv"
echo "  - signac_results/celltype_results/filtered/DA_filtering_summary.csv"
echo "  - signac_results/celltype_results/filtered/links_filtering_summary.csv"
echo "  - signac_results/celltype_results/filtered/DA_filtering_comparison.pdf"
echo "  - signac_results/celltype_results/filtered/correlation_distribution.pdf"
echo ""
echo "Note: DEG files are NOT filtered by this script (only DA peaks and peak-gene links)"
echo ""
echo "Recommendation: Use filtered/ results for publication"
echo ""
echo "================================================================================"
