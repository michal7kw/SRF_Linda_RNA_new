#!/bin/bash
#SBATCH --job-name=deg_cov_enh_sep
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=logs/run_viz2_deg_coverage_enhanced_separate_L1.out
#SBATCH --error=logs/run_viz2_deg_coverage_enhanced_separate_L1.err

################################################################################
# SLURM Job Script: Enhanced Coverage Plots for DEGs (Separate Genotypes)
#
# This script creates comprehensive visualizations including:
#   - Per-cell fragment abundance plots by cell type (4 groups)
#   - Coverage plots with ATAC + RNA
#   - Peak-gene linkage computation (genotype-stratified)
#   - Genomic links visualization
#
# Resource requirements:
#   - 8 CPUs for parallel processing
#   - 64GB RAM for Seurat + linkage computation
#   - 8 hours (linkage computation is time-consuming)
#
# Usage:
#   sbatch run_deg_coverage_enhanced_separate.sh
################################################################################

set -e
set -u
set -o pipefail

mkdir -p logs

echo "================================================================================"
echo "ENHANCED DEG COVERAGE PLOTS (SEPARATE GENOTYPES)"
echo "================================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "Running on: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 64GB"
echo "================================================================================"
echo ""

# Activate conda environment
echo "Activating conda environment: sc-chromatin2"
source ~/.bashrc
conda activate sc-chromatin2

echo ""
echo "Environment check:"
echo "  R version: $(R --version | head -1)"
echo "  Conda env: $CONDA_DEFAULT_ENV"
echo ""

# Set reticulate Python
export RETICULATE_PYTHON="/beegfs/scratch/ric.broccoli/kubacki.michal/conda/envs/sc-chromatin2/bin/python"
export RETICULATE_CONDA="/beegfs/scratch/ric.broccoli/kubacki.michal/conda/condabin/conda"

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline

echo "================================================================================"
echo "CREATING ENHANCED COVERAGE PLOTS WITH GENOMIC LINKS"
echo "================================================================================"
echo ""

Rscript scripts/separate_version/create_coverage_plots_from_DEGs_enhanced_separate_L1.R

if [ $? -ne 0 ]; then
    echo ""
    echo "✗ Error: Enhanced coverage plot generation failed!"
    exit 1
fi

echo ""
echo "✓ Enhanced coverage plot generation completed successfully"
echo ""

echo "================================================================================"
echo "WORKFLOW COMPLETE"
echo "================================================================================"
echo ""
echo "Finished: $(date)"
echo ""

RESULTS_DIR="signac_results_L1/DEG_coverage_plots_enhanced_separate"

echo "Output directory: ${RESULTS_DIR}/"
echo ""

if [ -d "$RESULTS_DIR" ]; then
    echo "Fragment abundance plots:"
    n_fragment=$(ls ${RESULTS_DIR}/fragment_abundance/*.pdf 2>/dev/null | wc -l)
    echo "  Cell types: $n_fragment"

    echo ""
    echo "Peak-gene linkage files:"
    n_links=$(ls ${RESULTS_DIR}/peak_gene_links/*.csv 2>/dev/null | wc -l)
    echo "  Genes with links: $n_links"

    echo ""
    echo "Coverage plots by cell type:"
    for dir in ${RESULTS_DIR}/*/; do
        if [ -d "$dir" ] && [[ ! "$dir" =~ (fragment_abundance|peak_gene_links) ]]; then
            n_plots=$(ls ${dir}*.pdf 2>/dev/null | wc -l)
            dirname=$(basename "$dir")
            echo "  ${dirname}: $n_plots plots"
        fi
    done

    echo ""
    total_plots=$(find ${RESULTS_DIR} -name "*coverage*.pdf" 2>/dev/null | wc -l)
    echo "Total coverage plots: $total_plots"
fi

echo ""
echo "✓ All tasks completed successfully!"
echo "================================================================================"
