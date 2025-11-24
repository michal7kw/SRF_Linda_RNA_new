#!/bin/bash
#SBATCH --job-name=coverage_plots_lit_CREs
#SBATCH --output=logs/coverage_plots_literature_CREs.log
#SBATCH --error=logs/coverage_plots_literature_CREs.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Run Coverage Plot Generation with Literature CREs
#
# This script creates coverage plots for top GABA DEGs showing:
# - ATAC coverage tracks (Ctrl vs Mut)
# - Gene expression
# - Genomic links from literature CREs to gene TSS
# - Gene annotations
#
# INPUT FILES:
# - CREs_from_literature/create_coverage_plots_with_literature_CREs.R: R script to generate coverage plots
# - analyze_GABA_DEGs_with_literature_CREs.py: Python script (prerequisite)
# - extract_bigwig_signals_and_visualize.py: Python script (prerequisite)
# - BigWig files with ATAC-seq signal data
# - Gene expression data
# - Literature CRE data
#
# OUTPUT FILES:
# - Coverage plots for top GABA DEGs (PNG format)
# - logs/coverage_plots_literature_CREs.log: Execution log file
# - logs/coverage_plots_literature_CREs.err: Error log file
#
################################################################################

# Set environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate sc-chromatin2

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline

# Create logs directory
mkdir -p CREs_from_literature/logs

# Set R environment variables
export RETICULATE_PYTHON=/beegfs/scratch/ric.broccoli/kubacki.michal/conda/envs/sc-chromatin2/bin/python
export RETICULATE_CONDA=/beegfs/scratch/ric.broccoli/kubacki.michal/conda

# Run R script
echo "Starting coverage plot generation..."
echo "Time: $(date)"
echo "================================"

Rscript CREs_from_literature/create_coverage_plots_with_literature_CREs.R

EXIT_CODE=$?

echo "================================"
echo "Time: $(date)"

if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Coverage plots completed successfully"
else
    echo "✗ Coverage plots failed with exit code $EXIT_CODE"
fi

exit $EXIT_CODE
