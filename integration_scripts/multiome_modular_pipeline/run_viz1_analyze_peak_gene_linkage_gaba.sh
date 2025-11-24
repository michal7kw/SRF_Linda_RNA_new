#!/bin/bash
#SBATCH --job-name=viz1_gaba
#SBATCH --partition=workq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/run_viz1_analyze_peak_gene_linkage_gaba.out
#SBATCH --error=logs/run_viz1_analyze_peak_gene_linkage_gaba.err

################################################################################
# Run Peak-Gene Linkage Analysis for GABA Neurons (Combined Controls - UPDATED)
#
# This wrapper script runs the updated peak-gene linkage analysis with:
# - More permissive threshold (r >= 0.15)
# - Mitochondrial gene filtering
# - Coverage plots with genomic links
# - Fragment abundance visualization
################################################################################

cd "$(dirname "$0")/.." || exit 1

echo "================================================================================"
echo "PEAK-GENE LINKAGE ANALYSIS: GABA NEURONS (UPDATED VERSION)"
echo "================================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "Running on: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 64GB"
echo "================================================================================"
echo

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline

# Activate conda environment
echo "Activating conda environment: sc-chromatin2"
source ~/.bashrc
conda activate sc-chromatin2

# Verify ggforce is available
echo "Verifying ggforce package..."
Rscript -e "library(ggforce); cat('  ✓ ggforce version:', as.character(packageVersion('ggforce')), '\n')" || {
  echo "ERROR: ggforce package not found!"
  echo "Please install: install.packages('ggforce')"
  exit 1
}

echo
echo "Environment check:"
echo "  R version: $(R --version | head -1)"
echo "  Conda env: $CONDA_DEFAULT_ENV"
echo

echo "================================================================================"
echo "RUNNING PEAK-GENE LINKAGE ANALYSIS"
echo "================================================================================"
echo

# Run the analysis
Rscript scripts/analyze_peak_gene_linkage_gaba_UPDATED.R

EXIT_CODE=$?

echo
echo "================================================================================"
echo "WORKFLOW COMPLETE"
echo "================================================================================"
echo "Finished: $(date)"
echo "Exit code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
  echo
  echo "✓ Analysis completed successfully!"
  echo
  echo "Results saved to: signac_results_L1/peak_gene_linkage_analysis_UPDATED/GABA/"
  echo
  echo "Output files:"
  echo "  - Peak-gene links: peak_gene_links/*.csv"
  echo "  - Coverage plots: coverage_plots/*_coverage_with_links.pdf"
  echo "  - Fragment abundance: fragment_abundance/GABA_fragment_abundance.pdf"
  echo "  - Summary: peak_gene_linkage_summary.csv"
else
  echo
  echo "✗ Analysis failed with exit code: $EXIT_CODE"
  echo "Check error log for details"
fi

echo "================================================================================"

exit $EXIT_CODE
