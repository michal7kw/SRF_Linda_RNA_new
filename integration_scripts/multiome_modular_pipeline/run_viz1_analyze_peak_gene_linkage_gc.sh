#!/bin/bash
#SBATCH --job-name=viz1_gc
#SBATCH --partition=workq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/run_viz1_analyze_peak_gene_linkage_gc.out
#SBATCH --error=logs/run_viz1_analyze_peak_gene_linkage_gc.err

################################################################################
# Run Peak-Gene Linkage Analysis for Granule Cells (Combined Controls - UPDATED)
################################################################################

cd "$(dirname "$0")/.." || exit 1

echo "================================================================================"
echo "PEAK-GENE LINKAGE ANALYSIS: GRANULE CELLS (UPDATED VERSION)"
echo "================================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "================================================================================"
echo

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline

source ~/.bashrc
conda activate sc-chromatin2

echo "Running analysis..."
Rscript scripts/analyze_peak_gene_linkage_gc_UPDATED.R

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
  echo "✓ Analysis completed successfully!"
  echo "Results: signac_results_L1/peak_gene_linkage_analysis_UPDATED/GC/"
else
  echo "✗ Analysis failed with exit code: $EXIT_CODE"
fi

exit $EXIT_CODE
