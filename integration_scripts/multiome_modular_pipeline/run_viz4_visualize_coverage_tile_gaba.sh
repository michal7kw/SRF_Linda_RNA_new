#!/bin/bash
#SBATCH --job-name=viz4_gaba
#SBATCH --partition=workq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/run_viz4_visualize_coverage_tile_gaba.out
#SBATCH --error=logs/run_viz4_visualize_coverage_tile_gaba.err

################################################################################
# SLURM Wrapper: Coverage and Tile Visualization for GABA Neurons (Combined)
#
# This script runs visualize_coverage_tile_gaba.R via SLURM
# Creates CoveragePlot (peaks=FALSE) and TilePlot for all genes in genes_inter.txt
################################################################################

set -euo pipefail

# Source conda
source ~/.bashrc
conda activate sc-chromatin2

# Set R environment variables
export RETICULATE_PYTHON="/beegfs/scratch/ric.broccoli/kubacki.michal/conda/envs/sc-chromatin2/bin/python"
export RETICULATE_CONDA="/beegfs/scratch/ric.broccoli/kubacki.michal/conda/bin/conda"

# Log start
echo "========================================"
echo "Coverage/Tile Visualization: GABA (Combined Controls)"
echo "Started: $(date)"
echo "========================================"
echo ""

# Show environment
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "R version:"
Rscript --version
echo ""

# Verify required packages
echo "Checking R packages..."
Rscript -e "
  required_pkgs <- c('Signac', 'Seurat', 'GenomicRanges', 'BSgenome.Mmusculus.UCSC.mm10',
                     'EnsDb.Mmusculus.v79', 'ggplot2', 'dplyr', 'tidyr', 'patchwork')
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf('Required package %s is not installed', pkg))
    }
  }
  cat('✓ All required packages available\n')
"
echo ""

# Navigate to correct directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline

# Run the R script
echo "Running visualization script..."
Rscript scripts/visualize_coverage_tile_gaba.R

# Check exit status
if [ $? -eq 0 ]; then
  echo ""
  echo "========================================"
  echo "✓ Visualization completed successfully"
  echo "Completed: $(date)"
  echo "========================================"
else
  echo ""
  echo "========================================"
  echo "✗ Visualization failed"
  echo "Failed: $(date)"
  echo "========================================"
  exit 1
fi
