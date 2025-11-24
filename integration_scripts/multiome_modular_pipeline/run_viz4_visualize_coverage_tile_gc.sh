#!/bin/bash
#SBATCH --job-name=viz4_gc
#SBATCH --partition=workq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/run_viz4_visualize_coverage_tile_gc.out
#SBATCH --error=logs/run_viz4_visualize_coverage_tile_gc.err

################################################################################
# SLURM Wrapper: Coverage and Tile Visualization for Granule Cells (Combined)
#
# This script runs visualize_coverage_tile_gc.R via SLURM
# Creates CoveragePlot (peaks=FALSE) and TilePlot for all genes in genes_gc.txt
################################################################################

set -euo pipefail

source ~/.bashrc
conda activate sc-chromatin2

export RETICULATE_PYTHON="/beegfs/scratch/ric.broccoli/kubacki.michal/conda/envs/sc-chromatin2/bin/python"
export RETICULATE_CONDA="/beegfs/scratch/ric.broccoli/kubacki.michal/conda/bin/conda"

echo "========================================"
echo "Coverage/Tile Visualization: Granule Cells (Combined Controls)"
echo "Started: $(date)"
echo "========================================"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline

Rscript scripts/visualize_coverage_tile_gc.R

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
