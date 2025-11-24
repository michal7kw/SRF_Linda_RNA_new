#!/bin/bash
#SBATCH --job-name=viz_cov_gc_sep
#SBATCH --partition=workq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/run_viz4_visualize_coverage_tile_gc_separate.out
#SBATCH --error=logs/run_viz4_visualize_coverage_tile_gc_separate.err

################################################################################
# SLURM Wrapper: Coverage and Tile Visualization for Granule Cells (Separate)
#
# This script runs visualize_coverage_tile_gc_separate.R via SLURM
# Creates plots for 4 separate groups: Nestin-Ctrl, Nestin-Mut, Emx1-Ctrl, Emx1-Mut
################################################################################

set -euo pipefail

source ~/.bashrc
conda activate sc-chromatin2

export RETICULATE_PYTHON="/beegfs/scratch/ric.broccoli/kubacki.michal/conda/envs/sc-chromatin2/bin/python"
export RETICULATE_CONDA="/beegfs/scratch/ric.broccoli/kubacki.michal/conda/bin/conda"

echo "========================================"
echo "Coverage/Tile Visualization: Granule Cells (Separate Genotypes)"
echo "Started: $(date)"
echo "========================================"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline

Rscript scripts/separate_version/visualize_coverage_tile_gc_separate.R

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
