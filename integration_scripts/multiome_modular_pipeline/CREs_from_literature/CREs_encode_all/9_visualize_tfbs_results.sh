#!/bin/bash
# =============================================================================
# SCRIPT: 9_visualize_tfbs_results.sh
# PURPOSE: SLURM wrapper for TFBS visualization
# =============================================================================
#SBATCH --job-name=vis_tfbs
#SBATCH --output=logs/9_visualize_tfbs_%j.out
#SBATCH --error=logs/9_visualize_tfbs_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=workq

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

# Activate conda environment with matplotlib, seaborn, pandas
CONDA_BASE="/beegfs/scratch/ric.broccoli/kubacki.michal/conda"
CONDA_ENV="/beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/analysis3_env"

echo "[$(date)] Starting TFBS visualization..."
source "${CONDA_BASE}/bin/activate" "${CONDA_ENV}"

# Run the Python script
python "${SCRIPT_DIR}/9_visualize_tfbs_results.py"

echo "[$(date)] TFBS visualization complete."
