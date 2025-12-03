#!/bin/bash
# =============================================================================
# SCRIPT: 8_differential_motif_analysis.sh
# PURPOSE: SLURM wrapper for differential motif enrichment analysis
# =============================================================================
#SBATCH --job-name=diff_motif
#SBATCH --output=logs/8_differential_motif_%j.out
#SBATCH --error=logs/8_differential_motif_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=workq

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

# Activate conda environment with pandas, scipy, numpy
CONDA_BASE="/beegfs/scratch/ric.broccoli/kubacki.michal/conda"
CONDA_ENV="/beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/analysis3_env"

echo "[$(date)] Starting differential motif analysis..."
source "${CONDA_BASE}/bin/activate" "${CONDA_ENV}"

# Run the Python script
python "${SCRIPT_DIR}/8_differential_motif_analysis.py"

echo "[$(date)] Differential motif analysis complete."
