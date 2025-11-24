#!/bin/bash

################################################################################
# SLURM Submission Script for Signac Multiome Analysis Pipeline (L1 VERSION)
################################################################################

#SBATCH --job-name=run_step4
#SBATCH --output=logs/run_step4.out
#SBATCH --error=logs/run_step4.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=workq
#SBATCH --account=kubacki.michal

# It runs:
# - signac_04_visualization.R

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

# CRITICAL: Tell R's reticulate to use conda Python (not system Python)
export RETICULATE_PYTHON=$(which python)
export RETICULATE_CONDA=$(which conda)
echo "Python for R: ${RETICULATE_PYTHON}"
echo "Conda for R: ${RETICULATE_CONDA}"


# Set number of threads for R
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

################################################################################
# Step 4: Advanced Visualization
################################################################################

echo "================================================================================"
echo "STEP 4: Advanced Visualization (L2)"
echo "================================================================================"
echo ""

START_TIME=$(date +%s)

Rscript ${BASE_DIR}/signac_04_visualization.R 2>&1

EXIT_CODE=$?

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "✗ Step 4 failed with exit code ${EXIT_CODE}"
    echo "  Check logs for errors"
    exit 1
fi

echo ""
echo "✓ Step 4 complete (elapsed time: ${ELAPSED}s)"
echo ""