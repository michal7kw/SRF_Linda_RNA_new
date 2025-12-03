#!/bin/bash
#SBATCH --job-name=differential_CRE_analysis
#SBATCH --output=logs/4_differential_CRE_%j.out
#SBATCH --error=logs/4_differential_CRE_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=workq

# =============================================================================
# Step 4: Differential CRE Analysis
# =============================================================================
#
# This script identifies differential CREs between conditions by comparing
# ATAC signal levels. CREs are analyzed separately for enhancers and silencers.
#
# USAGE:
#   sbatch 4_differential_CRE_analysis.sh
#
# ENVIRONMENT VARIABLES (optional):
#   MIN_SIGNAL=2.0     Minimum signal threshold
#   MIN_FC=3.0         Minimum fold change threshold
#   SKIP_INDIVIDUAL=0  Set to 1 to skip individual plots
#   USE_GABA=0         Set to 1 to use GABA-specific CREs only
#
# =============================================================================

set -e

# Script directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_directional"
cd "${SCRIPT_DIR}"

# Create logs directory
mkdir -p logs

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate rna_seq_analysis_deep

# Parse environment variables with defaults
# NOTE: Default thresholds adjusted for normalized BigWig signal values (0.01-1.0 range)
MIN_SIGNAL="${MIN_SIGNAL:-0.05}"
MIN_FC="${MIN_FC:-1.5}"
SKIP_INDIVIDUAL="${SKIP_INDIVIDUAL:-0}"
MAX_INDIVIDUAL="${MAX_INDIVIDUAL:-20}"
USE_GABA="${USE_GABA:-0}"
DIAGNOSTICS_ONLY="${DIAGNOSTICS_ONLY:-0}"

echo "=============================================================================="
echo "STEP 4: DIFFERENTIAL CRE ANALYSIS"
echo "=============================================================================="
echo ""
echo "Date: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Working directory: ${SCRIPT_DIR}"
echo ""
echo "Parameters:"
echo "  MIN_SIGNAL: ${MIN_SIGNAL}"
echo "  MIN_FC: ${MIN_FC}"
echo "  SKIP_INDIVIDUAL: ${SKIP_INDIVIDUAL}"
echo "  MAX_INDIVIDUAL: ${MAX_INDIVIDUAL}"
echo "  USE_GABA: ${USE_GABA}"
echo ""

# Build arguments
ARGS=""
ARGS="${ARGS} --min-signal ${MIN_SIGNAL}"
ARGS="${ARGS} --min-fc ${MIN_FC}"
ARGS="${ARGS} --max-individual ${MAX_INDIVIDUAL}"

if [[ "${SKIP_INDIVIDUAL}" == "1" ]]; then
    ARGS="${ARGS} --skip-individual"
fi

if [[ "${USE_GABA}" == "1" ]]; then
    ARGS="${ARGS} --use-gaba"
fi

if [[ "${DIAGNOSTICS_ONLY}" == "1" ]]; then
    ARGS="${ARGS} --diagnostics-only"
fi

echo "Running: python 4_differential_CRE_analysis.py ${ARGS}"
echo ""

python "${SCRIPT_DIR}/4_differential_CRE_analysis.py" ${ARGS}

echo ""
echo "=============================================================================="
echo "STEP 4 COMPLETE"
echo "=============================================================================="
echo "End time: $(date)"
