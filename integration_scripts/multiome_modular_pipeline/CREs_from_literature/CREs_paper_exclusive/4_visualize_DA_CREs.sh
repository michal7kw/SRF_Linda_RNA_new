#!/bin/bash
#SBATCH --job-name=4_visualize_DA_CREs
#SBATCH --output=logs/4_visualize_DA_CREs.log
#SBATCH --error=logs/4_visualize_DA_CREs.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=workq

################################################################################
# Visualize Differentially Accessible CREs with minSig/minFC Filtering
#
# This script creates individual CRE plots for cell-type-specific CREs that
# show significant differential accessibility between conditions.
#
# Prerequisites:
# - output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/matrix_*.gz
# - output/GABA_specific_CREs_genes.tsv (from 3_link_CREs_to_genes.py)
#
# Performance options (set before sbatch):
# - SKIP_INDIVIDUAL=1: Skip individual plots (default: 0)
# - PARALLEL_JOBS=8: Number of parallel processes (default: 8)
# - INDIVIDUAL_DPI=150: DPI for individual plots (default: 150)
# - MIN_SIGNAL=2.0: Minimum signal threshold (default: 2.0)
# - MIN_FC=3.0: Minimum fold change threshold (default: 3.0)
#
# Usage:
#   sbatch 4_visualize_DA_CREs.sh                    # Full analysis
#   SKIP_INDIVIDUAL=1 sbatch 4_visualize_DA_CREs.sh  # Fast mode (metaprofiles only)
#   MIN_SIGNAL=2.0 MIN_FC=3.0 sbatch 4_visualize_DA_CREs.sh  # Stricter filtering
################################################################################

echo "========================================================================"
echo "VISUALIZE DIFFERENTIALLY ACCESSIBLE CREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Create logs directory
mkdir -p logs

# Set defaults
SKIP_INDIVIDUAL=${SKIP_INDIVIDUAL:-0}
PARALLEL_JOBS=${PARALLEL_JOBS:-8}
INDIVIDUAL_DPI=${INDIVIDUAL_DPI:-150}
MIN_SIGNAL=${MIN_SIGNAL:-2.0}
MIN_FC=${MIN_FC:-3.0}

echo "Configuration:"
echo "  SKIP_INDIVIDUAL: $SKIP_INDIVIDUAL"
echo "  PARALLEL_JOBS: $PARALLEL_JOBS"
echo "  INDIVIDUAL_DPI: $INDIVIDUAL_DPI"
echo "  MIN_SIGNAL: $MIN_SIGNAL"
echo "  MIN_FC: $MIN_FC"
echo ""

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

# Build command
CMD="python 4_visualize_DA_CREs.py"
CMD="$CMD --min-signal $MIN_SIGNAL --min-fc $MIN_FC"

if [ "$SKIP_INDIVIDUAL" == "1" ]; then
    CMD="$CMD --skip-individual"
else
    CMD="$CMD --parallel $PARALLEL_JOBS --individual-dpi $INDIVIDUAL_DPI"
fi

echo "Running: $CMD"
echo ""

# Run the Python script
$CMD

echo ""
echo "========================================================================"
echo "Completed: $(date)"
echo "========================================================================"
