#!/bin/bash
#SBATCH --job-name=7_visualize_DA_CREs
#SBATCH --output=logs/7_visualize_DA_CREs.log
#SBATCH --error=logs/7_visualize_DA_CREs.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=workq

################################################################################
# Visualize Differentially Accessible CREs with minSig/minFC Filtering
#
# Creates individual CRE plots for CREs that show significant differential
# accessibility between conditions.
#
# Prerequisites:
# - output/GABA_DEG_analysis/heatmaps_deeptools/matrix_*.gz
# - output/hippocampal_interneuron_CREs_genes.tsv
#
# Performance options (set before sbatch):
# - SKIP_INDIVIDUAL=1: Skip individual plots (default: 0)
# - PARALLEL_JOBS=8: Number of parallel processes (default: 8)
# - INDIVIDUAL_DPI=150: DPI for individual plots (default: 150)
# - min_signal=2.0: Minimum signal threshold (default: 1.0)
# - min_fc=2.0: Minimum fold change threshold (default: 1.5)
#
# Usage:
#   sbatch 7_visualize_DA_CREs.sh                    # Full analysis
#   SKIP_INDIVIDUAL=1 sbatch 7_visualize_DA_CREs.sh  # Fast mode
################################################################################

echo "========================================================================"
echo "VISUALIZE DIFFERENTIALLY ACCESSIBLE CREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

mkdir -p logs

SKIP_INDIVIDUAL=${SKIP_INDIVIDUAL:-0}
PARALLEL_JOBS=${PARALLEL_JOBS:-8}
INDIVIDUAL_DPI=${INDIVIDUAL_DPI:-150}
MIN_SIGNAL=${MIN_SIGNAL:-2.0}
MIN_FC=${MIN_FC:-2.0}

echo "Configuration:"
echo "  SKIP_INDIVIDUAL: $SKIP_INDIVIDUAL"
echo "  PARALLEL_JOBS: $PARALLEL_JOBS"
echo "  MIN_SIGNAL: $MIN_SIGNAL"
echo "  MIN_FC: $MIN_FC"
echo ""

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

CMD="python 7_visualize_DA_CREs.py"
CMD="$CMD --min-signal $MIN_SIGNAL --min-fc $MIN_FC"

if [ "$SKIP_INDIVIDUAL" == "1" ]; then
    CMD="$CMD --skip-individual"
else
    CMD="$CMD --parallel $PARALLEL_JOBS --individual-dpi $INDIVIDUAL_DPI"
fi

echo "Running: $CMD"
$CMD

echo ""
echo "========================================================================"
echo "Completed: $(date)"
echo "========================================================================"
