#!/bin/bash
#SBATCH --job-name=5_visualize_bigwig
#SBATCH --output=logs/5_visualize_bigwig.log
#SBATCH --error=logs/5_visualize_bigwig.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=workq

################################################################################
# Visualize BigWig Signal at Splicing Gene ENCODE cCREs
#
# This script directly reads BigWig files and creates publication-quality
# visualizations of ATAC-seq signal at ENCODE cCREs.
#
# PERFORMANCE OPTIONS:
# -----------------------------------------
# DEFAULT: Quick mode (metaprofiles only):
#   sbatch 5_visualize_bigwig_signal.sh
#
# Full mode (create individual plots):
#   SKIP_INDIVIDUAL=0 sbatch 5_visualize_bigwig_signal.sh
#
# Full mode with parallel processing:
#   SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 5_visualize_bigwig_signal.sh
#
# Prerequisites:
# - output/splicing_encode_cCREs_all.bed
# - BigWig files from Signac pipeline
#
# Output:
# - Publication-quality metaprofiles
# - Ctrl vs Mut comparison plots
# - Individual CRE profiles (optional)
# - Statistical summaries
################################################################################

echo "========================================================================"
echo "VISUALIZE BIGWIG SIGNAL AT SPLICING GENE ENCODE cCREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/splicing_encode_cCREs"

# Activate environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# ============================================================================
# Performance options
# ============================================================================

SKIP_FLAG="--skip-individual"
if [ "${SKIP_INDIVIDUAL}" = "0" ]; then
    echo "Full mode: Creating individual plots for each CRE"
    SKIP_FLAG=""
else
    echo "Fast mode (DEFAULT): Skipping individual CRE plots"
    echo "   (To create individual plots, set SKIP_INDIVIDUAL=0)"
fi

PARALLEL_FLAG=""
if [ -n "${PARALLEL_JOBS}" ]; then
    PARALLEL_FLAG="--parallel ${PARALLEL_JOBS}"
    echo "   Using ${PARALLEL_JOBS} parallel processes"
fi

DPI_FLAG=""
if [ -n "${INDIVIDUAL_DPI}" ]; then
    DPI_FLAG="--individual-dpi ${INDIVIDUAL_DPI}"
    echo "   Individual plot DPI: ${INDIVIDUAL_DPI}"
fi

FILTER_FLAG=""
if [ -n "${MIN_SIGNAL}" ]; then
    FILTER_FLAG="$FILTER_FLAG --min-signal ${MIN_SIGNAL}"
    echo "   Min Signal: ${MIN_SIGNAL}"
fi
if [ -n "${MIN_FC}" ]; then
    FILTER_FLAG="$FILTER_FLAG --min-fc ${MIN_FC}"
    echo "   Min FC: ${MIN_FC}"
fi

echo ""

# ============================================================================
# Run visualization
# ============================================================================

echo "Running BigWig visualization script..."
echo ""

python 5_visualize_bigwig_signal.py $SKIP_FLAG $PARALLEL_FLAG $DPI_FLAG $FILTER_FLAG

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Visualization failed!"
    exit 1
fi

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "========================================================================"
echo "VISUALIZATION COMPLETE!"
echo "========================================================================"
echo ""
echo "Output directory: output/bigwig_profiles_*/"
echo ""
echo "Generated files:"
echo "  - metaprofile_nestin_ctrl_vs_mut.png"
echo "  - metaprofile_emx1_ctrl_vs_mut.png"
echo "  - summary_nestin.tsv"
echo "  - summary_emx1.tsv"
if [ "${SKIP_INDIVIDUAL}" = "0" ]; then
    echo "  - individual_Nestin_*.png"
    echo "  - individual_Emx1_*.png"
fi
echo ""
echo "Completed: $(date)"
echo "========================================================================"
