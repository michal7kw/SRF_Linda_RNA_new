#!/bin/bash
#SBATCH --job-name=5_visualize_bigwig_signal
#SBATCH --output=logs/5_visualize_bigwig_signal.log
#SBATCH --error=logs/5_visualize_bigwig_signal.err
#SBATCH --time=30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=workq

################################################################################
# Visualize BigWig Signal at Splicing Gene CREs
#
# This script directly reads BigWig files and creates publication-quality
# visualizations showing ATAC-seq signal at CREs linked to splicing genes.
#
# PERFORMANCE OPTIONS:
# -----------------------------------------
# DEFAULT: Quick mode (metaprofiles only, skips individual plots):
#   sbatch 5_visualize_bigwig_signal.sh
#
# Full mode (create individual plots):
#   SKIP_INDIVIDUAL=0 sbatch 5_visualize_bigwig_signal.sh
#
# Full mode with parallel processing (8x faster):
#   SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 5_visualize_bigwig_signal.sh
#
# Full mode with lower DPI (2x faster):
#   SKIP_INDIVIDUAL=0 INDIVIDUAL_DPI=100 sbatch 5_visualize_bigwig_signal.sh
#
# Combined (parallel + low DPI):
#   SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 INDIVIDUAL_DPI=100 sbatch 5_visualize_bigwig_signal.sh
#
# Output:
# - Metaprofiles comparing Ctrl vs Mut (always 300 DPI)
# - Individual CRE plots with difference panels (optional, configurable DPI)
# - Summary statistics tables
################################################################################

echo "========================================================================"
echo "VISUALIZE BIGWIG SIGNAL AT SPLICING GENE CREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/splicing_genes"

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# Performance options (can be set via environment variables):
# - SKIP_INDIVIDUAL=0: Create individual CRE plots (default: skip to save 30-50 minutes)
# - PARALLEL_JOBS=N: Use N parallel processes for individual plots
# - INDIVIDUAL_DPI=N: DPI for individual plots (default: 150)

SKIP_FLAG="--skip-individual"
if [ "${SKIP_INDIVIDUAL}" = "0" ]; then
    echo "📊 Full mode: Creating individual plots for both genotypes"
    SKIP_FLAG=""
else
    echo "⚡ Fast mode (DEFAULT): Skipping individual CRE plots"
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

echo ""

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
# Run visualization
python 5_visualize_bigwig_signal.py $SKIP_FLAG $PARALLEL_FLAG $DPI_FLAG $FILTER_FLAG

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Visualization failed!"
    exit 1
fi

echo ""
echo "========================================================================"
echo "COMPLETE!"
echo "========================================================================"
echo ""
echo "Output directory: output/splicing_genes_analysis/bigwig_profiles/"
echo ""
echo "Generated files:"
echo ""
echo "  📊 Metaprofiles (always created, 300 DPI):"
echo "      ├─ metaprofile_nestin_ctrl_vs_mut.png"
echo "      └─ metaprofile_emx1_ctrl_vs_mut.png"
echo ""

if [ "${SKIP_INDIVIDUAL}" = "0" ]; then
    echo "  📈 Individual CRE plots (DPI: ${INDIVIDUAL_DPI:-150}):"
    echo "      ├─ individual_Nestin_<Gene>_<CRE_ID>.png"
    echo "      └─ individual_Emx1_<Gene>_<CRE_ID>.png"
    if [ -n "${PARALLEL_JOBS}" ]; then
        echo "      (Created using ${PARALLEL_JOBS} parallel processes)"
    fi
    echo ""
else
    echo "  📈 Individual CRE plots: SKIPPED (default fast mode)"
    echo "      To create, run with SKIP_INDIVIDUAL=0"
    echo ""
fi

echo "  📋 Summary statistics:"
echo "      ├─ summary_nestin.tsv"
echo "      └─ summary_emx1.tsv"
echo ""
echo "Completed: $(date)"
echo "========================================================================"