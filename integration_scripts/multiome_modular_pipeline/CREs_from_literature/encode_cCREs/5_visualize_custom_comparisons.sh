#!/bin/bash
#SBATCH --job-name=5_visualize_encode
#SBATCH --output=logs/5_visualize_custom_comparisons.log
#SBATCH --error=logs/5_visualize_custom_comparisons.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=workq

################################################################################
# Create Custom Comparison Visualizations for ENCODE cCREs
#
# This script creates publication-quality metaprofiles and comparison plots
# for ENCODE cCREs associated with splicing genes.
#
# Prerequisites:
# - output/heatmaps_deeptools/matrix_GABA.tab (or other matrices)
# - output/encode_cCREs_GABA.bed
# - output/encode_cCREs_GABA.tsv
#
# Output:
# - Metaprofiles comparing conditions (Ctrl vs Mut, genotypes)
# - Optional individual CRE plots
# - Statistical summaries
#
# Usage:
#   sbatch 5_visualize_custom_comparisons.sh              # Default (GABA, fast)
#   sbatch 5_visualize_custom_comparisons.sh --full       # With individual plots
#   sbatch 5_visualize_custom_comparisons.sh --all        # All cell types matrix
################################################################################

echo "========================================================================"
echo "CREATE CUSTOM COMPARISON VISUALIZATIONS FOR ENCODE cCREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/encode_cCREs"

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

# Check if Python is available
if ! command -v python &> /dev/null; then
    echo "ERROR: Python not found!"
    exit 1
fi

echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo ""

# Create logs directory if needed
mkdir -p logs

# Parse command line arguments
FULL_MODE=false
MATRIX="GABA"
PARALLEL=8

for arg in "$@"; do
    case $arg in
        --full)
            FULL_MODE=true
            shift
            ;;
        --all)
            MATRIX="all_celltypes"
            shift
            ;;
        --matrix=*)
            MATRIX="${arg#*=}"
            shift
            ;;
        *)
            ;;
    esac
done

echo "Configuration:"
echo "  Matrix: $MATRIX"
echo "  Full mode (with individual plots): $FULL_MODE"
echo "  Parallel processes: $PARALLEL"
echo ""

# Build command
if [ "$FULL_MODE" = true ]; then
    echo "Running in FULL MODE (with individual CRE plots)..."
    python 5_visualize_custom_comparisons.py \
        --matrix "$MATRIX" \
        --parallel $PARALLEL \
        --individual-dpi 150 \
        --max-individual 100
else
    echo "Running in FAST MODE (metaprofiles only)..."
    python 5_visualize_custom_comparisons.py \
        --matrix "$MATRIX" \
        --skip-individual
fi

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "VISUALIZATION COMPLETE!"
    echo "========================================================================"
    echo ""
    echo "Output files in: ./output/custom_comparisons/"
    echo ""
    ls -la ./output/custom_comparisons/ 2>/dev/null || echo "(directory not found)"
    echo ""
    if [ -d "./output/custom_comparisons/profiles" ]; then
        echo "Profile plots:"
        ls -la ./output/custom_comparisons/profiles/*.png 2>/dev/null | head -20
    fi
else
    echo ""
    echo "ERROR: Visualization script failed!"
    exit 1
fi

echo ""
echo "Completed: $(date)"
echo "========================================================================"
