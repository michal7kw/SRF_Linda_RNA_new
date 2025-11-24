#!/bin/bash
#SBATCH --job-name=2_fold_enrichment_exclusive_analysis
#SBATCH --output=logs/2_fold_enrichment_exclusive_analysis.log
#SBATCH --error=logs/2_fold_enrichment_exclusive_analysis.err
#SBATCH --time=10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=workq

################################################################################
# Fold-Enrichment Analysis from deepTools Matrices
#
# This script runs POST-PROCESSING analysis on existing deepTools matrix files
# to calculate fold-enrichment statistics and create publication-ready figures.
#
# Key Features:
# - Reads existing deepTools matrices (NO recomputation)
# - Calculates fold-enrichment: GABA-specific / Excitatory-specific
# - Generates quantitative statistics with biological interpretation
# - Creates publication-ready comparison figures
# - FAST: ~1 minute (vs 3-4 hours for full Python pipeline)
#
# Prerequisites:
# - deepTools matrices must exist from previous run:
#   output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/matrix_*.gz
#
# Usage:
#   sbatch run_fold_enrichment_analysis.sh
#
# Output Files:
# - fold_enrichment_statistics.txt: Quantitative statistics
# - comparison_fold_enrichment.png: 4-panel comparison figure
# - fold_enrichment_by_condition.png: Detailed fold-enrichment plot
################################################################################

echo "========================================================================"
echo "FOLD-ENRICHMENT ANALYSIS FROM DEEPTOOLS MATRICES"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Navigate to correct directory
# cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature

# ============================================================================
# Check Prerequisites
# ============================================================================

echo "========================================================================"
echo "Checking prerequisites..."
echo "========================================================================"
echo ""

OUTPUT_DIR="output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools"
MISSING=0

if [ ! -f "$MATRIX_DIR/matrix_GABA_specific.gz" ]; then
    echo "✗ Missing: matrix_GABA_specific.gz"
    MISSING=$((MISSING + 1))
else
    echo "✓ Found: matrix_GABA_specific.gz"
fi

if [ ! -f "$MATRIX_DIR/matrix_Excitatory_specific.gz" ]; then
    echo "✗ Missing: matrix_Excitatory_specific.gz"
    MISSING=$((MISSING + 1))
else
    echo "✓ Found: matrix_Excitatory_specific.gz"
fi

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "ERROR: Matrix files not found!"
    echo ""
    echo "Please run deepTools pipeline first:"
    echo "  sbatch 1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh"
    echo ""
    echo "Or use fast replot if matrices exist elsewhere:"
    echo "  sbatch REPLOT_ONLY.sh"
    exit 1
fi

echo ""
echo "✓ All required matrix files found"
echo ""

# ============================================================================
# Activate Conda Environment
# ============================================================================

echo "========================================================================"
echo "Activating conda environment..."
echo "========================================================================"
echo ""

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# Check Python packages
echo "Checking Python packages..."
python -c "import numpy, pandas, matplotlib, seaborn, scipy; print('  ✓ All packages available')"

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Required Python packages not available!"
    echo "Please install: numpy, pandas, matplotlib, seaborn, scipy"
    exit 1
fi

echo ""

# ============================================================================
# Run Fold-Enrichment Analysis
# ============================================================================

echo "========================================================================"
echo "Running fold-enrichment analysis..."
echo "========================================================================"
echo ""

python 2_fold_enrichment_exclusive_analysis.py

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Fold-enrichment analysis failed!"
    echo "Check error messages above for details"
    exit 1
fi

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "========================================================================"
echo "ANALYSIS COMPLETE!"
echo "========================================================================"
echo ""
echo "Output directory:"
echo "  $MATRIX_DIR/"
echo ""
echo "Generated files:"
echo "  ★ Statistics:"
echo "      - fold_enrichment_statistics.txt"
echo "  ★ Figures:"
echo "      - comparison_fold_enrichment.png (4-panel comparison)"
echo "      - fold_enrichment_by_condition.png (detailed bar chart)"
echo ""
echo "Key advantages of this approach:"
echo "  ✓ Fast: ~1 minute (vs 3-4 hours for full Python pipeline)"
echo "  ✓ Reuses deepTools matrices (no recomputation)"
echo "  ✓ Publication-ready statistics and figures"
echo "  ✓ Quantitative fold-enrichment values"
echo "  ✓ Biological interpretation provided"
echo ""
echo "Next steps:"
echo "  1. Review fold_enrichment_statistics.txt for quantitative results"
echo "  2. Check comparison_fold_enrichment.png for visual summary"
echo "  3. Use fold-enrichment values in publication"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
