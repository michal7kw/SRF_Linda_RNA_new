#!/bin/bash
#SBATCH --job-name=gsea_focused_heatmaps
#SBATCH --output=logs/gsea_focused_heatmaps.out
#SBATCH --error=logs/gsea_focused_heatmaps.err
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

# Publication-Ready GSEA Heatmap Generator
# Creates compact, publication-quality heatmaps from GSEA results
#
# Generate focused GSEA heatmaps for specific biological categories
#
# This script creates small, focused heatmaps summarizing GSEA results:
# 1. Cell death categories in GC cluster for both Emx1 and Nestin models
# 2. Splicing categories in IN cluster for both Emx1 and Nestin models
#
# NOTE: The script automatically generates TWO versions of each heatmap:
#   - Standard version WITH FDR significance stars (*, **, ***)
#   - Clean version WITHOUT FDR indications (filename suffix: _no_fdr)
#
# Usage:
#   ./gsea_focused_heatmaps.sh                    # Generate both heatmaps
#   ./gsea_focused_heatmaps.sh cell_death         # Generate only cell death heatmap
#   ./gsea_focused_heatmaps.sh splicing           # Generate only splicing heatmap
#

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

# Set working directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/ANALYSIS/Enrichment/GSEA/Conditions_from_degs"
cd "$SCRIPT_DIR"

# Set default focus area (both, cell_death, or splicing)
FOCUS="${1:-both}"

# Set max terms to display (default: 12, optimized for compact figures)
MAX_TERMS="${2:-12}"

# Set FDR threshold (default: 0.25)
FDR_THRESHOLD="${3:-0.25}"

echo "========================================================================"
echo "Focused GSEA Heatmap Generator"
echo "========================================================================"
echo "Focus area: $FOCUS"
echo "Max terms: $MAX_TERMS"
echo "FDR threshold: $FDR_THRESHOLD"
echo "========================================================================"

# Run the Python script
python3 "$SCRIPT_DIR/gsea_focused_heatmaps.py" \
    --focus "$FOCUS" \
    --max_terms "$MAX_TERMS" \
    --fdr_threshold "$FDR_THRESHOLD"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "✓ SUCCESS: Heatmaps generated successfully!"
    echo "✓ Generated BOTH versions for each category:"
    echo "    - Standard version with FDR significance stars"
    echo "    - Clean version without FDR (_no_fdr suffix)"
    echo "✓ Check output: results_from_raw/gsea_focused_heatmaps/"
    echo "========================================================================"
else
    echo ""
    echo "========================================================================"
    echo "✗ ERROR: Heatmap generation failed with exit code $EXIT_CODE"
    echo "========================================================================"
fi

exit $EXIT_CODE
