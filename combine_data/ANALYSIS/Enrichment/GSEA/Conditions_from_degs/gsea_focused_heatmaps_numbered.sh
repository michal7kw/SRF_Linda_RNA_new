#!/bin/bash
#SBATCH --job-name=gsea_numbered_heatmaps
#SBATCH --output=logs/gsea_numbered_heatmaps.out
#SBATCH --error=logs/gsea_numbered_heatmaps.err
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

# Numbered GSEA Heatmap Generator
# Creates ultra-compact heatmaps with numbered terms and separate legend table
#
# This design is ideal for:
# - Very limited figure space
# - Multi-panel figures with many subplots
# - Professional publications with space constraints
#
# NOTE: The script automatically generates TWO versions of each heatmap:
#   - Standard version WITH FDR significance stars (*, **, ***)
#   - Clean version WITHOUT FDR indications (filename suffix: _no_fdr)
#
# Usage:
#   ./gsea_focused_heatmaps_numbered.sh                    # Generate both heatmaps
#   ./gsea_focused_heatmaps_numbered.sh cell_death         # Generate only cell death
#   ./gsea_focused_heatmaps_numbered.sh splicing           # Generate only splicing
#   ./gsea_focused_heatmaps_numbered.sh both 10            # 10 terms instead of 12
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

# Set max terms to display (default: 12)
MAX_TERMS="${2:-12}"

# Set FDR threshold (default: 0.25)
FDR_THRESHOLD="${3:-0.25}"

echo "========================================================================"
echo "Numbered GSEA Heatmap Generator (Ultra-Compact)"
echo "========================================================================"
echo "Focus area: $FOCUS"
echo "Max terms: $MAX_TERMS"
echo "FDR threshold: $FDR_THRESHOLD"
echo "Format: Numbers in heatmap, names in side-by-side legend table"
echo "========================================================================"

# Run the Python script
python3 "$SCRIPT_DIR/gsea_focused_heatmaps_numbered.py" \
    --focus "$FOCUS" \
    --max_terms "$MAX_TERMS" \
    --fdr_threshold "$FDR_THRESHOLD"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "✓ SUCCESS: Numbered heatmaps generated successfully!"
    echo "✓ Generated BOTH versions for each category:"
    echo "    - Standard version with FDR significance stars"
    echo "    - Clean version without FDR (_no_fdr suffix)"
    echo "✓ Check output directory: results_from_raw/gsea_focused_heatmaps_numbered/"
    echo "✓ Files include: PNG, PDF, and CSV legend"
    echo "========================================================================"
else
    echo ""
    echo "========================================================================"
    echo "✗ ERROR: Heatmap generation failed with exit code $EXIT_CODE"
    echo "========================================================================"
fi

exit $EXIT_CODE
