#!/bin/bash
#SBATCH --job-name=gsea_selected_terms_heatmap
#SBATCH --output=logs/gsea_selected_terms_heatmap.out
#SBATCH --error=logs/gsea_selected_terms_heatmap.err
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

# Selected Terms GSEA Heatmap Generator
# Creates a heatmap with only user-specified terms from selected_terms.txt
# No FDR significance indicators are shown on the heatmap
#
# Usage:
#   ./gsea_selected_terms_heatmap.sh                              # Use default Combined_GC
#   ./gsea_selected_terms_heatmap.sh Mature_GC                    # Create Mature_GC plot
#   ./gsea_selected_terms_heatmap.sh Mature_GC GABA               # Create Mature_GC with GABA order/scale
#   ./gsea_selected_terms_heatmap.sh Combined_GC GABA             # Create Combined_GC with GABA order/scale
#

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment pydata-env..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/pydata-env

# Set working directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/ANALYSIS/Enrichment/GSEA"
cd "$SCRIPT_DIR"

# Parse arguments
CLUSTER="${1:-Combined_GC}"           # Default cluster to analyze
REFERENCE_CLUSTER="${2:-}"            # Optional reference cluster for order/scale
TERMS_FILE="${3:-$SCRIPT_DIR/selected_terms.txt}"

echo "========================================================================"
echo "Selected Terms GSEA Heatmap Generator"
echo "========================================================================"
echo "Cluster: $CLUSTER"
if [ -n "$REFERENCE_CLUSTER" ]; then
    echo "Reference cluster: $REFERENCE_CLUSTER (for term order and color scale)"
fi
echo "Terms file: $TERMS_FILE"
echo "Format: Numbers in heatmap, names in side-by-side legend table"
echo "Note: NO FDR significance indicators on heatmap"
echo "========================================================================"

# Build Python command
PYTHON_CMD="python3 $SCRIPT_DIR/gsea_selected_terms_heatmap.py --cluster $CLUSTER --terms_file $TERMS_FILE"

# Add reference cluster if provided
if [ -n "$REFERENCE_CLUSTER" ]; then
    PYTHON_CMD="$PYTHON_CMD --reference_cluster $REFERENCE_CLUSTER"
fi

# Run the Python script
eval $PYTHON_CMD


EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "✓ SUCCESS: Selected terms heatmap generated successfully!"
    echo "✓ Check output directory: results_from_raw/gsea_selected_terms/"
    echo "✓ Files include: PNG, PDF, and CSV legend"
    echo "========================================================================"
else
    echo ""
    echo "========================================================================"
    echo "✗ ERROR: Heatmap generation failed with exit code $EXIT_CODE"
    echo "========================================================================"
fi

exit $EXIT_CODE
