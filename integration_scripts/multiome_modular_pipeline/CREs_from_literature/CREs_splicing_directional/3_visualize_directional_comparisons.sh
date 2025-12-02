#!/bin/bash
#SBATCH --job-name=viz_splicing
#SBATCH --output=logs/3_visualize_directional_%j.log
#SBATCH --error=logs/3_visualize_directional_%j.log
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=workq

# =============================================================================
# Step 3: Visualize Directional Comparisons
# =============================================================================

set -e

# Absolute path to script directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_directional"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate rna_seq_analysis_deep

# Change to script directory
cd "${SCRIPT_DIR}"

echo "=============================================================================="
echo "VISUALIZING DIRECTIONAL COMPARISONS"
echo "=============================================================================="
echo "Date: $(date)"
echo ""

# Parse arguments from environment variables
SKIP_INDIVIDUAL="${SKIP_INDIVIDUAL:-0}"
MAX_INDIVIDUAL="${MAX_INDIVIDUAL:-50}"
MIN_SIGNAL="${MIN_SIGNAL:-1.0}"
MIN_FC="${MIN_FC:-2.0}"
INDIVIDUAL_DPI="${INDIVIDUAL_DPI:-150}"

echo "Settings:"
echo "  SKIP_INDIVIDUAL: ${SKIP_INDIVIDUAL}"
echo "  MAX_INDIVIDUAL: ${MAX_INDIVIDUAL}"
echo "  MIN_SIGNAL: ${MIN_SIGNAL}"
echo "  MIN_FC: ${MIN_FC}"
echo "  INDIVIDUAL_DPI: ${INDIVIDUAL_DPI}"
echo ""

# Build arguments
ARGS=""
if [[ "${SKIP_INDIVIDUAL}" == "1" ]]; then
    ARGS="${ARGS} --skip-individual"
fi
ARGS="${ARGS} --max-individual ${MAX_INDIVIDUAL}"
ARGS="${ARGS} --min-signal ${MIN_SIGNAL}"
ARGS="${ARGS} --min-fc ${MIN_FC}"
ARGS="${ARGS} --individual-dpi ${INDIVIDUAL_DPI}"

# Run visualization script
echo "Working directory: $(pwd)"
echo "Running: python 3_visualize_directional_comparisons.py ${ARGS}"
echo ""
python "${SCRIPT_DIR}/3_visualize_directional_comparisons.py" ${ARGS}

echo ""
echo "Step 3 complete. Date: $(date)"

