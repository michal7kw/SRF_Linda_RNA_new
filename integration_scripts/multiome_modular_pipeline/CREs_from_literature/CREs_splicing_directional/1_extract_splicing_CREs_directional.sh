#!/bin/bash
#SBATCH --job-name=extract_splicing_CREs
#SBATCH --output=logs/1_extract_splicing_CREs_%j.log
#SBATCH --error=logs/1_extract_splicing_CREs_%j.log
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=workq

# =============================================================================
# Step 1: Extract Splicing Gene CREs with PCC Directionality
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
echo "EXTRACTING SPLICING GENE CREs WITH DIRECTIONALITY"
echo "=============================================================================="
echo "Date: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Run extraction script
python "${SCRIPT_DIR}/1_extract_splicing_CREs_directional.py"

echo ""
echo "Step 1 complete. Date: $(date)"
