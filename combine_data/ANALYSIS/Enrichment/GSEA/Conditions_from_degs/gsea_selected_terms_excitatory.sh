#!/bin/bash
#SBATCH --job-name=gsea_selected_excitatory
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/gsea_selected_excitatory.out
#SBATCH --error=logs/gsea_selected_excitatory.err

# Selected Terms GSEA Heatmap for Excitatory L1 Neurons
# This script creates a heatmap with only user-specified terms from a selected terms list

echo "=========================================="
echo "GSEA Selected Terms Heatmap - Excitatory L1"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "=========================================="

# Create logs directory if it doesn't exist
mkdir -p logs

# Activate conda environment (adjust if needed)
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/bin/activate
conda activate scrna-analysis

# Print Python and package versions
echo "Python version:"
python --version
echo ""
echo "Key packages:"
python -c "import pandas; print(f'pandas: {pandas.__version__}')"
python -c "import matplotlib; print(f'matplotlib: {matplotlib.__version__}')"
python -c "import seaborn; print(f'seaborn: {seaborn.__version__}')"
echo ""

# Define paths
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/ANALYSIS/Enrichment/GSEA/Conditions_from_degs"
SELECTED_TERMS_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/ANALYSIS/Enrichment/GSEA/selected_terms.txt"

echo "Selected terms file: $SELECTED_TERMS_FILE"
echo ""

# Run the script
cd "$SCRIPT_DIR"
python gsea_selected_terms_excitatory.py \
    --terms_file "$SELECTED_TERMS_FILE"

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "=========================================="

exit $EXIT_CODE
