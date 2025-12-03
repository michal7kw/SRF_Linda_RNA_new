#!/bin/bash
#SBATCH --job-name=Striatum_EP_viz
#SBATCH --output=logs/5_visualize_comparisons.log
#SBATCH --error=logs/5_visualize_comparisons.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --partition=workq

# Visualize Comparisons for Striatum_EP CREs

echo "========================================================================"
echo "VISUALIZE COMPARISONS FOR STRIATUM_EP CREs"
echo "========================================================================"
echo "Started: $(date)"

# Activate environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate striatum_viz

# Configuration
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/Striatum_EP_analysis"
cd $BASE_DIR

# Run script
python 5_visualize_comparisons.py

echo "========================================================================"
echo "Finished: $(date)"
