#!/bin/bash
#SBATCH --job-name=extract_Striatum_EP
#SBATCH --output=logs/1_extract_Striatum_EP.log
#SBATCH --error=logs/1_extract_Striatum_EP.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --partition=workq

# Extract Striatum_EP CREs

echo "========================================================================"
echo "EXTRACT STRIATUM_EP CREs"
echo "========================================================================"
echo "Started: $(date)"

# Activate environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate striatum_viz

# Configuration
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/Striatum_EP_analysis"
cd $BASE_DIR

# Run script
python 1_extract_Striatum_EP.py

echo "========================================================================"
echo "Finished: $(date)"
