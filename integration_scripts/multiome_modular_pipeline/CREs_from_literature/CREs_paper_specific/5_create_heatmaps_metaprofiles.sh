#!/bin/bash
#SBATCH --job-name=5_create_heatmaps_metaprofiles
#SBATCH --output=logs/5_create_heatmaps_metaprofiles.log
#SBATCH --error=logs/5_create_heatmaps_metaprofiles.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=workq

# Activate conda environment with pyBigWig
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

echo ""
echo "========================================================================"
echo "Running Python analysis..."
echo "========================================================================"
echo ""

# Run the Python script
python 5_create_heatmaps_metaprofiles.py
