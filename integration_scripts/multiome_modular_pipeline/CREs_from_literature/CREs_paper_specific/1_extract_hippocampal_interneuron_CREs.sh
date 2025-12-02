#!/bin/bash
#SBATCH --job-name=1_extract_hippocampal_interneuron_CREs
#SBATCH --output=logs/1_extract_hippocampal_interneuron_CREs.log
#SBATCH --error=logs/1_extract_hippocampal_interneuron_CREs.err
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
python 1_extract_hippocampal_interneuron_CREs.py
