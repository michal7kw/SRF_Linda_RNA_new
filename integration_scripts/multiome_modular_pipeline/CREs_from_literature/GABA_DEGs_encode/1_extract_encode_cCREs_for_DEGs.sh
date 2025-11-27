#!/bin/bash
#SBATCH --job-name=1_extract_DEGs_encode
#SBATCH --output=logs/1_extract_encode_cCREs_for_DEGs.log
#SBATCH --error=logs/1_extract_encode_cCREs_for_DEGs.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=workq

echo "Starting Step 1: Extract ENCODE cCREs for GABA DEGs"
echo "Date: $(date)"

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/GABA_DEGs_encode"

# Activate conda environment with bedtools (sc-chromatin2 has both bedtools and pandas)
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate sc-chromatin2

# Verify bedtools is available
echo "Checking bedtools: $(which bedtools)"

# Run extraction script
python 1_extract_encode_cCREs_for_DEGs.py

echo "Step 1 completed: $(date)"
