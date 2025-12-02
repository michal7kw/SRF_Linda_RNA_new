#!/bin/bash
#SBATCH --job-name=1_extract_encode_cCREs
#SBATCH --output=logs/1_extract_encode_cCREs.log
#SBATCH --error=logs/1_extract_encode_cCREs.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --partition=workq

# Extract ENCODE cCREs Associated with Splicing Genes (OPTIMIZED VERSION)
# Uses bedtools intersect for fast overlap detection
# Chunked processing of large Table 16 file

echo "========================================================================"
echo "EXTRACT ENCODE cCREs ASSOCIATED WITH SPLICING GENES (OPTIMIZED)"
echo "========================================================================"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Activate conda environment with bedtools
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate sc-chromatin2

# Check if bedtools is available
if ! command -v bedtools &> /dev/null; then
    echo "bedtools not found in sc-chromatin2, trying rna_seq_analysis_deep..."
    conda activate rna_seq_analysis_deep
fi

# Verify bedtools is now available
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found in any environment!"
    exit 1
fi

echo "Using bedtools: $(which bedtools)"
echo "bedtools version: $(bedtools --version)"
echo ""

# Run extraction script
python 1_extract_encode_cCREs.py

echo ""
echo "========================================================================"
echo "EXTRACTION COMPLETE"
echo "========================================================================"
echo "Finished: $(date)"
