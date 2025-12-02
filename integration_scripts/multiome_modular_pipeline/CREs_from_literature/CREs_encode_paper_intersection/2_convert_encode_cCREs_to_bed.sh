#!/bin/bash
#SBATCH --job-name=convert_encode_bed
#SBATCH --output=logs/2_convert_encode_cCREs_to_bed.log
#SBATCH --error=logs/2_convert_encode_cCREs_to_bed.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH --partition=workq

# Convert ENCODE cCRE TSV Files to BED Format
# Creates BED files for all CREs, GABA CREs, and type-specific CREs

echo "========================================================================"
echo "CONVERT ENCODE cCRE TSV FILES TO BED FORMAT"
echo "========================================================================"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate sc-chromatin2

# Run conversion script
python 2_convert_encode_cCREs_to_bed.py

echo ""
echo "========================================================================"
echo "CONVERSION COMPLETE"
echo "========================================================================"
echo "Finished: $(date)"
