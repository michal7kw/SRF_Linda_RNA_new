#!/bin/bash
#SBATCH --job-name=1_extract_GABA_cCREs
#SBATCH --output=logs/1_extract_GABA_encode_cCREs.log
#SBATCH --error=logs/1_extract_GABA_encode_cCREs.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Extract GABA Cell Type Specific ENCODE cCREs
#
# This script identifies ENCODE cCREs that are specifically associated with
# GABA/hippocampal cell types using enhanced SubType filtering.
#
# Prerequisites:
# - ../data/table_16.txt (literature CRE-gene correlations)
# - ../data/mm10-cCREs.bed (ENCODE cCREs)
#
# Output (TWO CRE sets):
# Table 16-only (no ENCODE intersection):
# - output/GABA_specific_table16_cCREs.tsv
#
# ENCODE-intersected:
# - output/GABA_specific_encode_cCREs.tsv
# - output/GABA_specific_encode_cCREs_by_type.tsv
#
# Summary:
# - output/SUMMARY_GABA_specific_cCREs.txt
################################################################################

echo "========================================================================"
echo "EXTRACT GABA CELL TYPE SPECIFIC ENCODE cCREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

mkdir -p logs

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh

# Use sc-chromatin2 for bedtools (required for ENCODE intersection)
conda activate sc-chromatin2

# Verify bedtools is available
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found in sc-chromatin2!"
    exit 1
fi

echo "Using environment: sc-chromatin2"
echo "bedtools version: $(bedtools --version)"
echo ""

python 1_extract_GABA_encode_cCREs.py

echo ""
echo "========================================================================"
echo "Completed: $(date)"
echo "========================================================================"
