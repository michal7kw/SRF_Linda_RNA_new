#!/bin/bash
#SBATCH --job-name=2_convert_to_bed
#SBATCH --output=logs/2_convert_to_bed.log
#SBATCH --error=logs/2_convert_to_bed.err
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=workq

################################################################################
# Convert GABA-Specific ENCODE cCRE TSV Files to BED Format
#
# Prerequisites:
# - output/GABA_specific_encode_cCREs.tsv
#
# Output:
# - output/GABA_specific_encode_cCREs.bed
# - output/GABA_specific_encode_cCREs_{type}.bed
################################################################################

echo "========================================================================"
echo "CONVERT GABA-SPECIFIC ENCODE cCRE TSV FILES TO BED FORMAT"
echo "========================================================================"
echo "Started: $(date)"
echo ""

mkdir -p logs

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

python 2_convert_to_bed.py

echo ""
echo "========================================================================"
echo "Completed: $(date)"
echo "========================================================================"
