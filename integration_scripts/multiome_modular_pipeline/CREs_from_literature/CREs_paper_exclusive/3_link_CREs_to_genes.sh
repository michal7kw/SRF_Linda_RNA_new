#!/bin/bash
#SBATCH --job-name=3_link_CREs_to_genes
#SBATCH --output=logs/3_link_CREs_to_genes.log
#SBATCH --error=logs/3_link_CREs_to_genes.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Link Cell-Type-Specific CREs to Genes Using Table 16
#
# This script links GABA-specific and Excitatory-specific CREs to genes
# using published CRE-gene correlations from ENCODE Table 16.
#
# Prerequisites:
# - output/GABA_specific_CREs.bed (from 1a_extract_cell_type_specific_CREs.py)
# - output/Excitatory_specific_CREs.bed
# - ../data/table_16.txt
#
# Output:
# - output/GABA_specific_CREs_genes.tsv
# - output/Excitatory_specific_CREs_genes.tsv
# - output/gene_linkage_summary.txt
################################################################################

echo "========================================================================"
echo "LINK CELL-TYPE-SPECIFIC CREs TO GENES"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Create logs directory
mkdir -p logs

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

# Run the Python script
python 3_link_CREs_to_genes.py

echo ""
echo "========================================================================"
echo "Completed: $(date)"
echo "========================================================================"
