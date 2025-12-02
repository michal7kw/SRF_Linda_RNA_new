#!/bin/bash
#SBATCH --job-name=6_link_CREs_to_genes
#SBATCH --output=logs/6_link_CREs_to_genes.log
#SBATCH --error=logs/6_link_CREs_to_genes.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Link Cell-Type CREs to Genes Using Table 16
#
# This script links hippocampal interneuron and excitatory neuron CREs to genes
# using published CRE-gene correlations from ENCODE Table 16.
#
# Prerequisites:
# - output/hippocampal_interneuron_CREs.bed
# - output/excitatory_neuron_CREs.bed (optional)
# - ../data/table_16.txt
#
# Output:
# - output/hippocampal_interneuron_CREs_genes.tsv
# - output/excitatory_neuron_CREs_genes.tsv
# - output/gene_linkage_summary.txt
################################################################################

echo "========================================================================"
echo "LINK CELL-TYPE CREs TO GENES"
echo "========================================================================"
echo "Started: $(date)"
echo ""

mkdir -p logs

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

python 6_link_CREs_to_genes.py

echo ""
echo "========================================================================"
echo "Completed: $(date)"
echo "========================================================================"
