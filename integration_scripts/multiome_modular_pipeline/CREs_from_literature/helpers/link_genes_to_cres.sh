#!/bin/bash
#SBATCH --job-name=link_genes_to_cres
#SBATCH --output=logs/link_genes_to_cres.log
#SBATCH --error=logs/link_genes_to_cres.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Link Genes to CREs Analysis
#
# This script links genes to cis-regulatory elements (CREs) using pyBigWig.
#
# DEPENDENCIES:
# - Python script: link_genes_to_cres.py (main analysis script)
# - Conda environment: rna_seq_analysis_deep (contains pyBigWig and dependencies)
# - Input data: CRE coordinates and gene annotation files
#
# OUTPUT FILES:
# - logs/link_genes_to_cres.log: Execution log file
# - logs/link_genes_to_cres.err: Error log file
# - Analysis results from link_genes_to_cres.py
#
################################################################################

# Activate conda environment with pyBigWig
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo ""
echo "========================================================================"
echo "Running Python analysis..."
echo "========================================================================"
echo ""

# Run the Python script
python link_genes_to_cres.py
