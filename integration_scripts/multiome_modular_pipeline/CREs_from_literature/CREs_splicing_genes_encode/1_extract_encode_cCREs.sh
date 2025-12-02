#!/bin/bash
#SBATCH --job-name=1_extract_encode_cCREs
#SBATCH --output=logs/1_extract_encode_cCREs.log
#SBATCH --error=logs/1_extract_encode_cCREs.err
#SBATCH --time=30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=workq

################################################################################
# Extract ENCODE cCREs Associated with Splicing Genes
#
# This script identifies ENCODE cCREs that are proximal to splicing genes
# using bedtools window with a configurable window size.
#
# INPUT:
# - Splicing genes list (from Reactome/GO pathways)
# - ENCODE mm10-cCREs.bed
#
# OUTPUT:
# - CREs_splicing_genes_encode_all.tsv (all associations)
# - CREs_splicing_genes_encode_GABA.tsv (same, for GABA analysis)
# - CREs_splicing_genes_encode_by_type.tsv (with CRE type info)
# - SUMMARY_CREs_splicing_genes_encode.txt (summary report)
#
# Runtime: ~5 minutes (includes BioMart query)
#
# Usage:
#   sbatch 1_extract_encode_cCREs.sh
################################################################################

echo "========================================================================"
echo "EXTRACT ENCODE cCREs ASSOCIATED WITH SPLICING GENES"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_genes_encode"

# ============================================================================
# Activate Conda Environment
# ============================================================================

echo "Activating conda environment..."
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# Check bedtools is available
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found!"
    echo "Please ensure bedtools is in the environment."
    exit 1
fi
echo "bedtools found: $(bedtools --version)"
echo ""

# ============================================================================
# Run Analysis
# ============================================================================

echo "Running ENCODE cCRE extraction..."
echo ""

python 1_extract_encode_cCREs.py

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Analysis failed!"
    exit 1
fi

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "========================================================================"
echo "ANALYSIS COMPLETE!"
echo "========================================================================"
echo ""
echo "Output directory:"
echo "  output/"
echo ""
echo "Generated files:"
echo "  1. CREs_splicing_genes_encode_all.tsv"
echo "  2. CREs_splicing_genes_encode_GABA.tsv"
echo "  3. CREs_splicing_genes_encode_by_type.tsv"
echo "  4. SUMMARY_CREs_splicing_genes_encode.txt"
echo ""
echo "Review the summary file for statistics:"
echo "  cat output/SUMMARY_CREs_splicing_genes_encode.txt"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
