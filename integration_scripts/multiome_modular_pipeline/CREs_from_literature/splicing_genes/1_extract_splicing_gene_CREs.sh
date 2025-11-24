#!/bin/bash
#SBATCH --job-name=1_extract_splicing_gene_CREs
#SBATCH --output=logs/1_extract_splicing_gene_CREs.log
#SBATCH --error=logs/1_extract_splicing_gene_CREs.err
#SBATCH --time=30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=workq

################################################################################
# Extract CREs Associated with Splicing Genes
#
# This script extracts CREs (cis-regulatory elements) that are linked to
# splicing-related genes from Reactome pathways.
#
# INPUT:
# - Splicing genes list (293 genes from Reactome)
# - CRE-gene links (from ENCODE + co-accessibility analysis)
# - GABA CREs (hippocampal interneuron CREs)
# - GABA-specific CREs (unique to GABA)
#
# OUTPUT:
# - splicing_genes_CREs_all_celltypes.tsv (all cell types)
# - splicing_genes_CREs_GABA.tsv (GABA cell types)
# - splicing_genes_CREs_GABA_specific.tsv (GABA-specific)
# - SUMMARY_splicing_genes_CREs.txt (summary report)
#
# Runtime: ~1 minute
#
# Usage:
#   sbatch run_extract_splicing_gene_CREs.sh
################################################################################

echo "========================================================================"
echo "EXTRACT CREs ASSOCIATED WITH SPLICING GENES"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/splicing_genes"

# ============================================================================
# Activate Conda Environment
# ============================================================================

echo "Activating conda environment..."
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# ============================================================================
# Run Analysis
# ============================================================================

echo "Running splicing gene CRE extraction..."
echo ""

python 1_extract_splicing_gene_CREs.py

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
echo "  output/splicing_genes_analysis/"
echo ""
echo "Generated files:"
echo "  1. splicing_genes_CREs_all_celltypes.tsv"
echo "  2. splicing_genes_CREs_GABA.tsv"
echo "  3. splicing_genes_CREs_GABA_specific.tsv"
echo "  4. SUMMARY_splicing_genes_CREs.txt"
echo ""
echo "Review the summary file for statistics:"
echo "  cat output/splicing_genes_analysis/SUMMARY_splicing_genes_CREs.txt"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
