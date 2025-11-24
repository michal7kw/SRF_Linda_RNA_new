#!/bin/bash

################################################################################
# Quick Start Script - Genotype-Stratified Multiome Figure (L1 VERSION)
################################################################################
#
# This script runs both data preparation and plotting for genotype-stratified
# multiome heatmaps using L1 cell type annotations (9 broad categories).
#
# PREREQUISITES:
#   - Signac L1 pipeline must be completed (run_signac_pipeline_L1.sh)
#   - Output directory: signac_results_L1/celltype_results/
#   - Required files per genotype (Nestin, Emx1):
#       * DEG/{celltype}_{genotype}_DEGs.csv
#       * DA/{celltype}_{genotype}_DA_peaks.csv
#       * peak_gene_links/{celltype}_{genotype}_peak_gene_links.csv
#   - Integrated Seurat object: signac_results_L1/integrated_seurat_processed.rds
#
# OUTPUT:
#   - multiome_heatmap_data_L1_{Genotype}.rds (prepared data)
#   - output_genotypes_L1/multiome_heatmap_L1_{Genotype}.pdf
#   - output_genotypes_L1/multiome_heatmap_L1_{Genotype}.png
#   - output_genotypes_L1/chromatin_driven_degs_L1_{Genotype}.csv
#
# USAGE:
#   ./QUICK_START_GENOTYPES_L1.sh
#
################################################################################

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# Load conda environment (if needed)
# Uncomment if running interactively:
# source ~/.bashrc
# conda activate sc-chromatin2

echo "================================================================================"
echo "  MULTIOME FIGURE GENERATION: GENOTYPE-STRATIFIED (L1 VERSION)"
echo "================================================================================"
echo ""
echo "Cell type annotation level: L1 (9 broad categories)"
echo "  Categories: Excitatory, GABA, Oligo, Astrocytes, Immune,"
echo "              Vascular, Progenitors, Unknown, Ependymal"
echo ""
echo "Genotype-stratified analysis:"
echo "  - Separate heatmaps for Nestin and Emx1 samples"
echo "  - Shows chromatin-driven DEGs (DEG + DA peak + peak-gene linkage)"
echo ""
echo "Working directory: ${SCRIPT_DIR}"
echo "Start time: $(date)"
echo ""

################################################################################
# Step 1: Prepare Data for Each Genotype
################################################################################

echo "================================================================================"
echo "STEP 1: Prepare Multiome Heatmap Data (L1)"
echo "================================================================================"
echo ""

START_TIME=$(date +%s)

Rscript "${SCRIPT_DIR}/prepare_multiome_heatmap_data_genotypes_L1.R"

EXIT_CODE=$?
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "✗ Data preparation failed with exit code ${EXIT_CODE}"
    echo "Please check the error messages above."
    exit 1
fi

echo ""
echo "✓ Data preparation complete (elapsed time: ${ELAPSED}s)"
echo ""

################################################################################
# Step 2: Generate Heatmap Visualizations
################################################################################

echo "================================================================================"
echo "STEP 2: Generate Multiome Heatmaps (L1)"
echo "================================================================================"
echo ""

START_TIME=$(date +%s)

Rscript "${SCRIPT_DIR}/plot_multiome_heatmap_genotypes_L1.R"

EXIT_CODE=$?
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "✗ Heatmap generation failed with exit code ${EXIT_CODE}"
    echo "Please check the error messages above."
    exit 1
fi

echo ""
echo "✓ Heatmap generation complete (elapsed time: ${ELAPSED}s)"
echo ""

################################################################################
# Summary
################################################################################

echo "================================================================================"
echo "MULTIOME FIGURE GENERATION COMPLETE (L1 VERSION)!"
echo "================================================================================"
echo ""
echo "End time: $(date)"
echo ""
echo "Output directory: ${SCRIPT_DIR}/output_genotypes_L1/"
echo ""
echo "Generated files:"
echo "  Data preparation:"
echo "    - multiome_heatmap_data_L1_Nestin.rds"
echo "    - multiome_heatmap_data_L1_Emx1.rds"
echo "    - output_genotypes_L1/chromatin_driven_degs_L1_Nestin.csv"
echo "    - output_genotypes_L1/chromatin_driven_degs_L1_Emx1.csv"
echo ""
echo "  Visualizations:"
echo "    - output_genotypes_L1/multiome_heatmap_L1_Nestin.pdf"
echo "    - output_genotypes_L1/multiome_heatmap_L1_Nestin.png"
echo "    - output_genotypes_L1/multiome_heatmap_L1_Emx1.pdf"
echo "    - output_genotypes_L1/multiome_heatmap_L1_Emx1.png"
echo ""
echo "✓ All figure generation steps completed successfully!"
echo ""
echo "================================================================================"
