#!/bin/bash
#SBATCH --job-name=encode_cCREs_pipeline
#SBATCH --output=logs/0_RUN_ENCODE_CCRES_ANALYSIS.log
#SBATCH --error=logs/0_RUN_ENCODE_CCRES_ANALYSIS.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# ENCODE cCREs Analysis Pipeline - Complete Pipeline
#
# This script runs the complete pipeline for extracting and visualizing
# ENCODE cCREs linked to splicing genes.
#
# Pipeline Steps:
# 1. Extract ENCODE cCREs linked to splicing genes
# 2. Convert TSV to BED format
# 3. Create heatmaps and metaprofiles with deepTools
# 4. Create custom comparison matrices
# 5. Create custom comparison visualizations
#
# NOTE: Emx1-Ctrl is EXCLUDED (failed sample) - using Nestin-Ctrl as control
#
# Prerequisites:
# - ../data/table_16.txt
# - ../data/mm10-cCREs.bed
# - BigWig files from Signac pipeline
#
# Usage:
#   sbatch 0_RUN_ENCODE_CCRES_ANALYSIS.sh
################################################################################

echo "========================================================================"
echo "ENCODE cCREs ANALYSIS PIPELINE - COMPLETE PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Change to script directory (use absolute path for SLURM compatibility)
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_encode_paper_intersection"
cd "$SCRIPT_DIR"

mkdir -p logs

# Activate conda
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh

# ============================================================================
# Step 1: Extract ENCODE cCREs linked to splicing genes
# ============================================================================
echo ""
echo "========================================================================"
echo "STEP 1: Extracting ENCODE cCREs linked to splicing genes"
echo "========================================================================"
echo ""

# Use sc-chromatin2 for bedtools (required for intersection)
conda activate sc-chromatin2
echo "Using environment: sc-chromatin2 (with bedtools)"
echo "bedtools version: $(bedtools --version)"

python 1_extract_encode_cCREs.py

if [ $? -ne 0 ]; then
    echo "ERROR: Step 1 failed!"
    exit 1
fi

echo "Step 1 completed successfully"

# ============================================================================
# Step 2: Convert to BED format
# ============================================================================
echo ""
echo "========================================================================"
echo "STEP 2: Converting TSV to BED format"
echo "========================================================================"
echo ""

# Switch to rna_seq_analysis_deep for pandas
conda activate rna_seq_analysis_deep
echo "Using environment: rna_seq_analysis_deep"

python 2_convert_encode_cCREs_to_bed.py

if [ $? -ne 0 ]; then
    echo "ERROR: Step 2 failed!"
    exit 1
fi

echo "Step 2 completed successfully"

# ============================================================================
# Step 3: Create heatmaps with deepTools
# ============================================================================
echo ""
echo "========================================================================"
echo "STEP 3: Creating heatmaps with deepTools"
echo "========================================================================"
echo ""

# rna_seq_analysis_deep has deepTools
bash 4_create_heatmaps_encode_cCREs.sh

if [ $? -ne 0 ]; then
    echo "WARNING: Step 3 had issues (may be due to missing BigWig files)"
fi

echo "Step 3 completed"

# ============================================================================
# Step 4: Create custom comparison matrices
# ============================================================================
echo ""
echo "========================================================================"
echo "STEP 4: Creating custom comparison matrices"
echo "========================================================================"
echo ""

bash 5_create_custom_comparisons.sh

if [ $? -ne 0 ]; then
    echo "WARNING: Step 4 had issues"
fi

echo "Step 4 completed"

# ============================================================================
# Step 5: Create custom comparison visualizations
# ============================================================================
echo ""
echo "========================================================================"
echo "STEP 5: Creating custom comparison visualizations"
echo "========================================================================"
echo ""

python 5_visualize_custom_comparisons.py --skip-individual

if [ $? -ne 0 ]; then
    echo "WARNING: Step 5 had issues"
fi

echo "Step 5 completed"

# ============================================================================
# Final Summary
# ============================================================================
echo ""
echo "========================================================================"
echo "PIPELINE COMPLETE!"
echo "========================================================================"
echo ""
echo "NOTE: Emx1-Ctrl is excluded (failed sample) - Nestin-Ctrl used as control"
echo ""
echo "Output files:"
echo ""
echo "Step 1 - Extraction:"
echo "  output/encode_cCREs_all_celltypes.tsv"
echo "  output/encode_cCREs_GABA.tsv"
echo "  output/encode_cCREs_by_type.tsv"
echo "  output/SUMMARY_encode_cCREs.txt"
echo ""
echo "Step 2 - BED files:"
echo "  output/encode_cCREs_all.bed"
echo "  output/encode_cCREs_GABA.bed"
echo "  output/encode_cCREs_{type}.bed"
echo ""
echo "Step 3 - Heatmaps:"
echo "  output/heatmaps_deeptools/heatmap_*.png"
echo "  output/heatmaps_deeptools/metaprofile_*.png"
echo ""
echo "Steps 4-5 - Custom Comparisons:"
echo "  output/custom_comparisons/matrix_*.gz"
echo "  output/custom_comparisons/profiles/metaprofile_*.png"
echo "  output/custom_comparisons/overview_all_conditions_*.png"
echo ""
echo "Comparisons performed (3 total):"
echo "  1. Nestin-Ctrl vs Nestin-Mut (within-genotype)"
echo "  2. Nestin-Ctrl vs Emx1-Mut (cross-genotype)"
echo "  3. Nestin-Mut vs Emx1-Mut (mutant comparison)"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
