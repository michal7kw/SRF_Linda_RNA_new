#!/bin/bash
#SBATCH --job-name=0_CREs_encode_GABA_specific
#SBATCH --output=logs/0_CREs_encode_GABA_specific_%j.log
#SBATCH --error=logs/0_CREs_encode_GABA_specific_%j.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# GABA Cell Type Specific ENCODE cCREs - Complete Pipeline
#
# This script runs the complete pipeline for extracting and visualizing
# GABA cell type specific ENCODE cCREs.
#
# Pipeline Steps:
# 1. Extract GABA-specific CREs from Table 16 + mm10-cCREs.bed
# 2. Convert to BED format
# 3. Create heatmaps and metaprofiles with deepTools
# 4. Visualize differentially accessible CREs with minSig/minFC filtering
#
# Prerequisites:
# - ../data/table_16.txt
# - ../data/mm10-cCREs.bed
# - BigWig files from Signac pipeline
#
# Usage:
#   sbatch 0_RUN_COMPLETE_PIPELINE.sh
################################################################################

echo "========================================================================"
echo "GABA CELL TYPE SPECIFIC ENCODE cCREs - COMPLETE PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Change to script directory (use absolute path for SLURM compatibility)
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_encode_GABA_specific"
cd "$SCRIPT_DIR"

mkdir -p logs

# Activate conda
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh

# ============================================================================
# Step 1: Extract GABA-specific ENCODE cCREs
# ============================================================================
echo ""
echo "========================================================================"
echo "STEP 1: Extracting GABA-specific ENCODE cCREs"
echo "========================================================================"
echo ""

# Use sc-chromatin2 for bedtools (required for intersection)
conda activate sc-chromatin2
echo "Using environment: sc-chromatin2 (with bedtools)"
echo "bedtools version: $(bedtools --version)"

python 1_extract_GABA_encode_cCREs.py

if [ $? -ne 0 ]; then
    echo "ERROR: Step 1 failed!"
    exit 1
fi

# ============================================================================
# Step 2: Convert to BED format
# ============================================================================
echo ""
echo "========================================================================"
echo "STEP 2: Converting to BED format"
echo "========================================================================"
echo ""

# Switch to rna_seq_analysis_deep for pandas
conda activate rna_seq_analysis_deep
echo "Using environment: rna_seq_analysis_deep"

python 2_convert_to_bed.py

if [ $? -ne 0 ]; then
    echo "ERROR: Step 2 failed!"
    exit 1
fi

# ============================================================================
# Step 3: Create heatmaps with deepTools
# ============================================================================
echo ""
echo "========================================================================"
echo "STEP 3: Creating heatmaps with deepTools"
echo "========================================================================"
echo ""

bash 3_create_heatmaps.sh

if [ $? -ne 0 ]; then
    echo "WARNING: Step 3 had issues (may be due to missing BigWig files)"
fi

# ============================================================================
# Step 4: Visualize differentially accessible CREs
# ============================================================================
echo ""
echo "========================================================================"
echo "STEP 4: Visualizing differentially accessible CREs"
echo "========================================================================"
echo ""

# Run with minSig=2.0, minFC=2.0 WITH individual plots
echo "Running with thresholds: minSig=2.0, minFC=2.0 (WITH individual CRE plots)..."
python 4_visualize_DA_CREs.py --min-signal 2.0 --min-fc 2.0 --parallel 8 --max-individual 100

if [ $? -ne 0 ]; then
    echo "WARNING: Step 4 (minSig=2.0, minFC=2.0) had issues"
fi

# Run with stricter thresholds minSig=2.0, minFC=3.0 WITH individual plots
echo ""
echo "Running with stricter thresholds: minSig=2.0, minFC=3.0 (WITH individual CRE plots)..."
python 4_visualize_DA_CREs.py --min-signal 2.0 --min-fc 3.0 --parallel 8 --max-individual 100

if [ $? -ne 0 ]; then
    echo "WARNING: Step 4 (minSig=2.0, minFC=3.0) had issues"
fi

# ============================================================================
# Final Summary
# ============================================================================
echo ""
echo "========================================================================"
echo "PIPELINE COMPLETE!"
echo "========================================================================"
echo ""
echo "Output files:"
echo ""
echo "Step 1 - Extraction (TWO CRE SETS):"
echo "  Table 16-only (no ENCODE intersection):"
echo "    output/GABA_specific_table16_cCREs.tsv"
echo "  ENCODE-intersected:"
echo "    output/GABA_specific_encode_cCREs.tsv"
echo "    output/GABA_specific_encode_cCREs_by_type.tsv"
echo "  Summary:"
echo "    output/SUMMARY_GABA_specific_cCREs.txt"
echo ""
echo "Step 2 - BED files:"
echo "  Table 16-only:"
echo "    output/GABA_specific_table16_cCREs.bed"
echo "  ENCODE-intersected:"
echo "    output/GABA_specific_encode_cCREs.bed"
echo "    output/GABA_specific_encode_cCREs_{type}.bed"
echo ""
echo "Step 3 - Heatmaps:"
echo "  Table 16-only:"
echo "    output/heatmaps_deeptools/heatmap_table16_only.png"
echo "    output/heatmaps_deeptools/heatmap_table16_only_{nestin,emx1}.png"
echo "  ENCODE-intersected:"
echo "    output/heatmaps_deeptools/heatmap_GABA_specific.png"
echo "    output/heatmaps_deeptools/heatmap_GABA_specific_{nestin,emx1}.png"
echo ""
echo "Step 4 - DA Visualization (WITH individual CRE plots):"
echo "  output/DA_profiles_minSig2.0_minFC2.0/"
echo "    - metaprofiles for all comparisons"
echo "    - profiles/ subfolder with individual CRE plots"
echo "  output/DA_profiles_minSig2.0_minFC3.0/"
echo "    - metaprofiles for all comparisons"
echo "    - profiles/ subfolder with individual CRE plots"
echo ""
echo "RECOMMENDATION:"
echo "  Use TABLE 16-ONLY for maximum coverage (all GABA CREs)"
echo "  Use ENCODE-INTERSECTED for high-confidence validated CREs"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
