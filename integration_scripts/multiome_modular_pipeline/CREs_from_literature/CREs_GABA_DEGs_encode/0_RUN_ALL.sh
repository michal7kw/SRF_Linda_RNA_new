#!/bin/bash
#SBATCH --job-name=0_CREs_GABA_DEGs_encode
#SBATCH --output=logs/0_CREs_GABA_DEGs_encode_%j.log
#SBATCH --error=logs/0_CREs_GABA_DEGs_encode_%j.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# MASTER SCRIPT: GABA DEG ENCODE cCRE Analysis Pipeline
#
# This script runs the complete GABA DEG ENCODE cCRE analysis pipeline
# as a single SLURM job, executing all steps sequentially.
#
# PIPELINE OVERVIEW:
# 1. Extract ENCODE cCREs overlapping with CREs linked to GABA DEGs
# 2. Compute signal matrices using deepTools
# 3. Create custom comparison visualizations
#
# COMPARISONS (excluding Emx1-Ctrl due to quality issues):
# - Nestin-Ctrl vs Nestin-Mut (within-genotype mutation effect)
# - Nestin-Ctrl vs Emx1-Mut (cross-genotype mutation effect)
# - Nestin-Mut vs Emx1-Mut (mutant genotype comparison)
#
# NOTE: Emx1-Ctrl is a FAILED SAMPLE and is completely excluded.
#       Nestin-Ctrl serves as the reference control for all comparisons.
#
# USAGE:
#   sbatch 0_RUN_GABA_DEG_ENCODE_ANALYSIS.sh
#
# REQUIREMENTS:
# - GABA DEG lists (up/down regulated)
# - Table 16 CRE-gene correlations
# - ENCODE mm10-cCREs.bed
# - BigWig files from Signac pipeline (Nestin-Ctrl, Nestin-Mut, Emx1-Mut only)
# - bedtools, deepTools (via conda)
#
################################################################################

set -e  # Exit on error

echo "========================================================================"
echo "GABA DEG ENCODE cCRE ANALYSIS PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Change to script directory
cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_GABA_DEGs_encode"

# Create directories
mkdir -p logs
mkdir -p output

# ============================================================================
# Check prerequisites
# ============================================================================
echo "Checking prerequisites..."

# DEG files
DEG_UP="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/biomarkers/sig_deg_lists/GABA/Cluster_GABA_vs_Rest_up_significant.csv"
DEG_DOWN="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/biomarkers/sig_deg_lists/GABA/Cluster_GABA_vs_Rest_down_significant.csv"

if [ ! -f "$DEG_UP" ]; then
    echo "ERROR: Up-DEG file not found: $DEG_UP"
    exit 1
fi
if [ ! -f "$DEG_DOWN" ]; then
    echo "ERROR: Down-DEG file not found: $DEG_DOWN"
    exit 1
fi
echo "  DEG files found"

# Table 16
TABLE16="../data/table_16.txt"
if [ ! -f "$TABLE16" ]; then
    echo "ERROR: Table 16 not found: $TABLE16"
    exit 1
fi
echo "  Table 16 found"

# ENCODE cCREs
ENCODE_CCRES="../data/mm10-cCREs.bed"
if [ ! -f "$ENCODE_CCRES" ]; then
    echo "ERROR: ENCODE cCREs not found: $ENCODE_CCRES"
    exit 1
fi
echo "  ENCODE cCREs found"

# BigWig files (only 3 conditions - Emx1-Ctrl excluded)
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
for BW in GABA_Nestin-Ctrl.bw GABA_Nestin-Mut.bw GABA_Emx1-Mut.bw; do
    if [ ! -f "$BIGWIG_BASE/$BW" ]; then
        echo "ERROR: BigWig file not found: $BIGWIG_BASE/$BW"
        exit 1
    fi
done
echo "  BigWig files found (3 conditions: Nestin-Ctrl, Nestin-Mut, Emx1-Mut)"

echo ""
echo "All prerequisites satisfied!"
echo ""

# ============================================================================
# STEP 1: Extract ENCODE cCREs for GABA DEGs
# ============================================================================
echo "========================================================================"
echo "STEP 1: Extract ENCODE cCREs for GABA DEGs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Activate conda environment with bedtools and pandas
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate sc-chromatin2

echo "Conda environment: sc-chromatin2"
echo "Checking bedtools: $(which bedtools)"
echo ""

# Run extraction script
python 1_extract_encode_cCREs_for_DEGs.py

echo ""
echo "Step 1 completed: $(date)"
echo ""

# ============================================================================
# STEP 2: Compute Signal Matrices
# ============================================================================
echo "========================================================================"
echo "STEP 2: Compute Signal Matrices"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Switch to deepTools environment
conda activate rna_seq_analysis_deep

echo "Conda environment: rna_seq_analysis_deep"
echo "deepTools version: $(computeMatrix --version 2>&1 || echo 'checking...')"
echo ""

# Source and run the matrix computation script directly (without sbatch)
# We'll inline the key parts here to avoid environment issues

BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR="./output/heatmaps_deeptools"
BED_DIR="./output"

mkdir -p $OUTPUT_DIR

# Parameters
WINDOW_SIZE=2000
BIN_SIZE=50
N_PROCESSORS=16

# BigWig files (Emx1-Ctrl excluded)
BW_NESTIN_CTRL="$BIGWIG_BASE/GABA_Nestin-Ctrl.bw"
BW_NESTIN_MUT="$BIGWIG_BASE/GABA_Nestin-Mut.bw"
BW_EMX1_MUT="$BIGWIG_BASE/GABA_Emx1-Mut.bw"

# BED files
BED_UP="$BED_DIR/encode_cCREs_up_DEGs.bed"
BED_DOWN="$BED_DIR/encode_cCREs_down_DEGs.bed"

# Function to create analysis for a DEG set
create_analysis() {
    local BED_FILE=$1
    local DESCRIPTION=$2
    local OUTPUT_PREFIX=$3

    if [ ! -f "$BED_FILE" ]; then
        echo "  Skipping: $DESCRIPTION (BED file not found)"
        return
    fi

    N_CRES=$(wc -l < "$BED_FILE")

    if [ $N_CRES -eq 0 ]; then
        echo "  Skipping: $DESCRIPTION (no CREs)"
        return
    fi

    echo ""
    echo "Processing: $DESCRIPTION"
    echo "  CREs: $N_CRES"
    echo ""

    # All three conditions (Emx1-Ctrl excluded)
    BIGWIGS="$BW_NESTIN_CTRL $BW_NESTIN_MUT $BW_EMX1_MUT"
    LABELS="Nestin-Ctrl Nestin-Mut Emx1-Mut"

    # Compute matrix
    echo "  Computing matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b $WINDOW_SIZE -a $WINDOW_SIZE \
        -R "$BED_FILE" \
        -S $BIGWIGS \
        --samplesLabel $LABELS \
        --binSize $BIN_SIZE \
        --sortRegions keep \
        --missingDataAsZero \
        -o "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}.gz" \
        -p $N_PROCESSORS \
        --outFileNameMatrix "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}.tab" \
        2>&1

    if [ ! -f "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}.gz" ]; then
        echo "  ERROR: computeMatrix failed"
        return
    fi
    echo "  Matrix computed successfully"

    # Create heatmap
    echo "  Creating heatmap..."
    plotHeatmap \
        -m "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}.gz" \
        -o "$OUTPUT_DIR/heatmap_${OUTPUT_PREFIX}.png" \
        --colorMap Reds \
        --dpi 300 \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 \
        --refPointLabel "CRE Center" \
        --heatmapHeight 10 \
        --heatmapWidth 4 \
        --sortRegions descend \
        --sortUsing mean \
        --plotTitle "ATAC Signal: $DESCRIPTION (n=$N_CRES)" \
        --xAxisLabel "Distance from CRE Center (bp)" \
        --yAxisLabel "CREs" \
        2>&1

    # Create metaprofile
    echo "  Creating metaprofile..."
    plotProfile \
        -m "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}.gz" \
        -o "$OUTPUT_DIR/metaprofile_${OUTPUT_PREFIX}.png" \
        --plotTitle "ATAC Signal: $DESCRIPTION (n=$N_CRES)" \
        --refPointLabel "CRE Center" \
        --averageType mean \
        --plotHeight 7 \
        --plotWidth 10 \
        --colors '#2E86AB' '#A23B72' '#C73E1D' \
        --yAxisLabel "Mean ATAC Signal" \
        2>&1

    echo "  Saved: heatmap_${OUTPUT_PREFIX}.png, metaprofile_${OUTPUT_PREFIX}.png"
}

# Create Nestin-only analysis
create_nestin_analysis() {
    local BED_FILE=$1
    local DESCRIPTION=$2
    local OUTPUT_PREFIX=$3

    if [ ! -f "$BED_FILE" ]; then
        return
    fi

    N_CRES=$(wc -l < "$BED_FILE")
    if [ $N_CRES -eq 0 ]; then
        return
    fi

    echo ""
    echo "Creating Nestin-only analysis: $DESCRIPTION"

    NESTIN_BIGWIGS="$BW_NESTIN_CTRL $BW_NESTIN_MUT"
    NESTIN_LABELS="Nestin-Ctrl Nestin-Mut"

    computeMatrix reference-point \
        --referencePoint center \
        -b $WINDOW_SIZE -a $WINDOW_SIZE \
        -R "$BED_FILE" \
        -S $NESTIN_BIGWIGS \
        --samplesLabel $NESTIN_LABELS \
        --binSize $BIN_SIZE \
        --sortRegions descend \
        --sortUsing mean \
        --missingDataAsZero \
        -o "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}_nestin.gz" \
        -p $N_PROCESSORS \
        --outFileNameMatrix "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}_nestin.tab" \
        2>&1

    if [ -f "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}_nestin.gz" ]; then
        plotHeatmap \
            -m "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}_nestin.gz" \
            -o "$OUTPUT_DIR/heatmap_${OUTPUT_PREFIX}_nestin.png" \
            --colorMap Reds \
            --dpi 300 \
            --whatToShow 'heatmap and colorbar' \
            --zMin 0 \
            --refPointLabel "CRE Center" \
            --heatmapHeight 10 \
            --heatmapWidth 3 \
            --plotTitle "Nestin: $DESCRIPTION (n=$N_CRES)" \
            --xAxisLabel "Distance from CRE Center (bp)" \
            2>&1

        plotProfile \
            -m "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}_nestin.gz" \
            -o "$OUTPUT_DIR/metaprofile_${OUTPUT_PREFIX}_nestin.png" \
            --plotTitle "Nestin: $DESCRIPTION (n=$N_CRES)" \
            --refPointLabel "CRE Center" \
            --colors '#2E86AB' '#A23B72' \
            --plotHeight 6 \
            --plotWidth 8 \
            2>&1

        echo "  Saved: heatmap_${OUTPUT_PREFIX}_nestin.png, metaprofile_${OUTPUT_PREFIX}_nestin.png"
    fi
}

# Run analyses
create_analysis "$BED_UP" "Up-Regulated GABA DEGs (ENCODE cCREs)" "up_DEGs"
create_analysis "$BED_DOWN" "Down-Regulated GABA DEGs (ENCODE cCREs)" "down_DEGs"
create_nestin_analysis "$BED_UP" "Up-DEG ENCODE cCREs" "up_DEGs"
create_nestin_analysis "$BED_DOWN" "Down-DEG ENCODE cCREs" "down_DEGs"

echo ""
echo "Step 2 completed: $(date)"
echo ""

# ============================================================================
# STEP 3: Create Custom Comparison Visualizations
# ============================================================================
echo "========================================================================"
echo "STEP 3: Create Custom Comparison Visualizations"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Switch back to sc-chromatin2 for visualization
conda activate sc-chromatin2

echo "Conda environment: sc-chromatin2"
echo ""

# Run visualization script
python 3_visualize_deg_comparisons.py

echo ""
echo "Step 3 completed: $(date)"
echo ""

# ============================================================================
# Summary
# ============================================================================
echo "========================================================================"
echo "PIPELINE COMPLETE!"
echo "========================================================================"
echo ""
echo "Finished: $(date)"
echo ""
echo "Output files:"
echo "  ./output/                              - TSV and BED files"
echo "  ./output/heatmaps_deeptools/           - Heatmaps and metaprofiles"
echo "  ./output/custom_comparisons_*/         - Custom comparison plots"
echo ""
echo "Key files:"
ls -la output/*.bed output/*.tsv 2>/dev/null || echo "  (BED/TSV files)"
echo ""
echo "Visualizations:"
ls -la output/heatmaps_deeptools/*.png 2>/dev/null | head -10 || echo "  (PNG files in heatmaps_deeptools/)"
echo ""
echo "========================================================================"
