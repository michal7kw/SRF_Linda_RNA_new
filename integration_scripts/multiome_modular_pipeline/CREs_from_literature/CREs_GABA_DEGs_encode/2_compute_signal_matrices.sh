#!/bin/bash
#SBATCH --job-name=2_matrices_DEGs_encode
#SBATCH --output=logs/2_compute_signal_matrices.log
#SBATCH --error=logs/2_compute_signal_matrices.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Compute Signal Matrices for GABA DEG ENCODE cCREs
#
# This script uses deepTools to compute signal matrices at ENCODE cCREs
# associated with GABA up-regulated and down-regulated DEGs.
#
# Creates matrices for:
# - UP-DEGs: Nestin-Ctrl vs Nestin-Mut, Nestin-Ctrl vs Emx1-Mut, Nestin-Mut vs Emx1-Mut
# - DOWN-DEGs: Same comparisons
#
# NOTE: Emx1-Ctrl is excluded due to quality issues.
#
# Prerequisites:
# - output/encode_cCREs_up_DEGs.bed
# - output/encode_cCREs_down_DEGs.bed
# - BigWig files from Signac pipeline
#
# Output:
# - Matrix files for downstream visualization
# - Heatmaps and metaprofiles for each DEG set
################################################################################

echo "========================================================================"
echo "COMPUTE SIGNAL MATRICES FOR GABA DEG ENCODE cCREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_GABA_DEGs_encode"

# Configuration
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR="./output/heatmaps_deeptools"
BED_DIR="./output"

# Create directories
mkdir -p $OUTPUT_DIR
mkdir -p logs

# Parameters
WINDOW_SIZE=2000  # bp around CRE center (±2kb)
BIN_SIZE=50       # bp per bin
N_PROCESSORS=16

# Activate deepTools environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

# Check if deepTools is available
if ! command -v computeMatrix &> /dev/null; then
    echo "ERROR: deepTools not found!"
    echo "Install with: conda install -c bioconda deeptools"
    exit 1
fi

echo "deepTools found: $(computeMatrix --version)"
echo ""

# ============================================================================
# Step 1: Check BED files exist
# ============================================================================
echo "========================================================================"
echo "STEP 1: Checking BED files"
echo "========================================================================"
echo ""

BED_UP="$BED_DIR/encode_cCREs_up_DEGs.bed"
BED_DOWN="$BED_DIR/encode_cCREs_down_DEGs.bed"

if [ ! -f "$BED_UP" ]; then
    echo "ERROR: Up-DEG BED file not found: $BED_UP"
    echo "Please run 1_extract_encode_cCREs_for_DEGs.py first"
    exit 1
fi

if [ ! -f "$BED_DOWN" ]; then
    echo "ERROR: Down-DEG BED file not found: $BED_DOWN"
    echo "Please run 1_extract_encode_cCREs_for_DEGs.py first"
    exit 1
fi

N_UP=$(wc -l < "$BED_UP")
N_DOWN=$(wc -l < "$BED_DOWN")

echo "Up-DEG ENCODE cCREs: $N_UP CREs"
echo "Down-DEG ENCODE cCREs: $N_DOWN CREs"
echo ""

# ============================================================================
# Step 2: Check BigWig files
# ============================================================================
echo "========================================================================"
echo "STEP 2: Checking BigWig files"
echo "========================================================================"
echo ""

# Define BigWig files (excluding Emx1-Ctrl due to quality issues)
BW_NESTIN_CTRL="$BIGWIG_BASE/GABA_Nestin-Ctrl.bw"
BW_NESTIN_MUT="$BIGWIG_BASE/GABA_Nestin-Mut.bw"
BW_EMX1_MUT="$BIGWIG_BASE/GABA_Emx1-Mut.bw"

MISSING=0

for BW_FILE in "$BW_NESTIN_CTRL" "$BW_NESTIN_MUT" "$BW_EMX1_MUT"; do
    if [ -f "$BW_FILE" ]; then
        echo "  Found: $(basename $BW_FILE)"
    else
        echo "  MISSING: $BW_FILE"
        MISSING=$((MISSING + 1))
    fi
done

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "ERROR: Missing $MISSING BigWig files!"
    exit 1
fi

echo ""

# ============================================================================
# Step 3: Function to create analysis
# ============================================================================

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
    echo "=========================================="
    echo "Processing: $DESCRIPTION"
    echo "=========================================="
    echo "  CREs: $N_CRES"
    echo "  Output prefix: $OUTPUT_PREFIX"
    echo ""

    # All three conditions
    BIGWIGS="$BW_NESTIN_CTRL $BW_NESTIN_MUT $BW_EMX1_MUT"
    LABELS="Nestin-Ctrl Nestin-Mut Emx1-Mut"

    # Compute matrix
    echo "  Computing matrix (all 3 conditions)..."
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
        2>&1 | tee "$OUTPUT_DIR/computeMatrix_${OUTPUT_PREFIX}.log"

    if [ ! -f "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}.gz" ] || [ ! -s "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}.gz" ]; then
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

    if [ -f "$OUTPUT_DIR/heatmap_${OUTPUT_PREFIX}.png" ]; then
        echo "  Saved: heatmap_${OUTPUT_PREFIX}.png"
    else
        echo "  WARNING: Heatmap not created"
    fi

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

    if [ -f "$OUTPUT_DIR/metaprofile_${OUTPUT_PREFIX}.png" ]; then
        echo "  Saved: metaprofile_${OUTPUT_PREFIX}.png"
    else
        echo "  WARNING: Metaprofile not created"
    fi
    echo ""
}

# ============================================================================
# Step 4: Create analyses for UP-DEGs and DOWN-DEGs
# ============================================================================
echo "========================================================================"
echo "STEP 4: Creating analyses"
echo "========================================================================"

create_analysis "$BED_UP" "Up-Regulated GABA DEGs (ENCODE cCREs)" "up_DEGs"
create_analysis "$BED_DOWN" "Down-Regulated GABA DEGs (ENCODE cCREs)" "down_DEGs"

# ============================================================================
# Step 5: Create Nestin-only comparison matrices
# ============================================================================
echo "========================================================================"
echo "STEP 5: Creating Nestin-specific analyses"
echo "========================================================================"

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
    echo "Creating Nestin analysis: $DESCRIPTION"

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

create_nestin_analysis "$BED_UP" "Up-DEG ENCODE cCREs" "up_DEGs"
create_nestin_analysis "$BED_DOWN" "Down-DEG ENCODE cCREs" "down_DEGs"

# ============================================================================
# Step 6: Summary
# ============================================================================
echo ""
echo "========================================================================"
echo "ANALYSIS COMPLETE!"
echo "========================================================================"
echo ""
echo "Output directory: $OUTPUT_DIR/"
echo ""
echo "Generated files:"
ls -1 $OUTPUT_DIR/*.png 2>/dev/null | while read f; do
    echo "  - $(basename $f)"
done
echo ""
echo "Matrix files (for custom visualization):"
ls -1 $OUTPUT_DIR/matrix_*.tab 2>/dev/null | while read f; do
    echo "  - $(basename $f)"
done
echo ""
echo "Completed: $(date)"
echo "========================================================================"
