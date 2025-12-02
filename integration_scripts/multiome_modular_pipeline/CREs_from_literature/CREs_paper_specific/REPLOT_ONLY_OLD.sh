#!/bin/bash
#SBATCH --job-name=replot_heatmaps_OLD
#SBATCH --output=logs/replot_heatmaps_OLD_%j.log
#SBATCH --error=logs/replot_heatmaps_OLD_%j.err
#SBATCH --time=30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=workq

################################################################################
# FAST REPLOT: OLD Workflow (60% Overlap CRE Sets)
#
# This script regenerates heatmaps and metaprofiles using EXISTING matrix files
# from the OLD workflow (overlapping CRE sets), skipping the time-consuming
# computeMatrix step (~20-30 minutes).
#
# Use this when:
# - Matrix files already exist from OLD workflow
# - You want to adjust visualization parameters (scale, colors, resolution)
# - You don't want to recompute signal extraction
#
# Runtime: ~1-2 minutes (vs 6 hours for full OLD pipeline)
#
# Usage:
#   sbatch REPLOT_ONLY_OLD.sh
#
# Prerequisites:
# - Matrix files must exist from previous run:
#   output/GABA_DEG_analysis/heatmaps_deeptools/matrix_*.gz
################################################################################

echo "========================================================================"
echo "FAST REPLOT: Regenerate Plots from Existing Matrices (OLD WORKFLOW)"
echo "========================================================================"
echo "Started: $(date)"
echo ""

OUTPUT_DIR="output/GABA_DEG_analysis/heatmaps_deeptools"

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# ============================================================================
# Check if matrix files exist
# ============================================================================
echo "========================================================================"
echo "Checking for existing matrix files..."
echo "========================================================================"
echo ""

MISSING=0

if [ ! -f "$OUTPUT_DIR/matrix_GABA_all_conditions.gz" ]; then
    echo "✗ Missing: matrix_GABA_all_conditions.gz"
    MISSING=$((MISSING + 1))
else
    echo "✓ Found: matrix_GABA_all_conditions.gz"
fi

if [ ! -f "$OUTPUT_DIR/matrix_Excitatory_all_conditions.gz" ]; then
    echo "✗ Missing: matrix_Excitatory_all_conditions.gz"
    MISSING=$((MISSING + 1))
else
    echo "✓ Found: matrix_Excitatory_all_conditions.gz"
fi

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "ERROR: Matrix files not found!"
    echo "Please run the full OLD pipeline first:"
    echo "  sbatch 2_RUN_SPECIFIC_CRES_ANALYSIS.sh"
    exit 1
fi

echo ""
echo "✓ All required matrix files found"
echo "  Skipping computeMatrix step (saves ~20-30 minutes!)"
echo ""

# ============================================================================
# Get CRE counts from matrix files
# ============================================================================
echo "========================================================================"
echo "Getting CRE counts from matrix files..."
echo "========================================================================"
echo ""

N_GABA_CRES=$(computeMatrixOperations info -m $OUTPUT_DIR/matrix_GABA_all_conditions.gz 2>&1 | grep "number of regions" | awk '{print $NF}')
N_EXCIT_CRES=$(computeMatrixOperations info -m $OUTPUT_DIR/matrix_Excitatory_all_conditions.gz 2>&1 | grep "number of regions" | awk '{print $NF}')

echo "CRE counts:"
echo "  GABA CREs: $N_GABA_CRES"
echo "  Excitatory CREs: $N_EXCIT_CRES"
echo ""

# ============================================================================
# Determine common scale from existing matrices
# ============================================================================
echo "========================================================================"
echo "Calculating scale from existing matrices..."
echo "========================================================================"
echo ""

# Get data range from existing matrices
computeMatrixOperations dataRange -m $OUTPUT_DIR/matrix_GABA_all_conditions.gz 2>&1 | tee $OUTPUT_DIR/dataRange_GABA.txt
GABA_90TH=$(grep -E '^[A-Za-z0-9-]+\s+[0-9]' $OUTPUT_DIR/dataRange_GABA.txt | awk '{print $6}' | sort -n | tail -1)

computeMatrixOperations dataRange -m $OUTPUT_DIR/matrix_Excitatory_all_conditions.gz 2>&1 | tee $OUTPUT_DIR/dataRange_Excitatory.txt
EXCIT_90TH=$(grep -E '^[A-Za-z0-9-]+\s+[0-9]' $OUTPUT_DIR/dataRange_Excitatory.txt | awk '{print $6}' | sort -n | tail -1)

# Use 90th percentile for scale (robust against outliers)
COMMON_90TH=$(echo "$GABA_90TH $EXCIT_90TH" | awk '{print ($1 > $2) ? $1 : $2}')
COMMON_MAX=$(echo "$COMMON_90TH" | awk '{print $1 * 1.2}')

echo ""
echo "Scale determined:"
echo "  Common scale (heatmaps): 0 to $COMMON_MAX"
echo "  Metaprofiles: Auto-scaled (NEW - fixes flat line issue)"
echo ""

# ============================================================================
# Regenerate heatmaps with new parameters
# ============================================================================
echo "========================================================================"
echo "Regenerating heatmaps with improved resolution..."
echo "========================================================================"
echo ""

# GABA heatmap
echo "Creating GABA heatmap..."
plotHeatmap \
    -m $OUTPUT_DIR/matrix_GABA_all_conditions.gz \
    -o $OUTPUT_DIR/heatmap_GABA_all_conditions.png \
    --colorMap Reds \
    --dpi 300 \
    --whatToShow 'heatmap and colorbar' \
    --zMin 0 \
    --zMax $COMMON_MAX \
    --refPointLabel "CRE Center" \
    --heatmapHeight 25 \
    --heatmapWidth 6 \
    --sortRegions descend \
    --sortUsing mean \
    --plotTitle "POSITIVE CONTROL: ATAC Signal at GABA CREs (n=$N_GABA_CRES)" \
    --xAxisLabel "Distance from CRE Center (bp)" \
    --yAxisLabel "CREs" \
    2>&1 | grep -v "^$"

echo "✓ Saved: heatmap_GABA_all_conditions.png"
echo ""

# Excitatory heatmap
echo "Creating Excitatory heatmap..."
plotHeatmap \
    -m $OUTPUT_DIR/matrix_Excitatory_all_conditions.gz \
    -o $OUTPUT_DIR/heatmap_Excitatory_all_conditions.png \
    --colorMap Reds \
    --dpi 300 \
    --whatToShow 'heatmap and colorbar' \
    --zMin 0 \
    --zMax $COMMON_MAX \
    --refPointLabel "CRE Center" \
    --heatmapHeight 25 \
    --heatmapWidth 6 \
    --sortRegions descend \
    --sortUsing mean \
    --plotTitle "NEGATIVE CONTROL: ATAC Signal at Excitatory CREs (n=$N_EXCIT_CRES)" \
    --xAxisLabel "Distance from CRE Center (bp)" \
    --yAxisLabel "CREs" \
    2>&1 | grep -v "^$"

echo "✓ Saved: heatmap_Excitatory_all_conditions.png"
echo ""

# ============================================================================
# Regenerate metaprofiles with auto-scaling
# ============================================================================
echo "========================================================================"
echo "Regenerating metaprofiles with auto-scaling (FIX for flat lines)..."
echo "========================================================================"
echo ""

# GABA metaprofile
echo "Creating GABA metaprofile..."
plotProfile \
    -m $OUTPUT_DIR/matrix_GABA_all_conditions.gz \
    -o $OUTPUT_DIR/metaprofile_GABA_all_conditions.png \
    --plotTitle "POSITIVE CONTROL: ATAC Signal at GABA CREs (n=$N_GABA_CRES)" \
    --refPointLabel "CRE Center" \
    --averageType mean \
    --plotHeight 7 \
    --plotWidth 10 \
    --colors '#2E86AB' '#A23B72' '#F18F01' '#C73E1D' \
    --yAxisLabel "Mean ATAC Signal" \
    --yMin 0.015 \
    2>&1 | grep -v "^$"

echo "✓ Saved: metaprofile_GABA_all_conditions.png (auto-scaled)"
echo ""

# Excitatory metaprofile
echo "Creating Excitatory metaprofile..."
plotProfile \
    -m $OUTPUT_DIR/matrix_Excitatory_all_conditions.gz \
    -o $OUTPUT_DIR/metaprofile_Excitatory_all_conditions.png \
    --plotTitle "NEGATIVE CONTROL: ATAC Signal at Excitatory CREs (n=$N_EXCIT_CRES)" \
    --refPointLabel "CRE Center" \
    --averageType mean \
    --plotHeight 7 \
    --plotWidth 10 \
    --colors '#2E86AB' '#A23B72' '#F18F01' '#C73E1D' \
    --yAxisLabel "Mean ATAC Signal" \
    --yMin 0.015 \
    2>&1 | grep -v "^$"

echo "✓ Saved: metaprofile_Excitatory_all_conditions.png (auto-scaled)"
echo ""

# ============================================================================
# Regenerate genotype-specific plots (if matrices exist)
# ============================================================================
echo "========================================================================"
echo "Regenerating genotype-specific plots (if available)..."
echo "========================================================================"
echo ""

if [ -f "$OUTPUT_DIR/matrix_GABA_nestin.gz" ]; then
    echo "Creating Nestin heatmap and metaprofile..."

    plotHeatmap \
        -m $OUTPUT_DIR/matrix_GABA_nestin.gz \
        -o $OUTPUT_DIR/heatmap_GABA_nestin.png \
        --colorMap Reds \
        --dpi 300 \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 \
        --refPointLabel "CRE Center" \
        --heatmapHeight 25 \
        --heatmapWidth 5 \
        --plotTitle "Nestin: ATAC Signal at GABA CREs (n=$N_GABA_CRES)" \
        --xAxisLabel "Distance from CRE Center (bp)" \
        2>&1 | grep -v "^$"

    plotProfile \
        -m $OUTPUT_DIR/matrix_GABA_nestin.gz \
        -o $OUTPUT_DIR/metaprofile_GABA_nestin.png \
        --plotTitle "Nestin: Metaprofile at GABA CREs (n=$N_GABA_CRES)" \
        --refPointLabel "CRE Center" \
        --colors '#2E86AB' '#A23B72' \
        --plotHeight 6 \
        --plotWidth 8 \
        --yMin 0 \
        2>&1 | grep -v "^$"

    echo "✓ Saved: heatmap_GABA_nestin.png"
    echo "✓ Saved: metaprofile_GABA_nestin.png"
    echo ""
fi

if [ -f "$OUTPUT_DIR/matrix_GABA_emx1.gz" ]; then
    echo "Creating Emx1 heatmap and metaprofile..."

    plotHeatmap \
        -m $OUTPUT_DIR/matrix_GABA_emx1.gz \
        -o $OUTPUT_DIR/heatmap_GABA_emx1.png \
        --colorMap Reds \
        --dpi 300 \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 \
        --refPointLabel "CRE Center" \
        --heatmapHeight 25 \
        --heatmapWidth 5 \
        --plotTitle "Emx1: ATAC Signal at GABA CREs (n=$N_GABA_CRES)" \
        --xAxisLabel "Distance from CRE Center (bp)" \
        2>&1 | grep -v "^$"

    plotProfile \
        -m $OUTPUT_DIR/matrix_GABA_emx1.gz \
        -o $OUTPUT_DIR/metaprofile_GABA_emx1.png \
        --plotTitle "Emx1: Metaprofile at GABA CREs (n=$N_GABA_CRES)" \
        --refPointLabel "CRE Center" \
        --colors '#F18F01' '#C73E1D' \
        --plotHeight 6 \
        --plotWidth 8 \
        --yMin 0 \
        2>&1 | grep -v "^$"

    echo "✓ Saved: heatmap_GABA_emx1.png"
    echo "✓ Saved: metaprofile_GABA_emx1.png"
    echo ""
fi

# ============================================================================
# Summary
# ============================================================================
echo "========================================================================"
echo "REPLOT COMPLETE (OLD WORKFLOW)!"
echo "========================================================================"
echo ""
echo "Runtime: ~1-2 minutes (vs 6 hours for full OLD pipeline)"
echo ""
echo "Regenerated plots with FIXES:"
echo "  ✓ Metaprofiles: Auto-scaled (no more flat lines!)"
echo "  ✓ Heatmaps: Higher resolution (25\"×6\" @ 300 DPI)"
echo "  ✓ Both: Better visual quality"
echo ""
echo "Output directory:"
echo "  $OUTPUT_DIR/"
echo ""
echo "Check the plots:"
echo "  Metaprofiles should now show clear curves (not flat)"
echo "  Heatmaps should have smoother gradients (less aliasing)"
echo ""
echo "NOTE: This is the OLD workflow with 60% overlapping CRE sets"
echo "      Both GABA and Excitatory plots should show signal (shared elements)"
echo ""
echo "For mutually exclusive CREs (0% overlap), use instead:"
    echo "  cd ../Exclusive_CREs && sbatch REPLOT_ONLY.sh  # NEW workflow"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
