#!/bin/bash
#SBATCH --job-name=3_heatmaps_GABA
#SBATCH --output=logs/3_create_heatmaps.log
#SBATCH --error=logs/3_create_heatmaps.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Create Heatmaps and Metaprofiles for GABA-Specific ENCODE cCREs
#
# This script uses deepTools to create heatmaps and metaprofiles for
# GABA cell type specific ENCODE cCREs. Creates visualizations for:
# 1. All GABA-specific CREs combined
# 2. Each individual CRE type (dELS, pELS, CTCF-only, etc.)
# 3. Genotype-specific analyses (Nestin/Emx1)
#
# Prerequisites:
# - output/GABA_specific_encode_cCREs.bed
# - output/GABA_specific_encode_cCREs_*.bed (type-specific BED files)
# - BigWig files from Signac pipeline
#
# Output:
# - Heatmaps and metaprofiles for each CRE category
# - Genotype-specific analyses (Nestin/Emx1)
################################################################################

echo "========================================================================"
echo "CREATE HEATMAPS FOR GABA-SPECIFIC ENCODE cCREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Change to script directory (use absolute path for SLURM compatibility)
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_encode_GABA_specific"
cd "$SCRIPT_DIR"

# Configuration
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR="./output/heatmaps_deeptools"
BED_DIR="./output"

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

# Create output directory
mkdir -p $OUTPUT_DIR

# ============================================================================
# Step 1: Check BED files exist
# ============================================================================
echo "========================================================================"
echo "STEP 1: Checking BED files"
echo "========================================================================"
echo ""

# Find all GABA-specific BED files
echo "BED files found:"
for BED_FILE in $BED_DIR/GABA_specific_encode_cCREs*.bed; do
    if [ -f "$BED_FILE" ]; then
        BASENAME=$(basename "$BED_FILE")
        N_CRES=$(wc -l < "$BED_FILE")
        echo "  $BASENAME ($N_CRES CREs)"
    fi
done
echo ""

# ============================================================================
# Step 2: Check BigWig files
# ============================================================================
echo "========================================================================"
echo "STEP 2: Checking BigWig files"
echo "========================================================================"
echo ""

# Collect BigWig files
BIGWIGS=""
LABELS=""
MISSING=0

# Define the samples to use (excluding failed Emx1-Ctrl)
declare -a SAMPLES=("Nestin-Ctrl" "Nestin-Mut" "Emx1-Mut")

for SAMPLE in "${SAMPLES[@]}"; do
    BW_FILE="$BIGWIG_BASE/GABA_${SAMPLE}.bw"
    if [ -f "$BW_FILE" ]; then
        BIGWIGS="$BIGWIGS $BW_FILE"
        LABELS="$LABELS ${SAMPLE}"
        echo "  Found: GABA_${SAMPLE}.bw"
    else
        echo "  Missing: $BW_FILE"
        MISSING=$((MISSING + 1))
    fi
done

echo ""
echo "NOTE: Emx1-Ctrl is excluded (failed sample). Using Nestin-Ctrl as control."

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "ERROR: Missing $MISSING BigWig files!"
    exit 1
fi

echo ""

# ============================================================================
# Step 3: Function to create heatmap and metaprofile
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

    # Step 3.1: Compute matrix
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
        2>&1 | tee "$OUTPUT_DIR/computeMatrix_${OUTPUT_PREFIX}.log"

    # Check if matrix file was created
    if [ ! -f "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}.gz" ] || [ ! -s "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}.gz" ]; then
        echo "  ERROR: computeMatrix failed (matrix file not created)"
        return
    fi
    echo "  Matrix computed"

    # Step 3.2: Create heatmap
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
        --heatmapWidth 6 \
        --sortRegions descend \
        --sortUsing mean \
        --plotTitle "ATAC Signal: $DESCRIPTION (n=$N_CRES)" \
        --xAxisLabel "Distance from CRE Center (bp)" \
        --yAxisLabel "CREs" \
        2>&1

    if [ -f "$OUTPUT_DIR/heatmap_${OUTPUT_PREFIX}.png" ]; then
        echo "  Heatmap saved: heatmap_${OUTPUT_PREFIX}.png"
    else
        echo "  WARNING: Heatmap not created"
    fi

    # Step 3.3: Create metaprofile
    echo "  Creating metaprofile..."
    plotProfile \
        -m "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}.gz" \
        -o "$OUTPUT_DIR/metaprofile_${OUTPUT_PREFIX}.png" \
        --plotTitle "ATAC Signal: $DESCRIPTION (n=$N_CRES)" \
        --refPointLabel "CRE Center" \
        --averageType mean \
        --plotHeight 7 \
        --plotWidth 10 \
        --colors '#2E86AB' '#A23B72' '#F18F01' '#C73E1D' \
        --yAxisLabel "Mean ATAC Signal" \
        2>&1

    if [ -f "$OUTPUT_DIR/metaprofile_${OUTPUT_PREFIX}.png" ]; then
        echo "  Metaprofile saved: metaprofile_${OUTPUT_PREFIX}.png"
    else
        echo "  WARNING: Metaprofile not created"
    fi
    echo ""
}

# ============================================================================
# Step 4: Create analyses for TABLE 16-ONLY CREs (no ENCODE intersection)
# ============================================================================
echo "========================================================================"
echo "STEP 4: Creating analyses for TABLE 16-ONLY CREs"
echo "========================================================================"

create_analysis "$BED_DIR/GABA_specific_table16_cCREs.bed" "Table 16-Only GABA CREs" "table16_only"

# ============================================================================
# Step 5: Create analyses for ENCODE-INTERSECTED CREs
# ============================================================================
echo "========================================================================"
echo "STEP 5: Creating analyses for ENCODE-INTERSECTED CREs"
echo "========================================================================"

create_analysis "$BED_DIR/GABA_specific_encode_cCREs.bed" "ENCODE-Intersected GABA CREs" "GABA_specific"

# ============================================================================
# Step 6: Create analyses for each CRE type (ENCODE-intersected only)
# ============================================================================
echo "========================================================================"
echo "STEP 6: Creating analyses for each CRE type (ENCODE-intersected)"
echo "========================================================================"

# Find all type-specific BED files
for BED_FILE in $BED_DIR/GABA_specific_encode_cCREs_*.bed; do
    if [ -f "$BED_FILE" ]; then
        BASENAME=$(basename "$BED_FILE" .bed)

        # Skip the main file (no type suffix)
        if [[ "$BASENAME" == "GABA_specific_encode_cCREs" ]]; then
            continue
        fi

        # Extract CRE type from filename
        CRE_TYPE=${BASENAME#GABA_specific_encode_cCREs_}
        DESCRIPTION="GABA-Specific $CRE_TYPE CREs"

        create_analysis "$BED_FILE" "$DESCRIPTION" "GABA_specific_$CRE_TYPE"
    fi
done

# ============================================================================
# Step 7: Create genotype-specific analyses for Table 16-only CREs
# ============================================================================
echo "========================================================================"
echo "STEP 7: Creating genotype-specific analyses (Table 16-only CREs)"
echo "========================================================================"
echo ""
echo "NOTE: Emx1-Ctrl is excluded (failed sample)."
echo "      Using Nestin-Ctrl as control for Emx1-Mut comparison."
echo ""

BED_TABLE16="$BED_DIR/GABA_specific_table16_cCREs.bed"

if [ -f "$BED_TABLE16" ]; then
    N_TABLE16_CRES=$(wc -l < "$BED_TABLE16")

    # Nestin: Ctrl vs Mut
    if [ -f "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" ] && [ -f "$BIGWIG_BASE/GABA_Nestin-Mut.bw" ]; then
        echo "Creating Nestin analysis (Table 16-only CREs)..."
        NESTIN_BIGWIGS="$BIGWIG_BASE/GABA_Nestin-Ctrl.bw $BIGWIG_BASE/GABA_Nestin-Mut.bw"
        NESTIN_LABELS="Nestin-Ctrl Nestin-Mut"

        computeMatrix reference-point \
            --referencePoint center \
            -b $WINDOW_SIZE -a $WINDOW_SIZE \
            -R "$BED_TABLE16" \
            -S $NESTIN_BIGWIGS \
            --samplesLabel $NESTIN_LABELS \
            --binSize $BIN_SIZE \
            --sortRegions descend \
            --sortUsing mean \
            --missingDataAsZero \
            -o "$OUTPUT_DIR/matrix_table16_only_nestin.gz" \
            -p $N_PROCESSORS \
            --outFileNameMatrix "$OUTPUT_DIR/matrix_table16_only_nestin.tab" \
            2>&1

        plotHeatmap \
            -m "$OUTPUT_DIR/matrix_table16_only_nestin.gz" \
            -o "$OUTPUT_DIR/heatmap_table16_only_nestin.png" \
            --colorMap Reds \
            --dpi 300 \
            --whatToShow 'heatmap and colorbar' \
            --zMin 0 \
            --refPointLabel "CRE Center" \
            --heatmapHeight 10 \
            --heatmapWidth 3 \
            --plotTitle "Nestin: ATAC at Table16-Only CREs (n=$N_TABLE16_CRES)" \
            --xAxisLabel "Distance from CRE Center (bp)" \
            2>&1

        if [ -f "$OUTPUT_DIR/heatmap_table16_only_nestin.png" ]; then
            echo "  Saved: heatmap_table16_only_nestin.png"
        fi

        plotProfile \
            -m "$OUTPUT_DIR/matrix_table16_only_nestin.gz" \
            -o "$OUTPUT_DIR/metaprofile_table16_only_nestin.png" \
            --plotTitle "Nestin: ATAC at Table16-Only CREs (n=$N_TABLE16_CRES)" \
            --refPointLabel "CRE Center" \
            --colors '#2E86AB' '#A23B72' \
            --plotHeight 6 \
            --plotWidth 8 \
            2>&1

        if [ -f "$OUTPUT_DIR/metaprofile_table16_only_nestin.png" ]; then
            echo "  Saved: metaprofile_table16_only_nestin.png"
        fi
        echo ""
    fi

    # Emx1-Mut vs Nestin-Ctrl (using Nestin-Ctrl as control since Emx1-Ctrl failed)
    if [ -f "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" ] && [ -f "$BIGWIG_BASE/GABA_Emx1-Mut.bw" ]; then
        echo "Creating Emx1-Mut vs Nestin-Ctrl analysis (Table 16-only CREs)..."
        EMX1_BIGWIGS="$BIGWIG_BASE/GABA_Nestin-Ctrl.bw $BIGWIG_BASE/GABA_Emx1-Mut.bw"
        EMX1_LABELS="Nestin-Ctrl Emx1-Mut"

        computeMatrix reference-point \
            --referencePoint center \
            -b $WINDOW_SIZE -a $WINDOW_SIZE \
            -R "$BED_TABLE16" \
            -S $EMX1_BIGWIGS \
            --samplesLabel $EMX1_LABELS \
            --binSize $BIN_SIZE \
            --sortRegions descend \
            --sortUsing mean \
            --missingDataAsZero \
            -o "$OUTPUT_DIR/matrix_table16_only_emx1.gz" \
            -p $N_PROCESSORS \
            2>&1

        plotHeatmap \
            -m "$OUTPUT_DIR/matrix_table16_only_emx1.gz" \
            -o "$OUTPUT_DIR/heatmap_table16_only_emx1.png" \
            --colorMap Reds \
            --dpi 300 \
            --whatToShow 'heatmap and colorbar' \
            --zMin 0 \
            --refPointLabel "CRE Center" \
            --heatmapHeight 10 \
            --heatmapWidth 3 \
            --plotTitle "Emx1-Mut vs Nestin-Ctrl: Table16 CREs (n=$N_TABLE16_CRES)" \
            --xAxisLabel "Distance from CRE Center (bp)" \
            2>&1

        if [ -f "$OUTPUT_DIR/heatmap_table16_only_emx1.png" ]; then
            echo "  Saved: heatmap_table16_only_emx1.png"
        fi

        plotProfile \
            -m "$OUTPUT_DIR/matrix_table16_only_emx1.gz" \
            -o "$OUTPUT_DIR/metaprofile_table16_only_emx1.png" \
            --plotTitle "Emx1-Mut vs Nestin-Ctrl: Table16 CREs (n=$N_TABLE16_CRES)" \
            --refPointLabel "CRE Center" \
            --colors '#2E86AB' '#C73E1D' \
            --plotHeight 6 \
            --plotWidth 8 \
            2>&1

        if [ -f "$OUTPUT_DIR/metaprofile_table16_only_emx1.png" ]; then
            echo "  Saved: metaprofile_table16_only_emx1.png"
        fi
        echo ""
    fi
else
    echo "  Skipping Table 16-only genotype analyses (BED file not found)"
fi

# ============================================================================
# Step 8: Create genotype-specific analyses for ENCODE-intersected CREs
# ============================================================================
echo "========================================================================"
echo "STEP 8: Creating genotype-specific analyses (ENCODE-intersected CREs)"
echo "========================================================================"
echo ""
echo "NOTE: Emx1-Ctrl is excluded (failed sample)."
echo "      Using Nestin-Ctrl as control for Emx1-Mut comparison."
echo ""

BED_GABA="$BED_DIR/GABA_specific_encode_cCREs.bed"

if [ ! -f "$BED_GABA" ]; then
    echo "  Skipping genotype-specific analyses (GABA BED not found)"
else
    N_GABA_CRES=$(wc -l < "$BED_GABA")

    # Nestin: Ctrl vs Mut
    if [ -f "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" ] && [ -f "$BIGWIG_BASE/GABA_Nestin-Mut.bw" ]; then
        echo "Creating Nestin analysis (GABA-Specific CREs)..."
        NESTIN_BIGWIGS="$BIGWIG_BASE/GABA_Nestin-Ctrl.bw $BIGWIG_BASE/GABA_Nestin-Mut.bw"
        NESTIN_LABELS="Nestin-Ctrl Nestin-Mut"

        computeMatrix reference-point \
            --referencePoint center \
            -b $WINDOW_SIZE -a $WINDOW_SIZE \
            -R "$BED_GABA" \
            -S $NESTIN_BIGWIGS \
            --samplesLabel $NESTIN_LABELS \
            --binSize $BIN_SIZE \
            --sortRegions descend \
            --sortUsing mean \
            --missingDataAsZero \
            -o "$OUTPUT_DIR/matrix_GABA_specific_nestin.gz" \
            -p $N_PROCESSORS \
            --outFileNameMatrix "$OUTPUT_DIR/matrix_GABA_specific_nestin.tab" \
            2>&1

        plotHeatmap \
            -m "$OUTPUT_DIR/matrix_GABA_specific_nestin.gz" \
            -o "$OUTPUT_DIR/heatmap_GABA_specific_nestin.png" \
            --colorMap Reds \
            --dpi 300 \
            --whatToShow 'heatmap and colorbar' \
            --zMin 0 \
            --refPointLabel "CRE Center" \
            --heatmapHeight 10 \
            --heatmapWidth 3 \
            --plotTitle "Nestin: ATAC at GABA-Specific CREs (n=$N_GABA_CRES)" \
            --xAxisLabel "Distance from CRE Center (bp)" \
            2>&1

        if [ -f "$OUTPUT_DIR/heatmap_GABA_specific_nestin.png" ]; then
            echo "  Saved: heatmap_GABA_specific_nestin.png"
        fi

        plotProfile \
            -m "$OUTPUT_DIR/matrix_GABA_specific_nestin.gz" \
            -o "$OUTPUT_DIR/metaprofile_GABA_specific_nestin.png" \
            --plotTitle "Nestin: ATAC at GABA-Specific CREs (n=$N_GABA_CRES)" \
            --refPointLabel "CRE Center" \
            --colors '#2E86AB' '#A23B72' \
            --plotHeight 6 \
            --plotWidth 8 \
            2>&1

        if [ -f "$OUTPUT_DIR/metaprofile_GABA_specific_nestin.png" ]; then
            echo "  Saved: metaprofile_GABA_specific_nestin.png"
        fi
        echo ""
    fi

    # Emx1-Mut vs Nestin-Ctrl (using Nestin-Ctrl as control since Emx1-Ctrl failed)
    if [ -f "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" ] && [ -f "$BIGWIG_BASE/GABA_Emx1-Mut.bw" ]; then
        echo "Creating Emx1-Mut vs Nestin-Ctrl analysis (using Nestin-Ctrl as control)..."
        EMX1_BIGWIGS="$BIGWIG_BASE/GABA_Nestin-Ctrl.bw $BIGWIG_BASE/GABA_Emx1-Mut.bw"
        EMX1_LABELS="Nestin-Ctrl Emx1-Mut"

        computeMatrix reference-point \
            --referencePoint center \
            -b $WINDOW_SIZE -a $WINDOW_SIZE \
            -R "$BED_GABA" \
            -S $EMX1_BIGWIGS \
            --samplesLabel $EMX1_LABELS \
            --binSize $BIN_SIZE \
            --sortRegions descend \
            --sortUsing mean \
            --missingDataAsZero \
            -o "$OUTPUT_DIR/matrix_GABA_specific_emx1.gz" \
            -p $N_PROCESSORS \
            --outFileNameMatrix "$OUTPUT_DIR/matrix_GABA_specific_emx1.tab" \
            2>&1

        plotHeatmap \
            -m "$OUTPUT_DIR/matrix_GABA_specific_emx1.gz" \
            -o "$OUTPUT_DIR/heatmap_GABA_specific_emx1.png" \
            --colorMap Reds \
            --dpi 300 \
            --whatToShow 'heatmap and colorbar' \
            --zMin 0 \
            --refPointLabel "CRE Center" \
            --heatmapHeight 10 \
            --heatmapWidth 3 \
            --plotTitle "Emx1-Mut vs Nestin-Ctrl: ATAC at GABA CREs (n=$N_GABA_CRES)" \
            --xAxisLabel "Distance from CRE Center (bp)" \
            2>&1

        if [ -f "$OUTPUT_DIR/heatmap_GABA_specific_emx1.png" ]; then
            echo "  Saved: heatmap_GABA_specific_emx1.png"
        fi

        plotProfile \
            -m "$OUTPUT_DIR/matrix_GABA_specific_emx1.gz" \
            -o "$OUTPUT_DIR/metaprofile_GABA_specific_emx1.png" \
            --plotTitle "Emx1-Mut vs Nestin-Ctrl: ATAC at GABA CREs (n=$N_GABA_CRES)" \
            --refPointLabel "CRE Center" \
            --colors '#2E86AB' '#C73E1D' \
            --plotHeight 6 \
            --plotWidth 8 \
            2>&1

        if [ -f "$OUTPUT_DIR/metaprofile_GABA_specific_emx1.png" ]; then
            echo "  Saved: metaprofile_GABA_specific_emx1.png"
        fi
        echo ""
    fi
fi

# ============================================================================
# Step 9: Create summary report
# ============================================================================
echo "========================================================================"
echo "STEP 9: Creating summary report"
echo "========================================================================"
echo ""

cat > $OUTPUT_DIR/README.txt << EOFR
================================================================================
DEEPTOOLS HEATMAPS: GABA CELL TYPE SPECIFIC CREs
================================================================================

Analysis date: $(date)

PURPOSE:
--------
Visualize ATAC-seq signal at GABA cell type specific CREs to validate
chromatin accessibility patterns at regulatory elements active in GABA neurons.

TWO CRE SETS ANALYZED:
----------------------
1. TABLE 16-ONLY: All GABA-specific CREs from Table 16 (no ENCODE intersection)
   - Maximum coverage of GABA regulatory regions
   - Uses coordinates directly from literature correlations

2. ENCODE-INTERSECTED: Table 16 CREs that overlap with mm10-cCREs.bed
   - High-confidence ENCODE-validated CREs
   - Includes CRE type annotations (dELS, pELS, etc.)

SOURCE DATA:
-----------
ENCODE cCREs:
  ../data/mm10-cCREs.bed
  Mouse (mm10) candidate cis-regulatory elements from ENCODE

CRE-gene links:
  ../data/table_16.txt (literature correlations, filtered for GABA SubTypes)

BigWig files:
  ../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw

GABA SUBTYPE KEYWORDS:
----------------------
Hippocampal regions: CA1, CA2, CA3, CA4, DG, DGNBL, HPF
GABA markers: GABA, GABAergic, INH, Interneuron, Inhibitory
Interneuron types: PV, PVALB, SST, VIP, LAMP5, LAMP, SNCG
Combined markers: PVGA, SSTGA, VIPGA, LAMGA
Granule cells: GRC, GC, Granule

PARAMETERS:
-----------
Window size: ±2000 bp around CRE center
Bin size: 50 bp
Reference point: CRE center
Sorting: By mean signal (descending)
Processors: 16 (parallel processing)

CONDITIONS:
-----------
- Nestin-Ctrl vs Nestin-Mut
- Nestin-Ctrl vs Emx1-Mut (Emx1-Ctrl excluded)

OUTPUT FILES:
-------------
TABLE 16-ONLY CREs:
  - heatmap_table16_only.png / metaprofile_table16_only.png
  - heatmap_table16_only_nestin.png / metaprofile_table16_only_nestin.png
  - heatmap_table16_only_emx1.png / metaprofile_table16_only_emx1.png

ENCODE-INTERSECTED CREs:
  - heatmap_GABA_specific.png / metaprofile_GABA_specific.png
  - heatmap_GABA_specific_nestin.png / metaprofile_GABA_specific_nestin.png
  - heatmap_GABA_specific_emx1.png / metaprofile_GABA_specific_emx1.png

CRE TYPE-SPECIFIC (ENCODE-intersected only):
  - heatmap_GABA_specific_{type}.png / metaprofile_GABA_specific_{type}.png

INTERPRETATION:
---------------
Heatmaps:
  - Red color = higher ATAC-seq signal (accessible chromatin)
  - Signal concentrated at CRE centers indicates active regulatory elements

Metaprofiles:
  - Compare Ctrl vs Mut to identify condition-specific accessibility changes
  - High signal expected at GABA-specific CREs (positive control)

RECOMMENDATION:
---------------
  - Use TABLE 16-ONLY for maximum coverage (all GABA CREs)
  - Use ENCODE-INTERSECTED for high-confidence validated CREs

Generated by: 3_create_heatmaps.sh
================================================================================
EOFR

echo "Saved: README.txt"
echo ""

# ============================================================================
# Final summary
# ============================================================================
echo "========================================================================"
echo "ANALYSIS COMPLETE!"
echo "========================================================================"
echo ""
echo "Output directory: $OUTPUT_DIR/"
echo ""
echo "Generated heatmaps and metaprofiles for:"
ls -1 $OUTPUT_DIR/heatmap_*.png 2>/dev/null | while read f; do
    echo "  - $(basename $f)"
done
echo ""
echo "Review the heatmaps to identify ATAC signal patterns at"
echo "GABA-specific ENCODE cCREs across conditions and genotypes."
echo ""
echo "Completed: $(date)"
echo "========================================================================"
