#!/bin/bash
#SBATCH --job-name=4_heatmaps_encode
#SBATCH --output=logs/4_create_heatmaps_encode_cCREs.log
#SBATCH --error=logs/4_create_heatmaps_encode_cCREs.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Create Heatmaps and Metaprofiles for ENCODE cCREs
#
# This script uses deepTools to create heatmaps and metaprofiles for ENCODE
# cCREs associated with splicing genes. Creates visualizations for:
# 1. All CREs combined
# 2. GABA cell type CREs
# 3. Each individual CRE type (dELS, pELS, CTCF-only, etc.)
#
# Prerequisites:
# - output/encode_cCREs_all.bed
# - output/encode_cCREs_GABA.bed
# - output/encode_cCREs_*.bed (type-specific BED files)
# - BigWig files from Signac pipeline
#
# Output:
# - Heatmaps and metaprofiles for each CRE category
# - Genotype-specific analyses (Nestin/Emx1)
################################################################################

echo "========================================================================"
echo "CREATE HEATMAPS FOR ENCODE cCREs USING DEEPTOOLS"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Change to script directory (use absolute path for SLURM compatibility)
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_encode_paper_intersection"
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

echo "✓ deepTools found: $(computeMatrix --version)"
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

# Find all BED files
echo "BED files found:"
for BED_FILE in $BED_DIR/encode_cCREs_*.bed; do
    if [ -f "$BED_FILE" ]; then
        BASENAME=$(basename "$BED_FILE")
        N_CRES=$(wc -l < "$BED_FILE")
        echo "  ✓ $BASENAME ($N_CRES CREs)"
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

# NOTE: Emx1-Ctrl is excluded (failed sample) - using only Nestin-Ctrl as control
for SAMPLE in "Nestin-Ctrl" "Nestin-Mut" "Emx1-Mut"; do
    BW_FILE="$BIGWIG_BASE/GABA_${SAMPLE}.bw"
    if [ -f "$BW_FILE" ]; then
        BIGWIGS="$BIGWIGS $BW_FILE"
        LABELS="$LABELS ${SAMPLE}"
        echo "  ✓ Found: GABA_${SAMPLE}.bw"
    else
        echo "  ✗ Missing: $BW_FILE"
        MISSING=$((MISSING + 1))
    fi
done

echo "  ⊗ Skipped: GABA_Emx1-Ctrl.bw (failed sample - excluded from analysis)"

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
        echo "  ⊗ Skipping: $DESCRIPTION (BED file not found)"
        return
    fi

    N_CRES=$(wc -l < "$BED_FILE")

    if [ $N_CRES -eq 0 ]; then
        echo "  ⊗ Skipping: $DESCRIPTION (no CREs)"
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

    # Check if matrix file was created (more reliable than exit code with pipes)
    if [ ! -f "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}.gz" ] || [ ! -s "$OUTPUT_DIR/matrix_${OUTPUT_PREFIX}.gz" ]; then
        echo "  ✗ ERROR: computeMatrix failed (matrix file not created)"
        return
    fi
    echo "  ✓ Matrix computed"

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
        echo "  ✓ Heatmap saved: heatmap_${OUTPUT_PREFIX}.png"
    else
        echo "  ✗ WARNING: Heatmap not created"
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
        echo "  ✓ Metaprofile saved: metaprofile_${OUTPUT_PREFIX}.png"
    else
        echo "  ✗ WARNING: Metaprofile not created"
    fi
    echo ""
}

# ============================================================================
# Step 4: Create analyses for ALL CREs
# ============================================================================
echo "========================================================================"
echo "STEP 4: Creating analyses for ALL CREs"
echo "========================================================================"

create_analysis "$BED_DIR/encode_cCREs_all.bed" "All ENCODE cCREs (All Cell Types)" "all_celltypes"

# ============================================================================
# Step 5: Create analyses for GABA CREs
# ============================================================================
echo "========================================================================"
echo "STEP 5: Creating analyses for GABA CREs"
echo "========================================================================"

create_analysis "$BED_DIR/encode_cCREs_GABA.bed" "ENCODE cCREs (GABA Cell Types)" "GABA"

# ============================================================================
# Step 6: Create analyses for each CRE type
# ============================================================================
echo "========================================================================"
echo "STEP 6: Creating analyses for each CRE type"
echo "========================================================================"

# Find all type-specific BED files
for BED_FILE in $BED_DIR/encode_cCREs_*.bed; do
    BASENAME=$(basename "$BED_FILE" .bed)

    # Skip all.bed and GABA.bed (already processed)
    if [[ "$BASENAME" == "encode_cCREs_all" ]] || [[ "$BASENAME" == "encode_cCREs_GABA" ]]; then
        continue
    fi

    # Extract CRE type from filename
    CRE_TYPE=${BASENAME#encode_cCREs_}
    DESCRIPTION="$CRE_TYPE CREs"

    create_analysis "$BED_FILE" "$DESCRIPTION" "$CRE_TYPE"
done

# ============================================================================
# Step 7: Create genotype-specific analyses for GABA CREs
# ============================================================================
echo "========================================================================"
echo "STEP 7: Creating genotype-specific analyses (GABA CREs)"
echo "========================================================================"
echo ""

BED_GABA="$BED_DIR/encode_cCREs_GABA.bed"

if [ ! -f "$BED_GABA" ]; then
    echo "  ⊗ Skipping genotype-specific analyses (GABA BED not found)"
else
    N_GABA_CRES=$(wc -l < "$BED_GABA")

    # Nestin only
    if [ -f "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" ] && [ -f "$BIGWIG_BASE/GABA_Nestin-Mut.bw" ]; then
        echo "Creating Nestin analysis (GABA CREs)..."
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
            -o "$OUTPUT_DIR/matrix_GABA_nestin.gz" \
            -p $N_PROCESSORS \
            2>&1

        plotHeatmap \
            -m "$OUTPUT_DIR/matrix_GABA_nestin.gz" \
            -o "$OUTPUT_DIR/heatmap_GABA_nestin.png" \
            --colorMap Reds \
            --dpi 300 \
            --whatToShow 'heatmap and colorbar' \
            --zMin 0 \
            --refPointLabel "CRE Center" \
            --heatmapHeight 10 \
            --heatmapWidth 3 \
            --plotTitle "Nestin: ATAC Signal at ENCODE cCREs (n=$N_GABA_CRES)" \
            --xAxisLabel "Distance from CRE Center (bp)" \
            2>&1

        if [ -f "$OUTPUT_DIR/heatmap_GABA_nestin.png" ]; then
            echo "  ✓ Saved: heatmap_GABA_nestin.png"
        fi

        plotProfile \
            -m "$OUTPUT_DIR/matrix_GABA_nestin.gz" \
            -o "$OUTPUT_DIR/metaprofile_GABA_nestin.png" \
            --plotTitle "Nestin: ATAC Signal at ENCODE cCREs (n=$N_GABA_CRES)" \
            --refPointLabel "CRE Center" \
            --colors '#2E86AB' '#A23B72' \
            --plotHeight 6 \
            --plotWidth 8 \
            2>&1

        if [ -f "$OUTPUT_DIR/metaprofile_GABA_nestin.png" ]; then
            echo "  ✓ Saved: metaprofile_GABA_nestin.png"
        fi
        echo ""
    fi

    # Emx1-Mut vs Nestin-Ctrl (using Nestin-Ctrl as control since Emx1-Ctrl failed)
    if [ -f "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" ] && [ -f "$BIGWIG_BASE/GABA_Emx1-Mut.bw" ]; then
        echo "Creating Emx1-Mut analysis using Nestin-Ctrl as control (GABA CREs)..."
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
            -o "$OUTPUT_DIR/matrix_GABA_emx1_mut.gz" \
            -p $N_PROCESSORS \
            2>&1

        plotHeatmap \
            -m "$OUTPUT_DIR/matrix_GABA_emx1_mut.gz" \
            -o "$OUTPUT_DIR/heatmap_GABA_emx1_mut.png" \
            --colorMap Reds \
            --dpi 300 \
            --whatToShow 'heatmap and colorbar' \
            --zMin 0 \
            --refPointLabel "CRE Center" \
            --heatmapHeight 10 \
            --heatmapWidth 3 \
            --plotTitle "Emx1-Mut vs Nestin-Ctrl: ATAC at ENCODE cCREs (n=$N_GABA_CRES)" \
            --xAxisLabel "Distance from CRE Center (bp)" \
            2>&1

        if [ -f "$OUTPUT_DIR/heatmap_GABA_emx1_mut.png" ]; then
            echo "  ✓ Saved: heatmap_GABA_emx1_mut.png"
        fi

        plotProfile \
            -m "$OUTPUT_DIR/matrix_GABA_emx1_mut.gz" \
            -o "$OUTPUT_DIR/metaprofile_GABA_emx1_mut.png" \
            --plotTitle "Emx1-Mut vs Nestin-Ctrl: ATAC at ENCODE cCREs (n=$N_GABA_CRES)" \
            --refPointLabel "CRE Center" \
            --colors '#2E86AB' '#C73E1D' \
            --plotHeight 6 \
            --plotWidth 8 \
            2>&1

        if [ -f "$OUTPUT_DIR/metaprofile_GABA_emx1_mut.png" ]; then
            echo "  ✓ Saved: metaprofile_GABA_emx1_mut.png"
        fi
        echo ""
    fi
fi

# ============================================================================
# Step 8: Create summary report
# ============================================================================
echo "========================================================================"
echo "STEP 8: Creating summary report"
echo "========================================================================"
echo ""

cat > $OUTPUT_DIR/README.txt << 'EOFR'
================================================================================
DEEPTOOLS HEATMAPS: ENCODE cCREs ASSOCIATED WITH SPLICING GENES
================================================================================

Analysis date: $(date)

PURPOSE:
--------
Visualize ATAC-seq signal at ENCODE cCREs (Candidate Cis-Regulatory Elements)
linked to splicing-related genes to understand chromatin accessibility patterns
at high-confidence regulatory elements.

SOURCE DATA:
-----------
ENCODE cCREs:
  ../data/mm10-cCREs.bed
  Mouse (mm10) candidate cis-regulatory elements from ENCODE

Splicing genes list:
  /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/splicing_genes/extracted_genes_final.csv

CRE-gene links:
  ../data/table_16.txt (literature correlations)

BigWig files:
  ../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw

CRE TYPES ANALYZED:
-------------------
- dELS: Distal enhancer-like signatures
- pELS: Proximal enhancer-like signatures
- CTCF-only: CTCF-bound sites
- DNase-H3K4me3: Promoter-like regions
- Multiple annotations: CREs with combined features

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
- Emx1-Ctrl vs Emx1-Mut

OUTPUT FILES:
-------------
All CREs:
  - heatmap_all_celltypes.png
  - metaprofile_all_celltypes.png

GABA CREs:
  - heatmap_GABA.png
  - metaprofile_GABA.png
  - heatmap_GABA_nestin.png / metaprofile_GABA_nestin.png
  - heatmap_GABA_emx1.png / metaprofile_GABA_emx1.png

Type-specific (for each CRE type):
  - heatmap_{type}.png
  - metaprofile_{type}.png

INTERPRETATION:
---------------
Heatmaps:
  - Red color = higher ATAC-seq signal (accessible chromatin)
  - Signal concentrated at CRE centers indicates active regulatory elements

Metaprofiles:
  - Compare Ctrl vs Mut to identify condition-specific accessibility changes
  - Compare CRE types to identify functional differences

Generated by: 4_create_heatmaps_encode_cCREs.sh
================================================================================
EOFR

echo "✓ Saved: README.txt"
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
ls -1 $OUTPUT_DIR/heatmap_*.png | while read f; do
    echo "  - $(basename $f)"
done
echo ""
echo "Review the heatmaps to identify ATAC signal patterns at"
echo "ENCODE cCREs across conditions, genotypes, and CRE types."
echo ""
echo "Completed: $(date)"
echo "========================================================================"
