#!/bin/bash
#SBATCH --job-name=4_create_heatmaps
#SBATCH --output=logs/4_create_heatmaps.log
#SBATCH --error=logs/4_create_heatmaps.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=workq

################################################################################
# Create Heatmaps and Metaprofiles for Splicing Gene ENCODE cCREs
#
# This script uses deepTools to create heatmaps and metaprofiles for ENCODE cCREs
# associated with splicing-related genes.
#
# Prerequisites:
# - output/CREs_splicing_genes_encode_all.bed
# - output/CREs_splicing_genes_encode_GABA.bed
# - BigWig files from Signac pipeline (signac_results_L1/bigwig_tracks_L1/)
#
# Output:
# - Heatmaps showing ATAC signal at splicing gene ENCODE cCREs
# - Metaprofiles showing average signal profiles
# - Genotype-specific analyses (Nestin/Emx1 separate)
################################################################################

echo "========================================================================"
echo "CREATE HEATMAPS FOR SPLICING GENE ENCODE cCREs USING DEEPTOOLS"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_genes_encode"

# Configuration
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR="./output/heatmaps_deeptools"
BED_ALL="./output/CREs_splicing_genes_encode_all.bed"
BED_GABA="./output/CREs_splicing_genes_encode_GABA.bed"

# Parameters
WINDOW_SIZE=2000  # bp around CRE center (+/-2kb)
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

# Check all BED file
if [ ! -f "$BED_ALL" ]; then
    echo "ERROR: BED file not found: $BED_ALL"
    echo "Please run: sbatch 2_convert_to_bed.sh"
    exit 1
fi

N_ALL_CRES=$(wc -l < $BED_ALL)
echo "Found BED file: $BED_ALL"
echo "  Total CREs: $N_ALL_CRES"
echo ""

# Check GABA BED file (same as all for this pipeline)
if [ ! -f "$BED_GABA" ]; then
    echo "WARNING: GABA BED file not found: $BED_GABA"
    echo "Using all CREs BED file instead."
    BED_GABA=$BED_ALL
fi

N_GABA_CRES=$(wc -l < $BED_GABA)
echo "GABA BED file: $BED_GABA"
echo "  Total CREs: $N_GABA_CRES"
echo ""

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Copy BED files to output directory for reference
cp $BED_ALL $OUTPUT_DIR/CREs_splicing_genes_encode_all.bed
cp $BED_GABA $OUTPUT_DIR/CREs_splicing_genes_encode_GABA.bed
echo "Copied BED files to output directory"
echo ""

# Show CRE IDs (first 20)
echo "CRE IDs in analysis (first 20):"
head -n 20 $BED_ALL | awk '{print "  - " $4 " (" $1 ":" $2 "-" $3 ")"}'
if [ $N_ALL_CRES -gt 20 ]; then
    echo "  ... and $((N_ALL_CRES - 20)) more"
fi
echo ""

# ============================================================================
# Step 2: Check BigWig files (Emx1-Ctrl excluded - failed sample)
# ============================================================================
echo "========================================================================"
echo "STEP 2: Checking BigWig files"
echo "========================================================================"
echo ""
echo "NOTE: Emx1-Ctrl is excluded (failed sample)"
echo ""

# Collect BigWig files - only valid samples (excluding Emx1-Ctrl)
BIGWIGS=""
LABELS=""
MISSING=0

# Check required samples: Nestin-Ctrl, Nestin-Mut, Emx1-Mut
for SAMPLE in "Nestin-Ctrl" "Nestin-Mut" "Emx1-Mut"; do
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

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "ERROR: Missing $MISSING BigWig files!"
    exit 1
fi

echo ""

# ============================================================================
# Step 3: Run computeMatrix for ALL CREs
# ============================================================================
echo "========================================================================"
echo "STEP 3: Computing signal matrices for ALL CREs"
echo "========================================================================"
echo ""

echo "Running computeMatrix on all CREs..."
echo "  Window size: +/-${WINDOW_SIZE} bp"
echo "  Bin size: ${BIN_SIZE} bp"
echo "  Processors: $N_PROCESSORS"
echo "  CREs: $N_ALL_CRES"
echo ""

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_ALL \
    -S $BIGWIGS \
    --samplesLabel $LABELS \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_all.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_all.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_all.log

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for all CREs"
    exit 1
fi

echo ""
echo "All CREs matrix computed and saved"
echo ""

# ============================================================================
# Step 4: Run computeMatrix for GABA CREs
# ============================================================================
echo "========================================================================"
echo "STEP 4: Computing signal matrices for GABA CREs"
echo "========================================================================"
echo ""

echo "Running computeMatrix on GABA CREs..."
echo "  Window size: +/-${WINDOW_SIZE} bp"
echo "  Bin size: ${BIN_SIZE} bp"
echo "  Processors: $N_PROCESSORS"
echo "  CREs: $N_GABA_CRES"
echo ""

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_GABA \
    -S $BIGWIGS \
    --samplesLabel $LABELS \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_GABA.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_GABA.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_GABA.log

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for GABA CREs"
    exit 1
fi

echo ""
echo "GABA CREs matrix computed and saved"
echo ""

# ============================================================================
# Step 5: Create heatmaps
# ============================================================================
echo "========================================================================"
echo "STEP 5: Creating heatmaps"
echo "========================================================================"
echo ""

# All CREs heatmap
echo "Creating all CREs heatmap..."
plotHeatmap \
    -m $OUTPUT_DIR/matrix_all.gz \
    -o $OUTPUT_DIR/heatmap_all.png \
    --colorMap Reds \
    --dpi 300 \
    --whatToShow 'heatmap and colorbar' \
    --zMin 0 \
    --refPointLabel "CRE Center" \
    --heatmapHeight 10 \
    --heatmapWidth 6 \
    --sortRegions descend \
    --sortUsing mean \
    --plotTitle "ATAC Signal at Splicing Gene ENCODE cCREs - All (n=$N_ALL_CRES)" \
    --xAxisLabel "Distance from CRE Center (bp)" \
    --yAxisLabel "CREs" \
    2>&1 | grep -v "^$"

echo "Saved: heatmap_all.png"

# GABA CREs heatmap
echo "Creating GABA CREs heatmap..."
plotHeatmap \
    -m $OUTPUT_DIR/matrix_GABA.gz \
    -o $OUTPUT_DIR/heatmap_GABA.png \
    --colorMap Reds \
    --dpi 300 \
    --whatToShow 'heatmap and colorbar' \
    --zMin 0 \
    --refPointLabel "CRE Center" \
    --heatmapHeight 10 \
    --heatmapWidth 6 \
    --sortRegions descend \
    --sortUsing mean \
    --plotTitle "ATAC Signal at Splicing Gene ENCODE cCREs - GABA (n=$N_GABA_CRES)" \
    --xAxisLabel "Distance from CRE Center (bp)" \
    --yAxisLabel "CREs" \
    2>&1 | grep -v "^$"

echo "Saved: heatmap_GABA.png"
echo ""

# ============================================================================
# Step 6: Create metaprofiles
# ============================================================================
echo "========================================================================"
echo "STEP 6: Creating metaprofiles"
echo "========================================================================"
echo ""

# All CREs metaprofile
echo "Creating all CREs metaprofile..."
plotProfile \
    -m $OUTPUT_DIR/matrix_all.gz \
    -o $OUTPUT_DIR/metaprofile_all.png \
    --plotTitle "ATAC Signal at Splicing Gene ENCODE cCREs - All (n=$N_ALL_CRES)" \
    --refPointLabel "CRE Center" \
    --averageType mean \
    --plotHeight 7 \
    --plotWidth 10 \
    --colors '#2E86AB' '#A23B72' '#F18F01' '#C73E1D' \
    --yAxisLabel "Mean ATAC Signal" \
    2>&1 | grep -v "^$"

echo "Saved: metaprofile_all.png"

# GABA CREs metaprofile
echo "Creating GABA CREs metaprofile..."
plotProfile \
    -m $OUTPUT_DIR/matrix_GABA.gz \
    -o $OUTPUT_DIR/metaprofile_GABA.png \
    --plotTitle "ATAC Signal at Splicing Gene ENCODE cCREs - GABA (n=$N_GABA_CRES)" \
    --refPointLabel "CRE Center" \
    --averageType mean \
    --plotHeight 7 \
    --plotWidth 10 \
    --colors '#2E86AB' '#A23B72' '#F18F01' '#C73E1D' \
    --yAxisLabel "Mean ATAC Signal" \
    2>&1 | grep -v "^$"

echo "Saved: metaprofile_GABA.png"
echo ""

# ============================================================================
# Step 7: Create genotype-specific analyses
# ============================================================================
echo "========================================================================"
echo "STEP 7: Creating genotype-specific analyses"
echo "========================================================================"
echo ""

# Nestin only
if [ -f "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" ] && [ -f "$BIGWIG_BASE/GABA_Nestin-Mut.bw" ]; then
    echo "Creating Nestin analysis (GABA CREs)..."
    NESTIN_BIGWIGS="$BIGWIG_BASE/GABA_Nestin-Ctrl.bw $BIGWIG_BASE/GABA_Nestin-Mut.bw"
    NESTIN_LABELS="Nestin-Ctrl Nestin-Mut"

    computeMatrix reference-point \
        --referencePoint center \
        -b $WINDOW_SIZE -a $WINDOW_SIZE \
        -R $BED_GABA \
        -S $NESTIN_BIGWIGS \
        --samplesLabel $NESTIN_LABELS \
        --binSize $BIN_SIZE \
        --sortRegions descend \
        --sortUsing mean \
        --missingDataAsZero \
        -o $OUTPUT_DIR/matrix_GABA_nestin.gz \
        -p $N_PROCESSORS \
        --outFileNameMatrix $OUTPUT_DIR/matrix_GABA_nestin.tab \
        2>&1 | grep -v "^$"

    plotHeatmap \
        -m $OUTPUT_DIR/matrix_GABA_nestin.gz \
        -o $OUTPUT_DIR/heatmap_GABA_nestin.png \
        --colorMap Reds \
        --dpi 300 \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 \
        --refPointLabel "CRE Center" \
        --heatmapHeight 10 \
        --heatmapWidth 3 \
        --plotTitle "Nestin: ATAC Signal at Splicing ENCODE cCREs (n=$N_GABA_CRES)" \
        --xAxisLabel "Distance from CRE Center (bp)" \
        2>&1 | grep -v "^$"

    echo "Saved: heatmap_GABA_nestin.png"

    plotProfile \
        -m $OUTPUT_DIR/matrix_GABA_nestin.gz \
        -o $OUTPUT_DIR/metaprofile_GABA_nestin.png \
        --plotTitle "Nestin: ATAC Signal at Splicing ENCODE cCREs (n=$N_GABA_CRES)" \
        --refPointLabel "CRE Center" \
        --colors '#2E86AB' '#A23B72' \
        --plotHeight 6 \
        --plotWidth 8 \
        2>&1 | grep -v "^$"

    echo "Saved: metaprofile_GABA_nestin.png"
    echo ""
fi

# Emx1 only (using Nestin-Ctrl as baseline since Emx1-Ctrl is a failed sample)
if [ -f "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" ] && [ -f "$BIGWIG_BASE/GABA_Emx1-Mut.bw" ]; then
    echo "Creating Emx1 analysis (GABA CREs)..."
    echo "NOTE: Using Nestin-Ctrl as baseline (Emx1-Ctrl is a failed sample)"
    EMX1_BIGWIGS="$BIGWIG_BASE/GABA_Nestin-Ctrl.bw $BIGWIG_BASE/GABA_Emx1-Mut.bw"
    EMX1_LABELS="Nestin-Ctrl Emx1-Mut"

    computeMatrix reference-point \
        --referencePoint center \
        -b $WINDOW_SIZE -a $WINDOW_SIZE \
        -R $BED_GABA \
        -S $EMX1_BIGWIGS \
        --samplesLabel $EMX1_LABELS \
        --binSize $BIN_SIZE \
        --sortRegions descend \
        --sortUsing mean \
        --missingDataAsZero \
        -o $OUTPUT_DIR/matrix_GABA_emx1.gz \
        -p $N_PROCESSORS \
        --outFileNameMatrix $OUTPUT_DIR/matrix_GABA_emx1.tab \
        2>&1 | grep -v "^$"

    plotHeatmap \
        -m $OUTPUT_DIR/matrix_GABA_emx1.gz \
        -o $OUTPUT_DIR/heatmap_GABA_emx1.png \
        --colorMap Reds \
        --dpi 300 \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 \
        --refPointLabel "CRE Center" \
        --heatmapHeight 10 \
        --heatmapWidth 3 \
        --plotTitle "Nestin-Ctrl vs Emx1-Mut: ATAC Signal at Splicing ENCODE cCREs (n=$N_GABA_CRES)" \
        --xAxisLabel "Distance from CRE Center (bp)" \
        2>&1 | grep -v "^$"

    echo "Saved: heatmap_GABA_emx1.png"

    plotProfile \
        -m $OUTPUT_DIR/matrix_GABA_emx1.gz \
        -o $OUTPUT_DIR/metaprofile_GABA_emx1.png \
        --plotTitle "Nestin-Ctrl vs Emx1-Mut: ATAC Signal at Splicing ENCODE cCREs (n=$N_GABA_CRES)" \
        --refPointLabel "CRE Center" \
        --colors '#2E86AB' '#F18F01' \
        --plotHeight 6 \
        --plotWidth 8 \
        2>&1 | grep -v "^$"

    echo "Saved: metaprofile_GABA_emx1.png"
    echo ""
fi

# ============================================================================
# Final summary
# ============================================================================
echo "========================================================================"
echo "ANALYSIS COMPLETE!"
echo "========================================================================"
echo ""
echo "Output directory: $OUTPUT_DIR/"
echo ""
echo "CREs analyzed:"
echo "  - All CREs: $N_ALL_CRES"
echo "  - GABA CREs: $N_GABA_CRES"
echo ""
echo "Generated files:"
echo "  Heatmaps:"
echo "    - heatmap_all.png"
echo "    - heatmap_GABA.png"
echo "    - heatmap_GABA_nestin.png (Nestin-Ctrl vs Nestin-Mut)"
echo "    - heatmap_GABA_emx1.png (Nestin-Ctrl vs Emx1-Mut)"
echo "  Metaprofiles:"
echo "    - metaprofile_all.png"
echo "    - metaprofile_GABA.png"
echo "    - metaprofile_GABA_nestin.png"
echo "    - metaprofile_GABA_emx1.png"
echo ""
echo "NOTE: Emx1-Ctrl excluded (failed sample) - Nestin-Ctrl used as baseline"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
