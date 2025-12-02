#!/bin/bash
#SBATCH --job-name=4_create_heatmaps_splicing_genes
#SBATCH --output=logs/4_create_heatmaps_splicing_genes.log
#SBATCH --error=logs/4_create_heatmaps_splicing_genes.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=workq

################################################################################
# Create Heatmaps and Metaprofiles for Splicing Gene CREs
#
# This script uses deepTools to create heatmaps and metaprofiles for CREs
# associated with splicing-related genes from Reactome pathways.
#
# Prerequisites:
# - output/splicing_genes_analysis/splicing_genes_CREs_all.bed
# - output/splicing_genes_analysis/splicing_genes_CREs_GABA.bed
# - BigWig files from Signac pipeline (signac_results_L1/bigwig_tracks_L1/)
#
# Output:
# - Heatmaps showing ATAC signal at splicing gene CREs
# - Metaprofiles showing average signal profiles
# - Genotype-specific analyses (Nestin/Emx1 separate)
################################################################################

echo "========================================================================"
echo "CREATE HEATMAPS FOR SPLICING GENE CREs USING DEEPTOOLS"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_genes_paper"

# Configuration
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR="./output/heatmaps_deeptools"
BED_ALL="./output/splicing_genes_CREs_all.bed"
BED_GABA="./output/splicing_genes_CREs_GABA.bed"

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

# ============================================================================
# Step 1: Check BED files exist
# ============================================================================
echo "========================================================================"
echo "STEP 1: Checking BED files"
echo "========================================================================"
echo ""

# Check all cell types BED file
if [ ! -f "$BED_ALL" ]; then
    echo "ERROR: All cell types BED file not found: $BED_ALL"
    echo "Please run: sbatch run_convert_splicing_CREs_to_bed.sh"
    exit 1
fi

N_ALL_CRES=$(wc -l < $BED_ALL)
echo "✓ Found all cell types BED file: $BED_ALL"
echo "  Total CREs (all cell types): $N_ALL_CRES"
echo ""

# Check GABA BED file
if [ ! -f "$BED_GABA" ]; then
    echo "ERROR: GABA BED file not found: $BED_GABA"
    echo "Please run: sbatch run_convert_splicing_CREs_to_bed.sh"
    exit 1
fi

N_GABA_CRES=$(wc -l < $BED_GABA)
echo "✓ Found GABA BED file: $BED_GABA"
echo "  Total CREs (GABA cell types): $N_GABA_CRES"
echo ""

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Copy BED files to output directory for reference
cp $BED_ALL $OUTPUT_DIR/splicing_genes_CREs_all.bed
cp $BED_GABA $OUTPUT_DIR/splicing_genes_CREs_GABA.bed
echo "✓ Copied BED files to output directory"
echo ""

# Show CRE IDs
echo "CRE IDs in analysis:"
awk '{print "  - " $4 " (" $1 ":" $2 "-" $3 ")"}' $BED_ALL
echo ""

# ============================================================================
# Step 2: Check BigWig files
# ============================================================================
echo "========================================================================"
echo "STEP 2: Checking BigWig files"
echo "========================================================================"
echo ""

# Collect BigWig files (excluding Emx1-Ctrl which is a failed sample)
BIGWIGS=""
LABELS=""
MISSING=0

# Only check for: Nestin-Ctrl, Nestin-Mut, Emx1-Mut (NOT Emx1-Ctrl)
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

echo ""
echo "NOTE: Emx1-Ctrl excluded (failed sample)"

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "ERROR: Missing $MISSING BigWig files!"
    exit 1
fi

echo ""

# ============================================================================
# Step 3: Run computeMatrix for ALL cell types
# ============================================================================
echo "========================================================================"
echo "STEP 3: Computing signal matrices for ALL CELL TYPES"
echo "========================================================================"
echo ""

echo "Running computeMatrix on all cell types CREs..."
echo "  Window size: ±${WINDOW_SIZE} bp"
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
    -o $OUTPUT_DIR/matrix_all_celltypes.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_all_celltypes.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_all_celltypes.log

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for all cell types"
    exit 1
fi

echo ""
echo "✓ All cell types matrix computed and saved"
echo ""

# ============================================================================
# Step 4: Run computeMatrix for GABA cell types
# ============================================================================
echo "========================================================================"
echo "STEP 4: Computing signal matrices for GABA CELL TYPES"
echo "========================================================================"
echo ""

echo "Running computeMatrix on GABA cell types CREs..."
echo "  Window size: ±${WINDOW_SIZE} bp"
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
    echo "ERROR: computeMatrix failed for GABA cell types"
    exit 1
fi

echo ""
echo "✓ GABA cell types matrix computed and saved"
echo ""

# ============================================================================
# Step 5: Create heatmaps
# ============================================================================
echo "========================================================================"
echo "STEP 5: Creating heatmaps"
echo "========================================================================"
echo ""

# All cell types heatmap
echo "Creating all cell types heatmap..."
plotHeatmap \
    -m $OUTPUT_DIR/matrix_all_celltypes.gz \
    -o $OUTPUT_DIR/heatmap_all_celltypes.png \
    --colorMap Reds \
    --dpi 300 \
    --whatToShow 'heatmap and colorbar' \
    --zMin 0 \
    --refPointLabel "CRE Center" \
    --heatmapHeight 10 \
    --heatmapWidth 6 \
    --sortRegions descend \
    --sortUsing mean \
    --plotTitle "ATAC Signal at Splicing Gene CREs - All Cell Types (n=$N_ALL_CRES)" \
    --xAxisLabel "Distance from CRE Center (bp)" \
    --yAxisLabel "CREs" \
    2>&1 | grep -v "^$"

echo "✓ Saved: heatmap_all_celltypes.png"

# GABA cell types heatmap
echo "Creating GABA cell types heatmap..."
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
    --plotTitle "ATAC Signal at Splicing Gene CREs - GABA Cell Types (n=$N_GABA_CRES)" \
    --xAxisLabel "Distance from CRE Center (bp)" \
    --yAxisLabel "CREs" \
    2>&1 | grep -v "^$"

echo "✓ Saved: heatmap_GABA.png"
echo ""

# ============================================================================
# Step 6: Create metaprofiles
# ============================================================================
echo "========================================================================"
echo "STEP 6: Creating metaprofiles"
echo "========================================================================"
echo ""

# All cell types metaprofile
echo "Creating all cell types metaprofile..."
plotProfile \
    -m $OUTPUT_DIR/matrix_all_celltypes.gz \
    -o $OUTPUT_DIR/metaprofile_all_celltypes.png \
    --plotTitle "ATAC Signal at Splicing Gene CREs - All Cell Types (n=$N_ALL_CRES)" \
    --refPointLabel "CRE Center" \
    --averageType mean \
    --plotHeight 7 \
    --plotWidth 10 \
    --colors '#2E86AB' '#A23B72' '#F18F01' '#C73E1D' \
    --yAxisLabel "Mean ATAC Signal" \
    2>&1 | grep -v "^$"

echo "✓ Saved: metaprofile_all_celltypes.png"

# GABA cell types metaprofile
echo "Creating GABA cell types metaprofile..."
plotProfile \
    -m $OUTPUT_DIR/matrix_GABA.gz \
    -o $OUTPUT_DIR/metaprofile_GABA.png \
    --plotTitle "ATAC Signal at Splicing Gene CREs - GABA Cell Types (n=$N_GABA_CRES)" \
    --refPointLabel "CRE Center" \
    --averageType mean \
    --plotHeight 7 \
    --plotWidth 10 \
    --colors '#2E86AB' '#A23B72' '#F18F01' '#C73E1D' \
    --yAxisLabel "Mean ATAC Signal" \
    2>&1 | grep -v "^$"

echo "✓ Saved: metaprofile_GABA.png"
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
        --plotTitle "Nestin: ATAC Signal at Splicing Gene CREs (n=$N_GABA_CRES)" \
        --xAxisLabel "Distance from CRE Center (bp)" \
        2>&1 | grep -v "^$"

    echo "✓ Saved: heatmap_GABA_nestin.png"

    plotProfile \
        -m $OUTPUT_DIR/matrix_GABA_nestin.gz \
        -o $OUTPUT_DIR/metaprofile_GABA_nestin.png \
        --plotTitle "Nestin: ATAC Signal at Splicing Gene CREs (n=$N_GABA_CRES)" \
        --refPointLabel "CRE Center" \
        --colors '#2E86AB' '#A23B72' \
        --plotHeight 6 \
        --plotWidth 8 \
        2>&1 | grep -v "^$"

    echo "✓ Saved: metaprofile_GABA_nestin.png"
    echo ""
fi

# Emx1 analysis (using Nestin-Ctrl as control since Emx1-Ctrl is a failed sample)
if [ -f "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" ] && [ -f "$BIGWIG_BASE/GABA_Emx1-Mut.bw" ]; then
    echo "Creating Emx1 analysis (GABA CREs)..."
    echo "NOTE: Using Nestin-Ctrl as control (Emx1-Ctrl is a failed sample)"
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
        --plotTitle "Emx1: ATAC Signal at Splicing Gene CREs (n=$N_GABA_CRES)" \
        --xAxisLabel "Distance from CRE Center (bp)" \
        2>&1 | grep -v "^$"

    echo "✓ Saved: heatmap_GABA_emx1.png"

    plotProfile \
        -m $OUTPUT_DIR/matrix_GABA_emx1.gz \
        -o $OUTPUT_DIR/metaprofile_GABA_emx1.png \
        --plotTitle "Emx1: ATAC Signal at Splicing Gene CREs (n=$N_GABA_CRES)" \
        --refPointLabel "CRE Center" \
        --colors '#F18F01' '#C73E1D' \
        --plotHeight 6 \
        --plotWidth 8 \
        2>&1 | grep -v "^$"

    echo "✓ Saved: metaprofile_GABA_emx1.png"
    echo ""
fi

# ============================================================================
# Step 8: Create summary report
# ============================================================================
echo "========================================================================"
echo "STEP 8: Creating summary report"
echo "========================================================================"
echo ""

# Get list of genes
GENE_LIST=$(awk -F'\t' 'NR>1 {print $1}' ./output/splicing_genes_CREs_all_celltypes.tsv | sort -u | tr '\n' ', ' | sed 's/,$//')

cat > $OUTPUT_DIR/README.txt << EOFR
================================================================================
DEEPTOOLS HEATMAPS: CREs ASSOCIATED WITH SPLICING GENES
================================================================================

Analysis date: $(date)

PURPOSE:
--------
Visualize ATAC-seq signal at CREs linked to splicing-related genes from
Reactome pathways to understand chromatin accessibility patterns at regulatory
elements controlling splicing machinery.

SPLICING GENES ANALYZED:
-------------------------
Total genes: 5
Gene list: $GENE_LIST

These genes are part of:
- mRNA splicing pathways
- Spliceosome assembly
- Pre-mRNA processing

SOURCE DATA:
-----------
Splicing genes list:
  /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/CREs_splicing_genes_paper/extracted_genes_final.csv
  Total: Variable number of genes from Reactome/GO pathways

CRE-gene links:
  /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/helpers/output/CRE_gene_links_detailed.tsv
  Created by: Signac peak-gene linkage analysis

BED files (converted from TSV):
  output/splicing_genes_analysis/splicing_genes_CREs_all.bed ($N_ALL_CRES CREs)
  output/splicing_genes_analysis/splicing_genes_CREs_GABA.bed ($N_GABA_CRES CREs)

BigWig files:
  ../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw
  Created by: signac_06_create_bigwig_tracks.R

CREs ANALYZED:
--------------
All cell types: $N_ALL_CRES CREs
GABA cell types: $N_GABA_CRES CREs
GABA-specific: 0 CREs (all GABA CREs are shared with other cell types)

Coverage: 5 of 293 splicing genes (1.7%) have associated CREs
(Low coverage due to subset of links in CRE_gene_links_detailed.tsv)

PARAMETERS:
-----------
Window size: ±${WINDOW_SIZE} bp around CRE center
Bin size: ${BIN_SIZE} bp
Reference point: CRE center
Sorting: By mean signal (descending)
Processors: $N_PROCESSORS (parallel processing)

CONDITIONS:
-----------
- Nestin-Ctrl
- Nestin-Mut
- Emx1-Mut
NOTE: Emx1-Ctrl excluded (failed sample); Nestin-Ctrl used as control for Emx1

OUTPUT FILES:
-------------
All cell types:
  - heatmap_all_celltypes.png
  - metaprofile_all_celltypes.png
  - matrix_all_celltypes.gz (reusable compressed matrix)

GABA cell types:
  - heatmap_GABA.png
  - metaprofile_GABA.png
  - matrix_GABA.gz (reusable compressed matrix)

Genotype-specific (GABA CREs):
  - heatmap_GABA_nestin.png / metaprofile_GABA_nestin.png
  - heatmap_GABA_emx1.png / metaprofile_GABA_emx1.png

Logs and metadata:
  - computeMatrix_all_celltypes.log
  - computeMatrix_GABA.log
  - README.txt (this file)

INTERPRETATION:
---------------
Heatmaps:
  - Red color indicates higher ATAC-seq signal (more accessible chromatin)
  - Signal should be concentrated at CRE centers
  - Sorted by mean signal (highest signal at top)

Metaprofiles:
  - Show average ATAC signal across all CREs
  - Peak at center indicates CRE-centered accessibility
  - Compare Ctrl vs Mut to identify condition-specific changes
  - Compare Nestin vs Emx1 to identify genotype-specific patterns

KEY FINDINGS:
-------------
1. All 5 CREs linked to splicing genes are active in GABA neurons
2. None are GABA-specific (all shared with other cell types)
3. This suggests splicing machinery is widely expressed across cell types
   (expected, as splicing is a fundamental cellular process)

BIOLOGICAL CONTEXT:
-------------------
Splicing genes analyzed include:
- Fus: RNA-binding protein involved in splicing
- Rbm5: RNA-binding motif protein 5 (splicing regulator)
- Srsf1, Srsf2, Srsf5: Serine/arginine-rich splicing factors

These genes are essential for:
- mRNA splicing and processing
- Alternative splicing regulation
- Spliceosome assembly and function

NEXT STEPS:
-----------
1. Compare heatmaps across conditions to identify accessibility changes
2. Correlate CRE accessibility with gene expression
3. Investigate whether SRF mutations affect splicing gene CRE accessibility
4. Expand analysis to full CRE-gene link dataset if needed

Generated by: $(basename $0)
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
echo "CREs analyzed:"
echo "  - All cell types: $N_ALL_CRES CREs"
echo "  - GABA cell types: $N_GABA_CRES CREs"
echo "  - GABA-specific: 0 CREs (skipped)"
echo ""
echo "Genes: 5 splicing genes (Fus, Rbm5, Srsf1, Srsf2, Srsf5)"
echo ""
echo "Generated files:"
echo "  Heatmaps:"
echo "    - heatmap_all_celltypes.png"
echo "    - heatmap_GABA.png"
echo "    - heatmap_GABA_nestin.png"
echo "    - heatmap_GABA_emx1.png"
echo "  Metaprofiles:"
echo "    - metaprofile_all_celltypes.png"
echo "    - metaprofile_GABA.png"
echo "    - metaprofile_GABA_nestin.png"
echo "    - metaprofile_GABA_emx1.png"
echo ""
echo "Review the heatmaps to identify ATAC signal patterns at"
echo "splicing gene CREs across conditions and genotypes."
echo ""
echo "Completed: $(date)"
echo "========================================================================"
