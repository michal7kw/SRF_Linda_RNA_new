#!/bin/bash
#SBATCH --job-name=5_create_custom_comparisons
#SBATCH --output=logs/5_create_custom_comparisons.log
#SBATCH --error=logs/5_create_custom_comparisons.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Custom ATAC Signal Comparisons for ENCODE cCREs
#
# IMPORTANT: Emx1-Ctrl is a FAILED SAMPLE and is excluded from all analyses.
# Only Nestin-Ctrl is used as the reference control for both mutant samples.
#
# This script creates custom comparisons (3 total):
# 1. Nestin-Ctrl vs Nestin-Mut (within-genotype mutation effect)
# 2. Nestin-Ctrl vs Emx1-Mut (cross-genotype, using Nestin-Ctrl as reference)
# 3. Nestin-Mut vs Emx1-Mut (mutant comparison)
#
# Prerequisites:
# - output/encode_cCREs_all.bed
# - BigWig files from Signac pipeline (signac_results_L1/bigwig_tracks_L1/)
#
# Output:
# - Custom comparison matrices for visualization
# - Metaprofiles and difference plots
# - Statistical summaries
################################################################################

echo "========================================================================"
echo "CUSTOM ATAC SIGNAL COMPARISONS FOR ENCODE cCREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""
echo "NOTE: Emx1-Ctrl is EXCLUDED (failed sample)"
echo "      Using Nestin-Ctrl as reference for all comparisons"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_encode_all"

# Configuration
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR="./output/custom_comparisons"
BED_ALL="./output/encode_cCREs_all.bed"

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
# Step 1: Check BED file exists
# ============================================================================
echo "========================================================================"
echo "STEP 1: Checking BED file"
echo "========================================================================"
echo ""

if [ ! -f "$BED_ALL" ]; then
    echo "ERROR: BED file not found: $BED_ALL"
    echo "Please run: sbatch 2_convert_encode_cCREs_to_bed.sh"
    exit 1
fi

N_CRES=$(wc -l < $BED_ALL)
echo "Found BED file: $BED_ALL"
echo "  Total CREs: $N_CRES"
echo ""

# Create output directory
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/profiles

# Copy BED file to output directory for reference
cp $BED_ALL $OUTPUT_DIR/encode_cCREs_all.bed
echo "Copied BED file to output directory"
echo ""

# ============================================================================
# Step 2: Check BigWig files (excluding Emx1-Ctrl - failed sample)
# ============================================================================
echo "========================================================================"
echo "STEP 2: Checking BigWig files"
echo "========================================================================"
echo ""

MISSING=0

# Only check for 3 valid samples (Emx1-Ctrl is excluded)
for SAMPLE in "GABA_Nestin-Ctrl" "GABA_Nestin-Mut" "GABA_Emx1-Mut"; do
    BW_FILE="$BIGWIG_BASE/${SAMPLE}.bw"
    if [ -f "$BW_FILE" ]; then
        echo "  Found: ${SAMPLE}.bw"
    else
        echo "  Missing: $BW_FILE"
        MISSING=$((MISSING + 1))
    fi
done

echo "  EXCLUDED: GABA_Emx1-Ctrl.bw (failed sample)"

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "ERROR: Missing $MISSING BigWig files!"
    exit 1
fi

echo ""

# ============================================================================
# Step 3: Compute matrices for custom comparisons (3 comparisons)
# ============================================================================

# Comparison 1: Nestin-Ctrl vs Nestin-Mut
echo "========================================================================"
echo "STEP 3A: Computing matrix for NESTIN-CTRL vs NESTIN-MUT"
echo "========================================================================"
echo ""

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_ALL \
    -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Nestin-Mut.bw" \
    --samplesLabel "Nestin-Ctrl" "Nestin-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_nestin_ctrl_vs_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_nestin_ctrl_vs_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_nestin_ctrl_vs_mut.log

if [ ! -f "$OUTPUT_DIR/matrix_nestin_ctrl_vs_mut.gz" ]; then
    echo "ERROR: computeMatrix failed for Nestin-Ctrl vs Nestin-Mut"
    exit 1
fi

echo "Matrix computed: Nestin-Ctrl vs Nestin-Mut"
echo ""

# Comparison 2: Nestin-Ctrl vs Emx1-Mut (cross-genotype, Nestin-Ctrl as reference)
echo "========================================================================"
echo "STEP 3B: Computing matrix for NESTIN-CTRL vs EMX1-MUT"
echo "========================================================================"
echo "(Using Nestin-Ctrl as reference control for Emx1-Mut)"
echo ""

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_ALL \
    -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
    --samplesLabel "Nestin-Ctrl" "Emx1-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_nestin_ctrl_vs_emx1_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_nestin_ctrl_vs_emx1_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_nestin_ctrl_vs_emx1_mut.log

if [ ! -f "$OUTPUT_DIR/matrix_nestin_ctrl_vs_emx1_mut.gz" ]; then
    echo "ERROR: computeMatrix failed for Nestin-Ctrl vs Emx1-Mut"
    exit 1
fi

echo "Matrix computed: Nestin-Ctrl vs Emx1-Mut"
echo ""

# Comparison 3: Nestin-Mut vs Emx1-Mut
echo "========================================================================"
echo "STEP 3C: Computing matrix for NESTIN-MUT vs EMX1-MUT"
echo "========================================================================"
echo ""

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_ALL \
    -S "$BIGWIG_BASE/GABA_Nestin-Mut.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
    --samplesLabel "Nestin-Mut" "Emx1-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_nestin_mut_vs_emx1_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_nestin_mut_vs_emx1_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_nestin_mut_vs_emx1_mut.log

if [ ! -f "$OUTPUT_DIR/matrix_nestin_mut_vs_emx1_mut.gz" ]; then
    echo "ERROR: computeMatrix failed for Nestin-Mut vs Emx1-Mut"
    exit 1
fi

echo "Matrix computed: Nestin-Mut vs Emx1-Mut"
echo ""

# ============================================================================
# Step 4: Create visualizations using deepTools plotProfile
# ============================================================================
echo "========================================================================"
echo "STEP 4: Creating metaprofile visualizations"
echo "========================================================================"
echo ""

# Create metaprofiles for each comparison using plotProfile
echo "Creating metaprofile: Nestin-Ctrl vs Nestin-Mut"
plotProfile -m $OUTPUT_DIR/matrix_nestin_ctrl_vs_mut.gz \
    --perGroup \
    --colors "#2E86AB" "#A23B72" \
    --plotTitle "Nestin: Ctrl vs Mut - ATAC Signal at ENCODE cCREs" \
    --yAxisLabel "Mean ATAC Signal" \
    --regionsLabel "CREs" \
    --plotHeight 6 --plotWidth 8 \
    -o $OUTPUT_DIR/profiles/metaprofile_nestin_ctrl_vs_mut.png \
    --dpi 300

echo "Creating metaprofile: Nestin-Ctrl vs Emx1-Mut"
plotProfile -m $OUTPUT_DIR/matrix_nestin_ctrl_vs_emx1_mut.gz \
    --perGroup \
    --colors "#2E86AB" "#C73E1D" \
    --plotTitle "Cross-genotype: Nestin-Ctrl vs Emx1-Mut - ATAC Signal at ENCODE cCREs" \
    --yAxisLabel "Mean ATAC Signal" \
    --regionsLabel "CREs" \
    --plotHeight 6 --plotWidth 8 \
    -o $OUTPUT_DIR/profiles/metaprofile_nestin_ctrl_vs_emx1_mut.png \
    --dpi 300

echo "Creating metaprofile: Nestin-Mut vs Emx1-Mut"
plotProfile -m $OUTPUT_DIR/matrix_nestin_mut_vs_emx1_mut.gz \
    --perGroup \
    --colors "#A23B72" "#C73E1D" \
    --plotTitle "Mutants: Nestin-Mut vs Emx1-Mut - ATAC Signal at ENCODE cCREs" \
    --yAxisLabel "Mean ATAC Signal" \
    --regionsLabel "CREs" \
    --plotHeight 6 --plotWidth 8 \
    -o $OUTPUT_DIR/profiles/metaprofile_nestin_vs_emx1_mut.png \
    --dpi 300

echo ""
echo "Metaprofiles created"
echo ""

# ============================================================================
# Step 4b: Create combined 3-sample matrix and overview (excludes Emx1-Ctrl)
# ============================================================================
echo "========================================================================"
echo "STEP 4B: Creating combined overview (3 conditions, excludes Emx1-Ctrl)"
echo "========================================================================"
echo ""

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_ALL \
    -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Nestin-Mut.bw" \
       "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
    --samplesLabel "Nestin-Ctrl" "Nestin-Mut" "Emx1-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_all_conditions.gz \
    -p $N_PROCESSORS \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_all_conditions.log

plotProfile -m $OUTPUT_DIR/matrix_all_conditions.gz \
    --perGroup \
    --colors "#2E86AB" "#A23B72" "#C73E1D" \
    --plotTitle "All Conditions - ATAC Signal at ENCODE cCREs ($N_CRES regions)" \
    --yAxisLabel "Mean ATAC Signal" \
    --regionsLabel "CREs" \
    --plotHeight 6 --plotWidth 10 \
    -o $OUTPUT_DIR/profiles/overview_all_conditions.png \
    --dpi 300

echo ""
echo "Overview created"
echo ""

# ============================================================================
# Step 5: Create summary report
# ============================================================================
echo "========================================================================"
echo "STEP 5: Creating summary report"
echo "========================================================================"
echo ""

cat > $OUTPUT_DIR/README.txt << 'EOFR'
================================================================================
CUSTOM ATAC SIGNAL COMPARISONS: ENCODE cCREs
================================================================================

IMPORTANT: Emx1-Ctrl is a FAILED SAMPLE and is excluded from all analyses.
           Nestin-Ctrl is used as the reference control for both mutant samples.

PURPOSE:
--------
Custom comparisons of ATAC-seq signal at all ENCODE cCREs,
using Nestin-Ctrl as the reference control for both mutants.

COMPARISONS PERFORMED (3 total):
--------------------------------
1. Nestin-Ctrl vs Nestin-Mut
   -> Within-genotype comparison (effect of mutation in Nestin background)

2. Nestin-Ctrl vs Emx1-Mut
   -> Cross-genotype comparison (using Nestin-Ctrl as reference for Emx1-Mut)

3. Nestin-Mut vs Emx1-Mut
   -> Mutant genotype comparison (comparing mutation effects across genotypes)

PARAMETERS:
-----------
Window size: +/-2000 bp around CRE center
Bin size: 50 bp
Reference point: CRE center
Processors: 16 (parallel processing)

OUTPUT FILES:
-------------
Metaprofiles (comparison overlays):
  - profiles/metaprofile_nestin_ctrl_vs_mut.png
  - profiles/metaprofile_nestin_ctrl_vs_emx1_mut.png
  - profiles/metaprofile_nestin_vs_emx1_mut.png

Overview:
  - profiles/overview_all_conditions.png (3 conditions)

Signal matrices (reusable):
  - matrix_nestin_ctrl_vs_mut.gz / .tab
  - matrix_nestin_ctrl_vs_emx1_mut.gz / .tab
  - matrix_nestin_mut_vs_emx1_mut.gz / .tab
  - matrix_all_conditions.gz (combined 3-sample matrix)

INTERPRETATION GUIDE:
---------------------
Comparison 1 (Nestin-Ctrl vs Nestin-Mut):
  - Shows effect of mutation within Nestin genotype
  - Higher signal in Mut = Increased accessibility due to mutation

Comparison 2 (Nestin-Ctrl vs Emx1-Mut):
  - Shows combined effect of genotype + mutation
  - Since Emx1-Ctrl failed, we use Nestin-Ctrl as the baseline
  - Differences reflect both genotype and mutation effects

Comparison 3 (Nestin-Mut vs Emx1-Mut):
  - Shows genotype effect under mutation
  - Differences indicate genotype-specific responses to mutation

KEY QUESTIONS TO ADDRESS:
--------------------------
1. Does mutation affect chromatin accessibility in Nestin?
   -> Check Comparison 1

2. How does Emx1-Mut compare to a control baseline?
   -> Check Comparison 2 (using Nestin-Ctrl as reference)

3. Do mutants differ between genotypes?
   -> Check Comparison 3

Generated by: 5_create_custom_comparisons.sh
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
echo "CREs analyzed: $N_CRES ENCODE cCREs (all types)"
echo ""
echo "NOTE: Emx1-Ctrl excluded (failed sample)"
echo ""
echo "Generated files:"
echo ""
echo "  Matrices (4 files):"
echo "      - matrix_nestin_ctrl_vs_mut.gz"
echo "      - matrix_nestin_ctrl_vs_emx1_mut.gz"
echo "      - matrix_nestin_mut_vs_emx1_mut.gz"
echo "      - matrix_all_conditions.gz (3 samples)"
echo ""
echo "  Metaprofiles (4 plots):"
ls -1 $OUTPUT_DIR/profiles/*.png 2>/dev/null | while read f; do
    echo "      - $(basename $f)"
done
echo ""
echo "  Documentation:"
echo "      - README.txt (interpretation guide)"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
