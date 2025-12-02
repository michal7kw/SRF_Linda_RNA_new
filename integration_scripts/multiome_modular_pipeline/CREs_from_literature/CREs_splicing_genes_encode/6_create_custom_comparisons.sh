#!/bin/bash
#SBATCH --job-name=6_custom_comparisons
#SBATCH --output=logs/6_custom_comparisons.log
#SBATCH --error=logs/6_custom_comparisons.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=workq

################################################################################
# Custom ATAC Signal Comparisons for Splicing Gene ENCODE cCREs
#
# This script creates custom comparisons across genotypes and conditions:
# 1. Nestin-Ctrl vs Nestin-Mut
# 2. Nestin-Ctrl vs Emx1-Mut
# 3. Nestin-Mut vs Emx1-Mut
#
# NOTE: Emx1-Ctrl is excluded (failed sample)
#
# Prerequisites:
# - output/splicing_encode_cCREs_all.bed
# - BigWig files from Signac pipeline
#
# Output:
# - Custom comparison metaprofiles
# - Difference plots for each comparison
################################################################################

echo "========================================================================"
echo "CUSTOM ATAC SIGNAL COMPARISONS FOR SPLICING GENE ENCODE cCREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/splicing_encode_cCREs"

# Configuration
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR="./output/custom_comparisons"
BED_ALL="./output/splicing_encode_cCREs_all.bed"

# Parameters
WINDOW_SIZE=2000  # bp around CRE center
BIN_SIZE=50       # bp per bin
N_PROCESSORS=16

# Activate deepTools environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

# Check if deepTools is available
if ! command -v computeMatrix &> /dev/null; then
    echo "ERROR: deepTools not found!"
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
    echo "Please run: sbatch 2_convert_to_bed.sh"
    exit 1
fi

N_CRES=$(wc -l < $BED_ALL)
echo "Found BED file: $BED_ALL"
echo "  Total CREs: $N_CRES"
echo ""

# Create output directory
mkdir -p $OUTPUT_DIR

# Copy BED file
cp $BED_ALL $OUTPUT_DIR/splicing_encode_cCREs_all.bed
echo "Copied BED file to output directory"
echo ""

# ============================================================================
# Step 2: Check BigWig files (excluding Emx1-Ctrl)
# ============================================================================
echo "========================================================================"
echo "STEP 2: Checking BigWig files"
echo "========================================================================"
echo ""

MISSING=0

for SAMPLE in "GABA_Nestin-Ctrl" "GABA_Nestin-Mut" "GABA_Emx1-Mut"; do
    BW_FILE="$BIGWIG_BASE/${SAMPLE}.bw"
    if [ -f "$BW_FILE" ]; then
        echo "  Found: ${SAMPLE}.bw"
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
echo "NOTE: Emx1-Ctrl excluded (failed sample)"
echo ""

# ============================================================================
# Step 3: Compute matrices for custom comparisons
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

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for Nestin-Ctrl vs Nestin-Mut"
    exit 1
fi

echo "Matrix computed: Nestin-Ctrl vs Nestin-Mut"
echo ""

# Comparison 2: Nestin-Ctrl vs Emx1-Mut
echo "========================================================================"
echo "STEP 3B: Computing matrix for NESTIN-CTRL vs EMX1-MUT"
echo "========================================================================"
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

if [ $? -ne 0 ]; then
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

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for Nestin-Mut vs Emx1-Mut"
    exit 1
fi

echo "Matrix computed: Nestin-Mut vs Emx1-Mut"
echo ""

# ============================================================================
# Step 4: Create visualizations using Python
# ============================================================================
echo "========================================================================"
echo "STEP 4: Creating custom comparison visualizations"
echo "========================================================================"
echo ""

# Performance options
SKIP_FLAG="--skip-individual"
if [ "${SKIP_INDIVIDUAL}" = "0" ]; then
    echo "Full mode: Creating individual plots"
    SKIP_FLAG=""
else
    echo "Fast mode (DEFAULT): Skipping individual CRE plots"
fi

PARALLEL_FLAG=""
if [ -n "${PARALLEL_JOBS}" ]; then
    PARALLEL_FLAG="--parallel ${PARALLEL_JOBS}"
    echo "   Using ${PARALLEL_JOBS} parallel processes"
fi

DPI_FLAG=""
if [ -n "${INDIVIDUAL_DPI}" ]; then
    DPI_FLAG="--individual-dpi ${INDIVIDUAL_DPI}"
    echo "   Individual plot DPI: ${INDIVIDUAL_DPI}"
fi

FILTER_FLAG=""
if [ -n "${MIN_SIGNAL}" ]; then
    FILTER_FLAG="$FILTER_FLAG --min-signal ${MIN_SIGNAL}"
fi
if [ -n "${MIN_FC}" ]; then
    FILTER_FLAG="$FILTER_FLAG --min-fc ${MIN_FC}"
fi

echo ""
echo "Running Python visualization script..."
python 6_visualize_custom_comparisons.py $SKIP_FLAG $PARALLEL_FLAG $DPI_FLAG $FILTER_FLAG

if [ $? -ne 0 ]; then
    echo "ERROR: Visualization failed!"
    exit 1
fi

echo ""
echo "Visualizations created"
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
CUSTOM ATAC SIGNAL COMPARISONS: SPLICING GENE ENCODE cCREs
================================================================================

PURPOSE:
--------
Custom comparisons of ATAC-seq signal at ENCODE cCREs linked to splicing genes,
focusing on cross-genotype and cross-condition comparisons.

DATA SOURCE:
------------
CRE source: ENCODE cCREs (mm10-cCREs.bed)
Linkage method: Genomic proximity (+/- 500kb window)

COMPARISONS PERFORMED:
----------------------
1. Nestin-Ctrl vs Nestin-Mut
   -> Within-genotype comparison (effect of mutation in Nestin background)

2. Nestin-Ctrl vs Emx1-Mut
   -> Cross-genotype comparison (wild-type Nestin vs mutant Emx1)

3. Nestin-Mut vs Emx1-Mut
   -> Mutant-to-mutant comparison (genotype effect under mutation)

NOTE: Emx1-Ctrl excluded as it is a failed sample

OUTPUT FILES:
-------------
Metaprofiles:
  - profiles/metaprofile_nestin_ctrl_vs_nestin_mut.png
  - profiles/metaprofile_nestin_ctrl_vs_emx1_mut.png
  - profiles/metaprofile_nestin_mut_vs_emx1_mut.png

Signal matrices:
  - matrix_nestin_ctrl_vs_mut.gz / .tab
  - matrix_nestin_ctrl_vs_emx1_mut.gz / .tab
  - matrix_nestin_mut_vs_emx1_mut.gz / .tab

INTERPRETATION:
---------------
Comparison 1 (Nestin-Ctrl vs Nestin-Mut):
  - Shows effect of mutation within Nestin genotype

Comparison 2 (Nestin-Ctrl vs Emx1-Mut):
  - Shows cross-genotype effect (wild-type vs mutant)

Comparison 3 (Nestin-Mut vs Emx1-Mut):
  - Shows genotype effect under mutation

Generated by: 6_create_custom_comparisons.sh
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
echo "CREs analyzed: $N_CRES ENCODE cCREs"
echo ""
echo "Generated files:"
echo "  Metaprofiles (3 custom comparisons):"
echo "      - profiles/metaprofile_nestin_ctrl_vs_nestin_mut.png"
echo "      - profiles/metaprofile_nestin_ctrl_vs_emx1_mut.png"
echo "      - profiles/metaprofile_nestin_mut_vs_emx1_mut.png"
echo ""
echo "COMPARISONS:"
echo "  1. Nestin-Ctrl vs Nestin-Mut  (within-genotype mutation effect)"
echo "  2. Nestin-Ctrl vs Emx1-Mut    (cross-genotype comparison)"
echo "  3. Nestin-Mut vs Emx1-Mut     (mutant-to-mutant genotype effect)"
echo ""
echo "NOTE: Emx1-Ctrl excluded (failed sample)"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
