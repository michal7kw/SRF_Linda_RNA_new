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
# This script creates custom comparisons across genotypes and conditions:
# 1. Nestin-Ctrl vs Nestin-Mut (within-genotype mutation effect)
# 2. Nestin-Ctrl vs Emx1-Mut (cross-genotype, using Nestin-Ctrl as control)
# 3. Nestin-Mut vs Emx1-Mut (mutant comparison)
#
# NOTE: Emx1-Ctrl is EXCLUDED (failed sample) - using Nestin-Ctrl as control
#
# Prerequisites:
# - output/encode_cCREs_GABA.bed
# - BigWig files from Signac pipeline (signac_results_L1/bigwig_tracks_L1/)
#
# Output:
# - Custom comparison matrices for visualization
# - Metaprofiles and difference plots
# - Statistical summaries
#
# PERFORMANCE OPTIONS:
# -----------------------------------------
# DEFAULT: Quick mode (metaprofiles only, skips individual plots):
#   sbatch 5_create_custom_comparisons.sh
#
# Full mode (create individual plots):
#   SKIP_INDIVIDUAL=0 sbatch 5_create_custom_comparisons.sh
#
# Full mode with parallel processing (8x faster):
#   SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 5_create_custom_comparisons.sh
################################################################################

echo "========================================================================"
echo "CUSTOM ATAC SIGNAL COMPARISONS FOR ENCODE cCREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Change to script directory (use absolute path for SLURM compatibility)
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_encode_paper_intersection"
cd "$SCRIPT_DIR"

# Configuration
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR="./output/custom_comparisons"
BED_GABA="./output/encode_cCREs_GABA.bed"

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

if [ ! -f "$BED_GABA" ]; then
    echo "ERROR: BED file not found: $BED_GABA"
    echo "Please run: sbatch 2_convert_encode_cCREs_to_bed.sh"
    exit 1
fi

N_CRES=$(wc -l < $BED_GABA)
echo "Found BED file: $BED_GABA"
echo "  Total CREs: $N_CRES"
echo ""

# Create output directory
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/profiles

# Copy BED file to output directory for reference
cp $BED_GABA $OUTPUT_DIR/encode_cCREs_GABA.bed
echo "Copied BED file to output directory"
echo ""

# ============================================================================
# Step 2: Check BigWig files
# ============================================================================
echo "========================================================================"
echo "STEP 2: Checking BigWig files"
echo "========================================================================"
echo ""

MISSING=0

# NOTE: Emx1-Ctrl is excluded (failed sample) - using only Nestin-Ctrl as control
for SAMPLE in "GABA_Nestin-Ctrl" "GABA_Nestin-Mut" "GABA_Emx1-Mut"; do
    BW_FILE="$BIGWIG_BASE/${SAMPLE}.bw"
    if [ -f "$BW_FILE" ]; then
        echo "  Found: ${SAMPLE}.bw"
    else
        echo "  Missing: $BW_FILE"
        MISSING=$((MISSING + 1))
    fi
done

echo "  Skipped: GABA_Emx1-Ctrl.bw (failed sample - excluded from analysis)"

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "ERROR: Missing $MISSING BigWig files!"
    exit 1
fi

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
    -R $BED_GABA \
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

# Comparison 2: Nestin-Ctrl vs Emx1-Mut (using Nestin-Ctrl as control since Emx1-Ctrl failed)
echo "========================================================================"
echo "STEP 3B: Computing matrix for NESTIN-CTRL vs EMX1-MUT"
echo "========================================================================"
echo ""

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_GABA \
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
    -R $BED_GABA \
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
# Step 4: Create visualizations using Python
# ============================================================================
echo "========================================================================"
echo "STEP 4: Creating custom comparison visualizations"
echo "========================================================================"
echo ""

# Performance options (can be set via environment variables):
# - SKIP_INDIVIDUAL=0: Create individual CRE plots (default: skip to save time)
# - PARALLEL_JOBS=N: Use N parallel processes for individual plots
# - INDIVIDUAL_DPI=N: DPI for individual plots (default: 150)

SKIP_FLAG="--skip-individual"
if [ "${SKIP_INDIVIDUAL}" = "0" ]; then
    echo "Full mode: Creating individual plots for each comparison"
    SKIP_FLAG=""
else
    echo "Fast mode (DEFAULT): Skipping individual CRE plots"
    echo "   (To create individual plots, set SKIP_INDIVIDUAL=0)"
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

echo ""
echo "Running Python visualization script..."
python 5_visualize_custom_comparisons.py $SKIP_FLAG $PARALLEL_FLAG $DPI_FLAG

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
CUSTOM ATAC SIGNAL COMPARISONS: ENCODE cCREs
================================================================================

PURPOSE:
--------
Custom comparisons of ATAC-seq signal at ENCODE cCREs linked to splicing genes,
focusing on within-genotype and cross-genotype comparisons.

COMPARISONS PERFORMED:
----------------------
NOTE: Emx1-Ctrl is EXCLUDED (failed sample) - using Nestin-Ctrl as control

1. Nestin-Ctrl vs Nestin-Mut
   -> Within-genotype comparison (effect of mutation in Nestin background)

2. Nestin-Ctrl vs Emx1-Mut
   -> Cross-genotype comparison (using Nestin-Ctrl as control for Emx1-Mut)

3. Nestin-Mut vs Emx1-Mut
   -> Mutant genotype comparison (mutation effect across genotypes)

PARAMETERS:
-----------
Window size: +/-2000 bp around CRE center
Bin size: 50 bp
Reference point: CRE center
Processors: 16 (parallel processing)

OUTPUT FILES:
-------------
Metaprofiles (comparison overlays):
  - profiles/metaprofile_nestin_ctrl_vs_mut_GABA.png
  - profiles/metaprofile_nestin_ctrl_vs_emx1_mut_GABA.png
  - profiles/metaprofile_nestin_vs_emx1_mut_GABA.png

Overview:
  - overview_all_conditions_GABA.png (3 conditions: Nestin-Ctrl, Nestin-Mut, Emx1-Mut)

Statistics:
  - comparison_statistics_GABA.txt

Signal matrices (reusable):
  - matrix_nestin_ctrl_vs_mut.gz / .tab
  - matrix_nestin_ctrl_vs_emx1_mut.gz / .tab
  - matrix_nestin_mut_vs_emx1_mut.gz / .tab

INTERPRETATION GUIDE:
---------------------
Comparison 1 (Nestin-Ctrl vs Nestin-Mut):
  - Shows effect of mutation within Nestin genotype
  - Higher signal in Mut = Increased accessibility due to mutation

Comparison 2 (Nestin-Ctrl vs Emx1-Mut):
  - Uses Nestin-Ctrl as control for Emx1-Mut (since Emx1-Ctrl failed)
  - Shows combined effect of genotype + mutation

Comparison 3 (Nestin-Mut vs Emx1-Mut):
  - Shows genotype effect under mutation
  - Differences indicate genotype-specific mutation responses

KEY QUESTIONS TO ADDRESS:
--------------------------
1. Does mutation affect chromatin accessibility in Nestin?
   -> Check Comparison 1

2. How does Emx1-Mut compare to baseline Nestin-Ctrl?
   -> Check Comparison 2

3. Do Nestin-Mut and Emx1-Mut show similar patterns?
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
echo "CREs analyzed: $N_CRES GABA-specific ENCODE cCREs"
echo ""
echo "Generated files:"
echo ""
echo "  Matrices (3 comparisons - Emx1-Ctrl excluded):"
echo "      - matrix_nestin_ctrl_vs_mut.gz"
echo "      - matrix_nestin_ctrl_vs_emx1_mut.gz"
echo "      - matrix_nestin_mut_vs_emx1_mut.gz"
echo ""
echo "  Metaprofiles (3 comparison overlays):"
ls -1 $OUTPUT_DIR/profiles/metaprofile_*.png 2>/dev/null | while read f; do
    echo "      - $(basename $f)"
done
echo ""
echo "  Documentation:"
echo "      - README.txt (interpretation guide)"
echo "      - comparison_statistics_GABA.txt"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
