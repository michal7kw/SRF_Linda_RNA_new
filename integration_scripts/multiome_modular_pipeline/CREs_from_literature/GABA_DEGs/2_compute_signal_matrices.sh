#!/bin/bash
#SBATCH --job-name=2_compute_deg_matrices
#SBATCH --output=logs/2_compute_signal_matrices.log
#SBATCH --error=logs/2_compute_signal_matrices.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Compute ATAC Signal Matrices for GABA DEGs CREs
#
# This script uses deepTools computeMatrix to compute signal matrices
# for CREs linked to up-regulated and down-regulated GABA DEGs.
#
# Comparisons (same as splicing genes analysis):
# 1. Nestin-Ctrl vs Nestin-Mut
# 2. Nestin-Ctrl vs Emx1-Mut
# 3. Nestin-Mut vs Emx1-Mut
#
# NOTE: Emx1-Ctrl is excluded (failed sample)
# NOTE: Uses GABA BigWig files ONLY for signal analysis
#
# Prerequisites:
# - output/GABA_DEGs_up_CREs.bed
# - output/GABA_DEGs_down_CREs.bed
# - BigWig files from Signac pipeline (signac_results_L1/bigwig_tracks_L1/)
#
# Output:
# - Signal matrices for up-DEGs (3 comparisons)
# - Signal matrices for down-DEGs (3 comparisons)
################################################################################

echo "========================================================================"
echo "COMPUTE ATAC SIGNAL MATRICES FOR GABA DEGs CREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/GABA_DEGs"

# Configuration
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR="./output"
BED_UP="./output/GABA_DEGs_up_CREs.bed"
BED_DOWN="./output/GABA_DEGs_down_CREs.bed"

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

if [ ! -f "$BED_UP" ]; then
    echo "ERROR: Up-DEGs BED file not found: $BED_UP"
    echo "Please run: python 1_extract_gaba_deg_CREs.py"
    exit 1
fi

if [ ! -f "$BED_DOWN" ]; then
    echo "ERROR: Down-DEGs BED file not found: $BED_DOWN"
    echo "Please run: python 1_extract_gaba_deg_CREs.py"
    exit 1
fi

N_CRES_UP=$(wc -l < $BED_UP)
N_CRES_DOWN=$(wc -l < $BED_DOWN)
echo "Found BED files:"
echo "  Up-DEGs CREs: $N_CRES_UP"
echo "  Down-DEGs CREs: $N_CRES_DOWN"
echo ""

# ============================================================================
# Step 2: Check BigWig files
# ============================================================================
echo "========================================================================"
echo "STEP 2: Checking BigWig files"
echo "========================================================================"
echo ""

MISSING=0

# Only check files we need (no Emx1-Ctrl)
for SAMPLE in "GABA_Nestin-Ctrl" "GABA_Nestin-Mut" "GABA_Emx1-Mut"; do
    BW_FILE="$BIGWIG_BASE/${SAMPLE}.bw"
    if [ -f "$BW_FILE" ]; then
        echo "  Found: ${SAMPLE}.bw"
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
echo "NOTE: Emx1-Ctrl excluded (failed sample)"
echo ""

# ============================================================================
# Step 3: Compute matrices for UP-DEGs
# ============================================================================
echo "========================================================================"
echo "STEP 3: Computing matrices for UP-regulated DEGs ($N_CRES_UP CREs)"
echo "========================================================================"
echo ""

# 3A: Nestin-Ctrl vs Nestin-Mut
echo "--- Comparison 1: Nestin-Ctrl vs Nestin-Mut ---"
computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_UP \
    -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Nestin-Mut.bw" \
    --samplesLabel "Nestin-Ctrl" "Nestin-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_up_nestin_ctrl_vs_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_up_nestin_ctrl_vs_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_up_nestin_ctrl_vs_mut.log

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for Up-DEGs Nestin-Ctrl vs Nestin-Mut"
    exit 1
fi
echo "Matrix computed: Up-DEGs Nestin-Ctrl vs Nestin-Mut"
echo ""

# 3B: Nestin-Ctrl vs Emx1-Mut
echo "--- Comparison 2: Nestin-Ctrl vs Emx1-Mut ---"
computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_UP \
    -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
    --samplesLabel "Nestin-Ctrl" "Emx1-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_up_nestin_ctrl_vs_emx1_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_up_nestin_ctrl_vs_emx1_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_up_nestin_ctrl_vs_emx1_mut.log

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for Up-DEGs Nestin-Ctrl vs Emx1-Mut"
    exit 1
fi
echo "Matrix computed: Up-DEGs Nestin-Ctrl vs Emx1-Mut"
echo ""

# 3C: Nestin-Mut vs Emx1-Mut
echo "--- Comparison 3: Nestin-Mut vs Emx1-Mut ---"
computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_UP \
    -S "$BIGWIG_BASE/GABA_Nestin-Mut.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
    --samplesLabel "Nestin-Mut" "Emx1-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_up_nestin_mut_vs_emx1_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_up_nestin_mut_vs_emx1_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_up_nestin_mut_vs_emx1_mut.log

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for Up-DEGs Nestin-Mut vs Emx1-Mut"
    exit 1
fi
echo "Matrix computed: Up-DEGs Nestin-Mut vs Emx1-Mut"
echo ""

# ============================================================================
# Step 4: Compute matrices for DOWN-DEGs
# ============================================================================
echo "========================================================================"
echo "STEP 4: Computing matrices for DOWN-regulated DEGs ($N_CRES_DOWN CREs)"
echo "========================================================================"
echo ""

# 4A: Nestin-Ctrl vs Nestin-Mut
echo "--- Comparison 1: Nestin-Ctrl vs Nestin-Mut ---"
computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_DOWN \
    -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Nestin-Mut.bw" \
    --samplesLabel "Nestin-Ctrl" "Nestin-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_down_nestin_ctrl_vs_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_down_nestin_ctrl_vs_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_down_nestin_ctrl_vs_mut.log

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for Down-DEGs Nestin-Ctrl vs Nestin-Mut"
    exit 1
fi
echo "Matrix computed: Down-DEGs Nestin-Ctrl vs Nestin-Mut"
echo ""

# 4B: Nestin-Ctrl vs Emx1-Mut
echo "--- Comparison 2: Nestin-Ctrl vs Emx1-Mut ---"
computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_DOWN \
    -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
    --samplesLabel "Nestin-Ctrl" "Emx1-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_down_nestin_ctrl_vs_emx1_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_down_nestin_ctrl_vs_emx1_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_down_nestin_ctrl_vs_emx1_mut.log

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for Down-DEGs Nestin-Ctrl vs Emx1-Mut"
    exit 1
fi
echo "Matrix computed: Down-DEGs Nestin-Ctrl vs Emx1-Mut"
echo ""

# 4C: Nestin-Mut vs Emx1-Mut
echo "--- Comparison 3: Nestin-Mut vs Emx1-Mut ---"
computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_DOWN \
    -S "$BIGWIG_BASE/GABA_Nestin-Mut.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
    --samplesLabel "Nestin-Mut" "Emx1-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_down_nestin_mut_vs_emx1_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_down_nestin_mut_vs_emx1_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_down_nestin_mut_vs_emx1_mut.log

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for Down-DEGs Nestin-Mut vs Emx1-Mut"
    exit 1
fi
echo "Matrix computed: Down-DEGs Nestin-Mut vs Emx1-Mut"
echo ""

# ============================================================================
# Summary
# ============================================================================
echo "========================================================================"
echo "MATRIX COMPUTATION COMPLETE!"
echo "========================================================================"
echo ""
echo "Output directory: $OUTPUT_DIR/"
echo ""
echo "Generated matrices:"
echo ""
echo "  UP-REGULATED DEGs ($N_CRES_UP CREs):"
echo "    matrix_up_nestin_ctrl_vs_mut.gz"
echo "    matrix_up_nestin_ctrl_vs_emx1_mut.gz"
echo "    matrix_up_nestin_mut_vs_emx1_mut.gz"
echo ""
echo "  DOWN-REGULATED DEGs ($N_CRES_DOWN CREs):"
echo "    matrix_down_nestin_ctrl_vs_mut.gz"
echo "    matrix_down_nestin_ctrl_vs_emx1_mut.gz"
echo "    matrix_down_nestin_mut_vs_emx1_mut.gz"
echo ""
echo "Parameters:"
echo "  Window size: +/- ${WINDOW_SIZE} bp"
echo "  Bin size: ${BIN_SIZE} bp"
echo "  Reference point: CRE center"
echo ""
echo "Next step: Run 3_visualize_deg_comparisons.py"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
