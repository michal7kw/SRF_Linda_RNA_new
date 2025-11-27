#!/bin/bash
#SBATCH --job-name=GABA_DEG_analysis
#SBATCH --output=logs/0_RUN_GABA_DEG_ANALYSIS.log
#SBATCH --error=logs/0_RUN_GABA_DEG_ANALYSIS.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=workq

################################################################################
# GABA DEGs CRE Analysis Pipeline
#
# This master script runs the complete analysis pipeline for visualizing
# ATAC signal at CREs linked to GABA differentially expressed genes.
#
# PIPELINE STEPS:
# 1. Extract CREs linked to GABA DEGs from Table 16 literature data
# 2. Compute deepTools signal matrices for up/down DEGs
# 3. Create visualization profiles (metaprofiles and individual CRE plots)
#
# INPUT:
# - GABA DEG lists (up_significant.csv, down_significant.csv)
# - Table 16 CRE-gene correlations
# - GABA BigWig files (Nestin-Ctrl, Nestin-Mut, Emx1-Mut)
#
# OUTPUT:
# - CRE TSV and BED files for up/down DEGs
# - deepTools signal matrices
# - Metaprofile visualizations (6 plots: 2 DEG sets x 3 comparisons)
# - Statistical summaries
#
# PERFORMANCE OPTIONS:
# -----------------------------------------
# DEFAULT: Quick mode (metaprofiles only):
#   sbatch 0_RUN_GABA_DEG_ANALYSIS.sh
#
# Full mode (with individual CRE plots):
#   SKIP_INDIVIDUAL=0 sbatch 0_RUN_GABA_DEG_ANALYSIS.sh
#
# Full mode with parallel processing:
#   SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 0_RUN_GABA_DEG_ANALYSIS.sh
#
# USAGE:
#   sbatch 0_RUN_GABA_DEG_ANALYSIS.sh
################################################################################

echo "========================================================================"
echo "GABA DEG CRE ANALYSIS PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/GABA_DEGs"

# Create logs directory
mkdir -p logs

# ============================================================================
# STEP 1: Extract CREs for GABA DEGs
# ============================================================================
echo "========================================================================"
echo "STEP 1: Extracting CREs linked to GABA DEGs"
echo "========================================================================"
echo ""

# Activate environment for Table 16 parsing
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "Running CRE extraction script..."
python 1_extract_gaba_deg_CREs.py 2>&1 | tee logs/1_extract_gaba_deg_CREs.log

if [ $? -ne 0 ]; then
    echo "ERROR: CRE extraction failed!"
    exit 1
fi

echo ""
echo "Step 1 completed successfully"
echo ""

# Check output files
if [ ! -f "./output/GABA_DEGs_up_CREs.bed" ] || [ ! -f "./output/GABA_DEGs_down_CREs.bed" ]; then
    echo "ERROR: BED files not generated!"
    exit 1
fi

N_UP=$(wc -l < ./output/GABA_DEGs_up_CREs.bed)
N_DOWN=$(wc -l < ./output/GABA_DEGs_down_CREs.bed)
echo "CREs extracted:"
echo "  Up-regulated DEGs: $N_UP CREs"
echo "  Down-regulated DEGs: $N_DOWN CREs"
echo ""

# ============================================================================
# STEP 2: Compute Signal Matrices
# ============================================================================
echo "========================================================================"
echo "STEP 2: Computing deepTools signal matrices"
echo "========================================================================"
echo ""

# Check if deepTools is available
if ! command -v computeMatrix &> /dev/null; then
    echo "ERROR: deepTools not found!"
    exit 1
fi

echo "deepTools found: $(computeMatrix --version)"
echo ""

# Configuration
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR="./output"
WINDOW_SIZE=2000
BIN_SIZE=50
N_PROCESSORS=16

# Check BigWig files
for SAMPLE in "GABA_Nestin-Ctrl" "GABA_Nestin-Mut" "GABA_Emx1-Mut"; do
    if [ ! -f "$BIGWIG_BASE/${SAMPLE}.bw" ]; then
        echo "ERROR: BigWig file not found: $BIGWIG_BASE/${SAMPLE}.bw"
        exit 1
    fi
done
echo "All BigWig files found"
echo ""

# Compute matrices for UP-DEGs
echo "--- Computing matrices for UP-regulated DEGs ---"

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R ./output/GABA_DEGs_up_CREs.bed \
    -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Nestin-Mut.bw" \
    --samplesLabel "Nestin-Ctrl" "Nestin-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_up_nestin_ctrl_vs_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_up_nestin_ctrl_vs_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_up_nestin_ctrl_vs_mut.log

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R ./output/GABA_DEGs_up_CREs.bed \
    -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
    --samplesLabel "Nestin-Ctrl" "Emx1-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_up_nestin_ctrl_vs_emx1_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_up_nestin_ctrl_vs_emx1_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_up_nestin_ctrl_vs_emx1_mut.log

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R ./output/GABA_DEGs_up_CREs.bed \
    -S "$BIGWIG_BASE/GABA_Nestin-Mut.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
    --samplesLabel "Nestin-Mut" "Emx1-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_up_nestin_mut_vs_emx1_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_up_nestin_mut_vs_emx1_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_up_nestin_mut_vs_emx1_mut.log

echo "Up-DEGs matrices computed"
echo ""

# Compute matrices for DOWN-DEGs
echo "--- Computing matrices for DOWN-regulated DEGs ---"

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R ./output/GABA_DEGs_down_CREs.bed \
    -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Nestin-Mut.bw" \
    --samplesLabel "Nestin-Ctrl" "Nestin-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_down_nestin_ctrl_vs_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_down_nestin_ctrl_vs_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_down_nestin_ctrl_vs_mut.log

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R ./output/GABA_DEGs_down_CREs.bed \
    -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
    --samplesLabel "Nestin-Ctrl" "Emx1-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_down_nestin_ctrl_vs_emx1_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_down_nestin_ctrl_vs_emx1_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_down_nestin_ctrl_vs_emx1_mut.log

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R ./output/GABA_DEGs_down_CREs.bed \
    -S "$BIGWIG_BASE/GABA_Nestin-Mut.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
    --samplesLabel "Nestin-Mut" "Emx1-Mut" \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_down_nestin_mut_vs_emx1_mut.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_down_nestin_mut_vs_emx1_mut.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_down_nestin_mut_vs_emx1_mut.log

echo "Down-DEGs matrices computed"
echo ""
echo "Step 2 completed successfully"
echo ""

# ============================================================================
# STEP 3: Create Visualizations
# ============================================================================
echo "========================================================================"
echo "STEP 3: Creating visualizations"
echo "========================================================================"
echo ""

# Performance options (can be set via environment variables):
# - SKIP_INDIVIDUAL=0: Create individual CRE plots (default: skip)
# - PARALLEL_JOBS=N: Use N parallel processes for individual plots
# - INDIVIDUAL_DPI=N: DPI for individual plots (default: 150)

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
    echo "Using ${PARALLEL_JOBS} parallel processes"
fi

DPI_FLAG=""
if [ -n "${INDIVIDUAL_DPI}" ]; then
    DPI_FLAG="--individual-dpi ${INDIVIDUAL_DPI}"
    echo "Individual plot DPI: ${INDIVIDUAL_DPI}"
fi

echo ""
echo "Running visualization script..."
python 3_visualize_deg_comparisons.py $SKIP_FLAG $PARALLEL_FLAG $DPI_FLAG 2>&1 | tee logs/3_visualize_deg_comparisons.log

if [ $? -ne 0 ]; then
    echo "ERROR: Visualization failed!"
    exit 1
fi

echo ""
echo "Step 3 completed successfully"
echo ""

# ============================================================================
# Create README
# ============================================================================
echo "========================================================================"
echo "Creating documentation"
echo "========================================================================"
echo ""

cat > ./output/README.txt << 'EOFR'
================================================================================
GABA DEG CRE ANALYSIS RESULTS
================================================================================

PURPOSE:
--------
Visualization of ATAC-seq signal at CREs (cis-regulatory elements) linked to
GABA-specific differentially expressed genes.

DEGs analyzed:
- Up-regulated: Genes enriched in GABA neurons vs Rest (GABA markers)
- Down-regulated: Genes depleted in GABA neurons vs Rest

DATA SOURCES:
-------------
1. GABA DEGs:
   - Up: Cluster_GABA_vs_Rest_up_significant.csv
   - Down: Cluster_GABA_vs_Rest_down_significant.csv

2. CRE-gene correlations: Table 16 (literature data)
   - Filters: FDR < 0.05, |PCC| > 0.2
   - Cell types: GABA/Hippocampal subtypes

3. ATAC signal: GABA BigWig files from Signac pipeline

COMPARISONS:
------------
For each DEG set (up/down), three comparisons are made:

1. Nestin-Ctrl vs Nestin-Mut
   -> Within-genotype mutation effect

2. Nestin-Ctrl vs Emx1-Mut
   -> Cross-genotype comparison

3. Nestin-Mut vs Emx1-Mut
   -> Mutant-to-mutant genotype effect

NOTE: Emx1-Ctrl excluded (failed sample)

OUTPUT FILES:
-------------
profiles_up/
  metaprofile_nestin_ctrl_vs_mut.png
  metaprofile_nestin_ctrl_vs_emx1_mut.png
  metaprofile_nestin_mut_vs_emx1_mut.png

profiles_down/
  metaprofile_nestin_ctrl_vs_mut.png
  metaprofile_nestin_ctrl_vs_emx1_mut.png
  metaprofile_nestin_mut_vs_emx1_mut.png

GABA_DEGs_up_CREs.tsv / .bed - CREs for up-regulated DEGs
GABA_DEGs_down_CREs.tsv / .bed - CREs for down-regulated DEGs
comparison_statistics.txt - Statistical summaries
SUMMARY_GABA_DEGs_CREs.txt - Analysis summary

INTERPRETATION:
---------------
Up-regulated DEGs (GABA markers):
- CREs linked to genes highly expressed in GABA neurons
- Changes may indicate altered GABA identity or function

Down-regulated DEGs:
- CREs linked to genes suppressed in GABA neurons
- Changes may indicate de-repression or cell type switching

Metaprofiles:
- Top panel: Mean ATAC signal (+/- SEM) for each condition
- Bottom panel: Difference (Condition2 - Condition1)
- Positive difference = increased accessibility in Condition2
- Negative difference = decreased accessibility in Condition2

Generated by: 0_RUN_GABA_DEG_ANALYSIS.sh
================================================================================
EOFR

echo "README created"
echo ""

# ============================================================================
# Final Summary
# ============================================================================
echo "========================================================================"
echo "ANALYSIS COMPLETE!"
echo "========================================================================"
echo ""
echo "Output directory: ./output/"
echo ""
echo "CREs analyzed:"
echo "  Up-regulated DEGs: $N_UP CREs"
echo "  Down-regulated DEGs: $N_DOWN CREs"
echo ""
echo "Generated visualizations:"
echo "  Up-regulated DEGs: output/profiles_up/"
echo "    - metaprofile_nestin_ctrl_vs_mut.png"
echo "    - metaprofile_nestin_ctrl_vs_emx1_mut.png"
echo "    - metaprofile_nestin_mut_vs_emx1_mut.png"
echo ""
echo "  Down-regulated DEGs: output/profiles_down/"
echo "    - metaprofile_nestin_ctrl_vs_mut.png"
echo "    - metaprofile_nestin_ctrl_vs_emx1_mut.png"
echo "    - metaprofile_nestin_mut_vs_emx1_mut.png"
echo ""
echo "Statistics: output/comparison_statistics.txt"
echo "Summary: output/SUMMARY_GABA_DEGs_CREs.txt"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
