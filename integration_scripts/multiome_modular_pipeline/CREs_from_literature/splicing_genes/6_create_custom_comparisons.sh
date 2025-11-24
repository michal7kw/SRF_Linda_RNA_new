#!/bin/bash
#SBATCH --job-name=6_create_custom_comparisons
#SBATCH --output=logs/6_create_custom_comparisons.log
#SBATCH --error=logs/6_create_custom_comparisons.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=workq

################################################################################
# Custom ATAC Signal Comparisons for Splicing Gene CREs
#
# This script creates custom comparisons across genotypes and conditions:
# 1. Nestin-Ctrl vs Nestin-Mut
# 2. Nestin-Ctrl vs Emx1-Mut
# 3. Nestin-Mut vs Emx1-Mut
#
# NOTE: Emx1-Ctrl is excluded (failed sample)
# NOTE: Uses ALL CREs (not just GABA-specific) but with GABA BigWig files ONLY for signal analysis
#
# PERFORMANCE OPTIONS:
# -----------------------------------------
# DEFAULT: Quick mode (metaprofiles only, skips individual plots):
#   sbatch 6_create_custom_comparisons.sh
#
# Full mode (create individual plots):
#   SKIP_INDIVIDUAL=0 sbatch 6_create_custom_comparisons.sh
#
# Full mode with parallel processing (8x faster):
#   SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 6_create_custom_comparisons.sh
#
# Full mode with lower DPI (2x faster):
#   SKIP_INDIVIDUAL=0 INDIVIDUAL_DPI=100 sbatch 6_create_custom_comparisons.sh
#
# Combined (parallel + low DPI):
#   SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 INDIVIDUAL_DPI=100 sbatch 6_create_custom_comparisons.sh
#
# Prerequisites:
# - output/splicing_genes_analysis/splicing_genes_CREs_all.bed
# - BigWig files from Signac pipeline (signac_results_L1/bigwig_tracks_L1/)
#
# Output:
# - Custom comparison metaprofiles (always 300 DPI)
# - Difference plots for each comparison
# - Individual CRE profiles (optional, configurable DPI)
# - Statistical summaries
################################################################################

echo "========================================================================"
echo "CUSTOM ATAC SIGNAL COMPARISONS FOR SPLICING GENE CREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/splicing_genes"

# Configuration
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR="./output/custom_comparisons"
BED_ALL="./output/splicing_genes_CREs_all.bed"
# NOTE: Using ALL CREs (not just GABA-specific) but with GABA BigWig files for signal analysis

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
# Step 1: Check BED file exists
# ============================================================================
echo "========================================================================"
echo "STEP 1: Checking BED file (ALL CREs, not just GABA-specific)"
echo "========================================================================"
echo ""

if [ ! -f "$BED_ALL" ]; then
    echo "ERROR: BED file not found: $BED_ALL"
    echo "Please run: sbatch 2_convert_splicing_CREs_to_bed.sh"
    exit 1
fi

N_CRES=$(wc -l < $BED_ALL)
echo "✓ Found BED file: $BED_ALL"
echo "  Total CREs: $N_CRES (all cell types)"
echo ""

# Create output directory
mkdir -p $OUTPUT_DIR

# Copy BED file to output directory for reference
cp $BED_ALL $OUTPUT_DIR/splicing_genes_CREs_all.bed
echo "✓ Copied BED file to output directory"
echo ""

# Show CRE IDs
echo "CREs in analysis:"
awk '{print "  - " $4 " (" $1 ":" $2 "-" $3 ")"}' $BED_ALL
echo ""

# ============================================================================
# Step 2: Check BigWig files (excluding Emx1-Ctrl)
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
        echo "  ✓ Found: ${SAMPLE}.bw"
    else
        echo "  ✗ Missing: $BW_FILE"
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

echo "✓ Matrix computed: Nestin-Ctrl vs Nestin-Mut"
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

echo "✓ Matrix computed: Nestin-Ctrl vs Emx1-Mut"
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

echo "✓ Matrix computed: Nestin-Mut vs Emx1-Mut"
echo ""

# ============================================================================
# Step 4: Create visualizations using Python
# ============================================================================
echo "========================================================================"
echo "STEP 4: Creating custom comparison visualizations"
echo "========================================================================"
echo ""

# Performance options (can be set via environment variables):
# - SKIP_INDIVIDUAL=0: Create individual CRE plots (default: skip to save 30-50 minutes)
# - PARALLEL_JOBS=N: Use N parallel processes for individual plots
# - INDIVIDUAL_DPI=N: DPI for individual plots (default: 150)

SKIP_FLAG="--skip-individual"
if [ "${SKIP_INDIVIDUAL}" = "0" ]; then
    echo "📊 Full mode: Creating individual plots for each comparison"
    SKIP_FLAG=""
else
    echo "⚡ Fast mode (DEFAULT): Skipping individual CRE plots"
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
python 6_visualize_custom_comparisons.py $SKIP_FLAG $PARALLEL_FLAG $DPI_FLAG

if [ $? -ne 0 ]; then
    echo "ERROR: Visualization failed!"
    exit 1
fi

echo ""
echo "✓ Visualizations created"
echo ""

# ============================================================================
# Step 5: Create summary report
# ============================================================================
echo "========================================================================"
echo "STEP 5: Creating summary report"
echo "========================================================================"
echo ""

# Get gene list
GENE_LIST=$(awk -F'\t' 'NR>1 {print $1}' ../output/splicing_genes_analysis/splicing_genes_CREs_all_celltypes.tsv | sort -u | tr '\n' ', ' | sed 's/,$//')
N_GENES=$(awk -F'\t' 'NR>1 {print $1}' ../output/splicing_genes_analysis/splicing_genes_CREs_all_celltypes.tsv | sort -u | wc -l)

cat > $OUTPUT_DIR/README.txt << 'EOFR'
================================================================================
CUSTOM ATAC SIGNAL COMPARISONS: SPLICING GENE CREs
================================================================================

Analysis date: $(date)

PURPOSE:
--------
Custom comparisons of ATAC-seq signal at CREs linked to splicing genes,
focusing on cross-genotype and cross-condition comparisons while excluding
the failed Emx1-Ctrl sample.

COMPARISONS PERFORMED:
----------------------
1. Nestin-Ctrl vs Nestin-Mut
   → Within-genotype comparison (effect of mutation in Nestin background)

2. Nestin-Ctrl vs Emx1-Mut
   → Cross-genotype comparison (wild-type Nestin vs mutant Emx1)

3. Nestin-Mut vs Emx1-Mut
   → Mutant-to-mutant comparison (genotype effect under mutation)

NOTE: Emx1-Ctrl excluded as it is a failed sample

CREs ANALYZED:
--------------
Total CREs: 10 (ALL cell types, not restricted to GABA)
All CREs are linked to splicing machinery genes

SPLICING GENES:
---------------
Total genes: 10
Gene list: Celf4, Fus, Khdrbs3, Nova1, Rbfox2, Rbm5, Sfpq, Srsf1, Srsf2, Srsf5

PARAMETERS:
-----------
Window size: ±2000 bp around CRE center
Bin size: 50 bp
Reference point: CRE center
Processors: 16 (parallel processing)

OUTPUT FILES:
-------------
Metaprofiles (comparison overlays):
  - profiles/metaprofile_nestin_ctrl_vs_nestin_mut.png
      → Line 1 (blue): Nestin-Ctrl
      → Line 2 (red): Nestin-Mut
      → Bottom panel: Difference (Nestin-Mut - Nestin-Ctrl)

  - profiles/metaprofile_nestin_ctrl_vs_emx1_mut.png
      → Line 1 (blue): Nestin-Ctrl
      → Line 2 (orange): Emx1-Mut
      → Bottom panel: Difference (Emx1-Mut - Nestin-Ctrl)

  - profiles/metaprofile_nestin_mut_vs_emx1_mut.png
      → Line 1 (red): Nestin-Mut
      → Line 2 (orange): Emx1-Mut
      → Bottom panel: Difference (Emx1-Mut - Nestin-Mut)

Individual CRE profiles:
  - profiles/individual_comparison1_<Gene>_<CRE_ID>.png (10 plots)
  - profiles/individual_comparison2_<Gene>_<CRE_ID>.png (10 plots)
  - profiles/individual_comparison3_<Gene>_<CRE_ID>.png (10 plots)

Signal matrices (reusable):
  - matrix_nestin_ctrl_vs_mut.gz / .tab
  - matrix_nestin_ctrl_vs_emx1_mut.gz / .tab
  - matrix_nestin_mut_vs_emx1_mut.gz / .tab

INTERPRETATION GUIDE:
---------------------
Comparison 1 (Nestin-Ctrl vs Nestin-Mut):
  - Shows effect of mutation within Nestin genotype
  - Red > Blue = Increased accessibility in Nestin mutant
  - Blue > Red = Decreased accessibility in Nestin mutant

Comparison 2 (Nestin-Ctrl vs Emx1-Mut):
  - Shows cross-genotype effect (wild-type vs mutant)
  - Orange > Blue = Higher accessibility in Emx1-Mut than Nestin-Ctrl
  - Identifies genotype-specific baseline differences + mutation effects

Comparison 3 (Nestin-Mut vs Emx1-Mut):
  - Shows genotype effect under mutation
  - Orange > Red = Emx1-Mut has higher accessibility than Nestin-Mut
  - Reveals whether mutation affects different genotypes similarly

KEY QUESTIONS TO ADDRESS:
--------------------------
1. Does the mutation affect splicing gene CREs consistently?
   → Compare Comparison 1 (within Nestin)
   → If strong effect, should see clear difference

2. Are there baseline genotype differences?
   → Compare Comparison 2 (Nestin-Ctrl vs Emx1-Mut)
   → Difference could be: genotype baseline + mutation effect

3. Do mutations in different genotypes have similar effects?
   → Compare Comparison 3 (Nestin-Mut vs Emx1-Mut)
   → Similar patterns = mutation drives effect
   → Different patterns = genotype-specific responses

4. Which comparison shows the strongest differences?
   → Identifies primary driver: mutation effect vs genotype effect

BIOLOGICAL CONTEXT:
-------------------
Splicing machinery genes control:
- mRNA processing and maturation
- Alternative splicing regulation
- Gene expression fine-tuning

Changes at splicing gene CREs may indicate:
- Altered expression of splicing factors
- Compensatory responses to SRF mutations
- Genotype-specific regulatory rewiring

NEXT STEPS:
-----------
1. Identify comparisons with largest differences
2. Correlate CRE accessibility with splicing gene expression (RNA-seq DEGs)
3. Check for alternative splicing changes in RNA-seq data
4. Validate with targeted assays (ChIP-qPCR, ATAC-qPCR)

Generated by: 6_create_custom_comparisons.sh
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
echo "CREs analyzed: $N_CRES CREs (ALL cell types, not just GABA-specific)"
echo "Signal analysis: GABA BigWig files ONLY (Nestin/Emx1 Ctrl vs Mut)"
echo "Genes: $N_GENES splicing genes"
echo ""
echo "Generated files:"
echo ""
echo "  📊 Metaprofiles (3 custom comparisons):"
echo "      ├─ profiles/metaprofile_nestin_ctrl_vs_nestin_mut.png"
echo "      ├─ profiles/metaprofile_nestin_ctrl_vs_emx1_mut.png"
echo "      └─ profiles/metaprofile_nestin_mut_vs_emx1_mut.png"
echo ""
echo "  📈 Individual CRE plots (30 plots total):"
echo "      ├─ profiles/individual_comparison1_*.png (10 plots)"
echo "      ├─ profiles/individual_comparison2_*.png (10 plots)"
echo "      └─ profiles/individual_comparison3_*.png (10 plots)"
echo ""
echo "  📋 Documentation:"
echo "      └─ README.txt (interpretation guide)"
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
