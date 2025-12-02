#!/bin/bash
#SBATCH --job-name=3_create_splicing_profiles
#SBATCH --output=logs/3_create_splicing_profiles.log
#SBATCH --error=logs/3_create_splicing_profiles.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=workq

################################################################################
# Create ATAC Signal Profiles for Splicing Gene CREs
#
# This script creates publication-quality metaprofiles showing ATAC signal
# at CREs linked to splicing genes, with focus on Ctrl vs Mut comparisons.
#
# PERFORMANCE OPTIONS:
# -----------------------------------------
# DEFAULT: Quick mode (metaprofiles only, skips individual plots):
#   sbatch 3_create_splicing_profiles.sh
#
# Full mode (create individual plots):
#   SKIP_INDIVIDUAL=0 sbatch 3_create_splicing_profiles.sh
#
# Full mode with parallel processing (8x faster):
#   SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 3_create_splicing_profiles.sh
#
# Full mode with lower DPI (2x faster):
#   SKIP_INDIVIDUAL=0 INDIVIDUAL_DPI=100 sbatch 3_create_splicing_profiles.sh
#
# Combined (parallel + low DPI):
#   SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 INDIVIDUAL_DPI=100 sbatch 3_create_splicing_profiles.sh
#
# Prerequisites:
# - output/splicing_genes_analysis/splicing_genes_CREs_GABA.bed
# - BigWig files from Signac pipeline (signac_results_L1/bigwig_tracks_L1/)
#
# Output:
# - Publication-quality metaprofiles with Ctrl vs Mut comparisons (always 300 DPI)
# - Difference plots (Mut - Ctrl)
# - Individual CRE profiles (optional, configurable DPI)
# - Statistical summaries
################################################################################

echo "========================================================================"
echo "CREATE ATAC SIGNAL PROFILES FOR SPLICING GENE CREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_genes_paper"

# Configuration
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR="./output/heatmaps_deeptools"
BED_ALL="./output/splicing_genes_CREs_all.bed"
# NOTE: Using ALL CREs (not just GABA-specific) but with GABA BigWig files ONLY for signal analysis

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
echo "STEP 1: Checking BED file"
echo "========================================================================"
echo ""

if [ ! -f "$BED_ALL" ]; then
    echo "ERROR: ALL CREs BED file not found: $BED_ALL"
    echo "Please run: sbatch 2_convert_splicing_CREs_to_bed.sh"
    exit 1
fi

N_CRES=$(wc -l < $BED_ALL)
echo "✓ Found ALL CREs BED file: $BED_ALL"
echo "  Total CREs: $N_CRES"
echo ""

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Copy BED file to output directory for reference
cp $BED_ALL $OUTPUT_DIR/splicing_genes_CREs_all.bed
echo "✓ Copied BED file to output directory"
echo ""

# Show CRE IDs and associated genes
echo "CREs in analysis:"
awk '{print "  - " $4 " (" $1 ":" $2 "-" $3 ")"}' $BED_ALL
echo ""

# ============================================================================
# Step 2: Check BigWig files
# ============================================================================
echo "========================================================================"
echo "STEP 2: Checking BigWig files"
echo "========================================================================"
echo ""

MISSING=0

for GENOTYPE in Nestin Emx1; do
    for CONDITION in Ctrl Mut; do
        BW_FILE="$BIGWIG_BASE/GABA_${GENOTYPE}-${CONDITION}.bw"
        if [ -f "$BW_FILE" ]; then
            echo "  ✓ Found: GABA_${GENOTYPE}-${CONDITION}.bw"
        else
            echo "  ✗ Missing: $BW_FILE"
            MISSING=$((MISSING + 1))
        fi
    done
done

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "ERROR: Missing $MISSING BigWig files!"
    exit 1
fi

echo ""

# ============================================================================
# Step 3: Run computeMatrix for Nestin (Ctrl vs Mut)
# ============================================================================
echo "========================================================================"
echo "STEP 3: Computing signal matrices for NESTIN (Ctrl vs Mut)"
echo "========================================================================"
echo ""

NESTIN_BIGWIGS="$BIGWIG_BASE/GABA_Nestin-Ctrl.bw $BIGWIG_BASE/GABA_Nestin-Mut.bw"
NESTIN_LABELS="Nestin-Ctrl Nestin-Mut"

echo "Running computeMatrix for Nestin..."
echo "  Window size: ±${WINDOW_SIZE} bp"
echo "  Bin size: ${BIN_SIZE} bp"
echo "  Processors: $N_PROCESSORS"
echo "  CREs: $N_CRES"
echo ""

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_ALL \
    -S $NESTIN_BIGWIGS \
    --samplesLabel $NESTIN_LABELS \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_GABA_nestin.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_GABA_nestin.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_nestin.log

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for Nestin"
    exit 1
fi

echo ""
echo "✓ Nestin matrix computed and saved"
echo ""

# ============================================================================
# Step 4: Run computeMatrix for Emx1 (Ctrl vs Mut)
# ============================================================================
echo "========================================================================"
echo "STEP 4: Computing signal matrices for EMX1 (Ctrl vs Mut)"
echo "========================================================================"
echo ""

EMX1_BIGWIGS="$BIGWIG_BASE/GABA_Nestin-Ctrl.bw $BIGWIG_BASE/GABA_Emx1-Mut.bw"
EMX1_LABELS="Nestin-Ctrl Emx1-Mut"

echo "Running computeMatrix for Emx1..."
echo "  Window size: ±${WINDOW_SIZE} bp"
echo "  Bin size: ${BIN_SIZE} bp"
echo "  Processors: $N_PROCESSORS"
echo "  CREs: $N_CRES"
echo ""

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_ALL \
    -S $EMX1_BIGWIGS \
    --samplesLabel $EMX1_LABELS \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_GABA_emx1.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_GABA_emx1.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_emx1.log

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for Emx1"
    exit 1
fi

echo ""
echo "✓ Emx1 matrix computed and saved"
echo ""

# ============================================================================
# Step 5: Create publication-quality visualizations using Python
# ============================================================================
echo "========================================================================"
echo "STEP 5: Creating publication-quality visualizations"
echo "========================================================================"
echo ""

# Performance options (can be set via environment variables):
# - SKIP_INDIVIDUAL=0: Create individual CRE plots (default: skip to save 30-50 minutes)
# - PARALLEL_JOBS=N: Use N parallel processes for individual plots
# - INDIVIDUAL_DPI=N: DPI for individual plots (default: 150)

SKIP_FLAG="--skip-individual"
if [ "${SKIP_INDIVIDUAL}" = "0" ]; then
    echo "📊 Full mode: Creating individual plots for each CRE"
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
FILTER_FLAG=""
if [ -n "${MIN_SIGNAL}" ]; then
    FILTER_FLAG="$FILTER_FLAG --min-signal ${MIN_SIGNAL}"
    echo "   Min Signal: ${MIN_SIGNAL}"
fi
if [ -n "${MIN_FC}" ]; then
    FILTER_FLAG="$FILTER_FLAG --min-fc ${MIN_FC}"
    echo "   Min FC: ${MIN_FC}"
fi

echo ""
echo "Running Python visualization script..."
python 3_create_splicing_profiles.py $SKIP_FLAG $PARALLEL_FLAG $DPI_FLAG $FILTER_FLAG

if [ $? -ne 0 ]; then
    echo "ERROR: Visualization failed!"
    exit 1
fi

echo ""
echo "✓ Visualizations created"
echo ""

# ============================================================================
# Step 6: Create summary report
# ============================================================================
echo "========================================================================"
echo "STEP 6: Creating summary report"
echo "========================================================================"
echo ""

# Get list of genes
GENE_LIST=$(awk -F'\t' 'NR>1 {print $1}' ./output/splicing_genes_CREs_all_celltypes.tsv | sort -u | tr '\n' ', ' | sed 's/,$//')
N_GENES=$(awk -F'\t' 'NR>1 {print $1}' ./output/splicing_genes_CREs_all_celltypes.tsv | sort -u | wc -l)

cat > $OUTPUT_DIR/README.txt << EOFR
================================================================================
PUBLICATION-QUALITY ATAC SIGNAL PROFILES: SPLICING GENE CREs
================================================================================

Analysis date: $(date)

PURPOSE:
--------
Visualize ATAC-seq signal at CREs linked to splicing-related genes with
focus on identifying Ctrl vs Mut differences at regulatory elements controlling
splicing machinery.

SPLICING GENES ANALYZED:
-------------------------
Total genes: $N_GENES
Gene list: $GENE_LIST

These genes are part of:
- mRNA splicing pathways
- Spliceosome assembly
- Pre-mRNA processing

CREs ANALYZED:
--------------
Total CREs: $N_CRES (all in GABA neurons)
All CREs are shared between GABA and other cell types (not GABA-specific)

GENOTYPE COMPARISONS:
---------------------
1. Nestin: Nestin-Ctrl vs Nestin-Mut
2. Emx1: Nestin-Ctrl vs Emx1-Mut (Emx1-Ctrl is a failed sample)

PARAMETERS:
-----------
Window size: ±${WINDOW_SIZE} bp around CRE center
Bin size: ${BIN_SIZE} bp
Reference point: CRE center
Processors: $N_PROCESSORS (parallel processing)

OUTPUT FILES:
-------------
Metaprofiles (Ctrl vs Mut comparison):
  - profiles/metaprofile_nestin_ctrl_vs_mut.png
      → Overlaid Ctrl vs Mut signals with SEM bands
      → Difference plot (Mut - Ctrl) below
      → Statistical annotations

  - profiles/metaprofile_emx1_ctrl_vs_mut.png
      → Same format as Nestin

Individual CRE profiles:
  - profiles/individual_Nestin_<Gene>_<CRE_ID>.png
      → One plot per CRE showing Nestin-Ctrl vs Nestin-Mut

  - profiles/individual_Emx1_<Gene>_<CRE_ID>.png
      → One plot per CRE showing Nestin-Ctrl vs Emx1-Mut

Signal matrices (reusable):
  - matrix_GABA_nestin.gz / matrix_GABA_nestin.tab
  - matrix_GABA_emx1.gz / matrix_GABA_emx1.tab

INTERPRETATION GUIDE:
---------------------
Metaprofiles:
  - Blue line: Ctrl samples (mean ± SEM)
  - Red line: Mut samples (mean ± SEM)
  - Shaded bands: Standard error of the mean
  - Bottom panel: Difference (Mut - Ctrl)
      • Red shading: Regions where Mut > Ctrl
      • Blue shading: Regions where Ctrl > Mut

Individual CRE plots:
  - Each CRE plotted separately
  - Allows identification of CRE-specific effects
  - Useful for follow-up validation

KEY QUESTIONS TO ADDRESS:
--------------------------
1. Do Mut samples show increased or decreased accessibility?
   → Check if red line is above (increased) or below (decreased) blue line

2. Are changes CRE-centered or broad?
   → Narrow peak at center = CRE-specific
   → Broad changes = regional accessibility effects

3. Are changes consistent across genotypes?
   → Compare Nestin vs Emx1 difference plots
   → Consistent pattern = robust effect
   → Genotype-specific = context-dependent regulation

4. Which CREs show strongest differences?
   → Check individual CRE plots
   → Identify candidates for validation

BIOLOGICAL CONTEXT:
-------------------
Splicing machinery is essential for:
- mRNA processing and maturation
- Alternative splicing regulation
- Gene expression control

Changes in CRE accessibility at splicing genes may indicate:
- Altered splicing factor expression
- Dysregulated mRNA processing
- Compensatory responses to mutations

NEXT STEPS:
-----------
1. Identify CREs with largest Ctrl vs Mut differences
2. Correlate CRE accessibility with splicing gene expression (RNA-seq)
3. Look for changes in alternative splicing patterns
4. Validate findings with targeted assays (ChIP, ATAC-qPCR)

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
echo "CREs analyzed: $N_CRES CREs (ALL cell types) linked to $N_GENES splicing genes"
echo "Signal analysis: GABA BigWig files ONLY (Nestin/Emx1 Ctrl vs Mut)"
echo "Analysis approach: ALL CREs with GABA BigWig signal"
echo ""
echo "Generated files:"
echo ""
echo "  📊 Metaprofiles (Ctrl vs Mut comparison):"
echo "      ├─ profiles/metaprofile_nestin_ctrl_vs_mut.png"
echo "      │   → Nestin: Ctrl vs Mut with difference plot"
echo "      └─ profiles/metaprofile_emx1_ctrl_vs_mut.png"
echo "          → Emx1: Ctrl vs Mut with difference plot"
echo ""
echo "  📈 Individual CRE plots:"
echo "      ├─ profiles/individual_Nestin_*.png ($N_CRES plots)"
echo "      └─ profiles/individual_Emx1_*.png ($N_CRES plots)"
echo ""
echo "  📋 Documentation:"
echo "      └─ README.txt (interpretation guide)"
echo ""
echo "WHAT TO LOOK FOR:"
echo "  ✓ Red line above blue = Increased accessibility in Mut"
echo "  ✓ Blue line above red = Decreased accessibility in Mut"
echo "  ✓ Bottom panel shows net difference (Mut - Ctrl)"
echo "  ✓ Individual plots identify CRE-specific effects"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
