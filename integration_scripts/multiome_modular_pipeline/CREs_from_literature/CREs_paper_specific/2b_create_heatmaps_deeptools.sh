#!/bin/bash
#SBATCH --job-name=2b_create_heatmaps_deeptools
#SBATCH --output=logs/2b_create_heatmaps_deeptools.log
#SBATCH --error=logs/2b_create_heatmaps_deeptools.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Create Heatmaps and Metaprofiles Using deepTools
#
# This script uses deepTools (computeMatrix + plotHeatmap + plotProfile) to
# create publication-ready heatmaps and metaprofiles of ATAC signal at:
# 1. GABA-specific CREs (POSITIVE CONTROL - expect HIGH signal)
# 2. Excitatory neuron CREs (NEGATIVE CONTROL - expect LOW signal)
#
# This demonstrates cell-type specificity of ATAC signal in GABA samples.
#
# Uses existing BED files from:
# - 1_extract_hippocampal_interneuron_CREs.py (GABA)
# - 2a_extract_excitatory_neuron_CREs.py (Excitatory)
#
# Requirements:
# - deepTools installed (conda install -c bioconda deeptools)
# - BED files of GABA and Excitatory CREs (already created)
# - BigWig files (already available)
#
# INPUT FILES:
# ------------
# 1. BED files:
#    a) GABA CREs (~239,691 CREs):
#       - output/hippocampal_interneuron_CREs.bed
#       - Created by: 1_extract_hippocampal_interneuron_CREs.py
#       - Format: chr, start, end, cCRE_ID, num_celltypes, celltype_list
#
#    b) Excitatory CREs (expected ~similar number):
#       - output/excitatory_neuron_CREs.bed
#       - Created by: 2a_extract_excitatory_neuron_CREs.py
#       - Format: chr, start, end, cCRE_ID, num_celltypes, celltype_list
#
# 2. BigWig files (ATAC signal, 4 conditions):
#    - ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Nestin-Ctrl.bw
#    - ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Nestin-Mut.bw
#    - ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Emx1-Ctrl.bw
#    - ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Emx1-Mut.bw
#    - Created by: signac_06_create_bigwig_tracks.R
#
# OUTPUT FILES:
# -------------
# All outputs in: output/GABA_DEG_analysis/heatmaps_deeptools/
#
# 1. BED files (copies for reference):
#    - GABA_CREs.bed
#    - Excitatory_CREs.bed
#
# 2. Signal matrices (compressed, reusable):
#    GABA CREs (positive control):
#      - matrix_GABA_all_conditions.gz
#      - matrix_GABA_nestin.gz
#      - matrix_GABA_emx1.gz
#
#    Excitatory CREs (negative control):
#      - matrix_Excitatory_all_conditions.gz
#      - matrix_Excitatory_nestin.gz
#      - matrix_Excitatory_emx1.gz
#
# 3. Heatmaps (PNG):
#    - heatmap_GABA_vs_Excitatory_all_conditions.png (side-by-side comparison)
#    - heatmap_GABA_all_conditions.png
#    - heatmap_Excitatory_all_conditions.png
#    - heatmap_GABA_nestin.png / heatmap_Excitatory_nestin.png
#    - heatmap_GABA_emx1.png / heatmap_Excitatory_emx1.png
#
# 4. Metaprofiles (PNG):
#    - metaprofile_GABA_vs_Excitatory_all_conditions.png (overlay comparison)
#    - metaprofile_GABA_all_conditions.png
#    - metaprofile_Excitatory_all_conditions.png
#    - metaprofile_GABA_nestin.png / metaprofile_Excitatory_nestin.png
#    - metaprofile_GABA_emx1.png / metaprofile_Excitatory_emx1.png
#
# 5. Logs and metadata:
#    - computeMatrix_GABA.log / computeMatrix_Excitatory.log
#    - plotHeatmap_*.log / plotProfile_*.log
#    - README.txt (analysis summary and parameters)
#
# 6. SLURM logs:
#    - logs/heatmaps_deeptools_<jobid>.log (standard output)
#    - logs/heatmaps_deeptools_<jobid>.err (standard error)
#
################################################################################

echo "========================================================================"
echo "CREATE HEATMAPS AND METAPROFILES USING DEEPTOOLS"
echo "POSITIVE CONTROL: GABA CREs (expect HIGH signal)"
echo "NEGATIVE CONTROL: Excitatory CREs (expect LOW signal)"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Configuration
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR="output/GABA_DEG_analysis/heatmaps_deeptools"
BED_GABA="output/hippocampal_interneuron_CREs.bed"
BED_EXCITATORY="output/excitatory_neuron_CREs.bed"

# Parameters
WINDOW_SIZE=2000  # bp around CRE center (±2kb)
BIN_SIZE=50       # bp per bin
N_PROCESSORS=8

# Create output directory
mkdir -p $OUTPUT_DIR
mkdir -p logs

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
echo "------------------------------------------------------------------------"

# Check GABA BED file
if [ ! -f "$BED_GABA" ]; then
    echo "ERROR: GABA BED file not found: $BED_GABA"
    echo "Expected file created by 1_extract_hippocampal_interneuron_CREs.py"
    exit 1
fi

N_GABA_CRES=$(wc -l < $BED_GABA)
echo "✓ Found GABA BED file: $BED_GABA"
echo "  Total GABA CREs: $N_GABA_CRES"
echo ""

# Check Excitatory BED file
if [ ! -f "$BED_EXCITATORY" ]; then
    echo "ERROR: Excitatory BED file not found: $BED_EXCITATORY"
    echo "Expected file created by 2a_extract_excitatory_neuron_CREs.py"
    echo ""
    echo "Please run first:"
    echo "  cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/Specific_CREs"
    echo "  conda activate bigwig"
    echo "  python 2a_extract_excitatory_neuron_CREs.py"
    exit 1
fi

N_EXCIT_CRES=$(wc -l < $BED_EXCITATORY)
echo "✓ Found Excitatory BED file: $BED_EXCITATORY"
echo "  Total Excitatory CREs: $N_EXCIT_CRES"
echo ""

# Copy BED files to output directory for reference
cp $BED_GABA $OUTPUT_DIR/GABA_CREs.bed
cp $BED_EXCITATORY $OUTPUT_DIR/Excitatory_CREs.bed
echo "✓ Copied BED files to output directory"
echo ""

# ============================================================================
# Step 2: Check BigWig files
# ============================================================================
echo "========================================================================"
echo "STEP 2: Checking BigWig files"
echo "------------------------------------------------------------------------"

# Collect BigWig files
# NOTE: Emx1-Ctrl is excluded (failed sample) - using only Nestin-Ctrl as control
BIGWIGS=""
LABELS=""
MISSING=0

# Define samples to include (excluding Emx1-Ctrl which is a failed sample)
SAMPLES_TO_INCLUDE="Nestin-Ctrl Nestin-Mut Emx1-Mut"

for SAMPLE in $SAMPLES_TO_INCLUDE; do
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
    echo "WARNING: Missing $MISSING BigWig files"
    echo "Will only analyze available conditions"
fi

echo ""

# ============================================================================
# Step 3: Run computeMatrix for GABA CREs (positive control)
# ============================================================================
echo "========================================================================"
echo "STEP 3: Computing signal matrices for GABA CREs (POSITIVE CONTROL)"
echo "------------------------------------------------------------------------"

echo "Running computeMatrix on GABA CREs..."
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
    -o $OUTPUT_DIR/matrix_GABA_all_conditions.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_GABA_all_conditions.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_GABA.log

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for GABA CREs"
    exit 1
fi

echo ""
echo "✓ GABA matrix computed and saved"
echo ""

# ============================================================================
# Step 4: Run computeMatrix for Excitatory CREs (negative control)
# ============================================================================
echo "========================================================================"
echo "STEP 4: Computing signal matrices for Excitatory CREs (NEGATIVE CONTROL)"
echo "------------------------------------------------------------------------"

echo "Running computeMatrix on Excitatory CREs..."
echo "  Window size: ±${WINDOW_SIZE} bp"
echo "  Bin size: ${BIN_SIZE} bp"
echo "  Processors: $N_PROCESSORS"
echo "  CREs: $N_EXCIT_CRES"
echo ""

computeMatrix reference-point \
    --referencePoint center \
    -b $WINDOW_SIZE -a $WINDOW_SIZE \
    -R $BED_EXCITATORY \
    -S $BIGWIGS \
    --samplesLabel $LABELS \
    --binSize $BIN_SIZE \
    --sortRegions keep \
    --missingDataAsZero \
    -o $OUTPUT_DIR/matrix_Excitatory_all_conditions.gz \
    -p $N_PROCESSORS \
    --outFileNameMatrix $OUTPUT_DIR/matrix_Excitatory_all_conditions.tab \
    2>&1 | tee $OUTPUT_DIR/computeMatrix_Excitatory.log

if [ $? -ne 0 ]; then
    echo "ERROR: computeMatrix failed for Excitatory CREs"
    exit 1
fi

echo ""
echo "✓ Excitatory matrix computed and saved"
echo ""

# ============================================================================
# Step 5: Determine common scale for comparison
# ============================================================================
echo "========================================================================"
echo "STEP 5: Determining common scale for GABA vs Excitatory comparison"
echo "------------------------------------------------------------------------"

echo "Calculating data ranges for both matrices..."

# Get data range from GABA matrix
computeMatrixOperations dataRange -m $OUTPUT_DIR/matrix_GABA_all_conditions.gz 2>&1 | tee $OUTPUT_DIR/dataRange_GABA.txt
# Extract 90th percentile (6th column) and max (3rd column) - skip header and non-numeric lines
GABA_90TH=$(grep -E '^[A-Za-z0-9-]+\s+[0-9]' $OUTPUT_DIR/dataRange_GABA.txt | awk '{print $6}' | sort -n | tail -1)
GABA_MAX=$(grep -E '^[A-Za-z0-9-]+\s+[0-9]' $OUTPUT_DIR/dataRange_GABA.txt | awk '{print $3}' | sort -n | tail -1)

# Get data range from Excitatory matrix
computeMatrixOperations dataRange -m $OUTPUT_DIR/matrix_Excitatory_all_conditions.gz 2>&1 | tee $OUTPUT_DIR/dataRange_Excitatory.txt
# Extract 90th percentile (6th column) and max (3rd column)
EXCIT_90TH=$(grep -E '^[A-Za-z0-9-]+\s+[0-9]' $OUTPUT_DIR/dataRange_Excitatory.txt | awk '{print $6}' | sort -n | tail -1)
EXCIT_MAX=$(grep -E '^[A-Za-z0-9-]+\s+[0-9]' $OUTPUT_DIR/dataRange_Excitatory.txt | awk '{print $3}' | sort -n | tail -1)

# Use 90th percentile for HEATMAP scale (more robust than max which can be outlier)
# Add 20% buffer for visibility
COMMON_90TH=$(echo "$GABA_90TH $EXCIT_90TH" | awk '{print ($1 > $2) ? $1 : $2}')
COMMON_MAX=$(echo "$COMMON_90TH" | awk '{print $1 * 1.2}')

# Use LOWER SCALE for METAPROFILE (metaprofiles show averaged signal, much lower than 90th percentile)
# If 90th percentile is very low (<0.5), use 20% of max instead; otherwise use 20% of 90th percentile
if (( $(echo "$COMMON_90TH < 0.5" | bc -l) )); then
    COMMON_MAX_VAL=$(echo "$GABA_MAX $EXCIT_MAX" | awk '{print ($1 > $2) ? $1 : $2}')
    METAPROFILE_MAX=$(echo "$COMMON_MAX_VAL" | awk '{print $1 * 0.2}')
    echo "Note: Using 20% of max ($COMMON_MAX_VAL) for metaprofile scale (sparse data)"
else
    METAPROFILE_MAX=$(echo "$COMMON_90TH" | awk '{print $1 * 0.2}')
    echo "Note: Using 20% of 90th percentile for metaprofile scale"
fi

echo ""
echo "Data ranges:"
echo "  GABA:"
echo "    90th percentile: $GABA_90TH"
echo "    Absolute max: $GABA_MAX"
echo "  Excitatory:"
echo "    90th percentile: $EXCIT_90TH"
echo "    Absolute max: $EXCIT_MAX"
echo ""
echo "Scale settings:"
echo "  Heatmaps: 90th %ile + 20% buffer = 0 to $COMMON_MAX"
echo "    (Shows individual CRE intensities, robust against outliers)"
echo "  Metaprofiles: 0 to $METAPROFILE_MAX"
echo "    (Adaptive scale optimized for visibility of mean signal curves)"
echo ""

echo "✓ Common scale determined for comparison"
echo ""

# ============================================================================
# Step 6: Create individual heatmaps with COMMON SCALE
# ============================================================================
echo "========================================================================"
echo "STEP 6: Creating individual heatmaps with common scale"
echo "------------------------------------------------------------------------"

# GABA heatmap (with common scale)
echo "Creating GABA heatmap (scale: 0-$COMMON_MAX)..."
plotHeatmap \
    -m $OUTPUT_DIR/matrix_GABA_all_conditions.gz \
    -o $OUTPUT_DIR/heatmap_GABA_all_conditions.png \
    --colorMap Reds \
    --dpi 300 \
    --whatToShow 'heatmap and colorbar' \
    --zMin 0 \
    --zMax $COMMON_MAX \
    --refPointLabel "CRE Center" \
    --heatmapHeight 25 \
    --heatmapWidth 6 \
    --sortRegions descend \
    --sortUsing mean \
    --plotTitle "POSITIVE CONTROL: ATAC Signal at GABA CREs (n=$N_GABA_CRES)" \
    --xAxisLabel "Distance from CRE Center (bp)" \
    --yAxisLabel "CREs" \
    2>&1 | grep -v "^$"

echo "✓ Saved: heatmap_GABA_all_conditions.png"

# Excitatory heatmap (with SAME common scale)
echo "Creating Excitatory heatmap (scale: 0-$COMMON_MAX)..."
plotHeatmap \
    -m $OUTPUT_DIR/matrix_Excitatory_all_conditions.gz \
    -o $OUTPUT_DIR/heatmap_Excitatory_all_conditions.png \
    --colorMap Reds \
    --dpi 300 \
    --whatToShow 'heatmap and colorbar' \
    --zMin 0 \
    --zMax $COMMON_MAX \
    --refPointLabel "CRE Center" \
    --heatmapHeight 25 \
    --heatmapWidth 6 \
    --sortRegions descend \
    --sortUsing mean \
    --plotTitle "NEGATIVE CONTROL: ATAC Signal at Excitatory CREs (n=$N_EXCIT_CRES)" \
    --xAxisLabel "Distance from CRE Center (bp)" \
    --yAxisLabel "CREs" \
    2>&1 | grep -v "^$"

echo "✓ Saved: heatmap_Excitatory_all_conditions.png"
echo ""
echo "NOTE: Both heatmaps use IDENTICAL scale (0-$COMMON_MAX) for direct comparison"
echo ""

# ============================================================================
# Step 7: Create metaprofiles with COMMON SCALE
# ============================================================================
echo "========================================================================"
echo "STEP 7: Creating metaprofiles with common scale"
echo "------------------------------------------------------------------------"

# Individual metaprofiles with common Y-axis scale
echo "Creating individual metaprofiles with common scale..."

# GABA metaprofile
# NOTE: 3 samples (Emx1-Ctrl excluded)
echo "Creating GABA metaprofile (Y-axis: auto-scaled)..."
plotProfile \
    -m $OUTPUT_DIR/matrix_GABA_all_conditions.gz \
    -o $OUTPUT_DIR/metaprofile_GABA_all_conditions.png \
    --plotTitle "POSITIVE CONTROL: ATAC Signal at GABA CREs (n=$N_GABA_CRES)" \
    --refPointLabel "CRE Center" \
    --averageType mean \
    --plotHeight 7 \
    --plotWidth 10 \
    --colors '#2E86AB' '#A23B72' '#C73E1D' \
    --yAxisLabel "Mean ATAC Signal" \
    --yMin 0 \
    2>&1 | grep -v "^$"

echo "✓ Saved: metaprofile_GABA_all_conditions.png"

# Excitatory metaprofile (auto-scaled independently)
# NOTE: 3 samples (Emx1-Ctrl excluded)
echo "Creating Excitatory metaprofile (Y-axis: auto-scaled)..."
plotProfile \
    -m $OUTPUT_DIR/matrix_Excitatory_all_conditions.gz \
    -o $OUTPUT_DIR/metaprofile_Excitatory_all_conditions.png \
    --plotTitle "NEGATIVE CONTROL: ATAC Signal at Excitatory CREs (n=$N_EXCIT_CRES)" \
    --refPointLabel "CRE Center" \
    --averageType mean \
    --plotHeight 7 \
    --plotWidth 10 \
    --colors '#2E86AB' '#A23B72' '#C73E1D' \
    --yAxisLabel "Mean ATAC Signal" \
    --yMin 0 \
    2>&1 | grep -v "^$"

echo "✓ Saved: metaprofile_Excitatory_all_conditions.png"
echo ""
echo "NOTE: Both metaprofiles auto-scale independently to their data range"
echo "      This ensures optimal visibility of signal curves in each plot"
echo ""

# ============================================================================
# Step 8: Create genotype-specific analyses (GABA CREs only)
# ============================================================================
echo "========================================================================"
echo "STEP 8: Creating genotype-specific analyses"
echo "------------------------------------------------------------------------"

# Nestin only (GABA CREs)
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
        --heatmapHeight 25 \
        --heatmapWidth 5 \
        --plotTitle "Nestin: ATAC Signal at GABA CREs (n=$N_GABA_CRES)" \
        --xAxisLabel "Distance from CRE Center (bp)" \
        2>&1 | grep -v "^$"

    echo "✓ Saved: heatmap_GABA_nestin.png"

    plotProfile \
        -m $OUTPUT_DIR/matrix_GABA_nestin.gz \
        -o $OUTPUT_DIR/metaprofile_GABA_nestin.png \
        --plotTitle "Nestin: Metaprofile at GABA CREs (n=$N_GABA_CRES)" \
        --refPointLabel "CRE Center" \
        --colors '#2E86AB' '#A23B72' \
        --plotHeight 6 \
        --plotWidth 8 \
        --yAxisLabel "Mean ATAC Signal" \
        --yMin 0 \
        2>&1 | grep -v "^$"

    echo "✓ Saved: metaprofile_GABA_nestin.png"
    echo ""
fi

# Emx1 only (GABA CREs)
# NOTE: Emx1-Ctrl is excluded (failed sample) - only Emx1-Mut is analyzed
if [ -f "$BIGWIG_BASE/GABA_Emx1-Mut.bw" ]; then
    echo "Creating Emx1 analysis (GABA CREs) - Emx1-Mut only (Emx1-Ctrl excluded)..."
    EMX1_BIGWIGS="$BIGWIG_BASE/GABA_Emx1-Mut.bw"
    EMX1_LABELS="Emx1-Mut"

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
        --heatmapHeight 25 \
        --heatmapWidth 5 \
        --plotTitle "Emx1-Mut: ATAC Signal at GABA CREs (n=$N_GABA_CRES)" \
        --xAxisLabel "Distance from CRE Center (bp)" \
        2>&1 | grep -v "^$"

    echo "✓ Saved: heatmap_GABA_emx1.png"

    plotProfile \
        -m $OUTPUT_DIR/matrix_GABA_emx1.gz \
        -o $OUTPUT_DIR/metaprofile_GABA_emx1.png \
        --plotTitle "Emx1-Mut: Metaprofile at GABA CREs (n=$N_GABA_CRES)" \
        --refPointLabel "CRE Center" \
        --colors '#C73E1D' \
        --plotHeight 6 \
        --plotWidth 8 \
        --yAxisLabel "Mean ATAC Signal" \
        --yMin 0 \
        2>&1 | grep -v "^$"

    echo "✓ Saved: metaprofile_GABA_emx1.png"
    echo ""
fi

# ============================================================================
# Step 9: Create summary report
# ============================================================================
echo "========================================================================"
echo "STEP 9: Creating summary report"
echo "------------------------------------------------------------------------"

cat > $OUTPUT_DIR/README.txt << EOFR
================================================================================
DEEPTOOLS HEATMAPS AND METAPROFILES: CELL-TYPE SPECIFICITY ANALYSIS
================================================================================

Analysis date: $(date)

PURPOSE:
--------
Demonstrate cell-type specificity of ATAC-seq signal in GABA samples by
comparing accessibility at:
- GABA-specific CREs (POSITIVE CONTROL - expect HIGH signal)
- Excitatory neuron CREs (NEGATIVE CONTROL - expect LOW signal)

CREs ANALYZED:
--------------
GABA CREs: $N_GABA_CRES (hippocampal interneurons)
Excitatory CREs: $N_EXCIT_CRES (glutamatergic neurons)

SOURCE DATA:
-----------
GABA BED file: output/hippocampal_interneuron_CREs.bed
  Created by: 1_extract_hippocampal_interneuron_CREs.py

Excitatory BED file: output/excitatory_neuron_CREs.bed
  Created by: 2a_extract_excitatory_neuron_CREs.py

BigWig files:
  ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw
  Created by: signac_06_create_bigwig_tracks.R

PARAMETERS:
-----------
Window size: ±${WINDOW_SIZE} bp around CRE center
Bin size: ${BIN_SIZE} bp
Reference point: CRE center
Sorting: By mean signal (descending)

CONDITIONS:
-----------
- Nestin-Ctrl
- Nestin-Mut
- Emx1-Ctrl
- Emx1-Mut

OUTPUT FILES:
-------------
BED files (copies for reference):
  - GABA_CREs.bed
  - Excitatory_CREs.bed

Signal matrices (compressed, reusable):
  GABA CREs:
    - matrix_GABA_all_conditions.gz (all 4 conditions)
    - matrix_GABA_nestin.gz (Nestin only)
    - matrix_GABA_emx1.gz (Emx1 only)

  Excitatory CREs:
    - matrix_Excitatory_all_conditions.gz
    - matrix_Excitatory_nestin.gz
    - matrix_Excitatory_emx1.gz

★ KEY COMPARISON PLOTS (use these for publication):
  - heatmap_GABA_all_conditions.png (POSITIVE CONTROL)
  - heatmap_Excitatory_all_conditions.png (NEGATIVE CONTROL)
  - metaprofile_GABA_all_conditions.png (POSITIVE CONTROL)
  - metaprofile_Excitatory_all_conditions.png (NEGATIVE CONTROL)

  NOTE: All plots use IDENTICAL scales for direct visual comparison

Individual heatmaps:
  - heatmap_GABA_all_conditions.png
  - heatmap_Excitatory_all_conditions.png
  - heatmap_GABA_nestin.png / heatmap_GABA_emx1.png

Individual metaprofiles:
  - metaprofile_GABA_all_conditions.png
  - metaprofile_Excitatory_all_conditions.png
  - metaprofile_GABA_nestin.png / metaprofile_GABA_emx1.png

Logs:
  - computeMatrix_GABA.log / computeMatrix_Excitatory.log
  - plotHeatmap_*.log
  - plotProfile_*.log

INTERPRETATION:
---------------
POSITIVE CONTROL (GABA CREs):
  - Should show HIGH ATAC signal in GABA samples
  - Signal concentrated at CRE centers
  - Validates that GABA samples capture GABAergic chromatin accessibility

NEGATIVE CONTROL (Excitatory CREs):
  - Should show LOW ATAC signal in GABA samples
  - Demonstrates cell-type specificity
  - If signal is high here, suggests contamination or non-specificity

Heatmaps:
  - Each row = one CRE
  - Color intensity = ATAC signal strength
  - CREs sorted by mean signal (highest at top)

Metaprofiles:
  - Average signal across all CREs
  - Line = mean, shaded area = SEM
  - Direct comparison of signal levels

EXPECTED RESULTS:
-----------------
1. GABA CREs should show prominent signal peaks at centers
2. Excitatory CREs should show minimal/no signal
3. Fold difference between GABA and Excitatory CREs validates specificity
4. Both Ctrl and Mut should show similar patterns (cell-type preserved)

QUALITY METRICS:
----------------
- ~7% of CREs may show zero signal (normal - not all CREs active simultaneously)
- "Skipping" warnings in logs are EXPECTED and NOT errors
- Check metaprofile_GABA_vs_Excitatory_*.png for clear separation

VISUALIZATION APPROACH:
-----------------------
Due to deepTools limitations, we create SEPARATE plots with IDENTICAL scales
instead of side-by-side combined plots. This ensures:
- Direct visual comparison (same color scale)
- Clear demonstration of signal differences
- Easy side-by-side placement in figures

NEXT STEPS:
-----------
1. Compare heatmaps side-by-side:
   → heatmap_GABA_all_conditions.png (expect STRONG red signal)
   → heatmap_Excitatory_all_conditions.png (expect WEAK signal)
   → Both use SAME color scale for fair comparison

2. Compare metaprofiles side-by-side:
   → metaprofile_GABA_all_conditions.png (expect HIGH curves)
   → metaprofile_Excitatory_all_conditions.png (expect LOW curves)
   → Both use SAME Y-axis scale for fair comparison

3. Calculate fold-enrichment (GABA/Excitatory signal ratio)
   → Use matrix files for quantitative analysis
   → Report mean signal difference

4. Compare Ctrl vs Mut within each genotype
   → Check if mutations alter accessibility patterns

5. For publication:
   → Place GABA and Excitatory plots side-by-side in figure panel
   → Highlight identical scales in figure legend
   → Use to demonstrate cell-type specificity

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
echo "  - GABA CREs (POSITIVE): $N_GABA_CRES"
echo "  - Excitatory CREs (NEGATIVE): $N_EXCIT_CRES"
echo ""
echo "★ KEY COMPARISON PLOTS (same scale for direct comparison):"
echo "  Heatmaps:"
echo "    - heatmap_GABA_all_conditions.png (POSITIVE CONTROL)"
echo "    - heatmap_Excitatory_all_conditions.png (NEGATIVE CONTROL)"
echo "  Metaprofiles:"
echo "    - metaprofile_GABA_all_conditions.png (POSITIVE CONTROL)"
echo "    - metaprofile_Excitatory_all_conditions.png (NEGATIVE CONTROL)"
echo ""
echo "  Common scale: 0 to $COMMON_MAX"
echo ""
echo "Genotype-specific (GABA only):"
echo "  - heatmap_GABA_nestin.png / metaprofile_GABA_nestin.png"
echo "  - heatmap_GABA_emx1.png / metaprofile_GABA_emx1.png"
echo ""
echo "EXPECTED RESULT:"
echo "  HIGH signal at GABA CREs (cell-type specific)"
echo "  LOW signal at Excitatory CREs (negative control)"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
