#!/bin/bash
#SBATCH --job-name=3_compare_cell_type_specific_CREs
#SBATCH --output=logs/3_compare_cell_type_specific_CREs.log
#SBATCH --error=logs/3_compare_cell_type_specific_CREs.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=workq

################################################################################
# Create Heatmaps and Metaprofiles: Cell-Type-Specific CREs
#
# This script analyzes ATAC signal at MUTUALLY EXCLUSIVE CRE sets:
# - GABA-specific CREs (POSITIVE CONTROL - expect HIGH signal)
# - Excitatory-specific CREs (NEGATIVE CONTROL - expect LOW signal)
#
# Key difference from previous analysis:
# - Previous: 60% overlap between CRE sets → identical-looking plots
# - Current: 0% overlap (mutually exclusive) → clear visual difference
#
# DEPENDENCIES:
# - Python script called: 3_compare_cell_type_specific_CREs.py (main analysis script)
# - Prerequisite script: 1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh (creates mutually exclusive CREs)
# - Conda environment: bigwig (contains pyBigWig and dependencies)
# - BigWig files from Signac pipeline: ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw
# - Signac pipeline script: ../signac_results_L1/run_signac_step6_bigwig_L1.sh (generates BigWig files)
#
# INPUT FILES:
# - output/GABA_specific_CREs.bed: BED file with GABA-specific CREs (mutually exclusive)
# - output/Excitatory_specific_CREs.bed: BED file with excitatory-specific CREs (mutually exclusive)
# - 3_compare_cell_type_specific_CREs.py: Python script to generate heatmaps and metaprofiles
# - ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw: BigWig files with ATAC-seq signal data
#
# OUTPUT FILES:
# - output/GABA_DEG_analysis/heatmaps_specific_CREs/heatmap_GABA_specific.png: Heatmap for GABA-specific CREs
# - output/GABA_DEG_analysis/heatmaps_specific_CREs/heatmap_Excitatory_specific.png: Heatmap for excitatory-specific CREs
# - output/GABA_DEG_analysis/heatmaps_specific_CREs/metaprofile_GABA_specific.png: Metaprofile for GABA-specific CREs
# - output/GABA_DEG_analysis/heatmaps_specific_CREs/metaprofile_Excitatory_specific.png: Metaprofile for excitatory-specific CREs
# - output/GABA_DEG_analysis/heatmaps_specific_CREs/comparison_cell_type_specific.png: Comparison plot
# - output/GABA_DEG_analysis/heatmaps_specific_CREs/cell_type_specificity_statistics.txt: Statistics file
# - logs/plot_specific_CREs.log: Execution log file
# - logs/plot_specific_CREs.err: Error log file
#
################################################################################

echo "========================================================================"
echo "CREATE HEATMAPS AND METAPROFILES: CELL-TYPE-SPECIFIC CREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Create logs directory
mkdir -p logs

# Activate conda environment with pyBigWig
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# ============================================================================
# Check prerequisites
# ============================================================================
echo "Checking prerequisites..."
echo ""

MISSING=0

# Check BED files
if [ ! -f "output/GABA_specific_CREs.bed" ]; then
    echo "✗ Missing: output/GABA_specific_CREs.bed"
    echo "  Please run: sbatch 1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh"
    MISSING=$((MISSING + 1))
else
    N_GABA=$(wc -l < output/GABA_specific_CREs.bed)
    echo "✓ Found GABA-specific CREs: $N_GABA"
fi

if [ ! -f "output/Excitatory_specific_CREs.bed" ]; then
    echo "✗ Missing: output/Excitatory_specific_CREs.bed"
    echo "  Please run: bash RUN_EXTRACT_SPECIFIC_CRES.sh"
    MISSING=$((MISSING + 1))
else
    N_EXCIT=$(wc -l < output/Excitatory_specific_CREs.bed)
    echo "✓ Found Excitatory-specific CREs: $N_EXCIT"
fi

echo ""

# Verify mutual exclusivity
if [ -f "output/GABA_specific_CREs.bed" ] && [ -f "output/Excitatory_specific_CREs.bed" ]; then
    echo "Verifying mutual exclusivity..."
    OVERLAP=$(comm -12 \
        <(awk '{print $1":"$2"-"$3}' output/GABA_specific_CREs.bed | sort) \
        <(awk '{print $1":"$2"-"$3}' output/Excitatory_specific_CREs.bed | sort) \
        | wc -l)

    if [ $OVERLAP -eq 0 ]; then
        echo "  ✓ VERIFIED: 0 overlapping coordinates (mutually exclusive)"
    else
        echo "  ✗ WARNING: $OVERLAP overlapping coordinates found!"
        echo "     Re-run: sbatch 1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh"
        MISSING=$((MISSING + 1))
    fi
    echo ""
fi

# Check BigWig files
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"

echo "Checking BigWig files..."
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

echo ""

if [ $MISSING -gt 0 ]; then
    echo "ERROR: Missing $MISSING required files!"
    echo ""
    echo "To fix:"
    echo "  1. Extract cell-type-specific CREs:"
    echo "     sbatch 1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh"
    echo ""
    echo "  2. Generate BigWig tracks (if missing):"
    echo "     cd ../signac_results_L1/"
    echo "     sbatch run_signac_step6_bigwig_L1.sh"
    exit 1
fi

# ============================================================================
# Run analysis
# ============================================================================
echo "========================================================================"
echo "Running Python analysis..."
echo "========================================================================"
echo ""

python 3_compare_cell_type_specific_CREs.py

EXIT_CODE=$?

echo ""
echo "========================================================================"

if [ $EXIT_CODE -eq 0 ]; then
    echo "ANALYSIS COMPLETE!"
    echo "========================================================================"
    echo ""
    echo "Output directory: output/GABA_DEG_analysis/heatmaps_specific_CREs/"
    echo ""
    echo "Key files created:"
    echo ""
    echo "  ★ Heatmaps (same color scale):"
    echo "      - heatmap_GABA_specific.png (POSITIVE CONTROL)"
    echo "      - heatmap_Excitatory_specific.png (NEGATIVE CONTROL)"
    echo ""
    echo "  ★ Metaprofiles (same Y-axis):"
    echo "      - metaprofile_GABA_specific.png (POSITIVE CONTROL)"
    echo "      - metaprofile_Excitatory_specific.png (NEGATIVE CONTROL)"
    echo ""
    echo "  ★ Comparison and statistics:"
    echo "      - comparison_cell_type_specific.png"
    echo "      - cell_type_specificity_statistics.txt"
    echo ""
    echo "========================================================================"
    echo "EXPECTED RESULTS"
    echo "========================================================================"
    echo ""
    echo "✓ GABA-specific CREs should show:"
    echo "    - Strong red signal in heatmaps"
    echo "    - High curves in metaprofiles"
    echo "    - Signal concentrated at CRE centers"
    echo ""
    echo "✓ Excitatory-specific CREs should show:"
    echo "    - Pale/white signal in heatmaps (LOW signal)"
    echo "    - Low/flat curves in metaprofiles"
    echo "    - Minimal signal throughout"
    echo ""
    echo "✓ Fold enrichment should be ≥3x"
    echo "    (demonstrates strong cell-type specificity)"
    echo ""
    echo "KEY IMPROVEMENT:"
    echo "  Previous analysis: 60% overlap → identical plots"
    echo "  Current analysis: 0% overlap → CLEAR DIFFERENCE!"
    echo ""
    echo "========================================================================"
    echo "NEXT STEPS"
    echo "========================================================================"
    echo ""
    echo "1. Review statistics:"
    echo "   cat output/GABA_DEG_analysis/heatmaps_specific_CREs/cell_type_specificity_statistics.txt"
    echo ""
    echo "2. View comparison plot:"
    echo "   Check: comparison_cell_type_specific.png"
    echo "   → Should show clear fold enrichment (≥3x)"
    echo ""
    echo "3. Compare heatmaps side-by-side:"
    echo "   - heatmap_GABA_specific.png (expect STRONG red)"
    echo "   - heatmap_Excitatory_specific.png (expect PALE/white)"
    echo ""
    echo "4. For publication:"
    echo "   → Place GABA and Excitatory heatmaps side-by-side in figure panel"
    echo "   → Highlight mutually exclusive CRE sets in methods"
    echo "   → Report fold enrichment as measure of specificity"
    echo ""
    echo "Completed: $(date)"
    echo "========================================================================"
else
    echo "ANALYSIS FAILED!"
    echo "========================================================================"
    echo ""
    echo "Check the error log for details:"
    echo "  logs/plot_specific_CREs_${SLURM_JOB_ID}.err"
    echo ""
    exit 1
fi
