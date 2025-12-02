#!/bin/bash
#SBATCH --job-name=1_RUN_EXCLUSIVE_CRES_ANALYSIS
#SBATCH --output=logs/1_RUN_EXCLUSIVE_CRES_ANALYSIS.log
#SBATCH --error=logs/1_RUN_EXCLUSIVE_CRES_ANALYSIS.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# COMPLETE CELL-TYPE SPECIFICITY ANALYSIS (Step 1 of Pipeline)
#
# This script runs the first part of the Exclusive CREs pipeline:
# 1. Extract mutually exclusive CRE sets (GABA-specific vs Excitatory-specific)
# 2. Create deepTools heatmaps and metaprofiles with direct comparison
#
# SCRIPTS IN THIS FOLDER:
# ────────────────────────────────────────────────────────────────────────────
# Master Scripts:
#   0_RUN_COMPLETE_EXCLUSIVE_PIPELINE.sh  - Submits all steps with dependencies
#   1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh      - This script (Step 1)
#
# Step 1 - Extract & Visualize (this script):
#   1a_extract_cell_type_specific_CREs.py/.sh  - Extract mutually exclusive CREs
#   1b_create_heatmaps_specific_CREs_deeptools.sh - deepTools heatmaps
#
# Step 2 - Fold Enrichment:
#   2_fold_enrichment_exclusive_analysis.py/.sh - Calculate fold enrichment
#
# Step 3 - Comparison & Gene Linkage:
#   3_compare_cell_type_specific_CREs.py/.sh   - pyBigWig comparison plots
#   3_link_CREs_to_genes.py/.sh                - Link CREs to genes (Table 16)
#
# Step 4 - DA Visualization:
#   4_visualize_DA_CREs.py/.sh                 - Visualize DA CREs with filtering
#
# Utilities:
#   REPLOT_ONLY.sh                             - Re-run plotting only
#   reproduce_scaling.py                        - Reproduce scaling parameters
# ────────────────────────────────────────────────────────────────────────────
#
# Key improvement over previous approach:
# - Previous: Used CRE sets with 60% overlap → identical-looking plots
# - Current: Uses mutually exclusive CRE sets → clear visual difference
#
# INPUT FILES:
# - ../data/table_1.xlsx: ENCODE data table 1
# - ../data/table_2.tsv: ENCODE data table 2
# - ../data/table_8.txt: ENCODE data table 8
# - 1a_extract_cell_type_specific_CREs.py: Python script to extract mutually exclusive CRE sets
# - 1b_create_heatmaps_specific_CREs_deeptools.sh: Script to generate heatmaps for specific CREs
# - ../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw: BigWig files with ATAC-seq signal data
#
# OUTPUT FILES:
# - output/GABA_specific_CREs.bed: BED file with GABA-specific CREs (mutually exclusive)
# - output/Excitatory_specific_CREs.bed: BED file with excitatory-specific CREs (mutually exclusive)
# - output/cell_type_specific_CREs_summary.txt: Summary report of CRE extraction
# - output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/heatmap_GABA_specific.png: Heatmap for GABA-specific CREs
# - output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/heatmap_Excitatory_specific.png: Heatmap for excitatory-specific CREs
# - output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/metaprofile_GABA_specific.png: Metaprofile for GABA-specific CREs
# - output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/metaprofile_Excitatory_specific.png: Metaprofile for excitatory-specific CREs
# - output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/comparison_cell_type_specific.png: Comparison plot
# - output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/cell_type_specificity_statistics.txt: Statistics file
# - logs/specificity_analysis.log: Execution log file
# - logs/specificity_analysis.err: Error log file
#
################################################################################

echo "========================================================================"
echo "COMPLETE CELL-TYPE SPECIFICITY ANALYSIS"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Create directories
mkdir -p logs
mkdir -p output

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# ============================================================================
# STEP 1: Extract Cell-Type-Specific CREs
# ============================================================================
echo "========================================================================"
echo "STEP 1/2: Extracting mutually exclusive CRE sets"
echo "========================================================================"
echo ""

if [ -f "output/GABA_specific_CREs.bed" ] && [ -f "output/Excitatory_specific_CREs.bed" ]; then
    echo "Cell-type-specific CRE files already exist:"
    echo "  - output/GABA_specific_CREs.bed ($(wc -l < output/GABA_specific_CREs.bed) CREs)"
    echo "  - output/Excitatory_specific_CREs.bed ($(wc -l < output/Excitatory_specific_CREs.bed) CREs)"
    echo ""
    echo "Using existing CRE files (set REEXTRACT=1 to force re-extraction)..."
    echo ""

    # Check if user wants to force re-extraction via environment variable
    if [ "${REEXTRACT}" = "1" ]; then
        echo "REEXTRACT=1 set, re-extracting CREs..."
        python 1a_extract_cell_type_specific_CREs.py

        if [ $? -ne 0 ]; then
            echo "ERROR: CRE extraction failed!"
            exit 1
        fi
    fi
else
    echo "Extracting cell-type-specific CREs..."
    python 1a_extract_cell_type_specific_CREs.py

    if [ $? -ne 0 ]; then
        echo "ERROR: CRE extraction failed!"
        exit 1
    fi
fi

echo ""

# Verify mutual exclusivity
echo "Verifying mutual exclusivity..."
OVERLAP=$(comm -12 \
    <(awk '{print $1":"$2"-"$3}' output/GABA_specific_CREs.bed | sort) \
    <(awk '{print $1":"$2"-"$3}' output/Excitatory_specific_CREs.bed | sort) \
    | wc -l)

if [ $OVERLAP -eq 0 ]; then
    echo "  ✓ VERIFIED: 0 overlapping coordinates"
    N_GABA=$(wc -l < output/GABA_specific_CREs.bed)
    N_EXCIT=$(wc -l < output/Excitatory_specific_CREs.bed)
    echo "  ✓ GABA-specific CREs: $N_GABA"
    echo "  ✓ Excitatory-specific CREs: $N_EXCIT"
else
    echo "  ✗ ERROR: $OVERLAP overlapping coordinates found!"
    echo "     CRE sets are not mutually exclusive!"
    exit 1
fi

echo ""

# ============================================================================
# STEP 2: Create Heatmaps and Metaprofiles (using fast deepTools)
# ============================================================================
echo "========================================================================"
echo "STEP 2/2: Creating heatmaps and metaprofiles (using deepTools)"
echo "========================================================================"
echo ""

bash 1b_create_heatmaps_specific_CREs_deeptools.sh

if [ $? -ne 0 ]; then
    echo "ERROR: Plotting failed!"
    exit 1
fi

echo ""

# ============================================================================
# FINAL SUMMARY
# ============================================================================
echo "========================================================================"
echo "ANALYSIS COMPLETE!"
echo "========================================================================"
echo ""

OUTPUT_DIR="output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools"

echo "CRE SETS (MUTUALLY EXCLUSIVE):"
echo "  - GABA-specific: $N_GABA CREs"
echo "  - Excitatory-specific: $N_EXCIT CREs"
echo "  - Overlap: 0 (verified)"
echo ""

echo "OUTPUT DIRECTORY:"
echo "  $OUTPUT_DIR/"
echo ""

echo "KEY FILES:"
echo ""
echo "  📊 Heatmaps (identical color scale):"
echo "      ├─ heatmap_GABA_specific.png"
echo "      │   → POSITIVE CONTROL: Expect STRONG red signal"
echo "      └─ heatmap_Excitatory_specific.png"
echo "          → NEGATIVE CONTROL: Expect PALE/white signal"
echo ""
echo "  📈 Metaprofiles (identical Y-axis):"
echo "      ├─ metaprofile_GABA_specific.png"
echo "      │   → POSITIVE CONTROL: Expect HIGH curves"
echo "      └─ metaprofile_Excitatory_specific.png"
echo "          → NEGATIVE CONTROL: Expect LOW/flat curves"
echo ""
echo "  📉 Comparison:"
echo "      ├─ comparison_cell_type_specific.png"
echo "      │   → Side-by-side comparison with fold enrichment"
echo "      └─ cell_type_specificity_statistics.txt"
echo "          → Quantitative metrics and interpretation"
echo ""

echo "========================================================================"
echo "EXPECTED RESULTS"
echo "========================================================================"
echo ""
echo "If analysis worked correctly, you should see:"
echo ""
echo "  ✓ GABA-specific heatmap: Strong red signal throughout"
echo "  ✓ Excitatory-specific heatmap: Pale/white (minimal signal)"
echo "  ✓ CLEAR VISUAL DIFFERENCE between the two"
echo "  ✓ Fold enrichment ≥3x in statistics file"
echo ""
echo "This demonstrates that GABA ATAC samples specifically capture"
echo "GABAergic chromatin accessibility!"
echo ""

echo "========================================================================"
echo "QUICK CHECK"
echo "========================================================================"
echo ""

if [ -f "$OUTPUT_DIR/cell_type_specificity_statistics.txt" ]; then
    echo "Fold enrichment values:"
    grep "Fold enrichment:" "$OUTPUT_DIR/cell_type_specificity_statistics.txt" | head -4
    echo ""
fi

echo "========================================================================"
echo "NEXT STEPS"
echo "========================================================================"
echo ""
echo "1. Review statistics file:"
echo "   cat $OUTPUT_DIR/cell_type_specificity_statistics.txt"
echo ""
echo "2. View plots in output directory:"
echo "   ls -lh $OUTPUT_DIR/*.png"
echo ""
echo "3. Compare with previous analysis (60% overlap):"
echo "   - Old: output/GABA_DEG_analysis/heatmaps_metaprofiles/"
echo "   - New: $OUTPUT_DIR/"
echo "   → New version should show CLEAR difference!"
echo ""
echo "4. For publication:"
echo "   → Use heatmap_GABA_specific.png vs heatmap_Excitatory_specific.png"
echo "   → Highlight mutually exclusive CRE sets in methods"
echo "   → Report fold enrichment values from statistics file"
echo ""

echo "Completed: $(date)"
echo "========================================================================"
