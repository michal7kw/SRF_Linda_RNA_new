#!/bin/bash
#SBATCH --job-name=2_RUN_SPECIFIC_CRES_ANALYSIS
#SBATCH --output=logs/2_RUN_SPECIFIC_CRES_ANALYSIS.log
#SBATCH --error=logs/2_RUN_SPECIFIC_CRES_ANALYSIS.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Cell-Type Specificity Analysis Pipeline - OLD VERSION
#
# ⚠️  NOTE: This is the OLDER version using CRE sets with 60% overlap
#
# Purpose: Demonstrate ATAC signal enrichment in OVERLAPPING CREs
#          (CREs that are active in both GABA and Excitatory neurons)
#
# This script uses:
# - hippocampal_interneuron_CREs.bed (GABA CREs - 239,691 CREs)
# - excitatory_neuron_CREs.bed (Excitatory CREs - 237,835 CREs)
# - 60% OVERLAP (144,015 shared CREs)
#
# Expected result: Both heatmaps show signal because many CREs are shared
#
# For MUTUALLY EXCLUSIVE CREs (0% overlap), use instead:
#   sbatch ../Exclusive_CREs/1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh
#
# Steps:
# 1. Extract Excitatory neuron CREs (if not already done)
# 2. Run deepTools heatmap and metaprofile analysis (directly, not as separate job)
# 3. Generate comparison plots
#
# Usage:
#   sbatch 2_RUN_SPECIFIC_CRES_ANALYSIS.sh
#
# INPUT FILES:
# - output/hippocampal_interneuron_CREs.bed: BED file containing GABA neuron CRE coordinates
# - 2a_extract_excitatory_neuron_CREs.py: Python script to extract excitatory neuron CREs
# - 2b_create_heatmaps_deeptools.sh: Script to generate heatmaps and metaprofiles
# - ../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw: BigWig files with ATAC-seq signal data
#
# OUTPUT FILES:
# - output/excitatory_neuron_CREs.bed: BED file containing excitatory neuron CRE coordinates
# - output/GABA_DEG_analysis/heatmaps_deeptools/heatmap_GABA_vs_Excitatory_all_conditions.png: Combined heatmap comparison
# - output/GABA_DEG_analysis/heatmaps_deeptools/metaprofile_GABA_vs_Excitatory_all_conditions.png: Combined metaprofile comparison
# - logs/heatmaps_deeptools_${JOB_ID}.log: Execution log file
# - logs/heatmaps_deeptools_${JOB_ID}.err: Error log file
#
################################################################################

set -e  # Exit on error

echo "========================================================================"
echo "CELL-TYPE SPECIFICITY ANALYSIS PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Navigate to correct directory
# cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature

# ============================================================================
# Step 1: Check if GABA CREs exist
# ============================================================================
echo "========================================================================"
echo "STEP 1: Checking GABA CREs (positive control)"
echo "------------------------------------------------------------------------"

if [ ! -f "output/hippocampal_interneuron_CREs.bed" ]; then
    echo "ERROR: GABA CREs not found!"
    echo "Please run 1_extract_hippocampal_interneuron_CREs.sh first"
    exit 1
fi

N_GABA=$(wc -l < output/hippocampal_interneuron_CREs.bed)
echo "✓ Found GABA CREs: $N_GABA"
echo ""

# ============================================================================
# Step 2: Extract Excitatory neuron CREs (if needed)
# ============================================================================
echo "========================================================================"
echo "STEP 2: Checking Excitatory CREs (negative control)"
echo "------------------------------------------------------------------------"

if [ -f "output/excitatory_neuron_CREs.bed" ]; then
    N_EXCIT=$(wc -l < output/excitatory_neuron_CREs.bed)
    echo "✓ Excitatory CREs already exist: $N_EXCIT"
    echo "Skipping extraction..."
    echo ""
else
    echo "Excitatory CREs not found. Extracting now..."
    echo ""

    # Activate conda environment with pandas
    source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
    conda activate rna_seq_analysis_deep

    # Run extraction script
    echo "Running 2a_extract_excitatory_neuron_CREs.py..."
    python 2a_extract_excitatory_neuron_CREs.py

    if [ $? -ne 0 ]; then
        echo "ERROR: Excitatory CRE extraction failed!"
        exit 1
    fi

    if [ ! -f "output/excitatory_neuron_CREs.bed" ]; then
        echo "ERROR: Expected output file not created!"
        exit 1
    fi

    N_EXCIT=$(wc -l < output/excitatory_neuron_CREs.bed)
    echo ""
    echo "✓ Successfully extracted Excitatory CREs: $N_EXCIT"
    echo ""
fi

# ============================================================================
# Step 3: Verify BigWig files
# ============================================================================
echo "========================================================================"
echo "STEP 3: Checking BigWig files"
echo "------------------------------------------------------------------------"

BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
MISSING=0

for GENOTYPE in Nestin Emx1; do
    for CONDITION in Ctrl Mut; do
        BW_FILE="$BIGWIG_BASE/GABA_${GENOTYPE}-${CONDITION}.bw"
        if [ -f "$BW_FILE" ]; then
            echo "✓ Found: GABA_${GENOTYPE}-${CONDITION}.bw"
        else
            echo "✗ Missing: $BW_FILE"
            MISSING=$((MISSING + 1))
        fi
    done
done

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "ERROR: Missing $MISSING BigWig files!"
    echo "Please run 2b_create_heatmaps_deeptools.sh first"
    exit 1
fi

echo ""
echo "✓ All BigWig files found"
echo ""

# ============================================================================
# Step 4: Run heatmap analysis directly (logged to this job's output)
# ============================================================================
echo "========================================================================"
echo "STEP 4: Running heatmap analysis with deepTools"
echo "------------------------------------------------------------------------"
echo ""

if [ ! -f "2b_create_heatmaps_deeptools.sh" ]; then
    echo "ERROR: 2b_create_heatmaps_deeptools.sh not found!"
    exit 1
fi

# Run the heatmap script directly (not as separate SLURM job)
# All output will be captured in this job's log files
bash 2b_create_heatmaps_deeptools.sh

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Heatmap analysis failed!"
    echo "Check the error messages above for details"
    exit 1
fi

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "========================================================================"
echo "ANALYSIS COMPLETE!"
echo "========================================================================"
echo ""
echo "Analysis components:"
echo "  - GABA CREs: $N_GABA"
echo "  - Excitatory CREs: $N_EXCIT"
echo "  - CRE overlap: ~60% (144,015 shared CREs)"
echo "  - BigWig samples: 3 (Nestin-Ctrl, Nestin-Mut, Emx1-Mut) - Emx1-Ctrl excluded"
echo ""
echo "Output directory:"
echo "  output/GABA_DEG_analysis/heatmaps_deeptools/"
echo ""
echo "★ KEY OUTPUTS:"
echo "  Heatmaps:"
echo "    - heatmap_GABA_all_conditions.png"
echo "    - heatmap_Excitatory_all_conditions.png"
echo "    - heatmap_GABA_vs_Excitatory_all_conditions.png (side-by-side)"
echo ""
echo "  Metaprofiles:"
echo "    - metaprofile_GABA_all_CREs.png"
echo "    - metaprofile_Excitatory_all_CREs.png"
echo "    - metaprofile_GABA_vs_Excitatory_all_conditions.png (overlay)"
echo ""
echo "  Genotype-specific (Nestin):"
echo "    - heatmap_GABA_nestin.png"
echo "    - metaprofile_GABA_nestin.png"
echo ""
echo "  Genotype-specific (Emx1):"
echo "    - heatmap_GABA_emx1.png"
echo "    - metaprofile_GABA_emx1.png"
echo ""
echo "EXPECTED RESULT:"
echo "  ⚠️  Both GABA and Excitatory heatmaps should show SIMILAR signal"
echo "     because 60% of CREs are shared between cell types"
echo ""
echo "  This demonstrates:"
echo "    - ATAC signal is present and enriched at CREs ✓"
echo "    - Many regulatory elements are shared (biology, not bug) ✓"
echo ""
echo "  For cell-type SPECIFICITY, compare with:"
echo "    output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/"
echo "#    (from ../Exclusive_CREs/1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh)"
echo ""
echo "Logs for this analysis:"
echo "  logs/2_RUN_SPECIFIC_CRES_ANALYSIS_${SLURM_JOB_ID}.log"
echo "#  logs/2_RUN_SPECIFIC_CRES_ANALYSIS_${SLURM_JOB_ID}.err"
echo ""
echo "Completed: $(date)"
echo "========================================================================"

