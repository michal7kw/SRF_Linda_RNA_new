#!/bin/bash
#SBATCH --job-name=4_compare_GABA_vs_Excitatory
#SBATCH --output=logs/4_compare_GABA_vs_Excitatory.log
#SBATCH --error=logs/4_compare_GABA_vs_Excitatory.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Create Heatmaps and Metaprofiles: Cell-Type Specificity Analysis
#
# This script runs the Python analysis comparing ATAC signal at:
# - GABA-specific CREs (POSITIVE CONTROL - expect HIGH signal)
# - Excitatory neuron CREs (NEGATIVE CONTROL - expect LOW signal)
#
# DEPENDENCIES:
# - Python script called: compare_GABA_vs_Excitatory.py (main analysis script)
# - Conda environment: rna_seq_analysis_deep (contains pyBigWig and dependencies)
# - Prerequisite scripts that generate input BED files:
#   * 1_extract_hippocampal_interneuron_CREs.py (creates GABA CREs)
#   * 2a_extract_excitatory_neuron_CREs.py (creates Excitatory CREs)
# - BigWig files from Signac pipeline: ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw
#
# INPUT FILES:
# - output/hippocampal_interneuron_CREs.bed: BED file containing GABA neuron CRE coordinates
# - output/excitatory_neuron_CREs.bed: BED file containing excitatory neuron CRE coordinates
# - ../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw: BigWig files with ATAC-seq signal data
#
# OUTPUT FILES:
# - Heatmaps with identical scales for both CRE types
# - Metaprofiles with identical Y-axes
# - Comparison plots and quantitative statistics
# - logs/create_heatmaps_metaprofiles.log: Execution log file
# - logs/create_heatmaps_metaprofiles.err: Error log file
#
################################################################################

echo "========================================================================"
echo "CREATE HEATMAPS AND METAPROFILES: GABA vs Excitatory"
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

# Check if required BED files exist
echo "Checking prerequisites..."
echo ""

if [ ! -f "output/hippocampal_interneuron_CREs.bed" ]; then
    echo "ERROR: GABA CREs not found!"
    echo "Please run: python 1_extract_hippocampal_interneuron_CREs.py"
    exit 1
fi

if [ ! -f "output/excitatory_neuron_CREs.bed" ]; then
    echo "ERROR: Excitatory CREs not found!"
    echo "Please run: python 2a_extract_excitatory_neuron_CREs.py"
    exit 1
fi

N_GABA=$(wc -l < output/hippocampal_interneuron_CREs.bed)
N_EXCIT=$(wc -l < output/excitatory_neuron_CREs.bed)

echo "✓ GABA CREs: $N_GABA"
echo "✓ Excitatory CREs: $N_EXCIT"
echo ""

# Check if BigWig files exist
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
    exit 1
fi

echo ""
echo "========================================================================"
echo "Running Python analysis..."
echo "========================================================================"
echo ""

# Run the Python script
python 4_compare_GABA_vs_Excitatory.py

EXIT_CODE=$?

echo ""
echo "========================================================================"

if [ $EXIT_CODE -eq 0 ]; then
    echo "ANALYSIS COMPLETE!"
    echo "========================================================================"
    echo ""
    echo "Output directory: output/GABA_DEG_analysis/heatmaps_metaprofiles/"
    echo ""
    echo "Key files created:"
    echo "  ★ Heatmaps (same color scale):"
    echo "      - heatmap_GABA_all_conditions.png (POSITIVE)"
    echo "      - heatmap_Excitatory_all_conditions.png (NEGATIVE)"
    echo ""
    echo "  ★ Metaprofiles (same Y-axis):"
    echo "      - metaprofile_GABA_all_CREs.png (POSITIVE)"
    echo "      - metaprofile_Excitatory_all_CREs.png (NEGATIVE)"
    echo ""
    echo "  ★ Comparison:"
    echo "      - comparison_GABA_vs_Excitatory.png"
    echo "      - enrichment_statistics.txt"
    echo ""
    echo "Completed: $(date)"
    echo "========================================================================"
else
    echo "ANALYSIS FAILED!"
    echo "========================================================================"
    echo ""
    echo "Check the error log for details:"
    echo "  logs/create_heatmaps_metaprofiles_${SLURM_JOB_ID}.err"
    echo ""
    exit 1
fi
