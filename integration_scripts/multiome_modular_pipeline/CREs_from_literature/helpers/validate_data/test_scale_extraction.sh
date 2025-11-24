#!/bin/bash
#SBATCH --job-name=test_scale_debug
#SBATCH --output=logs/test_scale_debug.log
#SBATCH --error=logs/test_scale_debug.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=workq

################################################################################
# Test Script: Debug Identical Plots Issue
#
# This script runs the Python comparison with debug output to investigate
# why GABA and Excitatory plots appear identical.
#
# DEPENDENCIES:
# - Python script called: compare_GABA_vs_Excitatory.py (main analysis script with debug output)
# - Conda environment: bigwig (contains pyBigWig and dependencies)
# - BigWig files from Signac pipeline: ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw
# - Prerequisite scripts that generate input BED files:
#   * extract_hippocampal_interneuron_CREs.py (creates GABA CREs)
#   * extract_excitatory_neuron_CREs.py (creates Excitatory CREs)
#
# INPUT FILES:
# - output/hippocampal_interneuron_CREs.bed: BED file containing GABA neuron CRE coordinates
# - output/excitatory_neuron_CREs.bed: BED file containing excitatory neuron CRE coordinates
# - compare_GABA_vs_Excitatory.py: Python script with debug output
# - ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw: BigWig files with ATAC-seq signal data
#
# OUTPUT FILES:
# - Debug output showing different BED files loaded
# - Debug output showing different CRE coordinates
# - Debug output showing different signal matrices
# - Heatmap and metaprofile files (same as compare_GABA_vs_Excitatory.py)
# - logs/test_scale_debug.log: Execution log file
# - logs/test_scale_debug.err: Error log file
#
################################################################################

echo "========================================================================"
echo "TEST: Debug Identical Plots Issue"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Create logs directory
mkdir -p logs

# Activate conda environment with pyBigWig
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate bigwig

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# Check if required BED files exist
echo "Checking BED files..."
echo ""

if [ ! -f "output/hippocampal_interneuron_CREs.bed" ]; then
    echo "ERROR: GABA CREs not found!"
    exit 1
fi

if [ ! -f "output/excitatory_neuron_CREs.bed" ]; then
    echo "ERROR: Excitatory CREs not found!"
    exit 1
fi

echo "✓ GABA CREs: $(wc -l < output/hippocampal_interneuron_CREs.bed) lines"
echo "✓ Excitatory CREs: $(wc -l < output/excitatory_neuron_CREs.bed) lines"
echo ""

# Show first few lines of each BED file
echo "First 3 lines of GABA CREs:"
head -3 output/hippocampal_interneuron_CREs.bed
echo ""

echo "First 3 lines of Excitatory CREs:"
head -3 output/excitatory_neuron_CREs.bed
echo ""

echo "========================================================================"
echo "Running Python script with debug output..."
echo "========================================================================"
echo ""

# Run the Python script
python compare_GABA_vs_Excitatory.py

EXIT_CODE=$?

echo ""
echo "========================================================================"

if [ $EXIT_CODE -eq 0 ]; then
    echo "SCRIPT COMPLETED!"
    echo "========================================================================"
    echo ""
    echo "Check the output above for DEBUG messages"
    echo ""
    echo "Key things to look for:"
    echo "  1. Are different BED files being loaded?"
    echo "  2. Do the CRE coordinates differ between GABA and Excitatory?"
    echo "  3. Do the signal matrices differ?"
    echo ""
    echo "Completed: $(date)"
    echo "========================================================================"
else
    echo "SCRIPT FAILED!"
    echo "========================================================================"
    echo ""
    exit 1
fi
