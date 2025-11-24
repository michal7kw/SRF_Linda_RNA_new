#!/bin/bash
#SBATCH --job-name=cre_length_analysis
#SBATCH --output=logs/cre_length_analysis.out
#SBATCH --error=logs/cre_length_analysis.err
#SBATCH --partition=workq
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00

################################################################################
# Analysis of hippocampal interneuron CRE lengths
# Calculates statistics and generates visualizations
#
# INPUT FILES:
# - analyze_cre_lengths.py: Python script to analyze CRE lengths
# - output/hippocampal_interneuron_CREs_with_header.bed: BED file containing CRE coordinates with header
#
# OUTPUT FILES:
# - output/CRE_length_analysis.png: Comprehensive 4-panel visualization
# - output/CRE_length_histogram_simple.png: Simple histogram visualization
# - output/CRE_length_statistics.txt: Text file with detailed length statistics
# - logs/cre_length_analysis.out: Execution log file
# - logs/cre_length_analysis.err: Error log file
#
################################################################################

set -euo pipefail

# Activate deepTools environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "=================================================="
echo "CRE Length Analysis"
echo "=================================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo ""

# Run analysis
python analyze_cre_lengths.py

echo ""
echo "=================================================="
echo "Analysis completed: $(date)"
echo "=================================================="
