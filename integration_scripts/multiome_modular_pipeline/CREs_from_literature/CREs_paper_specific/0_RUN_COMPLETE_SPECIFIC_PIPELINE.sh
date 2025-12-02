#!/bin/bash
#SBATCH --job-name=0_RUN_COMPLETE_SPECIFIC_PIPELINE
#SBATCH --output=logs/0_RUN_COMPLETE_SPECIFIC_PIPELINE_%j.log
#SBATCH --error=logs/0_RUN_COMPLETE_SPECIFIC_PIPELINE_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=workq

################################################################################
# Master script to run the complete Specific CREs pipeline
#
# This script submits the following jobs to SLURM with dependencies:
# 1. 1_extract_hippocampal_interneuron_CREs.sh (Extracts GABA CREs)
# 2. 2_RUN_SPECIFIC_CRES_ANALYSIS.sh (Extracts Excitatory CREs and runs deepTools analysis)
# 3. 3_fold_enrichment_specific_analysis.sh (Calculates fold enrichment)
# 4. 4_compare_GABA_vs_Excitatory.sh (Compares GABA vs Excitatory signal)
# 5. 5_create_heatmaps_metaprofiles.sh (Creates additional visualizations)
# 6. 6_link_CREs_to_genes.sh (Links CREs to genes using Table 16)
# 7. 7_visualize_DA_CREs.sh (Visualizes Differentially Accessible CREs)
#
# Samples analyzed (Emx1-Ctrl excluded - failed sample):
# - Nestin-Ctrl
# - Nestin-Mut
# - Emx1-Mut
#
# Utility scripts (not in main pipeline):
# - REPLOT_ONLY_OLD.sh (Fast replotting from existing matrices)
#
# Usage: sbatch 0_RUN_COMPLETE_SPECIFIC_PIPELINE.sh
################################################################################

set -e  # Exit on error

# Use SLURM_SUBMIT_DIR (the directory where sbatch was called)
# This is more reliable than BASH_SOURCE when running under SLURM
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

# Ensure logs directory exists
mkdir -p logs

echo "========================================================================"
echo "SUBMITTING SPECIFIC CRES PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo "Working directory: $SCRIPT_DIR"
echo ""

# Submit Step 1
JOB1=$(sbatch --parsable 1_extract_hippocampal_interneuron_CREs.sh)
echo "Submitted Step 1: 1_extract_hippocampal_interneuron_CREs.sh"
echo "  Job ID: $JOB1"
echo ""

# Submit Step 2 (depends on Step 1)
# This script needs the GABA CREs generated in Step 1
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 2_RUN_SPECIFIC_CRES_ANALYSIS.sh)
echo "Submitted Step 2: 2_RUN_SPECIFIC_CRES_ANALYSIS.sh"
echo "  Job ID: $JOB2"
echo "  Dependency: afterok:$JOB1"
echo ""

# Submit Step 3 (depends on Step 2)
# This script needs the matrices generated in Step 2
JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 3_fold_enrichment_specific_analysis.sh)
echo "Submitted Step 3: 3_fold_enrichment_specific_analysis.sh"
echo "  Job ID: $JOB3"
echo "  Dependency: afterok:$JOB2"
echo ""

# Submit Step 4 (depends on Step 2)
# This script needs the BED files generated in Step 1 & 2
JOB4=$(sbatch --parsable --dependency=afterok:$JOB2 4_compare_GABA_vs_Excitatory.sh)
echo "Submitted Step 4: 4_compare_GABA_vs_Excitatory.sh"
echo "  Job ID: $JOB4"
echo "  Dependency: afterok:$JOB2"
echo ""

# Submit Step 5 (depends on Step 2)
# This script creates additional visualizations
JOB5=$(sbatch --parsable --dependency=afterok:$JOB2 5_create_heatmaps_metaprofiles.sh)
echo "Submitted Step 5: 5_create_heatmaps_metaprofiles.sh"
echo "  Job ID: $JOB5"
echo "  Dependency: afterok:$JOB2"
echo ""

# Submit Step 6 (depends on Step 2)
# This script links CREs to genes using Table 16
# Needs: hippocampal_interneuron_CREs.bed (Step 1), excitatory_neuron_CREs.bed (Step 2)
JOB6=$(sbatch --parsable --dependency=afterok:$JOB2 6_link_CREs_to_genes.sh)
echo "Submitted Step 6: 6_link_CREs_to_genes.sh"
echo "  Job ID: $JOB6"
echo "  Dependency: afterok:$JOB2"
echo ""

# Submit Step 7 (depends on Step 2 and Step 6)
# This script visualizes Differentially Accessible CREs
# Needs: matrices from Step 2, gene linkages from Step 6
JOB7=$(sbatch --parsable --dependency=afterok:$JOB2:$JOB6 7_visualize_DA_CREs.sh)
echo "Submitted Step 7: 7_visualize_DA_CREs.sh"
echo "  Job ID: $JOB7"
echo "  Dependency: afterok:$JOB2:$JOB6"
echo ""

echo "========================================================================"
echo "PIPELINE SUBMITTED SUCCESSFULLY"
echo "========================================================================"
echo ""
echo "Samples analyzed:"
echo "  - Nestin-Ctrl (control)"
echo "  - Nestin-Mut"
echo "  - Emx1-Mut"
echo "  NOTE: Emx1-Ctrl excluded (failed sample)"
echo ""
echo "Job dependency graph:"
echo ""
echo "  Step 1 (Extract GABA CREs)"
echo "       |"
echo "       v"
echo "  Step 2 (Main Analysis) ----+---> Step 3 (Fold Enrichment)"
echo "       |                     |"
echo "       |                     +---> Step 4 (GABA vs Excitatory)"
echo "       |                     |"
echo "       |                     +---> Step 5 (Heatmaps)"
echo "       |                     |"
echo "       +-------------------> Step 6 (CRE-Gene Links)"
echo "                                  |"
echo "                                  v"
echo "                             Step 7 (DA CRE Visualization)"
echo ""
echo "Job IDs:"
echo "  Step 1: $JOB1"
echo "  Step 2: $JOB2"
echo "  Step 3: $JOB3"
echo "  Step 4: $JOB4"
echo "  Step 5: $JOB5"
echo "  Step 6: $JOB6"
echo "  Step 7: $JOB7"
echo ""
echo "Monitor status with: squeue -u $USER"
echo "View logs in: $SCRIPT_DIR/logs/"
echo ""
echo "Completed submission: $(date)"
echo "========================================================================"
