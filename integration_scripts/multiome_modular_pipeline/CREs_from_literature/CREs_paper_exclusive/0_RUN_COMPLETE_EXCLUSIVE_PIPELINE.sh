#!/bin/bash
#SBATCH --job-name=0_RUN_COMPLETE_EXCLUSIVE_PIPELINE
#SBATCH --output=logs/0_RUN_COMPLETE_EXCLUSIVE_PIPELINE.log
#SBATCH --error=logs/0_RUN_COMPLETE_EXCLUSIVE_PIPELINE.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=workq

################################################################################
# Master script to run the complete Exclusive CREs pipeline
#
# This script submits the following jobs to SLURM with dependencies:
# 1. 1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh (Extracts CREs and creates deepTools heatmaps)
#    - 1a_extract_cell_type_specific_CREs.py: Extract mutually exclusive CRE sets
#    - 1b_create_heatmaps_specific_CREs_deeptools.sh: Create deepTools heatmaps
# 2. 2_fold_enrichment_exclusive_analysis.sh (Calculates fold enrichment from matrices)
# 3. 3_compare_cell_type_specific_CREs.sh (Generates pyBigWig comparison plots)
# 4. 3_link_CREs_to_genes.sh (Links CREs to genes using Table 16)
# 5. 4_visualize_DA_CREs.sh (Visualizes differentially accessible CREs)
#
# Pipeline Order with Dependencies:
#   Step 1 ───┬─→ Step 2 (fold enrichment)
#             ├─→ Step 3 (pyBigWig comparison)
#             └─→ Step 4 (link CREs to genes) ──→ Step 5 (visualize DA CREs)
#
# Samples analyzed (Emx1-Ctrl EXCLUDED - failed sample):
#   - Nestin-Ctrl (reference control)
#   - Nestin-Mut
#   - Emx1-Mut
#
# Usage: sbatch 0_RUN_COMPLETE_EXCLUSIVE_PIPELINE.sh
#
# Utility Scripts (not in main pipeline):
# - REPLOT_ONLY.sh: Re-run plotting without recomputing matrices
################################################################################

# Change to script directory (use absolute path for SLURM compatibility)
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_paper_exclusive"
cd "$SCRIPT_DIR" || { echo "ERROR: Cannot cd to $SCRIPT_DIR"; exit 1; }

echo "========================================================================"
echo "SUBMITTING EXCLUSIVE CRES PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Create logs directory
mkdir -p logs

echo "Pipeline Steps:"
echo "  1. Extract CREs + deepTools heatmaps"
echo "  2. Fold enrichment analysis"
echo "  3. pyBigWig comparison plots"
echo "  4. Link CREs to genes"
echo "  5. Visualize DA CREs"
echo ""
echo "Samples: Nestin-Ctrl, Nestin-Mut, Emx1-Mut"
echo "NOTE: Emx1-Ctrl EXCLUDED (failed sample)"
echo ""

# Submit Step 1
# This script runs 1a (extraction) and 1b (heatmaps)
JOB1=$(sbatch --parsable 1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh)
echo "Submitted Step 1: 1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh"
echo "  Job ID: $JOB1"
echo "  Runs: 1a_extract_cell_type_specific_CREs.py"
echo "        1b_create_heatmaps_specific_CREs_deeptools.sh"
echo ""

# Submit Step 2 (depends on Step 1)
# This script needs the matrices generated in Step 1
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 2_fold_enrichment_exclusive_analysis.sh)
echo "Submitted Step 2: 2_fold_enrichment_exclusive_analysis.sh"
echo "  Job ID: $JOB2"
echo "  Dependency: afterok:$JOB1"
echo ""

# Submit Step 3 (depends on Step 1)
# This script needs the BED files generated in Step 1
JOB3=$(sbatch --parsable --dependency=afterok:$JOB1 3_compare_cell_type_specific_CREs.sh)
echo "Submitted Step 3: 3_compare_cell_type_specific_CREs.sh"
echo "  Job ID: $JOB3"
echo "  Dependency: afterok:$JOB1"
echo ""

# Submit Step 4 (depends on Step 1)
# Link CREs to genes using Table 16
JOB4=$(sbatch --parsable --dependency=afterok:$JOB1 3_link_CREs_to_genes.sh)
echo "Submitted Step 4: 3_link_CREs_to_genes.sh"
echo "  Job ID: $JOB4"
echo "  Dependency: afterok:$JOB1"
echo ""

# Submit Step 5 (depends on Step 1 AND Step 4)
# Visualize DA CREs (needs matrices from Step 1 and gene linkage from Step 4)
JOB5=$(sbatch --parsable --dependency=afterok:$JOB1:$JOB4 4_visualize_DA_CREs.sh)
echo "Submitted Step 5: 4_visualize_DA_CREs.sh"
echo "  Job ID: $JOB5"
echo "  Dependency: afterok:$JOB1:$JOB4"
echo ""

echo "========================================================================"
echo "PIPELINE SUBMITTED SUCCESSFULLY"
echo "========================================================================"
echo ""
echo "Job Summary:"
echo "  Step 1 (JOB1): $JOB1 - Extract CREs + heatmaps"
echo "  Step 2 (JOB2): $JOB2 - Fold enrichment (after JOB1)"
echo "  Step 3 (JOB3): $JOB3 - pyBigWig comparison (after JOB1)"
echo "  Step 4 (JOB4): $JOB4 - Link CREs to genes (after JOB1)"
echo "  Step 5 (JOB5): $JOB5 - Visualize DA CREs (after JOB1+JOB4)"
echo ""
echo "Monitor status with: squeue -u $USER"
echo ""
echo "Output directories:"
echo "  - output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/"
echo "  - output/GABA_specific_CREs_genes.tsv"
echo "  - output/Excitatory_specific_CREs_genes.tsv"
echo ""
echo "Completed submission: $(date)"
echo "========================================================================"
