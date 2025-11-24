#!/bin/bash
# Master script to run the complete Exclusive CREs pipeline
#
# This script submits the following jobs to SLURM with dependencies:
# 1. 1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh (Extracts CREs and creates heatmaps)
# 2. 2_fold_enrichment_exclusive_analysis.sh (Calculates fold enrichment from matrices)
# 3. 3_compare_cell_type_specific_CREs.sh (Generates comparison plots)
#
# Usage: ./0_RUN_COMPLETE_EXCLUSIVE_PIPELINE.sh

echo "========================================================================"
echo "SUBMITTING EXCLUSIVE CRES PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Submit Step 1
# This script runs 1a (extraction) and 1b (heatmaps)
JOB1=$(sbatch --parsable 1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh)
echo "Submitted Step 1: 1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh"
echo "  Job ID: $JOB1"
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

echo "========================================================================"
echo "PIPELINE SUBMITTED SUCCESSFULLY"
echo "========================================================================"
echo "Monitor status with: squeue -u $USER"
echo ""
