#!/bin/bash
# Master script to run the complete Specific CREs pipeline (OLD workflow)
#
# This script submits the following jobs to SLURM with dependencies:
# 1. 1_extract_hippocampal_interneuron_CREs.sh (Extracts GABA CREs)
# 2. 2_RUN_SPECIFIC_CRES_ANALYSIS.sh (Extracts Excitatory CREs and runs deepTools analysis)
# 3. 3_fold_enrichment_specific_analysis.sh (Calculates fold enrichment)
# 4. 4_compare_GABA_vs_Excitatory.sh (Compares GABA vs Excitatory signal)
# 5. 5_create_heatmaps_metaprofiles.sh (Creates additional visualizations)
#
# Usage: ./0_RUN_COMPLETE_SPECIFIC_PIPELINE.sh

echo "========================================================================"
echo "SUBMITTING SPECIFIC CRES PIPELINE (OLD WORKFLOW)"
echo "========================================================================"
echo "Started: $(date)"
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

echo "========================================================================"
echo "PIPELINE SUBMITTED SUCCESSFULLY"
echo "========================================================================"
echo "Monitor status with: squeue -u $USER"
echo ""
