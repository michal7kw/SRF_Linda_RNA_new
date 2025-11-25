#!/bin/bash
# Master script to run the complete ENCODE cCREs Analysis Pipeline
#
# This script submits the following jobs to SLURM with dependencies:
# 1. 1_extract_encode_cCREs.sh - Extracts ENCODE cCREs linked to splicing genes
# 2. 2_convert_encode_cCREs_to_bed.sh - Converts TSV to BED format (all + by type)
# 4. 4_create_heatmaps_encode_cCREs.sh - Creates heatmaps and metaprofiles for all CRE types
# 5. 5_visualize_custom_comparisons.sh - Creates custom comparison visualizations
#
# Usage: ./0_RUN_ENCODE_CCRES_ANALYSIS.sh

echo "========================================================================"
echo "SUBMITTING ENCODE cCREs ANALYSIS PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Submit Step 1: Extract ENCODE cCREs
JOB1=$(sbatch --parsable 1_extract_encode_cCREs.sh)
echo "Submitted Step 1: 1_extract_encode_cCREs.sh"
echo "  Job ID: $JOB1"
echo "  Description: Extract ENCODE cCREs linked to splicing genes"
echo ""

# Submit Step 2: Convert to BED format (depends on Step 1)
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 2_convert_encode_cCREs_to_bed.sh)
echo "Submitted Step 2: 2_convert_encode_cCREs_to_bed.sh"
echo "  Job ID: $JOB2"
echo "  Dependency: afterok:$JOB1"
echo "  Description: Convert TSV to BED format (all + by CRE type)"
echo ""

# Submit Step 4: Create heatmaps (depends on Step 2)
JOB4=$(sbatch --parsable --dependency=afterok:$JOB2 4_create_heatmaps_encode_cCREs.sh)
echo "Submitted Step 4: 4_create_heatmaps_encode_cCREs.sh"
echo "  Job ID: $JOB4"
echo "  Dependency: afterok:$JOB2"
echo "  Description: Create heatmaps and metaprofiles for all CRE types"
echo ""

# Submit Step 5: Custom comparison visualizations (depends on Step 4)
JOB5=$(sbatch --parsable --dependency=afterok:$JOB4 5_visualize_custom_comparisons.sh)
echo "Submitted Step 5: 5_visualize_custom_comparisons.sh"
echo "  Job ID: $JOB5"
echo "  Dependency: afterok:$JOB4"
echo "  Description: Create custom comparison visualizations (metaprofiles, statistics)"
echo ""

echo "========================================================================"
echo "PIPELINE SUBMITTED SUCCESSFULLY"
echo "========================================================================"
echo ""
echo "Monitor status with: squeue -u $USER"
echo ""
echo "Job chain:"
echo "  $JOB1 (Extract) → $JOB2 (Convert) → $JOB4 (Heatmaps) → $JOB5 (Comparisons)"
echo ""
echo "Expected runtime:"
echo "  Step 1: ~30-45 minutes (extracting overlaps)"
echo "  Step 2: ~5-10 minutes (converting to BED)"
echo "  Step 4: ~2-3 hours (creating heatmaps for all CRE types)"
echo "  Step 5: ~15-30 minutes (custom comparison visualizations)"
echo "  Total: ~3-4 hours"
echo ""
echo "Output will be in: ./output/"
echo "  - Heatmaps: ./output/heatmaps_deeptools/"
echo "  - Comparisons: ./output/custom_comparisons/"
echo "========================================================================"
