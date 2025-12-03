#!/bin/bash
#SBATCH --job-name=0_CREs_encode_all
#SBATCH --output=logs/0_CREs_encode_all_%j.log
#SBATCH --error=logs/0_CREs_encode_all_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=workq

################################################################################
# Master script to run the complete ENCODE cCREs Analysis Pipeline
#
# IMPORTANT: Emx1-Ctrl is a FAILED SAMPLE and is excluded from all analyses.
# Only 3 conditions are analyzed: Nestin-Ctrl, Nestin-Mut, Emx1-Mut.
#
# This script submits the following jobs to SLURM with dependencies:
# 1. 1_extract_encode_cCREs.sh - Extracts ENCODE cCREs linked to splicing genes
# 2. 2_convert_encode_cCREs_to_bed.sh - Converts TSV to BED format (all + by type)
# 4. 4_create_heatmaps_encode_cCREs.sh - Creates heatmaps and metaprofiles for all CRE types
# 5a. 5_create_custom_comparisons.sh - Computes matrices for custom comparisons
# 5b. 5_visualize_custom_comparisons.sh - Creates custom comparison visualizations
#      - Run with MIN_SIGNAL=2.0, MIN_FC=2.0 (individual plots)
#      - Run with MIN_SIGNAL=2.0, MIN_FC=3.0 (individual plots)
#
# Usage: sbatch 0_RUN_ENCODE_CCRES_ANALYSIS.sh
################################################################################

# Change to script directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_encode_all

# Create logs directory if it doesn't exist
mkdir -p logs

echo "========================================================================"
echo "SUBMITTING ENCODE cCREs ANALYSIS PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo "Working directory: $(pwd)"
echo ""
echo "NOTE: Emx1-Ctrl is EXCLUDED (failed sample)"
echo "      Only 3 conditions analyzed: Nestin-Ctrl, Nestin-Mut, Emx1-Mut"
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

# Submit Step 5a: Create custom comparisons (depends on Step 4)
JOB5A=$(sbatch --parsable --dependency=afterok:$JOB4 5_create_custom_comparisons.sh)
echo "Submitted Step 5a: 5_create_custom_comparisons.sh"
echo "  Job ID: $JOB5A"
echo "  Dependency: afterok:$JOB4"
echo "  Description: Compute matrices for custom comparisons"
echo ""

# Submit Step 5b-1: Custom comparison visualizations with MIN_SIGNAL=2.0, MIN_FC=2.0 (depends on Step 5a)
JOB5B1=$(sbatch --parsable --dependency=afterok:$JOB5A --export=ALL,MIN_SIGNAL=2.0,MIN_FC=2.0 5_visualize_custom_comparisons.sh --full)
echo "Submitted Step 5b-1: 5_visualize_custom_comparisons.sh (MIN_SIGNAL=2.0, MIN_FC=2.0)"
echo "  Job ID: $JOB5B1"
echo "  Dependency: afterok:$JOB5A"
echo "  Description: Create individual CRE plots (minSig=2.0, minFC=2.0)"
echo ""

# Submit Step 5b-2: Custom comparison visualizations with MIN_SIGNAL=2.0, MIN_FC=3.0 (depends on Step 5a, runs in parallel with 5b-1)
JOB5B2=$(sbatch --parsable --dependency=afterok:$JOB5A --export=ALL,MIN_SIGNAL=2.0,MIN_FC=3.0 5_visualize_custom_comparisons.sh --full)
echo "Submitted Step 5b-2: 5_visualize_custom_comparisons.sh (MIN_SIGNAL=2.0, MIN_FC=3.0)"
echo "  Job ID: $JOB5B2"
echo "  Dependency: afterok:$JOB5A"
echo "  Description: Create individual CRE plots (minSig=2.0, minFC=3.0)"
echo ""

echo "========================================================================"
echo "PIPELINE SUBMITTED SUCCESSFULLY"
echo "========================================================================"
echo ""
echo "Monitor status with: squeue -u $USER"
echo ""
echo "Job chain:"
echo "  $JOB1 (Extract) -> $JOB2 (Convert) -> $JOB4 (Heatmaps) -> $JOB5A (Custom matrices)"
echo "                                                              |"
echo "                                                    +---------+---------+"
echo "                                                    |                   |"
echo "                                              $JOB5B1            $JOB5B2"
echo "                                         (FC>=2.0)            (FC>=3.0)"
echo ""
echo "Expected runtime:"
echo "  Step 1: ~30-45 minutes (extracting overlaps)"
echo "  Step 2: ~5-10 minutes (converting to BED)"
echo "  Step 4: ~2-3 hours (creating heatmaps for all CRE types)"
echo "  Step 5a: ~1-2 hours (computing custom comparison matrices)"
echo "  Step 5b-1/2: ~1-2 hours each (individual CRE plots, run in parallel)"
echo "  Total: ~5-7 hours"
echo ""
echo "Output will be in: ./output/"
echo "  - Heatmaps: ./output/heatmaps_deeptools/"
echo "  - Comparisons: ./output/custom_comparisons/"
echo "  - Individual plots (FC>=2.0): ./output/custom_comparisons_minSig2.0_minFC2.0/"
echo "  - Individual plots (FC>=3.0): ./output/custom_comparisons_minSig2.0_minFC3.0/"
echo "========================================================================"
echo ""
echo "Completed: $(date)"
