#!/bin/bash
#SBATCH --job-name=0_RUN_SPLICING_GENES_ANALYSIS
#SBATCH --output=logs/0_RUN_SPLICING_GENES_ANALYSIS.log
#SBATCH --error=logs/0_RUN_SPLICING_GENES_ANALYSIS.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=workq

################################################################################
# Master script to run the complete Splicing Genes CRE Analysis Pipeline
#
# This script submits the following jobs to SLURM with dependencies:
# 1. 1_extract_splicing_gene_CREs.sh (Extracts CREs linked to splicing genes)
# 2. 2_convert_splicing_CREs_to_bed.sh (Converts TSV to BED format)
# 3. 3_create_splicing_profiles.sh (Creates profiles with deepTools)
# 4. 4_create_heatmaps_splicing_genes.sh (Creates heatmaps with deepTools)
# 5. 5_visualize_bigwig_signal.sh (Visualizes BigWig signal directly)
# 6. 6_create_custom_comparisons.sh (Creates custom comparison plots)
#
# NOTE: Emx1-Ctrl is a failed sample and excluded from all analyses.
#       Nestin-Ctrl is used as control for Emx1 comparisons.
#
# Usage: sbatch 0_RUN_SPLICING_GENES_ANALYSIS.sh
################################################################################

# Change to script directory
cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_genes_paper"

# Create logs directory if it doesn't exist
mkdir -p logs

echo "========================================================================"
echo "SUBMITTING SPLICING GENES PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Submit Step 1
JOB1=$(sbatch --parsable ./1_extract_splicing_gene_CREs.sh)
echo "Submitted Step 1: 1_extract_splicing_gene_CREs.sh"
echo "  Job ID: $JOB1"
echo ""

# Submit Step 2 (depends on Step 1)
# This script needs the TSV files generated in Step 1
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 ./2_convert_splicing_CREs_to_bed.sh)
echo "Submitted Step 2: 2_convert_splicing_CREs_to_bed.sh"
echo "  Job ID: $JOB2"
echo "  Dependency: afterok:$JOB1"
echo ""

# Submit Step 3 (depends on Step 2)
# This script needs the BED files generated in Step 2
JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 ./3_create_splicing_profiles.sh)
echo "Submitted Step 3: 3_create_splicing_profiles.sh"
echo "  Job ID: $JOB3"
echo "  Dependency: afterok:$JOB2"
echo ""

# Submit Step 4 (depends on Step 2)
# This script needs the BED files generated in Step 2
JOB4=$(sbatch --parsable --dependency=afterok:$JOB2 ./4_create_heatmaps_splicing_genes.sh)
echo "Submitted Step 4: 4_create_heatmaps_splicing_genes.sh"
echo "  Job ID: $JOB4"
echo "  Dependency: afterok:$JOB2"
echo ""

# Submit Step 5 (depends on Step 2)
# This script needs the BED files generated in Step 2
JOB5=$(sbatch --parsable --dependency=afterok:$JOB2 ./5_visualize_bigwig_signal.sh)
echo "Submitted Step 5: 5_visualize_bigwig_signal.sh"
echo "  Job ID: $JOB5"
echo "  Dependency: afterok:$JOB2"
echo ""

# Submit Step 6 (depends on Step 2)
# This script needs the BED files generated in Step 2
JOB6=$(sbatch --parsable --dependency=afterok:$JOB2 ./6_create_custom_comparisons.sh)
echo "Submitted Step 6: 6_create_custom_comparisons.sh"
echo "  Job ID: $JOB6"
echo "  Dependency: afterok:$JOB2"
echo ""

echo "========================================================================"
echo "PIPELINE SUBMITTED SUCCESSFULLY"
echo "========================================================================"
echo ""
echo "Job Summary:"
echo "  Step 1 (Extract CREs):      $JOB1"
echo "  Step 2 (Convert to BED):    $JOB2 (depends on $JOB1)"
echo "  Step 3 (Signal Profiles):   $JOB3 (depends on $JOB2)"
echo "  Step 4 (Heatmaps):          $JOB4 (depends on $JOB2)"
echo "  Step 5 (BigWig Viz):        $JOB5 (depends on $JOB2)"
echo "  Step 6 (Custom Comparisons):$JOB6 (depends on $JOB2)"
echo ""
echo "Monitor status with: squeue -u $USER"
echo "View logs in: logs/"
echo ""
echo "Completed submission: $(date)"
echo "========================================================================"
