#!/bin/bash
# Master script to run the complete Splicing Genes ENCODE cCRE Analysis Pipeline
#
# This script submits the following jobs to SLURM with dependencies:
# 1. 1_extract_encode_cCREs.sh (Extracts ENCODE cCREs linked to splicing genes)
# 2. 2_convert_to_bed.sh (Converts TSV to BED format)
# 3. 3_create_profiles.sh (Creates profiles with deepTools)
# 4. 4_create_heatmaps.sh (Creates heatmaps with deepTools)
# 5. 5_visualize_bigwig_signal.sh (Visualizes BigWig signal directly)
# 6. 6_create_custom_comparisons.sh (Creates custom comparison plots)
#
# DATA SOURCE: ENCODE cCREs (mm10-cCREs.bed) linked by genomic proximity
#
# Usage: ./0_RUN_CREs_splicing_genes_encode_ANALYSIS.sh

echo "========================================================================"
echo "SUBMITTING SPLICING GENES ENCODE cCRE PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo ""
echo "DATA SOURCE: ENCODE cCREs (mm10-cCREs.bed)"
echo "LINKAGE METHOD: Genomic proximity (+/- 500kb)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_genes_encode"

# Make all scripts executable
chmod +x *.sh *.py

# Submit Step 1
JOB1=$(sbatch --parsable 1_extract_encode_cCREs.sh)
echo "Submitted Step 1: 1_extract_encode_cCREs.sh"
echo "  Job ID: $JOB1"
echo ""

# Submit Step 2 (depends on Step 1)
# This script needs the TSV files generated in Step 1
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 2_convert_to_bed.sh)
echo "Submitted Step 2: 2_convert_to_bed.sh"
echo "  Job ID: $JOB2"
echo "  Dependency: afterok:$JOB1"
echo ""

# Submit Step 3 (depends on Step 2)
# This script needs the BED files generated in Step 2
JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 3_create_profiles.sh)
echo "Submitted Step 3: 3_create_profiles.sh"
echo "  Job ID: $JOB3"
echo "  Dependency: afterok:$JOB2"
echo ""

# Submit Step 4 (depends on Step 2)
# This script needs the BED files generated in Step 2
JOB4=$(sbatch --parsable --dependency=afterok:$JOB2 4_create_heatmaps.sh)
echo "Submitted Step 4: 4_create_heatmaps.sh"
echo "  Job ID: $JOB4"
echo "  Dependency: afterok:$JOB2"
echo ""

# Submit Step 5 (depends on Step 2)
# This script needs the BED files generated in Step 2
JOB5=$(sbatch --parsable --dependency=afterok:$JOB2 5_visualize_bigwig_signal.sh)
echo "Submitted Step 5: 5_visualize_bigwig_signal.sh"
echo "  Job ID: $JOB5"
echo "  Dependency: afterok:$JOB2"
echo ""

# Submit Step 6 (depends on Step 2)
# This script needs the BED files generated in Step 2
JOB6=$(sbatch --parsable --dependency=afterok:$JOB2 6_create_custom_comparisons.sh)
echo "Submitted Step 6: 6_create_custom_comparisons.sh"
echo "  Job ID: $JOB6"
echo "  Dependency: afterok:$JOB2"
echo ""

echo "========================================================================"
echo "PIPELINE SUBMITTED SUCCESSFULLY"
echo "========================================================================"
echo ""
echo "Pipeline Structure:"
echo "  Step 1: Extract ENCODE cCREs -> Step 2: Convert to BED"
echo "                                        |"
echo "                                        +-> Step 3: Profiles"
echo "                                        +-> Step 4: Heatmaps"
echo "                                        +-> Step 5: BigWig Visualization"
echo "                                        +-> Step 6: Custom Comparisons"
echo ""
echo "Monitor status with: squeue -u $USER"
echo ""
echo "Output will be in: ./output/"
echo ""
