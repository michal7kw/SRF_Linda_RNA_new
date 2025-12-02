#!/bin/bash
#SBATCH --job-name=2_convert_splicing_CREs_to_bed
#SBATCH --output=logs/2_convert_splicing_CREs_to_bed.log
#SBATCH --error=logs/2_convert_splicing_CREs_to_bed.err
#SBATCH --time=10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=workq

################################################################################
# Convert Splicing Gene CRE TSV Files to BED Format
#
# This script converts the CRE-gene linkage TSV files to BED format for use
# with deepTools visualization tools.
#
# INPUT:
# - output/splicing_genes_analysis/splicing_genes_CREs_all_celltypes.tsv
# - output/splicing_genes_analysis/splicing_genes_CREs_GABA.tsv
#
# OUTPUT:
# - output/splicing_genes_analysis/splicing_genes_CREs_all.bed
# - output/splicing_genes_analysis/splicing_genes_CREs_GABA.bed
#
# Runtime: ~1 minute
#
# Usage:
#   sbatch python 2_convert_splicing_CREs_to_bed.sh
################################################################################

echo "========================================================================"
echo "CONVERT SPLICING GENE CRE TSV FILES TO BED FORMAT"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Navigate to CREs_from_literature directory
cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_genes_paper"
OUTPUT_DIR="./output"

# ============================================================================
# Activate Conda Environment
# ============================================================================

echo "Activating conda environment..."
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# ============================================================================
# Run Conversion
# ============================================================================

echo "Converting TSV files to BED format..."
echo ""

python 2_convert_splicing_CREs_to_bed.py

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Conversion failed!"
    exit 1
fi

# ============================================================================
# Verify Output Files
# ============================================================================

echo ""
echo "========================================================================"
echo "VERIFYING OUTPUT FILES"
echo "========================================================================"
echo ""

for FILE in "splicing_genes_CREs_all.bed" "splicing_genes_CREs_GABA.bed"; do
    FULL_PATH="${OUTPUT_DIR}/${FILE}"
    if [ -f "$FULL_PATH" ]; then
        NUM_LINES=$(wc -l < "$FULL_PATH")
        echo "✓ $FILE: $NUM_LINES CREs"
        echo "  First 3 lines:"
        head -n 3 "$FULL_PATH" | awk '{print "    " $0}'
    else
        echo "✗ $FILE: NOT FOUND"
    fi
    echo ""
done

# ============================================================================
# Summary
# ============================================================================

echo "========================================================================"
echo "CONVERSION COMPLETE!"
echo "========================================================================"
echo ""
echo "Output directory:"
echo "  output/splicing_genes_analysis/"
echo ""
echo "Generated BED files:"
echo "  1. splicing_genes_CREs_all.bed (all cell types)"
echo "  2. splicing_genes_CREs_GABA.bed (GABA cell types)"
echo ""
echo "These files are ready for use with deepTools:"
echo "  - computeMatrix"
echo "  - plotHeatmap"
echo "  - plotProfile"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
