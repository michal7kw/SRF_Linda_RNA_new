#!/bin/bash
#SBATCH --job-name=2_convert_to_bed
#SBATCH --output=logs/2_convert_to_bed.log
#SBATCH --error=logs/2_convert_to_bed.err
#SBATCH --time=10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=workq

################################################################################
# Convert Splicing ENCODE cCRE TSV Files to BED Format
#
# This script converts the CRE-gene linkage TSV files to BED format for use
# with deepTools visualization tools.
#
# INPUT:
# - output/splicing_encode_cCREs_all.tsv
# - output/splicing_encode_cCREs_by_type.tsv
#
# OUTPUT:
# - output/splicing_encode_cCREs_all.bed
# - output/splicing_encode_cCREs_GABA.bed
# - output/splicing_encode_cCREs_{type}.bed (per CRE type)
#
# Runtime: ~1 minute
#
# Usage:
#   sbatch 2_convert_to_bed.sh
################################################################################

echo "========================================================================"
echo "CONVERT SPLICING ENCODE cCRE TSV FILES TO BED FORMAT"
echo "========================================================================"
echo "Started: $(date)"
echo ""

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/splicing_encode_cCREs"
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

python 2_convert_to_bed.py

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

for FILE in "splicing_encode_cCREs_all.bed" "splicing_encode_cCREs_GABA.bed"; do
    FULL_PATH="${OUTPUT_DIR}/${FILE}"
    if [ -f "$FULL_PATH" ]; then
        NUM_LINES=$(wc -l < "$FULL_PATH")
        echo "OK $FILE: $NUM_LINES CREs"
        echo "  First 3 lines:"
        head -n 3 "$FULL_PATH" | awk '{print "    " $0}'
    else
        echo "MISSING $FILE: NOT FOUND"
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
echo "  output/"
echo ""
echo "Generated BED files:"
echo "  1. splicing_encode_cCREs_all.bed (all CREs)"
echo "  2. splicing_encode_cCREs_GABA.bed (for GABA analysis)"
echo "  3. splicing_encode_cCREs_{type}.bed (per CRE type)"
echo ""
echo "These files are ready for use with deepTools:"
echo "  - computeMatrix"
echo "  - plotHeatmap"
echo "  - plotProfile"
echo ""
echo "Completed: $(date)"
echo "========================================================================"
