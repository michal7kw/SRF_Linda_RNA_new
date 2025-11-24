#!/bin/bash

################################################################################
# Add Sample Prefixes to Fragment Files
################################################################################
#
# This script creates new fragment files with sample prefixes added to barcodes.
# This is necessary because:
# - Fragment files have barcodes: "AAACAGCCACATTAAC-1-2"
# - Seurat cells have barcodes: "Nestin-Ctrl-adult_AAACAGCCACATTAAC-1-2"
#
# We need fragment files with prefixed barcodes to match the Seurat object.
#
# Usage:
#   bash scripts/add_prefix_to_fragments.sh <INPUT_DIR> <OUTPUT_DIR>
#
################################################################################

set -e
set -u
set -o pipefail

# Check arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <INPUT_DIR> <OUTPUT_DIR>"
    echo ""
    echo "Example:"
    echo "  $0 signac_results/translated_fragments signac_results/fragments_with_prefix"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

echo "================================================================================"
echo "ADDING SAMPLE PREFIXES TO FRAGMENT FILES"
echo "================================================================================"
echo ""
echo "Input directory:  $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Define sample mappings
declare -A SAMPLES
SAMPLES["R26-Nestin-Ctrl-adult"]="Nestin-Ctrl-adult_"
SAMPLES["R26-Nestin-Mut-adult"]="Nestin-Mut-adult_"
SAMPLES["R26-Emx1-Ctrl-adult"]="Emx1-Ctrl-adult_"
SAMPLES["R26-Emx1-Mut-adult"]="Emx1-Mut-adult_"

# Process each sample
for SAMPLE in "${!SAMPLES[@]}"; do
    PREFIX="${SAMPLES[$SAMPLE]}"
    INPUT_FILE="${INPUT_DIR}/${SAMPLE}_fragments_translated.tsv.gz"
    OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE}_fragments_with_prefix.tsv.gz"

    echo "Processing: $SAMPLE"
    echo "  Prefix: $PREFIX"
    echo "  Input:  $(basename $INPUT_FILE)"
    echo "  Output: $(basename $OUTPUT_FILE)"

    # Check input file exists
    if [ ! -f "$INPUT_FILE" ]; then
        echo "  ✗ ERROR: Input file not found!"
        exit 1
    fi

    # Add prefix to barcode column (column 4)
    # Fragment file format: chr start end barcode count
    echo "  Adding prefix to barcodes..."

    zcat "$INPUT_FILE" | awk -v prefix="$PREFIX" 'BEGIN{OFS="\t"} {$4=prefix$4; print}' | \
        bgzip -c > "$OUTPUT_FILE"

    # Index with tabix
    echo "  Indexing with tabix..."
    tabix -p bed "$OUTPUT_FILE"

    # Verify output
    OUTPUT_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
    INDEX_SIZE=$(du -h "${OUTPUT_FILE}.tbi" | cut -f1)

    echo "  ✓ Created: $OUTPUT_SIZE (data) + $INDEX_SIZE (index)"
    echo ""
done

echo "================================================================================"
echo "COMPLETE"
echo "================================================================================"
echo ""
echo "Created fragment files with prefixed barcodes in:"
echo "  $OUTPUT_DIR"
echo ""
echo "These files are now ready to use with fix_fragment_links.R scripts."
echo ""
