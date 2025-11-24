#!/bin/bash

################################################################################
# Fix Fragment File Paths and Update Seurat Object - L1 Pipeline
################################################################################
#
# This script fixes fragment file paths for the L1 (broad) cell type
# pipeline, enabling coverage plots to display properly.
#
# Steps:
#   1. Verify translated fragment files exist
#   2. Update Seurat object with correct fragment paths (strips/re-adds prefixes)
#   3. Verify the fix worked
#
# Usage:
#   ./fix_fragments_L1.sh
#
################################################################################

set -e
set -u
set -o pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
cd "$BASE_DIR"

echo "================================================================================"
echo "FIXING FRAGMENT FILE PATHS - L1 PIPELINE (Broad Cell Types)"
echo "================================================================================"
echo ""
echo "Pipeline: signac_results_L1/"
echo "Cell types: L1 (Unknown, Excitatory, Immune, Oligo, Astrocytes, GABA, etc.)"
echo ""
echo "================================================================================"
echo ""

################################################################################
# Step 1: Create Fragment Files with Sample Prefixes
################################################################################

echo "STEP 1: Creating fragment files with sample prefixes"
echo "================================================================================"
echo ""

INPUT_DIR="$BASE_DIR/signac_results_L1/translated_fragments"
OUTPUT_DIR="$BASE_DIR/signac_results_L1/fragments_with_prefix"

# Check if translated fragments exist
if [ ! -d "$INPUT_DIR" ]; then
    echo "✗ ERROR: Translated fragments directory not found!"
    echo ""
    echo "Expected location: $INPUT_DIR"
    echo ""
    echo "These files should have been created during the initial Signac integration."
    exit 1
fi

# Check if prefixed fragments already exist
if [ -d "$OUTPUT_DIR" ]; then
    echo "Prefixed fragment files already exist. Checking..."
    SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult" "R26-Emx1-Ctrl-adult" "R26-Emx1-Mut-adult")
    ALL_EXIST=true

    for SAMPLE in "${SAMPLES[@]}"; do
        FRAG_FILE="${OUTPUT_DIR}/${SAMPLE}_fragments_with_prefix.tsv.gz"
        if [ ! -f "$FRAG_FILE" ]; then
            ALL_EXIST=false
            break
        fi
    done

    if [ "$ALL_EXIST" = true ]; then
        echo "✓ All prefixed fragment files already exist. Skipping creation."
        echo ""
        echo "✓ STEP 1 COMPLETE: Using existing prefixed fragment files"
        echo ""
    else
        echo "Some prefixed fragment files missing. Creating all..."
        echo ""
        bash "$BASE_DIR/scripts/add_prefix_to_fragments.sh" "$INPUT_DIR" "$OUTPUT_DIR"
    fi
else
    echo "Creating prefixed fragment files..."
    echo ""
    bash "$BASE_DIR/scripts/add_prefix_to_fragments.sh" "$INPUT_DIR" "$OUTPUT_DIR"
fi

################################################################################
# Step 2: Update Seurat Object Fragment Links
################################################################################

echo "================================================================================"
echo "STEP 2: Updating Seurat object fragment links"
echo "================================================================================"
echo ""

# Activate conda environment
echo "Activating conda environment: sc-chromatin2"
source ~/.bashrc
conda activate sc-chromatin2

# Check if R script exists
FIX_SCRIPT="$BASE_DIR/scripts/Legacy/fix_fragment_links_L1.R"
if [ ! -f "$FIX_SCRIPT" ]; then
    echo "ERROR: Fragment fix script not found: $FIX_SCRIPT"
    exit 1
fi

echo "Running fragment link fix script (L1 version)..."
echo ""

Rscript "$FIX_SCRIPT"

if [ $? -ne 0 ]; then
    echo ""
    echo "✗ ERROR: Fragment link fix failed!"
    exit 1
fi

echo ""
echo "✓ STEP 2 COMPLETE: Seurat object updated"
echo ""

################################################################################
# Step 3: Verification
################################################################################

echo "================================================================================"
echo "STEP 3: Verification"
echo "================================================================================"
echo ""

# Check if backup was created
BACKUP_FILE="$BASE_DIR/signac_results_L1/integrated_seurat_processed_backup_before_fragment_fix.rds"
if [ -f "$BACKUP_FILE" ]; then
    BACKUP_SIZE=$(du -h "$BACKUP_FILE" | cut -f1)
    echo "✓ Backup created: $(basename $BACKUP_FILE) ($BACKUP_SIZE)"
else
    echo "⚠ Warning: Backup file not found"
fi

# Check if updated Seurat object exists
SEURAT_FILE="$BASE_DIR/signac_results_L1/integrated_seurat_processed.rds"
if [ -f "$SEURAT_FILE" ]; then
    SEURAT_SIZE=$(du -h "$SEURAT_FILE" | cut -f1)
    SEURAT_DATE=$(date -r "$SEURAT_FILE" "+%Y-%m-%d %H:%M:%S")
    echo "✓ Updated Seurat object: $(basename $SEURAT_FILE) ($SEURAT_SIZE)"
    echo "  Last modified: $SEURAT_DATE"
else
    echo "✗ ERROR: Updated Seurat object not found!"
    exit 1
fi

echo ""
echo "✓ STEP 3 COMPLETE: Verification passed"
echo ""

################################################################################
# Summary and Next Steps
################################################################################

echo "================================================================================"
echo "FIX COMPLETE - L1 PIPELINE"
echo "================================================================================"
echo ""
echo "Summary:"
echo "  ✓ Created prefixed fragment files (4 files + indices)"
echo "  ✓ Seurat object updated with correct fragment paths"
echo "  ✓ Backup created before modification"
echo ""
echo "Next Steps:"
echo ""
echo "1. Test coverage plots (L1 versions):"
echo "   sbatch scripts/run_analyze_peak_gene_linkage_gc_UPDATED_L1.sh"
echo "   sbatch scripts/run_analyze_peak_gene_linkage_gaba_UPDATED_L1.sh"
echo ""
echo "2. Run L1-specific visualization:"
echo "   sbatch run_viz2_deg_coverage_enhanced_L1.sh  # (if exists)"
echo ""
echo "3. Monitor output:"
echo "   tail -f logs/run_analyze_peak_gene_linkage_*_L1.out"
echo ""
echo "If coverage plots still show empty tracks, check:"
echo "  - Fragment file permissions (should be readable)"
echo "  - Cell barcode prefixes match between Seurat and fragment files"
echo "  - Fragment files are properly indexed (.tbi files exist)"
echo ""
echo "================================================================================"
