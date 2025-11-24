#!/bin/bash

################################################################################
# Fix Fragment File Paths - Master Script for Both Pipelines
################################################################################
#
# This master script fixes fragment file paths for both L1 and L2 pipelines.
#
# Usage:
#   ./fix_fragments_all.sh [L1|L2|both]
#
# Examples:
#   ./fix_fragments_all.sh L2       # Fix L2 (fine-grained) only
#   ./fix_fragments_all.sh L1       # Fix L1 (broad) only
#   ./fix_fragments_all.sh both     # Fix both (default)
#   ./fix_fragments_all.sh          # Fix both (default)
#
################################################################################

set -e
set -u
set -o pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
cd "$BASE_DIR"

# Parse argument
MODE="${1:-both}"

# Validate mode
if [[ ! "$MODE" =~ ^(L1|L2|both)$ ]]; then
    echo "ERROR: Invalid mode '$MODE'"
    echo ""
    echo "Usage: $0 [L1|L2|both]"
    echo ""
    echo "  L1    - Fix L1 pipeline only (broad cell types)"
    echo "  L2    - Fix L2 pipeline only (fine-grained cell types)"
    echo "  both  - Fix both pipelines (default)"
    echo ""
    exit 1
fi

echo "================================================================================"
echo "FRAGMENT FILE PATH FIX - MASTER SCRIPT"
echo "================================================================================"
echo ""
echo "Mode: $MODE"
echo ""
echo "This will fix broken fragment file paths that cause empty coverage plots."
echo ""
echo "================================================================================"
echo ""

START_TIME=$(date +%s)

################################################################################
# Fix L2 Pipeline
################################################################################

if [[ "$MODE" == "L2" || "$MODE" == "both" ]]; then
    echo "================================================================================"
    echo "FIXING L2 PIPELINE"
    echo "================================================================================"
    echo ""

    if [ -f "$BASE_DIR/fix_fragments_L2.sh" ]; then
        bash "$BASE_DIR/fix_fragments_L2.sh"

        if [ $? -ne 0 ]; then
            echo ""
            echo "✗ ERROR: L2 pipeline fix failed!"
            exit 1
        fi
    else
        echo "✗ ERROR: L2 fix script not found: fix_fragments_L2.sh"
        exit 1
    fi

    echo ""
fi

################################################################################
# Fix L1 Pipeline
################################################################################

if [[ "$MODE" == "L1" || "$MODE" == "both" ]]; then
    if [[ "$MODE" == "both" ]]; then
        echo ""
        echo "================================================================================"
        echo "FIXING L1 PIPELINE"
        echo "================================================================================"
        echo ""
    fi

    if [ -f "$BASE_DIR/fix_fragments_L1.sh" ]; then
        bash "$BASE_DIR/fix_fragments_L1.sh"

        if [ $? -ne 0 ]; then
            echo ""
            echo "✗ ERROR: L1 pipeline fix failed!"
            exit 1
        fi
    else
        echo "✗ ERROR: L1 fix script not found: fix_fragments_L1.sh"
        exit 1
    fi

    echo ""
fi

################################################################################
# Summary
################################################################################

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
ELAPSED_MIN=$((ELAPSED / 60))
ELAPSED_SEC=$((ELAPSED % 60))

echo "================================================================================"
echo "ALL FIXES COMPLETE"
echo "================================================================================"
echo ""
echo "Mode:          $MODE"
echo "Total time:    ${ELAPSED_MIN}m ${ELAPSED_SEC}s"
echo ""

case $MODE in
    L2)
        echo "Fixed pipeline:"
        echo "  ✓ L2 (signac_results/) - Fine-grained cell types"
        echo ""
        echo "Next steps:"
        echo "  sbatch run_viz2_deg_coverage_enhanced.sh"
        echo "  sbatch scripts/run_analyze_peak_gene_linkage_gc_UPDATED.sh"
        echo "  sbatch scripts/run_analyze_peak_gene_linkage_gaba_UPDATED.sh"
        ;;
    L1)
        echo "Fixed pipeline:"
        echo "  ✓ L1 (signac_results_L1/) - Broad cell types"
        echo ""
        echo "Next steps:"
        echo "  sbatch scripts/run_analyze_peak_gene_linkage_gc_UPDATED_L1.sh"
        echo "  sbatch scripts/run_analyze_peak_gene_linkage_gaba_UPDATED_L1.sh"
        ;;
    both)
        echo "Fixed pipelines:"
        echo "  ✓ L2 (signac_results/) - Fine-grained cell types"
        echo "  ✓ L1 (signac_results_L1/) - Broad cell types"
        echo ""
        echo "Next steps (L2):"
        echo "  sbatch run_viz2_deg_coverage_enhanced.sh"
        echo "  sbatch scripts/run_analyze_peak_gene_linkage_gc_UPDATED.sh"
        echo "  sbatch scripts/run_analyze_peak_gene_linkage_gaba_UPDATED.sh"
        echo ""
        echo "Next steps (L1):"
        echo "  sbatch scripts/run_analyze_peak_gene_linkage_gc_UPDATED_L1.sh"
        echo "  sbatch scripts/run_analyze_peak_gene_linkage_gaba_UPDATED_L1.sh"
        ;;
esac

echo ""
echo "Monitor visualization jobs:"
echo "  squeue -u \$USER"
echo "  tail -f logs/run_viz2_deg_coverage_enhanced.out"
echo ""
echo "================================================================================"
