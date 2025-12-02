#!/bin/bash
# =============================================================================
# SPLICING GENE CRE ANALYSIS PIPELINE - WITH DIRECTIONALITY
# =============================================================================
#
# This pipeline analyzes chromatin accessibility at splicing gene regulatory
# elements, specifically designed to identify NESTIN-SPECIFIC LOSS of
# accessibility in GABA neurons.
#
# KEY INNOVATION: Separates CREs by correlation direction
#   - Enhancer-like CREs (PCC > 0): accessibility ↑ = expression ↑
#   - Silencer-like CREs (PCC < 0): accessibility ↑ = expression ↓
#
# RESEARCH QUESTION:
#   Identify loss of chromatin accessibility at splicing gene regulatory
#   elements specifically in Nestin-MUT GABA neurons compared to:
#   1. Nestin-Ctrl (within-genotype mutation effect)
#   2. Emx1-Mut (genotype-specific response)
#
# BIOLOGICAL HYPOTHESIS:
#   If splicing disturbances are caused by reduced splicing gene expression
#   in Nestin-MUT, we expect to see:
#   - Loss of accessibility at ENHANCER-like CREs
#   - Pattern: Nestin-Ctrl HIGH → Nestin-Mut LOW, Emx1-Mut HIGH/NORMAL
#
# PIPELINE STEPS:
#   1. Extract CREs with PCC directionality (enhancers vs silencers)
#   2. Compute signal matrices with deepTools
#   3. Visualize comparisons and identify Nestin-specific loss
#
# USAGE:
#   ./0_RUN_SPLICING_DIRECTIONAL_PIPELINE.sh [OPTIONS]
#
# OPTIONS:
#   --dry-run           Show what would be executed without running
#   --skip-individual   Skip individual CRE plots (faster)
#   --step N            Run only step N (1, 2, or 3)
#   --help              Show this help message
#
# =============================================================================

set -e

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# Create logs directory
mkdir -p logs

# Default options
DRY_RUN=0
SKIP_INDIVIDUAL=0
SINGLE_STEP=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --skip-individual)
            SKIP_INDIVIDUAL=1
            shift
            ;;
        --step)
            SINGLE_STEP="$2"
            shift 2
            ;;
        --help|-h)
            head -50 "$0" | grep "^#" | sed 's/^# //' | sed 's/^#//'
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Export options for child scripts
export SKIP_INDIVIDUAL

echo "=============================================================================="
echo "SPLICING GENE CRE ANALYSIS PIPELINE - WITH DIRECTIONALITY"
echo "=============================================================================="
echo ""
echo "Date: $(date)"
echo "Working directory: ${SCRIPT_DIR}"
echo ""
echo "Options:"
echo "  Dry run: ${DRY_RUN}"
echo "  Skip individual plots: ${SKIP_INDIVIDUAL}"
echo "  Single step: ${SINGLE_STEP:-'all'}"
echo ""

# =============================================================================
# STEP 1: Extract CREs with PCC directionality
# =============================================================================

run_step1() {
    echo "=============================================================================="
    echo "STEP 1: Extract Splicing Gene CREs with Directionality"
    echo "=============================================================================="

    if [[ ${DRY_RUN} -eq 1 ]]; then
        echo "[DRY RUN] Would execute: sbatch 1_extract_splicing_CREs_directional.sh"
        echo "step1_job_id=12345"
        return
    fi

    step1_job=$(sbatch --parsable 1_extract_splicing_CREs_directional.sh)
    echo "Submitted Step 1: Job ID ${step1_job}"
    echo "${step1_job}"
}

# =============================================================================
# STEP 2: Compute signal matrices
# =============================================================================

run_step2() {
    local depend="$1"

    echo ""
    echo "=============================================================================="
    echo "STEP 2: Compute Signal Matrices"
    echo "=============================================================================="

    if [[ ${DRY_RUN} -eq 1 ]]; then
        echo "[DRY RUN] Would execute: sbatch --dependency=afterok:${depend} 2_compute_signal_matrices.sh"
        echo "step2_job_id=12346"
        return
    fi

    if [[ -n "${depend}" ]]; then
        step2_job=$(sbatch --parsable --dependency=afterok:${depend} 2_compute_signal_matrices.sh)
    else
        step2_job=$(sbatch --parsable 2_compute_signal_matrices.sh)
    fi
    echo "Submitted Step 2: Job ID ${step2_job}"
    echo "${step2_job}"
}

# =============================================================================
# STEP 3: Visualize and identify Nestin-specific loss
# =============================================================================

run_step3() {
    local depend="$1"

    echo ""
    echo "=============================================================================="
    echo "STEP 3: Visualize Comparisons and Identify Nestin-Specific Loss"
    echo "=============================================================================="

    if [[ ${DRY_RUN} -eq 1 ]]; then
        echo "[DRY RUN] Would execute: sbatch --dependency=afterok:${depend} 3_visualize_directional_comparisons.sh"
        echo "step3_job_id=12347"
        return
    fi

    if [[ -n "${depend}" ]]; then
        step3_job=$(sbatch --parsable --dependency=afterok:${depend} 3_visualize_directional_comparisons.sh)
    else
        step3_job=$(sbatch --parsable 3_visualize_directional_comparisons.sh)
    fi
    echo "Submitted Step 3: Job ID ${step3_job}"
    echo "${step3_job}"
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if [[ -n "${SINGLE_STEP}" ]]; then
    # Run single step
    case ${SINGLE_STEP} in
        1)
            run_step1
            ;;
        2)
            run_step2 ""
            ;;
        3)
            run_step3 ""
            ;;
        *)
            echo "ERROR: Invalid step number: ${SINGLE_STEP}"
            echo "Valid steps: 1, 2, 3"
            exit 1
            ;;
    esac
else
    # Run full pipeline with dependencies
    step1_job_id=$(run_step1 | tail -1)
    step2_job_id=$(run_step2 "${step1_job_id}" | tail -1)
    step3_job_id=$(run_step3 "${step2_job_id}" | tail -1)

    echo ""
    echo "=============================================================================="
    echo "PIPELINE SUBMITTED"
    echo "=============================================================================="
    echo ""
    echo "Job chain:"
    echo "  Step 1 (Extract CREs):     ${step1_job_id}"
    echo "  Step 2 (Compute matrices): ${step2_job_id} (depends on ${step1_job_id})"
    echo "  Step 3 (Visualize):        ${step3_job_id} (depends on ${step2_job_id})"
    echo ""
    echo "Monitor progress:"
    echo "  squeue -u \$USER"
    echo ""
    echo "View logs:"
    echo "  tail -f logs/1_extract_splicing_CREs_*.log"
    echo "  tail -f logs/2_compute_signal_matrices_*.log"
    echo "  tail -f logs/3_visualize_directional_*.log"
    echo ""
    echo "Expected outputs:"
    echo "  output/enhancer_CREs_GABA.bed"
    echo "  output/silencer_CREs_GABA.bed"
    echo "  output/heatmaps_deeptools/*.png"
    echo "  output/visualizations/*.png"
    echo "  output/nestin_specific_loss/*.tsv"
    echo ""
    echo "Estimated total runtime: ~2-3 hours"
fi

echo ""
echo "Date: $(date)"
