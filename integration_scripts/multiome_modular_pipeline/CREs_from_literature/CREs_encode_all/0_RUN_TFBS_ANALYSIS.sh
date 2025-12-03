#!/bin/bash
# =============================================================================
# SCRIPT: 0_RUN_TFBS_ANALYSIS.sh
# PURPOSE: Master script to run the complete TFBS analysis pipeline
#
# DESCRIPTION:
# This script orchestrates the entire Transcription Factor Binding Site (TFBS)
# analysis pipeline with SLURM job dependencies. It performs:
#
#   Step 6: Intersect ATAC peaks with PLS (Promoter-Like Signature) elements
#   Step 7: Run HOMER motif analysis on peaks overlapping PLS
#   Step 8: Differential motif enrichment comparison between conditions
#   Step 9: Visualization and summary figure generation
#
# USAGE:
#   ./0_RUN_TFBS_ANALYSIS.sh           # Run complete pipeline
#   ./0_RUN_TFBS_ANALYSIS.sh --dry-run # Show what would be run without executing
#   ./0_RUN_TFBS_ANALYSIS.sh --from 7  # Start from step 7 (skip intersection)
#
# REQUIREMENTS:
# - SLURM job scheduler
# - Conda environments: peak_calling_new (bedtools), annotation_enrichment (HOMER),
#   analysis3_env (Python: pandas, matplotlib, seaborn)
# - ENCODE cCREs already extracted (Steps 1-2 of main pipeline)
#
# OUTPUT:
# - ATAC peaks overlapping PLS: output/tfbs_analysis/peaks_on_PLS/
# - HOMER motif results: output/tfbs_analysis/homer_results/
# - Differential analysis: output/tfbs_analysis/differential_analysis/
# - Figures: output/tfbs_analysis/figures/
#
# AUTHOR: Claude Code
# DATE: December 2024
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/logs"
OUTPUT_DIR="${SCRIPT_DIR}/output/tfbs_analysis"

# Parse arguments
DRY_RUN=false
START_STEP=6

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --from)
            START_STEP=$2
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--dry-run] [--from STEP]"
            exit 1
            ;;
    esac
done

# =============================================================================
# FUNCTIONS
# =============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $1"
}

log_step() {
    echo ""
    echo "=========================================="
    echo "  STEP $1: $2"
    echo "=========================================="
}

submit_job() {
    local script="$1"
    local dependency="$2"
    local job_name="$3"

    if [[ ! -f "$script" ]]; then
        log_info "ERROR: Script not found: $script"
        exit 1
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        if [[ -n "$dependency" ]]; then
            echo "[DRY-RUN] sbatch --dependency=afterok:$dependency $script"
        else
            echo "[DRY-RUN] sbatch $script"
        fi
        echo "JOB_ID_PLACEHOLDER"
    else
        if [[ -n "$dependency" ]]; then
            job_id=$(sbatch --dependency=afterok:$dependency "$script" | awk '{print $4}')
        else
            job_id=$(sbatch "$script" | awk '{print $4}')
        fi
        log_info "Submitted $job_name: Job ID $job_id"
        echo "$job_id"
    fi
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

echo ""
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║          TFBS ANALYSIS PIPELINE FOR PLS ELEMENTS                     ║"
echo "║                                                                      ║"
echo "║  Transcription Factor Binding Site Search within ATAC Peaks          ║"
echo "║  at Promoter-Like Signature (PLS) Elements                           ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

log_info "Script directory: $SCRIPT_DIR"
log_info "Output directory: $OUTPUT_DIR"
log_info "Starting from step: $START_STEP"
log_info "Dry run: $DRY_RUN"

# Create directories
mkdir -p "${LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Verify prerequisite files exist
PLS_FILE="${SCRIPT_DIR}/output/encode_cCREs_PLS.bed"
if [[ ! -f "$PLS_FILE" && $START_STEP -le 6 ]]; then
    log_info "ERROR: PLS BED file not found: $PLS_FILE"
    log_info "Please run the ENCODE cCREs extraction pipeline first (Steps 1-2)"
    exit 1
fi

# Track job IDs for dependencies
JOB_ID_STEP6=""
JOB_ID_STEP7=""
JOB_ID_STEP8=""
JOB_ID_STEP9=""

# =============================================================================
# STEP 6: Intersect ATAC peaks with PLS elements
# =============================================================================

if [[ $START_STEP -le 6 ]]; then
    log_step 6 "Intersect ATAC peaks with PLS elements"
    log_info "Input: ENCODE PLS BED files + ATAC peak files"
    log_info "Output: Condition-specific peaks overlapping PLS"

    JOB_ID_STEP6=$(submit_job "${SCRIPT_DIR}/6_intersect_peaks_with_PLS.sh" "" "Step6_Intersect")
fi

# =============================================================================
# STEP 7: HOMER motif analysis (array job for 3 conditions)
# =============================================================================

if [[ $START_STEP -le 7 ]]; then
    log_step 7 "HOMER TFBS motif analysis"
    log_info "Input: HOMER-formatted BED files from Step 6"
    log_info "Output: De novo and known motif enrichment for each condition"
    log_info "Note: This is an array job processing 3 conditions in parallel"

    if [[ -n "$JOB_ID_STEP6" ]]; then
        JOB_ID_STEP7=$(submit_job "${SCRIPT_DIR}/7_run_homer_tfbs_analysis.sh" "$JOB_ID_STEP6" "Step7_HOMER")
    else
        JOB_ID_STEP7=$(submit_job "${SCRIPT_DIR}/7_run_homer_tfbs_analysis.sh" "" "Step7_HOMER")
    fi
fi

# =============================================================================
# STEP 8: Differential motif enrichment analysis
# =============================================================================

if [[ $START_STEP -le 8 ]]; then
    log_step 8 "Differential motif enrichment analysis"
    log_info "Input: HOMER knownResults.txt files from Step 7"
    log_info "Output: Differential enrichment tables and comparisons"

    if [[ -n "$JOB_ID_STEP7" ]]; then
        JOB_ID_STEP8=$(submit_job "${SCRIPT_DIR}/8_differential_motif_analysis.sh" "$JOB_ID_STEP7" "Step8_Differential")
    else
        JOB_ID_STEP8=$(submit_job "${SCRIPT_DIR}/8_differential_motif_analysis.sh" "" "Step8_Differential")
    fi
fi

# =============================================================================
# STEP 9: Visualization
# =============================================================================

if [[ $START_STEP -le 9 ]]; then
    log_step 9 "Visualization and summary figures"
    log_info "Input: Differential analysis results from Step 8"
    log_info "Output: Publication-quality figures (PNG/PDF)"

    if [[ -n "$JOB_ID_STEP8" ]]; then
        JOB_ID_STEP9=$(submit_job "${SCRIPT_DIR}/9_visualize_tfbs_results.sh" "$JOB_ID_STEP8" "Step9_Visualization")
    else
        JOB_ID_STEP9=$(submit_job "${SCRIPT_DIR}/9_visualize_tfbs_results.sh" "" "Step9_Visualization")
    fi
fi

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "=========================================="
echo "  PIPELINE SUBMITTED"
echo "=========================================="
echo ""
echo "Job Chain:"
if [[ -n "$JOB_ID_STEP6" ]]; then echo "  Step 6 (Intersection):  $JOB_ID_STEP6"; fi
if [[ -n "$JOB_ID_STEP7" ]]; then echo "  Step 7 (HOMER):         $JOB_ID_STEP7"; fi
if [[ -n "$JOB_ID_STEP8" ]]; then echo "  Step 8 (Differential):  $JOB_ID_STEP8"; fi
if [[ -n "$JOB_ID_STEP9" ]]; then echo "  Step 9 (Visualization): $JOB_ID_STEP9"; fi
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo "  tail -f ${LOG_DIR}/*.out"
echo ""
echo "Expected runtime: ~6-8 hours total"
echo "  - Step 6: ~5 minutes"
echo "  - Step 7: ~4-6 hours (HOMER is computationally intensive)"
echo "  - Step 8: ~5 minutes"
echo "  - Step 9: ~5 minutes"
echo ""
echo "Output directories:"
echo "  - Peaks on PLS: ${OUTPUT_DIR}/peaks_on_PLS/"
echo "  - HOMER results: ${OUTPUT_DIR}/homer_results/"
echo "  - Differential: ${OUTPUT_DIR}/differential_analysis/"
echo "  - Figures: ${OUTPUT_DIR}/figures/"
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "NOTE: This was a dry run. No jobs were actually submitted."
fi
