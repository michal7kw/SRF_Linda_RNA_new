#!/bin/bash

################################################################################
# Submit Peak-Gene Linkage Jobs in Batches with Delays
################################################################################
#
# Submits array jobs in smaller batches with delays to avoid concurrent
# conda environment access issues on shared BeeGFS filesystem
#
# Usage:
#   ./submit_in_batches.sh [L1|L2|both] [batch_size] [delay_seconds]
#
# Examples:
#   ./submit_in_batches.sh L1 3 60      # L1: 3 tasks at a time, 60s delay
#   ./submit_in_batches.sh L2 5 90      # L2: 5 tasks at a time, 90s delay
#   ./submit_in_batches.sh both 4 60    # Both: 4 tasks at a time, 60s delay
#
################################################################################

# Parse arguments
MODE=${1:-both}           # L1, L2, or both
BATCH_SIZE=${2:-4}        # Tasks per batch (default: 4)
DELAY=${3:-60}            # Delay between batches in seconds (default: 60)

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
cd "$BASE_DIR"

echo "================================================================================"
echo "BATCH SUBMISSION FOR PEAK-GENE LINKAGE"
echo "================================================================================"
echo ""
echo "Mode:          $MODE"
echo "Batch size:    $BATCH_SIZE tasks"
echo "Delay:         ${DELAY}s between batches"
echo ""
echo "This reduces concurrent conda environment access on BeeGFS"
echo ""
echo "================================================================================"
echo ""

################################################################################
# Function: Submit L2 Analysis in Batches
################################################################################

submit_L2() {
    echo "=== Submitting L2 Analysis (20 tasks) ==="
    echo ""

    TOTAL_TASKS=20
    CURRENT=0
    JOB_IDS=()

    while [ $CURRENT -lt $TOTAL_TASKS ]; do
        END=$((CURRENT + BATCH_SIZE - 1))
        if [ $END -ge $TOTAL_TASKS ]; then
            END=$((TOTAL_TASKS - 1))
        fi

        echo "Submitting batch: tasks $CURRENT-$END"

        JOB_ID=$(sbatch --array=$CURRENT-$END run_peak_gene_linkage_all.sh | awk '{print $4}')
        JOB_IDS+=($JOB_ID)

        echo "  Job ID: $JOB_ID"
        echo "  Status: $(squeue -j $JOB_ID -h | wc -l) tasks queued/running"

        CURRENT=$((END + 1))

        if [ $CURRENT -lt $TOTAL_TASKS ]; then
            echo "  Waiting ${DELAY}s before next batch..."
            sleep $DELAY
            echo ""
        fi
    done

    echo ""
    echo "✓ L2 submission complete"
    echo "  Job IDs: ${JOB_IDS[@]}"
    echo ""
}

################################################################################
# Function: Submit L1 Analysis in Batches
################################################################################

submit_L1() {
    echo "=== Submitting L1 Analysis (13 tasks) ==="
    echo ""

    TOTAL_TASKS=13
    CURRENT=0
    JOB_IDS=()

    while [ $CURRENT -lt $TOTAL_TASKS ]; do
        END=$((CURRENT + BATCH_SIZE - 1))
        if [ $END -ge $TOTAL_TASKS ]; then
            END=$((TOTAL_TASKS - 1))
        fi

        echo "Submitting batch: tasks $CURRENT-$END"

        JOB_ID=$(sbatch --array=$CURRENT-$END run_peak_gene_linkage_all_L1.sh | awk '{print $4}')
        JOB_IDS+=($JOB_ID)

        echo "  Job ID: $JOB_ID"
        echo "  Status: $(squeue -j $JOB_ID -h | wc -l) tasks queued/running"

        CURRENT=$((END + 1))

        if [ $CURRENT -lt $TOTAL_TASKS ]; then
            echo "  Waiting ${DELAY}s before next batch..."
            sleep $DELAY
            echo ""
        fi
    done

    echo ""
    echo "✓ L1 submission complete"
    echo "  Job IDs: ${JOB_IDS[@]}"
    echo ""
}

################################################################################
# Main Execution
################################################################################

START_TIME=$(date +%s)

case $MODE in
    L1)
        submit_L1
        ;;
    L2)
        submit_L2
        ;;
    both)
        submit_L2
        echo "Waiting ${DELAY}s before starting L1 submission..."
        sleep $DELAY
        echo ""
        submit_L1
        ;;
    *)
        echo "ERROR: Invalid mode '$MODE'"
        echo "Usage: $0 [L1|L2|both] [batch_size] [delay_seconds]"
        exit 1
        ;;
esac

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
ELAPSED_MIN=$((ELAPSED / 60))

echo "================================================================================"
echo "BATCH SUBMISSION COMPLETE"
echo "================================================================================"
echo ""
echo "Total time: ${ELAPSED}s (${ELAPSED_MIN} min)"
echo ""
echo "Monitor progress:"
echo "  squeue -u $USER"
echo ""
echo "Check for validation failures:"
echo "  grep -l 'ERROR: R environment validation failed' logs/run_peak_gene_linkage_all*.out"
echo ""
echo "================================================================================"
