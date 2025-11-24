#!/bin/bash

################################################################################
# Check Node Health Test Results
################################################################################
#
# Analyzes results from test_all_nodes.sh and provides a summary
#
# Usage:
#   ./check_node_health.sh
#
################################################################################

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"

echo "================================================================================"
echo "NODE HEALTH TEST RESULTS"
echo "================================================================================"
echo ""

# Check if any test files exist
if ! ls ${BASE_DIR}/logs/node_test_*.out >/dev/null 2>&1; then
    echo "ERROR: No test result files found"
    echo "Please run: ./test_all_nodes.sh first"
    exit 1
fi

# Arrays to store results
HEALTHY_NODES=()
FAILED_NODES=()
PENDING_NODES=()

# Analyze each node test result
for LOG_FILE in ${BASE_DIR}/logs/node_test_*.out; do
    NODE=$(basename "$LOG_FILE" | sed 's/node_test_\(.*\)\.out/\1/')

    # Check if job completed
    if ! grep -q "Exit code:" "$LOG_FILE" 2>/dev/null; then
        PENDING_NODES+=($NODE)
        continue
    fi

    # Check if passed or failed
    if grep -q "PASSED" "$LOG_FILE"; then
        HEALTHY_NODES+=($NODE)
    else
        FAILED_NODES+=($NODE)
    fi
done

# Print summary
NUM_HEALTHY=${#HEALTHY_NODES[@]}
NUM_FAILED=${#FAILED_NODES[@]}
NUM_PENDING=${#PENDING_NODES[@]}
NUM_TOTAL=$((NUM_HEALTHY + NUM_FAILED + NUM_PENDING))

echo "Summary:"
echo "  Total nodes tested:    $NUM_TOTAL"
echo "  ✓ Healthy nodes:       $NUM_HEALTHY"
echo "  ✗ Failed nodes:        $NUM_FAILED"
echo "  ⏳ Pending/Running:     $NUM_PENDING"
echo ""

if [ $NUM_HEALTHY -gt 0 ]; then
    echo "================================================================================"
    echo "✓ HEALTHY NODES (Environment OK)"
    echo "================================================================================"
    for NODE in "${HEALTHY_NODES[@]}"; do
        echo "  $NODE"
    done
    echo ""
fi

if [ $NUM_FAILED -gt 0 ]; then
    echo "================================================================================"
    echo "✗ FAILED NODES (Corrupted Environment)"
    echo "================================================================================"
    for NODE in "${FAILED_NODES[@]}"; do
        echo "  $NODE"

        # Show error details
        LOG_FILE="${BASE_DIR}/logs/node_test_${NODE}.out"
        if grep -q "Error" "$LOG_FILE"; then
            echo "    Error: $(grep "Error:" "$LOG_FILE" | head -1 | sed 's/^[[:space:]]*//')"
        fi
    done
    echo ""
fi

if [ $NUM_PENDING -gt 0 ]; then
    echo "================================================================================"
    echo "⏳ PENDING NODES (Still Running or Not Started)"
    echo "================================================================================"
    for NODE in "${PENDING_NODES[@]}"; do
        echo "  $NODE"
    done
    echo ""
    echo "Note: Run this script again after jobs complete"
    echo ""
fi

# Generate SLURM nodelist configuration
if [ $NUM_HEALTHY -gt 0 ]; then
    echo "================================================================================"
    echo "RECOMMENDED SLURM CONFIGURATION"
    echo "================================================================================"
    echo ""

    # Create comma-separated nodelist
    NODELIST=$(IFS=,; echo "${HEALTHY_NODES[*]}")

    echo "Add this line to your SLURM scripts:"
    echo ""
    echo "  #SBATCH --nodelist=$NODELIST"
    echo ""

    # Save to file for easy reference
    echo "$NODELIST" > ${BASE_DIR}/healthy_nodes.txt
    echo "Saved healthy node list to: healthy_nodes.txt"
    echo ""
fi

if [ $NUM_FAILED -gt 0 ]; then
    echo "================================================================================"
    echo "CORRUPTED NODES - NEED ATTENTION"
    echo "================================================================================"
    echo ""
    echo "The following nodes have corrupted R packages:"
    for NODE in "${FAILED_NODES[@]}"; do
        echo "  - $NODE"
    done
    echo ""
    echo "To fix, contact cluster admin or exclude these nodes:"
    EXCLUDE_LIST=$(IFS=,; echo "${FAILED_NODES[*]}")
    echo "  #SBATCH --exclude=$EXCLUDE_LIST"
    echo ""
fi

echo "================================================================================"
echo ""
echo "For detailed logs, see: logs/node_test_*.out"
echo ""
