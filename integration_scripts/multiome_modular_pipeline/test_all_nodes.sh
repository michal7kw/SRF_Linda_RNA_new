#!/bin/bash

################################################################################
# Test R Environment on ALL Available Nodes
################################################################################
#
# This script identifies all available nodes in the workq partition and
# submits a test job to each one individually.
#
# Usage:
#   ./test_all_nodes.sh
#
################################################################################

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"

# Create logs directory
mkdir -p ${BASE_DIR}/logs

echo "================================================================================"
echo "TESTING R ENVIRONMENT ON ALL NODES"
echo "================================================================================"
echo ""

# Get list of all available nodes in workq partition
echo "Querying available nodes in workq partition..."
NODES=$(sinfo -p workq -N -h -o "%N" | sort -u)

if [ -z "$NODES" ]; then
    echo "ERROR: No nodes found in workq partition"
    exit 1
fi

NODE_ARRAY=($NODES)
NUM_NODES=${#NODE_ARRAY[@]}

echo "Found $NUM_NODES nodes: ${NODE_ARRAY[@]}"
echo ""

# Create a temporary test script for direct node submission
cat > ${BASE_DIR}/test_single_node.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=test_node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:05:00
#SBATCH --partition=workq
#SBATCH --account=kubacki.michal

source ~/.bashrc
conda activate sc-chromatin2

echo "================================================================================"
echo "Node: $SLURMD_NODENAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "================================================================================"
echo ""
echo "Testing R packages..."
echo ""

Rscript -e "
# Load all packages that peak-gene linkage analysis needs
library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(rtracklayer)
library(MatrixGenerics)
library(S4Arrays)
library(DelayedArray)
library(dplyr)
cat('Node', Sys.info()['nodename'], 'is healthy\n')
" 2>&1

EXIT_CODE=$?

echo ""
echo "Exit code: $EXIT_CODE"
echo "================================================================================"

if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Node $SLURMD_NODENAME PASSED"
else
    echo "✗ Node $SLURMD_NODENAME FAILED"
fi

echo "================================================================================"
EOF

chmod +x ${BASE_DIR}/test_single_node.sh

# Submit one job per node
echo "Submitting test jobs..."
echo ""

SUBMITTED_JOBS=()

for NODE in "${NODE_ARRAY[@]}"; do
    echo "  Submitting test for node: $NODE"

    JOB_ID=$(sbatch --nodelist=$NODE \
                   --output=${BASE_DIR}/logs/node_test_${NODE}.out \
                   --error=${BASE_DIR}/logs/node_test_${NODE}.err \
                   ${BASE_DIR}/test_single_node.sh | awk '{print $4}')

    SUBMITTED_JOBS+=($JOB_ID)

    # Small delay to avoid overwhelming scheduler
    sleep 0.5
done

echo ""
echo "================================================================================"
echo "SUBMITTED $NUM_NODES TEST JOBS"
echo "================================================================================"
echo ""
echo "Job IDs: ${SUBMITTED_JOBS[@]}"
echo ""
echo "Monitor progress with:"
echo "  squeue -u $USER"
echo ""
echo "Check results when complete with:"
echo "  ${BASE_DIR}/check_node_health.sh"
echo ""
echo "================================================================================"
