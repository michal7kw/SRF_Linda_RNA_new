#!/bin/bash

################################################################################
# Submit Peak-Gene Linkage Jobs with Node Preference
################################################################################
#
# This wrapper submits jobs with preferred nodes but falls back to any node
# if preferred ones are busy
#
################################################################################

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"

echo "================================================================================"
echo "SUBMITTING PEAK-GENE LINKAGE JOBS"
echo "================================================================================"
echo ""

# Check if preferred nodes are available
PREFERRED_NODES="srcn03,srcn04"
AVAILABLE=$(sinfo -n $PREFERRED_NODES -h -t idle,mix -o "%N" 2>/dev/null | wc -l)

if [ $AVAILABLE -gt 0 ]; then
    echo "✓ Preferred nodes available: $PREFERRED_NODES"
    echo "  Submitting with nodelist constraint..."

    # Submit with nodelist
    JOB_ID=$(sbatch --nodelist=$PREFERRED_NODES run_peak_gene_linkage_all.sh | awk '{print $4}')
    echo "  Job submitted: $JOB_ID"
else
    echo "⚠ Preferred nodes busy: $PREFERRED_NODES"
    echo "  Submitting with flexible scheduling (will use any available node)..."
    echo "  Environment validation will ensure node health."

    # Submit without nodelist constraint
    JOB_ID=$(sbatch run_peak_gene_linkage_all.sh | awk '{print $4}')
    echo "  Job submitted: $JOB_ID"
fi

echo ""
echo "Monitor with: squeue -j $JOB_ID"
echo "================================================================================"
