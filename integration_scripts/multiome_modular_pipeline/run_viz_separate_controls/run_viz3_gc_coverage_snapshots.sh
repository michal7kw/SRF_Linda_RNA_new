#!/bin/bash

################################################################################
# SLURM Submission Script for Granule Cell Coverage Snapshots
################################################################################

#SBATCH --job-name=gc_coverage
#SBATCH --output=logs/run_viz3_gc_coverage_snapshots.out
#SBATCH --error=logs/run_viz3_gc_coverage_snapshots.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=workq
#SBATCH --account=kubacki.michal

################################################################################
# Configuration
################################################################################

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"

# Create logs directory
mkdir -p ${BASE_DIR}/logs

# Load conda environment
source ~/.bashrc
conda activate sc-chromatin2

# CRITICAL: Tell R's reticulate to use conda Python
export RETICULATE_PYTHON=$(which python)
export RETICULATE_CONDA=$(which conda)

# Set number of threads for R
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "================================================================================"
echo "  GRANULE CELL COVERAGE SNAPSHOTS - SPECIFIC GENES"
echo "================================================================================"
echo ""
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo "Start time: $(date)"
echo ""
echo "This script generates pseudobulk chromatin accessibility snapshots for:"
echo "  - Cell types: Immature GC + Mature GC (combined)"
echo "  - Genes: Rbfox1, Gli3, Cdk5"
echo "  - Comparisons: Nestin Ctrl-Mut, Emx1 Ctrl-Mut (separate)"
echo "  - Region: Gene body ± 100kb"
echo "  - Shows: ATAC accessibility + RNA expression"
echo ""
echo "================================================================================"
echo ""

################################################################################
# Run Coverage Snapshot Generation
################################################################################

START_TIME=$(date +%s)

Rscript ${BASE_DIR}/scripts/separate_version/create_gc_coverage_snapshots.R 2>&1

EXIT_CODE=$?

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "✗ Coverage snapshot generation failed with exit code ${EXIT_CODE}"
    echo "  Check logs for errors"
    exit 1
fi

echo ""
echo "✓ Coverage snapshot generation complete (elapsed time: ${ELAPSED}s)"
echo ""

################################################################################
# Final Summary
################################################################################

echo "================================================================================"
echo "COVERAGE SNAPSHOTS COMPLETE!"
echo "================================================================================"
echo ""
echo "End time: $(date)"
echo ""
echo "Results saved to: ${BASE_DIR}/signac_results/plots/gc_coverage_snapshots/"
echo ""
echo "Output files:"
echo "  - GC_coverage_snapshots_specific_genes.pdf (main output - all genes)"
echo "  - GC_coverage_Rbfox1.pdf (Rbfox1 only)"
echo "  - GC_coverage_Gli3.pdf (Gli3 only)"
echo "  - GC_coverage_Cdk5.pdf (Cdk5 only)"
echo "  - GC_coverage_summary_all_genes.pdf (summary by genotype)"
echo ""
echo "Each plot shows:"
echo "  • ATAC-seq coverage tracks (pseudobulk by condition)"
echo "  • RNA expression levels"
echo "  • Peak annotations"
echo "  • Gene structure"
echo "  • Separate comparisons for Nestin and Emx1 genotypes"
echo ""
echo "================================================================================"
