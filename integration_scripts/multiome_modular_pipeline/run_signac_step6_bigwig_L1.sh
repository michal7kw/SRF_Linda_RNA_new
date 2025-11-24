#!/bin/bash

################################################################################
# SLURM Submission Script for BigWig Track Generation (L1 VERSION)
################################################################################

#SBATCH --job-name=signac_bigwig_L1
#SBATCH --output=logs/signac_step6_bigwig_L1.out
#SBATCH --error=logs/signac_step6_bigwig_L1.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
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
echo "  SIGNAC MULTIOME ANALYSIS - STEP 6: BIGWIG TRACK GENERATION (L1 VERSION)"
echo "================================================================================"
echo ""
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo "Start time: $(date)"
echo ""
echo "================================================================================"
echo ""

################################################################################
# Run BigWig Generation
################################################################################

START_TIME=$(date +%s)

Rscript ${BASE_DIR}/signac_06_create_bigwig_tracks_L1.R 2>&1

EXIT_CODE=$?

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "✗ BigWig generation failed with exit code ${EXIT_CODE}"
    echo "  Check logs for errors"
    exit 1
fi

echo ""
echo "✓ BigWig generation complete (elapsed time: ${ELAPSED}s)"
echo ""

################################################################################
# Final Summary
################################################################################

echo "================================================================================"
echo "BIGWIG GENERATION COMPLETE (L1 VERSION)!"
echo "================================================================================"
echo ""
echo "End time: $(date)"
echo ""
echo "Results saved to: ${BASE_DIR}/signac_results/bigwig_tracks_L1/"
echo ""
echo "Key output files:"
echo "  - signac_results/bigwig_tracks_L1/by_sample/*.bw"
echo "  - signac_results/bigwig_tracks_L1/by_celltype/*.bw  (9 L1 lineages)"
echo "  - signac_results/bigwig_tracks_L1/igv_session.xml"
echo "  - signac_results/bigwig_tracks_L1/README.md"
echo ""
echo "L1 lineages: Excitatory, GABA, Oligo, Astrocytes, Immune, Vascular,"
echo "             Progenitors, Unknown, Ependymal"
echo ""
echo "To visualize in IGV:"
echo "  1. Download IGV from: https://software.broadinstitute.org/software/igv/"
echo "  2. Open IGV and load the session file (igv_session.xml)"
echo "  3. Navigate to your genes of interest"
echo ""
echo "================================================================================"
