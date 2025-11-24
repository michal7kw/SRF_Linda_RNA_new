#!/bin/bash

################################################################################
# SLURM Submission Script for Signac Multiome Analysis Pipeline
################################################################################

#SBATCH --job-name=run_signac_pipeline
#SBATCH --output=logs/run_signac_pipeline.out
#SBATCH --error=logs/run_signac_pipeline.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=workq
#SBATCH --account=kubacki.michal

# It runs:
# - signac_01_load_and_integrate.R
# - signac_02_process_and_reduce.R
# - signac_03_celltype_analysis.R
# - signac_04_visualization.R
# - signac_07_POST_PROCESSING_FILTERS.R


################################################################################
# Configuration
################################################################################

# Base directory
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"

# Create logs directory
mkdir -p ${BASE_DIR}/logs

# Load conda environment
source ~/.bashrc

conda activate sc-chromatin2

# CRITICAL: Tell R's reticulate to use conda Python (not system Python)
export RETICULATE_PYTHON=$(which python)
export RETICULATE_CONDA=$(which conda)
echo "Python for R: ${RETICULATE_PYTHON}"
echo "Conda for R: ${RETICULATE_CONDA}"


# Set number of threads for R
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "=================================================================================="
echo "  SIGNAC MULTIOME ANALYSIS PIPELINE"
echo "=================================================================================="
echo ""
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo "Start time: $(date)"
echo ""
echo "=================================================================================="
echo ""

################################################################################
# Step 1: Load and Integrate Data
################################################################################

echo "================================================================================"
echo "STEP 1: Loading and Integrating RNA and ATAC Data"
echo "================================================================================"
echo ""

START_TIME=$(date +%s)

Rscript ${BASE_DIR}/signac_01_load_and_integrate.R 2>&1

EXIT_CODE=$?

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "✗ Step 1 failed with exit code ${EXIT_CODE}"
    echo "  Check logs for errors"
    exit 1
fi

echo ""
echo "✓ Step 1 complete (elapsed time: ${ELAPSED}s)"
echo ""

################################################################################
# Step 2: Processing and Dimensionality Reduction
################################################################################

echo "================================================================================"
echo "STEP 2: Processing and Dimensionality Reduction"
echo "================================================================================"
echo ""

START_TIME=$(date +%s)

Rscript ${BASE_DIR}/signac_02_process_and_reduce.R 2>&1

EXIT_CODE=$?

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "✗ Step 2 failed with exit code ${EXIT_CODE}"
    echo "  Check logs for errors"
    exit 1
fi

echo ""
echo "✓ Step 2 complete (elapsed time: ${ELAPSED}s)"
echo ""

################################################################################
# Step 3: Cell-Type-Specific Differential Analysis
################################################################################

echo "================================================================================"
echo "STEP 3: Cell-Type-Specific Differential Analysis"
echo "================================================================================"
echo ""

START_TIME=$(date +%s)

Rscript ${BASE_DIR}/signac_03_celltype_analysis.R 2>&1

EXIT_CODE=$?

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "✗ Step 3 failed with exit code ${EXIT_CODE}"
    echo "  Check logs for errors"
    exit 1
fi

echo ""
echo "✓ Step 3 complete (elapsed time: ${ELAPSED}s)"
echo ""

################################################################################
# Step 4: Advanced Visualization
################################################################################

echo "================================================================================"
echo "STEP 4: Advanced Visualization"
echo "================================================================================"
echo ""

START_TIME=$(date +%s)

Rscript ${BASE_DIR}/signac_04_visualization.R 2>&1

EXIT_CODE=$?

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "✗ Step 4 failed with exit code ${EXIT_CODE}"
    echo "  Check logs for errors"
    exit 1
fi

echo ""
echo "✓ Step 4 complete (elapsed time: ${ELAPSED}s)"
echo ""

################################################################################
# Step 5: Post-Processing Filters
################################################################################

echo "================================================================================"
echo "STEP 5: Post-Processing Filters (L1)"
echo "================================================================================"
echo ""

START_TIME=$(date +%s)

Rscript ${BASE_DIR}/signac_07_POST_PROCESSING_FILTERS.R 2>&1

EXIT_CODE=$?

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "✗ Step 5 failed with exit code ${EXIT_CODE}"
    echo "  Check logs for errors"
    exit 1
fi

echo ""
echo "✓ Step 5 complete (elapsed time: ${ELAPSED}s)"
echo ""

################################################################################
# Final Summary
################################################################################

echo "================================================================================"
echo "PIPELINE COMPLETE!"
echo "================================================================================"
echo ""
echo "End time: $(date)"
echo ""
echo "Results saved to: ${BASE_DIR}/signac_results/"
echo ""
echo "Key output files:"
echo "  - signac_results/integrated_seurat_processed.rds"
echo "  - signac_results/celltype_results/summary_statistics.csv"
echo "  - signac_results/celltype_results/DEG/*.csv"
echo "  - signac_results/celltype_results/DA/*.csv"
echo "  - signac_results/celltype_results/peak_gene_links/*.csv"
echo "  - signac_results/celltype_results/filtered/ (post-processed high-confidence results)"
echo "  - signac_results/plots/*.pdf"
echo ""
echo "L1 Cell Categories Analyzed:"
echo "  - Excitatory, GABA, Oligo, Astrocytes, Immune,"
echo "  - Vascular, Progenitors, Unknown, Ependymal"
echo ""
echo "=================================================================================="
