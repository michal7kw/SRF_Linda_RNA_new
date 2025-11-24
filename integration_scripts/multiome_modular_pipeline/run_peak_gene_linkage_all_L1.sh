#!/bin/bash

################################################################################
# SLURM Array Job: Peak-Gene Linkage for All L1 Cell Types
################################################################################
#
# This script runs peak-gene linkage analysis for ALL available L1 cell types
# in parallel using SLURM array jobs.
#
# Includes all cell type × genotype combinations that have sufficient cells.
#
# USAGE:
#   sbatch run_peak_gene_linkage_all_L1.sh
#
# OUTPUT:
#   signac_results_L1/celltype_results/peak_gene_links/{CELLTYPE}_{GENOTYPE}_peak_gene_links.csv
#
# TIME ESTIMATE: 3-6 hours (all jobs in parallel)
#
################################################################################

#SBATCH --job-name=all_L1_linkage
#SBATCH --output=logs/run_peak_gene_linkage_all_L1_%a.out
#SBATCH --error=logs/run_peak_gene_linkage_all_L1_%a.err
#SBATCH --array=0-12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --partition=workq
#SBATCH --account=kubacki.michal
# #SBATCH --nodelist=srcn02,srcn04,srcn05,srcn10,srcn14,srcn17  # Only use specific nodes (can cause scheduling delays)
# Note: Comment out nodelist to allow flexible scheduling. Environment validation below will catch issues.

################################################################################
# Cell Type × Genotype Combinations
################################################################################
#
# Array Task ID Mapping (13 L1 broad categories):
#
# Cell types with EMX1 genotype (5 combinations):
#   0: Unknown Emx1
#   1: Excitatory Emx1
#   2: Oligo Emx1
#   3: Astrocytes Emx1
#   4: GABA Emx1
#
# Cell types with NESTIN genotype (8 combinations):
#   5: Unknown Nestin
#   6: Excitatory Nestin
#   7: Immune Nestin
#   8: Oligo Nestin
#   9: Astrocytes Nestin
#  10: GABA Nestin
#  11: Vascular Nestin
#  12: Progenitors Nestin
#
# NOTE: Ependymal excluded due to low cell counts
#
declare -a COMBINATIONS=(
    "Unknown Emx1"           # 0
    "Excitatory Emx1"        # 1
    "Oligo Emx1"             # 2
    "Astrocytes Emx1"        # 3
    "GABA Emx1"              # 4
    "Unknown Nestin"         # 5
    "Excitatory Nestin"      # 6
    "Immune Nestin"          # 7
    "Oligo Nestin"           # 8
    "Astrocytes Nestin"      # 9
    "GABA Nestin"            # 10
    "Vascular Nestin"        # 11
    "Progenitors Nestin"     # 12
)

################################################################################
# Parse Current Task
################################################################################

TASK_INFO="${COMBINATIONS[$SLURM_ARRAY_TASK_ID]}"
read -r CELLTYPE GENOTYPE <<< "$TASK_INFO"

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

# Set number of threads
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Validate critical R packages before starting (comprehensive check)
echo "Validating R environment..."
Rscript -e "library(Signac); library(Seurat); library(BSgenome.Mmusculus.UCSC.mm10); library(GenomicRanges); library(rtracklayer); library(MatrixGenerics); library(S4Arrays); library(DelayedArray); library(dplyr); cat('✓ All required packages loaded successfully\n')" 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: R environment validation failed"
    echo "Please reinstall corrupted packages:"
    echo "  conda activate sc-chromatin2"
    echo "  conda install -c bioconda --force-reinstall bioconductor-s4arrays bioconductor-rtracklayer bioconductor-matrixgenerics bioconductor-delayedarray"
    exit 1
fi

echo ""

################################################################################
# Job Information
################################################################################

echo "================================================================================"
echo "  PEAK-GENE LINKAGE: ALL L1 CELL TYPES"
echo "================================================================================"
echo ""
echo "Array Job ID:   ${SLURM_ARRAY_JOB_ID}"
echo "Task ID:        ${SLURM_ARRAY_TASK_ID}"
echo "Node:           ${SLURMD_NODENAME}"
echo "CPUs:           ${SLURM_CPUS_PER_TASK}"
echo "Memory:         ${SLURM_MEM_PER_NODE}MB"
echo "Start time:     $(date)"
echo ""
echo "Target:"
echo "  Cell type:    ${CELLTYPE}"
echo "  Genotype:     ${GENOTYPE}"
echo ""
echo "================================================================================"
echo ""

################################################################################
# Run Analysis
################################################################################

START_TIME=$(date +%s)

cd ${BASE_DIR}

Rscript signac_03b_peak_gene_linkage_single_celltype_L1.R "${CELLTYPE}" "${GENOTYPE}" 2>&1

EXIT_CODE=$?

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

################################################################################
# Summary
################################################################################

echo ""
echo "================================================================================"

if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ TASK ${SLURM_ARRAY_TASK_ID} COMPLETE"
    echo "================================================================================"
    echo ""
    echo "Cell type:      ${CELLTYPE}"
    echo "Genotype:       ${GENOTYPE}"
    echo "Elapsed time:   ${HOURS}h ${MINUTES}m ${SECONDS}s"
    echo ""
    echo "Output file:"
    echo "  signac_results_L1/celltype_results/peak_gene_links/${CELLTYPE}_${GENOTYPE}_peak_gene_links.csv"
else
    echo "✗ TASK ${SLURM_ARRAY_TASK_ID} FAILED"
    echo "================================================================================"
    echo ""
    echo "Exit code:      ${EXIT_CODE}"
    echo "Check log file for details"
    exit $EXIT_CODE
fi

echo ""
echo "End time:       $(date)"
echo "================================================================================"
