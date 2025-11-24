#!/bin/bash

################################################################################
# SLURM Array Job: Peak-Gene Linkage for All L2 Cell Types
################################################################################
#
# This script runs peak-gene linkage analysis for ALL available L2 cell types
# in parallel using SLURM array jobs.
#
# Includes all cell type × genotype combinations that have sufficient cells.
#
# USAGE:
#   sbatch run_peak_gene_linkage_all.sh
#
# OUTPUT:
#   signac_results/celltype_results/peak_gene_links/{CELLTYPE}_{GENOTYPE}_peak_gene_links.csv
#
# TIME ESTIMATE: 4-8 hours (all jobs in parallel)
#
################################################################################

#SBATCH --job-name=all_linkage
#SBATCH --output=logs/run_peak_gene_linkage_all_%a.out
#SBATCH --error=logs/run_peak_gene_linkage_all_%a.err
#SBATCH --array=0-19
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
# Array Task ID Mapping:
#
# Cell types with BOTH genotypes (12 combinations):
#   0-1:   Astrocytes (Nestin, Emx1)
#   2-3:   CA1 (Nestin, Emx1)
#   4-5:   CTX (Nestin, Emx1)
#   6-7:   ENT CTX (Nestin, Emx1)
#   8-9:   Mature GC (Nestin, Emx1)
#   10-11: rest (Nestin, Emx1)
#
# Cell types with NESTIN only (8 combinations):
#   12: CA3 Nestin
#   13: Immature GC Nestin
#   14: Microglia Nestin
#   15: MOL Nestin
#   16: NFOL Nestin
#   17: OPC Nestin
#   18: PV Nestin
#   19: SST Nestin
#
# NOTE: Some cell types like Subiculum are excluded due to low cell counts
#
declare -a COMBINATIONS=(
    # Both genotypes
    "Astrocytes|Nestin|450|30000"        # 0
    "Astrocytes|Emx1|280|30000"          # 1
    "CA1|Nestin|1200|35000"              # 2
    "CA1|Emx1|650|35000"                 # 3
    "CTX|Nestin|2500|40000"              # 4
    "CTX|Emx1|1100|40000"                # 5
    "ENT CTX|Nestin|800|32000"           # 6
    "ENT CTX|Emx1|400|32000"             # 7
    "Mature GC|Nestin|3500|45000"        # 8
    "Mature GC|Emx1|1800|45000"          # 9
    "rest|Nestin|200|25000"              # 10
    "rest|Emx1|150|25000"                # 11
    # Nestin only
    "CA3|Nestin|600|33000"               # 12
    "Immature GC|Nestin|1500|38000"      # 13
    "Microglia|Nestin|250|28000"         # 14
    "MOL|Nestin|1800|35000"              # 15 - Myelinating oligodendrocytes
    "NFOL|Nestin|900|32000"              # 16 - Newly formed oligodendrocytes
    "OPC|Nestin|400|30000"               # 17 - Oligodendrocyte precursor cells
    "PV|Nestin|500|31000"                # 18 - Parvalbumin interneurons
    "SST|Nestin|350|29000"               # 19 - Somatostatin interneurons
)

################################################################################
# Parse Current Task
################################################################################

TASK_INFO="${COMBINATIONS[$SLURM_ARRAY_TASK_ID]}"
IFS='|' read -r CELLTYPE GENOTYPE N_CELLS PEAKS_ESTIMATE <<< "$TASK_INFO"

################################################################################
# Configuration
################################################################################

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"

# Create logs directory
mkdir -p ${BASE_DIR}/logs

################################################################################
# Environment Setup
################################################################################

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
echo "  PEAK-GENE LINKAGE: ALL L2 CELL TYPES"
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
echo "  Expected cells: ~${N_CELLS}"
echo "  Expected peaks: ~${PEAKS_ESTIMATE}"
echo ""
echo "================================================================================"
echo ""

################################################################################
# Run Analysis
################################################################################

START_TIME=$(date +%s)

cd ${BASE_DIR}

Rscript signac_03b_peak_gene_linkage_single_celltype.R "${CELLTYPE}" "${GENOTYPE}" 2>&1

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
    echo "  signac_results/celltype_results/peak_gene_links/${CELLTYPE}_${GENOTYPE}_peak_gene_links.csv"
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
