#!/bin/bash
#SBATCH --job-name=test_env
#SBATCH --output=logs/test_env_%a_%N.out
#SBATCH --error=logs/test_env_%a_%N.err
#SBATCH --array=1-50%10  # Test up to 50 tasks, max 10 concurrent (spreads across nodes)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:05:00
#SBATCH --partition=workq
#SBATCH --account=kubacki.michal
#SBATCH --spread-job  # Force distribution across different nodes

source ~/.bashrc
conda activate sc-chromatin2

echo "================================================================================"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
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
    echo "✓ Task $SLURM_ARRAY_TASK_ID on $SLURMD_NODENAME PASSED"
else
    echo "✗ Task $SLURM_ARRAY_TASK_ID on $SLURMD_NODENAME FAILED"
fi

echo "================================================================================"
