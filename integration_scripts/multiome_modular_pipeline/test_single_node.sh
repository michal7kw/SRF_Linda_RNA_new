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
