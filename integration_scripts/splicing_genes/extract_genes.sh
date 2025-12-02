#!/bin/bash
#SBATCH --job-name=extract_genes
#SBATCH --output=./logs/extract_genes.out
#SBATCH --error=./logs/extract_genes.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=32G
#SBATCH --time=00:04:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Load conda environment
source ~/.bashrc

conda activate bioinf

# Set working directory to script location
WORKDIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/CREs_splicing_genes_paper"
cd "$WORKDIR"

python extract_genes.py