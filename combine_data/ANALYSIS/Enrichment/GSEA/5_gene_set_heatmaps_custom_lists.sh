#!/bin/bash
#SBATCH --job-name=gene_set_heatmaps_custom_lists
#SBATCH --output=logs/gene_set_heatmaps_custom_lists.out
#SBATCH --error=logs/gene_set_heatmaps_custom_lists.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Gene Set Scoring and Spatial Visualization Pipeline
# This script runs the gene set analysis for spatial transcriptomics data

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

# Set the working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/ANALYSIS/Enrichment/GSEA


# Run the analysis
echo "Running gene set scoring"
python gene_set_heatmaps_custom_lists.py