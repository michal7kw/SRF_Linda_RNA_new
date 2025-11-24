#!/bin/bash
#SBATCH --job-name=create_heatmap_collages
#SBATCH --output=logs/create_heatmap_collages.out
#SBATCH --error=logs/create_heatmap_collages.err
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
echo "Running create_heatmap_collages"
python create_heatmap_collages.py