#!/bin/bash
#SBATCH --job-name=3_visualize_DEGs_encode
#SBATCH --output=logs/3_visualize_deg_comparisons.log
#SBATCH --error=logs/3_visualize_deg_comparisons.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=workq

echo "Starting Step 3: Create custom comparison visualizations"
echo "Date: $(date)"

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_GABA_DEGs_encode"

# Activate conda environment (sc-chromatin2 has matplotlib, seaborn, scipy)
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate sc-chromatin2

# Run visualization script
python 3_visualize_deg_comparisons.py  

echo "Step 3 completed: $(date)"
