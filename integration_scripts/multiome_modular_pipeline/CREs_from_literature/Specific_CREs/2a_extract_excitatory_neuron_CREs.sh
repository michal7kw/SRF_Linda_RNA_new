#!/bin/bash
#SBATCH --job-name=2a_extract_excitatory_neuron_CREs
#SBATCH --output=logs/2a_extract_excitatory_neuron_CREs.log
#SBATCH --error=logs/2a_extract_excitatory_neuron_CREs.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=workq

################################################################################
# Extract Excitatory Neuron CREs
#
# This script extracts CREs specific to Excitatory neurons.
# Part of the Specific CREs workflow.
################################################################################

# Activate conda environment with pandas
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

echo "========================================================================"
echo "EXTRACT EXCITATORY NEURON CREs"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Run the Python script
python 2a_extract_excitatory_neuron_CREs.py

EXIT_CODE=$?

echo ""
echo "========================================================================"

if [ $EXIT_CODE -eq 0 ]; then
    echo "EXTRACTION COMPLETE!"
    if [ -f "output/excitatory_neuron_CREs.bed" ]; then
        N_CRES=$(wc -l < output/excitatory_neuron_CREs.bed)
        echo "Created output/excitatory_neuron_CREs.bed with $N_CRES CREs"
    fi
else
    echo "EXTRACTION FAILED!"
    exit 1
fi
