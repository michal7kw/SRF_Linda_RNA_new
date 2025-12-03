#!/bin/bash
#SBATCH --job-name=0_CREs_splicing_differential_all_celltypes
#SBATCH --output=logs/0_CREs_splicing_differential_all_celltypes_%j.log
#SBATCH --error=logs/0_CREs_splicing_differential_all_celltypes_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=workq

# =============================================================================
# MASTER SCRIPT: Differential Splicing CRE Analysis (All Cell Types)
# =============================================================================

set -e

SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_differential_all_celltypes"
cd "${SCRIPT_DIR}"

mkdir -p logs

echo "=============================================================================="
echo "STARTING PIPELINE: Differential Splicing CREs (All Cell Types)"
echo "=============================================================================="
echo "Date: $(date)"
echo ""

# Activate environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate rna_seq_analysis_deep

# -----------------------------------------------------------------------------
# Step 1: Extract CREs
# -----------------------------------------------------------------------------
echo "Step 1: Extracting CREs..."
python 1_extract_differential_candidates.py

# -----------------------------------------------------------------------------
# Step 2: Compute Matrices
# -----------------------------------------------------------------------------
echo "Step 2: Computing Signal Matrices..."
# We run this directly (not sbatch) to ensure sequential execution for this master script
# But we use the sbatch script as a source
bash 2_compute_signal_matrices.sh

# -----------------------------------------------------------------------------
# Step 3: Analyze Differential CREs
# -----------------------------------------------------------------------------
echo "Step 3: Analyzing Differential CREs..."
# Using min_fc=3.0 as requested
python 3_analyze_differential_cres.py --min-fc 1.5 --min-signal 2.0

echo ""
echo "=============================================================================="
echo "PIPELINE COMPLETE"
echo "=============================================================================="
