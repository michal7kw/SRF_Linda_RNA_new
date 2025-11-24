#!/bin/bash
#SBATCH --job-name=2_gsea_heatmap_plot
#SBATCH --output=logs/2_gsea_heatmap_plot.out
#SBATCH --error=logs/2_gsea_heatmap_plot.err
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

# GSEA Heatmap Plot Generator
#
# NOTE: The script automatically processes BOTH configurations:
#   - Default GSEA processing (gsea_analysis_between_conditions_default)
#   - Non-default GSEA processing (gsea_analysis_between_conditions)
#
# For EACH configuration, it generates TWO versions of each heatmap:
#   - Standard version WITH FDR annotations (NES; FDR values)
#   - Clean version WITHOUT FDR indications (filename suffix: _no_fdr)
#
# This script assumes that 'mamba' is available and the shell is configured
# to use it (e.g., by running 'conda init bash' or 'mamba init bash').

# Fail on any error
set -e

# Get the directory where the script is located
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/ANALYSIS/Enrichment/GSEA"

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

# Since the scripts are in the same directory as this .sh file, we'll cd there.
cd "$SCRIPT_DIR"

echo "NOTE: Processing BOTH configurations (default & non-default)"
echo "NOTE: Both versions (with and without FDR) will be generated for each"
echo ""

run_heatmap_script() {
    local cluster_type=$1
    local input_file=$2
    echo "  - Starting 2_gsea_heatmap_plot.py for cluster_type=${cluster_type}"
    python 2_gsea_heatmap_plot.py --cluster_type "${cluster_type}" --input_file "${input_file}"
}

# Run for both cluster types
run_heatmap_script "Mature_GC" "GSEA_Focus_Neuro_CellDeath.xlsx"
run_heatmap_script "Combined_GC" "GSEA_Focus_Neuro_CellDeath.xlsx"
run_heatmap_script "GABA" "GSEA_Focus_RNA_Splicing.xlsx"

echo "Waiting for all heatmap scripts to complete..."
wait

echo "==============================================="
echo "All GSEA heatmap scripts have completed!"
echo "✓ Processed BOTH configurations:"
echo "  - Default: results_from_raw/gsea_summary_visuals_default/"
echo "  - Non-default: results_from_raw/gsea_summary_visuals/"
echo "✓ Generated both FDR and no-FDR versions for each"
echo "==============================================="