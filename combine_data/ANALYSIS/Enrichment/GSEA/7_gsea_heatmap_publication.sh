#!/bin/bash
#SBATCH --job-name=gsea_pub_heatmap
#SBATCH --output=logs/gsea_pub_heatmap.out
#SBATCH --error=logs/gsea_pub_heatmap.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

# Publication-Ready GSEA Heatmap Generator
# Creates compact, publication-quality heatmaps from GSEA results
#
# NOTE: The script automatically generates TWO versions of each heatmap:
#   - Standard version WITH FDR significance stars (*, **, ***)
#   - Clean version WITHOUT FDR indications (filename suffix: _no_fdr)

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

# Set working directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/ANALYSIS/Enrichment/GSEA"
cd "$SCRIPT_DIR"

# Create logs directory if it doesn't exist
mkdir -p logs

echo "==============================================="
echo "Publication-Ready GSEA Heatmap Generator"
echo "==============================================="
echo "Start time: $(date)"
echo ""

# Configuration
MAX_TERMS=15
FDR_THRESHOLD=0.25

# Cluster types to process
CLUSTER_TYPES=("Mature_GC" "Combined_GC" "GABA")

# Option 1: Create standard heatmaps with specified max_terms
echo "Creating standard publication heatmaps..."
echo "Input files will be auto-detected based on cluster type:"
echo "  - Mature_GC/Combined_GC: GSEA_Focus_Neuro_CellDeath.csv"
echo "  - GABA: GSEA_Focus_RNA_Splicing.csv"
echo ""
echo "NOTE: Both versions (with and without FDR) will be generated for each cluster"
echo ""

for CLUSTER in "${CLUSTER_TYPES[@]}"; do
    echo "Processing cluster: $CLUSTER"
    python 7_gsea_heatmap_publication.py \
        --cluster_type "$CLUSTER" \
        --max_terms "$MAX_TERMS" \
        --fdr_threshold "$FDR_THRESHOLD"

    if [ $? -eq 0 ]; then
        echo "  ✓ Successfully created heatmaps for $CLUSTER (both with and without FDR)"
    else
        echo "  ✗ Failed to create heatmap for $CLUSTER"
    fi
    echo ""
done

# Option 2: Create multiple variants (uncomment to use)
# This creates top10, balanced, top15, and consistent versions
# NOTE: When using --create_variants, BOTH versions (with/without FDR) are generated for ALL variants
#
# echo "Creating variant heatmaps..."
# echo ""
#
# for CLUSTER in "${CLUSTER_TYPES[@]}"; do
#     echo "Processing cluster variants: $CLUSTER"
#     python 7_gsea_heatmap_publication.py \
#         --cluster_type "$CLUSTER" \
#         --create_variants
#
#     if [ $? -eq 0 ]; then
#         echo "  ✓ Successfully created variant heatmaps for $CLUSTER (all with and without FDR)"
#     else
#         echo "  ✗ Failed to create variant heatmaps for $CLUSTER"
#     fi
#     echo ""
# done

echo "==============================================="
echo "Heatmap generation completed!"
echo "✓ Generated both FDR and no-FDR versions"
echo "✓ Check output: results_from_raw/gsea_summary_visuals*/"
echo "End time: $(date)"
echo "==============================================="

# Deactivate conda environment
conda deactivate
