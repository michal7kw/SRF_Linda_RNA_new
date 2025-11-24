#!/bin/bash
#SBATCH --job-name=gsea_excitatory_L1
#SBATCH --output=logs/gsea_excitatory_L1.out
#SBATCH --error=logs/gsea_excitatory_L1.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16

# Complete GSEA Pipeline for Excitatory Neurons (cell_type_L1 level)
#
# This script performs a comprehensive analysis:
# 1. Loads scRNA-seq data (annotation_final.h5ad)
# 2. Subsets by cell_type_L1 == 'Excitatory' (ALL excitatory subtypes combined)
# 3. Performs DEG analysis (Mutant vs Control) for each genotype separately
# 4. Runs focused GSEA with custom pathway keywords
# 5. Creates numbered heatmaps for visualization
#
# This is a true L1-level analysis treating all excitatory neurons as one group,
# not analyzing L2 subtypes separately.
#
# Pathway categories analyzed:
#   - cell_death: Neurodegeneration, apoptosis, autophagy
#   - synaptic: Synaptic function, glutamatergic transmission
#   - development: Neurogenesis, axonogenesis, dendritogenesis
#   - metabolism: Mitochondrial function, energy metabolism
#
# Output directories:
#   - GSEA results: results_from_raw/gsea_analysis_excitatory_L1/
#   - Heatmaps: results_from_raw/gsea_focused_heatmaps_excitatory_L1/
#
# Usage:
#   sbatch gsea_excitatory_L1_complete_pipeline.sh     # Submit to SLURM
#   ./gsea_excitatory_L1_complete_pipeline.sh          # Interactive execution
#

# Initialize conda
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/scrna-analysis

# Set working directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/ANALYSIS/Enrichment/GSEA/Conditions_from_degs"
cd "$SCRIPT_DIR"

echo "========================================================================"
echo "Complete GSEA Pipeline - Excitatory Neurons (L1 Level)"
echo "========================================================================"
echo "Cell type filter: cell_type_L1 == 'Excitatory'"
echo "Analysis: Combines ALL excitatory subtypes as one group"
echo "Genotypes: Emx1 and Nestin (analyzed separately)"
echo "Comparison: Mutant vs Control (within each genotype)"
echo "Pathway focus: Cell death, synaptic, development, metabolism"
echo "========================================================================"
echo ""
echo "Pipeline steps:"
echo "  1. Load scRNA-seq data"
echo "  2. Subset to Excitatory neurons (L1 level)"
echo "  3. Perform DEG analysis for each genotype"
echo "  4. Run GSEA with multiple gene set libraries"
echo "  5. Apply focused keyword filtering"
echo "  6. Create numbered heatmaps"
echo "========================================================================"

# Run the Python script
python3 "$SCRIPT_DIR/gsea_focused_heatmaps_numbered_excitatory.py"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "✓ SUCCESS: Complete pipeline finished!"
    echo "========================================================================"
    echo "Results generated:"
    echo "  ✓ DEG results for each genotype"
    echo "  ✓ GSEA analysis with 4 gene set libraries:"
    echo "      - KEGG_2019_Mouse"
    echo "      - GO_Biological_Process_2023"
    echo "      - Reactome_2022"
    echo "      - MSigDB_Hallmark_2020"
    echo "  ✓ Focused pathway results (cell_death, synaptic, development, metabolism)"
    echo "  ✓ Numbered heatmaps (with and without FDR annotations)"
    echo ""
    echo "Output directories:"
    echo "  - GSEA results: results_from_raw/gsea_analysis_excitatory_L1/"
    echo "  - Heatmaps: results_from_raw/gsea_focused_heatmaps_excitatory_L1/"
    echo ""
    echo "Files include:"
    echo "  - DEG_results.csv (differential expression)"
    echo "  - GSEA_Summary_All_Genesets.csv (significant terms)"
    echo "  - GSEA_Full_Report_All_Genesets.csv (all terms)"
    echo "  - GSEA_Focus_*.csv (focused pathway results)"
    echo "  - PNG/PDF heatmaps (600 DPI, publication-ready)"
    echo "  - CSV legend files (pathway number → name mapping)"
    echo "========================================================================"
else
    echo ""
    echo "========================================================================"
    echo "✗ ERROR: Pipeline failed with exit code $EXIT_CODE"
    echo "Check logs for details: logs/gsea_excitatory_L1.err"
    echo "========================================================================"
fi

exit $EXIT_CODE
