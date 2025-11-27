#!/bin/bash
#SBATCH --job-name=Striatum_EP_heatmaps
#SBATCH --output=logs/4_create_heatmaps.log
#SBATCH --error=logs/4_create_heatmaps.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --partition=workq

# Create Heatmaps for Striatum_EP CREs

echo "========================================================================"
echo "CREATE HEATMAPS FOR STRIATUM_EP CREs"
echo "========================================================================"
echo "Started: $(date)"

# Activate environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate sc-chromatin4

# Configuration
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/Striatum_EP_analysis"
cd $BASE_DIR

OUTPUT_DIR="./output/heatmaps_deeptools"
mkdir -p $OUTPUT_DIR

BED_FILE="./output/Striatum_EP_splicing_genes.bed"

# BigWig files (GABA)
BW_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/bigwig_tracks_L1/by_celltype"
BW_FILES="$BW_DIR/GABA_Nestin-Ctrl.bw $BW_DIR/GABA_Nestin-Mut.bw $BW_DIR/GABA_Emx1-Ctrl.bw $BW_DIR/GABA_Emx1-Mut.bw"
LABELS="Nestin-Ctrl Nestin-Mut Emx1-Ctrl Emx1-Mut"

# deepTools parameters
REGION_BODY_LENGTH=2000
BIN_SIZE=50
PROCESSORS=16

echo "Using BED file: $BED_FILE"
echo "Output directory: $OUTPUT_DIR"

# 1. computeMatrix
echo "Running computeMatrix..."
computeMatrix reference-point \
    --referencePoint center \
    --regionsFileName $BED_FILE \
    --scoreFileName $BW_FILES \
    --upstream 2000 \
    --downstream 2000 \
    --binSize $BIN_SIZE \
    --numberOfProcessors $PROCESSORS \
    --outFileName $OUTPUT_DIR/matrix_Striatum_EP.gz \
    --outFileNameMatrix $OUTPUT_DIR/matrix_Striatum_EP.tab \
    --missingDataAsZero

# 2. plotHeatmap
echo "Running plotHeatmap..."
plotHeatmap \
    --matrixFile $OUTPUT_DIR/matrix_Striatum_EP.gz \
    --outFileName $OUTPUT_DIR/heatmap_Striatum_EP.png \
    --samplesLabel $LABELS \
    --colorMap Reds \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 12 \
    --heatmapWidth 4

# 3. plotProfile
echo "Running plotProfile..."
plotProfile \
    --matrixFile $OUTPUT_DIR/matrix_Striatum_EP.gz \
    --outFileName $OUTPUT_DIR/metaprofile_Striatum_EP.png \
    --samplesLabel $LABELS \
    --perGroup \
    --colors '#808080' '#FF0000' '#808080' '#FF0000' \
    --plotHeight 8 \
    --plotWidth 10

echo "========================================================================"
echo "Finished: $(date)"
