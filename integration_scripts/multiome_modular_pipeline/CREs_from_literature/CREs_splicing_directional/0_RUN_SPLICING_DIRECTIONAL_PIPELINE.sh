#!/bin/bash
#SBATCH --job-name=splicing_directional_pipeline
#SBATCH --output=logs/0_pipeline_%j.out
#SBATCH --error=logs/0_pipeline_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=workq

# =============================================================================
# SPLICING GENE CRE ANALYSIS PIPELINE - WITH DIRECTIONALITY
# =============================================================================
#
# This pipeline analyzes chromatin accessibility at splicing gene regulatory
# elements, specifically designed to identify NESTIN-SPECIFIC LOSS of
# accessibility in GABA neurons.
#
# KEY INNOVATION: Separates CREs by correlation direction
#   - Enhancer-like CREs (PCC > 0): accessibility ↑ = expression ↑
#   - Silencer-like CREs (PCC < 0): accessibility ↑ = expression ↓
#
# RESEARCH QUESTION:
#   Identify loss of chromatin accessibility at splicing gene regulatory
#   elements specifically in Nestin-MUT GABA neurons compared to:
#   1. Nestin-Ctrl (within-genotype mutation effect)
#   2. Emx1-Mut (genotype-specific response)
#
# USAGE:
#   sbatch 0_RUN_SPLICING_DIRECTIONAL_PIPELINE.sh
#
# ENVIRONMENT VARIABLES (optional):
#   SKIP_INDIVIDUAL=1    Skip individual CRE plots (faster)
#   MIN_SIGNAL=2.0       Minimum signal threshold for individual plots
#   MIN_FC=2.0           Minimum fold change for individual plots
#   MAX_INDIVIDUAL=50    Maximum number of individual plots
#
# =============================================================================

set -e

# Script directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_directional"
cd "${SCRIPT_DIR}"

# Create logs directory
mkdir -p logs

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate rna_seq_analysis_deep

# Parse environment variables with defaults
SKIP_INDIVIDUAL="${SKIP_INDIVIDUAL:-0}"
MIN_SIGNAL="${MIN_SIGNAL:-2.0}"
MIN_FC="${MIN_FC:-2.0}"
MAX_INDIVIDUAL="${MAX_INDIVIDUAL:-50}"
INDIVIDUAL_DPI="${INDIVIDUAL_DPI:-150}"

echo "=============================================================================="
echo "SPLICING GENE CRE ANALYSIS PIPELINE - WITH DIRECTIONALITY"
echo "=============================================================================="
echo ""
echo "Date: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Working directory: ${SCRIPT_DIR}"
echo ""
echo "Settings:"
echo "  SKIP_INDIVIDUAL: ${SKIP_INDIVIDUAL}"
echo "  MIN_SIGNAL: ${MIN_SIGNAL}"
echo "  MIN_FC: ${MIN_FC}"
echo "  MAX_INDIVIDUAL: ${MAX_INDIVIDUAL}"
echo "  INDIVIDUAL_DPI: ${INDIVIDUAL_DPI}"
echo ""

# =============================================================================
# STEP 1: Extract CREs with PCC directionality
# =============================================================================

echo "=============================================================================="
echo "STEP 1: Extract Splicing Gene CREs with Directionality"
echo "=============================================================================="
echo "Start time: $(date)"
echo ""

python "${SCRIPT_DIR}/1_extract_splicing_CREs_directional.py"

echo ""
echo "Step 1 complete: $(date)"

# =============================================================================
# STEP 2: Compute signal matrices
# =============================================================================

echo ""
echo "=============================================================================="
echo "STEP 2: Compute Signal Matrices"
echo "=============================================================================="
echo "Start time: $(date)"
echo ""

# Paths
OUTPUT_DIR="${SCRIPT_DIR}/output"
MATRIX_DIR="${OUTPUT_DIR}/matrices"
HEATMAP_DIR="${OUTPUT_DIR}/heatmaps_deeptools"

# BigWig files (GABA cell type) - NOTE: Emx1-Ctrl EXCLUDED (failed sample)
BIGWIG_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/bigwig_tracks_L1/by_celltype"
NESTIN_CTRL="${BIGWIG_DIR}/GABA_Nestin-Ctrl.bw"
NESTIN_MUT="${BIGWIG_DIR}/GABA_Nestin-Mut.bw"
EMX1_MUT="${BIGWIG_DIR}/GABA_Emx1-Mut.bw"

# Parameters
WINDOW_SIZE=2000
BIN_SIZE=50
PROCESSORS=16

# Create directories
mkdir -p "${MATRIX_DIR}" "${HEATMAP_DIR}"

echo "BigWig files:"
echo "  Nestin-Ctrl: ${NESTIN_CTRL}"
echo "  Nestin-Mut:  ${NESTIN_MUT}"
echo "  Emx1-Mut:    ${EMX1_MUT}"
echo ""
echo "NOTE: Emx1-Ctrl is EXCLUDED (failed sample)"
echo ""

# Verify files exist
for bw in "${NESTIN_CTRL}" "${NESTIN_MUT}" "${EMX1_MUT}"; do
    if [[ ! -f "${bw}" ]]; then
        echo "ERROR: BigWig file not found: ${bw}"
        exit 1
    fi
done

# Process each CRE type (enhancers and silencers)
for CRE_TYPE in "enhancer" "silencer"; do
    BED_FILE="${OUTPUT_DIR}/${CRE_TYPE}_CREs_GABA.bed"

    if [[ ! -f "${BED_FILE}" ]]; then
        echo "WARNING: BED file not found: ${BED_FILE}"
        echo "         Skipping ${CRE_TYPE} analysis"
        continue
    fi

    CRE_COUNT=$(wc -l < "${BED_FILE}")
    echo ""
    echo "=============================================================================="
    echo "Processing ${CRE_TYPE^^} CREs (${CRE_COUNT} regions)"
    echo "=============================================================================="

    # Output files
    MATRIX_FILE="${MATRIX_DIR}/matrix_${CRE_TYPE}_GABA.gz"
    MATRIX_TAB="${MATRIX_DIR}/matrix_${CRE_TYPE}_GABA.tab"
    HEATMAP_FILE="${HEATMAP_DIR}/heatmap_${CRE_TYPE}_all_conditions.png"
    PROFILE_FILE="${HEATMAP_DIR}/metaprofile_${CRE_TYPE}_all_conditions.png"

    # Compute matrix with all 3 conditions
    echo ""
    echo "Computing signal matrix..."
    echo "  Window: +/- ${WINDOW_SIZE} bp"
    echo "  Bin size: ${BIN_SIZE} bp"

    computeMatrix reference-point \
        --referencePoint center \
        --beforeRegionStartLength ${WINDOW_SIZE} \
        --afterRegionStartLength ${WINDOW_SIZE} \
        --binSize ${BIN_SIZE} \
        --regionsFileName "${BED_FILE}" \
        --scoreFileName "${NESTIN_CTRL}" "${NESTIN_MUT}" "${EMX1_MUT}" \
        --samplesLabel "Nestin-Ctrl" "Nestin-Mut" "Emx1-Mut" \
        --outFileName "${MATRIX_FILE}" \
        --outFileNameMatrix "${MATRIX_TAB}" \
        --numberOfProcessors ${PROCESSORS} \
        --missingDataAsZero \
        --skipZeros

    echo "  Saved: ${MATRIX_FILE}"

    # Create heatmap
    echo ""
    echo "Creating heatmap..."

    plotHeatmap \
        --matrixFile "${MATRIX_FILE}" \
        --outFileName "${HEATMAP_FILE}" \
        --colorMap "RdYlBu_r" \
        --whatToShow "heatmap and colorbar" \
        --sortUsing mean \
        --sortUsingSamples 1 \
        --zMin 0 \
        --heatmapHeight 15 \
        --heatmapWidth 5 \
        --dpi 300 \
        --plotTitle "${CRE_TYPE^} CREs at Splicing Genes (GABA)"

    echo "  Saved: ${HEATMAP_FILE}"

    # Create metaprofile
    echo ""
    echo "Creating metaprofile..."

    plotProfile \
        --matrixFile "${MATRIX_FILE}" \
        --outFileName "${PROFILE_FILE}" \
        --colors "#3182bd" "#de2d26" "#ff7f00" \
        --plotType lines \
        --perGroup \
        --dpi 300 \
        --plotTitle "${CRE_TYPE^} CREs at Splicing Genes (GABA)" \
        --yAxisLabel "ATAC Signal" \
        --regionsLabel "${CRE_TYPE^} CREs"

    echo "  Saved: ${PROFILE_FILE}"

    # Create pairwise comparison matrices
    echo ""
    echo "Creating pairwise comparison matrices..."

    # Comparison 1: Nestin-Ctrl vs Nestin-Mut
    MATRIX_NESTIN="${MATRIX_DIR}/matrix_${CRE_TYPE}_nestin_comparison.gz"

    computeMatrix reference-point \
        --referencePoint center \
        --beforeRegionStartLength ${WINDOW_SIZE} \
        --afterRegionStartLength ${WINDOW_SIZE} \
        --binSize ${BIN_SIZE} \
        --regionsFileName "${BED_FILE}" \
        --scoreFileName "${NESTIN_CTRL}" "${NESTIN_MUT}" \
        --samplesLabel "Nestin-Ctrl" "Nestin-Mut" \
        --outFileName "${MATRIX_NESTIN}" \
        --numberOfProcessors ${PROCESSORS} \
        --missingDataAsZero \
        --skipZeros

    plotProfile \
        --matrixFile "${MATRIX_NESTIN}" \
        --outFileName "${HEATMAP_DIR}/metaprofile_${CRE_TYPE}_nestin_ctrl_vs_mut.png" \
        --colors "#3182bd" "#de2d26" \
        --plotType lines \
        --perGroup \
        --dpi 300 \
        --plotTitle "Nestin Ctrl vs Mut: ${CRE_TYPE^} CREs" \
        --yAxisLabel "ATAC Signal"

    echo "  Saved: metaprofile_${CRE_TYPE}_nestin_ctrl_vs_mut.png"

    # Comparison 2: Nestin-Ctrl vs Emx1-Mut
    MATRIX_CROSS="${MATRIX_DIR}/matrix_${CRE_TYPE}_cross_comparison.gz"

    computeMatrix reference-point \
        --referencePoint center \
        --beforeRegionStartLength ${WINDOW_SIZE} \
        --afterRegionStartLength ${WINDOW_SIZE} \
        --binSize ${BIN_SIZE} \
        --regionsFileName "${BED_FILE}" \
        --scoreFileName "${NESTIN_CTRL}" "${EMX1_MUT}" \
        --samplesLabel "Nestin-Ctrl" "Emx1-Mut" \
        --outFileName "${MATRIX_CROSS}" \
        --numberOfProcessors ${PROCESSORS} \
        --missingDataAsZero \
        --skipZeros

    plotProfile \
        --matrixFile "${MATRIX_CROSS}" \
        --outFileName "${HEATMAP_DIR}/metaprofile_${CRE_TYPE}_nestin_ctrl_vs_emx1_mut.png" \
        --colors "#3182bd" "#ff7f00" \
        --plotType lines \
        --perGroup \
        --dpi 300 \
        --plotTitle "Nestin-Ctrl vs Emx1-Mut: ${CRE_TYPE^} CREs" \
        --yAxisLabel "ATAC Signal"

    echo "  Saved: metaprofile_${CRE_TYPE}_nestin_ctrl_vs_emx1_mut.png"

    # Comparison 3: Nestin-Mut vs Emx1-Mut
    MATRIX_MUT="${MATRIX_DIR}/matrix_${CRE_TYPE}_mutant_comparison.gz"

    computeMatrix reference-point \
        --referencePoint center \
        --beforeRegionStartLength ${WINDOW_SIZE} \
        --afterRegionStartLength ${WINDOW_SIZE} \
        --binSize ${BIN_SIZE} \
        --regionsFileName "${BED_FILE}" \
        --scoreFileName "${NESTIN_MUT}" "${EMX1_MUT}" \
        --samplesLabel "Nestin-Mut" "Emx1-Mut" \
        --outFileName "${MATRIX_MUT}" \
        --numberOfProcessors ${PROCESSORS} \
        --missingDataAsZero \
        --skipZeros

    plotProfile \
        --matrixFile "${MATRIX_MUT}" \
        --outFileName "${HEATMAP_DIR}/metaprofile_${CRE_TYPE}_nestin_mut_vs_emx1_mut.png" \
        --colors "#de2d26" "#ff7f00" \
        --plotType lines \
        --perGroup \
        --dpi 300 \
        --plotTitle "Nestin-Mut vs Emx1-Mut: ${CRE_TYPE^} CREs" \
        --yAxisLabel "ATAC Signal"

    echo "  Saved: metaprofile_${CRE_TYPE}_nestin_mut_vs_emx1_mut.png"

done

# Create combined heatmap comparing enhancers vs silencers
echo ""
echo "=============================================================================="
echo "Creating combined enhancer vs silencer comparison"
echo "=============================================================================="

ENHANCER_BED="${OUTPUT_DIR}/enhancer_CREs_GABA.bed"
SILENCER_BED="${OUTPUT_DIR}/silencer_CREs_GABA.bed"

if [[ -f "${ENHANCER_BED}" && -f "${SILENCER_BED}" ]]; then
    COMBINED_MATRIX="${MATRIX_DIR}/matrix_combined_enhancer_silencer.gz"

    computeMatrix reference-point \
        --referencePoint center \
        --beforeRegionStartLength ${WINDOW_SIZE} \
        --afterRegionStartLength ${WINDOW_SIZE} \
        --binSize ${BIN_SIZE} \
        --regionsFileName "${ENHANCER_BED}" "${SILENCER_BED}" \
        --scoreFileName "${NESTIN_CTRL}" "${NESTIN_MUT}" "${EMX1_MUT}" \
        --samplesLabel "Nestin-Ctrl" "Nestin-Mut" "Emx1-Mut" \
        --outFileName "${COMBINED_MATRIX}" \
        --numberOfProcessors ${PROCESSORS} \
        --missingDataAsZero \
        --skipZeros

    plotHeatmap \
        --matrixFile "${COMBINED_MATRIX}" \
        --outFileName "${HEATMAP_DIR}/heatmap_combined_enhancer_vs_silencer.png" \
        --colorMap "RdYlBu_r" \
        --whatToShow "heatmap and colorbar" \
        --sortUsing mean \
        --sortUsingSamples 1 \
        --zMin 0 \
        --heatmapHeight 15 \
        --heatmapWidth 5 \
        --dpi 300 \
        --regionsLabel "Enhancer CREs" "Silencer CREs" \
        --plotTitle "Enhancer vs Silencer CREs at Splicing Genes"

    plotProfile \
        --matrixFile "${COMBINED_MATRIX}" \
        --outFileName "${HEATMAP_DIR}/metaprofile_combined_enhancer_vs_silencer.png" \
        --colors "#3182bd" "#de2d26" "#ff7f00" \
        --plotType lines \
        --perGroup \
        --dpi 300 \
        --plotTitle "Enhancer vs Silencer CREs at Splicing Genes" \
        --yAxisLabel "ATAC Signal" \
        --regionsLabel "Enhancer CREs" "Silencer CREs"

    echo "  Saved: heatmap_combined_enhancer_vs_silencer.png"
    echo "  Saved: metaprofile_combined_enhancer_vs_silencer.png"
fi

echo ""
echo "Step 2 complete: $(date)"

# =============================================================================
# STEP 3: Visualize and identify Nestin-specific loss
# =============================================================================

echo ""
echo "=============================================================================="
echo "STEP 3: Visualize Comparisons and Identify Nestin-Specific Loss"
echo "=============================================================================="
echo "Start time: $(date)"
echo ""

# Build arguments for visualization script
VIZ_ARGS=""
if [[ "${SKIP_INDIVIDUAL}" == "1" ]]; then
    VIZ_ARGS="${VIZ_ARGS} --skip-individual"
fi
VIZ_ARGS="${VIZ_ARGS} --max-individual ${MAX_INDIVIDUAL}"
VIZ_ARGS="${VIZ_ARGS} --min-signal ${MIN_SIGNAL}"
VIZ_ARGS="${VIZ_ARGS} --min-fc ${MIN_FC}"
VIZ_ARGS="${VIZ_ARGS} --individual-dpi ${INDIVIDUAL_DPI}"

echo "Running visualization with arguments: ${VIZ_ARGS}"
echo ""

python "${SCRIPT_DIR}/3_visualize_directional_comparisons.py" ${VIZ_ARGS}

echo ""
echo "Step 3 complete: $(date)"

# =============================================================================
# PIPELINE COMPLETE
# =============================================================================

echo ""
echo "=============================================================================="
echo "PIPELINE COMPLETE"
echo "=============================================================================="
echo ""
echo "End time: $(date)"
echo ""
echo "Output directories:"
echo "  ${OUTPUT_DIR}/matrices/"
echo "  ${OUTPUT_DIR}/heatmaps_deeptools/"
echo "  ${OUTPUT_DIR}/visualizations/"
echo "  ${OUTPUT_DIR}/nestin_specific_loss/"
echo ""
echo "Key output files:"
echo "  output/enhancer_CREs_GABA.bed"
echo "  output/silencer_CREs_GABA.bed"
echo "  output/heatmaps_deeptools/*.png"
echo "  output/visualizations/*.png"
echo "  output/nestin_specific_loss/*.tsv"
echo ""
