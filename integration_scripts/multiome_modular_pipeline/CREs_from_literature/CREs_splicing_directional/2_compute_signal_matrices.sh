#!/bin/bash
#SBATCH --job-name=splicing_matrices
#SBATCH --output=logs/2_compute_signal_matrices_%j.log
#SBATCH --error=logs/2_compute_signal_matrices_%j.log
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=workq

# =============================================================================
# Step 2: Compute Signal Matrices for Splicing Gene CREs
# =============================================================================
#
# This script computes ATAC-seq signal matrices using deepTools for:
# - Enhancer-like CREs (positive PCC correlation)
# - Silencer-like CREs (negative PCC correlation)
#
# Comparisons (3 total, since Emx1-Ctrl is a FAILED sample):
# 1. Nestin-Ctrl vs Nestin-Mut (within-genotype mutation effect)
# 2. Nestin-Ctrl vs Emx1-Mut (cross-genotype mutation effect)
# 3. Nestin-Mut vs Emx1-Mut (mutant genotype comparison)
#
# =============================================================================

set -e  # Exit on error

# Absolute path to script directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_directional"

# Change to script directory
cd "${SCRIPT_DIR}"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate rna_seq_analysis_deep

# Paths
OUTPUT_DIR="${SCRIPT_DIR}/output"
MATRIX_DIR="${OUTPUT_DIR}/matrices"
HEATMAP_DIR="${OUTPUT_DIR}/heatmaps_deeptools"

# BigWig files (GABA cell type)
BIGWIG_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/bigwig_tracks_L1/by_celltype"

# NOTE: Emx1-Ctrl is EXCLUDED (failed sample with only 2.7 MB vs 12-28 MB for others)
NESTIN_CTRL="${BIGWIG_DIR}/GABA_Nestin-Ctrl.bw"
NESTIN_MUT="${BIGWIG_DIR}/GABA_Nestin-Mut.bw"
EMX1_MUT="${BIGWIG_DIR}/GABA_Emx1-Mut.bw"

# Parameters
WINDOW_SIZE=2000
BIN_SIZE=50
PROCESSORS=16

# Create directories
mkdir -p "${MATRIX_DIR}" "${HEATMAP_DIR}"

echo "=============================================================================="
echo "COMPUTING SIGNAL MATRICES FOR SPLICING GENE CREs"
echo "=============================================================================="
echo ""
echo "Date: $(date)"
echo "Output directory: ${OUTPUT_DIR}"
echo ""
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

# =============================================================================
# Process each CRE type (enhancers and silencers)
# =============================================================================

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

    # -------------------------------------------------------------------------
    # Compute matrix with all 3 conditions
    # -------------------------------------------------------------------------
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

    # -------------------------------------------------------------------------
    # Create heatmap
    # -------------------------------------------------------------------------
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

    # -------------------------------------------------------------------------
    # Create metaprofile
    # -------------------------------------------------------------------------
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

    # -------------------------------------------------------------------------
    # Create pairwise comparison matrices
    # -------------------------------------------------------------------------
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

# =============================================================================
# Create combined heatmap comparing enhancers vs silencers
# =============================================================================

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

# =============================================================================
# Summary
# =============================================================================

echo ""
echo "=============================================================================="
echo "MATRIX COMPUTATION COMPLETE"
echo "=============================================================================="
echo ""
echo "Output directories:"
echo "  Matrices: ${MATRIX_DIR}"
echo "  Heatmaps: ${HEATMAP_DIR}"
echo ""
echo "Key output files:"

for f in "${HEATMAP_DIR}"/*.png; do
    if [[ -f "$f" ]]; then
        echo "  $(basename "$f")"
    fi
done

echo ""
echo "Next: Run 3_visualize_directional_comparisons.py for detailed analysis"
echo ""
echo "Date: $(date)"
