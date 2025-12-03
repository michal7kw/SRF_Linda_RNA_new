#!/bin/bash
#SBATCH --job-name=diff_matrices
#SBATCH --output=logs/2_compute_signal_matrices_%j.log
#SBATCH --error=logs/2_compute_signal_matrices_%j.log
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=workq

# =============================================================================
# Step 2: Compute Signal Matrices for Differential CRE Analysis (All Cell Types)
# =============================================================================

set -e  # Exit on error

# Absolute path to script directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_differential_all_celltypes"

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
# Note: Even though we look at "all cell types" CREs, we are still analyzing the signal
# in our GABA ATAC-seq data to see if these CREs are active/differential in our samples.
BIGWIG_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/bigwig_tracks_L1/by_celltype"

# NOTE: Emx1-Ctrl is EXCLUDED (failed sample)
NESTIN_CTRL="${BIGWIG_DIR}/GABA_Nestin-Ctrl.bw"
NESTIN_MUT="${BIGWIG_DIR}/GABA_Nestin-Mut.bw"
EMX1_MUT="${BIGWIG_DIR}/GABA_Emx1-Mut.bw"

# Parameters
WINDOW_SIZE=2000
BIN_SIZE=50
PROCESSORS=16

# Create directories
mkdir -p "${MATRIX_DIR}" "${HEATMAP_DIR}" "logs"

echo "=============================================================================="
echo "COMPUTING SIGNAL MATRICES FOR DIFFERENTIAL CREs (ALL CELL TYPES)"
echo "=============================================================================="
echo ""
echo "Date: $(date)"
echo "Output directory: ${OUTPUT_DIR}"
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
    # Note: Using the "all_celltypes" BED files
    BED_FILE="${OUTPUT_DIR}/${CRE_TYPE}_CREs_all_celltypes.bed"

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
    MATRIX_FILE="${MATRIX_DIR}/matrix_${CRE_TYPE}_all.gz"
    MATRIX_TAB="${MATRIX_DIR}/matrix_${CRE_TYPE}_all.tab"
    HEATMAP_FILE="${HEATMAP_DIR}/heatmap_${CRE_TYPE}_all_conditions.png"
    PROFILE_FILE="${HEATMAP_DIR}/metaprofile_${CRE_TYPE}_all_conditions.png"

    # -------------------------------------------------------------------------
    # Compute matrix with all 3 conditions
    # -------------------------------------------------------------------------
    echo ""
    echo "Computing signal matrix..."
    
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
    # Create heatmap and profile
    # -------------------------------------------------------------------------
    echo "Creating visualizations..."

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
        --plotTitle "${CRE_TYPE^} CREs (All Cell Types)"

    plotProfile \
        --matrixFile "${MATRIX_FILE}" \
        --outFileName "${PROFILE_FILE}" \
        --colors "#3182bd" "#de2d26" "#ff7f00" \
        --plotType lines \
        --perGroup \
        --dpi 300 \
        --plotTitle "${CRE_TYPE^} CREs (All Cell Types)" \
        --yAxisLabel "ATAC Signal" \
        --regionsLabel "${CRE_TYPE^} CREs"

    echo "  Saved: ${HEATMAP_FILE}"
    echo "  Saved: ${PROFILE_FILE}"

done

echo ""
echo "=============================================================================="
echo "MATRIX COMPUTATION COMPLETE"
echo "=============================================================================="
echo "Next: Run 3_analyze_differential_cres.py"
