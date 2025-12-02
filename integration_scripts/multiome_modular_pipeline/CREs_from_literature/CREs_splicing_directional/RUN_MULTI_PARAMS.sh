#!/bin/bash
# =============================================================================
# RUN PIPELINE WITH MULTIPLE PARAMETER COMBINATIONS
# =============================================================================
#
# This script runs the splicing directional pipeline with multiple
# MIN_SIGNAL and MIN_FC parameter combinations.
#
# Steps 1-2 run once (extract CREs and compute matrices - same for all params)
# Step 3 runs multiple times with different thresholds for individual plots
#
# USAGE:
#   ./RUN_MULTI_PARAMS.sh           # Submit all jobs
#   ./RUN_MULTI_PARAMS.sh --dry-run # Show what would be submitted
#
# =============================================================================

set -e

SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_directional"
cd "${SCRIPT_DIR}"

# Create logs directory
mkdir -p logs

# Check for dry-run mode
DRY_RUN=0
if [[ "$1" == "--dry-run" ]]; then
    DRY_RUN=1
    echo "DRY-RUN MODE - showing commands without executing"
    echo ""
fi

# Define parameter combinations
# Format: "MIN_SIGNAL MIN_FC"
PARAM_COMBOS=(
    "3.0 3.0"
    "2.0 3.0"
    "2.0 2.0"
)

echo "=============================================================================="
echo "MULTI-PARAMETER PIPELINE SUBMISSION"
echo "=============================================================================="
echo ""
echo "Parameter combinations to run:"
for combo in "${PARAM_COMBOS[@]}"; do
    read -r SIG FC <<< "$combo"
    echo "  - MIN_SIGNAL=${SIG}, MIN_FC=${FC}"
done
echo ""

# =============================================================================
# STEP 1: Submit Steps 1-2 job (same for all parameter combinations)
# =============================================================================

echo "=============================================================================="
echo "Submitting Steps 1-2 (Extract CREs + Compute Matrices)"
echo "=============================================================================="

STEPS12_SCRIPT="${SCRIPT_DIR}/run_steps_1_2.sh"

# Create Steps 1-2 script if it doesn't exist
cat > "${STEPS12_SCRIPT}" << 'STEPS12_EOF'
#!/bin/bash
#SBATCH --job-name=splicing_steps12
#SBATCH --output=logs/steps_1_2_%j.log
#SBATCH --error=logs/steps_1_2_%j.log
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=workq

set -e

SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_directional"
cd "${SCRIPT_DIR}"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate rna_seq_analysis_deep

echo "=============================================================================="
echo "STEPS 1-2: Extract CREs and Compute Signal Matrices"
echo "=============================================================================="
echo "Date: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
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
echo "=============================================================================="
echo "STEPS 1-2 COMPLETE"
echo "=============================================================================="
echo "End time: $(date)"
STEPS12_EOF

chmod +x "${STEPS12_SCRIPT}"

if [[ $DRY_RUN -eq 1 ]]; then
    echo "  [DRY-RUN] Would submit: sbatch ${STEPS12_SCRIPT}"
    STEPS12_JOB="12345"
else
    STEPS12_JOB=$(sbatch --parsable "${STEPS12_SCRIPT}")
    echo "  Submitted job ${STEPS12_JOB}: Steps 1-2 (Extract CREs + Compute Matrices)"
fi

echo ""

# =============================================================================
# STEP 2: Submit Step 3 jobs for each parameter combination
# =============================================================================

echo "=============================================================================="
echo "Submitting Step 3 Jobs (Visualization with different thresholds)"
echo "=============================================================================="

STEP3_JOBS=()

for combo in "${PARAM_COMBOS[@]}"; do
    read -r MIN_SIG MIN_FC <<< "$combo"

    JOB_NAME="splicing_step3_sig${MIN_SIG}_fc${MIN_FC}"
    LOG_FILE="logs/step_3_sig${MIN_SIG}_fc${MIN_FC}_%j.log"

    echo ""
    echo "  MIN_SIGNAL=${MIN_SIG}, MIN_FC=${MIN_FC}"

    if [[ $DRY_RUN -eq 1 ]]; then
        echo "    [DRY-RUN] Would submit: sbatch --dependency=afterok:${STEPS12_JOB} --export=MIN_SIGNAL=${MIN_SIG},MIN_FC=${MIN_FC} ..."
    else
        JOB_ID=$(sbatch --parsable \
            --dependency=afterok:${STEPS12_JOB} \
            --job-name="${JOB_NAME}" \
            --output="${LOG_FILE}" \
            --error="${LOG_FILE}" \
            --export=MIN_SIGNAL=${MIN_SIG},MIN_FC=${MIN_FC} \
            "${SCRIPT_DIR}/3_visualize_directional_comparisons.sh")

        echo "    Submitted job ${JOB_ID} (depends on ${STEPS12_JOB})"
        STEP3_JOBS+=("${JOB_ID}")
    fi
done

echo ""
echo "=============================================================================="
echo "SUBMISSION COMPLETE"
echo "=============================================================================="
echo ""

if [[ $DRY_RUN -eq 1 ]]; then
    echo "DRY-RUN complete. No jobs were submitted."
else
    echo "Submitted jobs:"
    echo "  Steps 1-2: ${STEPS12_JOB}"
    echo "  Step 3 jobs: ${STEP3_JOBS[*]}"
    echo ""
    echo "Monitor progress with:"
    echo "  squeue -u \$USER"
    echo ""
    echo "View logs:"
    echo "  tail -f logs/steps_1_2_${STEPS12_JOB}.log"
    for combo in "${PARAM_COMBOS[@]}"; do
        read -r MIN_SIG MIN_FC <<< "$combo"
        echo "  tail -f logs/step_3_sig${MIN_SIG}_fc${MIN_FC}_*.log"
    done
fi
echo ""
