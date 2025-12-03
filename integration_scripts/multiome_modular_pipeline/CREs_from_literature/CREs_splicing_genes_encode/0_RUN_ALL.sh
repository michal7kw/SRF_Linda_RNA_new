#!/bin/bash
#SBATCH --job-name=0_CREs_splicing_genes_encode
#SBATCH --output=logs/0_CREs_splicing_genes_encode_%j.log
#SBATCH --error=logs/0_CREs_splicing_genes_encode_%j.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=workq

################################################################################
# Master Script: Splicing Genes ENCODE cCRE Analysis Pipeline
#
# This script runs the complete pipeline sequentially:
# 1. 1_extract_encode_cCREs (Extracts ENCODE cCREs linked to splicing genes)
# 2. 2_convert_to_bed (Converts TSV to BED format)
# 3. 3_create_profiles (Creates profiles with deepTools)
# 4. 4_create_heatmaps (Creates heatmaps with deepTools)
# 5. 5_visualize_bigwig_signal (Visualizes BigWig signal directly)
# 6. 6_create_custom_comparisons (Creates custom comparison plots)
#
# DATA SOURCE: ENCODE cCREs (mm10-cCREs.bed) linked by genomic proximity
# NOTE: Emx1-Ctrl is excluded (failed sample) - Nestin-Ctrl used as baseline
#
# Usage: sbatch 0_RUN_SPLICING_ENCODE_CCRES_ANALYSIS.sh
################################################################################

echo "========================================================================"
echo "SPLICING GENES ENCODE cCRE ANALYSIS PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""
echo "DATA SOURCE: ENCODE cCREs (mm10-cCREs.bed)"
echo "LINKAGE METHOD: Genomic proximity (+/- 500kb)"
echo "NOTE: Emx1-Ctrl excluded (failed sample)"
echo ""

# Change to pipeline directory
cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_genes_encode"

# Create logs directory if it doesn't exist
mkdir -p logs

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis_deep

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# Make all scripts executable
chmod +x *.sh *.py

# Track overall success
PIPELINE_SUCCESS=true

################################################################################
# STEP 1: Extract ENCODE cCREs
################################################################################
echo "========================================================================"
echo "STEP 1/6: Extracting ENCODE cCREs linked to splicing genes"
echo "========================================================================"
echo "Started: $(date)"
echo ""

python 1_extract_encode_cCREs.py

if [ $? -ne 0 ]; then
    echo "ERROR: Step 1 failed!"
    PIPELINE_SUCCESS=false
else
    echo ""
    echo "Step 1 completed successfully"
fi
echo ""

################################################################################
# STEP 2: Convert to BED format
################################################################################
if [ "$PIPELINE_SUCCESS" = true ]; then
    echo "========================================================================"
    echo "STEP 2/6: Converting TSV files to BED format"
    echo "========================================================================"
    echo "Started: $(date)"
    echo ""

    python 2_convert_to_bed.py

    if [ $? -ne 0 ]; then
        echo "ERROR: Step 2 failed!"
        PIPELINE_SUCCESS=false
    else
        echo ""
        echo "Step 2 completed successfully"
    fi
    echo ""
fi

################################################################################
# STEP 3: Create profiles with deepTools
################################################################################
if [ "$PIPELINE_SUCCESS" = true ]; then
    echo "========================================================================"
    echo "STEP 3/6: Creating ATAC signal profiles"
    echo "========================================================================"
    echo "Started: $(date)"
    echo ""

    # Source the shell script content directly (excluding SBATCH headers)
    # Run computeMatrix and Python visualization

    BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
    OUTPUT_DIR="./output/heatmaps_deeptools"
    BED_ALL="./output/CREs_splicing_genes_encode_all.bed"
    WINDOW_SIZE=2000
    BIN_SIZE=50
    N_PROCESSORS=16

    mkdir -p $OUTPUT_DIR

    # Nestin matrix
    echo "Computing Nestin signal matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b $WINDOW_SIZE -a $WINDOW_SIZE \
        -R $BED_ALL \
        -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Nestin-Mut.bw" \
        --samplesLabel "Nestin-Ctrl" "Nestin-Mut" \
        --binSize $BIN_SIZE \
        --sortRegions keep \
        --missingDataAsZero \
        -o $OUTPUT_DIR/matrix_GABA_nestin.gz \
        -p $N_PROCESSORS \
        --outFileNameMatrix $OUTPUT_DIR/matrix_GABA_nestin.tab

    # Emx1 matrix (using Nestin-Ctrl as baseline)
    echo "Computing Emx1 signal matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b $WINDOW_SIZE -a $WINDOW_SIZE \
        -R $BED_ALL \
        -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
        --samplesLabel "Nestin-Ctrl" "Emx1-Mut" \
        --binSize $BIN_SIZE \
        --sortRegions keep \
        --missingDataAsZero \
        -o $OUTPUT_DIR/matrix_GABA_emx1.gz \
        -p $N_PROCESSORS \
        --outFileNameMatrix $OUTPUT_DIR/matrix_GABA_emx1.tab

    # Run Python visualization
    echo "Creating profile visualizations..."
    python 3_create_profiles.py --skip-individual

    if [ $? -ne 0 ]; then
        echo "ERROR: Step 3 failed!"
        PIPELINE_SUCCESS=false
    else
        echo ""
        echo "Step 3 completed successfully"
    fi
    echo ""
fi

################################################################################
# STEP 4: Create heatmaps with deepTools
################################################################################
if [ "$PIPELINE_SUCCESS" = true ]; then
    echo "========================================================================"
    echo "STEP 4/6: Creating heatmaps and metaprofiles"
    echo "========================================================================"
    echo "Started: $(date)"
    echo ""

    OUTPUT_DIR="./output/heatmaps_deeptools"
    BED_ALL="./output/CREs_splicing_genes_encode_all.bed"
    N_CRES=$(wc -l < $BED_ALL)

    # All samples matrix
    echo "Computing all samples matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 2000 -a 2000 \
        -R $BED_ALL \
        -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Nestin-Mut.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
        --samplesLabel "Nestin-Ctrl" "Nestin-Mut" "Emx1-Mut" \
        --binSize 50 \
        --sortRegions keep \
        --missingDataAsZero \
        -o $OUTPUT_DIR/matrix_all.gz \
        -p 16 \
        --outFileNameMatrix $OUTPUT_DIR/matrix_all.tab

    # Create heatmaps
    echo "Creating heatmaps..."
    plotHeatmap \
        -m $OUTPUT_DIR/matrix_all.gz \
        -o $OUTPUT_DIR/heatmap_all.png \
        --colorMap Reds \
        --dpi 300 \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 \
        --refPointLabel "CRE Center" \
        --heatmapHeight 10 \
        --heatmapWidth 6 \
        --sortRegions descend \
        --sortUsing mean \
        --plotTitle "ATAC Signal at Splicing Gene ENCODE cCREs (n=$N_CRES)"

    # Nestin heatmap
    plotHeatmap \
        -m $OUTPUT_DIR/matrix_GABA_nestin.gz \
        -o $OUTPUT_DIR/heatmap_GABA_nestin.png \
        --colorMap Reds \
        --dpi 300 \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 \
        --refPointLabel "CRE Center" \
        --heatmapHeight 10 \
        --heatmapWidth 3 \
        --plotTitle "Nestin: Ctrl vs Mut at Splicing ENCODE cCREs (n=$N_CRES)"

    # Emx1 heatmap
    plotHeatmap \
        -m $OUTPUT_DIR/matrix_GABA_emx1.gz \
        -o $OUTPUT_DIR/heatmap_GABA_emx1.png \
        --colorMap Reds \
        --dpi 300 \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 \
        --refPointLabel "CRE Center" \
        --heatmapHeight 10 \
        --heatmapWidth 3 \
        --plotTitle "Nestin-Ctrl vs Emx1-Mut at Splicing ENCODE cCREs (n=$N_CRES)"

    # Create metaprofiles
    echo "Creating metaprofiles..."
    plotProfile \
        -m $OUTPUT_DIR/matrix_all.gz \
        -o $OUTPUT_DIR/metaprofile_all.png \
        --plotTitle "ATAC Signal at Splicing Gene ENCODE cCREs (n=$N_CRES)" \
        --refPointLabel "CRE Center" \
        --colors '#2E86AB' '#A23B72' '#F18F01' \
        --plotHeight 7 \
        --plotWidth 10

    plotProfile \
        -m $OUTPUT_DIR/matrix_GABA_nestin.gz \
        -o $OUTPUT_DIR/metaprofile_GABA_nestin.png \
        --plotTitle "Nestin: Ctrl vs Mut at Splicing ENCODE cCREs (n=$N_CRES)" \
        --refPointLabel "CRE Center" \
        --colors '#2E86AB' '#A23B72' \
        --plotHeight 6 \
        --plotWidth 8

    plotProfile \
        -m $OUTPUT_DIR/matrix_GABA_emx1.gz \
        -o $OUTPUT_DIR/metaprofile_GABA_emx1.png \
        --plotTitle "Nestin-Ctrl vs Emx1-Mut at Splicing ENCODE cCREs (n=$N_CRES)" \
        --refPointLabel "CRE Center" \
        --colors '#2E86AB' '#F18F01' \
        --plotHeight 6 \
        --plotWidth 8

    echo ""
    echo "Step 4 completed successfully"
    echo ""
fi

################################################################################
# STEP 5: Visualize BigWig signal directly
################################################################################
if [ "$PIPELINE_SUCCESS" = true ]; then
    echo "========================================================================"
    echo "STEP 5/6: Direct BigWig visualization"
    echo "========================================================================"
    echo "Started: $(date)"
    echo ""

    python 5_visualize_bigwig_signal.py --skip-individual

    if [ $? -ne 0 ]; then
        echo "ERROR: Step 5 failed!"
        PIPELINE_SUCCESS=false
    else
        echo ""
        echo "Step 5 completed successfully"
    fi
    echo ""
fi

################################################################################
# STEP 6: Create custom comparisons
################################################################################
if [ "$PIPELINE_SUCCESS" = true ]; then
    echo "========================================================================"
    echo "STEP 6/6: Creating custom cross-genotype comparisons"
    echo "========================================================================"
    echo "Started: $(date)"
    echo ""

    OUTPUT_DIR="./output/custom_comparisons"
    BED_ALL="./output/CREs_splicing_genes_encode_all.bed"
    mkdir -p $OUTPUT_DIR

    # Comparison 1: Nestin-Ctrl vs Nestin-Mut
    echo "Computing matrix: Nestin-Ctrl vs Nestin-Mut..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 2000 -a 2000 \
        -R $BED_ALL \
        -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Nestin-Mut.bw" \
        --samplesLabel "Nestin-Ctrl" "Nestin-Mut" \
        --binSize 50 \
        --sortRegions keep \
        --missingDataAsZero \
        -o $OUTPUT_DIR/matrix_nestin_ctrl_vs_mut.gz \
        -p 16 \
        --outFileNameMatrix $OUTPUT_DIR/matrix_nestin_ctrl_vs_mut.tab

    # Comparison 2: Nestin-Ctrl vs Emx1-Mut
    echo "Computing matrix: Nestin-Ctrl vs Emx1-Mut..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 2000 -a 2000 \
        -R $BED_ALL \
        -S "$BIGWIG_BASE/GABA_Nestin-Ctrl.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
        --samplesLabel "Nestin-Ctrl" "Emx1-Mut" \
        --binSize 50 \
        --sortRegions keep \
        --missingDataAsZero \
        -o $OUTPUT_DIR/matrix_nestin_ctrl_vs_emx1_mut.gz \
        -p 16 \
        --outFileNameMatrix $OUTPUT_DIR/matrix_nestin_ctrl_vs_emx1_mut.tab

    # Comparison 3: Nestin-Mut vs Emx1-Mut
    echo "Computing matrix: Nestin-Mut vs Emx1-Mut..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 2000 -a 2000 \
        -R $BED_ALL \
        -S "$BIGWIG_BASE/GABA_Nestin-Mut.bw" "$BIGWIG_BASE/GABA_Emx1-Mut.bw" \
        --samplesLabel "Nestin-Mut" "Emx1-Mut" \
        --binSize 50 \
        --sortRegions keep \
        --missingDataAsZero \
        -o $OUTPUT_DIR/matrix_nestin_mut_vs_emx1_mut.gz \
        -p 16 \
        --outFileNameMatrix $OUTPUT_DIR/matrix_nestin_mut_vs_emx1_mut.tab

    # Run Python visualization
    echo "Creating custom comparison visualizations..."
    python 6_visualize_custom_comparisons.py --skip-individual

    if [ $? -ne 0 ]; then
        echo "ERROR: Step 6 failed!"
        PIPELINE_SUCCESS=false
    else
        echo ""
        echo "Step 6 completed successfully"
    fi
    echo ""
fi

################################################################################
# FINAL SUMMARY
################################################################################
echo "========================================================================"
echo "PIPELINE COMPLETE"
echo "========================================================================"
echo ""
echo "Finished: $(date)"
echo ""

if [ "$PIPELINE_SUCCESS" = true ]; then
    echo "STATUS: SUCCESS - All steps completed successfully"
    echo ""
    echo "Output directories:"
    echo "  - output/                          (TSV and BED files)"
    echo "  - output/heatmaps_deeptools/       (heatmaps, metaprofiles, matrices)"
    echo "  - output/bigwig_profiles_*/        (BigWig visualizations)"
    echo "  - output/custom_comparisons/       (cross-genotype comparisons)"
    echo ""
    echo "Key output files:"
    echo "  - CREs_splicing_genes_encode_all.bed"
    echo "  - heatmap_all.png, heatmap_GABA_nestin.png, heatmap_GABA_emx1.png"
    echo "  - metaprofile_*.png"
    echo "  - custom_comparisons/profiles/*.png"
else
    echo "STATUS: FAILED - One or more steps encountered errors"
    echo "Check logs/0_pipeline_master.err for details"
fi

echo ""
echo "NOTE: Emx1-Ctrl excluded (failed sample) - Nestin-Ctrl used as baseline"
echo "========================================================================"
