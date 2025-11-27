#!/bin/bash
################################################################################
# MASTER SCRIPT: GABA DEG ENCODE cCRE Analysis Pipeline
#
# This script runs the complete GABA DEG ENCODE cCRE analysis pipeline using
# SLURM job dependencies to ensure proper execution order.
#
# PIPELINE OVERVIEW:
# 1. Extract ENCODE cCREs overlapping with CREs linked to GABA DEGs
# 2. Compute signal matrices using deepTools
# 3. Create custom comparison visualizations
#
# COMPARISONS (excluding Emx1-Ctrl due to quality issues):
# - Nestin-Ctrl vs Nestin-Mut (within-genotype mutation effect)
# - Nestin-Ctrl vs Emx1-Mut (cross-genotype mutation effect)
# - Nestin-Mut vs Emx1-Mut (mutant genotype comparison)
#
# USAGE:
#   ./0_RUN_GABA_DEG_ENCODE_ANALYSIS.sh [options]
#
# OPTIONS:
#   --skip-individual    Skip individual CRE plots in visualization
#   --parallel N         Use N parallel processes for individual plots
#   --dry-run            Show what would be run without submitting jobs
#
# REQUIREMENTS:
# - GABA DEG lists (up/down regulated)
# - Table 16 CRE-gene correlations
# - ENCODE mm10-cCREs.bed
# - BigWig files from Signac pipeline
# - bedtools, deepTools (via conda)
#
################################################################################

set -e  # Exit on error

echo "========================================================================"
echo "GABA DEG ENCODE cCRE ANALYSIS PIPELINE"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Change to script directory
cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/GABA_DEGs_encode"

# Parse command line arguments
SKIP_INDIVIDUAL=""
PARALLEL_ARG=""
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-individual)
            SKIP_INDIVIDUAL="--skip-individual"
            shift
            ;;
        --parallel)
            PARALLEL_ARG="--parallel $2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create directories
mkdir -p logs
mkdir -p output

echo "Configuration:"
echo "  Skip individual plots: ${SKIP_INDIVIDUAL:-no}"
echo "  Parallel processes: ${PARALLEL_ARG:-default}"
echo "  Dry run: $DRY_RUN"
echo ""

# ============================================================================
# Check prerequisites
# ============================================================================
echo "Checking prerequisites..."

# DEG files
DEG_UP="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/biomarkers/sig_deg_lists/GABA/Cluster_GABA_vs_Rest_up_significant.csv"
DEG_DOWN="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/biomarkers/sig_deg_lists/GABA/Cluster_GABA_vs_Rest_down_significant.csv"

if [ ! -f "$DEG_UP" ]; then
    echo "ERROR: Up-DEG file not found: $DEG_UP"
    exit 1
fi
if [ ! -f "$DEG_DOWN" ]; then
    echo "ERROR: Down-DEG file not found: $DEG_DOWN"
    exit 1
fi
echo "  DEG files found"

# Table 16
TABLE16="../data/table_16.txt"
if [ ! -f "$TABLE16" ]; then
    echo "ERROR: Table 16 not found: $TABLE16"
    exit 1
fi
echo "  Table 16 found"

# ENCODE cCREs
ENCODE_CCRES="../data/mm10-cCREs.bed"
if [ ! -f "$ENCODE_CCRES" ]; then
    echo "ERROR: ENCODE cCREs not found: $ENCODE_CCRES"
    exit 1
fi
echo "  ENCODE cCREs found"

# BigWig files
BIGWIG_BASE="../../signac_results_L1/bigwig_tracks_L1/by_celltype"
for BW in GABA_Nestin-Ctrl.bw GABA_Nestin-Mut.bw GABA_Emx1-Mut.bw; do
    if [ ! -f "$BIGWIG_BASE/$BW" ]; then
        echo "ERROR: BigWig file not found: $BIGWIG_BASE/$BW"
        exit 1
    fi
done
echo "  BigWig files found"

echo ""
echo "All prerequisites satisfied!"
echo ""

# ============================================================================
# Create SLURM job scripts
# ============================================================================

# Step 1: Extract ENCODE cCREs
cat > 1_extract_encode_cCREs_for_DEGs.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=1_extract_DEGs_encode
#SBATCH --output=logs/1_extract_encode_cCREs_for_DEGs.log
#SBATCH --error=logs/1_extract_encode_cCREs_for_DEGs.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=workq

echo "Starting Step 1: Extract ENCODE cCREs for GABA DEGs"
echo "Date: $(date)"

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/GABA_DEGs_encode"

# Activate conda environment with bedtools (sc-chromatin2 has both bedtools and pandas)
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate sc-chromatin2

# Verify bedtools is available
echo "Checking bedtools: $(which bedtools)"

# Run extraction script
python 1_extract_encode_cCREs_for_DEGs.py

echo "Step 1 completed: $(date)"
EOF

# Step 3 wrapper (visualization)
cat > 3_visualize_deg_comparisons.sh << EOF
#!/bin/bash
#SBATCH --job-name=3_visualize_DEGs_encode
#SBATCH --output=logs/3_visualize_deg_comparisons.log
#SBATCH --error=logs/3_visualize_deg_comparisons.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=workq

echo "Starting Step 3: Create custom comparison visualizations"
echo "Date: \$(date)"

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/GABA_DEGs_encode"

# Activate conda environment (sc-chromatin2 has matplotlib, seaborn, scipy)
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate sc-chromatin2

# Run visualization script
python 3_visualize_deg_comparisons.py $SKIP_INDIVIDUAL $PARALLEL_ARG

echo "Step 3 completed: \$(date)"
EOF

# ============================================================================
# Submit jobs with dependencies
# ============================================================================

if [ "$DRY_RUN" = true ]; then
    echo "DRY RUN - Would submit the following jobs:"
    echo ""
    echo "Step 1: sbatch 1_extract_encode_cCREs_for_DEGs.sh"
    echo "Step 2: sbatch --dependency=afterok:\$JOB1 2_compute_signal_matrices.sh"
    echo "Step 3: sbatch --dependency=afterok:\$JOB2 3_visualize_deg_comparisons.sh"
    echo ""
    exit 0
fi

echo "========================================================================"
echo "SUBMITTING JOBS"
echo "========================================================================"
echo ""

# Step 1: Extract ENCODE cCREs
echo "Submitting Step 1: Extract ENCODE cCREs..."
JOB1=$(sbatch 1_extract_encode_cCREs_for_DEGs.sh | awk '{print $4}')
echo "  Job ID: $JOB1"

# Step 2: Compute signal matrices (depends on Step 1)
echo "Submitting Step 2: Compute signal matrices..."
JOB2=$(sbatch --dependency=afterok:$JOB1 2_compute_signal_matrices.sh | awk '{print $4}')
echo "  Job ID: $JOB2 (depends on $JOB1)"

# Step 3: Visualization (depends on Step 2)
echo "Submitting Step 3: Create visualizations..."
JOB3=$(sbatch --dependency=afterok:$JOB2 3_visualize_deg_comparisons.sh | awk '{print $4}')
echo "  Job ID: $JOB3 (depends on $JOB2)"

echo ""
echo "========================================================================"
echo "PIPELINE SUBMITTED!"
echo "========================================================================"
echo ""
echo "Job chain:"
echo "  Step 1 (Extract): $JOB1"
echo "  Step 2 (Matrix):  $JOB2 (after $JOB1)"
echo "  Step 3 (Viz):     $JOB3 (after $JOB2)"
echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo "  tail -f logs/*.log"
echo ""
echo "Output will be in:"
echo "  ./output/                    - TSV and BED files"
echo "  ./output/heatmaps_deeptools/ - Heatmaps and metaprofiles"
echo "  ./output/custom_comparisons/ - Custom comparison plots"
echo ""
echo "========================================================================"
