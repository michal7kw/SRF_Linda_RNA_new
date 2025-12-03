#!/bin/bash
# =============================================================================
# SCRIPT: 7_run_homer_tfbs_analysis.sh
# PURPOSE: Run HOMER motif analysis on ATAC peaks overlapping PLS elements
#
# DESCRIPTION:
# This script uses HOMER (findMotifsGenome.pl) to identify enriched transcription
# factor binding site (TFBS) motifs in the promoter-like ATAC peaks. It performs:
# 1. De novo motif discovery to find novel enriched motifs
# 2. Known motif enrichment analysis against HOMER's TF motif database
# 3. Separate analysis for each condition (Nestin-Ctrl, Nestin-Mut, Emx1-Mut)
#
# The analysis reveals which TFs are likely binding to promoter regions of
# splicing-related genes and how SRF mutation affects TF accessibility.
#
# INPUT:
# - HOMER-formatted BED files from Step 6 (peaks_on_PLS_homer.bed)
# - mm10 genome (installed via HOMER or provided as FASTA)
#
# OUTPUT:
# - De novo motifs with statistical enrichment
# - Known TF motif enrichment tables
# - HTML reports for interactive exploration
#
# DEPENDENCIES:
# - HOMER (findMotifsGenome.pl) with mm10 genome
# - Conda environment: annotation_enrichment
#
# NOTE: First run may take longer if HOMER needs to install mm10 genome
# =============================================================================
#SBATCH --job-name=homer_tfbs
#SBATCH --output=logs/7_homer_tfbs_%A_%a.out
#SBATCH --error=logs/7_homer_tfbs_%A_%a.err
#SBATCH --array=0-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --partition=workq

set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

# Base directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${SCRIPT_DIR}/output/tfbs_analysis/peaks_on_PLS"
OUTPUT_DIR="${SCRIPT_DIR}/output/tfbs_analysis/homer_results"
LOG_DIR="${SCRIPT_DIR}/logs"

# Conda environment with HOMER
CONDA_BASE="/beegfs/scratch/ric.broccoli/kubacki.michal/conda"
CONDA_ENV="/beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/annotation_enrichment"
HOMER_PATH="/beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/annotation_enrichment/share/homer"

# Genome configuration
GENOME="mm10"
# Alternative: use FASTA path if HOMER genome not installed
# GENOME_FASTA="/path/to/mm10.fa"

# Sample configuration (array job)
SAMPLES=("Nestin_Ctrl" "Nestin_Mut" "Emx1_Mut")
SAMPLE_IDX=${SLURM_ARRAY_TASK_ID:-0}
CURRENT_SAMPLE="${SAMPLES[$SAMPLE_IDX]}"

# HOMER parameters
MOTIF_SIZE="200"          # Size of region around peak center to analyze
MOTIF_LENGTHS="8,10,12"   # De novo motif lengths to search
NUM_MOTIFS="25"           # Number of de novo motifs to find
NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

# =============================================================================
# FUNCTIONS
# =============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$CURRENT_SAMPLE] INFO: $1"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$CURRENT_SAMPLE] ERROR: $1" >&2
}

log_warning() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$CURRENT_SAMPLE] WARNING: $1"
}

check_homer_genome() {
    local genome="$1"
    log_info "Checking HOMER genome installation for $genome..."

    # Check if genome is installed
    if perl "${HOMER_PATH}/configureHomer.pl" -list 2>&1 | grep -q "^+.*${genome}"; then
        log_info "HOMER genome $genome is installed"
        return 0
    else
        log_warning "HOMER genome $genome not installed. Attempting to install..."
        perl "${HOMER_PATH}/configureHomer.pl" -install "$genome"
        if [[ $? -eq 0 ]]; then
            log_info "Successfully installed HOMER genome $genome"
            return 0
        else
            log_error "Failed to install HOMER genome $genome"
            return 1
        fi
    fi
}

validate_file() {
    local file="$1"
    local desc="$2"
    if [[ ! -f "$file" ]]; then
        log_error "$desc not found: $file"
        exit 1
    fi
    local count=$(wc -l < "$file")
    log_info "Validated $desc: $file ($count regions)"
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

log_info "=========================================="
log_info "Starting HOMER TFBS Analysis"
log_info "Sample: $CURRENT_SAMPLE (Array index: $SAMPLE_IDX)"
log_info "=========================================="

# Create output directories
mkdir -p "${OUTPUT_DIR}/${CURRENT_SAMPLE}"
mkdir -p "${LOG_DIR}"

# Activate conda environment
log_info "Activating conda environment..."
source "${CONDA_BASE}/bin/activate" "${CONDA_ENV}"

# Verify HOMER installation
if ! command -v findMotifsGenome.pl &> /dev/null; then
    log_error "HOMER findMotifsGenome.pl not found"
    exit 1
fi
log_info "Using HOMER: $(which findMotifsGenome.pl)"

# Check/install HOMER genome
if ! check_homer_genome "$GENOME"; then
    log_error "Cannot proceed without HOMER genome"
    log_info "Alternative: Provide mm10 FASTA path via GENOME_FASTA variable"
    exit 1
fi

# Define input file
INPUT_BED="${INPUT_DIR}/${CURRENT_SAMPLE}_peaks_on_PLS_homer.bed"
validate_file "$INPUT_BED" "Input BED file"

# Get peak count for reporting
PEAK_COUNT=$(wc -l < "$INPUT_BED")
log_info "Running HOMER on $PEAK_COUNT peaks overlapping PLS elements"

# Define output directory for this sample
SAMPLE_OUTPUT="${OUTPUT_DIR}/${CURRENT_SAMPLE}"

# Run HOMER findMotifsGenome.pl
# Parameters:
#   -size: Fragment size around peak center
#   -len: Motif lengths for de novo discovery
#   -S: Number of motifs to optimize
#   -p: Number of parallel threads
#   -mask: Mask repeats
#   -preparsedDir: Use cached preparsed data for faster runs
log_info "Running HOMER findMotifsGenome.pl..."
log_info "  Size: $MOTIF_SIZE bp"
log_info "  Motif lengths: $MOTIF_LENGTHS"
log_info "  De novo motifs: $NUM_MOTIFS"
log_info "  Threads: $NUM_THREADS"

# Create preparsed directory for HOMER cache
PREPARSED_DIR="${OUTPUT_DIR}/preparsed_${GENOME}"
mkdir -p "$PREPARSED_DIR"

findMotifsGenome.pl \
    "$INPUT_BED" \
    "$GENOME" \
    "$SAMPLE_OUTPUT" \
    -size "$MOTIF_SIZE" \
    -len "$MOTIF_LENGTHS" \
    -S "$NUM_MOTIFS" \
    -p "$NUM_THREADS" \
    -mask \
    -preparsedDir "$PREPARSED_DIR" \
    2>&1 | tee "${LOG_DIR}/7_homer_${CURRENT_SAMPLE}_detailed.log"

# Check if HOMER completed successfully
if [[ ! -f "${SAMPLE_OUTPUT}/homerResults.html" ]]; then
    log_error "HOMER did not complete successfully - no output HTML found"
    exit 1
fi

log_info "HOMER analysis completed successfully"

# =============================================================================
# POST-PROCESSING: Extract key results
# =============================================================================

log_info "Extracting key results..."

# Create summary of top known motifs
if [[ -f "${SAMPLE_OUTPUT}/knownResults.txt" ]]; then
    log_info "Extracting top 20 known motifs..."
    head -21 "${SAMPLE_OUTPUT}/knownResults.txt" > "${SAMPLE_OUTPUT}/top20_known_motifs.txt"

    # Count significantly enriched motifs (p < 0.01)
    sig_motifs=$(awk 'NR>1 && $3 < 0.01 {count++} END {print count+0}' "${SAMPLE_OUTPUT}/knownResults.txt")
    log_info "Significantly enriched known motifs (p < 0.01): $sig_motifs"
fi

# Create summary of de novo motifs
if [[ -f "${SAMPLE_OUTPUT}/homerResults.html" ]]; then
    denovo_count=$(ls -1 "${SAMPLE_OUTPUT}"/homerResults/motif*.motif 2>/dev/null | wc -l || echo "0")
    log_info "De novo motifs found: $denovo_count"
fi

# Create a quick summary file
SUMMARY_FILE="${SAMPLE_OUTPUT}/analysis_summary.txt"
cat > "$SUMMARY_FILE" << EOF
# HOMER TFBS Analysis Summary
# Sample: $CURRENT_SAMPLE
# Generated: $(date)

## Input
- Peak file: $INPUT_BED
- Number of peaks: $PEAK_COUNT
- Genome: $GENOME

## Parameters
- Analysis size: $MOTIF_SIZE bp
- Motif lengths: $MOTIF_LENGTHS
- De novo motifs: $NUM_MOTIFS

## Results
- Known motifs analyzed: $(wc -l < "${SAMPLE_OUTPUT}/knownResults.txt" 2>/dev/null || echo "N/A")
- Significantly enriched (p < 0.01): ${sig_motifs:-N/A}
- De novo motifs found: ${denovo_count:-N/A}

## Output Files
- HTML report: ${SAMPLE_OUTPUT}/homerResults.html
- Known motif results: ${SAMPLE_OUTPUT}/knownResults.txt
- Top 20 known motifs: ${SAMPLE_OUTPUT}/top20_known_motifs.txt
- De novo motifs: ${SAMPLE_OUTPUT}/homerResults/

EOF

log_info "=========================================="
log_info "Analysis Complete for $CURRENT_SAMPLE"
log_info "=========================================="
log_info "Output directory: $SAMPLE_OUTPUT"
log_info "HTML report: ${SAMPLE_OUTPUT}/homerResults.html"
log_info "Known motifs: ${SAMPLE_OUTPUT}/knownResults.txt"

log_info "Done!"
