#!/bin/bash
# =============================================================================
# SCRIPT: 6_intersect_peaks_with_PLS.sh
# PURPOSE: Intersect ATAC-seq peaks with PLS (Promoter-Like Signature) elements
#
# DESCRIPTION:
# This script identifies ATAC-seq peaks that overlap with ENCODE PLS elements
# (promoter-like cis-regulatory elements) for each experimental condition.
# The resulting peak sets will be used for TFBS (Transcription Factor Binding Site)
# motif analysis to identify which TFs may regulate splicing gene expression.
#
# INPUT:
# - PLS BED files from ENCODE cCRE extraction (Step 2 output)
# - ATAC-seq peak files (MACS2 narrowPeak format)
#
# OUTPUT:
# - Condition-specific BED files of peaks overlapping PLS elements
# - Summary statistics of peak-PLS overlaps
#
# NOTE: Emx1-Ctrl is excluded as a failed sample (per project documentation)
# =============================================================================
#SBATCH --job-name=intersect_PLS
#SBATCH --output=logs/6_intersect_PLS_%j.out
#SBATCH --error=logs/6_intersect_PLS_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=workq

set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

# Base directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/output/tfbs_analysis"
LOG_DIR="${SCRIPT_DIR}/logs"

# Input files
PLS_BED="${SCRIPT_DIR}/output/encode_cCREs_PLS.bed"
PLS_CTCF_BED="${SCRIPT_DIR}/output/encode_cCREs_PLS_CTCF_bound.bed"
ATAC_PEAK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_ATAC_data/chromap_final_output"

# Sample configuration (Emx1-Ctrl excluded - failed sample)
declare -A SAMPLES=(
    ["Nestin_Ctrl"]="R26-Nestin-Ctrl-adult_peaks_sorted.bed"
    ["Nestin_Mut"]="R26-Nestin-Mut-adult_peaks_sorted.bed"
    ["Emx1_Mut"]="R26-Emx1-Mut-adult_peaks_sorted.bed"
)

# Conda environment with bedtools
CONDA_BASE="/beegfs/scratch/ric.broccoli/kubacki.michal/conda"
CONDA_ENV="/beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/peak_calling_new"

# =============================================================================
# FUNCTIONS
# =============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $1"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2
}

validate_file() {
    local file="$1"
    local desc="$2"
    if [[ ! -f "$file" ]]; then
        log_error "$desc not found: $file"
        exit 1
    fi
    log_info "Validated $desc: $file ($(wc -l < "$file") lines)"
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

log_info "=========================================="
log_info "Starting ATAC Peak - PLS Intersection Analysis"
log_info "=========================================="

# Create output directories
mkdir -p "${OUTPUT_DIR}/peaks_on_PLS"
mkdir -p "${OUTPUT_DIR}/statistics"
mkdir -p "${LOG_DIR}"

# Activate conda environment
log_info "Activating conda environment..."
source "${CONDA_BASE}/bin/activate" "${CONDA_ENV}"

# Validate bedtools
if ! command -v bedtools &> /dev/null; then
    log_error "bedtools not found in conda environment"
    exit 1
fi
log_info "Using bedtools: $(which bedtools)"

# Validate input files
validate_file "$PLS_BED" "PLS BED file"
validate_file "$PLS_CTCF_BED" "PLS+CTCF BED file"

# Combine PLS files for complete promoter-like element set
COMBINED_PLS="${OUTPUT_DIR}/PLS_combined.bed"
log_info "Combining PLS BED files..."
cat "$PLS_BED" "$PLS_CTCF_BED" | sort -k1,1 -k2,2n | bedtools merge -i - > "$COMBINED_PLS"
log_info "Combined PLS elements: $(wc -l < "$COMBINED_PLS") unique regions"

# Initialize summary file
SUMMARY_FILE="${OUTPUT_DIR}/statistics/intersection_summary.txt"
echo "# ATAC Peak - PLS Intersection Summary" > "$SUMMARY_FILE"
echo "# Generated: $(date)" >> "$SUMMARY_FILE"
echo "# PLS elements (combined): $(wc -l < "$COMBINED_PLS")" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
printf "%-15s %12s %12s %12s %10s\n" "Sample" "Total_Peaks" "PLS_Overlaps" "Unique_PLS" "Pct_Peaks" >> "$SUMMARY_FILE"
echo "------------------------------------------------------------------------" >> "$SUMMARY_FILE"

# Process each sample
for sample_name in "${!SAMPLES[@]}"; do
    peak_file="${ATAC_PEAK_DIR}/${SAMPLES[$sample_name]}"

    log_info "Processing $sample_name..."

    # Validate peak file
    if [[ ! -f "$peak_file" ]]; then
        log_error "Peak file not found for $sample_name: $peak_file"
        continue
    fi

    total_peaks=$(wc -l < "$peak_file")
    log_info "  Total ATAC peaks: $total_peaks"

    # Output files
    output_bed="${OUTPUT_DIR}/peaks_on_PLS/${sample_name}_peaks_on_PLS.bed"
    output_homer="${OUTPUT_DIR}/peaks_on_PLS/${sample_name}_peaks_on_PLS_homer.bed"

    # Find peaks overlapping with PLS elements
    # Using -u to report each peak only once even if it overlaps multiple PLS
    bedtools intersect -a "$peak_file" -b "$COMBINED_PLS" -u > "$output_bed"

    overlap_count=$(wc -l < "$output_bed")
    pct_overlap=$(awk "BEGIN {printf \"%.2f\", ($overlap_count / $total_peaks) * 100}")

    log_info "  Peaks overlapping PLS: $overlap_count ($pct_overlap%)"

    # Count unique PLS elements covered
    unique_pls=$(bedtools intersect -a "$COMBINED_PLS" -b "$peak_file" -u | wc -l)
    log_info "  Unique PLS elements covered: $unique_pls"

    # Create HOMER-compatible format (chr start end name score strand)
    # HOMER expects: ID in column 4, strand in column 6
    awk -v sample="$sample_name" 'BEGIN{OFS="\t"} {
        print $1, $2, $3, sample"_peak_"NR, $5, ($6 != "" ? $6 : "+")
    }' "$output_bed" > "$output_homer"

    log_info "  Created HOMER-compatible BED: $output_homer"

    # Write to summary
    printf "%-15s %12d %12d %12d %9.2f%%\n" "$sample_name" "$total_peaks" "$overlap_count" "$unique_pls" "$pct_overlap" >> "$SUMMARY_FILE"
done

# Create merged peak set for background comparison
log_info "Creating merged peak set from all samples..."
MERGED_PEAKS="${OUTPUT_DIR}/peaks_on_PLS/All_samples_merged_peaks_on_PLS.bed"
cat "${OUTPUT_DIR}/peaks_on_PLS/"*_peaks_on_PLS.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - > "$MERGED_PEAKS"
log_info "Merged peak set: $(wc -l < "$MERGED_PEAKS") unique regions"

# Create HOMER-compatible merged file
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "merged_peak_"NR, 1000, "+"}' "$MERGED_PEAKS" > "${OUTPUT_DIR}/peaks_on_PLS/All_samples_merged_peaks_on_PLS_homer.bed"

# Add footer to summary
echo "" >> "$SUMMARY_FILE"
echo "------------------------------------------------------------------------" >> "$SUMMARY_FILE"
echo "Merged unique regions: $(wc -l < "$MERGED_PEAKS")" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
echo "Output files:" >> "$SUMMARY_FILE"
echo "  - Individual sample BED files: ${OUTPUT_DIR}/peaks_on_PLS/{sample}_peaks_on_PLS.bed" >> "$SUMMARY_FILE"
echo "  - HOMER-compatible files: ${OUTPUT_DIR}/peaks_on_PLS/{sample}_peaks_on_PLS_homer.bed" >> "$SUMMARY_FILE"
echo "  - Merged peaks: $MERGED_PEAKS" >> "$SUMMARY_FILE"

log_info "=========================================="
log_info "Intersection Analysis Complete"
log_info "=========================================="
log_info "Summary saved to: $SUMMARY_FILE"
cat "$SUMMARY_FILE"

log_info "Done!"
