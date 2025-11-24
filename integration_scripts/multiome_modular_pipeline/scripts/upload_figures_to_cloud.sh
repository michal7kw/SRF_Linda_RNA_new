#!/bin/bash
#
# upload_figures_to_cloud.sh
#
# Uploads analysis figures to cloud storage (S3/GCS) for easy stakeholder access
# Run after pipeline completion to sync latest results
#
# Usage:
#   ./upload_figures_to_cloud.sh [s3|gcs|dry-run]
#
# Prerequisites:
#   - AWS CLI: pip install awscli && aws configure
#   - OR Google Cloud: pip install google-cloud-storage && gcloud auth login
#

set -euo pipefail

# Configuration
PIPELINE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
RESULTS_DIR="${PIPELINE_DIR}/signac_results"
DATE=$(date +%Y%m%d_%H%M%S)
STORAGE_TYPE="${1:-s3}"  # s3, gcs, or dry-run

# Storage configuration (UPDATE THESE!)
S3_BUCKET="s3://your-bucket-name/srf-multiome"
GCS_BUCKET="gs://your-bucket-name/srf-multiome"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to upload to AWS S3
upload_to_s3() {
    local bucket="$1"

    log_info "Checking AWS CLI availability..."
    if ! command -v aws &> /dev/null; then
        log_error "AWS CLI not found. Install with: pip install awscli"
        return 1
    fi

    log_info "Uploading to S3: ${bucket}"

    # Upload main plots
    log_info "Uploading main plots..."
    aws s3 sync "${RESULTS_DIR}/plots/" "${bucket}/plots/latest/" \
        --exclude "*" --include "*.pdf" --include "*.png" \
        --acl public-read \
        --metadata last_updated="${DATE}"

    # Upload cell type results
    log_info "Uploading cell type results..."
    aws s3 sync "${RESULTS_DIR}/celltype_results/filtered/" "${bucket}/celltype_results/latest/" \
        --exclude "*" --include "*.pdf" --include "*.png" \
        --acl public-read

    # Upload to dated archive
    log_info "Creating dated archive..."
    aws s3 sync "${RESULTS_DIR}/plots/" "${bucket}/plots/${DATE}/" \
        --exclude "*" --include "*.pdf" --include "*.png" \
        --acl public-read

    # Generate and upload index file
    generate_index_html "${bucket}"

    log_info "Upload complete!"
    log_info "View figures at: https://$(echo ${bucket} | sed 's|s3://||').s3.amazonaws.com/plots/latest/"
}

# Function to upload to Google Cloud Storage
upload_to_gcs() {
    local bucket="$1"

    log_info "Checking gsutil availability..."
    if ! command -v gsutil &> /dev/null; then
        log_error "gsutil not found. Install Google Cloud SDK"
        return 1
    fi

    log_info "Uploading to GCS: ${bucket}"

    # Upload main plots
    gsutil -m rsync -r -x ".*\.(?!pdf$|png$)" \
        "${RESULTS_DIR}/plots/" "${bucket}/plots/latest/"

    # Set public read permissions
    gsutil -m acl ch -u AllUsers:R "${bucket}/plots/latest/*"

    log_info "Upload complete!"
    log_info "View figures at: https://storage.googleapis.com/$(echo ${bucket} | sed 's|gs://||')/plots/latest/"
}

# Generate simple HTML index
generate_index_html() {
    local bucket="$1"
    local index_file="/tmp/index.html"

    cat > "${index_file}" <<'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>SRF Multiome Analysis - Latest Figures</title>
    <style>
        body { font-family: Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; background: #f5f5f5; }
        h1 { color: #333; border-bottom: 3px solid #4CAF50; padding-bottom: 10px; }
        .update-info { background: #e3f2fd; padding: 10px; border-radius: 5px; margin: 20px 0; }
        .figure-section { background: white; margin: 20px 0; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        .figure-section h2 { color: #1976D2; margin-top: 0; }
        .figure-link { display: inline-block; margin: 10px 10px 10px 0; padding: 10px 20px; background: #4CAF50; color: white; text-decoration: none; border-radius: 5px; }
        .figure-link:hover { background: #45a049; }
        .description { color: #666; margin: 10px 0; }
    </style>
</head>
<body>
    <h1>SRF Multiome Integration - Analysis Figures</h1>

    <div class="update-info">
        <strong>Last Updated:</strong> DATE_PLACEHOLDER<br>
        <strong>Pipeline:</strong> Signac-based multiome integration (RNA + ATAC)<br>
        <strong>Samples:</strong> Nestin-Ctrl, Nestin-Mut, Emx1-Ctrl, Emx1-Mut
    </div>

    <div class="figure-section">
        <h2>Differential Expression Analysis</h2>
        <p class="description">Cell-type-specific DEGs comparing Ctrl vs Mut, stratified by genotype (Nestin/Emx1)</p>
        <a href="03_DEG_volcano_plots.pdf" class="figure-link">View Volcano Plots (PDF)</a>
    </div>

    <div class="figure-section">
        <h2>Differential Accessibility Analysis</h2>
        <p class="description">Cell-type-specific DA peaks comparing Ctrl vs Mut, stratified by genotype</p>
        <a href="04_DA_volcano_plots.pdf" class="figure-link">View DA Plots (PDF)</a>
    </div>

    <div class="figure-section">
        <h2>Peak-Gene Regulatory Linkages</h2>
        <p class="description">Summary of chromatin accessibility peaks linked to gene expression (within 200kb, r>0.2)</p>
        <a href="06_peak_gene_link_summary.pdf" class="figure-link">View Linkage Summary (PDF)</a>
    </div>

    <div class="figure-section">
        <h2>Summary Statistics</h2>
        <p class="description">Overall analysis summary across all cell types</p>
        <a href="07_summary_barplot.pdf" class="figure-link">View Summary (PDF)</a>
    </div>

    <hr style="margin: 40px 0;">

    <h2>Cell Type-Specific Results</h2>
    <p>Detailed results available in <a href="../celltype_results/latest/">celltype_results</a> directory</p>

    <footer style="margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd; color: #666; font-size: 0.9em;">
        Generated automatically by upload_figures_to_cloud.sh<br>
        Repository: SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline
    </footer>
</body>
</html>
EOF

    # Replace date placeholder
    sed -i "s/DATE_PLACEHOLDER/${DATE}/" "${index_file}"

    # Upload index
    if [[ "$bucket" == s3://* ]]; then
        aws s3 cp "${index_file}" "${bucket}/plots/latest/index.html" --acl public-read --content-type "text/html"
    elif [[ "$bucket" == gs://* ]]; then
        gsutil -h "Content-Type:text/html" cp "${index_file}" "${bucket}/plots/latest/index.html"
        gsutil acl ch -u AllUsers:R "${bucket}/plots/latest/index.html"
    fi

    rm "${index_file}"
}

# Dry run - just list what would be uploaded
dry_run() {
    log_info "DRY RUN - Files that would be uploaded:"
    echo ""
    echo "=== Main Plots ==="
    find "${RESULTS_DIR}/plots/" -type f \( -name "*.pdf" -o -name "*.png" \) -exec ls -lh {} \;
    echo ""
    echo "=== Cell Type Results ==="
    find "${RESULTS_DIR}/celltype_results/filtered/" -type f \( -name "*.pdf" -o -name "*.png" \) 2>/dev/null -exec ls -lh {} \; || echo "No filtered results found"
    echo ""

    # Calculate total size
    local total_size=$(find "${RESULTS_DIR}/plots/" -type f \( -name "*.pdf" -o -name "*.png" \) -exec du -ch {} + | tail -1 | cut -f1)
    log_info "Total size: ${total_size}"
}

# Main execution
main() {
    cd "${PIPELINE_DIR}"

    log_info "SRF Multiome Figure Upload Utility"
    log_info "Storage type: ${STORAGE_TYPE}"

    # Check if results exist
    if [ ! -d "${RESULTS_DIR}/plots" ]; then
        log_error "Results directory not found: ${RESULTS_DIR}/plots"
        log_error "Have you run the pipeline yet?"
        exit 1
    fi

    case "${STORAGE_TYPE}" in
        s3)
            if [ "${S3_BUCKET}" = "s3://your-bucket-name/srf-multiome" ]; then
                log_error "Please update S3_BUCKET variable in script with your actual bucket name"
                log_warn "Run 'aws s3 mb s3://your-unique-bucket-name' to create a bucket first"
                exit 1
            fi
            upload_to_s3 "${S3_BUCKET}"
            ;;
        gcs)
            if [ "${GCS_BUCKET}" = "gs://your-bucket-name/srf-multiome" ]; then
                log_error "Please update GCS_BUCKET variable in script with your actual bucket name"
                exit 1
            fi
            upload_to_gcs "${GCS_BUCKET}"
            ;;
        dry-run)
            dry_run
            ;;
        *)
            log_error "Unknown storage type: ${STORAGE_TYPE}"
            echo "Usage: $0 [s3|gcs|dry-run]"
            exit 1
            ;;
    esac
}

main "$@"
