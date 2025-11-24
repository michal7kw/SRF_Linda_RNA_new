#!/bin/bash
#
# upload_figures_gcs.sh
#
# Automated Google Cloud Storage upload script for SRF project figures
# Uploads analysis figures to GCS for easy stakeholder access
#
# Usage:
#   ./upload_figures_gcs.sh [command]
#
# Commands:
#   upload          - Full upload of all figures (default)
#   quick           - Quick sync (only new/changed files)
#   dry-run         - Show what would be uploaded without uploading
#   list            - List uploaded files
#   get-url         - Get public URL for dashboard
#   stats           - Show storage usage and costs
#   download <path> - Download all figures to local path
#   clean-old       - Remove old archives (keep last N)
#   delete-all      - Delete all uploaded content (CAREFUL!)
#   update-index    - Regenerate HTML index page
#   test            - Test GCS connection and permissions
#   help            - Show this help message
#

set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source configuration
if [[ ! -f "${SCRIPT_DIR}/gcs_config.sh" ]]; then
    echo "ERROR: Configuration file not found: ${SCRIPT_DIR}/gcs_config.sh"
    exit 1
fi

source "${SCRIPT_DIR}/gcs_config.sh"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log() {
    local level="$1"
    shift
    local message="$*"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')

    # Console output with color
    case "$level" in
        INFO)  echo -e "${GREEN}[INFO]${NC} ${message}" ;;
        WARN)  echo -e "${YELLOW}[WARN]${NC} ${message}" ;;
        ERROR) echo -e "${RED}[ERROR]${NC} ${message}" ;;
        DEBUG) [[ "${LOG_LEVEL}" == "DEBUG" ]] && echo -e "${BLUE}[DEBUG]${NC} ${message}" ;;
    esac

    # File logging
    echo "${timestamp} [${level}] ${message}" >> "${LOG_FILE}"
}

log_info() { log INFO "$@"; }
log_warn() { log WARN "$@"; }
log_error() { log ERROR "$@"; }
log_debug() { log DEBUG "$@"; }

# Progress bar function
show_progress() {
    local current=$1
    local total=$2
    local prefix=$3
    local bar_width=50

    local progress=$((current * 100 / total))
    local filled=$((current * bar_width / total))
    local empty=$((bar_width - filled))

    printf "\r${prefix} ["
    printf "%${filled}s" | tr ' ' '='
    printf "%${empty}s" | tr ' ' '-'
    printf "] %d%% (%d/%d)" "$progress" "$current" "$total"
}

# Validate configuration
validate_setup() {
    log_info "Validating configuration..."

    if ! validate_config; then
        log_error "Configuration validation failed"
        exit 1
    fi

    # Test gcloud authentication
    if ! gcloud auth list --filter=status:ACTIVE --format="value(account)" &>/dev/null; then
        log_error "Not authenticated with gcloud. Run: gcloud init"
        exit 1
    fi

    # Test bucket access
    if ! gsutil ls "gs://${GCS_BUCKET}" &>/dev/null; then
        log_error "Cannot access bucket: gs://${GCS_BUCKET}"
        log_error "Check: 1) Bucket exists, 2) You have permissions, 3) gcloud project is set"
        exit 1
    fi

    log_info "Configuration validated successfully"
}

# Get current timestamp for archives
get_timestamp() {
    date '+%Y%m%d_%H%M%S'
}

# Upload single file with metadata
upload_file() {
    local local_file="$1"
    local gcs_path="$2"

    local file_size=$(stat -c%s "$local_file" 2>/dev/null || stat -f%z "$local_file" 2>/dev/null)

    # Validate file size
    if [[ $file_size -lt $MIN_FILE_SIZE ]]; then
        log_warn "Skipping small file (${file_size} bytes): $(basename "$local_file")"
        return 1
    fi

    if [[ $file_size -gt $MAX_FILE_SIZE ]]; then
        log_warn "Skipping large file (${file_size} bytes): $(basename "$local_file")"
        return 1
    fi

    # Upload with metadata
    gsutil -h "Cache-Control:${CACHE_CONTROL}" \
           -h "Content-Type:application/pdf" \
           cp "$local_file" "$gcs_path" 2>&1 | tee -a "${LOG_FILE}" | grep -v "^$" || true

    return 0
}

# Full upload of all figures
do_upload() {
    local create_archive="${1:-true}"

    log_info "========================================="
    log_info "Starting full upload to GCS"
    log_info "Bucket: gs://${GCS_BUCKET}"
    log_info "========================================="

    local timestamp=$(get_timestamp)
    local upload_count=0
    local error_count=0

    # Check if results exist
    if [[ ! -d "${RESULTS_DIR}/plots" ]]; then
        log_error "Results directory not found: ${RESULTS_DIR}/plots"
        log_error "Have you run the pipeline yet?"
        exit 1
    fi

    # Upload main plots (L2 cell types)
    log_info "Uploading L2 cell type plots..."
    if gsutil -m rsync -r -d \
        -x ".*\.(log|tmp|bak)$" \
        "${RESULTS_DIR}/plots/" \
        "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/plots/latest/"; then
        log_info "L2 plots uploaded successfully"
        ((upload_count++))
    else
        log_error "Failed to upload L2 plots"
        ((error_count++))
    fi

    # Upload L1 plots if they exist
    if [[ -d "${RESULTS_DIR_L1}/plots" ]]; then
        log_info "Uploading L1 cell type plots..."
        if gsutil -m rsync -r -d \
            "${RESULTS_DIR_L1}/plots/" \
            "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/plots_L1/latest/"; then
            log_info "L1 plots uploaded successfully"
            ((upload_count++))
        else
            log_warn "Failed to upload L1 plots (continuing anyway)"
        fi
    fi

    # Upload cell type results (filtered)
    if [[ -d "${RESULTS_DIR}/celltype_results/filtered" ]]; then
        log_info "Uploading filtered cell type results..."
        if gsutil -m rsync -r \
            -x ".*\.(rds|RData)$" \
            "${RESULTS_DIR}/celltype_results/filtered/" \
            "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/celltype_results/latest/filtered/"; then
            log_info "Cell type results uploaded successfully"
            ((upload_count++))
        else
            log_warn "Failed to upload cell type results (continuing anyway)"
        fi
    fi

    # Upload peak-gene linkage results
    if [[ -d "${RESULTS_DIR}/peak_gene_linkage_analysis_UPDATED" ]]; then
        log_info "Uploading peak-gene linkage results..."
        if gsutil -m rsync -r \
            "${RESULTS_DIR}/peak_gene_linkage_analysis_UPDATED/" \
            "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/peak_gene_linkage/latest/"; then
            log_info "Peak-gene linkage results uploaded successfully"
            ((upload_count++))
        else
            log_warn "Failed to upload peak-gene linkage results (continuing anyway)"
        fi
    fi

    # Upload DEG coverage plots
    if [[ -d "${RESULTS_DIR}/DEG_coverage_enhanced" ]]; then
        log_info "Uploading DEG coverage plots..."
        if gsutil -m rsync -r \
            "${RESULTS_DIR}/DEG_coverage_enhanced/" \
            "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/DEG_coverage/latest/"; then
            log_info "DEG coverage plots uploaded successfully"
            ((upload_count++))
        else
            log_warn "Failed to upload DEG coverage plots (continuing anyway)"
        fi
    fi

    # Create dated archive if enabled
    if [[ "${ENABLE_ARCHIVES}" == "true" ]] && [[ "${create_archive}" == "true" ]]; then
        log_info "Creating dated archive: ${timestamp}"
        if gsutil -m rsync -r \
            "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/plots/latest/" \
            "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/plots/archive_${timestamp}/"; then
            log_info "Archive created successfully"
        else
            log_warn "Failed to create archive (continuing anyway)"
        fi
    fi

    # Set public access permissions
    log_info "Setting public access permissions..."
    if gsutil -m acl ch -r -u AllUsers:R "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/" 2>/dev/null; then
        log_info "Permissions set successfully"
    else
        log_warn "Could not set public permissions (files may not be accessible)"
    fi

    # Generate and upload index page
    generate_index_page

    log_info "========================================="
    log_info "Upload complete!"
    log_info "Uploaded: ${upload_count} directories"
    log_info "Errors: ${error_count}"
    log_info "========================================="

    # Send notification if configured
    send_notification "Upload completed" "Uploaded ${upload_count} directories with ${error_count} errors"

    # Show public URL
    get_public_url

    return 0
}

# Quick sync (only changed files)
do_quick_sync() {
    log_info "Quick sync: uploading only new/changed files..."

    # Use rsync with checksum comparison
    gsutil -m rsync -r -c \
        "${RESULTS_DIR}/plots/" \
        "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/plots/latest/"

    log_info "Quick sync complete"
    get_public_url
}

# Dry run - show what would be uploaded
do_dry_run() {
    log_info "========================================="
    log_info "DRY RUN - Files that would be uploaded:"
    log_info "========================================="

    echo ""
    echo "=== L2 Cell Type Plots ==="
    if [[ -d "${RESULTS_DIR}/plots" ]]; then
        find "${RESULTS_DIR}/plots/" -type f \( -name "*.pdf" -o -name "*.png" \) -exec ls -lh {} \;
    else
        echo "  Directory not found: ${RESULTS_DIR}/plots"
    fi

    echo ""
    echo "=== L1 Cell Type Plots ==="
    if [[ -d "${RESULTS_DIR_L1}/plots" ]]; then
        find "${RESULTS_DIR_L1}/plots/" -type f \( -name "*.pdf" -o -name "*.png" \) -exec ls -lh {} \;
    else
        echo "  Directory not found: ${RESULTS_DIR_L1}/plots"
    fi

    echo ""
    echo "=== Filtered Cell Type Results ==="
    if [[ -d "${RESULTS_DIR}/celltype_results/filtered" ]]; then
        find "${RESULTS_DIR}/celltype_results/filtered/" -type f \( -name "*.pdf" -o -name "*.png" -o -name "*.csv" \) -exec ls -lh {} \;
    else
        echo "  Directory not found: ${RESULTS_DIR}/celltype_results/filtered"
    fi

    echo ""
    echo "=== Storage Summary ==="
    if [[ -d "${RESULTS_DIR}/plots" ]]; then
        local total_size=$(find "${RESULTS_DIR}/plots/" -type f \( -name "*.pdf" -o -name "*.png" \) -exec du -ch {} + 2>/dev/null | tail -1 | cut -f1)
        echo "  Total size: ${total_size}"
    fi

    echo ""
    log_info "Run './upload_figures_gcs.sh upload' to perform actual upload"
}

# List uploaded files
do_list() {
    log_info "Listing uploaded files in gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/"
    echo ""

    gsutil ls -lh "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/plots/latest/" 2>/dev/null || {
        log_warn "No files found or bucket not accessible"
        return 1
    }
}

# Get public URL
get_public_url() {
    local base_url="https://storage.googleapis.com/${GCS_BUCKET}"

    echo ""
    log_info "========================================="
    log_info "Public URLs for sharing:"
    log_info "========================================="
    echo ""
    echo "Dashboard (index page):"
    echo "  ${base_url}/index.html"
    echo ""
    echo "Latest plots directory:"
    echo "  ${base_url}/${GCS_PATH_MULTIOME}/plots/latest/"
    echo ""
    echo "Specific figures:"
    echo "  ${base_url}/${GCS_PATH_MULTIOME}/plots/latest/03_DEG_volcano_plots.pdf"
    echo "  ${base_url}/${GCS_PATH_MULTIOME}/plots/latest/04_DA_volcano_plots.pdf"
    echo "  ${base_url}/${GCS_PATH_MULTIOME}/plots/latest/06_peak_gene_link_summary.pdf"
    echo ""
}

# Show storage statistics and estimated costs
do_stats() {
    log_info "Calculating storage usage and costs..."
    echo ""

    # Get bucket size
    local bucket_info=$(gsutil du -s "gs://${GCS_BUCKET}" 2>/dev/null)
    local bucket_bytes=$(echo "$bucket_info" | awk '{print $1}')
    local bucket_gb=$(echo "scale=2; $bucket_bytes / 1024 / 1024 / 1024" | bc)

    echo "Storage Usage:"
    echo "  Total size: ${bucket_gb} GB ($(numfmt --to=iec-i --suffix=B $bucket_bytes))"
    echo ""

    # Estimate costs (GCS Standard storage)
    local storage_cost=$(echo "scale=2; $bucket_gb * 0.020" | bc)

    echo "Estimated Monthly Costs (GCS Standard, us-central1):"
    echo "  Storage: \$${storage_cost}/month (${bucket_gb} GB × \$0.020/GB)"
    echo "  Operations: ~\$0.00/month (first 5000 reads free)"
    echo "  Network egress: ~\$0.10-0.50/month (depends on usage)"
    echo "  ----------------------------------------"
    echo "  Estimated total: \$$(echo "scale=2; $storage_cost + 0.30" | bc)/month"
    echo ""

    # List archives
    echo "Archives:"
    gsutil ls "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/plots/" 2>/dev/null | grep "archive_" | wc -l | xargs -I {} echo "  {} dated archives"
    echo ""
}

# Download all figures to local path
do_download() {
    local dest_path="${1:-.}"

    log_info "Downloading figures to: ${dest_path}"

    mkdir -p "${dest_path}"

    gsutil -m rsync -r \
        "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/plots/latest/" \
        "${dest_path}/"

    log_info "Download complete"
}

# Clean old archives
do_clean_old() {
    if [[ "${MAX_ARCHIVES}" -eq 0 ]]; then
        log_info "MAX_ARCHIVES=0, keeping all archives"
        return 0
    fi

    log_info "Cleaning old archives (keeping last ${MAX_ARCHIVES})..."

    # Get list of archives sorted by date
    local archives=($(gsutil ls "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/plots/" | grep "archive_" | sort -r))
    local num_archives=${#archives[@]}

    if [[ $num_archives -le $MAX_ARCHIVES ]]; then
        log_info "Only ${num_archives} archives exist, nothing to clean"
        return 0
    fi

    # Delete old archives
    local to_delete=$((num_archives - MAX_ARCHIVES))
    log_info "Deleting ${to_delete} old archives..."

    for ((i=MAX_ARCHIVES; i<num_archives; i++)); do
        local archive="${archives[$i]}"
        log_info "Deleting: $(basename "$archive")"
        gsutil -m rm -r "$archive"
    done

    log_info "Cleanup complete"
}

# Delete all uploaded content (DANGEROUS!)
do_delete_all() {
    echo ""
    log_warn "========================================="
    log_warn "WARNING: This will DELETE ALL uploaded figures!"
    log_warn "Bucket: gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/"
    log_warn "========================================="
    echo ""
    read -p "Type 'DELETE' to confirm: " confirm

    if [[ "$confirm" != "DELETE" ]]; then
        log_info "Cancelled"
        return 0
    fi

    log_info "Deleting all content..."
    gsutil -m rm -r "gs://${GCS_BUCKET}/${GCS_PATH_MULTIOME}/"
    log_info "All content deleted"
}

# Generate HTML index page
generate_index_page() {
    log_info "Generating HTML index page..."

    local index_file="/tmp/srf_index_$$.html"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')

    cat > "${index_file}" <<EOF
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SRF Multiome Analysis - Figure Dashboard</title>
    <style>
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            line-height: 1.6;
            color: #333;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 12px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
        }
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }
        .header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.2);
        }
        .header p {
            font-size: 1.2em;
            opacity: 0.9;
        }
        .update-info {
            background: #e3f2fd;
            padding: 20px;
            border-left: 4px solid #2196F3;
            margin: 20px 40px;
            border-radius: 4px;
        }
        .update-info strong {
            color: #1976D2;
        }
        .content {
            padding: 0 40px 40px;
        }
        .section {
            margin: 30px 0;
            padding: 25px;
            background: #f8f9fa;
            border-radius: 8px;
            border: 1px solid #e0e0e0;
        }
        .section h2 {
            color: #1976D2;
            margin-bottom: 15px;
            font-size: 1.8em;
            border-bottom: 3px solid #2196F3;
            padding-bottom: 10px;
        }
        .section h3 {
            color: #555;
            margin: 20px 0 10px;
            font-size: 1.3em;
        }
        .description {
            color: #666;
            margin: 10px 0;
            font-size: 1.05em;
        }
        .figure-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }
        .figure-link {
            display: block;
            padding: 15px 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            text-decoration: none;
            border-radius: 6px;
            font-weight: 500;
            transition: all 0.3s ease;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        .figure-link:hover {
            transform: translateY(-2px);
            box-shadow: 0 6px 12px rgba(0,0,0,0.2);
            background: linear-gradient(135deg, #764ba2 0%, #667eea 100%);
        }
        .figure-link::before {
            content: "📊 ";
            margin-right: 8px;
        }
        .folder-link {
            background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
        }
        .folder-link::before {
            content: "📁 ";
        }
        .info-box {
            background: #fff3cd;
            border: 1px solid #ffc107;
            border-radius: 6px;
            padding: 15px;
            margin: 20px 0;
        }
        .info-box strong {
            color: #856404;
        }
        footer {
            margin-top: 40px;
            padding: 30px;
            text-align: center;
            background: #f1f1f1;
            color: #666;
            font-size: 0.9em;
        }
        .badge {
            display: inline-block;
            padding: 4px 12px;
            background: #4CAF50;
            color: white;
            border-radius: 12px;
            font-size: 0.85em;
            margin-left: 10px;
            font-weight: 600;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 SRF Multiome Analysis</h1>
            <p>Integrated scRNA-seq + scATAC-seq Figure Dashboard</p>
        </div>

        <div class="update-info">
            <strong>Last Updated:</strong> ${timestamp}<br>
            <strong>Pipeline:</strong> Signac-based multiome integration (RNA + ATAC)<br>
            <strong>Samples:</strong> R26-Nestin-Ctrl, R26-Nestin-Mut, R26-Emx1-Ctrl, R26-Emx1-Mut<br>
            <strong>Cell Type Resolution:</strong> L2 (fine-grained, 20+ subtypes) + L1 (broad, 5-10 categories)
        </div>

        <div class="content">
            <div class="section">
                <h2>📈 Main Analysis Figures (L2 Cell Types)</h2>
                <p class="description">
                    Comprehensive analysis of fine-grained cell types comparing Control vs Mutant within each genotype (Nestin, Emx1).
                </p>

                <h3>Differential Expression (DEGs)</h3>
                <p class="description">Cell-type-specific differentially expressed genes</p>
                <div class="figure-grid">
                    <a href="${GCS_PATH_MULTIOME}/plots/latest/03_DEG_volcano_plots.pdf" class="figure-link">DEG Volcano Plots</a>
                </div>

                <h3>Differential Accessibility (DA Peaks)</h3>
                <p class="description">Cell-type-specific chromatin accessibility changes</p>
                <div class="figure-grid">
                    <a href="${GCS_PATH_MULTIOME}/plots/latest/04_DA_volcano_plots.pdf" class="figure-link">DA Volcano Plots</a>
                </div>

                <h3>Peak-Gene Regulatory Linkages</h3>
                <p class="description">Chromatin peaks linked to gene expression (200kb window, correlation r>0.2)</p>
                <div class="figure-grid">
                    <a href="${GCS_PATH_MULTIOME}/plots/latest/06_peak_gene_link_summary.pdf" class="figure-link">Linkage Summary</a>
                </div>

                <h3>Summary Statistics</h3>
                <p class="description">Overall analysis overview across all cell types</p>
                <div class="figure-grid">
                    <a href="${GCS_PATH_MULTIOME}/plots/latest/07_summary_barplot.pdf" class="figure-link">Summary Barplot</a>
                </div>
            </div>

            <div class="section">
                <h2>🔬 Detailed Results by Cell Type</h2>
                <p class="description">
                    Publication-ready filtered results for each cell type (FDR < 0.05, stringent quality thresholds)
                </p>
                <div class="figure-grid">
                    <a href="${GCS_PATH_MULTIOME}/celltype_results/latest/filtered/" class="figure-link folder-link">Browse Cell Type Results</a>
                </div>

                <div class="info-box">
                    <strong>Available Data:</strong> DEGs (CSV), DA peaks (CSV), Peak-gene links (CSV), Quality metrics, Cell counts
                </div>
            </div>

            <div class="section">
                <h2>🎯 Specialized Analyses</h2>

                <h3>Peak-Gene Linkage Analysis</h3>
                <p class="description">Detailed regulatory linkage analysis with coverage plots</p>
                <div class="figure-grid">
                    <a href="${GCS_PATH_MULTIOME}/peak_gene_linkage/latest/" class="figure-link folder-link">Peak-Gene Linkage Results</a>
                </div>

                <h3>DEG Coverage Plots</h3>
                <p class="description">Chromatin accessibility tracks for differentially expressed genes</p>
                <div class="figure-grid">
                    <a href="${GCS_PATH_MULTIOME}/DEG_coverage/latest/" class="figure-link folder-link">DEG Coverage Plots</a>
                </div>
            </div>

            <div class="section">
                <h2>📚 Additional Resources</h2>
                <div class="figure-grid">
                    <a href="${GCS_PATH_MULTIOME}/plots_L1/latest/" class="figure-link folder-link">L1 Cell Types (Broad Categories)</a>
                    <a href="${GCS_PATH_MULTIOME}/plots/" class="figure-link folder-link">All Versions (Historical Archives)</a>
                </div>
            </div>

            <div class="info-box">
                <strong>💡 Tip:</strong> All PDFs can be viewed directly in your browser or downloaded.
                Right-click any link and select "Save link as..." to download.
            </div>
        </div>

        <footer>
            <p>
                <strong>SRF Multiome Integration Pipeline</strong><br>
                Generated automatically by upload_figures_gcs.sh<br>
                Repository: SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline<br>
                <br>
                Questions? Contact your bioinformatics team.
            </p>
        </footer>
    </div>
</body>
</html>
EOF

    # Upload index page
    gsutil -h "Content-Type:text/html" \
           -h "Cache-Control:max-age=300, public" \
           cp "${index_file}" "gs://${GCS_BUCKET}/index.html"

    gsutil acl ch -u AllUsers:R "gs://${GCS_BUCKET}/index.html" 2>/dev/null || true

    rm "${index_file}"
    log_info "Index page uploaded to: https://storage.googleapis.com/${GCS_BUCKET}/index.html"
}

# Test GCS connection and permissions
do_test() {
    log_info "Testing GCS connection and permissions..."
    echo ""

    # Test 1: gcloud auth
    echo "Test 1: gcloud authentication"
    if gcloud auth list --filter=status:ACTIVE --format="value(account)" &>/dev/null; then
        local account=$(gcloud auth list --filter=status:ACTIVE --format="value(account)")
        echo "  ✓ Authenticated as: ${account}"
    else
        echo "  ✗ Not authenticated. Run: gcloud init"
        return 1
    fi

    # Test 2: Project configuration
    echo "Test 2: Project configuration"
    local current_project=$(gcloud config get-value project 2>/dev/null)
    if [[ -n "$current_project" ]]; then
        echo "  ✓ Project: ${current_project}"
    else
        echo "  ✗ No project set. Run: gcloud config set project ${GCS_PROJECT}"
        return 1
    fi

    # Test 3: Bucket access
    echo "Test 3: Bucket access"
    if gsutil ls "gs://${GCS_BUCKET}" &>/dev/null; then
        echo "  ✓ Can access bucket: gs://${GCS_BUCKET}"
    else
        echo "  ✗ Cannot access bucket: gs://${GCS_BUCKET}"
        echo "     Create bucket: gsutil mb -p ${GCS_PROJECT} -l us-central1 gs://${GCS_BUCKET}"
        return 1
    fi

    # Test 4: Write permissions
    echo "Test 4: Write permissions"
    local test_file="/tmp/gcs_test_$$.txt"
    echo "test" > "${test_file}"
    if gsutil cp "${test_file}" "gs://${GCS_BUCKET}/test_upload.txt" &>/dev/null; then
        echo "  ✓ Can write to bucket"
        gsutil rm "gs://${GCS_BUCKET}/test_upload.txt" &>/dev/null
    else
        echo "  ✗ Cannot write to bucket. Check IAM permissions."
        return 1
    fi
    rm "${test_file}"

    # Test 5: Public access
    echo "Test 5: Public access configuration"
    local test_url="https://storage.googleapis.com/${GCS_BUCKET}/test.txt"
    echo "test" | gsutil -h "Content-Type:text/plain" cp - "gs://${GCS_BUCKET}/test.txt" 2>/dev/null
    gsutil acl ch -u AllUsers:R "gs://${GCS_BUCKET}/test.txt" 2>/dev/null
    if curl -sf "${test_url}" &>/dev/null; then
        echo "  ✓ Public access is working"
    else
        echo "  ⚠ Public access may not be configured"
        echo "     Run: gsutil iam ch allUsers:objectViewer gs://${GCS_BUCKET}"
    fi
    gsutil rm "gs://${GCS_BUCKET}/test.txt" &>/dev/null

    echo ""
    log_info "All tests passed! Ready to upload."
    echo ""
    get_public_url
}

# Send notification (email or Slack)
send_notification() {
    local title="$1"
    local message="$2"

    # Email notification
    if [[ -n "${NOTIFICATION_EMAIL}" ]]; then
        log_debug "Sending email notification to ${NOTIFICATION_EMAIL}"
        echo "${message}" | mail -s "SRF GCS Upload: ${title}" "${NOTIFICATION_EMAIL}" 2>/dev/null || true
    fi

    # Slack notification
    if [[ -n "${SLACK_WEBHOOK_URL}" ]]; then
        log_debug "Sending Slack notification"
        curl -X POST -H 'Content-type: application/json' \
            --data "{\"text\":\"*${title}*\n${message}\"}" \
            "${SLACK_WEBHOOK_URL}" 2>/dev/null || true
    fi
}

# Show help
show_help() {
    cat <<EOF
SRF Google Cloud Storage Figure Upload Tool

Usage: ./upload_figures_gcs.sh [command]

Commands:
  upload          Full upload of all figures (creates archive)
  quick           Quick sync (only new/changed files, no archive)
  dry-run         Show what would be uploaded without uploading
  list            List uploaded files
  get-url         Get public URLs for sharing
  stats           Show storage usage and estimated costs
  download <dir>  Download all figures to local directory
  clean-old       Remove old archives (keep last N, configured in gcs_config.sh)
  delete-all      Delete all uploaded content (DANGEROUS!)
  update-index    Regenerate HTML index page
  test            Test GCS connection and permissions
  help            Show this help message

Examples:
  # First time setup
  gcloud init
  ./upload_figures_gcs.sh test
  ./upload_figures_gcs.sh upload

  # Daily usage
  ./upload_figures_gcs.sh quick

  # Check what's uploaded
  ./upload_figures_gcs.sh list
  ./upload_figures_gcs.sh stats

  # Get URL for stakeholders
  ./upload_figures_gcs.sh get-url

Configuration:
  Edit scripts/gcs_config.sh to configure bucket name, paths, and options

Documentation:
  See GCS_SETUP_GUIDE.md for complete setup instructions
EOF
}

# Main execution
main() {
    local command="${1:-upload}"

    # Show header
    echo ""
    log_info "========================================="
    log_info "SRF GCS Figure Upload Tool"
    log_info "========================================="
    echo ""

    # Handle commands that don't need validation
    case "$command" in
        help)
            show_help
            exit 0
            ;;
    esac

    # Validate setup for all other commands
    validate_setup

    # Execute command
    case "$command" in
        upload)
            do_upload true
            ;;
        quick)
            do_quick_sync
            ;;
        dry-run)
            do_dry_run
            ;;
        list)
            do_list
            ;;
        get-url)
            get_public_url
            ;;
        stats)
            do_stats
            ;;
        download)
            do_download "${2:-.}"
            ;;
        clean-old)
            do_clean_old
            ;;
        delete-all)
            do_delete_all
            ;;
        update-index)
            generate_index_page
            ;;
        test)
            do_test
            ;;
        *)
            log_error "Unknown command: $command"
            echo ""
            show_help
            exit 1
            ;;
    esac

    echo ""
}

# Run main function
main "$@"
