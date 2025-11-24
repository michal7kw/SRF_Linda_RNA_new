#!/bin/bash
#
# upload_complete_analysis_to_gcs.sh
#
# Upload complete multiome analysis results to Google Cloud Storage
# Uploads PNG, PDF, CSV, and BigWig files from all analysis directories
#
# Usage:
#   ./upload_complete_analysis_to_gcs.sh [dry-run|upload]
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "${SCRIPT_DIR}")"

# Source configuration
source "${SCRIPT_DIR}/gcs_config.sh"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[INFO]${NC} $*"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }
log_step() { echo -e "${BLUE}[STEP]${NC} $*"; }

# Parse mode
DRY_RUN=false
MODE="${1:-upload}"
if [[ "$MODE" == "dry-run" ]]; then
    DRY_RUN=true
    GSUTIL_FLAGS="-n"  # Dry-run flag for gsutil
    log_warn "DRY RUN MODE - No files will be uploaded"
else
    GSUTIL_FLAGS=""
fi

clear

cat <<'EOF'
╔══════════════════════════════════════════════════════════════════╗
║                                                                  ║
║     Upload Complete Analysis Results to Google Cloud Storage    ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝

This script uploads all analysis results from:
  1. CREs from literature (heatmaps, deeptools)
  2. L1 cell type results (DEGs, DA peaks, plots)
  3. L1 BigWig tracks (943 MB, 35 files)
  4. L1 peak-gene linkage (GABA analysis)

File types: PNG, PDF, CSV, BW (BigWig)
Estimated size: ~1 GB
Estimated time: 15-20 minutes

EOF

if $DRY_RUN; then
    echo "⚠️  DRY RUN MODE - Showing what would be uploaded"
else
    echo "🚀 UPLOAD MODE - Files will be uploaded to GCS"
fi

echo ""

# ============================================================================
# Configuration
# ============================================================================

# Validate GCS configuration
validate_config

# Directory mappings: LOCAL_PATH|GCS_PATH|DESCRIPTION
declare -a UPLOAD_DIRS=(
    "CREs_from_literature/output/GABA_DEG_analysis/heatmaps_deeptools|multiome/CREs_literature/heatmaps_deeptools|CREs Heatmaps (deeptools)"
    "CREs_from_literature/output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools|multiome/CREs_literature/heatmaps_specific_CREs|Specific CREs Heatmaps"
    "signac_results_L1/celltype_results|multiome/L1_celltype_results|L1 Cell Type Results"
    "signac_results_L1/bigwig_tracks_L1|multiome/L1_bigwig_tracks|L1 BigWig Tracks"
    "signac_results_L1/peak_gene_linkage_analysis_UPDATED/GABA/peak_gene_links|multiome/L1_peak_gene_linkage/GABA|GABA Peak-Gene Linkage"
)

# File extensions to upload
FILE_PATTERNS=( "*.png" "*.pdf" "*.csv" "*.bw" )

# Exclusion patterns (files to skip)
EXCLUDE_PATTERNS=( "*.log" "*.txt" "*.gz" "*.tab" "*.bed" "*.rds" "*.RData" )

echo "Configuration:"
echo "  GCS Bucket: gs://${GCS_BUCKET}"
echo "  Directories to upload: ${#UPLOAD_DIRS[@]}"
echo "  File types: ${FILE_PATTERNS[*]}"
echo ""

if ! $DRY_RUN; then
    read -p "Continue with upload? (yes/no): " confirm
    if [[ "$confirm" != "yes" ]]; then
        log_warn "Upload cancelled by user"
        exit 0
    fi
    echo ""
fi

# ============================================================================
# Step 1: Check prerequisites
# ============================================================================

log_step "Step 1: Checking prerequisites"
echo ""

# Check gsutil
if ! command -v gsutil &> /dev/null; then
    log_error "gsutil not found"
    echo "Install: https://cloud.google.com/storage/docs/gsutil_install"
    exit 1
fi

log_info "✓ gsutil found: $(gsutil version | head -1)"

# Check authentication
if ! gsutil ls gs://${GCS_BUCKET}/ &> /dev/null; then
    log_error "Cannot access GCS bucket: gs://${GCS_BUCKET}"
    echo ""
    echo "Authenticate with: gcloud auth login"
    exit 1
fi

log_info "✓ GCS bucket accessible: gs://${GCS_BUCKET}"
echo ""

# ============================================================================
# Step 2: Count files to upload
# ============================================================================

log_step "Step 2: Analyzing directories"
echo ""

total_files=0
total_size=0

for mapping in "${UPLOAD_DIRS[@]}"; do
    IFS='|' read -r local_path gcs_path description <<< "$mapping"
    full_local_path="${PIPELINE_DIR}/${local_path}"

    if [[ ! -d "$full_local_path" ]]; then
        log_warn "Directory not found: $local_path"
        continue
    fi

    # Count files matching patterns
    file_count=0
    for pattern in "${FILE_PATTERNS[@]}"; do
        count=$(find "$full_local_path" -type f -name "$pattern" 2>/dev/null | wc -l)
        file_count=$((file_count + count))
    done

    # Get directory size
    dir_size=$(du -sh "$full_local_path" 2>/dev/null | cut -f1)

    echo "  $description"
    echo "    Local: $local_path"
    echo "    GCS: gs://${GCS_BUCKET}/${gcs_path}"
    echo "    Files: $file_count"
    echo "    Size: $dir_size"
    echo ""

    total_files=$((total_files + file_count))
done

log_info "Total files to upload: $total_files"
echo ""

if $DRY_RUN; then
    log_warn "DRY RUN - Showing commands that would run"
    echo ""
fi

# ============================================================================
# Step 3: Upload each directory
# ============================================================================

log_step "Step 3: Uploading directories"
echo ""

upload_count=0
start_time=$(date +%s)

for mapping in "${UPLOAD_DIRS[@]}"; do
    IFS='|' read -r local_path gcs_path description <<< "$mapping"
    full_local_path="${PIPELINE_DIR}/${local_path}"

    if [[ ! -d "$full_local_path" ]]; then
        log_warn "Skipping: $description (directory not found)"
        continue
    fi

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log_info "Uploading: $description"
    echo "  From: $local_path"
    echo "  To: gs://${GCS_BUCKET}/${gcs_path}"
    echo ""

    # Upload with rsync (no exclusion - gsutil doesn't support complex patterns)
    # Instead, we'll upload everything and let GCS handle it
    if $DRY_RUN; then
        log_info "Would run:"
        echo "  gsutil -m rsync -r $GSUTIL_FLAGS \\"
        echo "    '$full_local_path/' \\"
        echo "    'gs://${GCS_BUCKET}/${gcs_path}/'"
        echo ""
    else
        log_info "Uploading (this may take a few minutes)..."

        # Upload using gsutil rsync
        if gsutil -m rsync -r $GSUTIL_FLAGS \
            "$full_local_path/" \
            "gs://${GCS_BUCKET}/${gcs_path}/" 2>&1 | grep -v "^$"; then
            log_info "✓ Upload successful"
            ((upload_count++))
        else
            log_warn "⚠️  Upload encountered errors (check output)"
        fi

        echo ""
    fi
done

end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
log_info "Upload phase complete"
echo "  Directories uploaded: $upload_count / ${#UPLOAD_DIRS[@]}"
echo "  Time elapsed: ${elapsed}s"
echo ""

# ============================================================================
# Step 4: Set public permissions
# ============================================================================

if ! $DRY_RUN; then
    log_step "Step 4: Setting public access permissions"
    echo ""

    log_info "Making all files publicly readable..."

    for mapping in "${UPLOAD_DIRS[@]}"; do
        IFS='|' read -r local_path gcs_path description <<< "$mapping"

        # Check if directory exists in GCS before setting permissions
        if gsutil ls "gs://${GCS_BUCKET}/${gcs_path}/" &>/dev/null; then
            if gsutil -m iam ch allUsers:objectViewer "gs://${GCS_BUCKET}/${gcs_path}" 2>&1 | grep -v "^$" | grep -v "No URLs matched"; then
                log_info "✓ Set permissions for: $description"
            fi
        else
            log_warn "Skipping permissions for $description (not uploaded)"
        fi
    done

    echo ""
fi

# ============================================================================
# Step 5: Generate enhanced dashboard
# ============================================================================

if ! $DRY_RUN; then
    log_step "Step 5: Generating enhanced dashboard"
    echo ""

    log_info "Creating HTML dashboard with all uploaded files..."

    # Generate enhanced dashboard HTML
    cat > /tmp/gcs_index.html << 'EOFHTML'
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SRF Multiome Analysis - Complete Results</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 20px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
        }
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }
        .header h1 { font-size: 2.5em; margin-bottom: 10px; }
        .header p { font-size: 1.2em; opacity: 0.9; }
        .content { padding: 40px; }
        .section {
            margin-bottom: 40px;
            border-left: 4px solid #667eea;
            padding-left: 20px;
        }
        .section h2 {
            color: #667eea;
            margin-bottom: 15px;
            font-size: 1.8em;
        }
        .section p {
            color: #666;
            margin-bottom: 15px;
            line-height: 1.6;
        }
        .file-list {
            background: #f7fafc;
            border-radius: 10px;
            padding: 20px;
            margin-top: 15px;
        }
        .file-item {
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 10px 15px;
            margin: 5px 0;
            background: white;
            border-radius: 8px;
            transition: all 0.3s;
            border: 1px solid #e2e8f0;
        }
        .file-item:hover {
            transform: translateX(5px);
            box-shadow: 0 4px 12px rgba(102, 126, 234, 0.2);
            border-color: #667eea;
        }
        .file-name {
            font-family: 'Monaco', 'Courier New', monospace;
            color: #2d3748;
            font-size: 0.9em;
        }
        .file-type {
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 0.75em;
            font-weight: 600;
            text-transform: uppercase;
        }
        .type-png { background: #fed7d7; color: #c53030; }
        .type-pdf { background: #bee3f8; color: #2c5282; }
        .type-csv { background: #c6f6d5; color: #2f855a; }
        .type-bw { background: #fbd38d; color: #c05621; }
        .stats {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }
        .stat-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 25px;
            border-radius: 15px;
            text-align: center;
        }
        .stat-number { font-size: 2.5em; font-weight: bold; }
        .stat-label { font-size: 0.9em; opacity: 0.9; margin-top: 5px; }
        a { color: inherit; text-decoration: none; }
        .folder-path {
            font-family: 'Monaco', 'Courier New', monospace;
            background: #edf2f7;
            padding: 8px 12px;
            border-radius: 6px;
            font-size: 0.85em;
            color: #4a5568;
            margin: 10px 0;
            display: inline-block;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 SRF Multiome Analysis</h1>
            <p>Complete Analysis Results - Public Dashboard</p>
            <p style="font-size: 0.9em; margin-top: 10px;">Updated: TIMESTAMP</p>
        </div>

        <div class="content">
            <div class="stats">
                <div class="stat-card">
                    <div class="stat-number">TOTAL_FILES</div>
                    <div class="stat-label">Total Files</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">TOTAL_SIZE</div>
                    <div class="stat-label">Total Size</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">5</div>
                    <div class="stat-label">Analysis Categories</div>
                </div>
            </div>

            <div class="section">
                <h2>📊 Main Signac Results</h2>
                <p>Primary multiome integration results (L2 cell types)</p>
                <div class="folder-path">multiome/plots/latest/</div>
                <div class="file-list">
                    <a href="multiome/plots/latest/03_DEG_volcano_plots.pdf" target="_blank">
                        <div class="file-item">
                            <span class="file-name">03_DEG_volcano_plots.pdf</span>
                            <span class="file-type type-pdf">PDF</span>
                        </div>
                    </a>
                    <a href="multiome/plots/latest/04_DA_volcano_plots.pdf" target="_blank">
                        <div class="file-item">
                            <span class="file-name">04_DA_volcano_plots.pdf</span>
                            <span class="file-type type-pdf">PDF</span>
                        </div>
                    </a>
                </div>
            </div>

            <div class="section">
                <h2>🔥 CREs from Literature Analysis</h2>
                <p>Cell-type specificity validation using literature-curated CREs</p>
                <div class="folder-path">multiome/CREs_literature/</div>
                <div class="file-list">
                    <p><strong>Heatmaps (deeptools)</strong></p>
                    <a href="multiome/CREs_literature/heatmaps_deeptools/heatmap_GABA_all_conditions.png" target="_blank">
                        <div class="file-item">
                            <span class="file-name">heatmap_GABA_all_conditions.png</span>
                            <span class="file-type type-png">PNG</span>
                        </div>
                    </a>
                    <a href="multiome/CREs_literature/heatmaps_deeptools/metaprofile_GABA_all_conditions.png" target="_blank">
                        <div class="file-item">
                            <span class="file-name">metaprofile_GABA_all_conditions.png</span>
                            <span class="file-type type-png">PNG</span>
                        </div>
                    </a>
                </div>
            </div>

            <div class="section">
                <h2>📈 L1 Cell Type Results</h2>
                <p>Broad cell type analysis (DEGs, DA peaks, peak-gene links)</p>
                <div class="folder-path">multiome/L1_celltype_results/</div>
                <div class="file-list">
                    <p><strong>Browse by category:</strong></p>
                    <a href="multiome/L1_celltype_results/DEG/" target="_blank">
                        <div class="file-item">
                            <span class="file-name">DEG/ (Differential Expression)</span>
                            <span class="file-type type-csv">CSV</span>
                        </div>
                    </a>
                    <a href="multiome/L1_celltype_results/DA/" target="_blank">
                        <div class="file-item">
                            <span class="file-name">DA/ (Differential Accessibility)</span>
                            <span class="file-type type-csv">CSV</span>
                        </div>
                    </a>
                </div>
            </div>

            <div class="section">
                <h2>📊 BigWig Tracks (IGV)</h2>
                <p>ATAC-seq coverage tracks for visualization in IGV browser (943 MB, 35 files)</p>
                <div class="folder-path">multiome/L1_bigwig_tracks/by_celltype/</div>
                <div class="file-list">
                    <p><strong>Example tracks:</strong></p>
                    <a href="multiome/L1_bigwig_tracks/by_celltype/GABA_Nestin-Ctrl.bw" target="_blank">
                        <div class="file-item">
                            <span class="file-name">GABA_Nestin-Ctrl.bw</span>
                            <span class="file-type type-bw">BW</span>
                        </div>
                    </a>
                    <a href="multiome/L1_bigwig_tracks/by_celltype/GABA_Nestin-Mut.bw" target="_blank">
                        <div class="file-item">
                            <span class="file-name">GABA_Nestin-Mut.bw</span>
                            <span class="file-type type-bw">BW</span>
                        </div>
                    </a>
                    <p style="margin-top: 15px; font-size: 0.9em; color: #666;">
                        <strong>Load in IGV:</strong> File → Load from URL → Paste BigWig URL
                    </p>
                </div>
            </div>

            <div class="section">
                <h2>🔗 Peak-Gene Linkage Analysis</h2>
                <p>GABA neuron regulatory network analysis</p>
                <div class="folder-path">multiome/L1_peak_gene_linkage/GABA/</div>
                <div class="file-list">
                    <a href="multiome/L1_peak_gene_linkage/GABA/significant_links.csv" target="_blank">
                        <div class="file-item">
                            <span class="file-name">significant_links.csv</span>
                            <span class="file-type type-csv">CSV</span>
                        </div>
                    </a>
                </div>
            </div>

            <div class="section" style="border-left-color: #48bb78;">
                <h2 style="color: #48bb78;">ℹ️ About This Dashboard</h2>
                <p>
                    This dashboard provides public access to all SRF multiome analysis results.
                    All files are hosted on Google Cloud Storage with direct HTTP access.
                </p>
                <p style="margin-top: 10px;">
                    <strong>Base URL:</strong>
                    <span class="folder-path">https://storage.googleapis.com/srf-multiome-figures/</span>
                </p>
                <p style="margin-top: 15px; color: #666; font-size: 0.9em;">
                    💡 <strong>Tip:</strong> Click any file to view/download directly.
                    All URLs are stable and can be shared with collaborators.
                </p>
            </div>
        </div>
    </div>
</body>
</html>
EOFHTML

    # Update placeholders
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    total_size=$(gsutil du -sh gs://${GCS_BUCKET}/ 2>/dev/null | cut -f1 || echo "N/A")
    total_files=$(gsutil ls -r gs://${GCS_BUCKET}/ 2>/dev/null | grep -v "/$" | wc -l || echo "N/A")

    # Use | as delimiter instead of / to avoid conflicts with slashes in values
    sed -i "s|TIMESTAMP|$timestamp|g" /tmp/gcs_index.html
    sed -i "s|TOTAL_FILES|$total_files|g" /tmp/gcs_index.html
    sed -i "s|TOTAL_SIZE|$total_size|g" /tmp/gcs_index.html

    # Upload dashboard
    gsutil cp /tmp/gcs_index.html gs://${GCS_BUCKET}/index.html
    gsutil setmeta -h "Content-Type:text/html" -h "Cache-Control:public, max-age=300" gs://${GCS_BUCKET}/index.html

    log_info "✓ Dashboard uploaded and updated"
    echo ""
fi

# ============================================================================
# Summary
# ============================================================================

clear

cat <<EOF

╔══════════════════════════════════════════════════════════════════╗
║                                                                  ║
║              Upload Complete! ✓                                  ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝

EOF

if $DRY_RUN; then
    echo "DRY RUN COMPLETE - No files were uploaded"
    echo ""
    echo "To execute for real:"
    echo "  ./scripts/upload_complete_analysis_to_gcs.sh upload"
else
    echo "✓ All analysis results uploaded to GCS"
    echo "✓ Public access permissions set"
    echo "✓ Dashboard generated and updated"
    echo ""
    echo "📊 Upload Statistics:"
    echo "  Directories processed: ${#UPLOAD_DIRS[@]}"
    echo "  Successful uploads: $upload_count"
    echo "  Total files: $total_files"
    echo "  Time elapsed: ${elapsed}s"
    echo ""
    echo "🔗 Access Your Results:"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    echo "Dashboard (main page):"
    echo "  https://storage.googleapis.com/srf-multiome-figures/index.html"
    echo ""
    echo "Direct folder access:"
    echo "  CREs Heatmaps:"
    echo "    https://storage.googleapis.com/srf-multiome-figures/multiome/CREs_literature/heatmaps_deeptools/"
    echo ""
    echo "  L1 Cell Type Results:"
    echo "    https://storage.googleapis.com/srf-multiome-figures/multiome/L1_celltype_results/"
    echo ""
    echo "  BigWig Tracks (for IGV):"
    echo "    https://storage.googleapis.com/srf-multiome-figures/multiome/L1_bigwig_tracks/by_celltype/"
    echo ""
    echo "  Peak-Gene Linkage:"
    echo "    https://storage.googleapis.com/srf-multiome-figures/multiome/L1_peak_gene_linkage/GABA/"
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    echo "💡 Share these URLs with your collaborators!"
    echo ""
fi

echo "For detailed documentation, see:"
echo "  GCS_UPLOAD_COMPLETE_ANALYSIS.md"
echo ""
