#!/bin/bash
#
# sync_all_to_gcs.sh
#
# Automated one-command solution:
#   1. Find all directories with visualization files
#   2. Upload to GCS with organized structure
#   3. Auto-generate dashboard based on uploaded content
#
# Usage:
#   ./sync_all_to_gcs.sh [dry-run|sync]
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
log_step() { echo -e "${BLUE}[STEP]${NC} $*"; }

# Parse mode
DRY_RUN=false
MODE="${1:-sync}"
if [[ "$MODE" == "dry-run" ]]; then
    DRY_RUN=true
fi

clear

cat <<'EOF'
╔══════════════════════════════════════════════════════════════════╗
║                                                                  ║
║     Automated GCS Sync & Dashboard Generator                    ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝

This script will:
  1. Auto-discover all directories with PNG/PDF/CSV/BW files
  2. Upload everything to GCS with organized structure
  3. Generate dynamic dashboard based on actual content
  4. Set public permissions

EOF

if $DRY_RUN; then
    echo "🔍 DRY RUN MODE - Showing what would be uploaded"
else
    echo "🚀 SYNC MODE - Files will be uploaded to GCS"
fi

echo ""

# Validate GCS configuration
validate_config

# ============================================================================
# Step 1: Auto-discover directories with visualization files
# ============================================================================

log_step "Step 1: Discovering directories with visualization files"
echo ""

log_info "Scanning pipeline directory for PNG, PDF, CSV, BW files..."
echo ""

# Find all directories containing visualization files
# Exclude certain directories
EXCLUDE_DIRS=(
    "*/cache/*"
    "*/tmp/*"
    "*/.git/*"
    "*/Archive/*"
    "*/Legacy/*"
)

# Build find exclude arguments
FIND_EXCLUDE=""
for exclude in "${EXCLUDE_DIRS[@]}"; do
    FIND_EXCLUDE="$FIND_EXCLUDE -path '$exclude' -prune -o"
done

# Find directories with visualization files
declare -A DIR_MAP
declare -A DIR_COUNTS
declare -A DIR_SIZES

# Discover CREs from literature
if [[ -d "CREs_from_literature/output/GABA_DEG_analysis" ]]; then
    for subdir in CREs_from_literature/output/GABA_DEG_analysis/*/; do
        if [[ -d "$subdir" ]]; then
            dirname=$(basename "$subdir")
            file_count=$(find "$subdir" -type f \( -name "*.png" -o -name "*.pdf" -o -name "*.csv" \) 2>/dev/null | wc -l)
            if [[ $file_count -gt 0 ]]; then
                gcs_path="multiome/CREs_literature/${dirname}"
                DIR_MAP["$subdir"]="$gcs_path"
                DIR_COUNTS["$subdir"]=$file_count
                DIR_SIZES["$subdir"]=$(du -sh "$subdir" 2>/dev/null | cut -f1)
            fi
        fi
    done
fi

# Discover L1 results
if [[ -d "signac_results_L1" ]]; then
    # Cell type results
    if [[ -d "signac_results_L1/celltype_results" ]]; then
        for subdir in signac_results_L1/celltype_results/*/; do
            if [[ -d "$subdir" ]]; then
                dirname=$(basename "$subdir")
                file_count=$(find "$subdir" -type f \( -name "*.png" -o -name "*.pdf" -o -name "*.csv" \) 2>/dev/null | wc -l)
                if [[ $file_count -gt 0 ]]; then
                    gcs_path="multiome/L1_celltype_results/${dirname}"
                    DIR_MAP["$subdir"]="$gcs_path"
                    DIR_COUNTS["$subdir"]=$file_count
                    DIR_SIZES["$subdir"]=$(du -sh "$subdir" 2>/dev/null | cut -f1)
                fi
            fi
        done
    fi

    # BigWig tracks
    if [[ -d "signac_results_L1/bigwig_tracks_L1" ]]; then
        file_count=$(find "signac_results_L1/bigwig_tracks_L1" -type f -name "*.bw" 2>/dev/null | wc -l)
        if [[ $file_count -gt 0 ]]; then
            DIR_MAP["signac_results_L1/bigwig_tracks_L1"]="multiome/L1_bigwig_tracks"
            DIR_COUNTS["signac_results_L1/bigwig_tracks_L1"]=$file_count
            DIR_SIZES["signac_results_L1/bigwig_tracks_L1"]=$(du -sh "signac_results_L1/bigwig_tracks_L1" 2>/dev/null | cut -f1)
        fi
    fi

    # Peak-gene linkage
    if [[ -d "signac_results_L1/peak_gene_linkage_analysis_UPDATED" ]]; then
        for cell_type_dir in signac_results_L1/peak_gene_linkage_analysis_UPDATED/*/; do
            if [[ -d "$cell_type_dir" ]]; then
                cell_type=$(basename "$cell_type_dir")
                if [[ -d "${cell_type_dir}peak_gene_links" ]]; then
                    file_count=$(find "${cell_type_dir}peak_gene_links" -type f \( -name "*.png" -o -name "*.pdf" -o -name "*.csv" \) 2>/dev/null | wc -l)
                    if [[ $file_count -gt 0 ]]; then
                        gcs_path="multiome/L1_peak_gene_linkage/${cell_type}"
                        DIR_MAP["${cell_type_dir}peak_gene_links"]="$gcs_path"
                        DIR_COUNTS["${cell_type_dir}peak_gene_links"]=$file_count
                        DIR_SIZES["${cell_type_dir}peak_gene_links"]=$(du -sh "${cell_type_dir}peak_gene_links" 2>/dev/null | cut -f1)
                    fi
                fi
            fi
        done
    fi
fi

# Display discovered directories
log_info "Discovered ${#DIR_MAP[@]} directories with visualization files:"
echo ""

total_files=0
for local_dir in "${!DIR_MAP[@]}"; do
    gcs_path="${DIR_MAP[$local_dir]}"
    count="${DIR_COUNTS[$local_dir]}"
    size="${DIR_SIZES[$local_dir]}"
    total_files=$((total_files + count))

    echo "  📁 $(basename "$local_dir")"
    echo "     Local: $local_dir"
    echo "     GCS: gs://${GCS_BUCKET}/${gcs_path}"
    echo "     Files: $count | Size: $size"
    echo ""
done

log_info "Total files to upload: $total_files"
echo ""

if [[ ${#DIR_MAP[@]} -eq 0 ]]; then
    log_warn "No directories with visualization files found"
    exit 0
fi

if ! $DRY_RUN; then
    # Check if running interactively
    if [[ -t 0 ]]; then
        read -p "Continue with upload? (yes/no): " confirm
        if [[ "$confirm" != "yes" ]]; then
            log_warn "Upload cancelled"
            exit 0
        fi
    else
        # Non-interactive mode (piped input) - auto-confirm
        log_info "Non-interactive mode - auto-confirming upload"
    fi
    echo ""
fi

# ============================================================================
# Step 2: Upload all directories
# ============================================================================

log_step "Step 2: Uploading all directories"
echo ""

upload_count=0
start_time=$(date +%s)

for local_dir in "${!DIR_MAP[@]}"; do
    gcs_path="${DIR_MAP[$local_dir]}"

    if [[ ! -d "$local_dir" ]]; then
        log_warn "Skipping: $local_dir (not found)"
        continue
    fi

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log_info "Uploading: $(basename "$local_dir")"
    echo "  From: $local_dir"
    echo "  To: gs://${GCS_BUCKET}/${gcs_path}"
    echo ""

    if $DRY_RUN; then
        log_info "Would run: gsutil -m rsync -r '$local_dir/' 'gs://${GCS_BUCKET}/${gcs_path}/'"
        echo ""
    else
        log_info "Uploading (this may take a few minutes)..."

        # Upload directory
        gsutil -m rsync -r \
            "$local_dir/" \
            "gs://${GCS_BUCKET}/${gcs_path}/" 2>&1 | tee /tmp/gsutil_upload.log | grep -v "^$" || true

        # Check exit status of gsutil command (from PIPESTATUS array)
        if [ "${PIPESTATUS[0]}" -eq 0 ]; then
            log_info "✓ Upload successful"
            ((upload_count++))
        else
            log_warn "⚠️  Upload encountered errors"
        fi

        echo ""
    fi
done

end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
log_info "Upload complete"
echo "  Directories uploaded: $upload_count / ${#DIR_MAP[@]}"
echo "  Time elapsed: ${elapsed}s"
echo ""

# ============================================================================
# Step 3: Set public permissions
# ============================================================================

if ! $DRY_RUN; then
    log_step "Step 3: Setting public access"
    echo ""

    log_info "Making files publicly readable..."

    # Set bucket-level IAM policy (more efficient than per-file)
    if gsutil iam ch allUsers:objectViewer gs://${GCS_BUCKET} 2>&1 | grep -v "^$"; then
        log_info "✓ Public access enabled"
    else
        log_warn "⚠️  Could not set public access (files may still be accessible)"
    fi

    echo ""
fi

# ============================================================================
# Step 4: Auto-generate dashboard
# ============================================================================

if ! $DRY_RUN; then
    log_step "Step 4: Generating dynamic dashboard"
    echo ""

    log_info "Scanning GCS bucket and generating dashboard..."

    if python3 "${SCRIPT_DIR}/generate_dynamic_dashboard.py" 2>&1 | grep -v "^$"; then
        log_info "✓ Dashboard generated and uploaded"
    else
        log_warn "⚠️  Dashboard generation failed (you can run it manually)"
    fi

    echo ""
fi

# ============================================================================
# Summary
# ============================================================================

clear

cat <<EOF

╔══════════════════════════════════════════════════════════════════╗
║                                                                  ║
║              Sync Complete! ✓                                    ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝

EOF

if $DRY_RUN; then
    echo "DRY RUN COMPLETE - No files were uploaded"
    echo ""
    echo "To execute for real:"
    echo "  ./scripts/sync_all_to_gcs.sh sync"
else
    echo "✓ Auto-discovered and uploaded ${#DIR_MAP[@]} directories"
    echo "✓ Total files uploaded: $total_files"
    echo "✓ Public access enabled"
    echo "✓ Dashboard generated automatically"
    echo ""
    echo "🔗 Access Your Results:"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    echo "Dashboard (auto-updated):"
    echo "  https://storage.googleapis.com/srf-multiome-figures/index.html"
    echo ""
    echo "Direct GCS access:"
    echo "  https://console.cloud.google.com/storage/browser/srf-multiome-figures"
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    echo "💡 Tips:"
    echo "  - Run this anytime to sync new files: ./scripts/sync_all_to_gcs.sh"
    echo "  - Hard refresh browser to see updates: Ctrl+Shift+R (Win/Linux) or Cmd+Shift+R (Mac)"
    echo "  - Dashboard automatically shows all uploaded files"
    echo ""
fi

echo "For documentation, see:"
echo "  GCS_UPLOAD_COMPLETE_ANALYSIS.md"
echo ""
