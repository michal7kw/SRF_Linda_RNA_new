#!/bin/bash
#
# clean_gcs.sh
#
# Remove old/unwanted directories from GCS
#
# Usage:
#   ./clean_gcs.sh [list|delete]
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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

MODE="${1:-list}"

echo "════════════════════════════════════════════════════════════════"
echo "  GCS Cleanup Tool"
echo "════════════════════════════════════════════════════════════════"
echo ""

# Validate GCS configuration
validate_config

# Define directories to potentially remove
UNWANTED_DIRS=(
    "multiome/plots"
)

if [[ "$MODE" == "list" ]]; then
    log_info "Checking for old/unwanted directories in GCS..."
    echo ""

    for dir in "${UNWANTED_DIRS[@]}"; do
        if gsutil ls "gs://${GCS_BUCKET}/${dir}/" &>/dev/null; then
            size=$(gsutil du -sh "gs://${GCS_BUCKET}/${dir}/" 2>/dev/null | cut -f1)
            count=$(gsutil ls -r "gs://${GCS_BUCKET}/${dir}/**" 2>/dev/null | grep -v ":$" | grep -v "/$" | wc -l)

            echo "❌ Found: ${dir}/"
            echo "   Files: $count"
            echo "   Size: $size"
            echo "   URL: gs://${GCS_BUCKET}/${dir}/"
            echo ""
        fi
    done

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    echo "To delete these directories, run:"
    echo "  ./scripts/clean_gcs.sh delete"
    echo ""
    echo "⚠️  Warning: This cannot be undone!"
    echo ""

elif [[ "$MODE" == "delete" ]]; then
    log_warn "This will DELETE the following directories from GCS:"
    echo ""

    for dir in "${UNWANTED_DIRS[@]}"; do
        if gsutil ls "gs://${GCS_BUCKET}/${dir}/" &>/dev/null; then
            echo "  ❌ gs://${GCS_BUCKET}/${dir}/"
        fi
    done

    echo ""
    read -p "Are you sure? Type 'DELETE' to confirm: " confirm

    if [[ "$confirm" == "DELETE" ]]; then
        echo ""
        log_info "Deleting directories..."

        for dir in "${UNWANTED_DIRS[@]}"; do
            if gsutil ls "gs://${GCS_BUCKET}/${dir}/" &>/dev/null; then
                log_info "Deleting ${dir}/..."
                if gsutil -m rm -r "gs://${GCS_BUCKET}/${dir}/" 2>&1 | grep -v "^$"; then
                    log_info "✓ Deleted ${dir}/"
                else
                    log_error "Failed to delete ${dir}/"
                fi
            fi
        done

        echo ""
        log_info "Cleanup complete!"
        echo ""
        log_info "Regenerating dashboard..."
        python3 "${SCRIPT_DIR}/generate_dynamic_dashboard.py"

        echo ""
        log_info "✓ Done! Hard refresh your browser to see changes."
    else
        log_warn "Deletion cancelled"
    fi

else
    log_error "Invalid mode: $MODE"
    echo "Usage: ./scripts/clean_gcs.sh [list|delete]"
    exit 1
fi
