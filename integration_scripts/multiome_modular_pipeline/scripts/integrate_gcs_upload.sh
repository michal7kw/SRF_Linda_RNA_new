#!/bin/bash
#
# integrate_gcs_upload.sh
#
# Integrates GCS upload functionality into existing pipeline scripts
# This script modifies run_signac_pipeline.sh and run_signac_pipeline_L1.sh
# to automatically upload figures after successful completion
#
# Usage:
#   ./integrate_gcs_upload.sh [install|uninstall|status]
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "${SCRIPT_DIR}")"

# Color output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[INFO]${NC} $*"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }

# GCS upload code to be added to pipeline scripts
GCS_UPLOAD_CODE='
################################################################################
# Automatic Figure Upload to Google Cloud Storage
################################################################################

echo "================================================================================"
echo "STEP 6: Uploading Figures to Google Cloud Storage"
echo "================================================================================"
echo ""

# Check if GCS upload is configured
GCS_UPLOAD_SCRIPT="${BASE_DIR}/scripts/upload_figures_gcs.sh"
GCS_CONFIG_SCRIPT="${BASE_DIR}/scripts/gcs_config.sh"

if [[ -f "${GCS_UPLOAD_SCRIPT}" ]] && [[ -f "${GCS_CONFIG_SCRIPT}" ]]; then
    # Source config to check if configured
    source "${GCS_CONFIG_SCRIPT}"

    if [[ "${GCS_BUCKET}" != "YOUR-BUCKET-NAME" ]] && [[ "${GCS_PROJECT}" != "YOUR-PROJECT-ID" ]]; then
        log_info "GCS upload is configured, uploading figures..."

        START_TIME=$(date +%s)

        if "${GCS_UPLOAD_SCRIPT}" quick 2>&1; then
            END_TIME=$(date +%s)
            ELAPSED=$((END_TIME - START_TIME))
            echo ""
            echo "✓ Figures uploaded to GCS (elapsed time: ${ELAPSED}s)"
            echo ""

            # Get and display public URL
            "${GCS_UPLOAD_SCRIPT}" get-url
        else
            echo ""
            echo "⚠ Figure upload failed (pipeline results are still available locally)"
            echo "  You can manually upload later with: ${GCS_UPLOAD_SCRIPT} upload"
            echo ""
        fi
    else
        echo "ℹ GCS upload not configured (skipping)"
        echo "  To enable automatic uploads, edit: ${GCS_CONFIG_SCRIPT}"
        echo "  See: ${BASE_DIR}/GCS_SETUP_GUIDE.md"
        echo ""
    fi
else
    echo "ℹ GCS upload not available (scripts not found)"
    echo ""
fi

echo "================================================================================"
'

# Backup a file
backup_file() {
    local file="$1"
    local backup="${file}.backup.$(date +%Y%m%d_%H%M%S)"

    if [[ -f "$file" ]]; then
        cp "$file" "$backup"
        log_info "Created backup: $(basename "$backup")"
        return 0
    fi

    return 1
}

# Check if GCS upload is already integrated
is_integrated() {
    local file="$1"

    if [[ -f "$file" ]]; then
        if grep -q "Uploading Figures to Google Cloud Storage" "$file"; then
            return 0  # Already integrated
        fi
    fi

    return 1  # Not integrated
}

# Install GCS upload integration
install_integration() {
    log_info "Installing GCS upload integration..."
    echo ""

    local modified_count=0

    # List of pipeline scripts to modify
    local scripts=(
        "${PIPELINE_DIR}/run_signac_pipeline.sh"
        "${PIPELINE_DIR}/run_signac_pipeline_L1.sh"
    )

    for script in "${scripts[@]}"; do
        if [[ ! -f "$script" ]]; then
            log_warn "Script not found: $(basename "$script")"
            continue
        fi

        if is_integrated "$script"; then
            log_info "Already integrated: $(basename "$script")"
            continue
        fi

        log_info "Integrating: $(basename "$script")"

        # Create backup
        backup_file "$script"

        # Find insertion point (before final "echo" statements)
        # Insert before the last echo "================" line
        local temp_file=$(mktemp)

        # Use awk to insert code before the final summary
        awk -v code="$GCS_UPLOAD_CODE" '
            /^echo "================================================================================"$/ && !inserted && NR > 200 {
                print code
                inserted = 1
            }
            { print }
        ' "$script" > "$temp_file"

        # Replace original file
        mv "$temp_file" "$script"
        chmod +x "$script"

        log_info "✓ Integrated: $(basename "$script")"
        ((modified_count++))
    done

    echo ""
    if [[ $modified_count -gt 0 ]]; then
        log_info "Integration complete! Modified ${modified_count} script(s)."
        echo ""
        log_info "Next steps:"
        echo "  1. Configure GCS: nano scripts/gcs_config.sh"
        echo "  2. Test setup: ./scripts/upload_figures_gcs.sh test"
        echo "  3. Run pipeline: sbatch run_signac_pipeline.sh"
        echo ""
        log_info "Figures will now automatically upload after pipeline completion!"
    else
        log_warn "No scripts were modified (already integrated or not found)"
    fi
}

# Uninstall GCS upload integration
uninstall_integration() {
    log_info "Uninstalling GCS upload integration..."
    echo ""

    local modified_count=0

    local scripts=(
        "${PIPELINE_DIR}/run_signac_pipeline.sh"
        "${PIPELINE_DIR}/run_signac_pipeline_L1.sh"
    )

    for script in "${scripts[@]}"; do
        if [[ ! -f "$script" ]]; then
            log_warn "Script not found: $(basename "$script")"
            continue
        fi

        if ! is_integrated "$script"; then
            log_info "Not integrated: $(basename "$script")"
            continue
        fi

        # Check if backup exists
        local latest_backup=$(ls -t "${script}.backup."* 2>/dev/null | head -1)

        if [[ -z "$latest_backup" ]]; then
            log_error "No backup found for: $(basename "$script")"
            log_error "Cannot safely uninstall. Please manually remove GCS upload section."
            continue
        fi

        log_info "Restoring from backup: $(basename "$latest_backup")"
        cp "$latest_backup" "$script"
        chmod +x "$script"

        log_info "✓ Restored: $(basename "$script")"
        ((modified_count++))
    done

    echo ""
    if [[ $modified_count -gt 0 ]]; then
        log_info "Uninstall complete! Restored ${modified_count} script(s) from backup."
    else
        log_warn "No scripts were restored"
    fi
}

# Show integration status
show_status() {
    log_info "GCS Upload Integration Status"
    echo ""

    local scripts=(
        "${PIPELINE_DIR}/run_signac_pipeline.sh"
        "${PIPELINE_DIR}/run_signac_pipeline_L1.sh"
    )

    for script in "${scripts[@]}"; do
        local basename=$(basename "$script")

        if [[ ! -f "$script" ]]; then
            echo "  $basename: NOT FOUND"
            continue
        fi

        if is_integrated "$script"; then
            echo "  $basename: ✓ INTEGRATED"
        else
            echo "  $basename: ✗ NOT INTEGRATED"
        fi

        # Check for backups
        local backup_count=$(ls -1 "${script}.backup."* 2>/dev/null | wc -l)
        if [[ $backup_count -gt 0 ]]; then
            echo "    Backups: ${backup_count} available"
        fi
    done

    echo ""

    # Check GCS configuration
    if [[ -f "${SCRIPT_DIR}/gcs_config.sh" ]]; then
        source "${SCRIPT_DIR}/gcs_config.sh"

        echo "GCS Configuration:"
        if [[ "${GCS_BUCKET}" == "YOUR-BUCKET-NAME" ]]; then
            echo "  Bucket: ✗ NOT CONFIGURED"
        else
            echo "  Bucket: ✓ ${GCS_BUCKET}"
        fi

        if [[ "${GCS_PROJECT}" == "YOUR-PROJECT-ID" ]]; then
            echo "  Project: ✗ NOT CONFIGURED"
        else
            echo "  Project: ✓ ${GCS_PROJECT}"
        fi
    else
        echo "GCS Configuration: NOT FOUND"
    fi

    echo ""
}

# Show help
show_help() {
    cat <<EOF
GCS Upload Integration Tool

Usage: ./integrate_gcs_upload.sh [command]

Commands:
  install     Add GCS upload to pipeline scripts (creates backups)
  uninstall   Remove GCS upload from pipeline scripts (restores from backup)
  status      Show integration and configuration status
  help        Show this help message

Description:
  This tool automatically integrates Google Cloud Storage upload functionality
  into your Signac pipeline scripts. After integration, figures will be
  automatically uploaded to GCS when the pipeline completes successfully.

Examples:
  # Check current status
  ./integrate_gcs_upload.sh status

  # Install integration
  ./integrate_gcs_upload.sh install

  # Remove integration
  ./integrate_gcs_upload.sh uninstall

Notes:
  - Backups are created automatically before modification
  - You must still configure GCS (edit scripts/gcs_config.sh)
  - See GCS_SETUP_GUIDE.md for complete setup instructions

EOF
}

# Main execution
main() {
    local command="${1:-status}"

    echo ""
    log_info "GCS Upload Integration Tool"
    echo ""

    case "$command" in
        install)
            install_integration
            ;;
        uninstall)
            uninstall_integration
            ;;
        status)
            show_status
            ;;
        help|--help|-h)
            show_help
            ;;
        *)
            log_error "Unknown command: $command"
            echo ""
            show_help
            exit 1
            ;;
    esac
}

main "$@"
