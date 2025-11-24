#!/bin/bash
#
# gcs_config.sh
#
# Configuration file for Google Cloud Storage figure uploads
# This file is sourced by upload_figures_gcs.sh
#
# SETUP INSTRUCTIONS:
# 1. Create a GCS bucket: https://console.cloud.google.com/storage
# 2. Update GCS_BUCKET below with your bucket name
# 3. Update GCS_PROJECT with your Google Cloud project ID
# 4. Run: gcloud init (one-time authentication)
# 5. Test: ./upload_figures_gcs.sh dry-run
#

# ==================== REQUIRED CONFIGURATION ====================

# Your Google Cloud Storage bucket name (REQUIRED - update this!)
# Example: srf-multiome-figures-kubacki
# Get from: https://console.cloud.google.com/storage/browser
GCS_BUCKET="srf-multiome-figures"

# Your Google Cloud Project ID (REQUIRED - update this!)
# Example: srf-multiome-figures-123456
# Get from: gcloud config get-value project
GCS_PROJECT="srf-multiome-figures"

# ==================== PATH CONFIGURATION ====================

# Base directory for multiome pipeline
PIPELINE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"

# Results directory
RESULTS_DIR="${PIPELINE_DIR}/signac_results"

# L1 results directory (broad cell types)
RESULTS_DIR_L1="${PIPELINE_DIR}/signac_results_L1"

# Spatial analysis results (optional)
SPATIAL_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation"

# ==================== UPLOAD CONFIGURATION ====================

# Bucket paths (paths within your GCS bucket)
GCS_PATH_MULTIOME="multiome"
GCS_PATH_SPATIAL="spatial"
GCS_PATH_ATAC="atac"

# Whether to create dated archives (true/false)
ENABLE_ARCHIVES=true

# Number of historical archives to keep (0 = keep all)
MAX_ARCHIVES=10

# File types to upload (space-separated)
UPLOAD_FILETYPES="*.pdf *.png"

# Exclude patterns (space-separated)
EXCLUDE_PATTERNS="*_temp.* *_draft.*"

# ==================== NOTIFICATION CONFIGURATION ====================

# Email notification when uploads complete (leave empty to disable)
NOTIFICATION_EMAIL=""

# Slack webhook URL for notifications (leave empty to disable)
SLACK_WEBHOOK_URL=""

# ==================== ADVANCED CONFIGURATION ====================

# Use parallel uploads for speed (true/false)
# Recommended: true for >100 files, false for <100 files
PARALLEL_UPLOADS=true

# Number of parallel upload threads (if PARALLEL_UPLOADS=true)
PARALLEL_THREADS=8

# Compression for text files (gzip) before upload (true/false)
# Note: PDFs/PNGs are already compressed, so this mainly affects CSVs
COMPRESS_TEXT_FILES=false

# Storage class (STANDARD, NEARLINE, COLDLINE, ARCHIVE)
# STANDARD = frequently accessed (recommended for figures)
# NEARLINE = accessed <1/month, cheaper storage
STORAGE_CLASS="STANDARD"

# Public access control
# Options: publicRead, private, authenticatedRead
PUBLIC_ACCESS="publicRead"

# Cache control header (how long browsers cache files)
# 3600 = 1 hour, 86400 = 1 day, 604800 = 1 week
CACHE_CONTROL="max-age=3600, public"

# ==================== LOGGING CONFIGURATION ====================

# Log file location
LOG_DIR="${PIPELINE_DIR}/logs"
LOG_FILE="${LOG_DIR}/gcs_upload.log"

# Log level (DEBUG, INFO, WARN, ERROR)
LOG_LEVEL="INFO"

# Keep log files for N days
LOG_RETENTION_DAYS=30

# ==================== VALIDATION ====================

# Minimum file size to upload (bytes) - prevents uploading corrupted files
MIN_FILE_SIZE=1024  # 1 KB

# Maximum file size to upload (bytes) - prevents uploading huge temp files
MAX_FILE_SIZE=104857600  # 100 MB

# ==================== DO NOT EDIT BELOW THIS LINE ====================

# Validation function
validate_config() {
    local errors=0

    if [[ "${GCS_BUCKET}" == "YOUR-BUCKET-NAME" ]]; then
        echo "ERROR: Please update GCS_BUCKET in scripts/gcs_config.sh"
        echo "       Create a bucket at: https://console.cloud.google.com/storage"
        errors=$((errors + 1))
    fi

    if [[ "${GCS_PROJECT}" == "YOUR-PROJECT-ID" ]]; then
        echo "ERROR: Please update GCS_PROJECT in scripts/gcs_config.sh"
        echo "       Get your project ID: gcloud config get-value project"
        errors=$((errors + 1))
    fi

    if ! command -v gsutil &> /dev/null; then
        echo "ERROR: gsutil not found. Please install Google Cloud SDK:"
        echo "       https://cloud.google.com/sdk/docs/install"
        errors=$((errors + 1))
    fi

    if ! command -v gcloud &> /dev/null; then
        echo "ERROR: gcloud not found. Please install Google Cloud SDK:"
        echo "       https://cloud.google.com/sdk/docs/install"
        errors=$((errors + 1))
    fi

    if [[ $errors -gt 0 ]]; then
        return 1
    fi

    return 0
}

# Create log directory if it doesn't exist
mkdir -p "${LOG_DIR}"

# Export variables for use in other scripts
export GCS_BUCKET GCS_PROJECT PIPELINE_DIR RESULTS_DIR RESULTS_DIR_L1
export GCS_PATH_MULTIOME GCS_PATH_SPATIAL GCS_PATH_ATAC
export ENABLE_ARCHIVES MAX_ARCHIVES UPLOAD_FILETYPES EXCLUDE_PATTERNS
export NOTIFICATION_EMAIL SLACK_WEBHOOK_URL
export PARALLEL_UPLOADS PARALLEL_THREADS COMPRESS_TEXT_FILES
export STORAGE_CLASS PUBLIC_ACCESS CACHE_CONTROL
export LOG_DIR LOG_FILE LOG_LEVEL LOG_RETENTION_DAYS
export MIN_FILE_SIZE MAX_FILE_SIZE
