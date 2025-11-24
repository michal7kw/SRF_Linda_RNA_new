#!/bin/bash
#
# gcs_quick_setup.sh
#
# Quick setup wizard for Google Cloud Storage figure hosting
# Guides you through the entire setup process interactively
#
# Usage:
#   ./gcs_quick_setup.sh
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/gcs_config.sh"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m'

clear

cat <<'EOF'
╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║        Google Cloud Storage Setup Wizard                     ║
║        SRF Multiome Figure Hosting                           ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝

This wizard will guide you through setting up Google Cloud Storage
for hosting and sharing your analysis figures with stakeholders.

Estimated time: 10-15 minutes

EOF

read -p "Press Enter to continue..."

# Step 1: Check prerequisites
echo ""
echo "═══════════════════════════════════════════════════════════"
echo " Step 1: Checking Prerequisites"
echo "═══════════════════════════════════════════════════════════"
echo ""

# Check gcloud
if command -v gcloud &> /dev/null; then
    echo -e "${GREEN}✓${NC} gcloud SDK is installed"
    gcloud --version | head -1
else
    echo -e "${RED}✗${NC} gcloud SDK is not installed"
    echo ""
    echo "Please install Google Cloud SDK:"
    echo ""
    echo "  cd ~"
    echo "  curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-linux-x86_64.tar.gz"
    echo "  tar -xf google-cloud-cli-linux-x86_64.tar.gz"
    echo "  ./google-cloud-sdk/install.sh"
    echo "  source ~/.bashrc"
    echo ""
    read -p "After installation, press Enter to continue or Ctrl+C to exit..."
    exit 1
fi

# Check gsutil
if command -v gsutil &> /dev/null; then
    echo -e "${GREEN}✓${NC} gsutil is installed"
else
    echo -e "${RED}✗${NC} gsutil is not installed (should come with gcloud SDK)"
    exit 1
fi

echo ""
read -p "Press Enter to continue..."

# Step 2: Authentication
echo ""
echo "═══════════════════════════════════════════════════════════"
echo " Step 2: Google Cloud Authentication"
echo "═══════════════════════════════════════════════════════════"
echo ""

# Check if already authenticated
if gcloud auth list --filter=status:ACTIVE --format="value(account)" &>/dev/null; then
    CURRENT_ACCOUNT=$(gcloud auth list --filter=status:ACTIVE --format="value(account)")
    echo -e "${GREEN}✓${NC} Already authenticated as: ${CURRENT_ACCOUNT}"
    echo ""
    read -p "Do you want to re-authenticate? (y/N): " reauth
    if [[ "$reauth" =~ ^[Yy]$ ]]; then
        gcloud auth login
    fi
else
    echo "You need to authenticate with Google Cloud."
    echo ""
    read -p "Press Enter to open authentication in browser..."
    gcloud auth login
fi

echo ""
read -p "Press Enter to continue..."

# Step 3: Project setup
echo ""
echo "═══════════════════════════════════════════════════════════"
echo " Step 3: Google Cloud Project"
echo "═══════════════════════════════════════════════════════════"
echo ""

# Check current project
CURRENT_PROJECT=$(gcloud config get-value project 2>/dev/null || echo "")

if [[ -n "$CURRENT_PROJECT" ]]; then
    echo -e "${GREEN}✓${NC} Current project: ${CURRENT_PROJECT}"
    echo ""
    read -p "Do you want to use this project? (Y/n): " use_current

    if [[ "$use_current" =~ ^[Nn]$ ]]; then
        CURRENT_PROJECT=""
    fi
fi

if [[ -z "$CURRENT_PROJECT" ]]; then
    echo ""
    echo "Available options:"
    echo "  1. Create a new project"
    echo "  2. Select an existing project"
    echo ""
    read -p "Choose option (1 or 2): " project_option

    if [[ "$project_option" == "1" ]]; then
        echo ""
        read -p "Enter new project name (e.g., srf-multiome-figures): " project_name
        PROJECT_ID="${project_name}-$(date +%s | tail -c 6)"

        echo ""
        echo "Creating project: ${PROJECT_ID}"
        gcloud projects create "${PROJECT_ID}" --name="${project_name}"

        CURRENT_PROJECT="${PROJECT_ID}"
    else
        echo ""
        echo "Available projects:"
        gcloud projects list
        echo ""
        read -p "Enter project ID: " CURRENT_PROJECT
    fi

    gcloud config set project "${CURRENT_PROJECT}"
fi

PROJECT_ID="${CURRENT_PROJECT}"

echo ""
echo -e "${GREEN}✓${NC} Using project: ${PROJECT_ID}"

echo ""
read -p "Press Enter to continue..."

# Step 4: Enable APIs
echo ""
echo "═══════════════════════════════════════════════════════════"
echo " Step 4: Enabling Required APIs"
echo "═══════════════════════════════════════════════════════════"
echo ""

echo "Enabling Cloud Storage API..."
gcloud services enable storage-api.googleapis.com --project="${PROJECT_ID}"

echo ""
echo -e "${GREEN}✓${NC} APIs enabled"

echo ""
read -p "Press Enter to continue..."

# Step 5: Create bucket
echo ""
echo "═══════════════════════════════════════════════════════════"
echo " Step 5: Creating Storage Bucket"
echo "═══════════════════════════════════════════════════════════"
echo ""

# Suggest bucket name
DEFAULT_BUCKET="srf-multiome-figures-${USER}"
echo "Suggested bucket name: ${DEFAULT_BUCKET}"
echo "(Bucket names must be globally unique)"
echo ""
read -p "Enter bucket name [${DEFAULT_BUCKET}]: " BUCKET_NAME
BUCKET_NAME="${BUCKET_NAME:-$DEFAULT_BUCKET}"

# Ask for region
echo ""
echo "Available regions:"
echo "  us-central1     (Iowa, USA)"
echo "  us-east1        (South Carolina, USA)"
echo "  europe-west1    (Belgium)"
echo "  europe-west4    (Netherlands)"
echo "  asia-east1      (Taiwan)"
echo ""
read -p "Enter region [us-central1]: " REGION
REGION="${REGION:-us-central1}"

echo ""
echo "Creating bucket: gs://${BUCKET_NAME}"
echo "Region: ${REGION}"
echo ""

if gsutil mb -p "${PROJECT_ID}" -l "${REGION}" -c STANDARD "gs://${BUCKET_NAME}"; then
    echo -e "${GREEN}✓${NC} Bucket created successfully"
else
    echo -e "${RED}✗${NC} Failed to create bucket (may already exist or name taken)"
    exit 1
fi

# Set public access
echo ""
echo "Setting public access permissions..."
gsutil iam ch allUsers:objectViewer "gs://${BUCKET_NAME}"

echo -e "${GREEN}✓${NC} Public access configured"

echo ""
read -p "Press Enter to continue..."

# Step 6: Update configuration
echo ""
echo "═══════════════════════════════════════════════════════════"
echo " Step 6: Updating Configuration"
echo "═══════════════════════════════════════════════════════════"
echo ""

echo "Updating ${CONFIG_FILE}..."

# Create backup
cp "${CONFIG_FILE}" "${CONFIG_FILE}.backup.$(date +%Y%m%d_%H%M%S)"

# Update config
sed -i "s|GCS_BUCKET=\"YOUR-BUCKET-NAME\"|GCS_BUCKET=\"${BUCKET_NAME}\"|g" "${CONFIG_FILE}"
sed -i "s|GCS_PROJECT=\"YOUR-PROJECT-ID\"|GCS_PROJECT=\"${PROJECT_ID}\"|g" "${CONFIG_FILE}"

echo -e "${GREEN}✓${NC} Configuration updated"

echo ""
read -p "Press Enter to continue..."

# Step 7: Test setup
echo ""
echo "═══════════════════════════════════════════════════════════"
echo " Step 7: Testing Setup"
echo "═══════════════════════════════════════════════════════════"
echo ""

echo "Running tests..."
"${SCRIPT_DIR}/upload_figures_gcs.sh" test

echo ""
read -p "Press Enter to continue..."

# Step 8: Integration
echo ""
echo "═══════════════════════════════════════════════════════════"
echo " Step 8: Pipeline Integration"
echo "═══════════════════════════════════════════════════════════"
echo ""

echo "Do you want to automatically upload figures after pipeline completion?"
echo ""
read -p "Enable automatic uploads? (Y/n): " enable_auto

if [[ ! "$enable_auto" =~ ^[Nn]$ ]]; then
    echo ""
    echo "Integrating GCS upload into pipeline scripts..."
    "${SCRIPT_DIR}/integrate_gcs_upload.sh" install
else
    echo ""
    echo "Automatic upload disabled. You can enable it later with:"
    echo "  ./scripts/integrate_gcs_upload.sh install"
fi

echo ""
read -p "Press Enter to continue..."

# Step 9: First upload (optional)
echo ""
echo "═══════════════════════════════════════════════════════════"
echo " Step 9: Initial Upload"
echo "═══════════════════════════════════════════════════════════"
echo ""

echo "Do you want to upload existing figures now?"
echo ""
read -p "Upload figures? (y/N): " do_upload

if [[ "$do_upload" =~ ^[Yy]$ ]]; then
    echo ""
    "${SCRIPT_DIR}/upload_figures_gcs.sh" upload
fi

# Final summary
clear

cat <<EOF

╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║                   🎉 Setup Complete! 🎉                       ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝

Configuration Summary:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  GCS Bucket:  ${BUCKET_NAME}
  Project ID:  ${PROJECT_ID}
  Region:      ${REGION}

Public URLs:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Dashboard:
    https://storage.googleapis.com/${BUCKET_NAME}/index.html

  Latest figures:
    https://storage.googleapis.com/${BUCKET_NAME}/multiome/plots/latest/

Next Steps:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  1. Share the dashboard URL with your stakeholders
  2. Run your analysis pipeline: sbatch run_signac_pipeline.sh
  3. Figures will automatically upload to GCS

Manual Commands:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Upload figures:          ./scripts/upload_figures_gcs.sh upload
  Quick sync:              ./scripts/upload_figures_gcs.sh quick
  View uploaded files:     ./scripts/upload_figures_gcs.sh list
  Check storage usage:     ./scripts/upload_figures_gcs.sh stats
  Get public URL:          ./scripts/upload_figures_gcs.sh get-url

Documentation:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  See GCS_SETUP_GUIDE.md for detailed documentation

Estimated Costs:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Storage (10GB):          ~$0.20/month
  Operations:              ~$0.00/month (free tier)
  Network egress:          ~$0.20/month
  Total:                   ~$0.40/month

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

EOF

echo ""
