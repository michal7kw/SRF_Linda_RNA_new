#!/bin/bash
#
# update_gcs.sh - Simple one-liner to update everything
#
# This is the easiest way to sync all your analysis results to GCS
# and update the dashboard. Just run this anytime!
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "════════════════════════════════════════════════════════════════"
echo "  Quick GCS Update"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "This will:"
echo "  1. Find all directories with PNG/PDF/CSV/BW files"
echo "  2. Upload to GCS (incremental sync)"
echo "  3. Auto-generate dashboard"
echo ""
echo "Press Ctrl+C to cancel, or wait 3 seconds to continue..."
sleep 3

# Run the full sync
"${SCRIPT_DIR}/sync_all_to_gcs.sh" sync
