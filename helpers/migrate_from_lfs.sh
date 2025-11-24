#!/bin/bash
#
# migrate_from_lfs.sh
#
# Migrate from Git LFS to GCS for figure hosting
# Run from: SRF_Linda_RNA/ (the Git repository root)
#
# Usage:
#   ./migrate_from_lfs.sh [dry-run|execute]
#

set -euo pipefail

# Verify we're in a Git repository
if [[ ! -d .git ]]; then
    echo "ERROR: Not in a Git repository root directory"
    echo "Current directory: $(pwd)"
    echo "Please cd to SRF_Linda_RNA/ and try again"
    exit 1
fi

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

DRY_RUN=false
[[ "${1:-execute}" == "dry-run" ]] && DRY_RUN=true

log_info() { echo -e "${GREEN}[INFO]${NC} $*"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }
log_step() { echo -e "${BLUE}[STEP]${NC} $*"; }

clear

cat <<'EOF'
╔══════════════════════════════════════════════════════════════════╗
║                                                                  ║
║          Git LFS to GCS Migration for SRF_Linda_RNA             ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝

Current LFS Usage:
  - 5,553 total LFS files
  - 2,128 PNGs (figures → moving to GCS)
  - 2,589 CSVs (data files)
  - 758 Excel files (data files)
  - 78 text files

This script will:
  ✓ Stop LFS tracking for PDFs and PNGs (figures)
  ✓ Keep LFS tracking for CSVs, Excel, TXT (data files)
  ✓ Update .gitignore to exclude figures
  ✓ Commit changes
  ✓ Clean up local LFS cache

EOF

if $DRY_RUN; then
    log_warn "DRY RUN MODE - No changes will be made"
    echo ""
fi

read -p "Press Enter to continue or Ctrl+C to cancel..."
echo ""

# ============================================================================
# Step 1: Show Current LFS Status
# ============================================================================

log_step "Step 1: Current LFS Configuration"
echo ""

log_info "Files tracked by LFS in .gitattributes:"
cat .gitattributes
echo ""

log_info "LFS files in repository (by type):"
git lfs ls-files | grep -o '\.[^.]*$' | sort | uniq -c | sort -rn
echo ""

log_info "Local LFS cache size:"
if [[ -d .git/lfs ]]; then
    du -sh .git/lfs
else
    echo "  (no LFS cache)"
fi

echo ""
read -p "Press Enter to continue..."

# ============================================================================
# Step 2: Untrack Figures (PDFs and PNGs)
# ============================================================================

log_step "Step 2: Stopping LFS tracking for figures (PDFs, PNGs)"
echo ""

if $DRY_RUN; then
    log_warn "[DRY RUN] Would run:"
    echo "  git lfs untrack '*.pdf'"
    echo "  git lfs untrack '*.png'"
else
    log_info "Untracking PDFs from LFS..."
    git lfs untrack "*.pdf" || log_warn "PDF untracking returned error (may not be tracked)"

    log_info "Untracking PNGs from LFS..."
    git lfs untrack "*.png" || log_warn "PNG untracking returned error (may not be tracked)"

    log_info "✓ Figures (PDF, PNG) no longer tracked by LFS"
    echo ""
    log_info "Updated .gitattributes:"
    cat .gitattributes
fi

echo ""
log_info "Data files (CSV, Excel, TXT) remain in LFS"
log_info "If you want to untrack them too, run manually:"
log_info "  git lfs untrack '*.csv'"
log_info "  git lfs untrack '*.xlsx'"
log_info "  git lfs untrack '*.txt'"
echo ""

read -p "Press Enter to continue..."

# ============================================================================
# Step 3: Verify .gitignore
# ============================================================================

log_step "Step 3: Verifying .gitignore excludes figures"
echo ""

if grep -q "^\*.pdf" .gitignore && grep -q "^\*.png" .gitignore; then
    log_info "✓ .gitignore already excludes PDFs and PNGs"
else
    log_warn ".gitignore doesn't properly exclude figures"
    log_info "This should have been added during GCS setup"

    if $DRY_RUN; then
        log_warn "[DRY RUN] Would check/update .gitignore"
    else
        read -p "Should I add figure exclusions to .gitignore? (y/N): " add_gitignore
        if [[ "$add_gitignore" =~ ^[Yy]$ ]]; then
            if ! grep -q "Figure Hosting - Google Cloud Storage" .gitignore; then
                cat >> .gitignore <<'GITIGNORE'

################################################################################
# Figure Hosting - Google Cloud Storage
################################################################################
# All figures are hosted on GCS, not in Git

# PDF and PNG figures (hosted on GCS)
*.pdf
*.png
!docs/workflow_diagrams/*.png
!README_figures/*.png
GITIGNORE
                log_info "✓ Added figure exclusions to .gitignore"
            fi
        fi
    fi
fi

echo ""
read -p "Press Enter to continue..."

# ============================================================================
# Step 4: Commit Changes
# ============================================================================

log_step "Step 4: Committing changes"
echo ""

if git diff --quiet .gitattributes .gitignore 2>/dev/null; then
    log_info "No changes to commit"
else
    if $DRY_RUN; then
        log_warn "[DRY RUN] Would commit:"
        echo ""
        git diff .gitattributes .gitignore
    else
        log_info "Staging changes..."
        git add .gitattributes .gitignore

        log_info "Committing..."
        git commit -m "Migrate figure hosting from Git LFS to Google Cloud Storage

Stop Git LFS tracking for PDFs and PNGs - now using GCS:
- Dashboard: https://storage.googleapis.com/srf-multiome-figures/index.html
- Cost: ~\$0.40/month (vs \$5+/month for Git LFS)

Data files (CSV, Excel, TXT) remain in Git LFS.
Figures (PDF, PNG) are now in .gitignore and hosted on GCS.

See: ../MIGRATE_FROM_GIT_LFS.md"

        log_info "✓ Changes committed"

        echo ""
        read -p "Push to remote? (y/N): " push_confirm
        if [[ "$push_confirm" =~ ^[Yy]$ ]]; then
            log_info "Pushing..."
            git push || {
                log_error "Push failed. You may need to push manually"
                log_info "Run: git push"
            }
            log_info "✓ Pushed to remote"
        else
            log_warn "Changes committed locally but NOT pushed"
            log_warn "Remember to push: git push"
        fi
    fi
fi

echo ""
read -p "Press Enter to continue..."

# ============================================================================
# Step 5: Clean Up LFS Cache
# ============================================================================

log_step "Step 5: Cleaning up local LFS cache"
echo ""

if [[ ! -d .git/lfs ]]; then
    log_info "No LFS cache to clean"
else
    cache_before=$(du -sh .git/lfs 2>/dev/null | cut -f1)
    log_info "LFS cache size before: ${cache_before}"

    if $DRY_RUN; then
        log_warn "[DRY RUN] Would run: git lfs prune && git gc --prune=now"
    else
        log_info "Running git lfs prune..."
        git lfs prune 2>&1 || log_warn "git lfs prune failed (not critical)"

        log_info "Running git gc..."
        git gc --prune=now 2>&1 || log_warn "git gc failed (not critical)"

        cache_after=$(du -sh .git/lfs 2>/dev/null | cut -f1)
        log_info "LFS cache size after: ${cache_after}"
        log_info "✓ Cache cleaned"
    fi
fi

echo ""
read -p "Press Enter to continue..."

# ============================================================================
# Step 6: Verify GCS Works
# ============================================================================

log_step "Step 6: Verifying GCS is configured"
echo ""

GCS_SCRIPT="integration_scripts/multiome_modular_pipeline/scripts/upload_figures_gcs.sh"

if [[ ! -f "$GCS_SCRIPT" ]]; then
    log_error "GCS script not found: $GCS_SCRIPT"
else
    if $DRY_RUN; then
        log_warn "[DRY RUN] Would test GCS upload"
    else
        log_info "Testing GCS configuration..."
        "$GCS_SCRIPT" test || {
            log_error "GCS test failed"
            log_info "Run manually: $GCS_SCRIPT test"
        }
    fi
fi

echo ""

# ============================================================================
# Summary
# ============================================================================

clear

cat <<EOF

╔══════════════════════════════════════════════════════════════════╗
║                                                                  ║
║                 Migration Complete! ✓                            ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝

EOF

if $DRY_RUN; then
    echo "DRY RUN COMPLETE - No changes were made"
    echo ""
    echo "To execute for real:"
    echo "  ./migrate_from_lfs.sh execute"
else
    echo "✓ Figures (PDF, PNG) no longer tracked by LFS"
    echo "✓ Data files (CSV, Excel, TXT) still in LFS"
    echo "✓ .gitignore excludes figures"
    echo "✓ Changes committed $(git diff --quiet HEAD || echo '(not pushed)')"
    echo "✓ Local LFS cache cleaned"
fi

echo ""
echo "What Changed:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  BEFORE: New PDFs/PNGs → Git LFS → \$\$\$ costs"
echo "  AFTER:  New PDFs/PNGs → Ignored → GCS (\$0.40/month)"
echo ""
echo "  Data files (CSV, Excel) still in Git LFS (unchanged)"
echo ""

echo "Next Steps:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "1. Test new workflow:"
echo "   cd integration_scripts/multiome_modular_pipeline"
echo "   sbatch run_signac_pipeline.sh"
echo "   # PDFs/PNGs should upload to GCS, not Git"
echo ""
echo "2. Verify figures accessible:"
echo "   https://storage.googleapis.com/srf-multiome-figures/index.html"
echo ""
echo "3. (Optional) Migrate data files from LFS:"
echo "   # CSVs, Excel, TXT can stay in LFS or move to .gitignore"
echo "   # See: ../MIGRATE_FROM_GIT_LFS.md"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
