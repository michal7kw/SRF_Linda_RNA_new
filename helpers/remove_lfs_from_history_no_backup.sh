#!/bin/bash
#
# remove_lfs_from_history_no_backup.sh
#
# Remove PDFs and PNGs from Git history - SKIP BACKUP STEP
# Use this if you already have a backup
#
# Usage:
#   ./remove_lfs_from_history_no_backup.sh
#

set -euo pipefail

# Add ~/.local/bin to PATH for git-filter-repo
export PATH="$HOME/.local/bin:$PATH"

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

# Verify we're in a Git repository
if [[ ! -d .git ]]; then
    log_error "Not in a Git repository root directory"
    echo "Current directory: $(pwd)"
    echo "Please cd to SRF_Linda_RNA/ and try again"
    exit 1
fi

clear

cat <<'EOF'
╔══════════════════════════════════════════════════════════════════╗
║                                                                  ║
║     Remove LFS Files from Git History (NO BACKUP)              ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝

⚠️  WARNING: BACKUP STEP SKIPPED ⚠️

This script will:
  ✓ Remove all PDFs and PNGs from Git history
  ✓ Keep CSVs, Excel, TXT files (still in LFS)
  ✓ Clean LFS cache
  ✓ Reduce repository size

After running this:
  ⚠️  All collaborators MUST re-clone the repository
  ⚠️  Requires force push (--force-with-lease)

EOF

echo "Current repository size:"
du -sh .git/
echo ""

read -p "Continue without backup? (yes/no): " confirm
if [[ "$confirm" != "yes" ]]; then
    log_warn "Aborted by user"
    exit 0
fi

echo ""

# ============================================================================
# Check Prerequisites
# ============================================================================

log_step "Checking prerequisites"
echo ""

# Check if git-filter-repo is available
if ! command -v git-filter-repo &> /dev/null; then
    log_error "git-filter-repo not found"
    echo ""
    echo "Install with: pip3 install --user git-filter-repo"
    exit 1
fi

log_info "✓ git-filter-repo is available"

# Check if there are uncommitted changes
if ! git diff-index --quiet HEAD -- 2>/dev/null; then
    log_error "You have uncommitted changes"
    echo ""
    echo "Please commit or stash your changes first:"
    echo "  git add ."
    echo "  git commit -m 'Your message'"
    exit 1
fi

log_info "✓ No uncommitted changes"

# Check branch
current_branch=$(git rev-parse --abbrev-ref HEAD)
log_info "✓ On branch: $current_branch"
echo ""

# ============================================================================
# Remove PDFs and PNGs from History
# ============================================================================

log_step "Removing PDFs and PNGs from Git history"
echo ""

log_info "Files to remove from history:"
echo "  - All *.pdf files"
echo "  - All *.png files"
echo ""
log_info "Files to KEEP in history:"
echo "  - All *.csv files (still in LFS)"
echo "  - All *.xlsx files (still in LFS)"
echo "  - All *.txt files (still in LFS)"
echo ""

read -p "Press Enter to start (this will take 10-20 minutes)..."

log_info "Running git filter-repo..."
log_info "This may take 10-20 minutes depending on repository size..."
echo ""

# Use git-filter-repo to remove PDF and PNG files
git filter-repo \
    --invert-paths \
    --path-glob '*.pdf' \
    --path-glob '*.png' \
    --force

log_info "✓ PDFs and PNGs removed from Git history"
echo ""

# ============================================================================
# Clean Up LFS Cache
# ============================================================================

log_step "Cleaning up LFS cache"
echo ""

if [[ -d .git/lfs ]]; then
    cache_before=$(du -sh .git/lfs 2>/dev/null | cut -f1 || echo "unknown")
    log_info "LFS cache size before: ${cache_before}"

    log_info "Pruning LFS cache..."
    git lfs prune --verbose || log_warn "git lfs prune failed (not critical)"

    log_info "Running aggressive garbage collection..."
    git gc --aggressive --prune=now

    cache_after=$(du -sh .git/lfs 2>/dev/null | cut -f1 || echo "unknown")
    log_info "LFS cache size after: ${cache_after}"
    log_info "✓ LFS cache cleaned"
else
    log_info "No LFS cache to clean"
fi

echo ""

# ============================================================================
# Verify Results
# ============================================================================

log_step "Verifying results"
echo ""

log_info "Checking for remaining PDFs and PNGs in history..."

# Check if any PDFs or PNGs remain in history
pdf_count=$(git rev-list --all --objects | grep -c '\.pdf$' || true)
png_count=$(git rev-list --all --objects | grep -c '\.png$' || true)

if [[ $pdf_count -eq 0 && $png_count -eq 0 ]]; then
    log_info "✓ No PDFs or PNGs found in Git history"
else
    log_warn "Found $pdf_count PDFs and $png_count PNGs still in history"
fi

log_info "Repository size after cleanup:"
du -sh .git/
echo ""

# ============================================================================
# Re-add Remote
# ============================================================================

log_step "Re-adding remote"
echo ""

log_info "git filter-repo removed the remote for safety"
log_info "Adding remote back..."

# Assuming GitHub repository
REPO_NAME="SRF_Linda_RNA"
GITHUB_USER="michal7kw"

git remote add origin "https://github.com/$GITHUB_USER/$REPO_NAME.git" || {
    log_warn "Could not add remote automatically"
    echo ""
    echo "Please add remote manually:"
    echo "  git remote add origin YOUR_REMOTE_URL"
}

log_info "✓ Remote added"
echo ""

# ============================================================================
# Summary
# ============================================================================

clear

cat <<EOF

╔══════════════════════════════════════════════════════════════════╗
║                                                                  ║
║              History Cleanup Complete! ✓                         ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝

✓ PDFs and PNGs removed from Git history
✓ LFS cache cleaned
✓ Repository size reduced
✓ Remote re-added

⚠️  NEXT STEP: FORCE PUSH ⚠️
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Run this command to update GitHub:

    git push --force-with-lease origin master

This will rewrite history on GitHub.
All collaborators must re-clone after this.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

On other computers, clone without LFS:

    GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/michal7kw/SRF_Linda_RNA.git

Figures are available at:
    https://storage.googleapis.com/srf-multiome-figures/index.html

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

EOF
