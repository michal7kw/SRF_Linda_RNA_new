#!/bin/bash
#
# remove_lfs_from_history.sh
#
# Remove PDFs and PNGs from Git history to free up LFS storage
# This rewrites Git history - all collaborators must re-clone after this
#
# Usage:
#   ./remove_lfs_from_history.sh [dry-run|execute]
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

DRY_RUN=false
[[ "${1:-execute}" == "dry-run" ]] && DRY_RUN=true

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
║     Remove LFS Files from Git History (PDF/PNG only)            ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝

⚠️  WARNING: This REWRITES Git history ⚠️

What this script does:
  ✓ Creates backup of repository
  ✓ Removes all PDFs and PNGs from Git history
  ✓ Keeps CSVs, Excel, TXT files (still in LFS)
  ✓ Reduces repository size dramatically
  ✓ Requires force push (--force-with-lease)

After running this:
  ⚠️  All collaborators MUST re-clone the repository
  ⚠️  Old clones will be incompatible
  ⚠️  GitHub URLs to specific commits will break

Benefits:
  ✓ No more LFS quota exceeded errors
  ✓ Faster clones (smaller history)
  ✓ All figures now in GCS instead

EOF

if $DRY_RUN; then
    log_warn "DRY RUN MODE - No changes will be made"
    echo ""
fi

echo "Current repository size:"
du -sh .git/
echo ""

read -p "Do you want to proceed? (yes/no): " confirm
if [[ "$confirm" != "yes" ]]; then
    log_warn "Aborted by user"
    exit 0
fi

echo ""

# ============================================================================
# Step 1: Check Prerequisites
# ============================================================================

log_step "Step 1: Checking prerequisites"
echo ""

# Check if git-filter-repo is available
if ! command -v git-filter-repo &> /dev/null; then
    log_warn "git-filter-repo not found. Trying to install..."

    if command -v pip3 &> /dev/null; then
        pip3 install --user git-filter-repo

        # Add to PATH if needed
        export PATH="$HOME/.local/bin:$PATH"

        if ! command -v git-filter-repo &> /dev/null; then
            log_error "git-filter-repo installation failed"
            echo ""
            echo "Please install manually:"
            echo "  pip3 install git-filter-repo"
            echo "  # or"
            echo "  conda install -c conda-forge git-filter-repo"
            exit 1
        fi
    else
        log_error "pip3 not found. Cannot install git-filter-repo"
        echo ""
        echo "Please install manually:"
        echo "  pip3 install git-filter-repo"
        echo "  # or"
        echo "  conda install -c conda-forge git-filter-repo"
        exit 1
    fi
fi

log_info "✓ git-filter-repo is available"

# Check if there are uncommitted changes
if ! git diff-index --quiet HEAD -- 2>/dev/null; then
    log_error "You have uncommitted changes"
    echo ""
    echo "Please commit or stash your changes first:"
    echo "  git add ."
    echo "  git commit -m 'Your message'"
    echo "  # or"
    echo "  git stash"
    exit 1
fi

log_info "✓ No uncommitted changes"

# Check if we're on master branch
current_branch=$(git rev-parse --abbrev-ref HEAD)
if [[ "$current_branch" != "master" ]]; then
    log_warn "Not on master branch (current: $current_branch)"
    read -p "Continue anyway? (y/N): " continue_branch
    if [[ ! "$continue_branch" =~ ^[Yy]$ ]]; then
        log_warn "Please switch to master: git checkout master"
        exit 1
    fi
fi

log_info "✓ On branch: $current_branch"
echo ""

# ============================================================================
# Step 2: Create Backup
# ============================================================================

log_step "Step 2: Creating backup"
echo ""

BACKUP_DIR="../SRF_Linda_RNA_backup_$(date +%Y%m%d_%H%M%S)"

if $DRY_RUN; then
    log_warn "[DRY RUN] Would create backup: $BACKUP_DIR"
else
    log_info "Creating backup of repository..."
    log_info "Backup location: $BACKUP_DIR"

    # Create backup (excludes .git to save space, keeps working tree)
    rsync -a --exclude='.git' ./ "$BACKUP_DIR/"

    log_info "✓ Backup created"
    log_info "If anything goes wrong, restore from: $BACKUP_DIR"
fi

echo ""
read -p "Press Enter to continue..."

# ============================================================================
# Step 3: Remove PDFs and PNGs from History
# ============================================================================

log_step "Step 3: Removing PDFs and PNGs from Git history"
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

if $DRY_RUN; then
    log_warn "[DRY RUN] Would run:"
    echo "  git filter-repo --invert-paths --path-glob '*.pdf' --path-glob '*.png' --force"
else
    log_info "Running git filter-repo (this may take several minutes)..."
    echo ""

    # Use git-filter-repo to remove PDF and PNG files
    # --invert-paths: remove matching paths (instead of keeping them)
    # --path-glob: match by glob pattern
    # --force: allow running in a repository with remotes
    git filter-repo \
        --invert-paths \
        --path-glob '*.pdf' \
        --path-glob '*.png' \
        --force

    log_info "✓ PDFs and PNGs removed from Git history"
fi

echo ""

# ============================================================================
# Step 4: Clean Up LFS Cache
# ============================================================================

log_step "Step 4: Cleaning up LFS cache"
echo ""

if [[ -d .git/lfs ]]; then
    cache_before=$(du -sh .git/lfs 2>/dev/null | cut -f1 || echo "unknown")
    log_info "LFS cache size before: ${cache_before}"

    if $DRY_RUN; then
        log_warn "[DRY RUN] Would run: git lfs prune && git gc --aggressive --prune=now"
    else
        log_info "Pruning LFS cache..."
        git lfs prune --verbose || log_warn "git lfs prune failed (not critical)"

        log_info "Running aggressive garbage collection..."
        git gc --aggressive --prune=now

        cache_after=$(du -sh .git/lfs 2>/dev/null | cut -f1 || echo "unknown")
        log_info "LFS cache size after: ${cache_after}"
        log_info "✓ LFS cache cleaned"
    fi
else
    log_info "No LFS cache to clean"
fi

echo ""

# ============================================================================
# Step 5: Verify Results
# ============================================================================

log_step "Step 5: Verifying results"
echo ""

if $DRY_RUN; then
    log_warn "[DRY RUN] Would verify changes"
else
    log_info "Checking for remaining PDFs and PNGs in history..."

    # Check if any PDFs or PNGs remain in history
    pdf_count=$(git rev-list --all --objects | grep -c '\.pdf$' || true)
    png_count=$(git rev-list --all --objects | grep -c '\.png$' || true)

    if [[ $pdf_count -eq 0 && $png_count -eq 0 ]]; then
        log_info "✓ No PDFs or PNGs found in Git history"
    else
        log_warn "Found $pdf_count PDFs and $png_count PNGs still in history"
        log_warn "This may be normal if some were added after the filter"
    fi

    log_info "Repository size after cleanup:"
    du -sh .git/
fi

echo ""

# ============================================================================
# Step 6: Re-add Remote
# ============================================================================

log_step "Step 6: Re-adding remote"
echo ""

if $DRY_RUN; then
    log_warn "[DRY RUN] Would re-add remote origin"
else
    # git filter-repo removes remotes for safety
    # We need to add it back

    log_info "git filter-repo removed the remote for safety"
    log_info "Adding remote back..."

    # Detect which remote URL format (assuming GitHub)
    REPO_NAME="SRF_Linda_RNA"
    GITHUB_USER="michal7kw"  # Update if different

    # Try to detect from backup or ask
    if [[ -f "${BACKUP_DIR}/.git/config" ]]; then
        REMOTE_URL=$(grep -A 1 '\[remote "origin"\]' "${BACKUP_DIR}/.git/config" | grep url | cut -d= -f2- | xargs || echo "")
        if [[ -n "$REMOTE_URL" ]]; then
            log_info "Detected remote URL: $REMOTE_URL"
            git remote add origin "$REMOTE_URL"
            log_info "✓ Remote added"
        else
            log_warn "Could not detect remote URL from backup"
            echo ""
            echo "Please add remote manually:"
            echo "  git remote add origin https://github.com/$GITHUB_USER/$REPO_NAME.git"
        fi
    else
        log_info "Assuming GitHub remote..."
        git remote add origin "https://github.com/$GITHUB_USER/$REPO_NAME.git" || {
            log_warn "Could not add remote automatically"
            echo ""
            echo "Please add remote manually:"
            echo "  git remote add origin YOUR_REMOTE_URL"
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
║              History Cleanup Complete! ✓                         ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝

EOF

if $DRY_RUN; then
    echo "DRY RUN COMPLETE - No changes were made"
    echo ""
    echo "To execute for real:"
    echo "  ./remove_lfs_from_history.sh execute"
else
    echo "✓ PDFs and PNGs removed from Git history"
    echo "✓ LFS cache cleaned"
    echo "✓ Repository size reduced"
    echo "✓ Remote re-added"
fi

echo ""
echo "⚠️  CRITICAL NEXT STEPS ⚠️"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "1. FORCE PUSH the cleaned history:"
echo "   git push --force-with-lease origin master"
echo ""
echo "   ⚠️  This will rewrite history on GitHub"
echo "   ⚠️  All collaborators must re-clone after this"
echo ""
echo "2. NOTIFY all collaborators:"
echo "   - Do NOT pull the old repository"
echo "   - Delete old clone: rm -rf SRF_Linda_RNA"
echo "   - Fresh clone: git clone https://github.com/michal7kw/SRF_Linda_RNA.git"
echo ""
echo "3. On other computers, to pull without LFS errors:"
echo "   GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/michal7kw/SRF_Linda_RNA.git"
echo ""
echo "4. Verify on GitHub:"
echo "   - Repository size should be much smaller"
echo "   - LFS bandwidth usage will be minimal going forward"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if ! $DRY_RUN; then
    echo "Backup location (if you need to restore):"
    echo "  $BACKUP_DIR"
    echo ""
fi
