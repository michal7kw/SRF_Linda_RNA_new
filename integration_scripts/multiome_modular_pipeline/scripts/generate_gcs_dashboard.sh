#!/bin/bash
#
# generate_gcs_dashboard.sh
#
# Generate enhanced HTML dashboard for GCS bucket
# Scans bucket contents and creates organized file browser
#
# Usage:
#   ./generate_gcs_dashboard.sh
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source configuration
source "${SCRIPT_DIR}/gcs_config.sh"

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[INFO]${NC} $*"; }
log_step() { echo -e "${BLUE}[STEP]${NC} $*"; }

echo "════════════════════════════════════════════════════════════════"
echo "  Generate GCS Dashboard"
echo "════════════════════════════════════════════════════════════════"
echo ""

# Validate GCS configuration
validate_config

log_step "Scanning GCS bucket contents"
echo ""

# Get statistics
log_info "Counting files and calculating sizes..."
total_files=$(gsutil ls -r gs://${GCS_BUCKET}/** 2>/dev/null | grep -v "/$" | wc -l)
total_size=$(gsutil du -sh gs://${GCS_BUCKET}/ 2>/dev/null | awk '{print $1}')
timestamp=$(date '+%Y-%m-%d %H:%M:%S')

log_info "Total files: $total_files"
log_info "Total size: $total_size"
echo ""

log_step "Generating HTML dashboard"
echo ""

# Create enhanced dashboard
cat > /tmp/gcs_dashboard_new.html << 'EOFHTML'
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
            max-width: 1400px;
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
        .update-time {
            font-size: 0.9em;
            margin-top: 10px;
            opacity: 0.8;
        }
        .content { padding: 40px; }
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

        .file-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
            gap: 15px;
            margin-top: 20px;
        }

        .file-card {
            background: #f7fafc;
            border: 2px solid #e2e8f0;
            border-radius: 12px;
            padding: 20px;
            transition: all 0.3s;
            cursor: pointer;
        }
        .file-card:hover {
            border-color: #667eea;
            box-shadow: 0 4px 12px rgba(102, 126, 234, 0.2);
            transform: translateY(-2px);
        }

        .file-icon {
            width: 50px;
            height: 50px;
            border-radius: 10px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 1.5em;
            margin-bottom: 12px;
        }
        .icon-png { background: #fed7d7; }
        .icon-pdf { background: #bee3f8; }
        .icon-csv { background: #c6f6d5; }
        .icon-bw { background: #fbd38d; }

        .file-name {
            font-family: 'Monaco', 'Courier New', monospace;
            font-size: 0.85em;
            color: #2d3748;
            word-break: break-all;
            margin-bottom: 8px;
        }

        .file-type-badge {
            display: inline-block;
            padding: 4px 10px;
            border-radius: 12px;
            font-size: 0.7em;
            font-weight: 600;
            text-transform: uppercase;
        }
        .badge-png { background: #fed7d7; color: #c53030; }
        .badge-pdf { background: #bee3f8; color: #2c5282; }
        .badge-csv { background: #c6f6d5; color: #2f855a; }
        .badge-bw { background: #fbd38d; color: #c05621; }

        .browse-link {
            display: inline-block;
            background: #667eea;
            color: white;
            padding: 10px 20px;
            border-radius: 8px;
            text-decoration: none;
            margin: 10px 10px 10px 0;
            transition: all 0.3s;
        }
        .browse-link:hover {
            background: #764ba2;
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(102, 126, 234, 0.3);
        }

        a { color: inherit; text-decoration: none; }

        .info-box {
            background: #f7fafc;
            border-left: 4px solid #48bb78;
            padding: 20px;
            margin: 20px 0;
            border-radius: 8px;
        }
        .info-box h3 {
            color: #48bb78;
            margin-bottom: 10px;
        }
        .info-box p {
            color: #4a5568;
            line-height: 1.6;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 SRF Multiome Analysis</h1>
            <p>Complete Analysis Results Dashboard</p>
            <p class="update-time">Last updated: TIMESTAMP</p>
        </div>

        <div class="content">
            <div class="stats">
                <div class="stat-card">
                    <div class="stat-number">TOTAL_FILES</div>
                    <div class="stat-label">Total Files</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">TOTAL_SIZE</div>
                    <div class="stat-label">Storage Used</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">5+</div>
                    <div class="stat-label">Analysis Categories</div>
                </div>
            </div>

            <!-- Main Signac Results -->
            <div class="section">
                <h2>📊 Main Signac Results (L2)</h2>
                <p>Primary multiome integration with fine-grained cell types</p>
                <div class="folder-path">multiome/plots/latest/</div>
                <a href="multiome/plots/latest/" target="_blank" class="browse-link">📁 Browse All Files</a>
                <div class="file-grid">
                    <a href="multiome/plots/latest/03_DEG_volcano_plots.pdf" target="_blank">
                        <div class="file-card">
                            <div class="file-icon icon-pdf">📄</div>
                            <div class="file-name">DEG Volcano Plots</div>
                            <span class="file-type-badge badge-pdf">PDF</span>
                        </div>
                    </a>
                    <a href="multiome/plots/latest/04_DA_volcano_plots.pdf" target="_blank">
                        <div class="file-card">
                            <div class="file-icon icon-pdf">📄</div>
                            <div class="file-name">DA Peak Volcano Plots</div>
                            <span class="file-type-badge badge-pdf">PDF</span>
                        </div>
                    </a>
                </div>
            </div>

            <!-- CREs from Literature -->
            <div class="section">
                <h2>🔥 CREs from Literature</h2>
                <p>Cell-type specificity validation using literature-curated regulatory elements</p>
                <div class="folder-path">multiome/CREs_literature/heatmaps_deeptools/</div>
                <a href="multiome/CREs_literature/heatmaps_deeptools/" target="_blank" class="browse-link">📁 Browse Heatmaps</a>
                <div class="file-grid">
                    <a href="multiome/CREs_literature/heatmaps_deeptools/heatmap_GABA_all_conditions.png" target="_blank">
                        <div class="file-card">
                            <div class="file-icon icon-png">🖼️</div>
                            <div class="file-name">GABA CREs - All Conditions</div>
                            <span class="file-type-badge badge-png">PNG</span>
                        </div>
                    </a>
                    <a href="multiome/CREs_literature/heatmaps_deeptools/heatmap_Excitatory_all_conditions.png" target="_blank">
                        <div class="file-card">
                            <div class="file-icon icon-png">🖼️</div>
                            <div class="file-name">Excitatory CREs - All Conditions</div>
                            <span class="file-type-badge badge-png">PNG</span>
                        </div>
                    </a>
                    <a href="multiome/CREs_literature/heatmaps_deeptools/metaprofile_GABA_all_conditions.png" target="_blank">
                        <div class="file-card">
                            <div class="file-icon icon-png">🖼️</div>
                            <div class="file-name">GABA Metaprofile</div>
                            <span class="file-type-badge badge-png">PNG</span>
                        </div>
                    </a>
                </div>
            </div>

            <!-- L1 Cell Type Results -->
            <div class="section">
                <h2>📈 L1 Cell Type Analysis</h2>
                <p>Broad cell type results: DEGs, DA peaks, and peak-gene linkages</p>
                <div class="folder-path">multiome/L1_celltype_results/</div>
                <a href="multiome/L1_celltype_results/DEG/" target="_blank" class="browse-link">📊 DEGs (CSV)</a>
                <a href="multiome/L1_celltype_results/DA/" target="_blank" class="browse-link">📊 DA Peaks (CSV)</a>
                <a href="multiome/L1_celltype_results/peak_gene_links/" target="_blank" class="browse-link">🔗 Peak-Gene Links</a>
            </div>

            <!-- BigWig Tracks -->
            <div class="section">
                <h2>📊 BigWig Tracks (IGV)</h2>
                <p>ATAC-seq coverage tracks for genome browser visualization (943 MB, 35 files)</p>
                <div class="folder-path">multiome/L1_bigwig_tracks/by_celltype/</div>
                <a href="multiome/L1_bigwig_tracks/by_celltype/" target="_blank" class="browse-link">📁 Browse BigWig Files</a>
                <div class="info-box">
                    <h3>💡 Load in IGV Browser</h3>
                    <p>
                        <strong>Step 1:</strong> Open IGV (Integrative Genomics Viewer)<br>
                        <strong>Step 2:</strong> File → Load from URL<br>
                        <strong>Step 3:</strong> Paste BigWig URL (e.g.,
                        <code style="background: #edf2f7; padding: 2px 6px; border-radius: 4px;">
                            https://storage.googleapis.com/srf-multiome-figures/multiome/L1_bigwig_tracks/by_celltype/GABA_Nestin-Ctrl.bw
                        </code>)
                    </p>
                </div>
            </div>

            <!-- Peak-Gene Linkage -->
            <div class="section">
                <h2>🔗 Peak-Gene Linkage Analysis</h2>
                <p>Regulatory network analysis for GABA neurons</p>
                <div class="folder-path">multiome/L1_peak_gene_linkage/GABA/</div>
                <a href="multiome/L1_peak_gene_linkage/GABA/" target="_blank" class="browse-link">📁 Browse Analysis Results</a>
            </div>

            <!-- About Section -->
            <div class="info-box">
                <h3>ℹ️ About This Dashboard</h3>
                <p>
                    This dashboard provides public access to SRF multiome analysis results.
                    All files are hosted on Google Cloud Storage with direct HTTP access.
                    <br><br>
                    <strong>Base URL:</strong>
                    <code style="background: #edf2f7; padding: 4px 8px; border-radius: 4px; font-family: monospace;">
                        https://storage.googleapis.com/srf-multiome-figures/
                    </code>
                    <br><br>
                    💡 <strong>Tip:</strong> All URLs are stable and can be shared directly with collaborators.
                    Click any file card or "Browse" button to explore directories.
                </p>
            </div>
        </div>
    </div>
</body>
</html>
EOFHTML

# Replace placeholders
sed -i "s|TIMESTAMP|$timestamp|g" /tmp/gcs_dashboard_new.html
sed -i "s|TOTAL_FILES|$total_files|g" /tmp/gcs_dashboard_new.html
sed -i "s|TOTAL_SIZE|$total_size|g" /tmp/gcs_dashboard_new.html

log_info "Dashboard HTML generated"
echo ""

log_step "Uploading dashboard to GCS"
echo ""

# Upload with proper metadata
gsutil cp /tmp/gcs_dashboard_new.html gs://${GCS_BUCKET}/index.html
gsutil setmeta -h "Content-Type:text/html" -h "Cache-Control:no-cache, must-revalidate" gs://${GCS_BUCKET}/index.html

log_info "✓ Dashboard uploaded"
log_info "✓ Cache headers set (no-cache)"
echo ""

echo "════════════════════════════════════════════════════════════════"
echo "  Dashboard Updated! ✓"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "Access your dashboard:"
echo "  https://storage.googleapis.com/srf-multiome-figures/index.html"
echo ""
echo "Statistics:"
echo "  Files: $total_files"
echo "  Size: $total_size"
echo "  Updated: $timestamp"
echo ""
echo "💡 Tip: If you see old content, hard refresh your browser:"
echo "   - Chrome/Firefox: Ctrl+Shift+R (Windows/Linux) or Cmd+Shift+R (Mac)"
echo "   - Safari: Cmd+Option+R"
echo ""
