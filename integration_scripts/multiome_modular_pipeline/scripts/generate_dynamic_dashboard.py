#!/usr/bin/env python3
"""
Generate dynamic GCS dashboard based on actual bucket contents
Scans GCS bucket and creates HTML with all uploaded files
"""

import subprocess
import json
from datetime import datetime
from collections import defaultdict
import sys

# Configuration
GCS_BUCKET = "srf-multiome-figures"
BASE_URL = f"https://storage.googleapis.com/{GCS_BUCKET}"
CONSOLE_URL = f"https://console.cloud.google.com/storage/browser/{GCS_BUCKET}"

def run_command(cmd):
    """Run shell command and return output"""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip()

def get_file_list():
    """Get all files from GCS bucket"""
    print("Scanning GCS bucket...")
    cmd = f"gsutil ls -r gs://{GCS_BUCKET}/** 2>/dev/null"
    output = run_command(cmd)

    files = []
    for line in output.split('\n'):
        line = line.strip()
        if line and not line.endswith(':') and not line.endswith('/'):
            # Extract relative path
            path = line.replace(f"gs://{GCS_BUCKET}/", "")
            files.append(path)

    return files

def organize_files(files):
    """Organize files by directory"""
    structure = defaultdict(list)

    for file in files:
        # Skip non-visualization files
        ext = file.split('.')[-1].lower()
        if ext not in ['png', 'pdf', 'csv', 'bw']:
            continue

        # Get directory
        parts = file.split('/')
        if len(parts) > 1:
            directory = '/'.join(parts[:-1])
            filename = parts[-1]
            structure[directory].append({
                'name': filename,
                'path': file,
                'ext': ext,
                'url': f"{BASE_URL}/{file}"
            })

    return structure

def get_directory_description(directory):
    """Get human-readable description for directory"""
    descriptions = {
        'multiome/CREs_literature/heatmaps_deeptools': {
            'title': '🔥 CREs Heatmaps (deepTools)',
            'desc': 'Cell-type specificity validation: GABA vs Excitatory CREs'
        },
        'multiome/CREs_literature/heatmaps_specific_CREs': {
            'title': '🎯 Specific CREs Analysis',
            'desc': 'Targeted analysis of specific regulatory elements'
        },
        'multiome/L1_celltype_results/DEG': {
            'title': '📊 Differential Expression (DEGs)',
            'desc': 'L1 cell type differential gene expression results'
        },
        'multiome/L1_celltype_results/DA': {
            'title': '📈 Differential Accessibility (DA Peaks)',
            'desc': 'L1 cell type differential chromatin accessibility'
        },
        'multiome/L1_bigwig_tracks/by_celltype': {
            'title': '📊 BigWig Tracks (IGV)',
            'desc': 'ATAC-seq coverage tracks for genome browser (943 MB, 35 files)'
        },
        'multiome/L1_peak_gene_linkage/GABA': {
            'title': '🔗 Peak-Gene Linkage (GABA)',
            'desc': 'Regulatory network analysis for GABA neurons'
        }
    }

    return descriptions.get(directory, {
        'title': f'📁 {directory.split("/")[-1]}',
        'desc': directory
    })

def generate_html(structure, stats):
    """Generate HTML dashboard"""

    # Start HTML
    html = '''<!DOCTYPE html>
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
            max-width: 1600px;
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
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }
        .stat-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 15px;
            text-align: center;
        }
        .stat-number { font-size: 2.5em; font-weight: bold; }
        .stat-label { font-size: 0.9em; opacity: 0.9; margin-top: 5px; }

        .section {
            margin-bottom: 50px;
            border-left: 4px solid #667eea;
            padding-left: 25px;
        }
        .section h2 {
            color: #667eea;
            margin-bottom: 15px;
            font-size: 1.9em;
        }
        .section p {
            color: #666;
            margin-bottom: 20px;
            line-height: 1.6;
            font-size: 1.05em;
        }

        .folder-path {
            font-family: 'Monaco', 'Courier New', monospace;
            background: #edf2f7;
            padding: 10px 15px;
            border-radius: 8px;
            font-size: 0.9em;
            color: #4a5568;
            margin: 15px 0;
            display: inline-block;
        }

        .file-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
            gap: 15px;
            margin-top: 25px;
        }

        .file-card {
            background: #f7fafc;
            border: 2px solid #e2e8f0;
            border-radius: 12px;
            padding: 20px;
            transition: all 0.3s;
            cursor: pointer;
            text-decoration: none;
            display: block;
            color: inherit;
        }
        .file-card:hover {
            border-color: #667eea;
            box-shadow: 0 4px 12px rgba(102, 126, 234, 0.3);
            transform: translateY(-3px);
        }

        .file-icon {
            width: 55px;
            height: 55px;
            border-radius: 12px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 1.6em;
            margin-bottom: 15px;
        }
        .icon-png { background: #fed7d7; }
        .icon-pdf { background: #bee3f8; }
        .icon-csv { background: #c6f6d5; }
        .icon-bw { background: #fbd38d; }

        .file-name {
            font-family: 'Monaco', 'Courier New', monospace;
            font-size: 0.85em;
            color: #2d3748;
            word-break: break-word;
            margin-bottom: 10px;
            font-weight: 500;
        }

        .file-type-badge {
            display: inline-block;
            padding: 5px 12px;
            border-radius: 14px;
            font-size: 0.7em;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        .badge-png { background: #fed7d7; color: #c53030; }
        .badge-pdf { background: #bee3f8; color: #2c5282; }
        .badge-csv { background: #c6f6d5; color: #2f855a; }
        .badge-bw { background: #fbd38d; color: #c05621; }

        .browse-link {
            display: inline-block;
            background: #667eea;
            color: white;
            padding: 12px 24px;
            border-radius: 10px;
            text-decoration: none;
            margin: 10px 10px 10px 0;
            transition: all 0.3s;
            font-weight: 500;
        }
        .browse-link:hover {
            background: #764ba2;
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(102, 126, 234, 0.4);
        }

        .info-box {
            background: #f7fafc;
            border-left: 4px solid #48bb78;
            padding: 25px;
            margin: 25px 0;
            border-radius: 10px;
        }
        .info-box h3 {
            color: #48bb78;
            margin-bottom: 12px;
            font-size: 1.3em;
        }
        .info-box p {
            color: #4a5568;
            line-height: 1.7;
        }
        code {
            background: #edf2f7;
            padding: 3px 8px;
            border-radius: 5px;
            font-family: 'Monaco', 'Courier New', monospace;
            font-size: 0.9em;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 SRF Multiome Analysis</h1>
            <p>Complete Analysis Results Dashboard</p>
            <p class="update-time">Last updated: ''' + stats['timestamp'] + '''</p>
        </div>

        <div class="content">
            <div class="stats">
                <div class="stat-card">
                    <div class="stat-number">''' + str(stats['total_files']) + '''</div>
                    <div class="stat-label">Total Files</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">''' + stats['total_size'] + '''</div>
                    <div class="stat-label">Storage Used</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">''' + str(stats['total_dirs']) + '''</div>
                    <div class="stat-label">Categories</div>
                </div>
            </div>
'''

    # Add sections for each directory
    # Sort directories by priority
    priority_dirs = [
        'multiome/CREs_literature/heatmaps_deeptools',
        'multiome/CREs_literature/heatmaps_specific_CREs',
        'multiome/L1_celltype_results/DEG',
        'multiome/L1_celltype_results/DA',
        'multiome/L1_bigwig_tracks/by_celltype',
        'multiome/L1_peak_gene_linkage/GABA'
    ]

    # Process priority directories first
    processed = set()
    for directory in priority_dirs:
        if directory in structure and directory not in processed:
            info = get_directory_description(directory)
            files = structure[directory]

            html += f'''
            <div class="section">
                <h2>{info['title']}</h2>
                <p>{info['desc']}</p>
                <div class="folder-path">{directory}/</div>
                <a href="{CONSOLE_URL}/{directory}" target="_blank" class="browse-link">📁 Browse in GCS Console</a>
                <div class="file-grid">
'''

            # Add files
            for file in files[:20]:  # Limit to 20 per section
                icon_emoji = {'png': '🖼️', 'pdf': '📄', 'csv': '📊', 'bw': '📈'}
                emoji = icon_emoji.get(file['ext'], '📄')

                html += f'''
                    <a href="{file['url']}" target="_blank" class="file-card">
                        <div class="file-icon icon-{file['ext']}">{emoji}</div>
                        <div class="file-name">{file['name']}</div>
                        <span class="file-type-badge badge-{file['ext']}">{file['ext'].upper()}</span>
                    </a>
'''

            if len(files) > 20:
                html += f'''
                    <div style="grid-column: 1/-1; text-align: center; padding: 20px; color: #666;">
                        ... and {len(files) - 20} more files
                        <a href="{CONSOLE_URL}/{directory}" target="_blank" style="color: #667eea; font-weight: bold;"> → Browse all in GCS Console</a>
                    </div>
'''

            html += '''
                </div>
            </div>
'''
            processed.add(directory)

    # Add remaining directories
    for directory in sorted(structure.keys()):
        if directory not in processed:
            info = get_directory_description(directory)
            files = structure[directory]

            html += f'''
            <div class="section">
                <h2>{info['title']}</h2>
                <p>{info['desc']}</p>
                <div class="folder-path">{directory}/</div>
                <a href="{CONSOLE_URL}/{directory}" target="_blank" class="browse-link">📁 Browse in GCS Console</a>
                <div class="file-grid">
'''

            for file in files[:20]:
                icon_emoji = {'png': '🖼️', 'pdf': '📄', 'csv': '📊', 'bw': '📈'}
                emoji = icon_emoji.get(file['ext'], '📄')

                html += f'''
                    <a href="{file['url']}" target="_blank" class="file-card">
                        <div class="file-icon icon-{file['ext']}">{emoji}</div>
                        <div class="file-name">{file['name']}</div>
                        <span class="file-type-badge badge-{file['ext']}">{file['ext'].upper()}</span>
                    </a>
'''

            if len(files) > 20:
                html += f'''
                    <div style="grid-column: 1/-1; text-align: center; padding: 20px; color: #666;">
                        ... and {len(files) - 20} more files
                        <a href="{CONSOLE_URL}/{directory}" target="_blank" style="color: #667eea; font-weight: bold;"> → Browse all in GCS Console</a>
                    </div>
'''

            html += '''
                </div>
            </div>
'''

    # Add info box at bottom
    html += '''
            <div class="info-box">
                <h3>ℹ️ About This Dashboard</h3>
                <p>
                    This dashboard provides public access to SRF multiome analysis results.
                    All files are hosted on Google Cloud Storage with direct HTTP access.
                    <br><br>
                    <strong>Base URL:</strong>
                    <code>https://storage.googleapis.com/''' + GCS_BUCKET + '''</code>
                    <br><br>
                    💡 <strong>Tip:</strong> All URLs are stable and can be shared directly with collaborators.
                    Click any file card to view/download, or use "Browse Directory" to see all files.
                </p>
            </div>
        </div>
    </div>
</body>
</html>'''

    return html

def main():
    print("════════════════════════════════════════════════════════════════")
    print("  Dynamic GCS Dashboard Generator")
    print("════════════════════════════════════════════════════════════════")
    print()

    # Get files from GCS
    files = get_file_list()
    print(f"Found {len(files)} files in bucket")

    # Organize by directory
    structure = organize_files(files)
    print(f"Organized into {len(structure)} directories")

    # Count visualization files
    viz_files = sum(len(f) for f in structure.values())
    print(f"Visualization files (PNG/PDF/CSV/BW): {viz_files}")
    print()

    # Get stats
    total_size = run_command(f"gsutil du -sh gs://{GCS_BUCKET}/ 2>/dev/null | awk '{{print $1}}'")
    stats = {
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'total_files': viz_files,
        'total_size': total_size,
        'total_dirs': len(structure)
    }

    # Generate HTML
    print("Generating HTML dashboard...")
    html = generate_html(structure, stats)

    # Save locally
    with open('/tmp/gcs_dashboard_dynamic.html', 'w') as f:
        f.write(html)
    print("✓ HTML generated")

    # Upload to GCS
    print("Uploading to GCS...")
    subprocess.run(f"gsutil cp /tmp/gcs_dashboard_dynamic.html gs://{GCS_BUCKET}/index.html", shell=True)
    subprocess.run(f"gsutil setmeta -h 'Content-Type:text/html' -h 'Cache-Control:no-cache, must-revalidate' gs://{GCS_BUCKET}/index.html", shell=True)
    print("✓ Dashboard uploaded")
    print()

    print("════════════════════════════════════════════════════════════════")
    print("  Dashboard Updated! ✓")
    print("════════════════════════════════════════════════════════════════")
    print()
    print(f"Access your dashboard:")
    print(f"  https://storage.googleapis.com/{GCS_BUCKET}/index.html")
    print()
    print(f"Statistics:")
    print(f"  Files: {viz_files}")
    print(f"  Directories: {len(structure)}")
    print(f"  Size: {total_size}")
    print(f"  Updated: {stats['timestamp']}")
    print()
    print("💡 Hard refresh your browser to see updates:")
    print("   Ctrl+Shift+R (Windows/Linux) or Cmd+Shift+R (Mac)")
    print()

if __name__ == '__main__':
    main()
