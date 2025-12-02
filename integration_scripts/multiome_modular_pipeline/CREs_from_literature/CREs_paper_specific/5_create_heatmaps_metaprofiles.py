#!/usr/bin/env python3
"""
Create Heatmaps and Metaprofiles: Cell-Type Specificity Analysis

Working dir: "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature"

This script demonstrates cell-type specificity of ATAC-seq signal by comparing:
1. GABA-specific CREs (POSITIVE CONTROL - expect HIGH signal in GABA samples)
2. Excitatory neuron CREs (NEGATIVE CONTROL - expect LOW signal in GABA samples)

Workflow:
1. Loads both GABA and Excitatory CREs from pre-extracted BED files
2. Extracts ATAC signal from BigWig files for all 4 conditions
3. Creates heatmaps with common scale for both CRE sets
4. Creates metaprofiles (average signal ± SEM) with common Y-axis
5. Generates comparison plots demonstrating cell-type specificity

Conditions analyzed:
- Nestin-Ctrl
- Nestin-Mut
- Emx1-Ctrl
- Emx1-Mut

INPUT FILES:
- output/hippocampal_interneuron_CREs.bed: GABA-specific CREs (positive control)
- output/excitatory_neuron_CREs.bed: Excitatory neuron CREs (negative control)
- ../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw: ATAC signal files (4 conditions)

OUTPUT FILES:
- output/GABA_DEG_analysis/heatmaps_metaprofiles/heatmap_GABA_all_conditions.png: GABA CREs heatmap
- output/GABA_DEG_analysis/heatmaps_metaprofiles/heatmap_Excitatory_all_conditions.png: Excitatory CREs heatmap
- output/GABA_DEG_analysis/heatmaps_metaprofiles/metaprofile_GABA_all_CREs.png: GABA metaprofile
- output/GABA_DEG_analysis/heatmaps_metaprofiles/metaprofile_Excitatory_all_CREs.png: Excitatory metaprofile
- output/GABA_DEG_analysis/heatmaps_metaprofiles/comparison_GABA_vs_Excitatory.png: Side-by-side comparison
- output/GABA_DEG_analysis/heatmaps_metaprofiles/enrichment_statistics.txt: Quantitative comparison
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import os
from datetime import datetime

os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature")

# Try to import pyBigWig
try:
    import pyBigWig
    PYBIGWIG_AVAILABLE = True
except ImportError:
    print("⚠️  pyBigWig not available. Install with: pip install pyBigWig")
    PYBIGWIG_AVAILABLE = False
    exit(1)

print("="*80)
print("CREATE HEATMAPS AND METAPROFILES: ATAC SIGNAL AT GABA CREs")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# ============================================================================
# Configuration
# ============================================================================
BIGWIG_BASE = "../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR = "output/GABA_DEG_analysis/heatmaps_metaprofiles"
BED_GABA = "output/hippocampal_interneuron_CREs.bed"
BED_EXCITATORY = "output/excitatory_neuron_CREs.bed"

# Analysis parameters
WINDOW_SIZE = 2000      # bp around CRE center (±2kb = 4kb total)
BIN_SIZE = 50           # bp per bin (4kb / 50bp = 80 bins)
N_BINS = WINDOW_SIZE * 2 // BIN_SIZE

# Use 90th percentile for scale (robust against outliers)
USE_PERCENTILE_SCALE = True
SCALE_PERCENTILE = 90
SCALE_BUFFER = 1.2      # 20% buffer above percentile

# Filter parameters
MIN_SIGNAL_THRESHOLD = 0.05  # Minimum signal in any condition

# Conditions to analyze
# NOTE: Emx1-Ctrl is excluded (failed sample)
CONDITIONS = {
    "Nestin-Ctrl": f"{BIGWIG_BASE}/GABA_Nestin-Ctrl.bw",
    "Nestin-Mut": f"{BIGWIG_BASE}/GABA_Nestin-Mut.bw",
    "Emx1-Mut": f"{BIGWIG_BASE}/GABA_Emx1-Mut.bw"
}

CONDITION_COLORS = {
    "Nestin-Ctrl": "#2E86AB",
    "Nestin-Mut": "#A23B72",
    "Emx1-Mut": "#C73E1D"
}

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Plotting style
sns.set_style("white")
plt.rcParams['figure.dpi'] = 150

# ============================================================================
# STEP 1: Load CREs from BED files (GABA + Excitatory)
# ============================================================================
print("STEP 1: Loading CREs from BED files...")
print("-"*80)

def load_bed_file(bed_file, cre_type="Unknown"):
    """Load CREs from BED file"""
    print(f"\nLoading {cre_type} CREs from: {bed_file}")

    if not os.path.exists(bed_file):
        print(f"  ✗ File not found: {bed_file}")
        return None

    # Read BED file (chr, start, end, name, ...)
    cres = pd.read_csv(bed_file, sep='\t', header=None,
                       names=['cre_chr', 'cre_start', 'cre_end', 'cre_id', 'extra1', 'extra2'])

    # Calculate CRE centers
    cres['cre_center'] = ((cres['cre_start'] + cres['cre_end']) / 2).astype(int)

    # Keep only necessary columns
    cres = cres[['cre_chr', 'cre_start', 'cre_end', 'cre_center', 'cre_id']]

    # Sort by chromosome and position
    cres = cres.sort_values(['cre_chr', 'cre_center']).reset_index(drop=True)

    print(f"  ✓ Loaded {len(cres)} CREs")
    print(f"  Chromosomes: {sorted(cres['cre_chr'].unique())}")

    return cres

# Load GABA CREs (positive control)
cres_gaba = load_bed_file(BED_GABA, "GABA (POSITIVE CONTROL)")

# Load Excitatory CREs (negative control)
cres_excitatory = load_bed_file(BED_EXCITATORY, "Excitatory (NEGATIVE CONTROL)")

# Check if both loaded successfully
if cres_gaba is None or cres_excitatory is None:
    print("\n⚠️  ERROR: Could not load both BED files!")
    print("Make sure to run:")
    print("  1. 1_extract_hippocampal_interneuron_CREs.py")
    print("  2. 2a_extract_excitatory_neuron_CREs.py")
    exit(1)

print(f"\n✓ Successfully loaded both CRE sets:")
print(f"  GABA CREs: {len(cres_gaba)}")
print(f"  Excitatory CREs: {len(cres_excitatory)}")

# ============================================================================
# STEP 2: Extract signal matrices from BigWig files
# ============================================================================
print(f"\n{'='*80}")
print("STEP 2: Extracting signal matrices from BigWig files...")
print("-"*80)

def extract_signal_window(bw_file, chrom, center, window_size=2000, bin_size=50):
    """
    Extract signal in bins around a genomic position

    Returns array of length (window_size*2 // bin_size)
    """
    try:
        bw = pyBigWig.open(bw_file)

        if chrom not in bw.chroms():
            bw.close()
            return np.full(N_BINS, np.nan)

        # Define window
        start = max(0, center - window_size)
        end = center + window_size

        # Ensure we don't exceed chromosome length
        chr_len = bw.chroms()[chrom]
        end = min(end, chr_len)

        # Extract signal in bins
        n_bins = (end - start) // bin_size
        signals = []

        for i in range(n_bins):
            bin_start = start + (i * bin_size)
            bin_end = min(bin_start + bin_size, end)

            signal = bw.stats(chrom, bin_start, bin_end, type="mean")

            if signal and signal[0] is not None:
                signals.append(signal[0])
            else:
                signals.append(0.0)

        bw.close()

        # Pad to correct length if needed
        while len(signals) < N_BINS:
            signals.append(0.0)

        return np.array(signals[:N_BINS])

    except Exception as e:
        print(f"    Error extracting signal: {e}")
        return np.full(N_BINS, np.nan)

# Extract signal matrices for all conditions
signal_matrices = {}

for condition_name, bw_file in CONDITIONS.items():
    print(f"\n{condition_name}:")

    if not os.path.exists(bw_file):
        print(f"  ✗ File not found: {bw_file}")
        print(f"  Skipping this condition")
        continue

    print(f"  Extracting signals from {len(cres)} CREs...")

    # Initialize matrix: rows = CREs, columns = bins
    matrix = np.zeros((len(cres), N_BINS))

    for idx, row in cres.iterrows():
        if idx % 100 == 0:
            print(f"    Processing {idx}/{len(cres)}...", end='\r')

        signals = extract_signal_window(
            bw_file,
            row['cre_chr'],
            row['cre_center'],
            window_size=WINDOW_SIZE,
            bin_size=BIN_SIZE
        )

        matrix[idx, :] = signals

    print(f"    Processing {len(cres)}/{len(cres)}... Done!")

    signal_matrices[condition_name] = matrix

    # Summary statistics
    print(f"  Mean signal: {np.nanmean(matrix):.3f}")
    print(f"  Median signal: {np.nanmedian(matrix):.3f}")
    print(f"  Max signal: {np.nanmax(matrix):.3f}")

if len(signal_matrices) == 0:
    print("\n⚠️  ERROR: No BigWig files found!")
    exit(1)

print(f"\n✓ Extracted signals for {len(signal_matrices)} conditions")

# ============================================================================
# STEP 3: Filter CREs by signal threshold
# ============================================================================
print(f"\n{'='*80}")
print("STEP 3: Filtering CREs by signal threshold...")
print("-"*80)

# Calculate maximum signal across all conditions and bins for each CRE
max_signals = np.zeros(len(cres))

for condition_name, matrix in signal_matrices.items():
    max_signals = np.maximum(max_signals, np.nanmax(matrix, axis=1))

# Filter
keep_mask = max_signals >= MIN_SIGNAL_THRESHOLD
cres_filtered = cres[keep_mask].reset_index(drop=True)
signal_matrices_filtered = {
    cond: matrix[keep_mask] for cond, matrix in signal_matrices.items()
}

print(f"Original CREs: {len(cres)}")
print(f"Filtered CREs (signal >= {MIN_SIGNAL_THRESHOLD}): {len(cres_filtered)}")
print(f"Removed: {len(cres) - len(cres_filtered)} ({100*(1-len(cres_filtered)/len(cres)):.1f}%)")

# Update variables
cres = cres_filtered
signal_matrices = signal_matrices_filtered

# ============================================================================
# STEP 4: Cluster CREs by signal pattern
# ============================================================================
print(f"\n{'='*80}")
print("STEP 4: Clustering CREs by signal pattern...")
print("-"*80)

# Use one condition for clustering (e.g., Nestin-Ctrl)
if CLUSTER_ON in signal_matrices:
    print(f"Clustering based on {CLUSTER_ON} signal...")

    cluster_matrix = signal_matrices[CLUSTER_ON]

    # Remove CREs with too many NaNs
    valid_rows = ~np.isnan(cluster_matrix).any(axis=1)
    cluster_matrix_clean = cluster_matrix[valid_rows]

    # Standardize signals for clustering
    scaler = StandardScaler()
    cluster_matrix_scaled = scaler.fit_transform(cluster_matrix_clean)

    # K-means clustering
    print(f"  Running K-means with {N_CLUSTERS} clusters...")
    kmeans = KMeans(n_clusters=N_CLUSTERS, random_state=42, n_init=10)
    cluster_labels_clean = kmeans.fit_predict(cluster_matrix_scaled)

    # Map back to original CREs
    cluster_labels = np.full(len(cres), -1)
    cluster_labels[valid_rows] = cluster_labels_clean

    cres['cluster'] = cluster_labels

    print(f"  Cluster sizes:")
    for i in range(N_CLUSTERS):
        n_in_cluster = (cluster_labels == i).sum()
        print(f"    Cluster {i+1}: {n_in_cluster} CREs ({100*n_in_cluster/len(cres):.1f}%)")

    print(f"\n  Cluster characteristics (mean signal in {CLUSTER_ON}):")
    for i in range(N_CLUSTERS):
        cluster_mask = cluster_labels == i
        cluster_mean_signal = np.nanmean(cluster_matrix[cluster_mask])
        print(f"    Cluster {i+1}: {cluster_mean_signal:.3f}")
else:
    print(f"⚠️  {CLUSTER_ON} not available, skipping clustering")
    cres['cluster'] = 0
    cluster_labels = np.zeros(len(cres))

# ============================================================================
# STEP 5: Sort CREs by signal intensity
# ============================================================================
print(f"\n{'='*80}")
print("STEP 5: Sorting CREs by signal intensity...")
print("-"*80)

if CLUSTER_ON in signal_matrices:
    # Calculate mean signal across window for each CRE
    mean_signals = np.nanmean(signal_matrices[CLUSTER_ON], axis=1)
    sort_idx = np.argsort(mean_signals)[::-1]  # Descending order

    cres_sorted = cres.iloc[sort_idx].reset_index(drop=True)
    signal_matrices_sorted = {
        cond: matrix[sort_idx] for cond, matrix in signal_matrices.items()
    }

    print(f"✓ Sorted {len(cres)} CREs by {CLUSTER_ON} signal intensity")
else:
    cres_sorted = cres
    signal_matrices_sorted = signal_matrices
    print(f"⚠️  Skipping sorting (no reference condition)")

# ============================================================================
# STEP 6: Create heatmaps
# ============================================================================
print(f"\n{'='*80}")
print("STEP 6: Creating heatmaps...")
print("-"*80)

# Heatmap 1: All conditions side-by-side
print("\nCreating combined heatmap (all conditions)...")

# NOTE: 3 conditions (Emx1-Ctrl excluded)
fig, axes = plt.subplots(1, 3, figsize=(12, 12), sharey=True)
fig.suptitle(f'ATAC Signal at GABA-Specific CREs (n={len(cres_sorted)})',
             fontsize=16, fontweight='bold', y=0.995)

for idx, (condition_name, matrix) in enumerate(signal_matrices_sorted.items()):
    ax = axes[idx]

    # Plot heatmap
    im = ax.imshow(matrix, aspect='auto', cmap='Reds',
                   vmin=0, vmax=np.nanpercentile(matrix, 99),
                   interpolation='nearest')

    # X-axis: position relative to CRE center
    x_ticks = [0, N_BINS//4, N_BINS//2, 3*N_BINS//4, N_BINS]
    x_labels = [f'-{WINDOW_SIZE//1000}kb', f'-{WINDOW_SIZE//2000}kb',
                'Center', f'+{WINDOW_SIZE//2000}kb', f'+{WINDOW_SIZE//1000}kb']
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, rotation=45, ha='right')

    # Y-axis
    if idx == 0:
        ax.set_ylabel('CREs (sorted by signal)', fontsize=12)

    ax.set_xlabel('Position', fontsize=12)
    ax.set_title(condition_name, fontsize=14, fontweight='bold')

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('ATAC Signal', fontsize=10)

plt.tight_layout()
plot_file = os.path.join(OUTPUT_DIR, "heatmap_all_conditions.png")
plt.savefig(plot_file, dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: heatmap_all_conditions.png")

# Heatmap 2: By cluster (if clustering was performed)
if 'cluster' in cres_sorted.columns and len(signal_matrices_sorted) > 0:
    print("\nCreating clustered heatmaps...")

    for cluster_id in range(N_CLUSTERS):
        cluster_mask = cres_sorted['cluster'] == cluster_id
        n_in_cluster = cluster_mask.sum()

        if n_in_cluster < 10:
            print(f"  Skipping cluster {cluster_id+1} (too few CREs: {n_in_cluster})")
            continue

        # NOTE: 3 conditions (Emx1-Ctrl excluded)
        fig, axes = plt.subplots(1, 3, figsize=(12, 8), sharey=True)
        fig.suptitle(f'ATAC Signal: Cluster {cluster_id+1} (n={n_in_cluster} CREs)',
                     fontsize=16, fontweight='bold', y=0.995)

        for idx, (condition_name, matrix) in enumerate(signal_matrices_sorted.items()):
            ax = axes[idx]

            cluster_matrix = matrix[cluster_mask]

            im = ax.imshow(cluster_matrix, aspect='auto', cmap='Reds',
                          vmin=0, vmax=np.nanpercentile(cluster_matrix, 99),
                          interpolation='nearest')

            x_ticks = [0, N_BINS//4, N_BINS//2, 3*N_BINS//4, N_BINS]
            x_labels = [f'-{WINDOW_SIZE//1000}kb', f'-{WINDOW_SIZE//2000}kb',
                       'Center', f'+{WINDOW_SIZE//2000}kb', f'+{WINDOW_SIZE//1000}kb']
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_labels, rotation=45, ha='right')

            if idx == 0:
                ax.set_ylabel('CREs', fontsize=12)

            ax.set_xlabel('Position', fontsize=12)
            ax.set_title(condition_name, fontsize=14, fontweight='bold')

            cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label('ATAC Signal', fontsize=10)

        plt.tight_layout()
        plot_file = os.path.join(OUTPUT_DIR, f"heatmap_cluster_{cluster_id+1}.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Saved: heatmap_cluster_{cluster_id+1}.png")

# ============================================================================
# STEP 7: Create metaprofiles (average profiles)
# ============================================================================
print(f"\n{'='*80}")
print("STEP 7: Creating metaprofiles...")
print("-"*80)

# Metaprofile 1: All CREs
print("\nCreating metaprofile for all CREs...")

fig, ax = plt.subplots(figsize=(10, 6))

x_coords = np.linspace(-WINDOW_SIZE, WINDOW_SIZE, N_BINS)

for condition_name, matrix in signal_matrices_sorted.items():
    # Calculate mean and SEM
    mean_signal = np.nanmean(matrix, axis=0)
    sem_signal = stats.sem(matrix, axis=0, nan_policy='omit')

    # Plot
    color = CONDITION_COLORS.get(condition_name, 'black')
    ax.plot(x_coords, mean_signal, color=color, linewidth=2.5,
            label=condition_name)
    ax.fill_between(x_coords,
                     mean_signal - sem_signal,
                     mean_signal + sem_signal,
                     color=color, alpha=0.2)

# Add center line
ax.axvline(0, color='gray', linestyle='--', alpha=0.5, linewidth=1)

ax.set_xlabel('Distance from CRE Center (bp)', fontsize=12)
ax.set_ylabel('Mean ATAC Signal', fontsize=12)
ax.set_title(f'Metaprofile: ATAC Signal at GABA CREs (n={len(cres_sorted)})',
             fontsize=14, fontweight='bold')
ax.legend(fontsize=11, frameon=True)
ax.grid(True, alpha=0.3, linestyle=':')

plt.tight_layout()
plot_file = os.path.join(OUTPUT_DIR, "metaprofile_all_CREs.png")
plt.savefig(plot_file, dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: metaprofile_all_CREs.png")

# Metaprofile 2: By cluster
if 'cluster' in cres_sorted.columns:
    print("\nCreating metaprofiles by cluster...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for cluster_id in range(N_CLUSTERS):
        ax = axes[cluster_id]
        cluster_mask = cres_sorted['cluster'] == cluster_id
        n_in_cluster = cluster_mask.sum()

        if n_in_cluster < 10:
            ax.text(0.5, 0.5, f'Cluster {cluster_id+1}\n(too few CREs: {n_in_cluster})',
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        for condition_name, matrix in signal_matrices_sorted.items():
            cluster_matrix = matrix[cluster_mask]

            mean_signal = np.nanmean(cluster_matrix, axis=0)
            sem_signal = stats.sem(cluster_matrix, axis=0, nan_policy='omit')

            color = CONDITION_COLORS.get(condition_name, 'black')
            ax.plot(x_coords, mean_signal, color=color, linewidth=2,
                   label=condition_name)
            ax.fill_between(x_coords,
                           mean_signal - sem_signal,
                           mean_signal + sem_signal,
                           color=color, alpha=0.2)

        ax.axvline(0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        ax.set_xlabel('Distance from CRE Center (bp)', fontsize=10)
        ax.set_ylabel('Mean ATAC Signal', fontsize=10)
        ax.set_title(f'Cluster {cluster_id+1} (n={n_in_cluster})', fontsize=12)
        ax.legend(fontsize=9, frameon=True)
        ax.grid(True, alpha=0.3, linestyle=':')

    plt.suptitle('Metaprofiles by Cluster', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plot_file = os.path.join(OUTPUT_DIR, "metaprofile_by_cluster.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: metaprofile_by_cluster.png")

# ============================================================================
# STEP 8: Create comparison plots
# ============================================================================
print(f"\n{'='*80}")
print("STEP 8: Creating comparison plots...")
print("-"*80)

# Plot 1: Nestin Ctrl vs Mut
print("\nCreating Nestin comparison...")

if "Nestin-Ctrl" in signal_matrices_sorted and "Nestin-Mut" in signal_matrices_sorted:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Metaprofile comparison
    ax = axes[0]
    for condition_name in ["Nestin-Ctrl", "Nestin-Mut"]:
        matrix = signal_matrices_sorted[condition_name]
        mean_signal = np.nanmean(matrix, axis=0)
        sem_signal = stats.sem(matrix, axis=0, nan_policy='omit')

        color = CONDITION_COLORS[condition_name]
        ax.plot(x_coords, mean_signal, color=color, linewidth=3,
               label=condition_name)
        ax.fill_between(x_coords,
                       mean_signal - sem_signal,
                       mean_signal + sem_signal,
                       color=color, alpha=0.2)

    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Distance from CRE Center (bp)', fontsize=12)
    ax.set_ylabel('Mean ATAC Signal', fontsize=12)
    ax.set_title('Nestin: Control vs Mutant', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, linestyle=':')

    # Difference plot
    ax = axes[1]
    ctrl_mean = np.nanmean(signal_matrices_sorted["Nestin-Ctrl"], axis=0)
    mut_mean = np.nanmean(signal_matrices_sorted["Nestin-Mut"], axis=0)
    diff = mut_mean - ctrl_mean

    ax.plot(x_coords, diff, color='purple', linewidth=3)
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.fill_between(x_coords, 0, diff, where=(diff < 0),
                     color='blue', alpha=0.3, label='Decreased in Mut')
    ax.fill_between(x_coords, 0, diff, where=(diff >= 0),
                     color='red', alpha=0.3, label='Increased in Mut')

    ax.set_xlabel('Distance from CRE Center (bp)', fontsize=12)
    ax.set_ylabel('Signal Difference (Mut - Ctrl)', fontsize=12)
    ax.set_title('Nestin: Signal Change', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, linestyle=':')

    plt.tight_layout()
    plot_file = os.path.join(OUTPUT_DIR, "comparison_nestin_ctrl_vs_mut.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: comparison_nestin_ctrl_vs_mut.png")

# Plot 2: Emx1 Ctrl vs Mut - SKIPPED (Emx1-Ctrl is a failed sample)
print("\nSkipping Emx1 comparison (Emx1-Ctrl excluded - failed sample)")
print("  NOTE: Only Emx1-Mut is available for analysis")

# ============================================================================
# STEP 9: Save signal matrices and CRE annotations
# ============================================================================
print(f"\n{'='*80}")
print("STEP 9: Saving signal matrices and annotations...")
print("-"*80)

# Save CRE list with cluster assignments
cres_out = cres_sorted[['cre_chr', 'cre_start', 'cre_end', 'cre_center', 'cre_id', 'cluster']].copy()
output_file = os.path.join(OUTPUT_DIR, "GABA_CREs_with_clusters.tsv")
cres_out.to_csv(output_file, sep='\t', index=False)
print(f"✓ Saved CRE annotations: GABA_CREs_with_clusters.tsv")

# Save signal matrices (as numpy arrays)
for condition_name, matrix in signal_matrices_sorted.items():
    output_file = os.path.join(OUTPUT_DIR, f"signal_matrix_{condition_name.replace('-', '_')}.npy")
    np.save(output_file, matrix)
    print(f"✓ Saved signal matrix: signal_matrix_{condition_name.replace('-', '_')}.npy")

# ============================================================================
# STEP 10: Create summary report
# ============================================================================
print(f"\n{'='*80}")
print("STEP 10: Creating summary report...")
print("-"*80)

report_file = os.path.join(OUTPUT_DIR, "analysis_summary.txt")
with open(report_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("HEATMAPS AND METAPROFILES: ATAC SIGNAL AT GABA CREs\n")
    f.write("="*80 + "\n\n")
    f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("PARAMETERS:\n")
    f.write("-"*80 + "\n")
    f.write(f"Window size: ±{WINDOW_SIZE} bp around CRE center\n")
    f.write(f"Bin size: {BIN_SIZE} bp\n")
    f.write(f"Number of bins: {N_BINS}\n")
    f.write(f"Min signal threshold: {MIN_SIGNAL_THRESHOLD}\n")
    f.write(f"Number of clusters: {N_CLUSTERS}\n\n")

    f.write("DATA:\n")
    f.write("-"*80 + "\n")
    f.write(f"Total GABA-specific CREs (from Table 16): {len(cres_filtered) + len(cres) - len(cres_filtered)}\n")
    f.write(f"After filtering (signal >= {MIN_SIGNAL_THRESHOLD}): {len(cres_sorted)}\n")
    f.write(f"Filtering rate: {100*(1-len(cres_sorted)/(len(cres_filtered) + len(cres) - len(cres_filtered))):.1f}%\n\n")

    f.write("CONDITIONS ANALYZED:\n")
    f.write("-"*80 + "\n")
    for condition_name in signal_matrices_sorted.keys():
        matrix = signal_matrices_sorted[condition_name]
        f.write(f"{condition_name}:\n")
        f.write(f"  Mean signal: {np.nanmean(matrix):.4f}\n")
        f.write(f"  Median signal: {np.nanmedian(matrix):.4f}\n")
        f.write(f"  Std signal: {np.nanstd(matrix):.4f}\n\n")

    if 'cluster' in cres_sorted.columns:
        f.write("CLUSTERS:\n")
        f.write("-"*80 + "\n")
        for i in range(N_CLUSTERS):
            cluster_mask = cres_sorted['cluster'] == i
            n_in_cluster = cluster_mask.sum()
            f.write(f"Cluster {i+1}: {n_in_cluster} CREs ({100*n_in_cluster/len(cres_sorted):.1f}%)\n")

            if CLUSTER_ON in signal_matrices_sorted:
                cluster_matrix = signal_matrices_sorted[CLUSTER_ON][cluster_mask]
                f.write(f"  Mean {CLUSTER_ON} signal: {np.nanmean(cluster_matrix):.4f}\n")
        f.write("\n")

    f.write("COMPARISONS:\n")
    f.write("-"*80 + "\n")

    # Nestin comparison
    if "Nestin-Ctrl" in signal_matrices_sorted and "Nestin-Mut" in signal_matrices_sorted:
        ctrl_mean = np.nanmean(signal_matrices_sorted["Nestin-Ctrl"])
        mut_mean = np.nanmean(signal_matrices_sorted["Nestin-Mut"])
        pct_change = 100 * (mut_mean - ctrl_mean) / ctrl_mean
        f.write(f"Nestin:\n")
        f.write(f"  Ctrl mean: {ctrl_mean:.4f}\n")
        f.write(f"  Mut mean: {mut_mean:.4f}\n")
        f.write(f"  Change: {pct_change:+.1f}%\n\n")

    # Emx1 comparison - SKIPPED (Emx1-Ctrl excluded)
    f.write(f"Emx1:\n")
    f.write(f"  NOTE: Emx1-Ctrl excluded (failed sample)\n")
    if "Emx1-Mut" in signal_matrices_sorted:
        mut_mean = np.nanmean(signal_matrices_sorted["Emx1-Mut"])
        f.write(f"  Emx1-Mut mean: {mut_mean:.4f}\n")
    f.write(f"  Ctrl vs Mut comparison: Not available\n\n")

    f.write("OUTPUT FILES:\n")
    f.write("-"*80 + "\n")
    f.write("Heatmaps:\n")
    f.write("  - heatmap_all_conditions.png\n")
    f.write("  - heatmap_cluster_*.png (by cluster)\n\n")
    f.write("Metaprofiles:\n")
    f.write("  - metaprofile_all_CREs.png\n")
    f.write("  - metaprofile_by_cluster.png\n\n")
    f.write("Comparisons:\n")
    f.write("  - comparison_nestin_ctrl_vs_mut.png\n")
    f.write("  NOTE: Emx1 comparison not available (Emx1-Ctrl excluded)\n\n")
    f.write("Data:\n")
    f.write("  - GABA_CREs_with_clusters.tsv\n")
    f.write("  - signal_matrix_*.npy (numpy arrays)\n\n")

print(f"✓ Saved: analysis_summary.txt")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print(f"\n{'='*80}")
print("ANALYSIS COMPLETE!")
print("="*80)

print(f"\nProcessed {len(cres_sorted)} GABA-specific CREs across 4 conditions:")
for condition_name in signal_matrices_sorted.keys():
    matrix = signal_matrices_sorted[condition_name]
    mean_sig = np.nanmean(matrix)
    print(f"  {condition_name}: mean signal = {mean_sig:.4f}")

print(f"\nOutput directory: {OUTPUT_DIR}/")
print(f"\nKey visualizations:")
print(f"  - heatmap_all_conditions.png: Side-by-side heatmaps")
print(f"  - metaprofile_all_CREs.png: Average profiles with SEM")
print(f"  - comparison_*_ctrl_vs_mut.png: Genotype-specific comparisons")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
