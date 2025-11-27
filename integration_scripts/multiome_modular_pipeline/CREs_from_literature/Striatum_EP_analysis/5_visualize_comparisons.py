#!/usr/bin/env python3
"""
Visualize Custom Comparisons for Striatum_EP CREs

This script generates publication-quality metaprofiles and statistical comparisons
for the Striatum_EP analysis.

Usage:
    python 5_visualize_comparisons.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy import stats

# Configuration
BASE_DIR = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/Striatum_EP_analysis"
os.chdir(BASE_DIR)

MATRIX_FILE = "./output/heatmaps_deeptools/matrix_Striatum_EP.tab"
OUTPUT_DIR = "./output/custom_comparisons"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "profiles"), exist_ok=True)

# Note: Emx1-Ctrl is a failed sample - use Nestin-Ctrl as reference for all comparisons
# Matrix still contains all 4 samples in order: Nestin-Ctrl, Nestin-Mut, Emx1-Ctrl, Emx1-Mut
MATRIX_LABELS = ["Nestin-Ctrl", "Nestin-Mut", "Emx1-Ctrl", "Emx1-Mut"]

# For visualization, exclude failed Emx1-Ctrl sample
LABELS = ["Nestin-Ctrl", "Nestin-Mut", "Emx1-Mut"]
COLORS = {'Nestin-Ctrl': '#808080', 'Nestin-Mut': '#FF0000', 'Emx1-Mut': '#0000FF'}
MATRIX_INDEX = {'Nestin-Ctrl': 0, 'Nestin-Mut': 1, 'Emx1-Ctrl': 2, 'Emx1-Mut': 3}

def load_deeptools_matrix(matrix_file):
    """Load deepTools matrix (tab-delimited).

    deepTools computeMatrix output format (.tab):
    - Line 1: #genes:N (comment with region count)
    - Line 2: #downstream:X upstream:Y ... (parameters)
    - Line 3: genes:N sample1 sample1 ... (column labels, repeated per bin)
    - Line 4+: Data rows (numeric values only, no BED columns in .tab format)
    """
    print(f"Loading matrix: {matrix_file}")

    # Count header lines (lines starting with # or non-numeric content)
    skip_rows = 0
    sample_labels = []
    n_regions = 0

    with open(matrix_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                # Parse region count from first line
                if 'genes:' in line:
                    try:
                        n_regions = int(line.split(':')[1])
                    except (ValueError, IndexError):
                        pass
                skip_rows += 1
            elif line.startswith('genes:'):
                # This is the sample label row
                parts = line.split('\t')
                sample_labels = parts[1:]  # Skip 'genes:N' prefix
                skip_rows += 1
                break
            else:
                # Reached data rows
                break

    print(f"Skipping {skip_rows} header rows")
    if sample_labels:
        unique_samples = list(dict.fromkeys(sample_labels))  # Preserve order, remove duplicates
        print(f"Samples found: {unique_samples}")

    df = pd.read_csv(matrix_file, sep='\t', skiprows=skip_rows, header=None)

    # In .tab format from computeMatrix, there are no BED columns - just signal values
    # The data is all numeric (signal values per bin)
    print(f"Loaded matrix with shape: {df.shape}")

    # Return None for metadata (no BED columns in .tab format) and the data matrix
    return None, df

def plot_metaprofile(data, title, filename):
    """Plot metaprofile for selected conditions (excludes failed Emx1-Ctrl)."""
    plt.figure(figsize=(10, 6))

    # Calculate mean profile for each condition
    # Data is concatenated: [Cond1_Bin1...Cond1_BinN, Cond2_Bin1...Cond2_BinN, ...]
    # Matrix has 4 samples, but we only plot 3 (excluding Emx1-Ctrl)
    n_bins = data.shape[1] // 4
    x = np.linspace(-2000, 2000, n_bins)

    for label in LABELS:
        idx = MATRIX_INDEX[label]
        start = idx * n_bins
        end = (idx + 1) * n_bins
        profile = data.iloc[:, start:end].mean(axis=0)

        # Standard Error
        sem = data.iloc[:, start:end].sem(axis=0)

        color = COLORS[label]
        linestyle = '-' if 'Ctrl' in label else '--'

        plt.plot(x, profile, label=label, color=color, linestyle=linestyle)
        plt.fill_between(x, profile-sem, profile+sem, color=color, alpha=0.2)

    plt.xlabel("Distance from Center (bp)")
    plt.ylabel("Mean ATAC Signal")
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"Saved: {filename}")


def plot_pairwise_comparison(data, label1, label2, title, filename):
    """Plot pairwise comparison profile between two conditions."""
    plt.figure(figsize=(10, 6))

    n_bins = data.shape[1] // 4
    x = np.linspace(-2000, 2000, n_bins)

    for label in [label1, label2]:
        idx = MATRIX_INDEX[label]
        start = idx * n_bins
        end = (idx + 1) * n_bins
        profile = data.iloc[:, start:end].mean(axis=0)
        sem = data.iloc[:, start:end].sem(axis=0)

        color = COLORS[label]
        linestyle = '-' if 'Ctrl' in label else '--'

        plt.plot(x, profile, label=label, color=color, linestyle=linestyle, linewidth=2)
        plt.fill_between(x, profile-sem, profile+sem, color=color, alpha=0.2)

    plt.xlabel("Distance from Center (bp)")
    plt.ylabel("Mean ATAC Signal")
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"Saved: {filename}")


def main():
    if not os.path.exists(MATRIX_FILE):
        print(f"Matrix file not found: {MATRIX_FILE}")
        print("Please run step 4 first.")
        return

    _, data = load_deeptools_matrix(MATRIX_FILE)
    print(f"Loaded {len(data)} regions.")
    
    # Plot Overall Metaprofile
    plot_metaprofile(data, "Striatum_EP CREs (Splicing Genes)", os.path.join(OUTPUT_DIR, "overview_all_conditions.png"))

    # Plot Pairwise Comparisons (saved to profiles/ subdirectory)
    profiles_dir = os.path.join(OUTPUT_DIR, "profiles")

    # Nestin-Ctrl vs Nestin-Mut
    plot_pairwise_comparison(
        data, "Nestin-Ctrl", "Nestin-Mut",
        "Nestin: Ctrl vs Mut",
        os.path.join(profiles_dir, "profile_Nestin_Ctrl_vs_Mut.png")
    )

    # Nestin-Ctrl vs Emx1-Mut
    plot_pairwise_comparison(
        data, "Nestin-Ctrl", "Emx1-Mut",
        "Nestin-Ctrl vs Emx1-Mut",
        os.path.join(profiles_dir, "profile_Nestin_Ctrl_vs_Emx1_Mut.png")
    )

    # Nestin-Mut vs Emx1-Mut
    plot_pairwise_comparison(
        data, "Nestin-Mut", "Emx1-Mut",
        "Nestin-Mut vs Emx1-Mut",
        os.path.join(profiles_dir, "profile_Nestin_Mut_vs_Emx1_Mut.png")
    )

    # Calculate Statistics (Total Signal)
    # Use Nestin-Ctrl as reference for both mutants (Emx1-Ctrl is failed sample)
    n_bins = data.shape[1] // 4
    stats_data = {}

    for label in LABELS:
        idx = MATRIX_INDEX[label]
        start = idx * n_bins
        end = (idx + 1) * n_bins
        # Sum signal per region
        total_signal = data.iloc[:, start:end].sum(axis=1)
        stats_data[label] = total_signal

    stats_df = pd.DataFrame(stats_data)

    # T-tests using Nestin-Ctrl as reference for all comparisons
    with open(os.path.join(OUTPUT_DIR, "comparison_statistics.txt"), 'w') as f:
        f.write("Comparison Statistics (Paired T-test)\n")
        f.write("=====================================\n")
        f.write("Note: Nestin-Ctrl used as reference (Emx1-Ctrl is failed sample)\n\n")

        # Nestin-Ctrl vs Nestin-Mut
        t, p = stats.ttest_rel(stats_df['Nestin-Ctrl'], stats_df['Nestin-Mut'])
        f.write(f"Nestin-Ctrl vs Nestin-Mut:\n")
        f.write(f"  Mean Ctrl:       {stats_df['Nestin-Ctrl'].mean():.2f}\n")
        f.write(f"  Mean Nestin-Mut: {stats_df['Nestin-Mut'].mean():.2f}\n")
        f.write(f"  T-statistic:     {t:.4f}\n")
        f.write(f"  P-value:         {p:.2e}\n\n")

        # Nestin-Ctrl vs Emx1-Mut
        t, p = stats.ttest_rel(stats_df['Nestin-Ctrl'], stats_df['Emx1-Mut'])
        f.write(f"Nestin-Ctrl vs Emx1-Mut:\n")
        f.write(f"  Mean Ctrl:      {stats_df['Nestin-Ctrl'].mean():.2f}\n")
        f.write(f"  Mean Emx1-Mut:  {stats_df['Emx1-Mut'].mean():.2f}\n")
        f.write(f"  T-statistic:    {t:.4f}\n")
        f.write(f"  P-value:        {p:.2e}\n\n")

        # Also compare the two mutants
        t, p = stats.ttest_rel(stats_df['Nestin-Mut'], stats_df['Emx1-Mut'])
        f.write(f"Nestin-Mut vs Emx1-Mut:\n")
        f.write(f"  Mean Nestin-Mut: {stats_df['Nestin-Mut'].mean():.2f}\n")
        f.write(f"  Mean Emx1-Mut:   {stats_df['Emx1-Mut'].mean():.2f}\n")
        f.write(f"  T-statistic:     {t:.4f}\n")
        f.write(f"  P-value:         {p:.2e}\n")

    print("Statistics saved.")

if __name__ == "__main__":
    main()
