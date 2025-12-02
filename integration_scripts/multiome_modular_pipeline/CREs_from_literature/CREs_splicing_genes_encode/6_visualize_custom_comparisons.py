#!/usr/bin/env python3
"""
Custom ATAC Signal Comparisons for Splicing Gene ENCODE cCREs

This script creates custom comparisons across genotypes and conditions:
1. Nestin-Ctrl vs Nestin-Mut
2. Nestin-Ctrl vs Emx1-Mut
3. Nestin-Mut vs Emx1-Mut

NOTE: Emx1-Ctrl is excluded (failed sample)

PERFORMANCE OPTIONS:
- --skip-individual: Skip individual CRE plots
- --parallel N: Use N parallel processes for individual plots
- --individual-dpi N: DPI for individual plots (default: 150)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import argparse

# Set publication-quality plotting defaults
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 11
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
sns.set_style("ticks")

def read_deeptools_matrix(matrix_file):
    """
    Read deepTools matrix file and extract signal data.
    """
    print(f"Reading matrix: {matrix_file}")

    matrix_tab = matrix_file.replace('.gz', '.tab')

    if not os.path.exists(matrix_tab):
        print(f"  ERROR: Matrix tab file not found: {matrix_tab}")
        return None, None, None

    with open(matrix_tab, 'r') as f:
        lines = f.readlines()

    n_regions = int(lines[0].strip().split(':')[1])
    print(f"  Regions: {n_regions}")

    params = lines[1].strip().split('\t')
    bin_size = None
    downstream = None
    upstream = None
    for param in params:
        if 'bin size:' in param:
            bin_size = int(param.split(':')[1])
        elif 'downstream:' in param:
            downstream = int(param.split(':')[1])
        elif 'upstream:' in param:
            upstream = int(param.split(':')[1])

    print(f"  Window: +/-{upstream} bp, Bin size: {bin_size} bp")

    header = lines[2].strip().split('\t')
    sample_names_all = header[1:]

    sample_names_unique = []
    prev_name = None
    bin_count = 0
    temp_bins = []

    for name in sample_names_all:
        if name != prev_name:
            if prev_name is not None:
                sample_names_unique.append(prev_name)
                temp_bins.append(bin_count)
            prev_name = name
            bin_count = 1
        else:
            bin_count += 1

    if prev_name is not None:
        sample_names_unique.append(prev_name)
        temp_bins.append(bin_count)

    n_bins = temp_bins[0]
    print(f"  Samples: {sample_names_unique}")
    print(f"  Bins per sample: {n_bins}")

    data_lines = lines[3:]
    data_matrix = []

    for line in data_lines:
        if line.strip():
            parts = line.strip().split('\t')
            data_matrix.append([float(x) for x in parts])

    data_array = np.array(data_matrix)
    print(f"  Data shape: {data_array.shape}")

    regions = [f"Region_{i+1}" for i in range(len(data_array))]

    data = {}
    for i, sample in enumerate(sample_names_unique):
        start_col = i * n_bins
        end_col = (i + 1) * n_bins
        data[sample] = data_array[:, start_col:end_col]

    bin_labels = np.linspace(-upstream, downstream, n_bins)

    return data, regions, bin_labels

def plot_comparison(data, bin_labels, output_file, title, sample1, sample2, label1, label2, color1, color2):
    """
    Create comparison plot between two samples.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10),
                                     gridspec_kw={'height_ratios': [3, 1]})

    # Top panel: Metaprofiles with SEM
    for label, sample, color in [(label1, sample1, color1), (label2, sample2, color2)]:
        if sample not in data:
            continue

        signal = data[sample]
        mean_signal = np.mean(signal, axis=0)
        sem_signal = stats.sem(signal, axis=0)

        ax1.plot(bin_labels, mean_signal, label=f'{label} (n={len(signal)} CREs)',
                color=color, linewidth=2.5, alpha=0.9)
        ax1.fill_between(bin_labels,
                         mean_signal - sem_signal,
                         mean_signal + sem_signal,
                         color=color, alpha=0.2)

    ax1.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5,
                label='CRE Center')
    ax1.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Mean ATAC Signal', fontsize=12, fontweight='bold')
    ax1.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax1.legend(loc='best', frameon=True, fontsize=11)
    ax1.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax1)

    # Bottom panel: Difference
    if sample1 in data and sample2 in data:
        mean1 = np.mean(data[sample1], axis=0)
        mean2 = np.mean(data[sample2], axis=0)

        difference = mean2 - mean1

        ax2.plot(bin_labels, difference, color='#6C464F', linewidth=2, alpha=0.9)
        ax2.fill_between(bin_labels, 0, difference, where=(difference > 0),
                         color=color2, alpha=0.3, label=f'{label2} > {label1}')
        ax2.fill_between(bin_labels, 0, difference, where=(difference < 0),
                         color=color1, alpha=0.3, label=f'{label1} > {label2}')
        ax2.axhline(y=0, color='black', linestyle='-', linewidth=1.5)
        ax2.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

        ax2.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
        ax2.set_ylabel(f'Delta ATAC Signal\n({label2} - {label1})', fontsize=11, fontweight='bold')
        ax2.legend(loc='best', frameon=True, fontsize=10)
        ax2.grid(True, alpha=0.3, linestyle=':')
        sns.despine(ax=ax2)

        max_diff = np.max(np.abs(difference))
        mean_diff = np.mean(difference)
        fold_change = np.mean(mean2) / np.mean(mean1) if np.mean(mean1) > 0 else np.nan

        stats_text = f'Max |Delta|: {max_diff:.4f}\nMean Delta: {mean_diff:.4f}\nFold Change: {fold_change:.2f}x'
        ax2.text(0.02, 0.95, stats_text,
                transform=ax2.transAxes, fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()

def main():
    """Main execution function"""
    parser = argparse.ArgumentParser(
        description='Create custom comparisons for splicing gene ENCODE cCREs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--skip-individual', action='store_true',
                       help='Skip individual CRE plots')
    parser.add_argument('--parallel', type=int, default=1, metavar='N',
                       help='Number of parallel processes')
    parser.add_argument('--individual-dpi', type=int, default=150, metavar='N',
                       help='DPI for individual plots')
    parser.add_argument('--min-signal', type=float, default=1.0, metavar='F',
                       help='Minimum max signal required to plot')
    parser.add_argument('--min-fc', type=float, default=1.5, metavar='F',
                       help='Minimum fold change required to plot')

    args = parser.parse_args()

    print("="*80)
    print("CUSTOM COMPARISONS: SPLICING GENE ENCODE cCREs")
    print("="*80)
    print()

    os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_genes_encode")

    matrix_dir = "./output/custom_comparisons"
    output_dir = f"{matrix_dir}/profiles"
    os.makedirs(output_dir, exist_ok=True)

    # ========================================================================
    # Comparison 1: Nestin-Ctrl vs Nestin-Mut
    # ========================================================================
    print("\n" + "="*80)
    print("COMPARISON 1: Nestin-Ctrl vs Nestin-Mut")
    print("="*80)

    matrix_file = f"{matrix_dir}/matrix_nestin_ctrl_vs_mut.gz"
    if os.path.exists(matrix_file.replace('.gz', '.tab')):
        data, regions, bin_labels = read_deeptools_matrix(matrix_file)

        if data is not None:
            plot_comparison(
                data, bin_labels,
                f"{output_dir}/metaprofile_nestin_ctrl_vs_nestin_mut.png",
                "Nestin: Ctrl vs Mut at Splicing Gene ENCODE cCREs",
                "Nestin-Ctrl", "Nestin-Mut",
                "Nestin-Ctrl", "Nestin-Mut",
                '#2E86AB', '#C73E1D'
            )
    else:
        print(f"  Matrix not found: {matrix_file}")

    # ========================================================================
    # Comparison 2: Nestin-Ctrl vs Emx1-Mut
    # ========================================================================
    print("\n" + "="*80)
    print("COMPARISON 2: Nestin-Ctrl vs Emx1-Mut")
    print("="*80)

    matrix_file = f"{matrix_dir}/matrix_nestin_ctrl_vs_emx1_mut.gz"
    if os.path.exists(matrix_file.replace('.gz', '.tab')):
        data, regions, bin_labels = read_deeptools_matrix(matrix_file)

        if data is not None:
            plot_comparison(
                data, bin_labels,
                f"{output_dir}/metaprofile_nestin_ctrl_vs_emx1_mut.png",
                "Cross-Genotype: Nestin-Ctrl vs Emx1-Mut at Splicing Gene ENCODE cCREs",
                "Nestin-Ctrl", "Emx1-Mut",
                "Nestin-Ctrl", "Emx1-Mut",
                '#2E86AB', '#F18F01'
            )
    else:
        print(f"  Matrix not found: {matrix_file}")

    # ========================================================================
    # Comparison 3: Nestin-Mut vs Emx1-Mut
    # ========================================================================
    print("\n" + "="*80)
    print("COMPARISON 3: Nestin-Mut vs Emx1-Mut")
    print("="*80)

    matrix_file = f"{matrix_dir}/matrix_nestin_mut_vs_emx1_mut.gz"
    if os.path.exists(matrix_file.replace('.gz', '.tab')):
        data, regions, bin_labels = read_deeptools_matrix(matrix_file)

        if data is not None:
            plot_comparison(
                data, bin_labels,
                f"{output_dir}/metaprofile_nestin_mut_vs_emx1_mut.png",
                "Mutant Comparison: Nestin-Mut vs Emx1-Mut at Splicing Gene ENCODE cCREs",
                "Nestin-Mut", "Emx1-Mut",
                "Nestin-Mut", "Emx1-Mut",
                '#A23B72', '#F18F01'
            )
    else:
        print(f"  Matrix not found: {matrix_file}")

    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "="*80)
    print("VISUALIZATION COMPLETE!")
    print("="*80)
    print(f"\nOutput directory: {output_dir}/")
    print("\nGenerated files:")
    print("  - metaprofile_nestin_ctrl_vs_nestin_mut.png")
    print("  - metaprofile_nestin_ctrl_vs_emx1_mut.png")
    print("  - metaprofile_nestin_mut_vs_emx1_mut.png")
    print()

if __name__ == "__main__":
    main()
