#!/usr/bin/env python3
"""
Create Publication-Quality ATAC Signal Profiles for Splicing Gene ENCODE cCREs

This script creates enhanced visualizations of ATAC-seq signal at ENCODE cCREs linked
to splicing genes, with focus on Ctrl vs Mut comparisons.

PERFORMANCE OPTIONS:
- --skip-individual: Skip individual CRE plots (saves 30-50 minutes)
- --parallel N: Use N parallel processes for individual plots (default: 1)
- --individual-dpi N: DPI for individual plots (default: 150, metaprofiles always 300)

INPUT:
- deepTools matrix files (matrix_GABA.gz, matrix_GABA_nestin.gz, matrix_GABA_emx1.gz)
- BED file with CRE coordinates

OUTPUT:
- Publication-quality metaprofiles (always created, 300 DPI)
- Ctrl vs Mut comparison plots
- Difference plots (Mut - Ctrl)
- Individual CRE profiles (optional, configurable DPI)
- Statistical summary
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import gzip
import os
import argparse
from multiprocessing import Pool
from functools import partial

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

    Returns:
    --------
    data : dict
        Dictionary with sample names as keys and signal arrays as values
    regions : list
        List of region names
    bin_labels : array
        Array of genomic positions (bp from center)
    """
    print(f"Reading matrix: {matrix_file}")

    # Read the tabular version of the matrix
    matrix_tab = matrix_file.replace('.gz', '.tab')

    if not os.path.exists(matrix_tab):
        print(f"  ERROR: Matrix tab file not found: {matrix_tab}")
        return None, None, None

    # Read file to get parameters and header
    with open(matrix_tab, 'r') as f:
        lines = f.readlines()

    # Line 1: #genes:N
    n_regions = int(lines[0].strip().split(':')[1])
    print(f"  Regions: {n_regions}")

    # Line 2: #downstream:X upstream:Y ... bin size:Z
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

    # Line 3: Header with sample names
    header = lines[2].strip().split('\t')
    # First element is "genes:N", rest are sample names repeated for each bin
    sample_names_all = header[1:]  # Skip "genes:N"

    # Get unique sample names (they repeat for each bin)
    # Count consecutive occurrences to determine bins per sample
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

    # Add last sample
    if prev_name is not None:
        sample_names_unique.append(prev_name)
        temp_bins.append(bin_count)

    n_bins = temp_bins[0]
    print(f"  Samples: {sample_names_unique}")
    print(f"  Bins per sample: {n_bins}")

    # Read data (starts from line 4, index 3)
    data_lines = lines[3:]
    data_matrix = []
    regions = []

    for line in data_lines:
        if line.strip():
            parts = line.strip().split('\t')
            data_matrix.append([float(x) for x in parts])

    data_array = np.array(data_matrix)
    print(f"  Data shape: {data_array.shape}")

    # Create region names (just use indices since we don't have BED info in this format)
    regions = [f"Region_{i+1}" for i in range(len(data_array))]

    # Split signal into samples
    data = {}
    for i, sample in enumerate(sample_names_unique):
        start_col = i * n_bins
        end_col = (i + 1) * n_bins
        data[sample] = data_array[:, start_col:end_col]

    # Create bin labels
    bin_labels = np.linspace(-upstream, downstream, n_bins)

    return data, regions, bin_labels

def plot_metaprofile_comparison(data, bin_labels, output_file, title,
                                 ctrl_samples, mut_samples):
    """
    Create metaprofile comparing Ctrl vs Mut with confidence intervals.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10),
                                     gridspec_kw={'height_ratios': [3, 1]})

    colors = {'Ctrl': '#2E86AB', 'Mut': '#C73E1D'}

    # Top panel: Metaprofiles with SEM
    for label, samples in [('Ctrl', ctrl_samples), ('Mut', mut_samples)]:
        # Get data for all samples in this condition
        sample_data = [data[s] for s in samples if s in data]

        if len(sample_data) == 0:
            continue

        # Average across CREs for each sample
        sample_means = [np.mean(d, axis=0) for d in sample_data]

        # Calculate mean and SEM across samples
        mean_signal = np.mean(sample_means, axis=0)
        sem_signal = stats.sem(sample_means, axis=0)

        # Plot mean line
        ax1.plot(bin_labels, mean_signal, label=label, color=colors[label],
                linewidth=2.5, alpha=0.9)

        # Add confidence interval (mean +/- SEM)
        ax1.fill_between(bin_labels,
                         mean_signal - sem_signal,
                         mean_signal + sem_signal,
                         color=colors[label], alpha=0.2)

    ax1.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5,
                label='CRE Center')
    ax1.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Mean ATAC Signal', fontsize=12, fontweight='bold')
    ax1.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax1.legend(loc='best', frameon=True, fontsize=11)
    ax1.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax1)

    # Bottom panel: Difference (Mut - Ctrl)
    if len(ctrl_samples) > 0 and len(mut_samples) > 0:
        # Calculate difference
        ctrl_data = [data[s] for s in ctrl_samples if s in data]
        mut_data = [data[s] for s in mut_samples if s in data]

        ctrl_means = [np.mean(d, axis=0) for d in ctrl_data]
        mut_means = [np.mean(d, axis=0) for d in mut_data]

        ctrl_mean = np.mean(ctrl_means, axis=0)
        mut_mean = np.mean(mut_means, axis=0)

        difference = mut_mean - ctrl_mean

        # Plot difference
        ax2.plot(bin_labels, difference, color='#6C464F', linewidth=2, alpha=0.9)
        ax2.fill_between(bin_labels, 0, difference, where=(difference > 0),
                         color='#C73E1D', alpha=0.3, label='Mut > Ctrl')
        ax2.fill_between(bin_labels, 0, difference, where=(difference < 0),
                         color='#2E86AB', alpha=0.3, label='Ctrl > Mut')
        ax2.axhline(y=0, color='black', linestyle='-', linewidth=1.5)
        ax2.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

        ax2.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Delta ATAC Signal\n(Mut - Ctrl)', fontsize=11, fontweight='bold')
        ax2.legend(loc='best', frameon=True, fontsize=10)
        ax2.grid(True, alpha=0.3, linestyle=':')
        sns.despine(ax=ax2)

        # Add statistical annotation
        max_diff = np.max(np.abs(difference))
        mean_diff = np.mean(difference)
        ax2.text(0.02, 0.95, f'Max |Delta|: {max_diff:.4f}\nMean Delta: {mean_diff:.4f}',
                transform=ax2.transAxes, fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()

def _plot_single_cre(args):
    """
    Helper function to plot a single CRE (for parallel processing).
    """
    (i, cre_id, gene, chrom, start, end, data_dict, bin_labels,
     ctrl_samples, mut_samples, genotype, output_dir, dpi) = args

    fig, ax = plt.subplots(figsize=(8, 5))

    colors = {'Ctrl': '#2E86AB', 'Mut': '#C73E1D'}

    for label, samples in [('Ctrl', ctrl_samples), ('Mut', mut_samples)]:
        sample_data = [data_dict[s][i, :] for s in samples if s in data_dict]

        if len(sample_data) == 0:
            continue

        mean_signal = np.mean(sample_data, axis=0)
        sem_signal = stats.sem(sample_data, axis=0) if len(sample_data) > 1 else np.zeros_like(mean_signal)

        ax.plot(bin_labels, mean_signal, label=f'{label} (n={len(sample_data)})',
               color=colors[label], linewidth=2.5, alpha=0.9)
        ax.fill_between(bin_labels,
                       mean_signal - sem_signal,
                       mean_signal + sem_signal,
                       color=colors[label], alpha=0.2)

    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Distance from CRE Center (bp)', fontsize=11, fontweight='bold')
    ax.set_ylabel('ATAC Signal', fontsize=11, fontweight='bold')
    ax.set_title(f'{gene} - {cre_id}\n{chrom}:{start}-{end} ({genotype})',
                fontsize=12, fontweight='bold')
    ax.legend(loc='best', frameon=True)
    ax.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax)

    plt.tight_layout()
    output_file = f"{output_dir}/individual_{genotype}_{gene}_{cre_id}.png"
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()

    return output_file

def plot_individual_cres(data, regions, bin_labels, bed_file, output_dir,
                         ctrl_samples, mut_samples, genotype, n_processes=1, dpi=150,
                         min_signal=2.0, min_fc=2.0):
    """
    Create individual plots for each CRE showing Ctrl vs Mut.
    """
    # Read BED file to get CRE IDs and coordinates
    bed_df = pd.read_csv(bed_file, sep='\t', header=None,
                         names=['chr', 'start', 'end', 'cre_id', 'score', 'strand'])

    # Read gene mapping from TSV
    tsv_file = bed_file.replace('.bed', '.tsv')
    gene_map = {}
    if os.path.exists(tsv_file):
        tsv_df = pd.read_csv(tsv_file, sep='\t')
        gene_map = tsv_df.groupby('cCRE_id1')['Gene'].first().to_dict()

    # Prepare arguments for each CRE
    plot_args = []
    skipped_count = 0

    for i in range(len(bed_df)):
        # Check significance first
        max_val = 0
        mean_ctrl = 0
        mean_mut = 0
        n_ctrl = 0
        n_mut = 0

        # Calculate means for FC
        for s in ctrl_samples:
            if s in data:
                sig = data[s][i, :]
                max_val = max(max_val, np.max(sig))
                mean_ctrl += np.mean(sig)
                n_ctrl += 1

        for s in mut_samples:
            if s in data:
                sig = data[s][i, :]
                max_val = max(max_val, np.max(sig))
                mean_mut += np.mean(sig)
                n_mut += 1

        if n_ctrl > 0: mean_ctrl /= n_ctrl
        if n_mut > 0: mean_mut /= n_mut

        # Filter by min signal
        if max_val < min_signal:
            skipped_count += 1
            continue

        # Filter by FC
        if mean_ctrl < 0.01 and mean_mut < 0.01:
            fc = 1.0
        elif mean_ctrl < 0.01:
            fc = 100.0
        else:
            fc = mean_mut / mean_ctrl

        # Check if FC meets criteria (up or down)
        if not (fc >= min_fc or fc <= (1.0/min_fc)):
            skipped_count += 1
            continue

        cre_id = bed_df.iloc[i]['cre_id']
        gene = gene_map.get(cre_id, 'Unknown')
        chrom = bed_df.iloc[i]['chr']
        start = bed_df.iloc[i]['start']
        end = bed_df.iloc[i]['end']

        plot_args.append((i, cre_id, gene, chrom, start, end, data, bin_labels,
                         ctrl_samples, mut_samples, genotype, output_dir, dpi))

    print(f"  Filtered: Found {len(plot_args)} significant plots (checked {len(bed_df)} CREs)")

    # Create plots (parallel or sequential)
    if n_processes > 1:
        print(f"  Creating {len(plot_args)} plots using {n_processes} parallel processes...")
        with Pool(n_processes) as pool:
            output_files = pool.map(_plot_single_cre, plot_args)
        print(f"  Completed {len(output_files)} plots")
    else:
        print(f"  Creating {len(plot_args)} plots sequentially...")
        for i, args in enumerate(plot_args):
            output_file = _plot_single_cre(args)
            if (i + 1) % 50 == 0:
                print(f"    Progress: {i+1}/{len(plot_args)} plots")
        print(f"  Completed {len(plot_args)} plots")

def calculate_statistics(data, ctrl_samples, mut_samples):
    """
    Calculate statistical summary of Ctrl vs Mut differences.
    """
    stats_dict = {}

    # Get data for each condition
    ctrl_data = [data[s] for s in ctrl_samples if s in data]
    mut_data = [data[s] for s in mut_samples if s in data]

    if len(ctrl_data) == 0 or len(mut_data) == 0:
        return stats_dict

    # Calculate mean signal across all CREs and bins
    ctrl_means = [np.mean(d) for d in ctrl_data]
    mut_means = [np.mean(d) for d in mut_data]

    stats_dict['ctrl_mean'] = np.mean(ctrl_means)
    stats_dict['ctrl_std'] = np.std(ctrl_means)
    stats_dict['mut_mean'] = np.mean(mut_means)
    stats_dict['mut_std'] = np.std(mut_means)
    stats_dict['fold_change'] = stats_dict['mut_mean'] / stats_dict['ctrl_mean'] if stats_dict['ctrl_mean'] > 0 else np.nan
    stats_dict['difference'] = stats_dict['mut_mean'] - stats_dict['ctrl_mean']

    # T-test (if we have replicates)
    if len(ctrl_means) > 1 and len(mut_means) > 1:
        t_stat, p_val = stats.ttest_ind(ctrl_means, mut_means)
        stats_dict['t_statistic'] = t_stat
        stats_dict['p_value'] = p_val

    return stats_dict

def main():
    """Main execution function"""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Create ATAC signal profiles for splicing gene ENCODE cCREs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--skip-individual', action='store_true',
                       help='Skip individual CRE plots (saves 30-50 minutes)')
    parser.add_argument('--parallel', type=int, default=1, metavar='N',
                       help='Number of parallel processes for individual plots')
    parser.add_argument('--individual-dpi', type=int, default=150, metavar='N',
                       help='DPI for individual plots (metaprofiles always 300)')
    parser.add_argument('--min-signal', type=float, default=1.0, metavar='F',
                       help='Minimum max signal required to plot (default: 1.0)')
    parser.add_argument('--min-fc', type=float, default=1.5, metavar='F',
                       help='Minimum fold change required to plot (default: 1.5)')
    args = parser.parse_args()

    print("="*80)
    print("CREATE PUBLICATION-QUALITY ATAC SIGNAL PROFILES")
    print("ENCODE cCREs LINKED TO SPLICING GENES")
    print("="*80)
    print()
    print("Performance settings:")
    print(f"  Individual plots: {'SKIPPED' if args.skip_individual else 'ENABLED'}")
    if not args.skip_individual:
        print(f"  Parallel processes: {args.parallel}")
        print(f"  Individual plot DPI: {args.individual_dpi}")
    print()

    # Define paths
    matrix_dir = "./output/heatmaps_deeptools"
    output_dir = f"{matrix_dir}/profiles_minSig{args.min_signal}_minFC{args.min_fc}"
    bed_file = "./output/CREs_splicing_genes_encode_all.bed"

    os.makedirs(output_dir, exist_ok=True)

    # ========================================================================
    # 1. Nestin: Ctrl vs Mut
    # ========================================================================
    print("\n" + "="*80)
    print("NESTIN: Ctrl vs Mut Comparison")
    print("="*80)

    matrix_nestin = f"{matrix_dir}/matrix_GABA_nestin.gz"
    if os.path.exists(matrix_nestin.replace('.gz', '.tab')):
        data_nestin, regions, bin_labels = read_deeptools_matrix(matrix_nestin)

        if data_nestin is not None:
            ctrl_samples = [s for s in data_nestin.keys() if 'Ctrl' in s]
            mut_samples = [s for s in data_nestin.keys() if 'Mut' in s]

            print(f"\nCtrl samples: {ctrl_samples}")
            print(f"Mut samples: {mut_samples}")

            # Metaprofile with comparison
            plot_metaprofile_comparison(
                data_nestin, bin_labels,
                f"{output_dir}/metaprofile_nestin_ctrl_vs_mut.png",
                "Nestin: ATAC Signal at Splicing Gene ENCODE cCREs\nCtrl vs Mut Comparison",
                ctrl_samples, mut_samples
            )

            # Individual CRE plots
            if not args.skip_individual:
                print("\nCreating individual CRE plots (Nestin)...")
                plot_individual_cres(data_nestin, regions, bin_labels, bed_file,
                                    output_dir, ctrl_samples, mut_samples, 'Nestin',
                                    n_processes=args.parallel, dpi=args.individual_dpi,
                                    min_signal=args.min_signal, min_fc=args.min_fc)
            else:
                print("\nSkipping individual CRE plots (use --no-skip-individual to enable)...")

            # Statistics
            stats_nestin = calculate_statistics(data_nestin, ctrl_samples, mut_samples)
            print("\nNestin Statistics:")
            for key, val in stats_nestin.items():
                print(f"  {key}: {val:.6f}")

    # ========================================================================
    # 2. Emx1: Ctrl vs Mut
    # ========================================================================
    print("\n" + "="*80)
    print("EMX1: Ctrl vs Mut Comparison")
    print("="*80)

    matrix_emx1 = f"{matrix_dir}/matrix_GABA_emx1.gz"
    if os.path.exists(matrix_emx1.replace('.gz', '.tab')):
        data_emx1, regions, bin_labels = read_deeptools_matrix(matrix_emx1)

        if data_emx1 is not None:
            ctrl_samples = [s for s in data_emx1.keys() if 'Ctrl' in s]
            mut_samples = [s for s in data_emx1.keys() if 'Mut' in s]

            print(f"\nCtrl samples: {ctrl_samples}")
            print(f"Mut samples: {mut_samples}")

            # Metaprofile with comparison
            plot_metaprofile_comparison(
                data_emx1, bin_labels,
                f"{output_dir}/metaprofile_emx1_ctrl_vs_mut.png",
                "Emx1: ATAC Signal at Splicing Gene ENCODE cCREs\nNestin-Ctrl vs Emx1-Mut Comparison",
                ctrl_samples, mut_samples
            )

            # Individual CRE plots
            if not args.skip_individual:
                print("\nCreating individual CRE plots (Emx1)...")
                plot_individual_cres(data_emx1, regions, bin_labels, bed_file,
                                    output_dir, ctrl_samples, mut_samples, 'Emx1',
                                    n_processes=args.parallel, dpi=args.individual_dpi,
                                    min_signal=args.min_signal, min_fc=args.min_fc)
            else:
                print("\nSkipping individual CRE plots (use --no-skip-individual to enable)...")

            # Statistics
            stats_emx1 = calculate_statistics(data_emx1, ctrl_samples, mut_samples)
            print("\nEmx1 Statistics:")
            for key, val in stats_emx1.items():
                print(f"  {key}: {val:.6f}")

    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "="*80)
    print("VISUALIZATION COMPLETE!")
    print("="*80)
    print(f"\nOutput directory: {output_dir}/")
    print("\nGenerated files:")
    print("  Metaprofiles:")
    print("    - metaprofile_nestin_ctrl_vs_mut.png")
    print("    - metaprofile_emx1_ctrl_vs_mut.png")
    print("\n  Individual CRE plots:")
    print("    - individual_Nestin_<gene>_<cre_id>.png")
    print("    - individual_Emx1_<gene>_<cre_id>.png")
    print()

if __name__ == "__main__":
    main()
