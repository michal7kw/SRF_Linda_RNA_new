#!/usr/bin/env python3
"""
Create Custom Comparison Visualizations for GABA DEG ENCODE cCREs

This script creates visualizations comparing ATAC signal at ENCODE cCREs
associated with GABA up-regulated and down-regulated DEGs.

COMPARISONS:
For each DEG set (up/down):
1. Nestin-Ctrl vs Nestin-Mut (within-genotype mutation effect)
2. Nestin-Ctrl vs Emx1-Mut (cross-genotype mutation effect)
3. Nestin-Mut vs Emx1-Mut (mutant genotype comparison)

NOTE: Emx1-Ctrl is excluded due to quality issues.

INPUT:
- deepTools matrix files from Step 2
- BED files with CRE coordinates
- TSV files with CRE-gene associations

OUTPUT:
- Publication-quality metaprofiles for each comparison
- Difference plots
- Statistical summaries
- Optional individual CRE plots

PERFORMANCE OPTIONS:
- --skip-individual: Skip individual CRE plots (saves time)
- --parallel N: Use N parallel processes for individual plots
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import argparse
from multiprocessing import Pool
from functools import partial
from datetime import datetime

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
    sample_names_all = header[1:]  # Skip "genes:N"

    # Get unique sample names
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

    for line in data_lines:
        if line.strip():
            parts = line.strip().split('\t')
            data_matrix.append([float(x) for x in parts])

    data_array = np.array(data_matrix)
    print(f"  Data shape: {data_array.shape}")

    # Create region names
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


def plot_comparison_metaprofile(data, bin_labels, output_file, title,
                                sample1_name, sample2_name, color1, color2,
                                n_cres):
    """
    Create metaprofile comparing two samples with difference plot.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10),
                                    gridspec_kw={'height_ratios': [3, 1]})

    # Top panel: Metaprofiles
    sample1_data = data[sample1_name]
    sample2_data = data[sample2_name]

    # Calculate mean across CREs
    sample1_mean = np.mean(sample1_data, axis=0)
    sample2_mean = np.mean(sample2_data, axis=0)

    # Calculate SEM across CREs
    sample1_sem = stats.sem(sample1_data, axis=0)
    sample2_sem = stats.sem(sample2_data, axis=0)

    # Plot mean lines
    ax1.plot(bin_labels, sample1_mean, label=sample1_name, color=color1,
            linewidth=2.5, alpha=0.9)
    ax1.plot(bin_labels, sample2_mean, label=sample2_name, color=color2,
            linewidth=2.5, alpha=0.9)

    # Add confidence intervals (mean +/- SEM)
    ax1.fill_between(bin_labels,
                     sample1_mean - sample1_sem,
                     sample1_mean + sample1_sem,
                     color=color1, alpha=0.2)
    ax1.fill_between(bin_labels,
                     sample2_mean - sample2_sem,
                     sample2_mean + sample2_sem,
                     color=color2, alpha=0.2)

    ax1.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5,
                label='CRE Center')
    ax1.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Mean ATAC Signal', fontsize=12, fontweight='bold')
    ax1.set_title(f'{title}\n(n={n_cres} ENCODE cCREs)', fontsize=14, fontweight='bold', pad=20)
    ax1.legend(loc='best', frameon=True, fontsize=11)
    ax1.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax1)

    # Bottom panel: Difference (sample2 - sample1)
    difference = sample2_mean - sample1_mean

    # Plot difference
    ax2.plot(bin_labels, difference, color='#6C464F', linewidth=2, alpha=0.9)
    ax2.fill_between(bin_labels, 0, difference, where=(difference > 0),
                     color=color2, alpha=0.3, label=f'{sample2_name} > {sample1_name}')
    ax2.fill_between(bin_labels, 0, difference, where=(difference < 0),
                     color=color1, alpha=0.3, label=f'{sample1_name} > {sample2_name}')
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1.5)
    ax2.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

    ax2.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
    ax2.set_ylabel(f'Delta ATAC Signal\n({sample2_name} - {sample1_name})',
                   fontsize=11, fontweight='bold')
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


def plot_combined_overview(data, bin_labels, output_file, title, n_cres, deg_type):
    """
    Create a combined overview plot showing all 3 conditions.
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    # Colors for the 3 conditions (excluding Emx1-Ctrl)
    colors = {
        'Nestin-Ctrl': '#2E86AB',  # Blue
        'Nestin-Mut': '#A23B72',   # Magenta
        'Emx1-Mut': '#C73E1D'      # Red
    }

    for sample_name in ['Nestin-Ctrl', 'Nestin-Mut', 'Emx1-Mut']:
        if sample_name in data:
            sample_data = data[sample_name]
            sample_mean = np.mean(sample_data, axis=0)
            sample_sem = stats.sem(sample_data, axis=0)

            ax.plot(bin_labels, sample_mean, label=sample_name,
                   color=colors[sample_name], linewidth=2.5, alpha=0.9)
            ax.fill_between(bin_labels,
                           sample_mean - sample_sem,
                           sample_mean + sample_sem,
                           color=colors[sample_name], alpha=0.15)

    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mean ATAC Signal', fontsize=12, fontweight='bold')

    # Title with DEG type indicator
    deg_label = "Up-Regulated" if deg_type == "up" else "Down-Regulated"
    ax.set_title(f'{title}\n{deg_label} GABA DEGs (n={n_cres} ENCODE cCREs)',
                fontsize=14, fontweight='bold', pad=20)
    ax.legend(loc='best', frameon=True, fontsize=11)
    ax.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()


def calculate_comparison_statistics(data, sample1_name, sample2_name):
    """
    Calculate statistical summary for comparison.
    """
    stats_dict = {}

    sample1_data = data[sample1_name]
    sample2_data = data[sample2_name]

    # Calculate mean signal across all CREs and bins
    sample1_mean = np.mean(sample1_data)
    sample2_mean = np.mean(sample2_data)

    stats_dict['sample1_mean'] = sample1_mean
    stats_dict['sample1_std'] = np.std(sample1_data)
    stats_dict['sample2_mean'] = sample2_mean
    stats_dict['sample2_std'] = np.std(sample2_data)
    stats_dict['fold_change'] = sample2_mean / sample1_mean if sample1_mean > 0 else np.nan
    stats_dict['difference'] = sample2_mean - sample1_mean

    # Calculate per-CRE means for paired t-test
    sample1_cre_means = np.mean(sample1_data, axis=1)
    sample2_cre_means = np.mean(sample2_data, axis=1)

    # Paired t-test (same CREs in both conditions)
    t_stat, p_val = stats.ttest_rel(sample2_cre_means, sample1_cre_means)
    stats_dict['t_statistic'] = t_stat
    stats_dict['p_value'] = p_val

    return stats_dict


def plot_up_vs_down_comparison(up_data, down_data, bin_labels, output_dir,
                                sample_name, color_up, color_down, n_up, n_down):
    """
    Create comparison plot of up-DEGs vs down-DEGs for a given condition.
    """
    fig, ax = plt.subplots(figsize=(10, 7))

    # Calculate means
    up_mean = np.mean(up_data[sample_name], axis=0)
    down_mean = np.mean(down_data[sample_name], axis=0)

    # Calculate SEM
    up_sem = stats.sem(up_data[sample_name], axis=0)
    down_sem = stats.sem(down_data[sample_name], axis=0)

    # Plot
    ax.plot(bin_labels, up_mean, label=f'Up-DEGs (n={n_up})',
           color=color_up, linewidth=2.5, alpha=0.9)
    ax.plot(bin_labels, down_mean, label=f'Down-DEGs (n={n_down})',
           color=color_down, linewidth=2.5, alpha=0.9)

    ax.fill_between(bin_labels, up_mean - up_sem, up_mean + up_sem,
                   color=color_up, alpha=0.2)
    ax.fill_between(bin_labels, down_mean - down_sem, down_mean + down_sem,
                   color=color_down, alpha=0.2)

    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mean ATAC Signal', fontsize=12, fontweight='bold')
    ax.set_title(f'ATAC Signal at ENCODE cCREs: Up vs Down DEGs\n{sample_name}',
                fontsize=14, fontweight='bold', pad=20)
    ax.legend(loc='best', frameon=True, fontsize=11)
    ax.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax)

    plt.tight_layout()
    output_file = f"{output_dir}/profiles/metaprofile_up_vs_down_{sample_name.replace('-', '_')}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()


def _plot_single_cre_comparison(args):
    """
    Helper function for parallel processing of individual CRE plots.
    """
    (i, cre_id, gene, chrom, start, end, data, bin_labels,
     sample1_name, sample2_name, color1, color2, comparison_name,
     output_dir, dpi, deg_type) = args

    fig, ax = plt.subplots(figsize=(8, 5))

    # Get data for both samples
    sample1_signal = data[sample1_name][i, :]
    sample2_signal = data[sample2_name][i, :]

    ax.plot(bin_labels, sample1_signal, label=sample1_name,
           color=color1, linewidth=2.5, alpha=0.9)
    ax.plot(bin_labels, sample2_signal, label=sample2_name,
           color=color2, linewidth=2.5, alpha=0.9)

    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Distance from CRE Center (bp)', fontsize=11, fontweight='bold')
    ax.set_ylabel('ATAC Signal', fontsize=11, fontweight='bold')

    deg_label = "Up-DEG" if deg_type == "up" else "Down-DEG"
    ax.set_title(f'{gene} ({deg_label}) - {cre_id}\n{chrom}:{start}-{end}\n{comparison_name}',
                fontsize=12, fontweight='bold')
    ax.legend(loc='best', frameon=True)
    ax.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax)

    plt.tight_layout()
    safe_gene = str(gene).replace('/', '_').replace(' ', '_')
    safe_cre = str(cre_id).replace('/', '_').replace(' ', '_')
    output_file = f"{output_dir}/individual_{deg_type}_{comparison_name}_{safe_gene}_{safe_cre}.png"
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()

    return output_file


def plot_individual_cres_comparison(data, regions, bin_labels, bed_file, tsv_file,
                                    output_dir, sample1_name, sample2_name,
                                    color1, color2, comparison_name, deg_type,
                                    n_processes=1, dpi=150, max_plots=50,
                                    min_signal=2.0, min_fc=2.0):
    """
    Create individual plots for each CRE showing comparison.
    """
    # Read BED file to get CRE IDs and coordinates
    bed_df = pd.read_csv(bed_file, sep='\t', header=None,
                         names=['chr', 'start', 'end', 'cre_id', 'score', 'strand'])

    # Read gene mapping from TSV
    gene_map = {}
    if os.path.exists(tsv_file):
        tsv_df = pd.read_csv(tsv_file, sep='\t')
        # Group by cCRE_id1 and get first gene
        if 'cCRE_id1' in tsv_df.columns and 'Gene' in tsv_df.columns:
            gene_map = tsv_df.groupby('cCRE_id1')['Gene'].first().to_dict()

    # Limit number of plots
    # n_cres = min(len(bed_df), max_plots, len(regions))
    # if len(bed_df) > max_plots:
    #     print(f"  NOTE: Limiting to {max_plots} individual plots (out of {len(bed_df)} CREs)")

    # Prepare arguments for parallel processing
    plot_args = []
    
    for i in range(len(bed_df)):
        # Check filtering criteria first
        sample1_signal = data[sample1_name][i, :]
        sample2_signal = data[sample2_name][i, :]
        
        mean1 = np.mean(sample1_signal)
        mean2 = np.mean(sample2_signal)
        max_signal = max(np.max(sample1_signal), np.max(sample2_signal))
        
        # Calculate Fold Change (handle zeros)
        if mean1 > 0.01:
            fc = mean2 / mean1
        elif mean2 > 0.01:
            fc = 100.0 # High FC if 1 is 0 and 2 is positive
        else:
            fc = 1.0 # Both zero
            
        # Check thresholds
        # 1. Minimum signal check (at least one condition must have signal)
        if max_signal < min_signal:
            continue
            
        # 2. Fold change check (must be significant change up OR down)
        # FC >= min_fc OR FC <= 1/min_fc
        if not (fc >= min_fc or fc <= (1.0/min_fc)):
            continue
            
        # Stop if we have enough plots
        if len(plot_args) >= max_plots:
            break
            
        cre_id = bed_df.iloc[i]['cre_id']
        gene = gene_map.get(cre_id, 'Unknown')
        chrom = bed_df.iloc[i]['chr']
        start = bed_df.iloc[i]['start']
        end = bed_df.iloc[i]['end']

        args = (i, cre_id, gene, chrom, start, end, data, bin_labels,
                sample1_name, sample2_name, color1, color2, comparison_name,
                output_dir, dpi, deg_type)
        plot_args.append(args)

    print(f"  Filtered: Found {len(plot_args)} significant plots (checked {i+1}/{len(bed_df)} CREs)")

    # Create plots (parallel or sequential)
    if n_processes > 1:
        print(f"  Using {n_processes} parallel processes...")
        with Pool(n_processes) as pool:
            output_files = pool.map(_plot_single_cre_comparison, plot_args)
        print(f"  Created {len(output_files)} plots in parallel")
    else:
        print(f"  Creating {len(plot_args)} plots sequentially...")
        output_files = []
        for args in plot_args:
            output_file = _plot_single_cre_comparison(args)
            output_files.append(output_file)
        print(f"  Created {len(output_files)} plots")


def process_deg_set(deg_type, matrix_file, bed_file, tsv_file, output_dir, args):
    """
    Process one DEG set (up or down).
    """
    deg_label = "Up-Regulated" if deg_type == "up" else "Down-Regulated"
    print(f"\n{'='*80}")
    print(f"PROCESSING {deg_label.upper()} DEGs")
    print(f"{'='*80}\n")

    # Read matrix
    data, regions, bin_labels = read_deeptools_matrix(matrix_file)

    if data is None:
        print(f"ERROR: Failed to read matrix file for {deg_label} DEGs")
        return None, None, None

    n_cres = len(regions)
    print(f"\nLoaded {n_cres} CREs from matrix\n")

    # Create combined overview
    print("-" * 40)
    print("Creating combined overview...")
    plot_combined_overview(
        data, bin_labels,
        f"{output_dir}/overview_all_conditions_{deg_type}_DEGs.png",
        "ATAC Signal at ENCODE cCREs",
        n_cres, deg_type
    )

    # Define comparisons (3 comparisons, excluding Emx1-Ctrl)
    comparisons = [
        {
            'name': 'nestin_ctrl_vs_mut',
            'sample1': 'Nestin-Ctrl',
            'sample2': 'Nestin-Mut',
            'color1': '#2E86AB',  # Blue
            'color2': '#A23B72',  # Magenta
            'title': f'Nestin: Ctrl vs Mut - {deg_label} DEGs',
            'description': 'Within-genotype mutation effect (Nestin)'
        },
        {
            'name': 'nestin_ctrl_vs_emx1_mut',
            'sample1': 'Nestin-Ctrl',
            'sample2': 'Emx1-Mut',
            'color1': '#2E86AB',  # Blue
            'color2': '#C73E1D',  # Red
            'title': f'Nestin-Ctrl vs Emx1-Mut - {deg_label} DEGs',
            'description': 'Cross-genotype mutation effect'
        },
        {
            'name': 'nestin_mut_vs_emx1_mut',
            'sample1': 'Nestin-Mut',
            'sample2': 'Emx1-Mut',
            'color1': '#A23B72',  # Magenta
            'color2': '#C73E1D',  # Red
            'title': f'Mutant Comparison: Nestin vs Emx1 - {deg_label} DEGs',
            'description': 'Genotype comparison (Mutants only)'
        }
    ]

    # Process each comparison
    all_stats = {}

    for comp in comparisons:
        # Check if both samples exist in data
        if comp['sample1'] not in data or comp['sample2'] not in data:
            print(f"\n  Skipping: {comp['description']} (samples not in matrix)")
            continue

        print(f"\n{'-'*40}")
        print(f"COMPARISON: {comp['description']}")
        print(f"({comp['sample1']} vs {comp['sample2']})")
        print("-" * 40)

        # Metaprofile with comparison
        print(f"\nCreating metaprofile...")
        plot_comparison_metaprofile(
            data, bin_labels,
            f"{output_dir}/profiles/metaprofile_{comp['name']}_{deg_type}_DEGs.png",
            comp['title'],
            comp['sample1'], comp['sample2'],
            comp['color1'], comp['color2'],
            n_cres
        )

        # Individual CRE plots (optional)
        if not args.skip_individual and os.path.exists(bed_file):
            print(f"\nCreating individual CRE plots ({comp['name']})...")
            plot_individual_cres_comparison(
                data, regions, bin_labels, bed_file, tsv_file,
                f"{output_dir}/profiles",
                comp['sample1'], comp['sample2'],
                comp['color1'], comp['color2'],
                comp['name'], deg_type,
                n_processes=args.parallel,
                dpi=args.individual_dpi,
                max_plots=args.max_individual,
                min_signal=args.min_signal,
                min_fc=args.min_fc
            )
        elif args.skip_individual:
            print(f"\n  Skipping individual CRE plots for {comp['name']} (fast mode)")

        # Statistics
        stats_result = calculate_comparison_statistics(data, comp['sample1'], comp['sample2'])
        all_stats[comp['name']] = {
            'comparison': comp,
            'stats': stats_result
        }

        print(f"\nStatistics:")
        print(f"  {comp['sample1']} mean: {stats_result['sample1_mean']:.6f}")
        print(f"  {comp['sample2']} mean: {stats_result['sample2_mean']:.6f}")
        print(f"  Fold change ({comp['sample2']}/{comp['sample1']}): {stats_result['fold_change']:.4f}")
        print(f"  Difference: {stats_result['difference']:.6f}")
        print(f"  T-statistic: {stats_result['t_statistic']:.4f}")
        print(f"  P-value: {stats_result['p_value']:.4e}")

    return data, n_cres, all_stats


def main():
    """Main execution function"""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Create custom comparison visualizations for GABA DEG ENCODE cCREs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Performance options:
  Fast mode (metaprofiles only):
    python %(prog)s --skip-individual

  Parallel processing:
    python %(prog)s --parallel 8

Examples:
  # Quick overview (metaprofiles only):
  python %(prog)s --skip-individual

  # Full analysis with parallelization:
  python %(prog)s --parallel 8
        """
    )

    parser.add_argument('--skip-individual', action='store_true',
                       help='Skip individual CRE plots (saves time)')
    parser.add_argument('--parallel', type=int, default=1, metavar='N',
                       help='Number of parallel processes for individual plots (default: 1)')
    parser.add_argument('--individual-dpi', type=int, default=150, metavar='N',
                       help='DPI for individual plots (metaprofiles always 300, default: 150)')
    parser.add_argument('--max-individual', type=int, default=50, metavar='N',
                       help='Maximum number of individual plots per comparison (default: 50)')
    parser.add_argument('--min-signal', type=float, default=1.0, metavar='F',
                       help='Minimum max signal required to plot (default: 1.0)')
    parser.add_argument('--min-fc', type=float, default=1.5, metavar='F',
                       help='Minimum fold change required to plot (default: 1.5)')

    args = parser.parse_args()

    print("=" * 80)
    print("CUSTOM COMPARISON VISUALIZATIONS FOR GABA DEG ENCODE cCREs")
    print("=" * 80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # Change to script directory
    os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_GABA_DEGs_encode")

    # Display performance settings
    if args.skip_individual:
        print("FAST MODE: Skipping individual CRE plots")
    else:
        print(f"FULL MODE: Creating individual plots")
        if args.parallel > 1:
            print(f"   -> Using {args.parallel} parallel processes")
        print(f"   -> Individual plot DPI: {args.individual_dpi}")
        print(f"   -> Max individual plots: {args.max_individual}")
        print(f"   -> Filtering: Min Signal >= {args.min_signal}, Min FC >= {args.min_fc}")
    print()

    # Define paths
    matrix_dir = "./output/heatmaps_deeptools"
    output_dir = f"./output/custom_comparisons_minSig{args.min_signal}_minFC{args.min_fc}"

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/profiles", exist_ok=True)

    # Process UP-DEGs
    up_matrix = f"{matrix_dir}/matrix_up_DEGs.gz"
    up_bed = "./output/encode_cCREs_up_DEGs.bed"
    up_tsv = "./output/encode_cCREs_up_DEGs.tsv"

    up_data, n_up, up_stats = None, 0, {}
    if os.path.exists(up_matrix.replace('.gz', '.tab')):
        up_data, n_up, up_stats = process_deg_set("up", up_matrix, up_bed, up_tsv, output_dir, args)
    else:
        print(f"\nWARNING: Up-DEG matrix not found: {up_matrix}")
        print("Run step 2 (compute_signal_matrices.sh) first.")

    # Process DOWN-DEGs
    down_matrix = f"{matrix_dir}/matrix_down_DEGs.gz"
    down_bed = "./output/encode_cCREs_down_DEGs.bed"
    down_tsv = "./output/encode_cCREs_down_DEGs.tsv"

    down_data, n_down, down_stats = None, 0, {}
    if os.path.exists(down_matrix.replace('.gz', '.tab')):
        down_data, n_down, down_stats = process_deg_set("down", down_matrix, down_bed, down_tsv, output_dir, args)
    else:
        print(f"\nWARNING: Down-DEG matrix not found: {down_matrix}")
        print("Run step 2 (compute_signal_matrices.sh) first.")

    # Create up vs down comparisons for each condition
    if up_data is not None and down_data is not None:
        print("\n" + "=" * 80)
        print("CREATING UP vs DOWN COMPARISONS")
        print("=" * 80)

        # Read bin_labels from matrix file
        with open(up_matrix.replace('.gz', '.tab'), 'r') as f:
            lines = f.readlines()
        params = lines[1].strip().split('\t')
        upstream = downstream = 2000
        for param in params:
            if 'upstream:' in param:
                upstream = int(param.split(':')[1])
            if 'downstream:' in param:
                downstream = int(param.split(':')[1])

        for sample_name in ['Nestin-Ctrl', 'Nestin-Mut', 'Emx1-Mut']:
            if sample_name in up_data and sample_name in down_data:
                n_bins = up_data[sample_name].shape[1]
                bin_labels = np.linspace(-upstream, downstream, n_bins)

                plot_up_vs_down_comparison(
                    up_data, down_data, bin_labels, output_dir,
                    sample_name, '#E63946', '#457B9D', n_up, n_down
                )

    # Save statistics to file
    print("\n" + "=" * 80)
    print("SAVING STATISTICS")
    print("=" * 80)

    stats_file = f"{output_dir}/comparison_statistics.txt"
    with open(stats_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CUSTOM COMPARISON STATISTICS - GABA DEG ENCODE cCREs\n")
        f.write("=" * 80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("UP-REGULATED DEGs\n")
        f.write("=" * 40 + "\n")
        f.write(f"Number of ENCODE cCREs: {n_up}\n\n")

        for name, result in up_stats.items():
            comp = result['comparison']
            s = result['stats']
            f.write(f"{comp['description'].upper()}\n")
            f.write(f"({comp['sample1']} vs {comp['sample2']})\n")
            f.write("-" * 40 + "\n")
            f.write(f"  {comp['sample1']} mean signal: {s['sample1_mean']:.6f} +/- {s['sample1_std']:.6f}\n")
            f.write(f"  {comp['sample2']} mean signal: {s['sample2_mean']:.6f} +/- {s['sample2_std']:.6f}\n")
            f.write(f"  Fold change: {s['fold_change']:.4f}\n")
            f.write(f"  Difference: {s['difference']:.6f}\n")
            f.write(f"  T-statistic: {s['t_statistic']:.4f}\n")
            f.write(f"  P-value: {s['p_value']:.4e}\n")
            significance = "***" if s['p_value'] < 0.001 else "**" if s['p_value'] < 0.01 else "*" if s['p_value'] < 0.05 else "ns"
            f.write(f"  Significance: {significance}\n\n")

        f.write("\nDOWN-REGULATED DEGs\n")
        f.write("=" * 40 + "\n")
        f.write(f"Number of ENCODE cCREs: {n_down}\n\n")

        for name, result in down_stats.items():
            comp = result['comparison']
            s = result['stats']
            f.write(f"{comp['description'].upper()}\n")
            f.write(f"({comp['sample1']} vs {comp['sample2']})\n")
            f.write("-" * 40 + "\n")
            f.write(f"  {comp['sample1']} mean signal: {s['sample1_mean']:.6f} +/- {s['sample1_std']:.6f}\n")
            f.write(f"  {comp['sample2']} mean signal: {s['sample2_mean']:.6f} +/- {s['sample2_std']:.6f}\n")
            f.write(f"  Fold change: {s['fold_change']:.4f}\n")
            f.write(f"  Difference: {s['difference']:.6f}\n")
            f.write(f"  T-statistic: {s['t_statistic']:.4f}\n")
            f.write(f"  P-value: {s['p_value']:.4e}\n")
            significance = "***" if s['p_value'] < 0.001 else "**" if s['p_value'] < 0.01 else "*" if s['p_value'] < 0.05 else "ns"
            f.write(f"  Significance: {significance}\n\n")

    print(f"  Saved: {stats_file}")

    # Summary
    print("\n" + "=" * 80)
    print("VISUALIZATION COMPLETE!")
    print("=" * 80)
    print(f"\nOutput directory: {output_dir}/")
    print(f"\nGenerated files:")

    if n_up > 0:
        print(f"\n  Up-DEGs ({n_up} ENCODE cCREs):")
        print(f"    - overview_all_conditions_up_DEGs.png")
        print(f"    - profiles/metaprofile_nestin_ctrl_vs_mut_up_DEGs.png")
        print(f"    - profiles/metaprofile_nestin_ctrl_vs_emx1_mut_up_DEGs.png")
        print(f"    - profiles/metaprofile_nestin_mut_vs_emx1_mut_up_DEGs.png")

    if n_down > 0:
        print(f"\n  Down-DEGs ({n_down} ENCODE cCREs):")
        print(f"    - overview_all_conditions_down_DEGs.png")
        print(f"    - profiles/metaprofile_nestin_ctrl_vs_mut_down_DEGs.png")
        print(f"    - profiles/metaprofile_nestin_ctrl_vs_emx1_mut_down_DEGs.png")
        print(f"    - profiles/metaprofile_nestin_mut_vs_emx1_mut_down_DEGs.png")

    if n_up > 0 and n_down > 0:
        print(f"\n  Up vs Down comparisons:")
        print(f"    - profiles/metaprofile_up_vs_down_Nestin_Ctrl.png")
        print(f"    - profiles/metaprofile_up_vs_down_Nestin_Mut.png")
        print(f"    - profiles/metaprofile_up_vs_down_Emx1_Mut.png")

    print(f"\n  Statistics:")
    print(f"    - comparison_statistics.txt")

    print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)


if __name__ == "__main__":
    main()
