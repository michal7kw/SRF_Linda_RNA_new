#!/usr/bin/env python3
"""
Create Custom Comparison Visualizations for Splicing Gene CREs

This script creates visualizations for three custom comparisons:
1. Nestin-Ctrl vs Nestin-Mut
2. Nestin-Ctrl vs Emx1-Mut
3. Nestin-Mut vs Emx1-Mut

INPUT:
- deepTools matrix files (matrix_nestin_ctrl_vs_mut.gz, etc.)
- BED file with CRE coordinates

OUTPUT:
- Publication-quality metaprofiles for each comparison
- Difference plots
- Individual CRE profiles (optional, can be skipped for faster execution)
- Statistical summaries

PERFORMANCE OPTIONS:
- --skip-individual: Skip individual CRE plots (saves 30-50 minutes)
- --parallel N: Use N parallel processes for individual plots
- --individual-dpi N: DPI for individual plots (default: 150, metaprofiles always 300)
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

    print(f"  Window: ±{upstream} bp, Bin size: {bin_size} bp")

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
                                sample1_name, sample2_name, color1, color2):
    """
    Create metaprofile comparing two samples.

    Parameters:
    -----------
    data : dict
        Dictionary with sample names as keys
    bin_labels : array
        Genomic positions
    output_file : str
        Output file path
    title : str
        Plot title
    sample1_name : str
        Name of first sample
    sample2_name : str
        Name of second sample
    color1 : str
        Color for sample 1
    color2 : str
        Color for sample 2
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

    # Add confidence intervals (mean ± SEM)
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
    ax1.set_title(title, fontsize=14, fontweight='bold', pad=20)
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
    ax2.set_ylabel(f'Δ ATAC Signal\n({sample2_name} - {sample1_name})',
                   fontsize=11, fontweight='bold')
    ax2.legend(loc='best', frameon=True, fontsize=10)
    ax2.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax2)

    # Add statistical annotation
    max_diff = np.max(np.abs(difference))
    mean_diff = np.mean(difference)
    ax2.text(0.02, 0.95, f'Max |Δ|: {max_diff:.4f}\nMean Δ: {mean_diff:.4f}',
            transform=ax2.transAxes, fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()

def _plot_single_cre_comparison(args):
    """
    Helper function for parallel processing of individual CRE plots.
    This must be a top-level function for multiprocessing to work.
    """
    (i, cre_id, gene, chrom, start, end, data, bin_labels,
     sample1_name, sample2_name, color1, color2, comparison_name,
     output_dir, dpi) = args

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
    ax.set_title(f'{gene} - {cre_id}\n{chrom}:{start}-{end}\n{comparison_name}',
                fontsize=12, fontweight='bold')
    ax.legend(loc='best', frameon=True)
    ax.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax)

    plt.tight_layout()
    output_file = f"{output_dir}/individual_{comparison_name}_{gene}_{cre_id}.png"
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()

    return output_file

def plot_individual_cres_comparison(data, regions, bin_labels, bed_file, output_dir,
                                    sample1_name, sample2_name, color1, color2,
                                    comparison_name, n_processes=1, dpi=150,
                                    min_signal=1.0, min_fc=1.5):
    """
    Create individual plots for each CRE showing comparison.

    Parameters:
    -----------
    n_processes : int
        Number of parallel processes to use (default: 1 = sequential)
    dpi : int
        DPI for individual plots (default: 150)
    min_signal : float
        Minimum max signal required to plot (default: 1.0)
    min_fc : float
        Minimum fold change required to plot (default: 1.5)
    """
    # Read BED file to get CRE IDs and coordinates
    bed_df = pd.read_csv(bed_file, sep='\t', header=None,
                         names=['chr', 'start', 'end', 'cre_id', 'score', 'strand'])

    # Read gene mapping from TSV
    tsv_file = bed_file.replace('.bed', '_celltypes.tsv')
    gene_map = {}
    if os.path.exists(tsv_file):
        tsv_df = pd.read_csv(tsv_file, sep='\t')
        gene_map = tsv_df.groupby('cCRE1')['Gene'].first().to_dict()

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
            
        cre_id = bed_df.iloc[i]['cre_id']
        gene = gene_map.get(cre_id, 'Unknown')
        chrom = bed_df.iloc[i]['chr']
        start = bed_df.iloc[i]['start']
        end = bed_df.iloc[i]['end']

        args = (i, cre_id, gene, chrom, start, end, data, bin_labels,
                sample1_name, sample2_name, color1, color2, comparison_name,
                output_dir, dpi)
        plot_args.append(args)

    print(f"  Filtered: Found {len(plot_args)} significant plots (checked {len(bed_df)} CREs)")

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
            print(f"  Saved: {output_file}")
            output_files.append(output_file)

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

def main():
    """Main execution function"""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Create custom comparison visualizations for splicing gene CREs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Performance options:
  Fast mode (metaprofiles only):
    python %(prog)s --skip-individual

  Parallel processing (8x faster):
    python %(prog)s --parallel 8

  Lower DPI (2x faster):
    python %(prog)s --individual-dpi 100

  Combined (fastest):
    python %(prog)s --skip-individual
    # OR for full analysis:
    python %(prog)s --parallel 8 --individual-dpi 100

Examples:
  # Quick overview (metaprofiles only):
  python %(prog)s --skip-individual

  # Full analysis with parallelization:
  python %(prog)s --parallel 8

  # Publication quality with parallelization:
  python %(prog)s --parallel 8 --individual-dpi 300
        """
    )

    parser.add_argument('--skip-individual', action='store_true',
                       help='Skip individual CRE plots (saves 30-50 minutes)')
    parser.add_argument('--parallel', type=int, default=1, metavar='N',
                       help='Number of parallel processes for individual plots (default: 1)')
    parser.add_argument('--individual-dpi', type=int, default=150, metavar='N',
                       help='DPI for individual plots (metaprofiles always 300, default: 150)')
    parser.add_argument('--min-signal', type=float, default=1.0, metavar='F',
                       help='Minimum max signal required to plot (default: 1.0)')
    parser.add_argument('--min-fc', type=float, default=1.5, metavar='F',
                       help='Minimum fold change required to plot (default: 1.5)')

    args = parser.parse_args()

    print("="*80)
    print("CUSTOM COMPARISON VISUALIZATIONS FOR SPLICING GENE CREs")
    print("="*80)
    print()

    # Display performance settings
    if args.skip_individual:
        print("⚡ FAST MODE: Skipping individual CRE plots")
    else:
        print(f"📊 FULL MODE: Creating individual plots")
        if args.parallel > 1:
            print(f"   → Using {args.parallel} parallel processes")
        print(f"   → Individual plot DPI: {args.individual_dpi}")
        print(f"   → Filtering: Min Signal >= {args.min_signal}, Min FC >= {args.min_fc}")
        print(f"   → Metaprofile DPI: 300 (always high quality)")
    print()

    # Define paths
    matrix_dir = "./output/custom_comparisons"
    output_dir = f"{matrix_dir}_minSig{args.min_signal}_minFC{args.min_fc}/profiles"
    bed_file = "./output/splicing_genes_CREs_all.bed"
    # NOTE: Using ALL CREs (not just GABA-specific) but with GABA BigWig files ONLY for signal analysis

    os.makedirs(output_dir, exist_ok=True)

    # Define comparisons
    comparisons = [
        {
            'name': 'comparison1',
            'matrix': 'matrix_nestin_ctrl_vs_mut',
            'sample1': 'Nestin-Ctrl',
            'sample2': 'Nestin-Mut',
            'color1': '#2E86AB',  # Blue
            'color2': '#C73E1D',  # Red
            'title': 'Nestin: Ctrl vs Mut\nATAC Signal at Splicing Gene CREs',
            'description': 'Within-genotype mutation effect'
        },
        {
            'name': 'comparison2',
            'matrix': 'matrix_nestin_ctrl_vs_emx1_mut',
            'sample1': 'Nestin-Ctrl',
            'sample2': 'Emx1-Mut',
            'color1': '#2E86AB',  # Blue
            'color2': '#FF8C42',  # Orange
            'title': 'Nestin-Ctrl vs Emx1-Mut\nATAC Signal at Splicing Gene CREs',
            'description': 'Cross-genotype comparison'
        },
        {
            'name': 'comparison3',
            'matrix': 'matrix_nestin_mut_vs_emx1_mut',
            'sample1': 'Nestin-Mut',
            'sample2': 'Emx1-Mut',
            'color1': '#C73E1D',  # Red
            'color2': '#FF8C42',  # Orange
            'title': 'Nestin-Mut vs Emx1-Mut\nATAC Signal at Splicing Gene CREs',
            'description': 'Mutant-to-mutant genotype effect'
        }
    ]

    # Process each comparison
    all_stats = {}

    for comp in comparisons:
        print("\n" + "="*80)
        print(f"COMPARISON: {comp['description'].upper()}")
        print(f"({comp['sample1']} vs {comp['sample2']})")
        print("="*80)

        matrix_file = f"{matrix_dir}/{comp['matrix']}.gz"

        if not os.path.exists(matrix_file.replace('.gz', '.tab')):
            print(f"  WARNING: Matrix file not found, skipping: {matrix_file}")
            continue

        # Read matrix
        data, regions, bin_labels = read_deeptools_matrix(matrix_file)

        if data is None:
            continue

        # Metaprofile with comparison
        print(f"\nCreating metaprofile...")
        plot_comparison_metaprofile(
            data, bin_labels,
            f"{output_dir}/metaprofile_{comp['sample1'].lower().replace('-', '_')}_vs_{comp['sample2'].lower().replace('-', '_')}.png",
            comp['title'],
            comp['sample1'], comp['sample2'],
            comp['color1'], comp['color2']
        )

        # Individual CRE plots (optional)
        if not args.skip_individual:
            print(f"\nCreating individual CRE plots ({comp['name']})...")
            plot_individual_cres_comparison(
                data, regions, bin_labels, bed_file,
                output_dir,
                comp['sample1'], comp['sample2'],
                comp['color1'], comp['color2'],
                comp['name'],
                n_processes=args.parallel,
                dpi=args.individual_dpi,
                min_signal=args.min_signal,
                min_fc=args.min_fc
            )
        else:
            print(f"\n⚡ Skipping individual CRE plots for {comp['name']} (fast mode)")
            print(f"   To create individual plots, run without --skip-individual")

        # Statistics
        stats_result = calculate_comparison_statistics(data, comp['sample1'], comp['sample2'])
        all_stats[comp['name']] = stats_result

        print(f"\n{comp['description']} Statistics:")
        print(f"  {comp['sample1']} mean: {stats_result['sample1_mean']:.6f}")
        print(f"  {comp['sample2']} mean: {stats_result['sample2_mean']:.6f}")
        print(f"  Fold change ({comp['sample2']}/{comp['sample1']}): {stats_result['fold_change']:.4f}")
        print(f"  Difference ({comp['sample2']} - {comp['sample1']}): {stats_result['difference']:.6f}")
        print(f"  T-statistic: {stats_result['t_statistic']:.4f}")
        print(f"  P-value: {stats_result['p_value']:.4e}")

    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "="*80)
    print("VISUALIZATION COMPLETE!")
    print("="*80)
    print(f"\nOutput directory: {output_dir}/")
    print("\nGenerated files:")
    print("  Metaprofiles (always created):")
    print("    - metaprofile_nestin_ctrl_vs_nestin_mut.png")
    print("    - metaprofile_nestin_ctrl_vs_emx1_mut.png")
    print("    - metaprofile_nestin_mut_vs_emx1_mut.png")

    if not args.skip_individual:
        # Get CRE count from BED file
        bed_df = pd.read_csv(bed_file, sep='\t', header=None,
                            names=['chr', 'start', 'end', 'cre_id', 'score', 'strand'])
        n_cres = len(bed_df)
        print(f"\n  Individual CRE plots (DPI: {args.individual_dpi}):")
        print(f"    - individual_comparison1_<gene>_<cre_id>.png ({n_cres} plots)")
        print(f"    - individual_comparison2_<gene>_<cre_id>.png ({n_cres} plots)")
        print(f"    - individual_comparison3_<gene>_<cre_id>.png ({n_cres} plots)")
        print(f"    Total: {n_cres * 3} individual plots")
        if args.parallel > 1:
            print(f"    (Created using {args.parallel} parallel processes)")
    else:
        print("\n  Individual CRE plots: SKIPPED (fast mode)")
        print("    To create individual plots, run without --skip-individual")

    # Save statistics to file
    stats_file = f"{matrix_dir}/comparison_statistics.txt"
    with open(stats_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("CUSTOM COMPARISON STATISTICS\n")
        f.write("="*80 + "\n\n")

        for comp in comparisons:
            if comp['name'] in all_stats:
                s = all_stats[comp['name']]
                f.write(f"{comp['description'].upper()}\n")
                f.write(f"({comp['sample1']} vs {comp['sample2']})\n")
                f.write("-"*80 + "\n")
                f.write(f"  {comp['sample1']} mean signal: {s['sample1_mean']:.6f} ± {s['sample1_std']:.6f}\n")
                f.write(f"  {comp['sample2']} mean signal: {s['sample2_mean']:.6f} ± {s['sample2_std']:.6f}\n")
                f.write(f"  Fold change: {s['fold_change']:.4f}\n")
                f.write(f"  Difference: {s['difference']:.6f}\n")
                f.write(f"  T-statistic: {s['t_statistic']:.4f}\n")
                f.write(f"  P-value: {s['p_value']:.4e}\n")
                f.write("\n")

    print(f"\n✓ Saved statistics: {stats_file}")
    print()

if __name__ == "__main__":
    main()
