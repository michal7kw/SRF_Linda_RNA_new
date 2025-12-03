#!/usr/bin/env python3
"""
Visualize Differentially Accessible CREs with minSig/minFC Filtering

This script creates visualizations for GABA-specific ENCODE cCREs that show
significant differential accessibility between conditions. Includes:

1. Metaprofiles for each comparison (Ctrl vs Mut, Nestin vs Emx1)
2. Individual CRE plots for significant CREs (based on minSig/minFC thresholds)
3. Statistical summaries

Filtering Criteria:
- min_signal: Minimum max signal required across conditions (default: 1.0)
- min_fc: Minimum fold change required (FC >= min_fc OR FC <= 1/min_fc, default: 1.5)

Usage:
    python 4_visualize_DA_CREs.py
    python 4_visualize_DA_CREs.py --min-signal 2.0 --min-fc 2.0
    python 4_visualize_DA_CREs.py --skip-individual
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
        List of region names (CRE IDs)
    bin_labels : array
        Array of genomic positions (bp from center)
    """
    print(f"Reading matrix: {matrix_file}")

    matrix_tab = matrix_file.replace('.gz', '.tab')

    if not os.path.exists(matrix_tab):
        print(f"  ERROR: Matrix tab file not found: {matrix_tab}")
        return None, None, None

    with open(matrix_tab, 'r') as f:
        lines = f.readlines()

    # Line 1: #genes:N
    n_regions = int(lines[0].strip().split(':')[1])
    print(f"  Regions: {n_regions}")

    # Line 2: Parameters
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
    sample_names_all = header[1:]

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

    if prev_name is not None:
        sample_names_unique.append(prev_name)
        temp_bins.append(bin_count)

    n_bins = temp_bins[0]
    print(f"  Samples: {sample_names_unique}")
    print(f"  Bins per sample: {n_bins}")

    # Read data
    data_lines = lines[3:]
    data_matrix = []

    for line in data_lines:
        if line.strip():
            parts = line.strip().split('\t')
            data_matrix.append([float(x) for x in parts])

    data_array = np.array(data_matrix)
    print(f"  Data shape: {data_array.shape}")

    regions = [f"CRE_{i+1}" for i in range(len(data_array))]

    # Split signal into samples
    data = {}
    for i, sample in enumerate(sample_names_unique):
        start_col = i * n_bins
        end_col = (i + 1) * n_bins
        data[sample] = data_array[:, start_col:end_col]

    bin_labels = np.linspace(-upstream, downstream, n_bins)

    return data, regions, bin_labels


def load_cre_info(tsv_file, bed_file):
    """Load CRE information from TSV and BED files."""
    cre_info = {}

    # Read BED file for coordinates
    if os.path.exists(bed_file):
        bed_df = pd.read_csv(bed_file, sep='\t', header=None,
                             names=['chr', 'start', 'end', 'cre_id', 'score', 'strand'])
        for idx, row in bed_df.iterrows():
            cre_info[idx] = {
                'cre_id': row['cre_id'],
                'chr': row['chr'],
                'start': row['start'],
                'end': row['end']
            }

    # Read TSV file for gene associations
    if os.path.exists(tsv_file):
        tsv_df = pd.read_csv(tsv_file, sep='\t')
        gene_map = tsv_df.groupby('cCRE_id1')['Gene'].first().to_dict()

        for idx in cre_info:
            cre_id = cre_info[idx]['cre_id']
            cre_info[idx]['gene'] = gene_map.get(cre_id, 'Unknown')

    return cre_info


def filter_significant_cres(data, regions, cre_info, sample1_name, sample2_name,
                            min_signal=2.0, min_fc=2.0):
    """
    Filter CREs based on signal and fold change thresholds.

    Parameters:
    -----------
    min_signal : float
        Minimum max signal required across conditions
    min_fc : float
        Minimum fold change (FC >= min_fc OR FC <= 1/min_fc)

    Returns:
    --------
    significant_cres : list
        List of (index, cre_info, fc, max_signal) for significant CREs
    """
    sample1_data = data[sample1_name]
    sample2_data = data[sample2_name]

    significant_cres = []

    for i in range(len(regions)):
        # Calculate mean signal for each sample
        mean1 = np.mean(sample1_data[i, :])
        mean2 = np.mean(sample2_data[i, :])

        # Calculate max signal across both conditions
        max_signal = max(np.max(sample1_data[i, :]), np.max(sample2_data[i, :]))

        # Check min signal threshold
        if max_signal < min_signal:
            continue

        # Calculate fold change (add small epsilon to avoid division by zero)
        eps = 1e-6
        fc = (mean2 + eps) / (mean1 + eps)

        # Check fold change threshold (up OR down)
        if not (fc >= min_fc or fc <= (1.0 / min_fc)):
            continue

        # Get CRE info
        info = cre_info.get(i, {
            'cre_id': f'CRE_{i+1}',
            'chr': 'unknown',
            'start': 0,
            'end': 0,
            'gene': 'Unknown'
        })

        significant_cres.append((i, info, fc, max_signal, mean1, mean2))

    return significant_cres


def _plot_single_cre(args):
    """Helper function for parallel processing of individual CRE plots."""
    (i, info, fc, max_signal, mean1, mean2, data, bin_labels,
     sample1_name, sample2_name, color1, color2, comparison_name,
     output_dir, dpi) = args

    fig, ax = plt.subplots(figsize=(8, 5))

    sample1_signal = data[sample1_name][i, :]
    sample2_signal = data[sample2_name][i, :]

    ax.plot(bin_labels, sample1_signal, label=sample1_name,
           color=color1, linewidth=2.5, alpha=0.9)
    ax.plot(bin_labels, sample2_signal, label=sample2_name,
           color=color2, linewidth=2.5, alpha=0.9)

    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Distance from CRE Center (bp)', fontsize=11, fontweight='bold')
    ax.set_ylabel('ATAC Signal', fontsize=11, fontweight='bold')

    direction = "UP" if fc > 1 else "DOWN"
    title = f"{info.get('gene', 'Unknown')} - {info.get('cre_id', f'CRE_{i+1}')}\n"
    title += f"{info.get('chr', '')}:{info.get('start', '')}-{info.get('end', '')}\n"
    title += f"FC={fc:.2f} ({direction}) | {comparison_name}"

    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(loc='best', frameon=True)
    ax.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax)

    plt.tight_layout()

    safe_gene = str(info.get('gene', 'Unknown')).replace('/', '_').replace(' ', '_')
    safe_cre = str(info.get('cre_id', f'CRE_{i+1}')).replace('/', '_').replace(' ', '_')
    output_file = f"{output_dir}/individual_{direction.lower()}_{comparison_name}_{safe_gene}_{safe_cre}.png"
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()

    return output_file


def plot_comparison_metaprofile(data, bin_labels, significant_cres, output_file, title,
                                sample1_name, sample2_name, color1, color2,
                                use_all_cres=False, all_cres_count=None):
    """Create metaprofile comparing two samples for significant CREs."""

    if use_all_cres:
        # Use all CREs
        sample1_data = data[sample1_name]
        sample2_data = data[sample2_name]
        n_cres = all_cres_count
    else:
        # Use only significant CREs
        if len(significant_cres) == 0:
            print(f"  No significant CREs for metaprofile")
            return

        indices = [x[0] for x in significant_cres]
        sample1_data = data[sample1_name][indices, :]
        sample2_data = data[sample2_name][indices, :]
        n_cres = len(significant_cres)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10),
                                     gridspec_kw={'height_ratios': [3, 1]})

    # Calculate mean across CREs
    sample1_mean = np.mean(sample1_data, axis=0)
    sample2_mean = np.mean(sample2_data, axis=0)

    # Calculate SEM
    sample1_sem = stats.sem(sample1_data, axis=0)
    sample2_sem = stats.sem(sample2_data, axis=0)

    # Top panel: Metaprofiles
    ax1.plot(bin_labels, sample1_mean, label=sample1_name, color=color1,
            linewidth=2.5, alpha=0.9)
    ax1.plot(bin_labels, sample2_mean, label=sample2_name, color=color2,
            linewidth=2.5, alpha=0.9)

    ax1.fill_between(bin_labels, sample1_mean - sample1_sem, sample1_mean + sample1_sem,
                     color=color1, alpha=0.2)
    ax1.fill_between(bin_labels, sample2_mean - sample2_sem, sample2_mean + sample2_sem,
                     color=color2, alpha=0.2)

    ax1.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax1.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Mean ATAC Signal', fontsize=12, fontweight='bold')
    ax1.set_title(f'{title}\n(n={n_cres} CREs)', fontsize=14, fontweight='bold', pad=20)
    ax1.legend(loc='best', frameon=True, fontsize=11)
    ax1.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax1)

    # Bottom panel: Difference
    difference = sample2_mean - sample1_mean

    ax2.plot(bin_labels, difference, color='#6C464F', linewidth=2, alpha=0.9)
    ax2.fill_between(bin_labels, 0, difference, where=(difference > 0),
                     color=color2, alpha=0.3, label=f'{sample2_name} > {sample1_name}')
    ax2.fill_between(bin_labels, 0, difference, where=(difference < 0),
                     color=color1, alpha=0.3, label=f'{sample1_name} > {sample2_name}')
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1.5)
    ax2.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

    ax2.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
    ax2.set_ylabel(f'Delta Signal\n({sample2_name} - {sample1_name})',
                   fontsize=11, fontweight='bold')
    ax2.legend(loc='best', frameon=True, fontsize=10)
    ax2.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax2)

    # Add stats
    max_diff = np.max(np.abs(difference))
    mean_diff = np.mean(difference)
    ax2.text(0.02, 0.95, f'Max |Delta|: {max_diff:.4f}\nMean Delta: {mean_diff:.4f}',
            transform=ax2.transAxes, fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Visualize differentially accessible GABA-specific CREs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Filtering Criteria:
  --min-signal: Minimum max signal required (default: 1.0)
  --min-fc: Minimum fold change, applied as FC >= X OR FC <= 1/X (default: 1.5)

Examples:
  # Default thresholds (minSig=2.0, minFC=2.0)
  python 4_visualize_DA_CREs.py

  # Stricter thresholds
  python 4_visualize_DA_CREs.py --min-signal 2.0 --min-fc 2.0

  # Fast mode (metaprofiles only)
  python 4_visualize_DA_CREs.py --skip-individual
        """
    )

    parser.add_argument('--min-signal', type=float, default=1.0,
                       help='Minimum max signal threshold (default: 1.0)')
    parser.add_argument('--min-fc', type=float, default=1.5,
                       help='Minimum fold change threshold (default: 1.5)')
    parser.add_argument('--skip-individual', action='store_true',
                       help='Skip individual CRE plots')
    parser.add_argument('--parallel', type=int, default=8,
                       help='Number of parallel processes (default: 8)')
    parser.add_argument('--individual-dpi', type=int, default=150,
                       help='DPI for individual plots (default: 150)')
    parser.add_argument('--max-individual', type=int, default=100,
                       help='Maximum individual plots per comparison (default: 100)')

    args = parser.parse_args()

    print("="*80)
    print("VISUALIZE DIFFERENTIALLY ACCESSIBLE GABA-SPECIFIC CREs")
    print("="*80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")

    # Change to script directory
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    os.chdir(SCRIPT_DIR)

    # Configuration
    min_signal = args.min_signal
    min_fc = args.min_fc

    print(f"Filtering criteria:")
    print(f"  min_signal: {min_signal}")
    print(f"  min_fc: {min_fc} (up) or {1/min_fc:.2f} (down)")
    print("")

    # Define paths
    matrix_dir = "./output/heatmaps_deeptools"
    output_dir = f"./output/DA_profiles_minSig{min_signal}_minFC{min_fc}"
    bed_file = "./output/GABA_specific_encode_cCREs.bed"
    tsv_file = "./output/GABA_specific_encode_cCREs.tsv"

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/profiles", exist_ok=True)

    # Load CRE info
    cre_info = load_cre_info(tsv_file, bed_file)
    print(f"Loaded info for {len(cre_info)} CREs")
    print("")

    # Define comparisons
    # NOTE: Emx1-Ctrl is excluded as it is a failed sample.
    #       Using Nestin-Ctrl as control for Emx1-Mut comparison.
    comparisons = [
        {
            'name': 'nestin_ctrl_vs_mut',
            'matrix': 'matrix_GABA_specific_nestin.gz',
            'sample1': 'Nestin-Ctrl',
            'sample2': 'Nestin-Mut',
            'color1': '#2E86AB',
            'color2': '#A23B72',
            'title': 'Nestin: Ctrl vs Mut',
            'description': 'Within-genotype mutation effect (Nestin)'
        },
        {
            'name': 'nestin_ctrl_vs_emx1_mut',
            'matrix': 'matrix_GABA_specific_emx1.gz',
            'sample1': 'Nestin-Ctrl',
            'sample2': 'Emx1-Mut',
            'color1': '#2E86AB',
            'color2': '#C73E1D',
            'title': 'Nestin-Ctrl vs Emx1-Mut',
            'description': 'Cross-genotype comparison (Nestin-Ctrl as reference)'
        }
    ]

    # Cross-genotype comparisons from main matrix
    cross_comparisons = [
        {
            'name': 'nestin_mut_vs_emx1_mut',
            'sample1': 'Nestin-Mut',
            'sample2': 'Emx1-Mut',
            'color1': '#A23B72',
            'color2': '#C73E1D',
            'title': 'Nestin-Mut vs Emx1-Mut',
            'description': 'Mutant comparison'
        }
    ]

    all_stats = {}

    # Process within-genotype comparisons
    for comp in comparisons:
        matrix_file = os.path.join(matrix_dir, comp['matrix'])

        if not os.path.exists(matrix_file.replace('.gz', '.tab')):
            print(f"\nSkipping {comp['name']}: matrix not found")
            continue

        print("\n" + "="*80)
        print(f"COMPARISON: {comp['description'].upper()}")
        print(f"({comp['sample1']} vs {comp['sample2']})")
        print("="*80)

        # Read matrix
        data, regions, bin_labels = read_deeptools_matrix(matrix_file)

        if data is None:
            continue

        # Filter significant CREs
        significant_cres = filter_significant_cres(
            data, regions, cre_info,
            comp['sample1'], comp['sample2'],
            min_signal, min_fc
        )

        print(f"\nSignificant CREs: {len(significant_cres)} / {len(regions)}")

        # Split into up/down
        up_cres = [(i, info, fc, ms, m1, m2) for i, info, fc, ms, m1, m2 in significant_cres if fc > 1]
        down_cres = [(i, info, fc, ms, m1, m2) for i, info, fc, ms, m1, m2 in significant_cres if fc < 1]
        print(f"  UP: {len(up_cres)}")
        print(f"  DOWN: {len(down_cres)}")

        # Metaprofile for ALL CREs (overview)
        print(f"\nCreating metaprofile (all CREs)...")
        plot_comparison_metaprofile(
            data, bin_labels, None,
            f"{output_dir}/metaprofile_{comp['name']}_all.png",
            f"{comp['title']} (All CREs)",
            comp['sample1'], comp['sample2'],
            comp['color1'], comp['color2'],
            use_all_cres=True, all_cres_count=len(regions)
        )

        # Metaprofile for significant CREs
        # if len(significant_cres) > 0:
        #     print(f"Creating metaprofile (significant CREs)...")
        #     plot_comparison_metaprofile(
        #         data, bin_labels, significant_cres,
        #         f"{output_dir}/metaprofile_{comp['name']}_significant.png",
        #         f"{comp['title']} (Significant CREs: minSig={min_signal}, minFC={min_fc})",
        #         comp['sample1'], comp['sample2'],
        #         comp['color1'], comp['color2']
        #     )

        # Individual plots (if not skipped)
        if not args.skip_individual and len(significant_cres) > 0:
            n_plots = min(len(significant_cres), args.max_individual)
            print(f"\nCreating {n_plots} individual CRE plots...")

            # Sort by fold change magnitude
            sorted_cres = sorted(significant_cres, key=lambda x: abs(np.log2(x[2])), reverse=True)

            plot_args = []
            for i, info, fc, max_signal, mean1, mean2 in sorted_cres[:n_plots]:
                plot_args.append((
                    i, info, fc, max_signal, mean1, mean2, data, bin_labels,
                    comp['sample1'], comp['sample2'],
                    comp['color1'], comp['color2'],
                    comp['name'], f"{output_dir}/profiles", args.individual_dpi
                ))

            if args.parallel > 1:
                with Pool(args.parallel) as pool:
                    output_files = pool.map(_plot_single_cre, plot_args)
            else:
                output_files = [_plot_single_cre(a) for a in plot_args]

            print(f"  Created {len(output_files)} plots")

        # Save statistics
        all_stats[comp['name']] = {
            'total_cres': len(regions),
            'significant_cres': len(significant_cres),
            'up_cres': len(up_cres),
            'down_cres': len(down_cres)
        }

    # Process cross-genotype comparisons using main matrix
    main_matrix = os.path.join(matrix_dir, "matrix_GABA_specific.gz")
    if os.path.exists(main_matrix.replace('.gz', '.tab')):
        data, regions, bin_labels = read_deeptools_matrix(main_matrix)

        if data is not None:
            for comp in cross_comparisons:
                if comp['sample1'] not in data or comp['sample2'] not in data:
                    continue

                print("\n" + "="*80)
                print(f"COMPARISON: {comp['description'].upper()}")
                print(f"({comp['sample1']} vs {comp['sample2']})")
                print("="*80)

                significant_cres = filter_significant_cres(
                    data, regions, cre_info,
                    comp['sample1'], comp['sample2'],
                    min_signal, min_fc
                )

                print(f"\nSignificant CREs: {len(significant_cres)} / {len(regions)}")

                # Metaprofile
                print(f"\nCreating metaprofile (all CREs)...")
                plot_comparison_metaprofile(
                    data, bin_labels, None,
                    f"{output_dir}/metaprofile_{comp['name']}_all.png",
                    f"{comp['title']} (All CREs)",
                    comp['sample1'], comp['sample2'],
                    comp['color1'], comp['color2'],
                    use_all_cres=True, all_cres_count=len(regions)
                )

                # if len(significant_cres) > 0:
                #     plot_comparison_metaprofile(
                #         data, bin_labels, significant_cres,
                #         f"{output_dir}/metaprofile_{comp['name']}_significant.png",
                #         f"{comp['title']} (Significant CREs)",
                #         comp['sample1'], comp['sample2'],
                #         comp['color1'], comp['color2']
                #     )

                up_cres = [x for x in significant_cres if x[2] > 1]
                down_cres = [x for x in significant_cres if x[2] < 1]

                all_stats[comp['name']] = {
                    'total_cres': len(regions),
                    'significant_cres': len(significant_cres),
                    'up_cres': len(up_cres),
                    'down_cres': len(down_cres)
                }

    # Save statistics
    stats_file = f"{output_dir}/DA_statistics.txt"
    with open(stats_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("DIFFERENTIALLY ACCESSIBLE GABA-SPECIFIC CREs - STATISTICS\n")
        f.write("="*80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Filtering criteria:\n")
        f.write(f"  min_signal: {min_signal}\n")
        f.write(f"  min_fc: {min_fc}\n\n")

        for name, s in all_stats.items():
            f.write(f"{name.upper()}\n")
            f.write("-"*80 + "\n")
            f.write(f"  Total CREs: {s['total_cres']}\n")
            f.write(f"  Significant CREs: {s['significant_cres']} ({100*s['significant_cres']/s['total_cres']:.1f}%)\n")
            f.write(f"    UP: {s['up_cres']}\n")
            f.write(f"    DOWN: {s['down_cres']}\n\n")

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nOutput directory: {output_dir}/")
    print(f"Statistics saved: {stats_file}")
    print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
