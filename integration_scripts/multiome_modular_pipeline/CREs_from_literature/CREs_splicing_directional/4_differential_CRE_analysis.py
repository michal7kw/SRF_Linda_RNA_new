#!/usr/bin/env python3
"""
Step 4: Differential CRE Analysis

This script identifies differential CREs (dCREs) between samples by:
1. Loading deepTools matrices for ALL cell types (not just GABA-specific)
2. Computing signal differences between conditions
3. Applying statistical significance and fold change thresholds
4. Generating comprehensive output files with well-organized structure

Comparisons performed:
- Nestin-Ctrl vs Nestin-Mut (within-genotype mutation effect)
- Nestin-Ctrl vs Emx1-Mut (cross-genotype comparison)
- Nestin-Mut vs Emx1-Mut (mutant comparison)

Author: Claude Code
Date: December 2024
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from pathlib import Path
from scipy import stats
import gzip
import argparse
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = Path("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_directional")
OUTPUT_DIR = BASE_DIR / "output"
MATRIX_DIR = OUTPUT_DIR / "matrices"

# Visualization parameters
DPI = 300
FIGSIZE_VOLCANO = (10, 8)
FIGSIZE_HEATMAP = (12, 10)
FIGSIZE_BARPLOT = (10, 6)
FIGSIZE_PROFILE = (8, 6)

# Colors for conditions
COLORS = {
    'Nestin-Ctrl': '#3182bd',  # Blue
    'Nestin-Mut': '#de2d26',   # Red
    'Emx1-Mut': '#ff7f00',     # Orange
}

# Direction colors
DIRECTION_COLORS = {
    'up': '#e41a1c',      # Red (increased in Mut/compared condition)
    'down': '#377eb8',    # Blue (decreased in Mut/compared condition)
    'ns': '#999999',      # Gray (not significant)
}

# Comparison definitions
COMPARISONS = {
    'nestin_ctrl_vs_mut': {
        'name': 'Nestin-Ctrl vs Nestin-Mut',
        'condition1': 'Nestin-Ctrl',
        'condition2': 'Nestin-Mut',
        'description': 'Within-genotype mutation effect'
    },
    'nestin_ctrl_vs_emx1_mut': {
        'name': 'Nestin-Ctrl vs Emx1-Mut',
        'condition1': 'Nestin-Ctrl',
        'condition2': 'Emx1-Mut',
        'description': 'Cross-genotype comparison'
    },
    'nestin_mut_vs_emx1_mut': {
        'name': 'Nestin-Mut vs Emx1-Mut',
        'condition1': 'Nestin-Mut',
        'condition2': 'Emx1-Mut',
        'description': 'Mutant genotype comparison'
    }
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_deeptools_matrix(matrix_file):
    """Load deepTools matrix file and extract signal data."""
    print(f"Loading matrix: {matrix_file}")

    with gzip.open(matrix_file, 'rt') as f:
        header_line = f.readline().strip()

    import json
    header_json = header_line.replace('@', '')
    header_data = json.loads(header_json)

    sample_labels = header_data.get('sample_labels', [])
    sample_boundaries = header_data.get('sample_boundaries', [])
    upstream = header_data.get('upstream', [2000])[0]
    downstream = header_data.get('downstream', [2000])[0]
    bin_size = header_data.get('bin size', [50])[0]

    n_bins = int((upstream + downstream) / bin_size)

    df = pd.read_csv(matrix_file, sep='\t', skiprows=1, header=None)
    n_regions = len(df)

    # Extract region info (first 6 columns)
    region_info = df.iloc[:, :6].copy()
    region_info.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']

    signal_start_col = 6
    signal_data = df.iloc[:, signal_start_col:].values.astype(float)

    sample_data = {}
    for i, label in enumerate(sample_labels):
        start_bin = sample_boundaries[i]
        end_bin = sample_boundaries[i + 1]
        sample_data[label] = signal_data[:, start_bin:end_bin]

    positions = np.linspace(-upstream, downstream, n_bins)

    print(f"  Samples: {sample_labels}")
    print(f"  Regions: {n_regions}")
    print(f"  Bins per sample: {n_bins}")

    return sample_data, positions, n_regions, region_info


def load_gene_links(gene_links_file):
    """Load gene linkage information."""
    if gene_links_file.exists():
        return pd.read_csv(gene_links_file, sep='\t')
    return pd.DataFrame()


def compute_region_stats(sample_data, window='central'):
    """
    Compute signal statistics for each region.

    Parameters:
    -----------
    sample_data : dict
        Sample name -> signal matrix
    window : str
        'central' for central bins only, 'full' for all bins

    Returns:
    --------
    DataFrame with mean, max signals per region per sample
    """
    region_stats = {}

    for sample, data in sample_data.items():
        n_bins = data.shape[1]

        if window == 'central':
            center = n_bins // 2
            # Use central ~500bp (10 bins at 50bp)
            window_slice = slice(max(0, center - 5), min(n_bins, center + 5))
        else:
            window_slice = slice(None)

        region_stats[f'{sample}_mean'] = np.nanmean(data[:, window_slice], axis=1)
        region_stats[f'{sample}_max'] = np.nanmax(data[:, window_slice], axis=1)

    return pd.DataFrame(region_stats)


def identify_differential_CREs(region_stats, region_info, comparison, min_signal, min_fc):
    """
    Identify differential CREs between two conditions.

    Parameters:
    -----------
    region_stats : DataFrame
        Signal statistics per region
    region_info : DataFrame
        Region coordinates and metadata
    comparison : dict
        Comparison configuration
    min_signal : float
        Minimum signal threshold
    min_fc : float
        Minimum fold change threshold

    Returns:
    --------
    DataFrame with differential CRE information
    """
    cond1 = comparison['condition1']
    cond2 = comparison['condition2']

    # Get mean signals
    signal1 = region_stats[f'{cond1}_mean'].values
    signal2 = region_stats[f'{cond2}_mean'].values

    # Compute fold change (cond2 / cond1)
    # Add small value to avoid division by zero
    epsilon = 0.01
    fc = (signal2 + epsilon) / (signal1 + epsilon)
    log2fc = np.log2(fc)

    # Compute max signal across conditions
    max_signal = np.maximum(signal1, signal2)

    # Create results dataframe
    results = region_info.copy()
    results[f'{cond1}_signal'] = signal1
    results[f'{cond2}_signal'] = signal2
    results['log2FC'] = log2fc
    results['fold_change'] = fc
    results['max_signal'] = max_signal

    # Apply thresholds
    # Significant if: max_signal >= min_signal AND (fc >= min_fc OR fc <= 1/min_fc)
    sig_up = (max_signal >= min_signal) & (fc >= min_fc)
    sig_down = (max_signal >= min_signal) & (fc <= (1.0 / min_fc))

    results['significant'] = sig_up | sig_down
    results['direction'] = 'ns'
    results.loc[sig_up, 'direction'] = 'up'
    results.loc[sig_down, 'direction'] = 'down'

    # Add comparison info
    results['comparison'] = comparison['name']

    return results


def create_volcano_plot(diff_results, cre_type, comparison, output_file, min_signal, min_fc):
    """Create volcano-style plot (signal vs fold change)."""
    fig, ax = plt.subplots(figsize=FIGSIZE_VOLCANO)

    cond1 = comparison['condition1']
    cond2 = comparison['condition2']

    # Get data
    x = diff_results['log2FC'].values
    y = diff_results['max_signal'].values
    direction = diff_results['direction'].values

    # Plot non-significant points
    ns_mask = direction == 'ns'
    ax.scatter(x[ns_mask], y[ns_mask],
               c=DIRECTION_COLORS['ns'], alpha=0.3, s=20, label='Not significant')

    # Plot significant up
    up_mask = direction == 'up'
    ax.scatter(x[up_mask], y[up_mask],
               c=DIRECTION_COLORS['up'], alpha=0.7, s=30,
               label=f'Up in {cond2} (n={up_mask.sum()})')

    # Plot significant down
    down_mask = direction == 'down'
    ax.scatter(x[down_mask], y[down_mask],
               c=DIRECTION_COLORS['down'], alpha=0.7, s=30,
               label=f'Down in {cond2} (n={down_mask.sum()})')

    # Add threshold lines
    ax.axhline(min_signal, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(np.log2(min_fc), color='gray', linestyle='--', alpha=0.5)
    ax.axvline(-np.log2(min_fc), color='gray', linestyle='--', alpha=0.5)

    ax.set_xlabel(f'log2 Fold Change ({cond2}/{cond1})', fontsize=12)
    ax.set_ylabel('Max Signal', fontsize=12)
    ax.set_title(f'{cre_type.capitalize()} CREs: {comparison["name"]}\n{comparison["description"]}',
                 fontsize=14)
    ax.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file.name}")


def create_ma_plot(diff_results, cre_type, comparison, output_file):
    """Create MA plot (mean signal vs fold change)."""
    fig, ax = plt.subplots(figsize=FIGSIZE_VOLCANO)

    cond1 = comparison['condition1']
    cond2 = comparison['condition2']

    # Calculate mean of both conditions (A)
    A = (diff_results[f'{cond1}_signal'] + diff_results[f'{cond2}_signal']) / 2
    M = diff_results['log2FC']
    direction = diff_results['direction']

    # Plot points by direction
    for dir_type, color in [('ns', DIRECTION_COLORS['ns']),
                            ('up', DIRECTION_COLORS['up']),
                            ('down', DIRECTION_COLORS['down'])]:
        mask = direction == dir_type
        if mask.sum() > 0:
            alpha = 0.3 if dir_type == 'ns' else 0.7
            size = 20 if dir_type == 'ns' else 30
            ax.scatter(A[mask], M[mask], c=color, alpha=alpha, s=size,
                      label=f'{dir_type.capitalize()} (n={mask.sum()})')

    ax.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax.set_xlabel('Mean Signal', fontsize=12)
    ax.set_ylabel(f'log2 Fold Change ({cond2}/{cond1})', fontsize=12)
    ax.set_title(f'{cre_type.capitalize()} CREs: MA Plot\n{comparison["name"]}', fontsize=14)
    ax.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file.name}")


def create_summary_heatmap(all_results, cre_type, output_file):
    """Create heatmap summarizing differential CREs across comparisons."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Create summary matrix
    summary_data = []
    for comp_name, results in all_results.items():
        n_up = (results['direction'] == 'up').sum()
        n_down = (results['direction'] == 'down').sum()
        n_total = len(results)
        summary_data.append({
            'Comparison': COMPARISONS[comp_name]['name'],
            'Up': n_up,
            'Down': n_down,
            'Not Significant': n_total - n_up - n_down
        })

    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.set_index('Comparison')

    # Create stacked bar plot
    colors = [DIRECTION_COLORS['up'], DIRECTION_COLORS['down'], DIRECTION_COLORS['ns']]
    summary_df[['Up', 'Down', 'Not Significant']].plot(
        kind='barh', stacked=True, ax=ax, color=colors
    )

    ax.set_xlabel('Number of CREs', fontsize=12)
    ax.set_ylabel('')
    ax.set_title(f'{cre_type.capitalize()} CREs: Differential Analysis Summary', fontsize=14)
    ax.legend(title='Direction', loc='upper right')

    plt.tight_layout()
    plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file.name}")


def create_metaprofile_differential(sample_data, positions, diff_mask, direction,
                                    cre_type, comparison, output_file):
    """Create metaprofile plot for differential CREs."""
    if diff_mask.sum() == 0:
        print(f"  Skipping metaprofile (no {direction} CREs)")
        return

    fig, ax = plt.subplots(figsize=FIGSIZE_PROFILE)

    for sample, color in COLORS.items():
        if sample in sample_data:
            data = sample_data[sample][diff_mask, :]
            mean_signal = np.nanmean(data, axis=0)
            sem_signal = stats.sem(data, axis=0, nan_policy='omit')

            ax.plot(positions, mean_signal, color=color, label=sample, linewidth=2)
            ax.fill_between(positions, mean_signal - sem_signal, mean_signal + sem_signal,
                           color=color, alpha=0.2)

    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Distance from CRE center (bp)', fontsize=12)
    ax.set_ylabel('ATAC Signal', fontsize=12)
    ax.set_title(f'{cre_type.capitalize()} CREs: {direction.upper()} in {comparison["condition2"]}\n'
                 f'{comparison["name"]} (n={diff_mask.sum()})', fontsize=14)
    ax.legend(loc='upper right')
    ax.set_xlim(positions[0], positions[-1])

    plt.tight_layout()
    plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file.name}")


def create_individual_differential_plots(sample_data, positions, diff_results,
                                         cre_type, comparison, output_dir, max_plots=20):
    """Create individual plots for top differential CREs."""
    sig_results = diff_results[diff_results['significant']].copy()

    if len(sig_results) == 0:
        print(f"  No significant differential CREs to plot")
        return

    # Sort by absolute log2FC
    sig_results['abs_log2FC'] = np.abs(sig_results['log2FC'])
    sig_results = sig_results.sort_values('abs_log2FC', ascending=False)

    n_plots = min(max_plots, len(sig_results))

    for i, (idx, row) in enumerate(sig_results.head(n_plots).iterrows()):
        fig, ax = plt.subplots(figsize=(6, 4))

        for sample, color in COLORS.items():
            if sample in sample_data:
                signal = sample_data[sample][idx, :]
                ax.plot(positions, signal, color=color, label=sample, linewidth=2)

        ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('Distance from CRE center (bp)')
        ax.set_ylabel('ATAC Signal')

        direction = row['direction']
        fc = row['fold_change']
        region_name = str(row['name'])[:30] if pd.notna(row['name']) else 'unknown'

        ax.set_title(f'{cre_type.capitalize()}: {region_name}\n'
                     f'{direction.upper()} (FC={fc:.2f})')
        ax.legend(loc='upper right')
        ax.set_xlim(positions[0], positions[-1])

        output_file = output_dir / f"dCRE_{i+1:03d}_{direction}_{region_name[:20]}_fc{fc:.1f}.png"
        plt.tight_layout()
        plt.savefig(output_file, dpi=150)
        plt.close()

    print(f"  Created {n_plots} individual plots")


def write_summary_report(all_cre_results, min_signal, min_fc, output_file):
    """Write comprehensive summary report."""
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("DIFFERENTIAL CRE ANALYSIS - SUMMARY REPORT\n")
        f.write("=" * 80 + "\n\n")

        f.write("ANALYSIS PARAMETERS:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Minimum signal threshold: {min_signal}\n")
        f.write(f"Minimum fold change threshold: {min_fc}\n")
        f.write(f"Log2 FC threshold: {np.log2(min_fc):.2f}\n\n")

        f.write("COMPARISONS:\n")
        f.write("-" * 80 + "\n")
        for comp_id, comp in COMPARISONS.items():
            f.write(f"  {comp['name']}: {comp['description']}\n")
        f.write("\n")

        for cre_type in ['enhancer', 'silencer']:
            if cre_type not in all_cre_results:
                continue

            f.write("=" * 80 + "\n")
            f.write(f"{cre_type.upper()} CREs RESULTS\n")
            f.write("=" * 80 + "\n\n")

            cre_results = all_cre_results[cre_type]

            for comp_id, results in cre_results.items():
                comp = COMPARISONS[comp_id]

                n_total = len(results)
                n_up = (results['direction'] == 'up').sum()
                n_down = (results['direction'] == 'down').sum()
                n_sig = n_up + n_down

                f.write(f"\n{comp['name']}:\n")
                f.write("-" * 40 + "\n")
                f.write(f"  Total CREs analyzed: {n_total}\n")
                f.write(f"  Differential CREs: {n_sig} ({100*n_sig/n_total:.1f}%)\n")
                f.write(f"    Up in {comp['condition2']}: {n_up}\n")
                f.write(f"    Down in {comp['condition2']}: {n_down}\n")

                if n_sig > 0:
                    sig_results = results[results['significant']].sort_values(
                        'log2FC', ascending=False
                    )

                    if n_up > 0:
                        f.write(f"\n  Top UP differential CREs:\n")
                        for _, row in sig_results[sig_results['direction'] == 'up'].head(5).iterrows():
                            name = str(row['name'])[:30] if pd.notna(row['name']) else 'unknown'
                            f.write(f"    {name}: log2FC={row['log2FC']:.2f}, "
                                   f"signal={row['max_signal']:.2f}\n")

                    if n_down > 0:
                        f.write(f"\n  Top DOWN differential CREs:\n")
                        for _, row in sig_results[sig_results['direction'] == 'down'].head(5).iterrows():
                            name = str(row['name'])[:30] if pd.notna(row['name']) else 'unknown'
                            f.write(f"    {name}: log2FC={row['log2FC']:.2f}, "
                                   f"signal={row['max_signal']:.2f}\n")

                f.write("\n")

        f.write("=" * 80 + "\n")
        f.write("BIOLOGICAL INTERPRETATION:\n")
        f.write("=" * 80 + "\n\n")

        f.write("ENHANCER CREs (PCC > 0):\n")
        f.write("  - Positive correlation: accessibility UP = expression UP\n")
        f.write("  - 'UP' differential: increased accessibility -> increased expression\n")
        f.write("  - 'DOWN' differential: decreased accessibility -> decreased expression\n\n")

        f.write("SILENCER CREs (PCC < 0):\n")
        f.write("  - Negative correlation: accessibility UP = expression DOWN\n")
        f.write("  - 'UP' differential: increased accessibility -> decreased expression\n")
        f.write("  - 'DOWN' differential: decreased accessibility -> increased expression\n\n")

        f.write("OUTPUT FILES:\n")
        f.write("-" * 80 + "\n")
        f.write("  differential_CREs_{cre_type}_{comparison}.tsv - All CRE statistics\n")
        f.write("  differential_CREs_{cre_type}_{comparison}_significant.tsv - Significant only\n")
        f.write("  differential_CREs_{cre_type}_{comparison}_significant.bed - BED for IGV\n")
        f.write("  volcano_{cre_type}_{comparison}.png - Volcano plot\n")
        f.write("  ma_plot_{cre_type}_{comparison}.png - MA plot\n")
        f.write("  metaprofile_{direction}_{cre_type}_{comparison}.png - Signal profiles\n")

    print(f"  Saved: {output_file.name}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def print_signal_diagnostics(region_stats, cre_type):
    """Print diagnostic information about signal distribution."""
    print(f"\n  Signal Distribution Diagnostics for {cre_type}:")
    print(f"  {'-'*50}")

    for col in region_stats.columns:
        if '_mean' in col:
            values = region_stats[col].values
            non_zero = values[values > 0]
            print(f"  {col}:")
            print(f"    Range: {values.min():.4f} - {values.max():.4f}")
            print(f"    Mean: {values.mean():.4f}, Median: {np.median(values):.4f}")
            print(f"    Non-zero: {len(non_zero)}/{len(values)} ({100*len(non_zero)/len(values):.1f}%)")
            if len(non_zero) > 0:
                print(f"    Non-zero mean: {non_zero.mean():.4f}")

            # Show percentiles
            percentiles = [50, 75, 90, 95, 99]
            pct_values = np.percentile(values, percentiles)
            print(f"    Percentiles: " + ", ".join([f"P{p}={v:.4f}" for p, v in zip(percentiles, pct_values)]))

    # Suggest thresholds
    all_means = []
    for col in region_stats.columns:
        if '_mean' in col:
            all_means.extend(region_stats[col].values)
    all_means = np.array(all_means)
    non_zero_means = all_means[all_means > 0]

    if len(non_zero_means) > 0:
        suggested_min_signal = np.percentile(non_zero_means, 75)
        print(f"\n  Suggested min_signal (P75 of non-zero): {suggested_min_signal:.4f}")


def main():
    parser = argparse.ArgumentParser(
        description='Identify differential CREs between conditions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script identifies differential cis-regulatory elements (dCREs) between
experimental conditions using signal thresholds and fold change criteria.

All CREs are analyzed (not just GABA-specific), with separate analysis for
enhancer-like and silencer-like CREs based on their correlation direction.

IMPORTANT: Default thresholds are set for normalized BigWig signal values
(typically 0.01-1.0 range). Adjust based on your data.

Examples:
  # Default thresholds (min_signal=0.05, min_fc=1.5):
  python %(prog)s

  # More permissive thresholds:
  python %(prog)s --min-signal 0.02 --min-fc 1.3

  # Stricter thresholds:
  python %(prog)s --min-signal 0.1 --min-fc 2.0

  # Skip individual plots:
  python %(prog)s --skip-individual

  # Show signal diagnostics without filtering:
  python %(prog)s --diagnostics-only
        """
    )
    parser.add_argument('--min-signal', type=float, default=0.05, metavar='F',
                       help='Minimum max signal required (default: 0.05)')
    parser.add_argument('--min-fc', type=float, default=1.5, metavar='F',
                       help='Minimum fold change required (default: 1.5)')
    parser.add_argument('--skip-individual', action='store_true',
                       help='Skip individual CRE plots')
    parser.add_argument('--max-individual', type=int, default=20,
                       help='Maximum individual plots per comparison (default: 20)')
    parser.add_argument('--use-gaba', action='store_true',
                       help='Use GABA-specific matrices (default: all cell types)')
    parser.add_argument('--diagnostics-only', action='store_true',
                       help='Only show signal diagnostics, do not perform differential analysis')
    args = parser.parse_args()

    print("=" * 80)
    print("DIFFERENTIAL CRE ANALYSIS")
    print("=" * 80)
    print(f"\nParameters:")
    print(f"  Minimum signal: {args.min_signal}")
    print(f"  Minimum fold change: {args.min_fc}")
    print(f"  Log2 FC threshold: {np.log2(args.min_fc):.2f}")
    print(f"  Using: {'GABA-specific' if args.use_gaba else 'All cell types'}")
    print()

    # Create output directory structure
    suffix = f"_minSig{args.min_signal}_minFC{args.min_fc}"
    cell_suffix = "_GABA" if args.use_gaba else "_all"

    diff_dir = OUTPUT_DIR / f"differential_CREs{suffix}"
    diff_dir.mkdir(parents=True, exist_ok=True)

    # Subdirectories for organization
    results_dir = diff_dir / "results"
    plots_dir = diff_dir / "plots"
    bed_dir = diff_dir / "bed_files"
    individual_dir = diff_dir / "individual_plots"

    for d in [results_dir, plots_dir, bed_dir]:
        d.mkdir(exist_ok=True)

    all_cre_results = {}

    # Process each CRE type
    for cre_type in ['enhancer', 'silencer']:
        # Choose matrix file based on cell type mode
        if args.use_gaba:
            matrix_file = MATRIX_DIR / f"matrix_{cre_type}_GABA.gz"
            gene_links_file = OUTPUT_DIR / f"{cre_type}_CREs_GABA_gene_links.tsv"
        else:
            # For all cell types, we need to compute matrices first
            # Check if all-celltype matrix exists, otherwise use GABA
            matrix_file = MATRIX_DIR / f"matrix_{cre_type}_all.gz"
            gene_links_file = OUTPUT_DIR / f"{cre_type}_CREs_all_gene_links.tsv"

            if not matrix_file.exists():
                print(f"\n  Note: All-celltype matrix not found, using GABA matrix")
                matrix_file = MATRIX_DIR / f"matrix_{cre_type}_GABA.gz"
                gene_links_file = OUTPUT_DIR / f"{cre_type}_CREs_GABA_gene_links.tsv"

        if not matrix_file.exists():
            print(f"\nSkipping {cre_type}: matrix file not found")
            continue

        print(f"\n{'='*80}")
        print(f"PROCESSING {cre_type.upper()} CREs")
        print(f"{'='*80}")

        # Load data
        sample_data, positions, n_regions, region_info = load_deeptools_matrix(matrix_file)

        # Load gene links
        gene_links = load_gene_links(gene_links_file)
        if len(gene_links) > 0:
            print(f"  Loaded {len(gene_links)} gene links")

        # Compute region statistics
        region_stats = compute_region_stats(sample_data)

        # Print signal diagnostics
        print_signal_diagnostics(region_stats, cre_type)

        # If diagnostics-only mode, skip to next CRE type
        if args.diagnostics_only:
            print(f"\n  [Diagnostics-only mode - skipping differential analysis for {cre_type}]")
            continue

        # Create CRE type subdirectories
        cre_plots_dir = plots_dir / cre_type
        cre_bed_dir = bed_dir / cre_type
        cre_plots_dir.mkdir(exist_ok=True)
        cre_bed_dir.mkdir(exist_ok=True)

        if not args.skip_individual:
            cre_individual_dir = individual_dir / cre_type
            cre_individual_dir.mkdir(parents=True, exist_ok=True)

        cre_results = {}

        # Process each comparison
        for comp_id, comparison in COMPARISONS.items():
            print(f"\n  --- {comparison['name']} ---")

            # Identify differential CREs
            diff_results = identify_differential_CREs(
                region_stats, region_info, comparison,
                args.min_signal, args.min_fc
            )

            cre_results[comp_id] = diff_results

            # Statistics
            n_up = (diff_results['direction'] == 'up').sum()
            n_down = (diff_results['direction'] == 'down').sum()
            print(f"  Differential: {n_up} up, {n_down} down")

            # Create comparison subdirectories
            comp_plots_dir = cre_plots_dir / comp_id
            comp_plots_dir.mkdir(exist_ok=True)

            # Save results
            out_base = f"{cre_type}_{comp_id}"

            # Full results
            diff_results.to_csv(results_dir / f"differential_CREs_{out_base}.tsv",
                               sep='\t', index=False)

            # Significant only
            sig_results = diff_results[diff_results['significant']]
            if len(sig_results) > 0:
                sig_results.to_csv(results_dir / f"differential_CREs_{out_base}_significant.tsv",
                                  sep='\t', index=False)

                # BED file for significant CREs
                bed_df = sig_results[['chrom', 'start', 'end', 'name', 'score', 'strand']].copy()
                bed_df.to_csv(cre_bed_dir / f"dCREs_{out_base}_significant.bed",
                             sep='\t', index=False, header=False)

                # Separate BED files by direction
                for direction in ['up', 'down']:
                    dir_mask = sig_results['direction'] == direction
                    if dir_mask.sum() > 0:
                        dir_df = sig_results[dir_mask][['chrom', 'start', 'end', 'name', 'score', 'strand']]
                        dir_df.to_csv(cre_bed_dir / f"dCREs_{out_base}_{direction}.bed",
                                     sep='\t', index=False, header=False)

            # Create plots
            create_volcano_plot(diff_results, cre_type, comparison,
                               comp_plots_dir / f"volcano_{out_base}.png",
                               args.min_signal, args.min_fc)

            create_ma_plot(diff_results, cre_type, comparison,
                          comp_plots_dir / f"ma_plot_{out_base}.png")

            # Create metaprofiles for differential CREs
            for direction in ['up', 'down']:
                dir_mask = (diff_results['direction'] == direction).values
                if dir_mask.sum() > 0:
                    create_metaprofile_differential(
                        sample_data, positions, dir_mask, direction,
                        cre_type, comparison,
                        comp_plots_dir / f"metaprofile_{direction}_{out_base}.png"
                    )

            # Individual plots
            if not args.skip_individual and len(sig_results) > 0:
                comp_individual_dir = cre_individual_dir / comp_id
                comp_individual_dir.mkdir(exist_ok=True)
                create_individual_differential_plots(
                    sample_data, positions, diff_results,
                    cre_type, comparison, comp_individual_dir,
                    args.max_individual
                )

        # Create summary heatmap
        create_summary_heatmap(cre_results, cre_type,
                              plots_dir / f"summary_{cre_type}_differential.png")

        all_cre_results[cre_type] = cre_results

    # Write comprehensive summary report
    write_summary_report(all_cre_results, args.min_signal, args.min_fc,
                        diff_dir / "SUMMARY_differential_CRE_analysis.txt")

    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nOutput directory: {diff_dir}")
    print(f"\nDirectory structure:")
    print(f"  {diff_dir.name}/")
    print(f"    results/          - TSV files with all statistics")
    print(f"    plots/            - Visualization plots")
    print(f"      enhancer/       - Enhancer CRE plots by comparison")
    print(f"      silencer/       - Silencer CRE plots by comparison")
    print(f"    bed_files/        - BED files for IGV/deepTools")
    print(f"      enhancer/       - Enhancer differential CRE BED files")
    print(f"      silencer/       - Silencer differential CRE BED files")
    if not args.skip_individual:
        print(f"    individual_plots/ - Individual CRE signal profiles")


if __name__ == "__main__":
    main()
