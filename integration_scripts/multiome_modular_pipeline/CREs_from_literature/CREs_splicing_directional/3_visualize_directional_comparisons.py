#!/usr/bin/env python3
"""
Step 3: Visualize Directional Comparisons and Identify Nestin-Specific Loss

This script:
1. Loads deepTools matrices from Step 2
2. Creates publication-quality comparison plots
3. Identifies CREs with NESTIN-SPECIFIC LOSS of accessibility
4. Generates individual CRE plots for significant hits

Nestin-Specific Loss Pattern:
- Nestin-Ctrl: HIGH signal (baseline)
- Nestin-Mut: LOW signal (Nestin-specific loss)
- Emx1-Mut: HIGH/NORMAL signal (preserved in Emx1 background)

This pattern indicates regulatory elements specifically affected by
SRF mutation in the Nestin-Cre lineage but not in Emx1-Cre lineage.

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
VIZ_DIR = OUTPUT_DIR / "visualizations"
NESTIN_LOSS_DIR = OUTPUT_DIR / "nestin_specific_loss"

# Parameters for identifying Nestin-specific loss
MIN_CTRL_SIGNAL = 0.5       # Minimum signal in Nestin-Ctrl (baseline must be detectable)
MIN_FOLD_CHANGE = 1.5       # Minimum fold change (Ctrl/Mut)
MAX_EMX1_LOSS_RATIO = 0.7   # Emx1-Mut must retain at least 70% of Ctrl signal

# Visualization parameters
DPI = 300
FIGSIZE_PROFILE = (8, 6)
FIGSIZE_HEATMAP = (10, 8)
FIGSIZE_INDIVIDUAL = (6, 4)

# Colors
COLORS = {
    'Nestin-Ctrl': '#3182bd',  # Blue
    'Nestin-Mut': '#de2d26',   # Red
    'Emx1-Mut': '#ff7f00',     # Orange
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_deeptools_matrix(matrix_file):
    """Load deepTools matrix file and extract signal data."""
    print(f"Loading matrix: {matrix_file}")

    with gzip.open(matrix_file, 'rt') as f:
        # Read header line starting with @
        header_line = f.readline().strip()

    # Parse header for sample names and parameters
    import json
    header_json = header_line.replace('@', '')
    header_data = json.loads(header_json)

    sample_labels = header_data.get('sample_labels', [])
    sample_boundaries = header_data.get('sample_boundaries', [])
    upstream = header_data.get('upstream', [2000])[0]
    downstream = header_data.get('downstream', [2000])[0]
    bin_size = header_data.get('bin size', [50])[0]

    # Number of bins per sample
    n_bins = int((upstream + downstream) / bin_size)
    n_samples = len(sample_labels)

    # Load data using pandas to handle mixed types
    # deepTools matrix format: chr, start, end, name, score, strand, then signal values
    df = pd.read_csv(matrix_file, sep='\t', skiprows=1, header=None)

    n_regions = len(df)

    # First 6 columns are region info, rest are signal values
    # Column indices: 0=chr, 1=start, 2=end, 3=name, 4=score, 5=strand, 6+=signal
    signal_start_col = 6
    signal_data = df.iloc[:, signal_start_col:].values.astype(float)

    # Split signal data by sample using sample_boundaries
    sample_data = {}
    for i, label in enumerate(sample_labels):
        start_bin = sample_boundaries[i]
        end_bin = sample_boundaries[i + 1]
        sample_data[label] = signal_data[:, start_bin:end_bin]

    # Create position array
    positions = np.linspace(-upstream, downstream, n_bins)

    print(f"  Samples: {sample_labels}")
    print(f"  Regions: {n_regions}")
    print(f"  Bins per sample: {n_bins}")
    print(f"  Signal columns: {signal_data.shape[1]}")

    return sample_data, positions, n_regions


def load_bed_with_genes(bed_file, gene_links_file):
    """Load BED file and associate with gene information."""
    # Load BED
    bed_df = pd.read_csv(bed_file, sep='\t', header=None,
                         names=['chrom', 'start', 'end', 'cre_id', 'score', 'strand'])

    # Load gene links if available
    if gene_links_file.exists():
        gene_df = pd.read_csv(gene_links_file, sep='\t')
        # Create CRE -> gene mapping (one CRE may link to multiple genes)
        cre_genes = gene_df.groupby(['chrom', 'start', 'end'])['gene'].apply(
            lambda x: ', '.join(sorted(set(x)))
        ).reset_index()
        cre_genes.columns = ['chrom', 'start', 'end', 'genes']

        bed_df = bed_df.merge(cre_genes, on=['chrom', 'start', 'end'], how='left')
        bed_df['genes'] = bed_df['genes'].fillna('unknown')
    else:
        bed_df['genes'] = 'unknown'

    return bed_df


def compute_region_stats(sample_data):
    """Compute mean signal for each region across each sample."""
    region_stats = {}
    for sample, data in sample_data.items():
        # Mean signal across central bins (within 500bp of center)
        n_bins = data.shape[1]
        center = n_bins // 2
        central_bins = slice(center - 5, center + 5)  # ~500bp window

        region_stats[sample] = np.nanmean(data[:, central_bins], axis=1)

    return pd.DataFrame(region_stats)


def identify_nestin_specific_loss(region_stats, bed_df, cre_type):
    """Identify CREs with Nestin-specific loss of accessibility."""
    print(f"\nIdentifying Nestin-specific loss in {cre_type} CREs...")

    # Calculate metrics
    df = region_stats.copy()
    df['cre_idx'] = range(len(df))

    # Fold changes
    df['fc_nestin'] = df['Nestin-Ctrl'] / (df['Nestin-Mut'] + 0.01)  # Add small value to avoid division by zero
    df['fc_emx1'] = df['Nestin-Ctrl'] / (df['Emx1-Mut'] + 0.01)
    df['emx1_retention'] = df['Emx1-Mut'] / (df['Nestin-Ctrl'] + 0.01)

    # Criteria for Nestin-specific loss:
    # 1. Nestin-Ctrl has detectable signal
    # 2. Nestin-Mut shows significant loss (high fold change)
    # 3. Emx1-Mut retains signal (low fold change / high retention)

    nestin_loss_mask = (
        (df['Nestin-Ctrl'] >= MIN_CTRL_SIGNAL) &           # Baseline detectable
        (df['fc_nestin'] >= MIN_FOLD_CHANGE) &              # Loss in Nestin-Mut
        (df['emx1_retention'] >= MAX_EMX1_LOSS_RATIO)       # Preserved in Emx1-Mut
    )

    nestin_loss_df = df[nestin_loss_mask].copy()

    print(f"  Total CREs: {len(df)}")
    print(f"  With detectable signal: {(df['Nestin-Ctrl'] >= MIN_CTRL_SIGNAL).sum()}")
    print(f"  With Nestin-specific loss: {len(nestin_loss_df)}")

    if len(nestin_loss_df) > 0:
        # Add BED information
        for col in bed_df.columns:
            nestin_loss_df[col] = bed_df.iloc[nestin_loss_df['cre_idx']][col].values

        # Sort by fold change
        nestin_loss_df = nestin_loss_df.sort_values('fc_nestin', ascending=False)

        print(f"\n  Top 10 Nestin-specific loss CREs:")
        print(f"  {'Gene':<15} {'FC(Ctrl/Mut)':<15} {'Emx1 Retention':<15}")
        print(f"  {'-'*45}")
        for _, row in nestin_loss_df.head(10).iterrows():
            genes = row.get('genes', 'unknown')[:15]
            print(f"  {genes:<15} {row['fc_nestin']:<15.2f} {row['emx1_retention']:<15.2f}")

    return nestin_loss_df, df


def create_metaprofile_with_difference(sample_data, positions, title, output_file):
    """Create metaprofile plot with difference panel."""
    fig, axes = plt.subplots(2, 1, figsize=FIGSIZE_PROFILE, height_ratios=[3, 1],
                             gridspec_kw={'hspace': 0.1})

    # Top panel: Metaprofiles
    ax1 = axes[0]

    for sample, color in COLORS.items():
        if sample in sample_data:
            data = sample_data[sample]
            mean_signal = np.nanmean(data, axis=0)
            sem_signal = stats.sem(data, axis=0, nan_policy='omit')

            ax1.plot(positions, mean_signal, color=color, label=sample, linewidth=2)
            ax1.fill_between(positions, mean_signal - sem_signal, mean_signal + sem_signal,
                           color=color, alpha=0.2)

    ax1.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax1.set_ylabel('ATAC Signal', fontsize=12)
    ax1.set_title(title, fontsize=14)
    ax1.legend(loc='upper right')
    ax1.set_xlim(positions[0], positions[-1])
    ax1.tick_params(labelbottom=False)

    # Bottom panel: Difference (Nestin-Mut - Nestin-Ctrl)
    ax2 = axes[1]

    if 'Nestin-Ctrl' in sample_data and 'Nestin-Mut' in sample_data:
        ctrl_mean = np.nanmean(sample_data['Nestin-Ctrl'], axis=0)
        mut_mean = np.nanmean(sample_data['Nestin-Mut'], axis=0)
        diff = mut_mean - ctrl_mean

        ax2.fill_between(positions, 0, diff, where=(diff < 0),
                        color=COLORS['Nestin-Mut'], alpha=0.5, label='Loss in Mut')
        ax2.fill_between(positions, 0, diff, where=(diff >= 0),
                        color='gray', alpha=0.3, label='Gain in Mut')
        ax2.axhline(0, color='black', linestyle='-', linewidth=0.5)
        ax2.axvline(0, color='gray', linestyle='--', alpha=0.5)

    ax2.set_xlabel('Distance from CRE center (bp)', fontsize=12)
    ax2.set_ylabel('Δ Signal', fontsize=10)
    ax2.set_xlim(positions[0], positions[-1])

    plt.tight_layout()
    plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file.name}")


def create_three_way_comparison_plot(region_stats, title, output_file, highlight_nestin_loss=None):
    """Create scatter plot comparing all three conditions."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Plot 1: Nestin-Ctrl vs Nestin-Mut
    ax1 = axes[0]
    ax1.scatter(region_stats['Nestin-Ctrl'], region_stats['Nestin-Mut'],
               alpha=0.3, s=10, c='gray')
    if highlight_nestin_loss is not None and len(highlight_nestin_loss) > 0:
        ax1.scatter(highlight_nestin_loss['Nestin-Ctrl'], highlight_nestin_loss['Nestin-Mut'],
                   alpha=0.7, s=20, c='red', label='Nestin-specific loss')
    max_val = max(region_stats['Nestin-Ctrl'].max(), region_stats['Nestin-Mut'].max())
    ax1.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    ax1.set_xlabel('Nestin-Ctrl')
    ax1.set_ylabel('Nestin-Mut')
    ax1.set_title('Within-genotype\n(Mutation Effect)')
    ax1.legend()

    # Plot 2: Nestin-Ctrl vs Emx1-Mut
    ax2 = axes[1]
    ax2.scatter(region_stats['Nestin-Ctrl'], region_stats['Emx1-Mut'],
               alpha=0.3, s=10, c='gray')
    if highlight_nestin_loss is not None and len(highlight_nestin_loss) > 0:
        ax2.scatter(highlight_nestin_loss['Nestin-Ctrl'], highlight_nestin_loss['Emx1-Mut'],
                   alpha=0.7, s=20, c='red', label='Nestin-specific loss')
    max_val = max(region_stats['Nestin-Ctrl'].max(), region_stats['Emx1-Mut'].max())
    ax2.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    ax2.set_xlabel('Nestin-Ctrl')
    ax2.set_ylabel('Emx1-Mut')
    ax2.set_title('Cross-genotype\n(Ctrl vs Emx1-Mut)')
    ax2.legend()

    # Plot 3: Nestin-Mut vs Emx1-Mut
    ax3 = axes[2]
    ax3.scatter(region_stats['Nestin-Mut'], region_stats['Emx1-Mut'],
               alpha=0.3, s=10, c='gray')
    if highlight_nestin_loss is not None and len(highlight_nestin_loss) > 0:
        ax3.scatter(highlight_nestin_loss['Nestin-Mut'], highlight_nestin_loss['Emx1-Mut'],
                   alpha=0.7, s=20, c='red', label='Nestin-specific loss')
    max_val = max(region_stats['Nestin-Mut'].max(), region_stats['Emx1-Mut'].max())
    ax3.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    ax3.set_xlabel('Nestin-Mut')
    ax3.set_ylabel('Emx1-Mut')
    ax3.set_title('Mutant Comparison\n(Genotype Effect)')
    ax3.legend()

    fig.suptitle(title, fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file.name}")


def check_individual_plot_thresholds(sample_data, cre_idx, min_signal, min_fc):
    """
    Check if a CRE passes thresholds for individual plotting.

    Parameters:
    -----------
    sample_data : dict
        Dictionary with sample names as keys and signal arrays as values
    cre_idx : int
        Index of the CRE to check
    min_signal : float
        Minimum max signal required (at least one condition must exceed this)
    min_fc : float
        Minimum fold change required (Ctrl/Mut >= min_fc OR Mut/Ctrl >= min_fc)

    Returns:
    --------
    passes : bool
        True if CRE passes thresholds
    info : dict
        Dictionary with signal and fold change information
    """
    # Get signal for each condition
    ctrl_signal = sample_data['Nestin-Ctrl'][cre_idx, :]
    mut_signal = sample_data['Nestin-Mut'][cre_idx, :]
    emx1_signal = sample_data['Emx1-Mut'][cre_idx, :]

    # Calculate mean and max signals
    ctrl_mean = np.nanmean(ctrl_signal)
    mut_mean = np.nanmean(mut_signal)
    emx1_mean = np.nanmean(emx1_signal)

    max_signal = max(np.nanmax(ctrl_signal), np.nanmax(mut_signal), np.nanmax(emx1_signal))

    # Calculate fold change (handle zeros)
    if ctrl_mean > 0.01:
        fc = ctrl_mean / mut_mean if mut_mean > 0.01 else 100.0
    elif mut_mean > 0.01:
        fc = 0.01  # Very low FC if ctrl is 0 and mut is positive
    else:
        fc = 1.0  # Both near zero

    info = {
        'ctrl_mean': ctrl_mean,
        'mut_mean': mut_mean,
        'emx1_mean': emx1_mean,
        'max_signal': max_signal,
        'fc': fc
    }

    # Check thresholds
    # 1. Minimum signal check (at least one condition must have signal)
    if max_signal < min_signal:
        return False, info

    # 2. Fold change check (must be significant change up OR down)
    # FC >= min_fc OR FC <= 1/min_fc
    if not (fc >= min_fc or fc <= (1.0 / min_fc)):
        return False, info

    return True, info


def create_individual_cre_plot(sample_data, positions, cre_idx, cre_info, output_file, dpi=150):
    """Create individual CRE plot showing all three conditions."""
    fig, ax = plt.subplots(figsize=FIGSIZE_INDIVIDUAL)

    for sample, color in COLORS.items():
        if sample in sample_data:
            signal = sample_data[sample][cre_idx, :]
            ax.plot(positions, signal, color=color, label=sample, linewidth=2)

    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Distance from CRE center (bp)')
    ax.set_ylabel('ATAC Signal')

    genes = cre_info.get('genes', 'unknown')
    fc = cre_info.get('fc_nestin', 0)
    ax.set_title(f"Gene: {genes}\nFold Change (Ctrl/Nestin-Mut): {fc:.2f}")

    ax.legend(loc='upper right')
    ax.set_xlim(positions[0], positions[-1])

    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi)
    plt.close()


def create_summary_barplot(all_stats, output_file):
    """Create summary barplot of mean signals."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for idx, (cre_type, stats_df) in enumerate(all_stats.items()):
        ax = axes[idx]

        means = [stats_df[col].mean() for col in ['Nestin-Ctrl', 'Nestin-Mut', 'Emx1-Mut']]
        sems = [stats.sem(stats_df[col]) for col in ['Nestin-Ctrl', 'Nestin-Mut', 'Emx1-Mut']]

        x = np.arange(3)
        bars = ax.bar(x, means, yerr=sems, capsize=5,
                     color=[COLORS['Nestin-Ctrl'], COLORS['Nestin-Mut'], COLORS['Emx1-Mut']])

        ax.set_xticks(x)
        ax.set_xticklabels(['Nestin-Ctrl', 'Nestin-Mut', 'Emx1-Mut'])
        ax.set_ylabel('Mean ATAC Signal')
        ax.set_title(f'{cre_type.capitalize()} CREs')

        # Add significance annotation
        # T-test between Nestin-Ctrl and Nestin-Mut
        t_stat, p_val = stats.ttest_rel(stats_df['Nestin-Ctrl'], stats_df['Nestin-Mut'])
        if p_val < 0.001:
            sig = '***'
        elif p_val < 0.01:
            sig = '**'
        elif p_val < 0.05:
            sig = '*'
        else:
            sig = 'ns'

        y_max = max(means) + max(sems) * 1.5
        ax.plot([0, 1], [y_max, y_max], 'k-', linewidth=1)
        ax.text(0.5, y_max * 1.02, sig, ha='center', fontsize=12)

    plt.tight_layout()
    plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file.name}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Visualize directional CRE comparisons',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Filtering options for individual plots:
  --min-signal: Minimum max signal required (at least one condition must exceed)
  --min-fc: Minimum fold change required (Ctrl/Mut >= min_fc OR Mut/Ctrl >= min_fc)

Examples:
  # Quick overview (metaprofiles only):
  python %(prog)s --skip-individual

  # Full analysis with default thresholds (min_signal=2.0, min_fc=2.0):
  python %(prog)s

  # Even stricter thresholds (fewer, most significant plots):
  python %(prog)s --min-signal 3.0 --min-fc 3.0

  # More permissive thresholds (more plots):
  python %(prog)s --min-signal 1.0 --min-fc 1.5
        """
    )
    parser.add_argument('--skip-individual', action='store_true',
                       help='Skip individual CRE plots (faster)')
    parser.add_argument('--max-individual', type=int, default=50,
                       help='Maximum number of individual plots to create')
    parser.add_argument('--min-signal', type=float, default=2.0, metavar='F',
                       help='Minimum max signal required to plot (default: 2.0)')
    parser.add_argument('--min-fc', type=float, default=2.0, metavar='F',
                       help='Minimum fold change required to plot (default: 2.0)')
    parser.add_argument('--individual-dpi', type=int, default=150, metavar='N',
                       help='DPI for individual plots (default: 150)')
    args = parser.parse_args()

    print("=" * 80)
    print("VISUALIZING DIRECTIONAL COMPARISONS")
    print("=" * 80)

    # Display settings
    if args.skip_individual:
        print("\nFAST MODE: Skipping individual CRE plots")
    else:
        print(f"\nFULL MODE: Creating individual plots")
        print(f"   -> Filtering: Min Signal >= {args.min_signal}, Min FC >= {args.min_fc}")
        print(f"   -> Max individual plots per CRE type: {args.max_individual}")
        print(f"   -> Individual plot DPI: {args.individual_dpi}")
    print()

    # Create output directories with threshold suffix for individual plots
    suffix = f"_minSig{args.min_signal}_minFC{args.min_fc}"

    VIZ_DIR.mkdir(parents=True, exist_ok=True)
    NESTIN_LOSS_DIR.mkdir(parents=True, exist_ok=True)

    all_stats = {}
    all_nestin_loss = {}

    # Process each CRE type
    for cre_type in ['enhancer', 'silencer']:
        matrix_file = MATRIX_DIR / f"matrix_{cre_type}_GABA.gz"
        bed_file = OUTPUT_DIR / f"{cre_type}_CREs_GABA.bed"
        gene_links_file = OUTPUT_DIR / f"{cre_type}_CREs_GABA_gene_links.tsv"

        if not matrix_file.exists():
            print(f"\nSkipping {cre_type}: matrix file not found")
            continue

        print(f"\n{'='*80}")
        print(f"PROCESSING {cre_type.upper()} CREs")
        print(f"{'='*80}")

        # Load data
        sample_data, positions, n_regions = load_deeptools_matrix(matrix_file)
        bed_df = load_bed_with_genes(bed_file, gene_links_file)

        # Compute region statistics
        region_stats = compute_region_stats(sample_data)
        all_stats[cre_type] = region_stats

        # Create metaprofile with difference plot
        create_metaprofile_with_difference(
            sample_data, positions,
            f"{cre_type.capitalize()} CREs at Splicing Genes (GABA)",
            VIZ_DIR / f"metaprofile_{cre_type}_with_difference.png"
        )

        # Identify Nestin-specific loss
        nestin_loss_df, full_df = identify_nestin_specific_loss(region_stats, bed_df, cre_type)
        all_nestin_loss[cre_type] = nestin_loss_df

        # Create three-way comparison plot
        create_three_way_comparison_plot(
            region_stats,
            f"{cre_type.capitalize()} CREs: Three-way Comparison",
            VIZ_DIR / f"scatter_{cre_type}_three_way_comparison.png",
            highlight_nestin_loss=nestin_loss_df
        )

        # Save Nestin-specific loss results
        if len(nestin_loss_df) > 0:
            nestin_loss_df.to_csv(
                NESTIN_LOSS_DIR / f"nestin_specific_loss_{cre_type}.tsv",
                sep='\t', index=False
            )
            print(f"\n  Saved: nestin_specific_loss_{cre_type}.tsv")

            # Create individual CRE plots for top hits (with threshold filtering)
            if not args.skip_individual:
                individual_dir = NESTIN_LOSS_DIR / f"individual_{cre_type}{suffix}"
                individual_dir.mkdir(exist_ok=True)

                print(f"\n  Creating individual plots (minSig>={args.min_signal}, minFC>={args.min_fc})...")
                print(f"  Checking {len(nestin_loss_df)} Nestin-specific loss CREs...")

                plots_created = 0
                plots_skipped = 0

                for i, (_, row) in enumerate(nestin_loss_df.iterrows()):
                    if plots_created >= args.max_individual:
                        break

                    cre_idx = row['cre_idx']

                    # Check thresholds
                    passes, info = check_individual_plot_thresholds(
                        sample_data, cre_idx, args.min_signal, args.min_fc
                    )

                    if not passes:
                        plots_skipped += 1
                        continue

                    genes = str(row.get('genes', 'unknown')).replace(', ', '_')[:30]
                    fc = row['fc_nestin']

                    output_file = individual_dir / f"cre_{plots_created+1:03d}_{genes}_fc{fc:.1f}.png"
                    create_individual_cre_plot(sample_data, positions, cre_idx, row, output_file,
                                              dpi=args.individual_dpi)
                    plots_created += 1

                print(f"  Created {plots_created} individual plots (skipped {plots_skipped} below threshold)")

    # Create summary comparison
    if len(all_stats) > 0:
        print(f"\n{'='*80}")
        print("CREATING SUMMARY PLOTS")
        print(f"{'='*80}")

        create_summary_barplot(all_stats, VIZ_DIR / "summary_barplot_all_conditions.png")

    # Write summary report
    write_summary_report(all_stats, all_nestin_loss)

    print(f"\n{'='*80}")
    print("VISUALIZATION COMPLETE")
    print(f"{'='*80}")
    print(f"\nOutput directories:")
    print(f"  Visualizations: {VIZ_DIR}")
    print(f"  Nestin-specific loss: {NESTIN_LOSS_DIR}")


def write_summary_report(all_stats, all_nestin_loss):
    """Write summary report of the analysis."""
    report_path = NESTIN_LOSS_DIR / "SUMMARY_nestin_specific_loss.txt"

    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("NESTIN-SPECIFIC LOSS OF ACCESSIBILITY ANALYSIS\n")
        f.write("=" * 80 + "\n\n")

        f.write("ANALYSIS PARAMETERS:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Minimum Ctrl signal: {MIN_CTRL_SIGNAL}\n")
        f.write(f"Minimum fold change (Ctrl/Mut): {MIN_FOLD_CHANGE}\n")
        f.write(f"Minimum Emx1 retention ratio: {MAX_EMX1_LOSS_RATIO}\n\n")

        f.write("NESTIN-SPECIFIC LOSS PATTERN:\n")
        f.write("-" * 80 + "\n")
        f.write("CREs where:\n")
        f.write("  - Nestin-Ctrl has HIGH signal (baseline)\n")
        f.write("  - Nestin-Mut has LOW signal (accessibility loss)\n")
        f.write("  - Emx1-Mut has HIGH/NORMAL signal (preserved)\n\n")

        f.write("RESULTS SUMMARY:\n")
        f.write("-" * 80 + "\n\n")

        for cre_type in ['enhancer', 'silencer']:
            if cre_type in all_stats:
                stats_df = all_stats[cre_type]
                loss_df = all_nestin_loss.get(cre_type, pd.DataFrame())

                f.write(f"{cre_type.upper()} CREs:\n")
                f.write(f"  Total CREs: {len(stats_df)}\n")
                f.write(f"  With Nestin-specific loss: {len(loss_df)}\n")
                f.write(f"  Percentage: {100*len(loss_df)/len(stats_df):.1f}%\n\n")

                if len(loss_df) > 0:
                    f.write(f"  Mean signals:\n")
                    f.write(f"    Nestin-Ctrl: {stats_df['Nestin-Ctrl'].mean():.3f}\n")
                    f.write(f"    Nestin-Mut:  {stats_df['Nestin-Mut'].mean():.3f}\n")
                    f.write(f"    Emx1-Mut:    {stats_df['Emx1-Mut'].mean():.3f}\n\n")

                    f.write(f"  Top genes with Nestin-specific loss:\n")
                    for _, row in loss_df.head(20).iterrows():
                        genes = row.get('genes', 'unknown')
                        fc = row['fc_nestin']
                        retention = row['emx1_retention']
                        f.write(f"    {genes}: FC={fc:.2f}, Emx1_retention={retention:.2f}\n")
                    f.write("\n")

        f.write("=" * 80 + "\n")
        f.write("BIOLOGICAL INTERPRETATION:\n")
        f.write("=" * 80 + "\n\n")

        f.write("ENHANCER CREs with Nestin-specific loss:\n")
        f.write("  -> These regulatory elements normally activate splicing genes\n")
        f.write("  -> Loss of accessibility in Nestin-Mut -> reduced expression\n")
        f.write("  -> Preserved in Emx1-Mut -> genotype-specific effect\n")
        f.write("  -> May explain splicing defects specifically in Nestin lineage\n\n")

        f.write("SILENCER CREs with Nestin-specific loss:\n")
        f.write("  -> These regulatory elements normally repress splicing genes\n")
        f.write("  -> Loss of accessibility in Nestin-Mut -> increased expression\n")
        f.write("  -> May indicate compensatory or aberrant gene activation\n\n")

    print(f"  Saved: {report_path.name}")


if __name__ == "__main__":
    main()
