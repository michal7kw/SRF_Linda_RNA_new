#!/usr/bin/env python3
"""
Visualize BigWig Signal at Splicing Gene CREs

This script directly reads BigWig files and creates publication-quality
visualizations of ATAC-seq signal at CREs, with focus on Ctrl vs Mut comparisons.

INPUT:
- BigWig files (GABA_Nestin-Ctrl.bw, GABA_Nestin-Mut.bw, etc.)
- BED file with CRE coordinates

OUTPUT:
- Publication-quality signal profiles
- Ctrl vs Mut comparison plots
- Difference plots (Mut - Ctrl)
- Individual CRE profiles (optional, can be skipped for faster execution)
- Statistical summary

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
import pyBigWig
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

def extract_signal_from_bigwig(bw_file, chrom, start, end, window=2000, bins=80):
    """
    Extract signal from BigWig file around a genomic region.

    Parameters:
    -----------
    bw_file : str
        Path to BigWig file
    chrom : str
        Chromosome name
    start : int
        Region start coordinate
    end : int
        Region end coordinate
    window : int
        Window size around region center (±window bp)
    bins : int
        Number of bins for signal

    Returns:
    --------
    signal : array
        Signal values across the region
    positions : array
        Genomic positions for each bin
    """
    # Calculate region center
    center = (start + end) // 2

    # Define window around center
    window_start = max(0, center - window)
    window_end = center + window

    # Open BigWig file
    bw = pyBigWig.open(bw_file)

    # Extract signal
    try:
        signal = bw.stats(chrom, window_start, window_end, type="mean", nBins=bins)
        # Replace None values with 0
        signal = np.array([x if x is not None else 0.0 for x in signal])
    except:
        # If chromosome not found or other error, return zeros
        signal = np.zeros(bins)

    bw.close()

    # Create position array (relative to center)
    positions = np.linspace(-window, window, bins)

    return signal, positions

def plot_metaprofile_comparison(signals_ctrl, signals_mut, positions, output_file,
                                 title, cre_info):
    """
    Create metaprofile comparing Ctrl vs Mut with confidence intervals.

    Parameters:
    -----------
    signals_ctrl : list of arrays
        Signal arrays for control samples
    signals_mut : list of arrays
        Signal arrays for mutant samples
    positions : array
        Position array (bp from center)
    output_file : str
        Output file path
    title : str
        Plot title
    cre_info : list of dicts
        Information about each CRE (gene, cre_id, coordinates)
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10),
                                     gridspec_kw={'height_ratios': [3, 1]})

    colors = {'Ctrl': '#2E86AB', 'Mut': '#C73E1D'}

    # Top panel: Metaprofiles with SEM
    for label, signals in [('Ctrl', signals_ctrl), ('Mut', signals_mut)]:
        if len(signals) == 0:
            continue

        # Calculate mean and SEM across all CREs
        mean_signal = np.mean(signals, axis=0)
        sem_signal = stats.sem(signals, axis=0)

        # Plot mean line
        ax1.plot(positions, mean_signal, label=f'{label} (n={len(signals)} CREs)',
                color=colors[label], linewidth=2.5, alpha=0.9)

        # Add confidence interval (mean ± SEM)
        ax1.fill_between(positions,
                         mean_signal - sem_signal,
                         mean_signal + sem_signal,
                         color=colors[label], alpha=0.2)

    ax1.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5,
                label='CRE Center')
    ax1.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('ATAC Signal', fontsize=12, fontweight='bold')
    ax1.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax1.legend(loc='best', frameon=True, fontsize=11)
    ax1.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax1)

    # Bottom panel: Difference (Mut - Ctrl)
    if len(signals_ctrl) > 0 and len(signals_mut) > 0:
        mean_ctrl = np.mean(signals_ctrl, axis=0)
        mean_mut = np.mean(signals_mut, axis=0)

        difference = mean_mut - mean_ctrl

        # Plot difference
        ax2.plot(positions, difference, color='#6C464F', linewidth=2, alpha=0.9)
        ax2.fill_between(positions, 0, difference, where=(difference > 0),
                         color='#C73E1D', alpha=0.3, label='Mut > Ctrl')
        ax2.fill_between(positions, 0, difference, where=(difference < 0),
                         color='#2E86AB', alpha=0.3, label='Ctrl > Mut')
        ax2.axhline(y=0, color='black', linestyle='-', linewidth=1.5)
        ax2.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

        ax2.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Δ ATAC Signal\n(Mut - Ctrl)', fontsize=11, fontweight='bold')
        ax2.legend(loc='best', frameon=True, fontsize=10)
        ax2.grid(True, alpha=0.3, linestyle=':')
        sns.despine(ax=ax2)

        # Add statistical annotation
        max_diff = np.max(np.abs(difference))
        mean_diff = np.mean(difference)
        fold_change = np.mean(mean_mut) / np.mean(mean_ctrl) if np.mean(mean_ctrl) > 0 else np.nan

        stats_text = f'Max |Δ|: {max_diff:.4f}\nMean Δ: {mean_diff:.4f}\nFold Change: {fold_change:.2f}x'
        ax2.text(0.02, 0.95, stats_text,
                transform=ax2.transAxes, fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()

def _plot_single_cre_bigwig(args):
    """
    Helper function for parallel processing of individual CRE plots from BigWig.
    This must be a top-level function for multiprocessing to work.
    """
    (bw_ctrl, bw_mut, chrom, start, end, window, bins, gene, cre_id,
     genotype, output_file, dpi) = args

    # Extract signals
    signal_ctrl, positions = extract_signal_from_bigwig(bw_ctrl, chrom, start, end, window, bins)
    signal_mut, _ = extract_signal_from_bigwig(bw_mut, chrom, start, end, window, bins)

    # Create plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8),
                                     gridspec_kw={'height_ratios': [2, 1]})

    colors = {'Ctrl': '#2E86AB', 'Mut': '#C73E1D'}

    # Top panel: Signal profiles
    ax1.plot(positions, signal_ctrl, label='Ctrl',
            color=colors['Ctrl'], linewidth=2.5, alpha=0.9)
    ax1.plot(positions, signal_mut, label='Mut',
            color=colors['Mut'], linewidth=2.5, alpha=0.9)

    ax1.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax1.set_ylabel('ATAC Signal', fontsize=11, fontweight='bold')
    ax1.set_title(f'{gene} - {cre_id}\n{chrom}:{start:,}-{end:,} ({genotype})',
                 fontsize=12, fontweight='bold')
    ax1.legend(loc='best', frameon=True)
    ax1.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax1)

    # Bottom panel: Difference
    difference = signal_mut - signal_ctrl

    ax2.plot(positions, difference, color='#6C464F', linewidth=2, alpha=0.9)
    ax2.fill_between(positions, 0, difference, where=(difference > 0),
                     color='#C73E1D', alpha=0.3, label='Mut > Ctrl')
    ax2.fill_between(positions, 0, difference, where=(difference < 0),
                     color='#2E86AB', alpha=0.3, label='Ctrl > Mut')
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1.5)
    ax2.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

    ax2.set_xlabel('Distance from CRE Center (bp)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Δ Signal\n(Mut - Ctrl)', fontsize=10, fontweight='bold')
    ax2.legend(loc='best', frameon=True, fontsize=9)
    ax2.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax2)

    # Statistics
    max_diff = np.max(np.abs(difference))
    mean_diff = np.mean(difference)
    mean_ctrl = np.mean(signal_ctrl)
    mean_mut = np.mean(signal_mut)
    fold_change = mean_mut / mean_ctrl if mean_ctrl > 0 else np.nan

    stats_text = f'Mean Ctrl: {mean_ctrl:.4f}\nMean Mut: {mean_mut:.4f}\nFold Change: {fold_change:.2f}x'
    ax2.text(0.98, 0.95, stats_text,
            transform=ax2.transAxes, fontsize=8, verticalalignment='top',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()

    # Return statistics
    return {
        'mean_ctrl': mean_ctrl,
        'mean_mut': mean_mut,
        'mean_diff': mean_diff,
        'max_diff': max_diff,
        'fold_change': fold_change
    }

def plot_individual_cre(signal_ctrl, signal_mut, positions, gene, cre_id,
                        chrom, start, end, genotype, output_file, dpi=300):
    """
    Create individual plot for a single CRE showing Ctrl vs Mut.

    Parameters:
    -----------
    dpi : int
        DPI for the plot (default: 300)
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8),
                                     gridspec_kw={'height_ratios': [2, 1]})

    colors = {'Ctrl': '#2E86AB', 'Mut': '#C73E1D'}

    # Top panel: Signal profiles
    ax1.plot(positions, signal_ctrl, label='Ctrl',
            color=colors['Ctrl'], linewidth=2.5, alpha=0.9)
    ax1.plot(positions, signal_mut, label='Mut',
            color=colors['Mut'], linewidth=2.5, alpha=0.9)

    ax1.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax1.set_ylabel('ATAC Signal', fontsize=11, fontweight='bold')
    ax1.set_title(f'{gene} - {cre_id}\n{chrom}:{start:,}-{end:,} ({genotype})',
                 fontsize=12, fontweight='bold')
    ax1.legend(loc='best', frameon=True)
    ax1.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax1)

    # Bottom panel: Difference
    difference = signal_mut - signal_ctrl

    ax2.plot(positions, difference, color='#6C464F', linewidth=2, alpha=0.9)
    ax2.fill_between(positions, 0, difference, where=(difference > 0),
                     color='#C73E1D', alpha=0.3, label='Mut > Ctrl')
    ax2.fill_between(positions, 0, difference, where=(difference < 0),
                     color='#2E86AB', alpha=0.3, label='Ctrl > Mut')
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1.5)
    ax2.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

    ax2.set_xlabel('Distance from CRE Center (bp)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Δ Signal\n(Mut - Ctrl)', fontsize=10, fontweight='bold')
    ax2.legend(loc='best', frameon=True, fontsize=9)
    ax2.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax2)

    # Statistics
    max_diff = np.max(np.abs(difference))
    mean_diff = np.mean(difference)
    mean_ctrl = np.mean(signal_ctrl)
    mean_mut = np.mean(signal_mut)
    fold_change = mean_mut / mean_ctrl if mean_ctrl > 0 else np.nan

    stats_text = f'Mean Ctrl: {mean_ctrl:.4f}\nMean Mut: {mean_mut:.4f}\nFold Change: {fold_change:.2f}x'
    ax2.text(0.98, 0.95, stats_text,
            transform=ax2.transAxes, fontsize=8, verticalalignment='top',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()

def create_summary_table(cre_info_list, stats_list, output_file):
    """
    Create summary table with statistics for all CREs.
    """
    df = pd.DataFrame({
        'Gene': [info['gene'] for info in cre_info_list],
        'CRE_ID': [info['cre_id'] for info in cre_info_list],
        'Chromosome': [info['chrom'] for info in cre_info_list],
        'Start': [info['start'] for info in cre_info_list],
        'End': [info['end'] for info in cre_info_list],
        'Mean_Ctrl': [s['mean_ctrl'] for s in stats_list],
        'Mean_Mut': [s['mean_mut'] for s in stats_list],
        'Fold_Change': [s['fold_change'] for s in stats_list],
        'Mean_Difference': [s['mean_diff'] for s in stats_list],
        'Max_Abs_Difference': [s['max_diff'] for s in stats_list]
    })

    df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
    print(f"  Saved summary table: {output_file}")

    return df

def main():
    """Main execution function"""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Visualize BigWig signal at splicing gene CREs',
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
    print("VISUALIZE BIGWIG SIGNAL AT SPLICING GENE CREs")
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
        print(f"   → Metaprofile DPI: 300 (always high quality)")
    print()

    os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_genes_paper")

    # Define paths
    bigwig_base = "../../signac_results_L1/bigwig_tracks_L1/by_celltype"
    bed_file = "./output/splicing_genes_CREs_all.bed"
    # NOTE: Using ALL CREs (not just GABA-specific) but with GABA BigWig files ONLY for signal analysis
    tsv_file = "./output/splicing_genes_CREs_all_celltypes.tsv"
    output_dir = f"./output/bigwig_profiles_minSig{args.min_signal}_minFC{args.min_fc}"

    os.makedirs(output_dir, exist_ok=True)

    # Parameters
    window = 2000  # ±2kb around CRE center
    bins = 80      # Number of bins

    # Read BED file
    print("Reading CRE coordinates...")
    bed_df = pd.read_csv(bed_file, sep='\t', header=None,
                        names=['chr', 'start', 'end', 'cre_id', 'score', 'strand'])
    print(f"  Found {len(bed_df)} CREs")

    # Read gene mapping
    print("Reading gene mappings...")
    tsv_df = pd.read_csv(tsv_file, sep='\t')
    gene_map = tsv_df.groupby('cCRE1')['Gene'].first().to_dict()

    # ========================================================================
    # NESTIN Analysis
    # ========================================================================
    print("\n" + "="*80)
    print("NESTIN: Extracting BigWig signals")
    print("="*80)

    bw_nestin_ctrl = f"{bigwig_base}/GABA_Nestin-Ctrl.bw"
    bw_nestin_mut = f"{bigwig_base}/GABA_Nestin-Mut.bw"

    if not os.path.exists(bw_nestin_ctrl):
        print(f"ERROR: BigWig file not found: {bw_nestin_ctrl}")
        return

    signals_nestin_ctrl = []
    signals_nestin_mut = []
    cre_info_nestin = []
    stats_nestin = []
    positions = None

    for i, row in bed_df.iterrows():
        cre_id = row['cre_id']
        gene = gene_map.get(cre_id, 'Unknown')
        chrom = row['chr']
        start = row['start']
        end = row['end']

        print(f"\n  Processing: {gene} - {cre_id} ({chrom}:{start}-{end})")

        # Extract signals
        signal_ctrl, pos = extract_signal_from_bigwig(bw_nestin_ctrl, chrom, start, end, window, bins)
        signal_mut, _ = extract_signal_from_bigwig(bw_nestin_mut, chrom, start, end, window, bins)

        if positions is None:
            positions = pos

        signals_nestin_ctrl.append(signal_ctrl)
        signals_nestin_mut.append(signal_mut)

        # Store CRE info
        info = {
            'gene': gene,
            'cre_id': cre_id,
            'chrom': chrom,
            'start': start,
            'end': end
        }
        cre_info_nestin.append(info)

        # Calculate statistics
        mean_ctrl = np.mean(signal_ctrl)
        mean_mut = np.mean(signal_mut)
        mean_diff = mean_mut - mean_ctrl
        max_diff = np.max(np.abs(signal_mut - signal_ctrl))
        fold_change = mean_mut / mean_ctrl if mean_ctrl > 0 else np.nan

        stats = {
            'mean_ctrl': mean_ctrl,
            'mean_mut': mean_mut,
            'mean_diff': mean_diff,
            'max_diff': max_diff,
            'fold_change': fold_change
        }
        stats_nestin.append(stats)

        print(f"    Mean Ctrl: {mean_ctrl:.6f}, Mean Mut: {mean_mut:.6f}, FC: {fold_change:.2f}x")

        # Check significance for plotting
        is_significant = True
        
        # 1. Min signal check
        max_val = max(np.max(signal_ctrl), np.max(signal_mut))
        if max_val < args.min_signal:
            is_significant = False
            
        # 2. Fold change check
        # Handle zero/low values
        if mean_ctrl < 0.01 and mean_mut < 0.01:
            fc_check = 1.0
        elif mean_ctrl < 0.01:
            fc_check = 100.0
        else:
            fc_check = mean_mut / mean_ctrl
            
        if not (fc_check >= args.min_fc or fc_check <= (1.0/args.min_fc)):
            is_significant = False
            
        # Store significance
        info['is_significant'] = is_significant

    # Create individual plots (optional, parallel or sequential)
    if not args.skip_individual:
        print(f"\nCreating Nestin individual CRE plots...")

        if args.parallel > 1:
            # Parallel processing
            print(f"  Using {args.parallel} parallel processes...")
            plot_args = []
            for i, row in bed_df.iterrows():
                # Check if significant (using index to match cre_info_nestin)
                if not cre_info_nestin[i]['is_significant']:
                    continue
                    
                cre_id = row['cre_id']
                gene = gene_map.get(cre_id, 'Unknown')
                chrom = row['chr']
                start = row['start']
                end = row['end']
                output_file = f"{output_dir}/individual_Nestin_{gene}_{cre_id}.png"

                args_tuple = (bw_nestin_ctrl, bw_nestin_mut, chrom, start, end,
                             window, bins, gene, cre_id, 'Nestin', output_file,
                             args.individual_dpi)
                plot_args.append(args_tuple)

            with Pool(args.parallel) as pool:
                stats_parallel = pool.map(_plot_single_cre_bigwig, plot_args)

            print(f"  Created {len(stats_parallel)} Nestin plots in parallel")
        else:
            # Sequential processing
            for i, row in bed_df.iterrows():
                # Check if significant
                if not cre_info_nestin[i]['is_significant']:
                    continue

                cre_id = row['cre_id']
                gene = gene_map.get(cre_id, 'Unknown')
                chrom = row['chr']
                start = row['start']
                end = row['end']

                signal_ctrl = signals_nestin_ctrl[i]
                signal_mut = signals_nestin_mut[i]

                output_file = f"{output_dir}/individual_Nestin_{gene}_{cre_id}.png"
                plot_individual_cre(signal_ctrl, signal_mut, positions, gene, cre_id,
                                   chrom, start, end, 'Nestin', output_file,
                                   dpi=args.individual_dpi)
    else:
        print(f"\n⚡ Skipping Nestin individual CRE plots (fast mode)")

    # Create metaprofile
    print("\nCreating Nestin metaprofile...")
    plot_metaprofile_comparison(
        signals_nestin_ctrl, signals_nestin_mut, positions,
        f"{output_dir}/metaprofile_nestin_ctrl_vs_mut.png",
        "Nestin: ATAC Signal at Splicing Gene CREs\nCtrl vs Mut Comparison",
        cre_info_nestin
    )

    # Create summary table
    df_nestin = create_summary_table(cre_info_nestin, stats_nestin,
                                     f"{output_dir}/summary_nestin.tsv")

    # ========================================================================
    # EMX1 Analysis
    # ========================================================================
    print("\n" + "="*80)
    print("EMX1: Extracting BigWig signals")
    print("="*80)

    bw_emx1_ctrl = f"{bigwig_base}/GABA_Nestin-Ctrl.bw"
    bw_emx1_mut = f"{bigwig_base}/GABA_Emx1-Mut.bw"

    if not os.path.exists(bw_emx1_ctrl):
        print(f"ERROR: BigWig file not found: {bw_emx1_ctrl}")
        return

    signals_emx1_ctrl = []
    signals_emx1_mut = []
    cre_info_emx1 = []
    stats_emx1 = []

    for i, row in bed_df.iterrows():
        cre_id = row['cre_id']
        gene = gene_map.get(cre_id, 'Unknown')
        chrom = row['chr']
        start = row['start']
        end = row['end']

        print(f"\n  Processing: {gene} - {cre_id} ({chrom}:{start}-{end})")

        # Extract signals
        signal_ctrl, _ = extract_signal_from_bigwig(bw_emx1_ctrl, chrom, start, end, window, bins)
        signal_mut, _ = extract_signal_from_bigwig(bw_emx1_mut, chrom, start, end, window, bins)

        signals_emx1_ctrl.append(signal_ctrl)
        signals_emx1_mut.append(signal_mut)

        # Store CRE info
        info = {
            'gene': gene,
            'cre_id': cre_id,
            'chrom': chrom,
            'start': start,
            'end': end
        }
        cre_info_emx1.append(info)

        # Calculate statistics
        mean_ctrl = np.mean(signal_ctrl)
        mean_mut = np.mean(signal_mut)
        mean_diff = mean_mut - mean_ctrl
        max_diff = np.max(np.abs(signal_mut - signal_ctrl))
        fold_change = mean_mut / mean_ctrl if mean_ctrl > 0 else np.nan

        stats = {
            'mean_ctrl': mean_ctrl,
            'mean_mut': mean_mut,
            'mean_diff': mean_diff,
            'max_diff': max_diff,
            'fold_change': fold_change
        }
        stats_emx1.append(stats)

        print(f"    Mean Ctrl: {mean_ctrl:.6f}, Mean Mut: {mean_mut:.6f}, FC: {fold_change:.2f}x")

        # Check significance for plotting
        is_significant = True
        
        # 1. Min signal check
        max_val = max(np.max(signal_ctrl), np.max(signal_mut))
        if max_val < args.min_signal:
            is_significant = False
            
        # 2. Fold change check
        if mean_ctrl < 0.01 and mean_mut < 0.01:
            fc_check = 1.0
        elif mean_ctrl < 0.01:
            fc_check = 100.0
        else:
            fc_check = mean_mut / mean_ctrl
            
        if not (fc_check >= args.min_fc or fc_check <= (1.0/args.min_fc)):
            is_significant = False
            
        # Store significance
        info['is_significant'] = is_significant

    # Create individual plots (optional, parallel or sequential)
    if not args.skip_individual:
        print(f"\nCreating Emx1 individual CRE plots...")

        if args.parallel > 1:
            # Parallel processing
            print(f"  Using {args.parallel} parallel processes...")
            plot_args = []
            for i, row in bed_df.iterrows():
                # Check if significant
                if not cre_info_emx1[i]['is_significant']:
                    continue

                cre_id = row['cre_id']
                gene = gene_map.get(cre_id, 'Unknown')
                chrom = row['chr']
                start = row['start']
                end = row['end']
                output_file = f"{output_dir}/individual_Emx1_{gene}_{cre_id}.png"

                args_tuple = (bw_emx1_ctrl, bw_emx1_mut, chrom, start, end,
                             window, bins, gene, cre_id, 'Emx1', output_file,
                             args.individual_dpi)
                plot_args.append(args_tuple)

            with Pool(args.parallel) as pool:
                stats_parallel = pool.map(_plot_single_cre_bigwig, plot_args)

            print(f"  Created {len(stats_parallel)} Emx1 plots in parallel")
        else:
            # Sequential processing
            for i, row in bed_df.iterrows():
                # Check if significant
                if not cre_info_emx1[i]['is_significant']:
                    continue

                cre_id = row['cre_id']
                gene = gene_map.get(cre_id, 'Unknown')
                chrom = row['chr']
                start = row['start']
                end = row['end']

                signal_ctrl = signals_emx1_ctrl[i]
                signal_mut = signals_emx1_mut[i]

                output_file = f"{output_dir}/individual_Emx1_{gene}_{cre_id}.png"
                plot_individual_cre(signal_ctrl, signal_mut, positions, gene, cre_id,
                                   chrom, start, end, 'Emx1', output_file,
                                   dpi=args.individual_dpi)
    else:
        print(f"\n⚡ Skipping Emx1 individual CRE plots (fast mode)")

    # Create metaprofile
    print("\nCreating Emx1 metaprofile...")
    plot_metaprofile_comparison(
        signals_emx1_ctrl, signals_emx1_mut, positions,
        f"{output_dir}/metaprofile_emx1_ctrl_vs_mut.png",
        "Emx1: ATAC Signal at Splicing Gene CREs\nNestin-Ctrl vs Emx1-Mut Comparison",
        cre_info_emx1
    )

    # Create summary table
    df_emx1 = create_summary_table(cre_info_emx1, stats_emx1,
                                   f"{output_dir}/summary_emx1.tsv")

    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nOutput directory: {output_dir}/")
    print("\nGenerated files:")
    print("  Metaprofiles (always created, 300 DPI):")
    print("    - metaprofile_nestin_ctrl_vs_mut.png")
    print("    - metaprofile_emx1_ctrl_vs_mut.png")

    if not args.skip_individual:
        print(f"\n  Individual CRE plots (DPI: {args.individual_dpi}):")
        print(f"    - individual_Nestin_*.png ({len(bed_df)} plots)")
        print(f"    - individual_Emx1_*.png ({len(bed_df)} plots)")
        print(f"    Total: {len(bed_df) * 2} individual plots")
        if args.parallel > 1:
            print(f"    (Created using {args.parallel} parallel processes)")
    else:
        print("\n  Individual CRE plots: SKIPPED (fast mode)")
        print("    To create individual plots, run without --skip-individual")

    print("\n  Summary tables:")
    print("    - summary_nestin.tsv")
    print("    - summary_emx1.tsv")
    print()

    # Print summary statistics
    print("="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print("\nNestin:")
    print(df_nestin.to_string(index=False))
    print("\nEmx1:")
    print(df_emx1.to_string(index=False))
    print()

if __name__ == "__main__":
    main()
