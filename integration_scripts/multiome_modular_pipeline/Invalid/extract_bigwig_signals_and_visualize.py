#!/usr/bin/env python3
"""
Extract BigWig Signals at Literature CREs and Create Visualizations

Purpose:
- For genes with literature CREs, extract ATAC signal from BigWig files
- Compare Ctrl vs Mut signal at each CRE
- Create visual comparisons for hypothesis generation

This script should be run AFTER analyze_GABA_DEGs_with_literature_CREs.py

INPUT FILES:
- output/GABA_DEG_analysis/DEG_CRE_links_*.tsv: DEG-CRE linkage files for each genotype
- ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Nestin-Ctrl.bw: ATAC signal for Nestin control
- ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Nestin-Mut.bw: ATAC signal for Nestin mutant
- ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Emx1-Ctrl.bw: ATAC signal for Emx1 control
- ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Emx1-Mut.bw: ATAC signal for Emx1 mutant

OUTPUT FILES:
- output/GABA_DEG_analysis/bigwig_signals/CRE_signals_*.tsv: Signal data for each genotype
- output/GABA_DEG_analysis/plots/scatter_ctrl_vs_mut_*.png: Signal comparison scatter plots
- output/GABA_DEG_analysis/plots/signal_log2fc_dist_*.png: Log2FC distribution plots
- output/GABA_DEG_analysis/plots/gene_vs_cre_*.png: Gene expression vs CRE accessibility plots
- output/GABA_DEG_analysis/plots/top_genes_heatmap_*.png: Heatmaps of top DEGs
- output/GABA_DEG_analysis/bigwig_signals/bigwig_analysis_report.txt: Complete analysis report
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns

# Try to import pyBigWig
try:
    import pyBigWig
    PYBIGWIG_AVAILABLE = True
except ImportError:
    print("⚠️  pyBigWig not available. Install with: pip install pyBigWig")
    PYBIGWIG_AVAILABLE = False

print("="*80)
print("EXTRACT BIGWIG SIGNALS AND CREATE VISUALIZATIONS")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# ============================================================================
# Configuration
# ============================================================================
INPUT_DIR = "output/GABA_DEG_analysis"
BIGWIG_BASE = "../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR = "output/GABA_DEG_analysis/bigwig_signals"
PLOT_DIR = "output/GABA_DEG_analysis/plots"

GENOTYPES = ["Nestin", "Emx1"]
CELL_TYPE = "GABA"

# Create output directories
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)

# Plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 150

# ============================================================================
# STEP 1: Check prerequisites
# ============================================================================
print("STEP 1: Checking prerequisites...")
print("-"*80)

if not PYBIGWIG_AVAILABLE:
    print("\n⚠️  ERROR: pyBigWig not available!")
    print("Install with: pip install pyBigWig")
    print("or: conda install -c bioconda pybigwig")
    exit(1)

# Check input files
missing_files = []
for genotype in GENOTYPES:
    input_file = os.path.join(INPUT_DIR, f"DEG_CRE_links_{genotype}.tsv")
    if not os.path.exists(input_file):
        missing_files.append(input_file)

if missing_files:
    print("\n⚠️  ERROR: Missing input files:")
    for f in missing_files:
        print(f"  {f}")
    print("\nRun analyze_GABA_DEGs_with_literature_CREs.py first!")
    exit(1)

print("✓ All input files found")

# Check BigWig files
bigwig_files = {}
missing_bw = []

for genotype in GENOTYPES:
    for condition in ['Ctrl', 'Mut']:
        bw_file = os.path.join(BIGWIG_BASE, f"GABA_{genotype}-{condition}.bw")
        if os.path.exists(bw_file):
            bigwig_files[f"{genotype}_{condition}"] = bw_file
            print(f"✓ Found: GABA_{genotype}-{condition}.bw")
        else:
            missing_bw.append(f"GABA_{genotype}-{condition}.bw")

if missing_bw:
    print("\n⚠️  WARNING: Missing BigWig files:")
    for f in missing_bw:
        print(f"  {f}")
    print("\nWill only analyze genotypes with complete data.")

# ============================================================================
# STEP 2: Extract BigWig signals at CRE locations
# ============================================================================
print(f"\n{'='*80}")
print("STEP 2: Extracting BigWig signals at CRE locations...")
print("-"*80)

def extract_signal(bw_file, chrom, start, end):
    """Extract mean signal from BigWig file in given region"""
    try:
        bw = pyBigWig.open(bw_file)

        # Check if chromosome exists
        if chrom not in bw.chroms():
            bw.close()
            return np.nan

        # Get signal
        signal = bw.stats(chrom, int(start), int(end), type="mean")
        bw.close()

        if signal and signal[0] is not None:
            return signal[0]
        else:
            return 0.0
    except Exception as e:
        print(f"    Error extracting signal from {bw_file} at {chrom}:{start}-{end}: {e}")
        return np.nan

signal_results = {}

for genotype in GENOTYPES:
    print(f"\n{genotype}:")

    # Check if we have both BigWig files
    ctrl_key = f"{genotype}_Ctrl"
    mut_key = f"{genotype}_Mut"

    if ctrl_key not in bigwig_files or mut_key not in bigwig_files:
        print(f"  Skipping (missing BigWig files)")
        continue

    # Load DEG-CRE links
    input_file = os.path.join(INPUT_DIR, f"DEG_CRE_links_{genotype}.tsv")
    deg_cre = pd.read_csv(input_file, sep='\t')

    print(f"  Loading {len(deg_cre)} gene-CRE pairs...")

    # Extract signals
    signals = []

    for idx, row in deg_cre.iterrows():
        if idx % 10 == 0:
            print(f"    Processing {idx}/{len(deg_cre)}...", end='\r')

        chrom = row['cre_chr']
        start = row['cre_start']
        end = row['cre_end']

        ctrl_signal = extract_signal(bigwig_files[ctrl_key], chrom, start, end)
        mut_signal = extract_signal(bigwig_files[mut_key], chrom, start, end)

        signals.append({
            'gene': row['gene'],
            'cre_id': row['cre_id'],
            'cre_chr': chrom,
            'cre_start': start,
            'cre_end': end,
            'gene_log2fc': row['gene_log2fc'],
            'gene_direction': row['gene_direction'],
            'ctrl_signal': ctrl_signal,
            'mut_signal': mut_signal,
            'signal_log2fc': np.log2((mut_signal + 0.01) / (ctrl_signal + 0.01)),
            'has_da_peak': row.get('has_da_peak', False),
            'consistency': row.get('consistency', 'UNKNOWN')
        })

    print(f"    Processing {len(deg_cre)}/{len(deg_cre)}... Done!")

    signals_df = pd.DataFrame(signals)
    signal_results[genotype] = signals_df

    # Save signals
    output_file = os.path.join(OUTPUT_DIR, f"CRE_signals_{genotype}.tsv")
    signals_df.to_csv(output_file, sep='\t', index=False)
    print(f"  ✓ Saved signals: {output_file}")

    # Summary statistics
    valid_signals = signals_df.dropna(subset=['ctrl_signal', 'mut_signal'])
    print(f"  Valid signals: {len(valid_signals)} / {len(signals_df)}")
    print(f"    Mean Ctrl signal: {valid_signals['ctrl_signal'].mean():.3f}")
    print(f"    Mean Mut signal: {valid_signals['mut_signal'].mean():.3f}")

# ============================================================================
# STEP 3: Create visualizations
# ============================================================================
print(f"\n{'='*80}")
print("STEP 3: Creating visualizations...")
print("-"*80)

for genotype in GENOTYPES:
    if genotype not in signal_results:
        continue

    print(f"\n{genotype}:")
    signals_df = signal_results[genotype]

    # Filter for valid signals
    valid_signals = signals_df.dropna(subset=['ctrl_signal', 'mut_signal'])

    if len(valid_signals) == 0:
        print("  No valid signals to plot")
        continue

    # ========================================================================
    # Plot 1: Ctrl vs Mut signal scatter
    # ========================================================================
    fig, ax = plt.subplots(figsize=(8, 8))

    # Color by gene expression direction
    colors = valid_signals['gene_direction'].map({'Up': 'red', 'Down': 'blue'})

    ax.scatter(valid_signals['ctrl_signal'], valid_signals['mut_signal'],
               c=colors, alpha=0.6, s=50)

    # Add diagonal line
    max_val = max(valid_signals['ctrl_signal'].max(), valid_signals['mut_signal'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='No change')

    ax.set_xlabel('ATAC Signal (Control)', fontsize=12)
    ax.set_ylabel('ATAC Signal (Mutant)', fontsize=12)
    ax.set_title(f'CRE Accessibility: {genotype} Ctrl vs Mut\n'
                 f'(Red=Gene Up, Blue=Gene Down)', fontsize=14)
    ax.legend()

    plt.tight_layout()
    plot_file = os.path.join(PLOT_DIR, f"scatter_ctrl_vs_mut_{genotype}.png")
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved: scatter_ctrl_vs_mut_{genotype}.png")

    # ========================================================================
    # Plot 2: Signal log2FC distribution
    # ========================================================================
    fig, ax = plt.subplots(figsize=(10, 6))

    # Separate by gene expression direction
    up_genes = valid_signals[valid_signals['gene_direction'] == 'Up']['signal_log2fc']
    down_genes = valid_signals[valid_signals['gene_direction'] == 'Down']['signal_log2fc']

    ax.hist(up_genes, bins=30, alpha=0.6, color='red', label='Up-regulated genes')
    ax.hist(down_genes, bins=30, alpha=0.6, color='blue', label='Down-regulated genes')
    ax.axvline(0, color='black', linestyle='--', alpha=0.5, label='No change')

    ax.set_xlabel('CRE Accessibility log2FC (Mut/Ctrl)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'CRE Accessibility Changes: {genotype}', fontsize=14)
    ax.legend()

    plt.tight_layout()
    plot_file = os.path.join(PLOT_DIR, f"signal_log2fc_dist_{genotype}.png")
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved: signal_log2fc_dist_{genotype}.png")

    # ========================================================================
    # Plot 3: Gene expression vs CRE accessibility
    # ========================================================================
    fig, ax = plt.subplots(figsize=(10, 8))

    ax.scatter(valid_signals['gene_log2fc'], valid_signals['signal_log2fc'],
               c=colors, alpha=0.6, s=50)

    # Add quadrant lines
    ax.axhline(0, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(0, color='gray', linestyle='--', alpha=0.3)

    # Highlight consistent points
    consistent = valid_signals[valid_signals['consistency'] == 'CONSISTENT']
    has_consistent = len(consistent) > 0

    if has_consistent:
        ax.scatter(consistent['gene_log2fc'], consistent['signal_log2fc'],
                   edgecolors='green', facecolors='none', s=100, linewidths=2,
                   label='Consistent (DA peak)')

    ax.set_xlabel('Gene Expression log2FC (Mut/Ctrl)', fontsize=12)
    ax.set_ylabel('CRE Accessibility log2FC (Mut/Ctrl)', fontsize=12)

    if has_consistent:
        ax.set_title(f'Gene Expression vs CRE Accessibility: {genotype}\n'
                     f'(Green circles = Consistent with DA peaks)', fontsize=14)
        ax.legend()
    else:
        ax.set_title(f'Gene Expression vs CRE Accessibility: {genotype}\n'
                     f'(No DA peak overlaps found)', fontsize=14)

    plt.tight_layout()
    plot_file = os.path.join(PLOT_DIR, f"gene_vs_cre_{genotype}.png")
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved: gene_vs_cre_{genotype}.png")

    # ========================================================================
    # Plot 4: Top genes heatmap (if enough genes)
    # ========================================================================
    # Group by gene and take mean signal
    gene_signals = valid_signals.groupby('gene').agg({
        'ctrl_signal': 'mean',
        'mut_signal': 'mean',
        'gene_log2fc': 'first',
        'gene_direction': 'first'
    }).reset_index()

    # Sort by absolute gene log2FC
    gene_signals['abs_log2fc'] = gene_signals['gene_log2fc'].abs()
    gene_signals = gene_signals.sort_values('abs_log2fc', ascending=False)

    # Take top 20 genes
    top_genes = gene_signals.head(20)

    if len(top_genes) >= 5:
        fig, ax = plt.subplots(figsize=(10, max(8, len(top_genes) * 0.4)))

        # Create matrix
        data = top_genes[['ctrl_signal', 'mut_signal']].values

        # Plot heatmap
        im = ax.imshow(data, cmap='YlOrRd', aspect='auto')

        # Set ticks
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Control', 'Mutant'])
        ax.set_yticks(range(len(top_genes)))
        ax.set_yticklabels(top_genes['gene'].tolist())

        # Add colorbar
        plt.colorbar(im, ax=ax, label='Mean ATAC Signal')

        # Add gene expression direction
        for i, row in enumerate(top_genes.itertuples()):
            direction_color = 'red' if row.gene_direction == 'Up' else 'blue'
            ax.text(-0.5, i, '▲' if row.gene_direction == 'Up' else '▼',
                   ha='center', va='center', color=direction_color, fontsize=12)

        ax.set_title(f'Top {len(top_genes)} DEGs: Mean CRE Accessibility - {genotype}\n'
                     f'(▲=Up-regulated, ▼=Down-regulated)', fontsize=14)

        plt.tight_layout()
        plot_file = os.path.join(PLOT_DIR, f"top_genes_heatmap_{genotype}.png")
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: top_genes_heatmap_{genotype}.png")

# ============================================================================
# STEP 4: Create summary report
# ============================================================================
print(f"\n{'='*80}")
print("STEP 4: Creating summary report...")
print("-"*80)

report_file = os.path.join(OUTPUT_DIR, "bigwig_analysis_report.txt")
with open(report_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("BIGWIG SIGNAL EXTRACTION AND VISUALIZATION REPORT\n")
    f.write("="*80 + "\n\n")
    f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("RESULTS BY GENOTYPE:\n")
    f.write("-"*80 + "\n\n")

    for genotype in GENOTYPES:
        if genotype not in signal_results:
            f.write(f"{genotype}: No data\n\n")
            continue

        signals_df = signal_results[genotype]
        valid_signals = signals_df.dropna(subset=['ctrl_signal', 'mut_signal'])

        f.write(f"{genotype}:\n")
        f.write(f"  Total gene-CRE pairs: {len(signals_df)}\n")
        f.write(f"  Valid signals extracted: {len(valid_signals)}\n\n")

        if len(valid_signals) > 0:
            f.write(f"  Signal statistics:\n")
            f.write(f"    Mean Ctrl signal: {valid_signals['ctrl_signal'].mean():.3f}\n")
            f.write(f"    Mean Mut signal: {valid_signals['mut_signal'].mean():.3f}\n")
            f.write(f"    Mean log2FC: {valid_signals['signal_log2fc'].mean():.3f}\n\n")

            # Correlation with gene expression
            corr = valid_signals[['gene_log2fc', 'signal_log2fc']].corr().iloc[0, 1]
            f.write(f"  Correlation (gene expr vs CRE access): {corr:.3f}\n\n")

            # Consistency
            if 'consistency' in valid_signals.columns:
                n_consistent = (valid_signals['consistency'] == 'CONSISTENT').sum()
                n_inconsistent = (valid_signals['consistency'] == 'INCONSISTENT').sum()
                f.write(f"  DA peak consistency:\n")
                f.write(f"    Consistent: {n_consistent}\n")
                f.write(f"    Inconsistent: {n_inconsistent}\n\n")

        f.write("\n")

    f.write("="*80 + "\n")
    f.write("OUTPUT FILES:\n")
    f.write("="*80 + "\n")
    f.write(f"  Signals: {OUTPUT_DIR}/CRE_signals_{{genotype}}.tsv\n")
    f.write(f"  Plots: {PLOT_DIR}/\n")
    f.write(f"    - scatter_ctrl_vs_mut_{{genotype}}.png\n")
    f.write(f"    - signal_log2fc_dist_{{genotype}}.png\n")
    f.write(f"    - gene_vs_cre_{{genotype}}.png\n")
    f.write(f"    - top_genes_heatmap_{{genotype}}.png\n\n")

    f.write("INTERPRETATION GUIDE:\n")
    f.write("-"*80 + "\n")
    f.write("1. Scatter plot (Ctrl vs Mut): Points below diagonal = decreased accessibility\n")
    f.write("2. Distribution plot: Look for shift in up/down-regulated genes\n")
    f.write("3. Gene vs CRE plot: Points in quadrants 1 & 3 = consistent changes\n")
    f.write("4. Heatmap: Visual comparison of top DEGs' CRE accessibility\n\n")

    f.write("HYPOTHESIS GENERATION:\n")
    f.write("-"*80 + "\n")
    f.write("- Look for genes where CRE accessibility changes match expression changes\n")
    f.write("- Green circles in gene_vs_cre plot = supported by DA peaks\n")
    f.write("- Consider genes with multiple consistent CREs for follow-up\n")

print(f"✓ Saved: {report_file}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print(f"\n{'='*80}")
print("ANALYSIS COMPLETE!")
print("="*80)

for genotype in GENOTYPES:
    if genotype in signal_results:
        signals_df = signal_results[genotype]
        valid = signals_df.dropna(subset=['ctrl_signal', 'mut_signal'])

        print(f"\n{genotype}:")
        print(f"  Signals extracted: {len(valid)} / {len(signals_df)}")

        if len(valid) > 0:
            corr = valid[['gene_log2fc', 'signal_log2fc']].corr().iloc[0, 1]
            print(f"  Correlation (gene vs CRE): {corr:.3f}")

print(f"\nOutput files:")
print(f"  Signals: {OUTPUT_DIR}/")
print(f"  Plots: {PLOT_DIR}/")

print(f"\nNext steps:")
print("  1. Review plots for overall trends")
print("  2. Identify genes with strong CRE support")
print("  3. Prioritize candidates for experimental validation")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
