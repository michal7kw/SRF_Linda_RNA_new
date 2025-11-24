#!/usr/bin/env python3
"""
Extract BigWig Signals at Literature CREs and Create Visualizations (CORRECTED)

Purpose:
- For genes with literature CREs, extract ATAC signal from BigWig files
- Compare Ctrl vs Mut signal at each CRE
- Create visual comparisons for hypothesis generation
- CORRECTED: Proper handling of low-signal CREs and pseudocounts

Changes from original:
1. Adaptive pseudocount based on median signal
2. Filter out very low-signal CREs (likely noise)
3. Cap extreme log2FC values at biologically reasonable thresholds
4. Add diagnostic plots showing signal distributions
5. Add quality control metrics

This script should be run AFTER analyze_GABA_DEGs_with_literature_CREs.py

INPUT FILES:
- output/GABA_DEG_analysis/DEG_CRE_links_*.tsv: DEG-CRE linkage files for each genotype
- ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Nestin-Ctrl.bw: ATAC signal for Nestin control
- ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Nestin-Mut.bw: ATAC signal for Nestin mutant
- ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Emx1-Ctrl.bw: ATAC signal for Emx1 control
- ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Emx1-Mut.bw: ATAC signal for Emx1 mutant

OUTPUT FILES:
- output/GABA_DEG_analysis/bigwig_signals_corrected/CRE_signals_*.tsv: Filtered signal data for each genotype
- output/GABA_DEG_analysis/plots_corrected/diagnostic_qc_*.png: Quality control diagnostic plots
- output/GABA_DEG_analysis/plots_corrected/scatter_ctrl_vs_mut_*.png: Signal comparison scatter plots
- output/GABA_DEG_analysis/plots_corrected/signal_log2fc_dist_*.png: Log2FC distribution plots
- output/GABA_DEG_analysis/plots_corrected/gene_vs_cre_*.png: Gene expression vs CRE accessibility plots
- output/GABA_DEG_analysis/plots_corrected/top_genes_heatmap_*.png: Heatmaps of top DEGs
- output/GABA_DEG_analysis/bigwig_signals_corrected/bigwig_analysis_report_corrected.txt: Complete analysis report
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
print("EXTRACT BIGWIG SIGNALS AND CREATE VISUALIZATIONS (CORRECTED)")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# ============================================================================
# Configuration
# ============================================================================
INPUT_DIR = "output/GABA_DEG_analysis"
BIGWIG_BASE = "../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR = "output/GABA_DEG_analysis/bigwig_signals_corrected"
PLOT_DIR = "output/GABA_DEG_analysis/plots_corrected"

GENOTYPES = ["Nestin", "Emx1"]
CELL_TYPE = "GABA"

# Quality control parameters
MIN_SIGNAL_THRESHOLD = 0.05  # Minimum signal in either Ctrl or Mut to be considered
MAX_LOG2FC = 3.0             # Cap extreme fold-changes at ±3 (8-fold change)
PSEUDOCOUNT_PERCENTILE = 10  # Use 10th percentile of non-zero signals as pseudocount

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
raw_signal_results = {}  # Keep raw data for diagnostics

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
            'has_da_peak': row.get('has_da_peak', False),
            'consistency': row.get('consistency', 'UNKNOWN')
        })

    print(f"    Processing {len(deg_cre)}/{len(deg_cre)}... Done!")

    signals_df = pd.DataFrame(signals)
    raw_signal_results[genotype] = signals_df.copy()

    # ========================================================================
    # QUALITY CONTROL: Calculate adaptive pseudocount and filter
    # ========================================================================
    print(f"\n  Quality control:")

    valid_initial = signals_df.dropna(subset=['ctrl_signal', 'mut_signal'])
    print(f"    Initial valid signals: {len(valid_initial)} / {len(signals_df)}")

    # Calculate adaptive pseudocount from non-zero signals
    all_nonzero_signals = pd.concat([
        valid_initial[valid_initial['ctrl_signal'] > 0]['ctrl_signal'],
        valid_initial[valid_initial['mut_signal'] > 0]['mut_signal']
    ])

    if len(all_nonzero_signals) > 0:
        pseudocount = np.percentile(all_nonzero_signals, PSEUDOCOUNT_PERCENTILE)
        print(f"    Adaptive pseudocount (p{PSEUDOCOUNT_PERCENTILE}): {pseudocount:.4f}")
    else:
        pseudocount = 0.01
        print(f"    Using default pseudocount: {pseudocount}")

    # Filter for minimum signal threshold
    signals_df['max_signal'] = signals_df[['ctrl_signal', 'mut_signal']].max(axis=1)
    signals_filtered = signals_df[signals_df['max_signal'] >= MIN_SIGNAL_THRESHOLD].copy()

    print(f"    After filtering (signal >= {MIN_SIGNAL_THRESHOLD}): {len(signals_filtered)} / {len(signals_df)}")
    print(f"    Removed low-signal CREs: {len(signals_df) - len(signals_filtered)}")

    # Calculate log2FC with adaptive pseudocount
    signals_filtered['signal_log2fc_raw'] = np.log2(
        (signals_filtered['mut_signal'] + pseudocount) /
        (signals_filtered['ctrl_signal'] + pseudocount)
    )

    # Cap extreme values
    signals_filtered['signal_log2fc'] = signals_filtered['signal_log2fc_raw'].clip(
        lower=-MAX_LOG2FC,
        upper=MAX_LOG2FC
    )

    n_capped = (signals_filtered['signal_log2fc_raw'].abs() > MAX_LOG2FC).sum()
    print(f"    Capped extreme log2FC (|fc| > {MAX_LOG2FC}): {n_capped}")

    signal_results[genotype] = signals_filtered

    # Save filtered signals
    output_file = os.path.join(OUTPUT_DIR, f"CRE_signals_{genotype}.tsv")
    signals_filtered.to_csv(output_file, sep='\t', index=False)
    print(f"  ✓ Saved filtered signals: {output_file}")

    # Summary statistics
    print(f"\n  Signal statistics (filtered data):")
    print(f"    Mean Ctrl signal: {signals_filtered['ctrl_signal'].mean():.3f}")
    print(f"    Mean Mut signal: {signals_filtered['mut_signal'].mean():.3f}")
    print(f"    Median Ctrl signal: {signals_filtered['ctrl_signal'].median():.3f}")
    print(f"    Median Mut signal: {signals_filtered['mut_signal'].median():.3f}")
    print(f"    Mean log2FC: {signals_filtered['signal_log2fc'].mean():.3f}")
    print(f"    Median log2FC: {signals_filtered['signal_log2fc'].median():.3f}")

# ============================================================================
# STEP 3: Create diagnostic plots
# ============================================================================
print(f"\n{'='*80}")
print("STEP 3: Creating diagnostic plots...")
print("-"*80)

for genotype in GENOTYPES:
    if genotype not in signal_results:
        continue

    print(f"\n{genotype}:")
    raw_signals = raw_signal_results[genotype]
    filtered_signals = signal_results[genotype]

    # ========================================================================
    # Diagnostic Plot: Signal distribution before/after filtering
    # ========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Signal Distribution QC: {genotype}', fontsize=16, fontweight='bold')

    # Raw signal histograms
    ax = axes[0, 0]
    valid_raw = raw_signals.dropna(subset=['ctrl_signal', 'mut_signal'])
    ax.hist(valid_raw['ctrl_signal'], bins=50, alpha=0.7, color='blue', label='Control')
    ax.hist(valid_raw['mut_signal'], bins=50, alpha=0.7, color='red', label='Mutant')
    ax.axvline(MIN_SIGNAL_THRESHOLD, color='black', linestyle='--',
               label=f'Filter threshold ({MIN_SIGNAL_THRESHOLD})')
    ax.set_xlabel('ATAC Signal')
    ax.set_ylabel('Count')
    ax.set_title('Raw Signal Distribution (All CREs)')
    ax.set_xlim(0, 2)
    ax.legend()

    # Log-scale signal histograms
    ax = axes[0, 1]
    valid_nonzero = valid_raw[(valid_raw['ctrl_signal'] > 0) & (valid_raw['mut_signal'] > 0)]
    ax.hist(np.log10(valid_nonzero['ctrl_signal'] + 1e-10), bins=50, alpha=0.7,
            color='blue', label='Control')
    ax.hist(np.log10(valid_nonzero['mut_signal'] + 1e-10), bins=50, alpha=0.7,
            color='red', label='Mutant')
    ax.axvline(np.log10(MIN_SIGNAL_THRESHOLD), color='black', linestyle='--',
               label=f'Filter threshold')
    ax.set_xlabel('log10(ATAC Signal)')
    ax.set_ylabel('Count')
    ax.set_title('Log-Scale Signal Distribution')
    ax.legend()

    # Filtered signal scatter
    ax = axes[1, 0]
    ax.scatter(filtered_signals['ctrl_signal'], filtered_signals['mut_signal'],
               alpha=0.4, s=20, c='gray')
    max_val = max(filtered_signals['ctrl_signal'].max(),
                  filtered_signals['mut_signal'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='No change')
    ax.set_xlabel('Control Signal')
    ax.set_ylabel('Mutant Signal')
    ax.set_title(f'Filtered CREs (n={len(filtered_signals)})')
    ax.legend()

    # Log2FC distribution comparison
    ax = axes[1, 1]
    # Raw log2FC (from unfiltered data)
    valid_raw_with_fc = valid_raw.copy()
    valid_raw_with_fc['signal_log2fc_raw'] = np.log2(
        (valid_raw_with_fc['mut_signal'] + 0.01) /
        (valid_raw_with_fc['ctrl_signal'] + 0.01)
    )

    ax.hist(valid_raw_with_fc['signal_log2fc_raw'], bins=50, alpha=0.5,
            color='gray', label='Original (with artifacts)', range=(-10, 10))
    ax.hist(filtered_signals['signal_log2fc'], bins=50, alpha=0.7,
            color='green', label='Corrected (filtered & capped)')
    ax.axvline(0, color='black', linestyle='--', alpha=0.5)
    ax.axvline(-MAX_LOG2FC, color='red', linestyle=':', alpha=0.5, label=f'Cap at ±{MAX_LOG2FC}')
    ax.axvline(MAX_LOG2FC, color='red', linestyle=':', alpha=0.5)
    ax.set_xlabel('CRE Accessibility log2FC (Mut/Ctrl)')
    ax.set_ylabel('Count')
    ax.set_title('Log2FC Distribution: Before vs After Correction')
    ax.legend()

    plt.tight_layout()
    plot_file = os.path.join(PLOT_DIR, f"diagnostic_qc_{genotype}.png")
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved: diagnostic_qc_{genotype}.png")

# ============================================================================
# STEP 4: Create main visualizations (corrected data)
# ============================================================================
print(f"\n{'='*80}")
print("STEP 4: Creating main visualizations (corrected data)...")
print("-"*80)

for genotype in GENOTYPES:
    if genotype not in signal_results:
        continue

    print(f"\n{genotype}:")
    valid_signals = signal_results[genotype]

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
    ax.set_title(f'CRE Accessibility: {genotype} Ctrl vs Mut (Corrected)\n'
                 f'(Red=Gene Up, Blue=Gene Down) | n={len(valid_signals)} CREs', fontsize=14)
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

    ax.hist(up_genes, bins=40, alpha=0.6, color='red', label=f'Up-regulated genes (n={len(up_genes)})')
    ax.hist(down_genes, bins=40, alpha=0.6, color='blue', label=f'Down-regulated genes (n={len(down_genes)})')
    ax.axvline(0, color='black', linestyle='--', alpha=0.5, label='No change')
    ax.axvline(valid_signals['signal_log2fc'].median(), color='green',
               linestyle='-', alpha=0.7, linewidth=2,
               label=f'Median log2FC: {valid_signals["signal_log2fc"].median():.2f}')

    ax.set_xlabel('CRE Accessibility log2FC (Mut/Ctrl)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'CRE Accessibility Changes: {genotype} (Corrected)', fontsize=14)
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

    # Calculate correlation
    corr = valid_signals[['gene_log2fc', 'signal_log2fc']].corr().iloc[0, 1]

    # Add trend line
    from scipy import stats
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        valid_signals['gene_log2fc'], valid_signals['signal_log2fc']
    )
    x_trend = np.array([valid_signals['gene_log2fc'].min(),
                        valid_signals['gene_log2fc'].max()])
    y_trend = slope * x_trend + intercept
    ax.plot(x_trend, y_trend, 'g--', alpha=0.5, linewidth=2,
            label=f'Trend: r={corr:.3f}, p={p_value:.2e}')

    ax.set_xlabel('Gene Expression log2FC (Mut/Ctrl)', fontsize=12)
    ax.set_ylabel('CRE Accessibility log2FC (Mut/Ctrl)', fontsize=12)

    if has_consistent:
        ax.set_title(f'Gene Expression vs CRE Accessibility: {genotype} (Corrected)\n'
                     f'(Green circles = Consistent with DA peaks) | n={len(valid_signals)} CREs',
                     fontsize=14)
    else:
        ax.set_title(f'Gene Expression vs CRE Accessibility: {genotype} (Corrected)\n'
                     f'(No DA peak overlaps found) | n={len(valid_signals)} CREs',
                     fontsize=14)

    ax.legend()

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

        ax.set_title(f'Top {len(top_genes)} DEGs: Mean CRE Accessibility - {genotype} (Corrected)\n'
                     f'(▲=Up-regulated, ▼=Down-regulated)', fontsize=14)

        plt.tight_layout()
        plot_file = os.path.join(PLOT_DIR, f"top_genes_heatmap_{genotype}.png")
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: top_genes_heatmap_{genotype}.png")

# ============================================================================
# STEP 5: Create summary report
# ============================================================================
print(f"\n{'='*80}")
print("STEP 5: Creating summary report...")
print("-"*80)

report_file = os.path.join(OUTPUT_DIR, "bigwig_analysis_report_corrected.txt")
with open(report_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("BIGWIG SIGNAL EXTRACTION AND VISUALIZATION REPORT (CORRECTED)\n")
    f.write("="*80 + "\n\n")
    f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("CORRECTIONS APPLIED:\n")
    f.write("-"*80 + "\n")
    f.write(f"1. Minimum signal threshold: {MIN_SIGNAL_THRESHOLD}\n")
    f.write(f"   - Filters out low-signal CREs (likely noise/artifacts)\n\n")
    f.write(f"2. Adaptive pseudocount: {PSEUDOCOUNT_PERCENTILE}th percentile of non-zero signals\n")
    f.write(f"   - Avoids artifacts from fixed pseudocount dominating low signals\n\n")
    f.write(f"3. Maximum log2FC cap: ±{MAX_LOG2FC}\n")
    f.write(f"   - Caps extreme fold-changes at biologically reasonable levels\n")
    f.write(f"   - ±3 = 8-fold change\n\n")

    f.write("RESULTS BY GENOTYPE:\n")
    f.write("-"*80 + "\n\n")

    for genotype in GENOTYPES:
        if genotype not in signal_results:
            f.write(f"{genotype}: No data\n\n")
            continue

        raw_signals = raw_signal_results[genotype]
        filtered_signals = signal_results[genotype]

        f.write(f"{genotype}:\n")
        f.write(f"  Total gene-CRE pairs (raw): {len(raw_signals)}\n")
        f.write(f"  Filtered pairs (signal >= {MIN_SIGNAL_THRESHOLD}): {len(filtered_signals)}\n")
        f.write(f"  Filtering rate: {100 * (1 - len(filtered_signals)/len(raw_signals)):.1f}%\n\n")

        if len(filtered_signals) > 0:
            f.write(f"  Signal statistics (filtered data):\n")
            f.write(f"    Mean Ctrl signal: {filtered_signals['ctrl_signal'].mean():.3f}\n")
            f.write(f"    Mean Mut signal: {filtered_signals['mut_signal'].mean():.3f}\n")
            f.write(f"    Median Ctrl signal: {filtered_signals['ctrl_signal'].median():.3f}\n")
            f.write(f"    Median Mut signal: {filtered_signals['mut_signal'].median():.3f}\n\n")

            f.write(f"  Log2FC statistics:\n")
            f.write(f"    Mean log2FC: {filtered_signals['signal_log2fc'].mean():.3f}\n")
            f.write(f"    Median log2FC: {filtered_signals['signal_log2fc'].median():.3f}\n")
            f.write(f"    Std log2FC: {filtered_signals['signal_log2fc'].std():.3f}\n\n")

            # Correlation with gene expression
            corr = filtered_signals[['gene_log2fc', 'signal_log2fc']].corr().iloc[0, 1]
            f.write(f"  Correlation (gene expr vs CRE access): {corr:.3f}\n\n")

            # Direction analysis
            n_decreased = (filtered_signals['signal_log2fc'] < -0.5).sum()
            n_increased = (filtered_signals['signal_log2fc'] > 0.5).sum()
            n_unchanged = ((filtered_signals['signal_log2fc'] >= -0.5) &
                          (filtered_signals['signal_log2fc'] <= 0.5)).sum()

            f.write(f"  CRE accessibility changes:\n")
            f.write(f"    Decreased (log2FC < -0.5): {n_decreased} ({100*n_decreased/len(filtered_signals):.1f}%)\n")
            f.write(f"    Unchanged (-0.5 to 0.5): {n_unchanged} ({100*n_unchanged/len(filtered_signals):.1f}%)\n")
            f.write(f"    Increased (log2FC > 0.5): {n_increased} ({100*n_increased/len(filtered_signals):.1f}%)\n\n")

            # Consistency
            if 'consistency' in filtered_signals.columns:
                n_consistent = (filtered_signals['consistency'] == 'CONSISTENT').sum()
                n_inconsistent = (filtered_signals['consistency'] == 'INCONSISTENT').sum()
                f.write(f"  DA peak consistency:\n")
                f.write(f"    Consistent: {n_consistent}\n")
                f.write(f"    Inconsistent: {n_inconsistent}\n\n")

        f.write("\n")

    f.write("="*80 + "\n")
    f.write("OUTPUT FILES:\n")
    f.write("="*80 + "\n")
    f.write(f"  Filtered signals: {OUTPUT_DIR}/CRE_signals_{{genotype}}.tsv\n")
    f.write(f"  Plots: {PLOT_DIR}/\n")
    f.write(f"    - diagnostic_qc_{{genotype}}.png (QC diagnostics)\n")
    f.write(f"    - scatter_ctrl_vs_mut_{{genotype}}.png\n")
    f.write(f"    - signal_log2fc_dist_{{genotype}}.png\n")
    f.write(f"    - gene_vs_cre_{{genotype}}.png\n")
    f.write(f"    - top_genes_heatmap_{{genotype}}.png\n\n")

    f.write("INTERPRETATION GUIDE:\n")
    f.write("-"*80 + "\n")
    f.write("1. Diagnostic QC plot: Shows filtering effects and removal of artifacts\n")
    f.write("2. Scatter plot (Ctrl vs Mut): Points below diagonal = decreased accessibility\n")
    f.write("3. Distribution plot: Shift towards negative = chromatin compaction\n")
    f.write("4. Gene vs CRE plot: Points in quadrants 1 & 3 = consistent changes\n")
    f.write("5. Heatmap: Visual comparison of top DEGs' CRE accessibility\n\n")

    f.write("KEY IMPROVEMENTS OVER ORIGINAL:\n")
    f.write("-"*80 + "\n")
    f.write("- Removed low-signal artifacts (binary ON/OFF CREs)\n")
    f.write("- Used adaptive pseudocount to avoid dominating low signals\n")
    f.write("- Capped extreme fold-changes at biologically reasonable levels\n")
    f.write("- Added quality control metrics and diagnostic plots\n")
    f.write("- More accurate representation of gradual chromatin changes\n")

print(f"✓ Saved: {report_file}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print(f"\n{'='*80}")
print("ANALYSIS COMPLETE! (CORRECTED VERSION)")
print("="*80)

for genotype in GENOTYPES:
    if genotype in signal_results:
        raw_df = raw_signal_results[genotype]
        filtered_df = signal_results[genotype]

        print(f"\n{genotype}:")
        print(f"  Raw CRE-gene pairs: {len(raw_df)}")
        print(f"  Filtered pairs (signal >= {MIN_SIGNAL_THRESHOLD}): {len(filtered_df)}")
        print(f"  Filtering rate: {100 * (1 - len(filtered_df)/len(raw_df)):.1f}%")

        if len(filtered_df) > 0:
            corr = filtered_df[['gene_log2fc', 'signal_log2fc']].corr().iloc[0, 1]
            print(f"  Correlation (gene vs CRE): {corr:.3f}")
            print(f"  Median CRE log2FC: {filtered_df['signal_log2fc'].median():.3f}")

print(f"\nOutput files:")
print(f"  Filtered signals: {OUTPUT_DIR}/")
print(f"  Corrected plots: {PLOT_DIR}/")

print(f"\nKey improvements:")
print(f"  ✓ Removed {100 * (1 - len(signal_results['Emx1'])/len(raw_signal_results['Emx1'])):.1f}% low-signal artifacts")
print(f"  ✓ Adaptive pseudocount prevents artifacts")
print(f"  ✓ Capped extreme fold-changes at ±{MAX_LOG2FC}")
print(f"  ✓ Added diagnostic QC plots")

print(f"\nNext steps:")
print("  1. Compare diagnostic_qc plots (before/after filtering)")
print("  2. Review corrected gene_vs_cre plots")
print("  3. Check if correlation improves with filtered data")
print("  4. Focus on high-signal CREs for biological interpretation")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
