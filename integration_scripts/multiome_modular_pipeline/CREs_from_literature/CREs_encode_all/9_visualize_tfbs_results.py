#!/usr/bin/env python3
"""
9_visualize_tfbs_results.py

PURPOSE:
Create publication-quality visualizations of TFBS motif enrichment analysis
at promoter-like ATAC peaks.

DESCRIPTION:
This script generates:
1. Heatmaps of TF motif enrichment across conditions
2. Volcano plots of differential motif enrichment
3. Bar plots of top enriched TF families
4. Motif enrichment comparison plots
5. Summary figures for publication

INPUT:
- Differential analysis results from Step 8
- TF enrichment matrix

OUTPUT:
- PNG/PDF figures for publication
- Interactive HTML plots (optional)

AUTHOR: Claude Code
DATE: December 2024
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import logging

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for cluster
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm

# =============================================================================
# CONFIGURATION
# =============================================================================

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Default paths
SCRIPT_DIR = Path(__file__).parent.resolve()
INPUT_DIR = SCRIPT_DIR / "output" / "tfbs_analysis" / "differential_analysis"
OUTPUT_DIR = SCRIPT_DIR / "output" / "tfbs_analysis" / "figures"

# Plot styling
plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'sans-serif'
})

# Color palettes
CONDITION_COLORS = {
    'Nestin_Ctrl': '#2ecc71',   # Green
    'Nestin_Mut': '#e74c3c',    # Red
    'Emx1_Mut': '#9b59b6'       # Purple
}

# Thresholds for visualization
P_THRESHOLD = 0.01
LOG_P_THRESHOLD = 2.0  # -log10(0.01)
FC_THRESHOLD = 1.5


# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def plot_enrichment_heatmap(
    matrix_df: pd.DataFrame,
    output_dir: Path,
    top_n: int = 50
):
    """
    Create a heatmap of TF motif enrichment across conditions.
    """
    logger.info(f"Creating enrichment heatmap (top {top_n} motifs)...")

    # Get enrichment columns
    enr_cols = [c for c in matrix_df.columns if c.endswith('_LogP')]
    if len(enr_cols) == 0:
        logger.warning("No enrichment columns found")
        return

    # Calculate mean enrichment for ranking
    matrix_df['Mean_LogP'] = matrix_df[enr_cols].mean(axis=1)
    top_df = matrix_df.nlargest(top_n, 'Mean_LogP')

    # Prepare heatmap data
    heatmap_data = top_df[enr_cols].copy()
    heatmap_data.index = top_df['TF_Name']
    heatmap_data.columns = [c.replace('_LogP', '') for c in heatmap_data.columns]

    # Create figure
    fig, ax = plt.subplots(figsize=(8, max(12, top_n * 0.25)))

    # Create heatmap
    sns.heatmap(
        heatmap_data,
        cmap='YlOrRd',
        annot=False,
        fmt='.1f',
        linewidths=0.5,
        linecolor='white',
        cbar_kws={'label': '-log10(p-value)', 'shrink': 0.6},
        ax=ax
    )

    ax.set_title(f'TF Motif Enrichment at Promoter-Like ATAC Peaks\n(Top {top_n} by mean enrichment)', pad=15)
    ax.set_xlabel('Condition')
    ax.set_ylabel('Transcription Factor')

    # Rotate x labels
    plt.xticks(rotation=45, ha='right')

    # Save
    output_file = output_dir / 'heatmap_tf_enrichment.png'
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    logger.info(f"  Saved: {output_file}")

    # Also save as PDF for publication
    pdf_file = output_dir / 'heatmap_tf_enrichment.pdf'
    fig, ax = plt.subplots(figsize=(8, max(12, top_n * 0.25)))
    sns.heatmap(heatmap_data, cmap='YlOrRd', annot=False, linewidths=0.5, linecolor='white',
                cbar_kws={'label': '-log10(p-value)', 'shrink': 0.6}, ax=ax)
    ax.set_title(f'TF Motif Enrichment at Promoter-Like ATAC Peaks\n(Top {top_n} by mean enrichment)', pad=15)
    ax.set_xlabel('Condition')
    ax.set_ylabel('Transcription Factor')
    plt.xticks(rotation=45, ha='right')
    plt.savefig(pdf_file, bbox_inches='tight')
    plt.close()


def plot_volcano(
    diff_df: pd.DataFrame,
    comparison_name: str,
    output_dir: Path
):
    """
    Create a volcano plot of differential motif enrichment.
    """
    logger.info(f"Creating volcano plot for {comparison_name}...")

    if len(diff_df) == 0:
        logger.warning("No data for volcano plot")
        return

    # Parse comparison name
    parts = comparison_name.replace('differential_motifs_', '').replace('.tsv', '').split('_vs_')
    name1, name2 = parts[0], parts[1]

    # Get relevant columns
    x_col = 'Log2_Enrichment_Ratio'
    y_col = f'Log_P_value_{name2}'  # Use mutant p-value for y-axis

    if y_col not in diff_df.columns:
        y_col = [c for c in diff_df.columns if 'Log_P_value' in c][0]

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Color based on differential status
    colors = []
    for _, row in diff_df.iterrows():
        if row.get('Differential_Status', '') == f'Enriched_in_{name2}':
            colors.append('#e74c3c')  # Red for mutant-enriched
        elif row.get('Differential_Status', '') == f'Enriched_in_{name1}':
            colors.append('#2ecc71')  # Green for control-enriched
        else:
            colors.append('#95a5a6')  # Gray for no change

    # Scatter plot
    scatter = ax.scatter(
        diff_df[x_col],
        diff_df[y_col],
        c=colors,
        alpha=0.6,
        s=30,
        edgecolors='none'
    )

    # Add threshold lines
    ax.axhline(y=LOG_P_THRESHOLD, color='gray', linestyle='--', alpha=0.5, label=f'p=0.01')
    ax.axvline(x=np.log2(FC_THRESHOLD), color='gray', linestyle=':', alpha=0.5)
    ax.axvline(x=-np.log2(FC_THRESHOLD), color='gray', linestyle=':', alpha=0.5)

    # Label top differential motifs
    sig_df = diff_df[diff_df['Differential_Status'] != 'No_Change'].copy()
    if len(sig_df) > 0:
        top_sig = sig_df.nlargest(10, x_col, key=abs)
        for _, row in top_sig.iterrows():
            ax.annotate(
                row['TF_Name'],
                (row[x_col], row[y_col]),
                fontsize=8,
                alpha=0.8,
                xytext=(5, 5),
                textcoords='offset points'
            )

    ax.set_xlabel(f'Log2 Enrichment Ratio ({name2} / {name1})')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(f'Differential TF Motif Enrichment\n{name1} vs {name2}')

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#e74c3c', label=f'Enriched in {name2}'),
        Patch(facecolor='#2ecc71', label=f'Enriched in {name1}'),
        Patch(facecolor='#95a5a6', label='No significant change')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    # Save
    output_file = output_dir / f'volcano_{name1}_vs_{name2}.png'
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    logger.info(f"  Saved: {output_file}")


def plot_tf_family_enrichment(
    matrix_df: pd.DataFrame,
    output_dir: Path
):
    """
    Create a bar plot of TF family enrichment.
    """
    logger.info("Creating TF family enrichment plot...")

    # Get enrichment columns
    enr_cols = [c for c in matrix_df.columns if c.endswith('_LogP')]

    # Group by TF family
    family_data = []
    for family in matrix_df['TF_Family'].unique():
        family_df = matrix_df[matrix_df['TF_Family'] == family]
        row = {'TF_Family': family, 'Count': len(family_df)}
        for col in enr_cols:
            row[col.replace('_LogP', '_Mean_LogP')] = family_df[col].mean()
        family_data.append(row)

    family_df = pd.DataFrame(family_data)
    family_df = family_df.sort_values('Count', ascending=False).head(20)

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))

    # Left: Count of motifs per family
    ax1 = axes[0]
    bars = ax1.barh(family_df['TF_Family'], family_df['Count'], color='steelblue')
    ax1.set_xlabel('Number of Motifs')
    ax1.set_ylabel('TF Family')
    ax1.set_title('TF Family Representation')
    ax1.invert_yaxis()

    # Right: Mean enrichment per family by condition
    ax2 = axes[1]
    mean_cols = [c for c in family_df.columns if '_Mean_LogP' in c]
    x = np.arange(len(family_df))
    width = 0.25

    for i, col in enumerate(mean_cols):
        condition = col.replace('_Mean_LogP', '')
        color = CONDITION_COLORS.get(condition, 'gray')
        ax2.barh(x + i * width, family_df[col], width, label=condition, color=color, alpha=0.8)

    ax2.set_xlabel('Mean -log10(p-value)')
    ax2.set_ylabel('TF Family')
    ax2.set_title('Mean Enrichment by TF Family')
    ax2.set_yticks(x + width)
    ax2.set_yticklabels(family_df['TF_Family'])
    ax2.legend(loc='lower right')
    ax2.invert_yaxis()

    plt.tight_layout()

    # Save
    output_file = output_dir / 'barplot_tf_family_enrichment.png'
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    logger.info(f"  Saved: {output_file}")


def plot_condition_comparison(
    matrix_df: pd.DataFrame,
    output_dir: Path
):
    """
    Create scatter plots comparing enrichment between conditions.
    """
    logger.info("Creating condition comparison plots...")

    enr_cols = [c for c in matrix_df.columns if c.endswith('_LogP')]
    conditions = [c.replace('_LogP', '') for c in enr_cols]

    # Create pairwise comparisons
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    comparisons = [
        ('Nestin_Ctrl', 'Nestin_Mut', 'Nestin: Ctrl vs Mut'),
        ('Nestin_Ctrl', 'Emx1_Mut', 'Ctrl (Nestin) vs Mut (Emx1)'),
        ('Nestin_Mut', 'Emx1_Mut', 'Mut: Nestin vs Emx1')
    ]

    for idx, (cond1, cond2, title) in enumerate(comparisons):
        ax = axes[idx]

        col1 = f'{cond1}_LogP'
        col2 = f'{cond2}_LogP'

        if col1 not in matrix_df.columns or col2 not in matrix_df.columns:
            ax.set_visible(False)
            continue

        # Scatter plot
        ax.scatter(
            matrix_df[col1],
            matrix_df[col2],
            alpha=0.5,
            s=20,
            c='steelblue'
        )

        # Add diagonal line
        max_val = max(matrix_df[col1].max(), matrix_df[col2].max())
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='y=x')

        # Calculate correlation
        corr = matrix_df[[col1, col2]].corr().iloc[0, 1]
        ax.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax.transAxes,
                fontsize=10, verticalalignment='top')

        ax.set_xlabel(f'{cond1} -log10(p)')
        ax.set_ylabel(f'{cond2} -log10(p)')
        ax.set_title(title)

    plt.tight_layout()

    # Save
    output_file = output_dir / 'scatter_condition_comparisons.png'
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    logger.info(f"  Saved: {output_file}")


def plot_summary_figure(
    matrix_df: pd.DataFrame,
    diff_results: dict,
    output_dir: Path
):
    """
    Create a multi-panel summary figure for publication.
    """
    logger.info("Creating summary figure...")

    fig = plt.figure(figsize=(16, 12))

    # Panel A: Top enriched TFs heatmap (small version)
    ax1 = fig.add_subplot(2, 2, 1)
    enr_cols = [c for c in matrix_df.columns if c.endswith('_LogP')]
    matrix_df['Mean_LogP'] = matrix_df[enr_cols].mean(axis=1)
    top_df = matrix_df.nlargest(25, 'Mean_LogP')
    heatmap_data = top_df[enr_cols].copy()
    heatmap_data.index = top_df['TF_Name']
    heatmap_data.columns = [c.replace('_LogP', '') for c in heatmap_data.columns]

    sns.heatmap(heatmap_data, cmap='YlOrRd', annot=False, ax=ax1,
                cbar_kws={'shrink': 0.6, 'label': '-log10(p)'})
    ax1.set_title('A. Top 25 Enriched TF Motifs', fontweight='bold')
    ax1.set_xlabel('')
    ax1.set_ylabel('')

    # Panel B: TF family counts
    ax2 = fig.add_subplot(2, 2, 2)
    family_counts = matrix_df['TF_Family'].value_counts().head(15)
    family_counts.plot(kind='barh', ax=ax2, color='steelblue')
    ax2.set_xlabel('Number of Motifs')
    ax2.set_title('B. TF Family Distribution', fontweight='bold')
    ax2.invert_yaxis()

    # Panel C: Differential enrichment summary
    ax3 = fig.add_subplot(2, 2, 3)
    diff_summary = []
    for name, df in diff_results.items():
        if df is not None and 'Differential_Status' in df.columns:
            parts = name.split('_vs_')
            name1, name2 = parts[0], parts[1]
            up = (df['Differential_Status'] == f'Enriched_in_{name2}').sum()
            down = (df['Differential_Status'] == f'Enriched_in_{name1}').sum()
            diff_summary.append({
                'Comparison': f'{name1} vs {name2}',
                f'Up in {name2}': up,
                f'Down in {name2}': -down
            })

    if diff_summary:
        summary_df = pd.DataFrame(diff_summary)
        x = np.arange(len(summary_df))
        width = 0.35

        # Find the column names dynamically
        up_cols = [c for c in summary_df.columns if c.startswith('Up')]
        down_cols = [c for c in summary_df.columns if c.startswith('Down')]

        if up_cols and down_cols:
            for i, (up_col, down_col) in enumerate(zip(up_cols, down_cols)):
                ax3.barh(i, summary_df.iloc[i][up_col], color='#e74c3c', alpha=0.8, label='Up in Mut' if i == 0 else '')
                ax3.barh(i, summary_df.iloc[i][down_col], color='#2ecc71', alpha=0.8, label='Down in Mut' if i == 0 else '')

            ax3.set_yticks(range(len(summary_df)))
            ax3.set_yticklabels(summary_df['Comparison'])
            ax3.axvline(x=0, color='black', linewidth=0.5)
            ax3.set_xlabel('Number of Differential Motifs')
            ax3.legend(loc='lower right')

    ax3.set_title('C. Differential TF Motif Enrichment', fontweight='bold')

    # Panel D: Correlation between conditions
    ax4 = fig.add_subplot(2, 2, 4)
    if len(enr_cols) >= 2:
        corr_matrix = matrix_df[enr_cols].corr()
        corr_matrix.index = [c.replace('_LogP', '') for c in corr_matrix.index]
        corr_matrix.columns = [c.replace('_LogP', '') for c in corr_matrix.columns]
        sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm',
                    center=0, ax=ax4, vmin=-1, vmax=1,
                    cbar_kws={'shrink': 0.6, 'label': 'Correlation'})
    ax4.set_title('D. Enrichment Correlation Between Conditions', fontweight='bold')

    plt.tight_layout()

    # Save
    output_file = output_dir / 'summary_figure_tfbs_analysis.png'
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    logger.info(f"  Saved: {output_file}")

    # Also save as PDF
    fig = plt.figure(figsize=(16, 12))
    # Recreate the same figure for PDF
    ax1 = fig.add_subplot(2, 2, 1)
    sns.heatmap(heatmap_data, cmap='YlOrRd', annot=False, ax=ax1,
                cbar_kws={'shrink': 0.6, 'label': '-log10(p)'})
    ax1.set_title('A. Top 25 Enriched TF Motifs', fontweight='bold')

    ax2 = fig.add_subplot(2, 2, 2)
    family_counts.plot(kind='barh', ax=ax2, color='steelblue')
    ax2.set_xlabel('Number of Motifs')
    ax2.set_title('B. TF Family Distribution', fontweight='bold')
    ax2.invert_yaxis()

    ax3 = fig.add_subplot(2, 2, 3)
    ax3.set_title('C. Differential TF Motif Enrichment', fontweight='bold')

    ax4 = fig.add_subplot(2, 2, 4)
    if len(enr_cols) >= 2:
        sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm',
                    center=0, ax=ax4, vmin=-1, vmax=1)
    ax4.set_title('D. Enrichment Correlation Between Conditions', fontweight='bold')

    plt.tight_layout()
    pdf_file = output_dir / 'summary_figure_tfbs_analysis.pdf'
    plt.savefig(pdf_file, bbox_inches='tight')
    plt.close()


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Visualize TFBS motif enrichment results'
    )
    parser.add_argument(
        '--input-dir', '-i',
        type=Path,
        default=INPUT_DIR,
        help=f'Input directory with differential analysis results (default: {INPUT_DIR})'
    )
    parser.add_argument(
        '--output-dir', '-o',
        type=Path,
        default=OUTPUT_DIR,
        help=f'Output directory for figures (default: {OUTPUT_DIR})'
    )
    parser.add_argument(
        '--top-n',
        type=int,
        default=50,
        help='Number of top motifs to show in heatmap (default: 50)'
    )
    args = parser.parse_args()

    logger.info("=" * 60)
    logger.info("TFBS Motif Enrichment Visualization")
    logger.info("=" * 60)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load enrichment matrix
    matrix_file = args.input_dir / "tf_enrichment_matrix_all_samples.tsv"
    if not matrix_file.exists():
        logger.error(f"Enrichment matrix not found: {matrix_file}")
        sys.exit(1)

    matrix_df = pd.read_csv(matrix_file, sep='\t')
    logger.info(f"Loaded enrichment matrix: {len(matrix_df)} motifs")

    # Load differential analysis results
    diff_results = {}
    for diff_file in args.input_dir.glob("differential_motifs_*.tsv"):
        if 'significant' not in diff_file.name:
            name = diff_file.stem.replace('differential_motifs_', '')
            diff_results[name] = pd.read_csv(diff_file, sep='\t')
            logger.info(f"Loaded differential results: {name}")

    # Generate visualizations
    plot_enrichment_heatmap(matrix_df, args.output_dir, args.top_n)
    plot_tf_family_enrichment(matrix_df, args.output_dir)
    plot_condition_comparison(matrix_df, args.output_dir)

    # Volcano plots for each comparison
    for name, df in diff_results.items():
        plot_volcano(df, name, args.output_dir)

    # Summary figure
    plot_summary_figure(matrix_df, diff_results, args.output_dir)

    logger.info("=" * 60)
    logger.info("Visualization Complete")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
