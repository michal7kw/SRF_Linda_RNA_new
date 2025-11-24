#!/usr/bin/env python3
"""
Focused GSEA Heatmap Generator - Numbered Version

Creates ultra-compact heatmaps with numbered terms instead of names.
The term names are displayed in a separate table/legend alongside the heatmap.

This design is optimized for:
- Maximum space efficiency in multi-panel figures
- Very limited figure space
- Professional publication quality

Features:
- Numbered rows (1, 2, 3...) instead of long term names
- Separate term legend table
- Side-by-side layout (heatmap + legend)
- Extremely compact dimensions
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import os
import pathlib
import numpy as np
import argparse
from matplotlib import rcParams

# --- Configuration ---
PROJECT_DIR = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)

BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")
GSEA_BASE_DIR = os.path.join(BASE_RESULTS_DIR, 'gsea_analysis_between_conditions')
OUTPUT_DIR = os.path.join(BASE_RESULTS_DIR, 'gsea_focused_heatmaps_numbered')

pathlib.Path(OUTPUT_DIR).mkdir(exist_ok=True)

# Publication settings - ultra-compact
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
rcParams['font.size'] = 6.5
rcParams['axes.labelsize'] = 7
rcParams['axes.titlesize'] = 8
rcParams['xtick.labelsize'] = 6.5
rcParams['ytick.labelsize'] = 6.5
rcParams['legend.fontsize'] = 5.5
rcParams['figure.titlesize'] = 8
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

# Define genotypes to compare
GENOTYPES = {
    'Emx1': 'emx1',
    'Nestin': 'nestin'
}

# Define focus configurations
FOCUS_CONFIG = {
    'cell_death': {
        'keywords': [
            'neurodegeneration', 'neuronal death', 'axonopathy', 'neuron death',
            'apoptosis', 'apoptotic', 'necroptosis', 'pyroptosis', 'ferroptosis',
            'autophagy', 'autophagic', 'parthanatos', 'cell death'
        ],
        'file_suffix': 'Neuro_CellDeath',
        'title_name': 'Cell Death',
        'clusters': ['Mature_GC', 'Combined_GC']
    },
    'splicing': {
        'keywords': [
            'splicing', 'spliceosome', 'splice', 'rna processing',
            'alternative splicing', 'splice variant', 'isoform',
            'pre-mrna', 'intron', 'exon',
            'splice site', 'exon skipping', 'intron retention',
            'sr protein', 'hnrnp', 'splicing factor',
            'aberrant splicing', 'spliceopathy',
            'nonsense-mediated decay', 'nmd'
        ],
        'file_suffix': 'RNA_Splicing',
        'title_name': 'RNA Splicing',
        'clusters': ['GABA']
    }
}


def shorten_term_name(term, max_length=60):
    """
    Intelligently shorten pathway/GO term names for the legend table.

    More conservative than the compact version since these will be in a table.
    """
    # Remove common prefixes
    prefixes_to_remove = [
        'GOBP_', 'GOCC_', 'GOMF_', 'KEGG_', 'REACTOME_', 'HALLMARK_',
        'GO_', 'WP_', 'BIOCARTA_', 'PID_'
    ]

    for prefix in prefixes_to_remove:
        if term.startswith(prefix):
            term = term[len(prefix):]
            break

    # Replace underscores with spaces and title case
    term = term.replace('_', ' ').title()

    # Moderate abbreviations for table readability
    abbreviations = {
        'Regulation Of': 'Reg.',
        'Positive Regulation': '+Reg.',
        'Negative Regulation': '-Reg.',
        'Programmed Cell Death': 'PCD',
        'Processing': 'Proc.',
        'Independent': 'Indep.',
        'Dependent': 'Dep.',
        'Mediated': 'Med.',
        'Activation': 'Activ.',
        'Recognition': 'Recog.',
    }

    for full, abbrev in sorted(abbreviations.items(), key=lambda x: len(x[0]), reverse=True):
        term = term.replace(full, abbrev)

    # Remove extra spaces
    term = ' '.join(term.split())

    # Truncate if still too long
    if len(term) > max_length:
        term = term[:max_length-3] + '...'

    return term


def load_focused_gsea_results(genotype_dir, cluster_name, focus_key):
    """Load focused GSEA results for a given genotype and cluster."""
    focus_config = FOCUS_CONFIG[focus_key]
    file_suffix = focus_config['file_suffix']

    comparison_name = "Mutant_vs_Control"
    file_path = pathlib.Path(GSEA_BASE_DIR) / genotype_dir / f"{cluster_name}_{comparison_name}" / f"GSEA_Focus_{file_suffix}.csv"

    # Try .csv first, then .xlsx if not found
    if not file_path.exists():
        file_path = file_path.with_suffix('.xlsx')
        if not file_path.exists():
            print(f"Warning: File not found for {genotype_dir}/{cluster_name}: {file_path}")
            return None

    print(f"Loading results from {file_path}")

    try:
        if file_path.suffix == '.csv':
            df = pd.read_csv(file_path)
        else:
            try:
                df = pd.read_excel(file_path)
            except ImportError:
                csv_path = file_path.with_suffix('.csv')
                if csv_path.exists():
                    print(f"  openpyxl not available, using CSV version: {csv_path}")
                    df = pd.read_csv(csv_path)
                else:
                    raise ImportError("openpyxl not installed and no CSV version found.")
    except Exception as e:
        print(f"Error loading file {file_path}: {e}")
        return None

    if 'Term' not in df.columns or 'nes' not in df.columns or 'fdr' not in df.columns:
        print(f"Warning: Required columns not found in {file_path}")
        return None

    df = df[['Term', 'nes', 'fdr']].copy()

    if df['Term'].duplicated().any():
        print(f"  - Found duplicate terms. Resolving by max absolute NES.")
        df['abs_nes'] = df['nes'].abs()
        df = df.loc[df.groupby('Term')['abs_nes'].idxmax()]
        df = df.drop(columns=['abs_nes'])

    return df


def create_term_legend_table(ax, term_dict, fontsize=6):
    """
    Create a clean table showing term numbers and their names.
    Single column format with "1. Pathway Name" style.

    Parameters:
    -----------
    ax : matplotlib axis
        Axis to draw the table on
    term_dict : dict
        Dictionary mapping term numbers to term names
    fontsize : int
        Font size for the table text
    """
    ax.axis('off')

    # Create table data with single column (number + term combined)
    table_data = []
    for num, term in term_dict.items():
        # Combine number and term in single cell
        table_data.append([f"{num}. {term}"])

    # Create table with single column
    table = ax.table(
        cellText=table_data,
        cellLoc='left',
        loc='center',
        bbox=[0, 0, 1, 1]  # Full axis space
    )

    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)

    # Calculate row height to fit all rows
    row_height = 1.0 / len(table_data)

    # Style data cells
    for i in range(len(table_data)):
        cell = table[(i, 0)]
        # Alternate row colors for readability
        if i % 2 == 0:
            cell.set_facecolor('#F8F8F8')
        else:
            cell.set_facecolor('white')
        cell.set_text_props(fontsize=fontsize)
        cell.set_edgecolor('lightgray')
        cell.set_linewidth(0.3)
        cell.set_height(row_height)

    return table


def create_numbered_heatmap(focus_key, max_terms=12, fdr_threshold=0.25, no_fdr=False):
    """
    Create a numbered heatmap with separate term legend.

    Parameters:
    -----------
    focus_key : str
        Focus area key ('cell_death' or 'splicing')
    max_terms : int
        Maximum number of terms to display
    fdr_threshold : float
        FDR threshold for significance
    no_fdr : bool
        If True, create heatmap without FDR significance stars
    """
    focus_config = FOCUS_CONFIG[focus_key]
    title_name = focus_config['title_name']
    clusters = focus_config['clusters']

    print(f"\n{'='*70}")
    print(f"Creating numbered heatmap for {title_name}{' (no FDR)' if no_fdr else ''}")
    print(f"{'='*70}")

    # Load data
    all_results = []
    cluster_found = None

    for cluster_name in clusters:
        print(f"\nTrying cluster: {cluster_name}")
        temp_results = []

        for genotype_name, genotype_dir in GENOTYPES.items():
            df = load_focused_gsea_results(genotype_dir, cluster_name, focus_key)
            if df is not None:
                df_renamed = df.copy()
                df_renamed.rename(columns={
                    'nes': f"{genotype_name}_nes",
                    'fdr': f"{genotype_name}_fdr"
                }, inplace=True)
                df_renamed.set_index('Term', inplace=True)
                temp_results.append(df_renamed)

        if len(temp_results) == len(GENOTYPES):
            all_results = temp_results
            cluster_found = cluster_name
            print(f"  ✓ Found results for cluster: {cluster_name}")
            break

    if not all_results:
        print(f"\nError: No GSEA result files found for {title_name}.")
        print(f"Tried clusters: {', '.join(clusters)}")
        return

    merged_df = pd.concat(all_results, axis=1, join='outer')

    if merged_df.empty:
        print(f"No data found for {title_name}. Exiting.")
        return

    print(f"\nFound {len(merged_df)} total terms across both genotypes.")

    # Filter and select terms
    nes_cols = [f"{name}_nes" for name in GENOTYPES.keys()]
    fdr_cols = [f"{name}_fdr" for name in GENOTYPES.keys()]

    merged_df['min_fdr'] = merged_df[fdr_cols].min(axis=1)
    merged_df['mean_abs_nes'] = merged_df[nes_cols].abs().mean(axis=1)

    significant_df = merged_df[merged_df['min_fdr'] < fdr_threshold].copy()

    if len(significant_df) == 0:
        print(f"  Warning: No terms pass FDR threshold {fdr_threshold}. Using all terms.")
        significant_df = merged_df.copy()

    top_terms = significant_df.nlargest(max_terms, 'mean_abs_nes')

    print(f"Selected top {len(top_terms)} terms by mean absolute NES")

    # Create term number mapping
    term_dict = {}
    numbered_terms = []
    for i, term in enumerate(top_terms.index, 1):
        term_dict[i] = shorten_term_name(term)
        numbered_terms.append(str(i))

    # Prepare plot data with numbers
    plot_df = top_terms[nes_cols].copy()
    plot_df.columns = GENOTYPES.keys()
    plot_df.index = numbered_terms

    # Calculate dimensions
    n_terms = len(plot_df)
    n_cols = len(plot_df.columns)

    # Ultra-compact heatmap (since we're using numbers)
    heatmap_width = min(1.8 + (n_cols * 0.5), 3.5)
    heatmap_height = min(1.5 + (n_terms * 0.2), 6)

    # Legend table width (wider to fit full term names)
    legend_width = 4.5

    # Total figure width
    total_width = heatmap_width + legend_width + 0.5  # 0.5 for spacing

    # Create figure with GridSpec for layout
    fig = plt.figure(figsize=(total_width, heatmap_height))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[heatmap_width, legend_width], wspace=0.2)

    ax_heatmap = fig.add_subplot(gs[0])
    ax_legend = fig.add_subplot(gs[1])

    # Color scale
    vmax = np.nanmax(np.abs(plot_df.values))
    vmax = max(vmax, 0.1)
    vmin = -vmax

    # Create significance annotations (or empty if no_fdr)
    annot_df = pd.DataFrame(index=plot_df.index, columns=plot_df.columns, dtype=str)

    if not no_fdr:
        for idx, (num, orig_term) in enumerate(zip(numbered_terms, top_terms.index)):
            for genotype_name in GENOTYPES.keys():
                fdr_col = f"{genotype_name}_fdr"

                if orig_term in top_terms.index:
                    fdr_val = top_terms.loc[orig_term, fdr_col]

                    if pd.notna(fdr_val):
                        if fdr_val < 0.001:
                            sig = '***'
                        elif fdr_val < 0.01:
                            sig = '**'
                        elif fdr_val < fdr_threshold:
                            sig = '*'
                        else:
                            sig = ''
                        annot_df.loc[num, genotype_name] = sig
                    else:
                        annot_df.loc[num, genotype_name] = ''
    else:
        # Fill with empty strings for no FDR version
        annot_df = annot_df.fillna('')

    # Create heatmap
    sns.heatmap(
        plot_df,
        cmap='RdBu_r',
        vmin=vmin,
        vmax=vmax,
        center=0,
        annot=annot_df,
        fmt='',
        linewidths=0.3,
        linecolor='lightgray',
        cbar_kws={
            'label': 'NES',
            'shrink': 0.7,
            'aspect': 12,
            'pad': 0.02
        },
        ax=ax_heatmap
    )

    # Format heatmap
    cluster_display = cluster_found.replace('_', ' ')
    ax_heatmap.set_xlabel('Genotype', fontsize=7, fontweight='bold')
    ax_heatmap.set_ylabel('Pathway #', fontsize=7, fontweight='bold')
    ax_heatmap.set_title(f'{title_name} in {cluster_display}\n(Mutant vs Control)',
                         fontsize=8, fontweight='bold', pad=8)

    ax_heatmap.tick_params(axis='x', rotation=45, labelsize=6.5)
    ax_heatmap.tick_params(axis='y', rotation=0, labelsize=6.5)

    # Add significance legend below x-axis label (only if not no_fdr)
    if not no_fdr:
        legend_elements = [
            mpatches.Patch(facecolor='white', edgecolor='black', label='* FDR<0.25'),
            mpatches.Patch(facecolor='white', edgecolor='black', label='** FDR<0.01'),
            mpatches.Patch(facecolor='white', edgecolor='black', label='*** FDR<0.001')
        ]

        ax_heatmap.legend(
            handles=legend_elements,
            bbox_to_anchor=(0.5, -0.20),  # Further down to avoid overlap
            loc='upper center',
            frameon=True,
            fontsize=5.5,
            ncol=3,
            handlelength=0.8,
            handletextpad=0.4,
            columnspacing=0.8
        )

    # Create term legend table
    create_term_legend_table(ax_legend, term_dict, fontsize=6)

    # Adjust layout with more bottom padding for FDR legend (if present)
    if not no_fdr:
        plt.tight_layout(rect=[0, 0.08, 1, 1])  # Leave 8% space at bottom for FDR legend
    else:
        plt.tight_layout()

    # Save outputs with more padding
    fdr_suffix = '_no_fdr' if no_fdr else ''
    base_filename = f"GSEA_Focused_{focus_key}_{cluster_found}_numbered{fdr_suffix}"

    png_path = os.path.join(OUTPUT_DIR, f"{base_filename}.png")
    plt.savefig(png_path, dpi=600, bbox_inches='tight', facecolor='white', pad_inches=0.15)
    print(f"\n✓ Saved PNG (600 DPI): {png_path}")

    pdf_path = os.path.join(OUTPUT_DIR, f"{base_filename}.pdf")
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white', pad_inches=0.15)
    print(f"✓ Saved PDF (vector): {pdf_path}")

    # Also save term mapping as CSV for reference
    term_df = pd.DataFrame(list(term_dict.items()), columns=['Number', 'Pathway'])
    csv_path = os.path.join(OUTPUT_DIR, f"{base_filename}_legend.csv")
    term_df.to_csv(csv_path, index=False)
    print(f"✓ Saved legend CSV: {csv_path}")

    plt.close()


def main():
    """Main function to generate numbered GSEA heatmaps."""
    parser = argparse.ArgumentParser(
        description="Generate ultra-compact GSEA heatmaps with numbered terms."
    )
    parser.add_argument(
        '--focus',
        type=str,
        choices=['cell_death', 'splicing', 'both'],
        default='both',
        help="Focus area: 'cell_death', 'splicing', or 'both' (default: both)"
    )
    parser.add_argument(
        '--max_terms',
        type=int,
        default=12,
        help="Maximum number of terms to display (default: 12)"
    )
    parser.add_argument(
        '--fdr_threshold',
        type=float,
        default=0.25,
        help="FDR threshold for significance (default: 0.25)"
    )

    args = parser.parse_args()

    print("=" * 70)
    print("Numbered GSEA Heatmap Generator (Ultra-Compact Version)")
    print("=" * 70)
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Max terms per heatmap: {args.max_terms}")
    print(f"FDR threshold: {args.fdr_threshold}")
    print(f"Format: Numbers in heatmap, names in separate legend")

    if args.focus == 'both' or args.focus == 'cell_death':
        # Create version with FDR
        create_numbered_heatmap(
            'cell_death',
            max_terms=args.max_terms,
            fdr_threshold=args.fdr_threshold,
            no_fdr=False
        )
        # Create version without FDR
        create_numbered_heatmap(
            'cell_death',
            max_terms=args.max_terms,
            fdr_threshold=args.fdr_threshold,
            no_fdr=True
        )

    if args.focus == 'both' or args.focus == 'splicing':
        # Create version with FDR
        create_numbered_heatmap(
            'splicing',
            max_terms=args.max_terms,
            fdr_threshold=args.fdr_threshold,
            no_fdr=False
        )
        # Create version without FDR
        create_numbered_heatmap(
            'splicing',
            max_terms=args.max_terms,
            fdr_threshold=args.fdr_threshold,
            no_fdr=True
        )

    print("\n" + "=" * 70)
    print("✓ Numbered heatmap generation completed!")
    print(f"✓ Output directory: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == '__main__':
    main()
