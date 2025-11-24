#!/usr/bin/env python3
"""
Publication-Ready GSEA Heatmap Generator

Creates compact, publication-quality heatmaps of GSEA results suitable for
inclusion in scientific paper figures. Focuses on clarity, readability, and
space efficiency while maintaining scientific rigor.

Features:
- Compact layout optimized for multi-panel figures
- Smart term selection and abbreviation
- Publication-standard fonts and formatting
- Multiple export formats (PDF vector + high-res PNG)
- Customizable term limits and filtering
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import pathlib
import numpy as np
import argparse
from matplotlib import rcParams

# --- Configuration ---
PROJECT_DIR = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)

DEFAULT_GSEA_scRNA_PROCESSING = True

BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")

if DEFAULT_GSEA_scRNA_PROCESSING == True:
    GSEA_BASE_DIR = os.path.join(BASE_RESULTS_DIR, 'gsea_analysis_between_conditions_default')
    OUTPUT_DIR = os.path.join(BASE_RESULTS_DIR, 'gsea_summary_visuals_default')
else:
    GSEA_BASE_DIR = os.path.join(BASE_RESULTS_DIR, 'gsea_analysis_between_conditions')
    OUTPUT_DIR = os.path.join(BASE_RESULTS_DIR, 'gsea_summary_visuals')

pathlib.Path(OUTPUT_DIR).mkdir(exist_ok=True)

# Publication settings
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
rcParams['font.size'] = 7
rcParams['axes.labelsize'] = 8
rcParams['axes.titlesize'] = 9
rcParams['xtick.labelsize'] = 7
rcParams['ytick.labelsize'] = 7
rcParams['legend.fontsize'] = 7
rcParams['figure.titlesize'] = 9
rcParams['pdf.fonttype'] = 42  # TrueType fonts for editability
rcParams['ps.fonttype'] = 42
rcParams['svg.fonttype'] = 'none'

# Define the analyses to include in the heatmap
ANALYSES = {
    'Emx1': 'emx1',
    'Nestin': 'nestin',
    'Both': 'both_genotypes'
}

def shorten_term_name(term, max_length=50):
    """
    Intelligently shorten pathway/GO term names for better readability.

    Parameters:
    -----------
    term : str
        Original term name
    max_length : int
        Maximum length for the shortened term

    Returns:
    --------
    str
        Shortened term name
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

    # Common abbreviations for biological terms
    abbreviations = {
        'Regulation Of': 'Reg.',
        'Response To': 'Resp.',
        'Positive Regulation': '+Reg.',
        'Negative Regulation': '-Reg.',
        'Biological Process': 'Process',
        'Signaling Pathway': 'Signaling',
        'Pathway': 'Path.',
        'Activity': 'Act.',
        'Process': 'Proc.',
        'Development': 'Dev.',
        'Differentiation': 'Diff.',
        'Metabolism': 'Metab.',
        'Biosynthesis': 'Biosynth.',
        'Organization': 'Org.',
        'Assembly': 'Asm.',
        'Transport': 'Transp.',
        'Localization': 'Local.',
        'Proliferation': 'Prolif.',
        'Apoptotic': 'Apopt.',
        'Programmed Cell Death': 'PCD',
    }

    for full, abbrev in abbreviations.items():
        term = term.replace(full, abbrev)

    # Truncate if still too long
    if len(term) > max_length:
        term = term[:max_length-3] + '...'

    return term


def load_gsea_results(analysis_name, genotype_dir, cluster_type, input_file):
    """Loads filtered GSEA results for a given analysis."""
    sanitized_cluster_name = cluster_type.replace(' ', '_')
    comparison_name = "Mutant_vs_Control"

    file_path = pathlib.Path(GSEA_BASE_DIR) / genotype_dir / f"{sanitized_cluster_name}_{comparison_name}" / input_file

    # Try .xlsx first, then .csv if not found
    if not file_path.exists():
        # Try CSV version if Excel not found
        csv_path = file_path.with_suffix('.csv')
        if csv_path.exists():
            file_path = csv_path
        else:
            print(f"Warning: File not found for {analysis_name}: {file_path}")
            print(f"  Also tried: {csv_path}")
            return None

    print(f"Loading results for {analysis_name} from {file_path}")

    # Load based on file extension
    try:
        if file_path.suffix == '.csv':
            df = pd.read_csv(file_path)
        else:
            try:
                df = pd.read_excel(file_path)
            except ImportError:
                # If openpyxl not available, try CSV version
                csv_path = file_path.with_suffix('.csv')
                if csv_path.exists():
                    print(f"  openpyxl not available, using CSV version: {csv_path}")
                    df = pd.read_csv(csv_path)
                else:
                    raise ImportError("openpyxl not installed and no CSV version found. "
                                    "Please install openpyxl or provide CSV files.")
    except Exception as e:
        print(f"Error loading file {file_path}: {e}")
        return None

    if 'Term' not in df.columns or 'nes' not in df.columns or 'fdr' not in df.columns:
        print(f"Warning: 'Term', 'nes', or 'fdr' column not found in {file_path}")
        return None

    # Handle potential duplicate terms by keeping the one with the max absolute NES
    if df['Term'].duplicated().any():
        print(f"  - Found duplicate terms in {analysis_name}. Resolving by max absolute NES.")
        df['abs_nes'] = df['nes'].abs()
        df = df.loc[df.groupby('Term')['abs_nes'].idxmax()]
        df = df.drop(columns=['abs_nes'])

    # Keep only Term, nes, and fdr, and rename columns
    df = df[['Term', 'nes', 'fdr']].copy()
    df.rename(columns={
        'nes': f"{analysis_name}_nes",
        'fdr': f"{analysis_name}_fdr"
    }, inplace=True)
    df.set_index('Term', inplace=True)

    return df


def create_publication_heatmap(plot_df, merged_df, cluster_type, max_terms=20,
                               fdr_threshold=0.25, output_suffix='', no_fdr=False):
    """
    Create a compact, publication-ready heatmap.

    Parameters:
    -----------
    plot_df : DataFrame
        NES values to plot
    merged_df : DataFrame
        Full merged data including FDR values
    cluster_type : str
        Cluster type name
    max_terms : int
        Maximum number of terms to display
    fdr_threshold : float
        FDR threshold for significance stars
    output_suffix : str
        Optional suffix for output filename
    no_fdr : bool
        If True, create heatmap without FDR significance stars
    """

    # Limit number of terms
    if len(plot_df) > max_terms:
        print(f"  Limiting to top {max_terms} terms (from {len(plot_df)})")
        plot_df = plot_df.head(max_terms)

    # Shorten term names
    plot_df.index = [shorten_term_name(term) for term in plot_df.index]

    # Calculate figure dimensions based on number of terms
    n_terms = len(plot_df)
    n_cols = len(plot_df.columns)

    # Compact sizing for publication
    fig_width = min(2.5 + (n_cols * 0.6), 7)  # Max 7 inches (single column)
    fig_height = min(1.5 + (n_terms * 0.20), 10)  # 0.2 inches per term

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Color scale centered at 0
    vmax = np.nanmax(np.abs(plot_df.values))
    vmax = max(vmax, 0.1)
    vmin = -vmax

    # Create annotation with significance stars (or empty if no_fdr)
    annot_df = pd.DataFrame(index=plot_df.index, columns=plot_df.columns, dtype=str)

    original_terms = merged_df.index[:len(plot_df)]  # Get original term names

    if not no_fdr:
        for idx, (short_term, orig_term) in enumerate(zip(plot_df.index, original_terms)):
            for analysis in ANALYSES.keys():
                fdr_col = f"{analysis}_fdr"
                nes_col = f"{analysis}_nes"

                if orig_term in merged_df.index:
                    fdr_val = merged_df.loc[orig_term, fdr_col]
                    nes_val = merged_df.loc[orig_term, nes_col]

                    if pd.notna(fdr_val) and pd.notna(nes_val):
                        # Add significance stars
                        if fdr_val < 0.001:
                            sig = '***'
                        elif fdr_val < 0.01:
                            sig = '**'
                        elif fdr_val < fdr_threshold:
                            sig = '*'
                        else:
                            sig = ''

                        annot_df.loc[short_term, analysis] = sig
                    else:
                        annot_df.loc[short_term, analysis] = ''
                else:
                    annot_df.loc[short_term, analysis] = ''
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
        linewidths=0.5,
        linecolor='gray',
        cbar_kws={
            'label': 'NES',
            'shrink': 0.8,
            'aspect': 20,
            'pad': 0.02
        },
        ax=ax
    )

    # Formatting
    ax.set_xlabel('', fontweight='bold')
    ax.set_ylabel('', fontweight='bold')
    ax.set_title(f'{cluster_type.replace("_", " ")}',
                fontweight='bold', pad=10)

    # Rotate labels
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    # Add legend for significance (only if not no_fdr)
    if not no_fdr:
        legend_elements = [
            mpatches.Patch(facecolor='white', edgecolor='black', label='* FDR < 0.25'),
            mpatches.Patch(facecolor='white', edgecolor='black', label='** FDR < 0.01'),
            mpatches.Patch(facecolor='white', edgecolor='black', label='*** FDR < 0.001')
        ]

        # Position legend outside plot area
        ax.legend(handles=legend_elements,
                 bbox_to_anchor=(1.02, 1),
                 loc='upper left',
                 frameon=True,
                 fontsize=6)

    plt.tight_layout()

    # Save in multiple formats
    fdr_suffix = '_no_fdr' if no_fdr else ''
    base_filename = f"GSEA_publication_{cluster_type}{output_suffix}{fdr_suffix}"

    # High-resolution PNG
    png_path = os.path.join(OUTPUT_DIR, f"{base_filename}.png")
    plt.savefig(png_path, dpi=600, bbox_inches='tight', facecolor='white')
    print(f"  Saved PNG: {png_path}")

    # Vector PDF for editing (commented out for now)
    # pdf_path = os.path.join(OUTPUT_DIR, f"{base_filename}.pdf")
    # plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    # print(f"  Saved PDF: {pdf_path}")

    # SVG for maximum editability (commented out for now)
    # svg_path = os.path.join(OUTPUT_DIR, f"{base_filename}.svg")
    # plt.savefig(svg_path, bbox_inches='tight', facecolor='white')
    # print(f"  Saved SVG: {svg_path}")

    plt.close()


def create_compact_versions(merged_df, cluster_type, nes_cols, no_fdr=False):
    """
    Create multiple compact versions with different term selections.

    Parameters:
    -----------
    merged_df : DataFrame
        Merged GSEA results
    cluster_type : str
        Cluster type name
    nes_cols : list
        List of NES column names
    no_fdr : bool
        If True, create heatmaps without FDR significance stars
    """

    print(f"\n--- Creating Publication-Ready Versions{' (no FDR)' if no_fdr else ''} ---")

    # Version 1: Top 10 most significant (by mean absolute NES)
    print("\n1. Top 10 by mean absolute NES:")
    merged_df['mean_abs_nes'] = merged_df[nes_cols].abs().mean(axis=1)
    top_10 = merged_df.nlargest(10, 'mean_abs_nes')
    plot_df = top_10[nes_cols].copy()
    plot_df.columns = ANALYSES.keys()
    create_publication_heatmap(plot_df, top_10, cluster_type,
                               max_terms=10, output_suffix='_top10', no_fdr=no_fdr)

    # Version 2: Top 5 upregulated + Top 5 downregulated
    print("\n2. Top 5 up + Top 5 down:")
    merged_df['mean_nes'] = merged_df[nes_cols].mean(axis=1)
    top_up = merged_df[merged_df['mean_nes'] > 0].nlargest(5, 'mean_nes')
    top_down = merged_df[merged_df['mean_nes'] < 0].nsmallest(5, 'mean_nes')
    balanced = pd.concat([top_up, top_down]).sort_values('mean_nes', ascending=False)
    plot_df = balanced[nes_cols].copy()
    plot_df.columns = ANALYSES.keys()
    create_publication_heatmap(plot_df, balanced, cluster_type,
                               max_terms=10, output_suffix='_balanced', no_fdr=no_fdr)

    # Version 3: Top 15 (standard)
    print("\n3. Top 15 terms:")
    top_15 = merged_df.nlargest(15, 'mean_abs_nes')
    plot_df = top_15[nes_cols].copy()
    plot_df.columns = ANALYSES.keys()
    create_publication_heatmap(plot_df, top_15, cluster_type,
                               max_terms=15, output_suffix='_top15', no_fdr=no_fdr)

    # Version 4: All significant terms with consistent enrichment
    print("\n4. Consistently enriched terms (all analyses FDR < 0.25):")
    fdr_cols = [f"{name}_fdr" for name in ANALYSES.keys()]
    consistent = merged_df[merged_df[fdr_cols].max(axis=1) < 0.25]

    if len(consistent) > 0:
        # Limit to top 20 by mean abs NES
        if len(consistent) > 20:
            consistent = consistent.nlargest(20, 'mean_abs_nes')

        plot_df = consistent[nes_cols].copy()
        plot_df.columns = ANALYSES.keys()
        create_publication_heatmap(plot_df, consistent, cluster_type,
                                   max_terms=20, output_suffix='_consistent', no_fdr=no_fdr)
        print(f"  Found {len(consistent)} consistently enriched terms")
    else:
        print("  No consistently enriched terms found")


def main():
    """
    Generates publication-ready GSEA heatmaps in multiple compact formats.
    """
    parser = argparse.ArgumentParser(
        description="Generate compact, publication-ready GSEA heatmaps."
    )
    parser.add_argument(
        '--cluster_type',
        type=str,
        required=True,
        choices=['Mature_GC', 'Combined_GC', 'GABA'],
        help="The cluster type to generate the heatmap for"
    )
    parser.add_argument(
        '--input_file',
        type=str,
        default=None,
        help="Input GSEA summary file name (.xlsx or .csv). If not specified, auto-detects based on cluster type."
    )
    parser.add_argument(
        '--max_terms',
        type=int,
        default=15,
        help="Maximum number of terms to display in main plot (default: 15)"
    )
    parser.add_argument(
        '--fdr_threshold',
        type=float,
        default=0.25,
        help="FDR threshold for significance (default: 0.25)"
    )
    parser.add_argument(
        '--create_variants',
        action='store_true',
        help="Create multiple variant heatmaps with different term selections"
    )

    args = parser.parse_args()
    cluster_type = args.cluster_type
    input_file = args.input_file

    # Auto-detect input file based on cluster type if not specified
    if input_file is None:
        if cluster_type == 'GABA':
            input_file = 'GSEA_Focus_RNA_Splicing.csv'
            print(f"Auto-detected input file for GABA: {input_file}")
        else:  # Mature_GC or Combined_GC
            input_file = 'GSEA_Focus_Neuro_CellDeath.csv'
            print(f"Auto-detected input file for {cluster_type}: {input_file}")

    print("=" * 70)
    print(f"Publication-Ready GSEA Heatmap Generator")
    print(f"Cluster type: {cluster_type}")
    print(f"Input file: {input_file}")
    print("=" * 70)

    # Load data from all analyses
    all_results = []
    for analysis_name, genotype_dir in ANALYSES.items():
        results_df = load_gsea_results(analysis_name, genotype_dir, cluster_type, input_file)
        if results_df is not None:
            all_results.append(results_df)

    if not all_results:
        print(f"No GSEA result files found for cluster type '{cluster_type}'. Exiting.")
        return

    # Merge all dataframes on the 'Term' index
    merged_df = pd.concat(all_results, axis=1, join='outer')

    if merged_df.empty:
        print("No GSEA result files found or they were empty. Exiting.")
        return

    print(f"\nFound {len(merged_df)} total terms across all analyses.")

    # Identify NES columns
    nes_cols = [f"{name}_nes" for name in ANALYSES.keys()]

    if args.create_variants:
        # Create multiple compact versions (with FDR)
        create_compact_versions(merged_df, cluster_type, nes_cols, no_fdr=False)
        # Create multiple compact versions (without FDR)
        create_compact_versions(merged_df, cluster_type, nes_cols, no_fdr=True)
    else:
        # Create single standard version (with FDR)
        print("\n--- Creating Standard Publication Heatmap ---")

        # Calculate mean absolute NES for ranking
        merged_df['mean_abs_nes'] = merged_df[nes_cols].abs().mean(axis=1)

        # Select top terms
        top_terms = merged_df.nlargest(args.max_terms, 'mean_abs_nes')

        print(f"Selected top {len(top_terms)} terms by mean absolute NES")

        # Prepare plot data
        plot_df = top_terms[nes_cols].copy()
        plot_df.columns = ANALYSES.keys()

        # Create heatmap with FDR
        create_publication_heatmap(
            plot_df,
            top_terms,
            cluster_type,
            max_terms=args.max_terms,
            fdr_threshold=args.fdr_threshold,
            no_fdr=False
        )

        # Create heatmap without FDR
        print("\n--- Creating Standard Publication Heatmap (no FDR) ---")
        create_publication_heatmap(
            plot_df,
            top_terms,
            cluster_type,
            max_terms=args.max_terms,
            fdr_threshold=args.fdr_threshold,
            no_fdr=True
        )

    print("\n" + "=" * 70)
    print("Publication heatmap generation completed!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == '__main__':
    main()
