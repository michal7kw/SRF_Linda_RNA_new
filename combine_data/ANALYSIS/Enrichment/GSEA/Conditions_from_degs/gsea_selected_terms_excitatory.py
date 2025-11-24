#!/usr/bin/env python3
"""
Selected Terms GSEA Heatmap Generator - Excitatory L1 Neurons Only

Creates a heatmap with only user-specified terms from a selected terms list.
Analyzes excitatory neurons (cell_type_L1 == 'Excitatory') combining all subtypes.

Features:
- Plots only terms from a selected terms file
- Numbered rows with separate legend table
- No significance asterisks on the heatmap
- Clean, publication-ready format
- Genotype-stratified analysis (Emx1 vs Nestin)
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
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
GSEA_BASE_DIR = os.path.join(BASE_RESULTS_DIR, 'gsea_analysis_excitatory_L1')
OUTPUT_DIR = os.path.join(BASE_RESULTS_DIR, 'gsea_selected_terms_excitatory_L1')
SELECTED_TERMS_FILE = os.path.join(WORKING_DIR, "ANALYSIS/Enrichment/GSEA/selected_terms.txt")

pathlib.Path(OUTPUT_DIR).mkdir(exist_ok=True)

# Publication settings
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


def load_selected_terms(filepath):
    """
    Load selected terms from the text file.

    Returns:
    --------
    list : List of selected term names
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Selected terms file not found: {filepath}")

    with open(filepath, 'r') as f:
        terms = [line.strip() for line in f if line.strip()]

    print(f"Loaded {len(terms)} selected terms from {filepath}")
    return terms


def normalize_term_name(term):
    """
    Normalize term names for matching between different formats.
    Handles both GO/Reactome format variations.
    """
    # Remove common prefixes
    prefixes = ['GOBP_', 'GOCC_', 'GOMF_', 'GO_', 'REACTOME_', 'R-HSA-', 'KEGG_', 'HALLMARK_']
    for prefix in prefixes:
        if term.startswith(prefix):
            term = term[len(prefix):]
            break

    # Convert to uppercase and replace separators with spaces
    term = term.upper()
    term = term.replace('_', ' ')
    term = term.replace('-', ' ')

    # Remove extra spaces
    term = ' '.join(term.split())

    return term


def shorten_term_name(term, max_length=60):
    """
    Intelligently shorten pathway/GO term names for the legend table.
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
        'Mitochondrial': 'Mito.',
    }

    for full, abbrev in sorted(abbreviations.items(), key=lambda x: len(x[0]), reverse=True):
        term = term.replace(full, abbrev)

    # Remove extra spaces
    term = ' '.join(term.split())

    # Truncate if still too long
    if len(term) > max_length:
        term = term[:max_length-3] + '...'

    return term


def load_gsea_results(genotype_dir):
    """Load GSEA results for excitatory L1 neurons from a given genotype."""
    comparison_name = "Excitatory_L1_Mutant_vs_Control"

    # Try multiple possible file patterns
    possible_files = [
        pathlib.Path(GSEA_BASE_DIR) / genotype_dir / comparison_name / "GSEA_Full_Report_All_Genesets.csv",
        pathlib.Path(GSEA_BASE_DIR) / genotype_dir / comparison_name / "GSEA_Full_Report_All_Genesets.xlsx",
        pathlib.Path(GSEA_BASE_DIR) / genotype_dir / comparison_name / "GSEA_Summary_All_Genesets.csv",
        pathlib.Path(GSEA_BASE_DIR) / genotype_dir / comparison_name / "GSEA_Summary_All_Genesets.xlsx",
        pathlib.Path(GSEA_BASE_DIR) / genotype_dir / comparison_name / "GSEA_results.csv",
        pathlib.Path(GSEA_BASE_DIR) / genotype_dir / comparison_name / "GSEA_results.xlsx",
    ]

    file_path = None
    for fp in possible_files:
        if fp.exists():
            file_path = fp
            break

    if file_path is None:
        print(f"Warning: No GSEA results file found for {genotype_dir}")
        print(f"Tried locations:")
        for fp in possible_files:
            print(f"  - {fp}")
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
        print(f"Available columns: {df.columns.tolist()}")
        return None

    df = df[['Term', 'nes', 'fdr']].copy()

    # Handle duplicate terms by keeping the one with lowest FDR
    if df['Term'].duplicated().any():
        print(f"  - Found duplicate terms. Resolving by lowest FDR.")
        df = df.sort_values('fdr')
        df = df[~df['Term'].duplicated(keep='first')]

    return df


def match_selected_terms(gsea_df, selected_terms):
    """
    Match selected terms to GSEA results using fuzzy matching.

    Returns:
    --------
    DataFrame with matched terms
    """
    matched_terms = []

    # Normalize selected terms for matching
    normalized_selected = {normalize_term_name(t): t for t in selected_terms}

    # Normalize GSEA terms
    gsea_df['normalized_term'] = gsea_df['Term'].apply(normalize_term_name)

    # Match terms
    for norm_selected, orig_selected in normalized_selected.items():
        matching_rows = gsea_df[gsea_df['normalized_term'] == norm_selected]

        if len(matching_rows) > 0:
            matched_terms.append(matching_rows.iloc[0]['Term'])
            print(f"  ✓ Matched: {orig_selected}")
        else:
            print(f"  ✗ Not found: {orig_selected}")

    # Filter to matched terms
    result_df = gsea_df[gsea_df['Term'].isin(matched_terms)].copy()
    result_df = result_df.drop(columns=['normalized_term'])

    return result_df


def create_term_legend_table(ax, term_dict, fontsize=6):
    """
    Create a clean table showing term numbers and their names.
    Single column format with "1. Pathway Name" style.
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


def create_selected_terms_heatmap(selected_terms_file=SELECTED_TERMS_FILE, output_suffix=''):
    """
    Create a heatmap with only selected terms (no FDR significance indicators).

    Parameters:
    -----------
    selected_terms_file : str
        Path to the selected terms text file
    output_suffix : str
        Optional suffix for output files
    """
    print(f"\n{'='*70}")
    print(f"Creating heatmap for selected terms in Excitatory L1 neurons")
    print(f"{'='*70}")

    # Load selected terms
    selected_terms = load_selected_terms(selected_terms_file)

    # Load data for both genotypes
    all_results = []

    for genotype_name, genotype_dir in GENOTYPES.items():
        df = load_gsea_results(genotype_dir)
        if df is not None:
            # Match selected terms
            print(f"\nMatching selected terms in {genotype_name} data:")
            df_matched = match_selected_terms(df, selected_terms)

            if len(df_matched) == 0:
                print(f"Warning: No matching terms found for {genotype_name}")
                continue

            df_renamed = df_matched.copy()
            df_renamed.rename(columns={
                'nes': f"{genotype_name}_nes",
                'fdr': f"{genotype_name}_fdr"
            }, inplace=True)

            # Handle any remaining duplicates before setting index
            if df_renamed['Term'].duplicated().any():
                print(f"  - Removing duplicates in {genotype_name} matched results")
                df_renamed = df_renamed.sort_values(f"{genotype_name}_fdr")
                df_renamed = df_renamed[~df_renamed['Term'].duplicated(keep='first')]

            df_renamed.set_index('Term', inplace=True)
            all_results.append(df_renamed)

    if len(all_results) != len(GENOTYPES):
        print(f"\nError: Could not load results for all genotypes.")
        return

    # Merge results
    merged_df = pd.concat(all_results, axis=1, join='outer')

    if merged_df.empty:
        print(f"No data found. Exiting.")
        return

    print(f"\nFound {len(merged_df)} matched terms across both genotypes.")

    # Get NES columns
    nes_cols = [f"{name}_nes" for name in GENOTYPES.keys()]

    # Sort by mean absolute NES to order terms
    merged_df['mean_abs_nes'] = merged_df[nes_cols].abs().mean(axis=1)
    merged_df = merged_df.sort_values('mean_abs_nes', ascending=False)

    # Create term number mapping
    term_dict = {}
    numbered_terms = []
    for i, term in enumerate(merged_df.index, 1):
        term_dict[i] = shorten_term_name(term)
        numbered_terms.append(str(i))

    # Prepare plot data with numbers
    plot_df = merged_df[nes_cols].copy()
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

    # Create heatmap WITHOUT annotations (no FDR stars)
    sns.heatmap(
        plot_df,
        cmap='RdBu_r',
        vmin=vmin,
        vmax=vmax,
        center=0,
        annot=False,  # No annotations
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
    ax_heatmap.set_xlabel('Genotype', fontsize=7, fontweight='bold')
    ax_heatmap.set_ylabel('Pathway #', fontsize=7, fontweight='bold')
    ax_heatmap.set_title(f'Selected Terms in Excitatory L1 Neurons\n(Mutant vs Control)',
                         fontsize=8, fontweight='bold', pad=8)

    ax_heatmap.tick_params(axis='x', rotation=45, labelsize=6.5)
    ax_heatmap.tick_params(axis='y', rotation=0, labelsize=6.5)

    # Create term legend table
    create_term_legend_table(ax_legend, term_dict, fontsize=6)

    # Adjust layout
    plt.tight_layout()

    # Save outputs
    suffix_str = f"_{output_suffix}" if output_suffix else ""
    base_filename = f"GSEA_Selected_Terms_Excitatory_L1{suffix_str}"

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

    # Save the actual NES and FDR values
    data_export = merged_df.copy()
    data_export.index.name = 'Term'
    data_csv_path = os.path.join(OUTPUT_DIR, f"{base_filename}_data.csv")
    data_export.to_csv(data_csv_path)
    print(f"✓ Saved data CSV: {data_csv_path}")

    plt.close()


def main():
    """Main function to generate selected terms GSEA heatmap."""
    parser = argparse.ArgumentParser(
        description="Generate GSEA heatmap with only selected terms for Excitatory L1 neurons (no FDR indicators)."
    )
    parser.add_argument(
        '--terms_file',
        type=str,
        default=SELECTED_TERMS_FILE,
        help=f"Path to selected terms file (default: {SELECTED_TERMS_FILE})"
    )
    parser.add_argument(
        '--output_suffix',
        type=str,
        default='',
        help="Optional suffix for output files"
    )

    args = parser.parse_args()

    print("=" * 70)
    print("Selected Terms GSEA Heatmap Generator - Excitatory L1 Neurons")
    print("=" * 70)
    print(f"GSEA base directory: {GSEA_BASE_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Selected terms file: {args.terms_file}")
    print(f"Format: Numbers in heatmap, names in separate legend (NO FDR indicators)")

    create_selected_terms_heatmap(
        selected_terms_file=args.terms_file,
        output_suffix=args.output_suffix
    )

    print("\n" + "=" * 70)
    print("✓ Selected terms heatmap generation completed!")
    print(f"✓ Output directory: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == '__main__':
    main()
