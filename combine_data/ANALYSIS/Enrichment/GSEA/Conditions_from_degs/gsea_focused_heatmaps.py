#!/usr/bin/env python3
"""
Focused GSEA Heatmap Generator

Creates small, focused heatmaps summarizing GSEA results for specific
biological categories across different model pairs (genotypes).

This script generates:
1. Cell death categories in GC cluster for both Emx1 and Nestin models
2. Splicing categories in IN cluster for both Emx1 and Nestin models

Each heatmap shows:
- Rows: Gene sets (pathways/terms) related to the focus area
- Columns: Different genotype comparisons (Emx1, Nestin)
- Colors: Normalized Enrichment Score (NES)
- Annotations: Significance stars based on FDR
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

BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")
GSEA_BASE_DIR = os.path.join(BASE_RESULTS_DIR, 'gsea_analysis_between_conditions')
OUTPUT_DIR = os.path.join(BASE_RESULTS_DIR, 'gsea_focused_heatmaps')

pathlib.Path(OUTPUT_DIR).mkdir(exist_ok=True)

# Publication settings - compact for limited space
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
rcParams['font.size'] = 7
rcParams['axes.labelsize'] = 8
rcParams['axes.titlesize'] = 9
rcParams['xtick.labelsize'] = 7
rcParams['ytick.labelsize'] = 7
rcParams['legend.fontsize'] = 6
rcParams['figure.titlesize'] = 9
rcParams['pdf.fonttype'] = 42  # TrueType fonts for editability
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
        'clusters': ['Mature_GC', 'Combined_GC']  # Try both cluster names
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
        'clusters': ['GABA']  # IN cluster
    }
}


def shorten_term_name(term, max_length=45):
    """
    Intelligently shorten pathway/GO term names for better readability in compact figures.

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

    # Common abbreviations for biological terms (more aggressive)
    abbreviations = {
        'Regulation Of': 'Reg.',
        'Positive Regulation Of': '+Reg.',
        'Negative Regulation Of': '-Reg.',
        'Positive Regulation': '+Reg.',
        'Negative Regulation': '-Reg.',
        'Response To': 'Resp.',
        'Biological Process': 'Proc.',
        'Signaling Pathway': 'Signal.',
        'Signaling': 'Signal.',
        'Pathway': 'Path.',
        'Activity': 'Act.',
        'Process': 'Proc.',
        'Development': 'Dev.',
        'Differentiation': 'Diff.',
        'Metabolism': 'Metab.',
        'Metabolic': 'Metab.',
        'Biosynthesis': 'Biosynth.',
        'Biosynthetic': 'Biosynth.',
        'Organization': 'Org.',
        'Assembly': 'Asm.',
        'Transport': 'Transp.',
        'Localization': 'Local.',
        'Proliferation': 'Prolif.',
        'Apoptotic': 'Apopt.',
        'Programmed Cell Death': 'PCD',
        'Independent': 'Indep.',
        'Dependent': 'Dep.',
        'Mediated': 'Med.',
        'Processing': 'Proc.',
        'Enhanced By': 'Enh. By',
        'Activation': 'Activ.',
        'Receptor': 'Rec.',
        'Recognition': 'Recog.',
        'And Cytosol': '& Cytosol',
        'In Nucleolus And': 'In Nucleolus &',
        'Via Spliceosome': 'Via Spliceos.',
        'Complex': 'Cmplx',
    }

    # Apply abbreviations in order (longer phrases first)
    for full, abbrev in sorted(abbreviations.items(), key=lambda x: len(x[0]), reverse=True):
        term = term.replace(full, abbrev)

    # Remove extra spaces
    term = ' '.join(term.split())

    # Truncate if still too long
    if len(term) > max_length:
        term = term[:max_length-3] + '...'

    return term


def load_focused_gsea_results(genotype_dir, cluster_name, focus_key):
    """
    Load focused GSEA results for a given genotype and cluster.

    Parameters:
    -----------
    genotype_dir : str
        Directory name for the genotype (e.g., 'emx1', 'nestin')
    cluster_name : str
        Cluster name (e.g., 'Mature_GC', 'GABA')
    focus_key : str
        Focus area key ('cell_death' or 'splicing')

    Returns:
    --------
    DataFrame or None
        DataFrame with Term, nes, and fdr columns, or None if file not found
    """
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
                    raise ImportError("openpyxl not installed and no CSV version found.")
    except Exception as e:
        print(f"Error loading file {file_path}: {e}")
        return None

    if 'Term' not in df.columns or 'nes' not in df.columns or 'fdr' not in df.columns:
        print(f"Warning: Required columns not found in {file_path}")
        return None

    # Keep only necessary columns
    df = df[['Term', 'nes', 'fdr']].copy()

    # Handle potential duplicate terms by keeping the one with the max absolute NES
    if df['Term'].duplicated().any():
        print(f"  - Found duplicate terms. Resolving by max absolute NES.")
        df['abs_nes'] = df['nes'].abs()
        df = df.loc[df.groupby('Term')['abs_nes'].idxmax()]
        df = df.drop(columns=['abs_nes'])

    return df


def create_focused_heatmap(focus_key, max_terms=20, fdr_threshold=0.25, no_fdr=False):
    """
    Create a focused heatmap for a specific biological category.

    Parameters:
    -----------
    focus_key : str
        Focus area key ('cell_death' or 'splicing')
    max_terms : int
        Maximum number of terms to display
    fdr_threshold : float
        FDR threshold for significance stars
    no_fdr : bool
        If True, create heatmap without FDR significance stars
    """
    focus_config = FOCUS_CONFIG[focus_key]
    title_name = focus_config['title_name']
    clusters = focus_config['clusters']

    print(f"\n{'='*70}")
    print(f"Creating heatmap for {title_name}{' (no FDR)' if no_fdr else ''}")
    print(f"{'='*70}")

    # Try each cluster name until we find one that works
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

        # If we found results for both genotypes, use this cluster
        if len(temp_results) == len(GENOTYPES):
            all_results = temp_results
            cluster_found = cluster_name
            print(f"  ✓ Found results for cluster: {cluster_name}")
            break

    if not all_results:
        print(f"\nError: No GSEA result files found for {title_name}.")
        print(f"Tried clusters: {', '.join(clusters)}")
        print(f"Please run gsea_conditions_filter.py first for these configurations.")
        return

    # Merge all dataframes on the 'Term' index
    merged_df = pd.concat(all_results, axis=1, join='outer')

    if merged_df.empty:
        print(f"No data found for {title_name}. Exiting.")
        return

    print(f"\nFound {len(merged_df)} total terms across both genotypes.")

    # Calculate mean absolute NES for ranking
    nes_cols = [f"{name}_nes" for name in GENOTYPES.keys()]
    fdr_cols = [f"{name}_fdr" for name in GENOTYPES.keys()]

    # Filter to keep only terms significant in at least one genotype
    merged_df['min_fdr'] = merged_df[fdr_cols].min(axis=1)
    merged_df['mean_abs_nes'] = merged_df[nes_cols].abs().mean(axis=1)

    # Keep only terms with at least one significant result
    significant_df = merged_df[merged_df['min_fdr'] < fdr_threshold].copy()

    if len(significant_df) == 0:
        print(f"  Warning: No terms pass FDR threshold {fdr_threshold}. Using all terms.")
        significant_df = merged_df.copy()

    # Select top terms by mean absolute NES from significant terms
    top_terms = significant_df.nlargest(max_terms, 'mean_abs_nes')

    print(f"Selected top {len(top_terms)} terms by mean absolute NES")

    # Prepare plot data
    plot_df = top_terms[nes_cols].copy()
    plot_df.columns = GENOTYPES.keys()

    # Shorten term names
    original_terms = plot_df.index.tolist()
    plot_df.index = [shorten_term_name(term) for term in plot_df.index]

    # Calculate figure dimensions based on number of terms
    n_terms = len(plot_df)
    n_cols = len(plot_df.columns)

    # Very compact sizing for publication with limited space
    fig_width = min(2.5 + (n_cols * 0.6), 5)  # Narrower width
    fig_height = min(1.5 + (n_terms * 0.18), 8)  # Reduced height per term (0.18 instead of 0.25)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Color scale centered at 0
    vmax = np.nanmax(np.abs(plot_df.values))
    vmax = max(vmax, 0.1)
    vmin = -vmax

    # Create annotation with significance stars (or empty if no_fdr)
    annot_df = pd.DataFrame(index=plot_df.index, columns=plot_df.columns, dtype=str)

    if not no_fdr:
        for idx, (short_term, orig_term) in enumerate(zip(plot_df.index, original_terms)):
            for genotype_name in GENOTYPES.keys():
                fdr_col = f"{genotype_name}_fdr"
                nes_col = f"{genotype_name}_nes"

                if orig_term in top_terms.index:
                    fdr_val = top_terms.loc[orig_term, fdr_col]
                    nes_val = top_terms.loc[orig_term, nes_col]

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

                        annot_df.loc[short_term, genotype_name] = sig
                    else:
                        annot_df.loc[short_term, genotype_name] = ''
                else:
                    annot_df.loc[short_term, genotype_name] = ''
    else:
        # Fill with empty strings for no FDR version
        annot_df = annot_df.fillna('')

    # Create heatmap with compact styling
    sns.heatmap(
        plot_df,
        cmap='RdBu_r',
        vmin=vmin,
        vmax=vmax,
        center=0,
        annot=annot_df,
        fmt='',
        linewidths=0.3,  # Thinner lines
        linecolor='lightgray',
        cbar_kws={
            'label': 'NES',
            'shrink': 0.7,
            'aspect': 15,
            'pad': 0.01
        },
        ax=ax
    )

    # Formatting - more compact
    cluster_display = cluster_found.replace('_', ' ')
    ax.set_xlabel('Genotype', fontsize=8, fontweight='bold')
    ax.set_ylabel('', fontweight='bold')
    ax.set_title(f'{title_name} in {cluster_display}\n(Mutant vs Control)',
                fontsize=9, fontweight='bold', pad=8)

    # Rotate labels
    plt.xticks(rotation=45, ha='right', fontsize=7)
    plt.yticks(rotation=0, fontsize=6.5)  # Slightly smaller y-axis labels

    # Add compact legend for significance (only if not no_fdr)
    if not no_fdr:
        legend_elements = [
            mpatches.Patch(facecolor='white', edgecolor='black', label='* FDR<0.25'),
            mpatches.Patch(facecolor='white', edgecolor='black', label='** FDR<0.01'),
            mpatches.Patch(facecolor='white', edgecolor='black', label='*** FDR<0.001')
        ]

        # Position legend outside plot area - more compact
        ax.legend(handles=legend_elements,
                 bbox_to_anchor=(1.01, 1),
                 loc='upper left',
                 frameon=True,
                 fontsize=5.5,
                 handlelength=1,
                 handletextpad=0.5)

    plt.tight_layout()

    # Save in multiple formats
    fdr_suffix = '_no_fdr' if no_fdr else ''
    base_filename = f"GSEA_Focused_{focus_key}_{cluster_found}_compact{fdr_suffix}"

    # High-resolution PNG for publication
    png_path = os.path.join(OUTPUT_DIR, f"{base_filename}.png")
    plt.savefig(png_path, dpi=600, bbox_inches='tight', facecolor='white', pad_inches=0.05)
    print(f"\n✓ Saved PNG (600 DPI): {png_path}")

    # PDF for vector graphics
    pdf_path = os.path.join(OUTPUT_DIR, f"{base_filename}.pdf")
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white', pad_inches=0.05)
    print(f"✓ Saved PDF (vector): {pdf_path}")

    plt.close()


def main():
    """
    Main function to generate focused GSEA heatmaps.
    """
    parser = argparse.ArgumentParser(
        description="Generate focused GSEA heatmaps for specific biological categories."
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
        help="Maximum number of terms to display (default: 12, optimized for compact figures)"
    )
    parser.add_argument(
        '--fdr_threshold',
        type=float,
        default=0.25,
        help="FDR threshold for significance (default: 0.25)"
    )

    args = parser.parse_args()

    print("=" * 70)
    print("Focused GSEA Heatmap Generator (Compact Publication Version)")
    print("=" * 70)
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Max terms per heatmap: {args.max_terms} (compact for limited space)")
    print(f"FDR threshold: {args.fdr_threshold} (terms must be significant in ≥1 genotype)")

    # Generate heatmaps based on focus argument
    if args.focus == 'both' or args.focus == 'cell_death':
        # Create version with FDR
        create_focused_heatmap(
            'cell_death',
            max_terms=args.max_terms,
            fdr_threshold=args.fdr_threshold,
            no_fdr=False
        )
        # Create version without FDR
        create_focused_heatmap(
            'cell_death',
            max_terms=args.max_terms,
            fdr_threshold=args.fdr_threshold,
            no_fdr=True
        )

    if args.focus == 'both' or args.focus == 'splicing':
        # Create version with FDR
        create_focused_heatmap(
            'splicing',
            max_terms=args.max_terms,
            fdr_threshold=args.fdr_threshold,
            no_fdr=False
        )
        # Create version without FDR
        create_focused_heatmap(
            'splicing',
            max_terms=args.max_terms,
            fdr_threshold=args.fdr_threshold,
            no_fdr=True
        )

    print("\n" + "=" * 70)
    print("✓ Focused heatmap generation completed!")
    print(f"✓ Output directory: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == '__main__':
    main()
