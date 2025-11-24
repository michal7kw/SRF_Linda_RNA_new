#!/usr/bin/env python3
"""
Gene Set Scoring with Custom Gene Lists - Excitatory Neurons Only

This script performs gene set scoring using custom gene lists, specifically
focusing on excitatory neurons (cell_type_L1 == 'Excitatory'). It analyzes
gene set enrichment within excitatory neuron subtypes.

Filtering Strategy:
1. Filter cells to only include excitatory neurons (cell_type_L1 == 'Excitatory')
2. Analyze subtypes within excitatory neurons using cell_type_L2_new annotations
3. Create visualizations showing enrichment patterns across excitatory subtypes

Input:
    - annotation_final.h5ad: scRNA-seq data with cell type annotations
    - gene_lists.csv: Custom gene lists with up/down regulation indicators

Output:
    - Publication-ready heatmaps for excitatory neuron subtypes
    - Summary statistics and reports
"""

import scanpy as sc
import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# --- Configuration ---
WORKING_DIR = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data"
os.chdir(WORKING_DIR)

# Set up directories
BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")
INPUT_DIR = BASE_RESULTS_DIR
ADATA_PATH = os.path.join(INPUT_DIR, "annotation_final.h5ad")

# Gene sets path (custom gene lists)
GENE_LISTS_PATH = os.path.join("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/DATA/GENE_LIST/gene_lists.csv")

# Output directory (excitatory-specific)
OUTPUT_DIR = os.path.join(BASE_RESULTS_DIR, "gene_set_enrichment_heatmaps_custom_excitatory")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Cell type filtering
CELL_TYPE_L1_FILTER = 'Excitatory'  # Only analyze excitatory neurons

# Cell type columns to analyze (within excitatory neurons)
CELL_TYPE_COLUMNS = ['cell_type_L2_new']  # Detailed subtypes within excitatory neurons

# Minimum cells per cell type to include in analysis
MIN_CELLS_PER_TYPE = 50

# Publication-ready plot settings
plt.rcParams['font.size'] = 8
plt.rcParams['axes.labelsize'] = 9
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 10
plt.rcParams['pdf.fonttype'] = 42  # TrueType fonts for editability
plt.rcParams['ps.fonttype'] = 42


def load_custom_gene_lists(gene_lists_path):
    """
    Load custom gene lists from CSV file with up/down regulation indicators.

    The CSV format is:
    - Row 1: Gene set names (biological processes)
    - Row 2: up/down indicators (correlation direction)
    - Rows 3+: Gene symbols

    Parameters:
    -----------
    gene_lists_path : str
        Path to the gene lists CSV file

    Returns:
    --------
    dict
        Nested dictionary structure:
        {
            'gene_set_name': {
                'up': [list of positively correlated genes],
                'down': [list of negatively correlated genes]
            }
        }
    """
    print(f"Loading custom gene lists from: {gene_lists_path}")

    # Read the CSV file
    df = pd.read_csv(gene_lists_path, header=None)

    # First row is gene set names
    gene_set_names = df.iloc[0].values
    # Second row is up/down indicators
    correlation_dirs = df.iloc[1].values
    # Remaining rows are genes
    gene_rows = df.iloc[2:].values

    # Build the gene sets dictionary
    gene_sets = {}

    # Process each column
    for col_idx in range(len(gene_set_names)):
        gene_set_name = gene_set_names[col_idx]
        correlation_dir = correlation_dirs[col_idx]

        # Skip if gene set name is NaN
        if pd.isna(gene_set_name):
            continue

        # Initialize gene set if not exists
        if gene_set_name not in gene_sets:
            gene_sets[gene_set_name] = {'up': [], 'down': []}

        # Collect genes from this column
        genes = []
        for row_idx in range(len(gene_rows)):
            gene = gene_rows[row_idx][col_idx]
            if pd.notna(gene) and gene != '':
                genes.append(gene)

        # Add genes to appropriate list based on correlation direction
        if correlation_dir == 'up':
            gene_sets[gene_set_name]['up'].extend(genes)
        elif correlation_dir == 'down':
            gene_sets[gene_set_name]['down'].extend(genes)
        else:
            print(f"  Warning: Unknown correlation direction '{correlation_dir}' for {gene_set_name}")

    # Remove duplicates and clean up
    for gene_set_name in gene_sets:
        gene_sets[gene_set_name]['up'] = list(set(gene_sets[gene_set_name]['up']))
        gene_sets[gene_set_name]['down'] = list(set(gene_sets[gene_set_name]['down']))

    print(f"\nLoaded {len(gene_sets)} gene sets:")
    for gene_set_name, gene_dict in gene_sets.items():
        n_up = len(gene_dict['up'])
        n_down = len(gene_dict['down'])
        print(f"  - {gene_set_name}: {n_up} up-genes, {n_down} down-genes")

    return gene_sets


def calculate_gene_set_scores_with_direction(adata, gene_sets_dict, score_type='aggregate'):
    """
    Calculate gene set scores considering positive/negative correlation directions.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object containing the gene expression data
    gene_sets_dict : dict
        Dictionary with gene set names as keys and 'up'/'down' gene lists as values
    score_type : str
        Type of scoring: 'positive', 'negative', or 'aggregate'
        - 'positive': Only score using up-regulated genes
        - 'negative': Only score using down-regulated genes
        - 'aggregate': Combine both (up genes positive, down genes negative)

    Returns:
    --------
    AnnData
        Updated AnnData object with gene set scores in .obs
    """
    print(f"\nCalculating gene set scores (type: {score_type})...")

    # Make a copy to avoid modifying the original
    adata_copy = adata.copy()

    # Calculate scores for each gene set
    for gene_set_name, gene_dict in gene_sets_dict.items():
        print(f"  Scoring gene set: {gene_set_name}")

        up_genes = gene_dict['up']
        down_genes = gene_dict['down']

        # Filter genes to only those present in the dataset
        up_genes_in_data = [gene for gene in up_genes if gene in adata_copy.var.index]
        down_genes_in_data = [gene for gene in down_genes if gene in adata_copy.var.index]

        up_missing = len(up_genes) - len(up_genes_in_data)
        down_missing = len(down_genes) - len(down_genes_in_data)

        if up_missing > 0:
            print(f"    Up-genes: {len(up_genes_in_data)}/{len(up_genes)} found ({up_missing} missing)")
        if down_missing > 0:
            print(f"    Down-genes: {len(down_genes_in_data)}/{len(down_genes)} found ({down_missing} missing)")

        if score_type == 'positive':
            # Only score using up-regulated genes
            if len(up_genes_in_data) == 0:
                print(f"    Error: No up-regulated genes found. Skipping.")
                continue

            sc.tl.score_genes(adata_copy, gene_list=up_genes_in_data,
                             score_name=f"{gene_set_name}_positive", use_raw=False)

        elif score_type == 'negative':
            # Only score using down-regulated genes
            if len(down_genes_in_data) == 0:
                print(f"    Error: No down-regulated genes found. Skipping.")
                continue

            sc.tl.score_genes(adata_copy, gene_list=down_genes_in_data,
                             score_name=f"{gene_set_name}_negative", use_raw=False)

        elif score_type == 'aggregate':
            # Calculate both scores and combine them
            up_score = None
            down_score = None

            if len(up_genes_in_data) > 0:
                sc.tl.score_genes(adata_copy, gene_list=up_genes_in_data,
                                 score_name='_temp_up', use_raw=False)
                up_score = adata_copy.obs['_temp_up'].values

            if len(down_genes_in_data) > 0:
                sc.tl.score_genes(adata_copy, gene_list=down_genes_in_data,
                                 score_name='_temp_down', use_raw=False)
                down_score = adata_copy.obs['_temp_down'].values

            # Combine scores: up genes contribute positively, down genes negatively
            if up_score is not None and down_score is not None:
                # Aggregate = (up_score) - (down_score)
                # This means high expression of up-genes increases score
                # and high expression of down-genes decreases score
                aggregate_score = up_score - down_score
            elif up_score is not None:
                aggregate_score = up_score
            elif down_score is not None:
                aggregate_score = -down_score
            else:
                print(f"    Error: No genes found for {gene_set_name}. Skipping.")
                continue

            adata_copy.obs[f"{gene_set_name}_aggregate"] = aggregate_score

            # Clean up temporary columns
            if '_temp_up' in adata_copy.obs.columns:
                adata_copy.obs.drop('_temp_up', axis=1, inplace=True)
            if '_temp_down' in adata_copy.obs.columns:
                adata_copy.obs.drop('_temp_down', axis=1, inplace=True)

        else:
            raise ValueError(f"Unknown score_type: {score_type}")

    return adata_copy


def create_individual_geneset_heatmaps(adata, gene_sets_dict, cell_type_column,
                                       split_mode='condition', score_type='aggregate',
                                       output_dir=OUTPUT_DIR):
    """
    Create separate heatmaps for each gene set with side-by-side condition/genotype comparisons.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names as keys
    cell_type_column : str
        Column name containing cell type annotations
    split_mode : str
        'condition' for Control/Mutant split (2 columns)
        'genotype_condition' for Mutant-Nestin, Control-Nestin, Mutant-Emx1, Control-Emx1 (4 columns)
    score_type : str
        Type of score: 'positive', 'negative', or 'aggregate'
    output_dir : str
        Directory to save the plots
    """
    print(f"\nCreating individual gene set heatmaps for: {cell_type_column}")
    print(f"  Split mode: {split_mode}, Score type: {score_type}")

    # Check if column exists
    if cell_type_column not in adata.obs.columns:
        print(f"  Error: Column '{cell_type_column}' not found")
        return

    # Get cell types with sufficient cells
    cell_type_counts = adata.obs[cell_type_column].value_counts()
    valid_cell_types = cell_type_counts[cell_type_counts >= MIN_CELLS_PER_TYPE].index.tolist()

    print(f"  Including {len(valid_cell_types)} cell types (excitatory subtypes)")

    # Filter data
    adata_filtered = adata[adata.obs[cell_type_column].isin(valid_cell_types)].copy()

    if len(adata_filtered) == 0:
        print(f"  No data available")
        return

    # Create a separate heatmap for each gene set
    for gene_set_name in gene_sets_dict.keys():
        print(f"\n  Processing gene set: {gene_set_name}")
        score_column = f"{gene_set_name}_{score_type}"

        if score_column not in adata_filtered.obs.columns:
            print(f"    Score column '{score_column}' not found, skipping")
            continue

        # Prepare data based on split mode
        heatmap_data = []

        if split_mode == 'condition' and 'condition' in adata_filtered.obs.columns:
            # Two columns: Control and Mutant
            conditions = ['Control', 'Mutant']

            for cell_type in sorted(valid_cell_types):
                row_data = {'Cell_Type': cell_type}

                for condition in conditions:
                    mask = (adata_filtered.obs[cell_type_column] == cell_type) & \
                           (adata_filtered.obs['condition'] == condition)
                    cell_data = adata_filtered[mask]

                    if len(cell_data) > 0:
                        row_data[condition] = cell_data.obs[score_column].mean()
                        row_data[f'{condition}_n'] = len(cell_data)
                    else:
                        row_data[condition] = np.nan
                        row_data[f'{condition}_n'] = 0

                heatmap_data.append(row_data)

        elif split_mode == 'genotype_condition' and 'genotype' in adata_filtered.obs.columns and 'condition' in adata_filtered.obs.columns:
            # Four columns: Control-Emx1, Mutant-Emx1, Control-Nestin, Mutant-Nestin
            genotypes = sorted(adata_filtered.obs['genotype'].unique())
            conditions = ['Control', 'Mutant']

            for cell_type in sorted(valid_cell_types):
                row_data = {'Cell_Type': cell_type}

                for genotype in genotypes:
                    for condition in conditions:
                        col_name = f"{condition}-{genotype}"
                        mask = (adata_filtered.obs[cell_type_column] == cell_type) & \
                               (adata_filtered.obs['genotype'] == genotype) & \
                               (adata_filtered.obs['condition'] == condition)
                        cell_data = adata_filtered[mask]

                        if len(cell_data) > 0:
                            row_data[col_name] = cell_data.obs[score_column].mean()
                            row_data[f'{col_name}_n'] = len(cell_data)
                        else:
                            row_data[col_name] = np.nan
                            row_data[f'{col_name}_n'] = 0

                heatmap_data.append(row_data)
        else:
            print(f"    Error: Required columns not found for split mode '{split_mode}'")
            continue

        # Convert to DataFrame
        heatmap_df = pd.DataFrame(heatmap_data)

        if heatmap_df.empty:
            print(f"    No data available")
            continue

        # Get score columns (exclude Cell_Type and _n columns)
        score_columns = [col for col in heatmap_df.columns
                        if col != 'Cell_Type' and not col.endswith('_n')]

        if not score_columns:
            print(f"    No score columns found")
            continue

        # Set up the data for heatmap
        plot_data = heatmap_df.set_index('Cell_Type')[score_columns]

        print(f"    Original score range: [{plot_data.min().min():.3f}, {plot_data.max().max():.3f}]")

        # Calculate figure size based on number of conditions and cell types
        n_conditions = len(score_columns)
        n_celltypes = len(plot_data)

        # Width: ~0.6 inches per condition + 3 for labels and colorbar
        fig_width = min(n_conditions * 0.8 + 3, 10)
        # Height: ~0.35 inches per cell type + 1.5 for title/labels
        fig_height = min(n_celltypes * 0.35 + 1.5, 12)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        # Use actual data range but compress it to highlight differences
        # Calculate range based on actual data
        valid_scores = plot_data.values[~np.isnan(plot_data.values)]

        if len(valid_scores) > 0:
            data_min = np.min(valid_scores)
            data_max = np.max(valid_scores)
            data_range = data_max - data_min

            # Use the actual min/max but add small padding for better visualization
            vmin = data_min - 0.05 * data_range
            vmax = data_max + 0.05 * data_range

            print(f"    Color scale: [{vmin:.3f}, {vmax:.3f}]")
        else:
            vmin, vmax = -0.2, 0.2

        # Create heatmap with white-to-red color scale
        # 'Reds' goes from white (low) to dark red (high)
        sns.heatmap(plot_data, ax=ax, cmap='Reds',
                   vmin=vmin, vmax=vmax,
                   cbar_kws={'label': f'Mean Gene Set Score ({score_type})', 'shrink': 0.8},
                   linewidths=0.5, linecolor='lightgray',
                   annot=False, fmt='.3f',
                   cbar=True)

        # Improve labels
        ax.set_xlabel('Condition', fontweight='bold', fontsize=9)
        ax.set_ylabel('Excitatory Neuron Subtype', fontweight='bold', fontsize=9)

        plt.xticks(rotation=45, ha='right', fontsize=8)
        plt.yticks(rotation=0, fontsize=8)

        # Add title
        title = f'{gene_set_name} ({score_type})\nExcitatory Neurons - {cell_type_column.replace("_", " ").title()}'
        plt.title(title, fontsize=10, fontweight='bold', pad=15)

        plt.tight_layout()

        # Save files
        safe_gene_set = gene_set_name.replace('/', '_').replace(' ', '_').replace('(', '').replace(')', '')
        base_filename = f'{safe_gene_set}_{score_type}_excitatory_{cell_type_column}_{split_mode}'

        plt.savefig(os.path.join(output_dir, f'{base_filename}.png'),
                   bbox_inches='tight', dpi=300)
        plt.savefig(os.path.join(output_dir, f'{base_filename}.pdf'),
                   bbox_inches='tight', dpi=300)
        plt.close()

        print(f"    Saved: {base_filename}")

        # Also save the data
        save_df = heatmap_df.copy()
        save_df.to_csv(os.path.join(output_dir, f'{base_filename}_data.csv'), index=False)
        print(f"    Saved data: {base_filename}_data.csv")


def create_summary_report(adata, gene_sets_dict, score_type, output_dir):
    """
    Create a summary report of gene set scores.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names
    score_type : str
        Type of score calculated
    output_dir : str
        Directory to save the report
    """
    print("\nCreating summary report...")

    summary_data = []

    for gene_set_name in gene_sets_dict.keys():
        score_column = f"{gene_set_name}_{score_type}"

        if score_column not in adata.obs.columns:
            continue

        scores = adata.obs[score_column]

        summary_stats = {
            'Gene_Set': gene_set_name,
            'Score_Type': score_type,
            'N_Cells': len(scores),
            'Mean_Score': scores.mean(),
            'Std_Score': scores.std(),
            'Median_Score': scores.median(),
            'Min_Score': scores.min(),
            'Max_Score': scores.max(),
            'N_Up_Genes': len(gene_sets_dict[gene_set_name]['up']),
            'N_Down_Genes': len(gene_sets_dict[gene_set_name]['down'])
        }
        summary_data.append(summary_stats)

    # Convert to DataFrame and save
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(output_dir, f'gene_set_scores_summary_excitatory_{score_type}.csv'),
                     index=False)

    print(f"  Saved: gene_set_scores_summary_excitatory_{score_type}.csv")


def main():
    """Main function to run the gene set scoring and heatmap visualization pipeline."""

    print("=" * 70)
    print("Gene Set Scoring - EXCITATORY NEURONS ONLY")
    print("=" * 70)

    # Load the scRNA-seq data
    if not os.path.exists(ADATA_PATH):
        print(f"Error: Data file not found at {ADATA_PATH}")
        sys.exit(1)

    print(f"\nLoading scRNA data from: {ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)
    print(f"Loaded data: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"Available obs columns: {list(adata.obs.columns)}")

    # Check if cell_type_L1 column exists
    if 'cell_type_L1' not in adata.obs.columns:
        print(f"Error: 'cell_type_L1' column not found in data")
        sys.exit(1)

    print(f"\nAvailable cell_type_L1 categories: {list(adata.obs['cell_type_L1'].unique())}")

    # Filter for excitatory neurons only
    print(f"\nFiltering for excitatory neurons (cell_type_L1 == '{CELL_TYPE_L1_FILTER}')...")
    adata_excitatory = adata[adata.obs['cell_type_L1'] == CELL_TYPE_L1_FILTER].copy()
    print(f"Filtered data: {adata_excitatory.n_obs} excitatory neurons ({adata_excitatory.n_obs/adata.n_obs*100:.1f}% of total)")

    if adata_excitatory.n_obs == 0:
        print(f"Error: No excitatory neurons found in data")
        sys.exit(1)

    # Check available subtypes within excitatory neurons
    for col in CELL_TYPE_COLUMNS:
        if col in adata_excitatory.obs.columns:
            print(f"\nExcitatory neuron subtypes in '{col}':")
            subtype_counts = adata_excitatory.obs[col].value_counts()
            for subtype, count in subtype_counts.items():
                print(f"  - {subtype}: {count} cells")

    # Load custom gene sets
    if not os.path.exists(GENE_LISTS_PATH):
        print(f"Error: Gene lists file not found at {GENE_LISTS_PATH}")
        sys.exit(1)

    print("\n" + "=" * 70)
    gene_sets_dict = load_custom_gene_lists(GENE_LISTS_PATH)
    print("=" * 70)

    if not gene_sets_dict:
        print("No gene sets loaded. Exiting.")
        sys.exit(1)

    # Process all three score types
    score_types = ['positive', 'negative', 'aggregate']

    for score_type in score_types:
        print("\n" + "=" * 70)
        print(f"Processing Score Type: {score_type.upper()}")
        print("=" * 70)

        # Calculate gene set scores
        adata_scored = calculate_gene_set_scores_with_direction(
            adata_excitatory, gene_sets_dict, score_type=score_type
        )

        # Create individual gene set heatmaps for each cell type column
        print("\n" + "=" * 70)
        print(f"Creating Individual Gene Set Heatmaps ({score_type})")
        print("=" * 70)

        for cell_type_col in CELL_TYPE_COLUMNS:
            if cell_type_col in adata_scored.obs.columns:
                # Create heatmaps with 2-column split (Control/Mutant)
                create_individual_geneset_heatmaps(
                    adata_scored, gene_sets_dict, cell_type_col,
                    split_mode='condition', score_type=score_type,
                    output_dir=OUTPUT_DIR
                )

                # Create heatmaps with 4-column split (Control-Nestin, Mutant-Nestin, Control-Emx1, Mutant-Emx1)
                if 'genotype' in adata_scored.obs.columns:
                    create_individual_geneset_heatmaps(
                        adata_scored, gene_sets_dict, cell_type_col,
                        split_mode='genotype_condition', score_type=score_type,
                        output_dir=OUTPUT_DIR
                    )
            else:
                print(f"\nWarning: Column '{cell_type_col}' not found in data")

        # Create summary report
        print("\n" + "=" * 70)
        create_summary_report(adata_scored, gene_sets_dict, score_type, OUTPUT_DIR)

    print("\n" + "=" * 70)
    print("Pipeline completed successfully!")
    print(f"Results saved to: {OUTPUT_DIR}")
    print(f"Total excitatory neurons analyzed: {adata_excitatory.n_obs}")
    print("=" * 70)


if __name__ == "__main__":
    main()
