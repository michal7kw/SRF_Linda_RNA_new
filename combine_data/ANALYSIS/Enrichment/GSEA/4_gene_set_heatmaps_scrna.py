#!/usr/bin/env python3
"""
Gene Set Scoring and Heatmap Visualization for scRNA Data

This script performs gene set scoring at single cell level and creates publication-ready
heatmap visualizations showing gene set scores across different cell types (clusters).
Focuses on cell death and splicing-related gene sets.

Input:
    - annotation_final.h5ad: scRNA-seq data with cell type annotations
    - import_ready.csv: Gene sets for scoring

Output:
    - Publication-ready heatmaps for L1 and L2 cell type annotations
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

# Gene sets paths (MSigDB gene sets)
GENE_SETS_DIR = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA"
GENE_SETS_FILES = [
    os.path.join(GENE_SETS_DIR, "GOBP_RNA_SPLICING.v2025.1.Mm.grp"),
    os.path.join(GENE_SETS_DIR, "HALLMARK_APOPTOSIS.v2025.1.Mm.grp")
]

# Output directory
OUTPUT_DIR = os.path.join(BASE_RESULTS_DIR, "gene_set_enrichment_heatmaps")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Cell type columns to analyze
CELL_TYPE_COLUMNS = ['cell_type_L1', 'cell_type_L2_new']

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


def load_gene_sets_from_grp(gene_sets_files):
    """
    Load gene sets from MSigDB .grp files.

    Parameters:
    -----------
    gene_sets_files : list
        List of paths to .grp files

    Returns:
    --------
    dict
        Dictionary with gene set names as keys and lists of genes as values
    """
    print("Loading gene sets from MSigDB .grp files...")

    gene_sets_dict = {}

    for gene_set_path in gene_sets_files:
        if not os.path.exists(gene_set_path):
            print(f"  Warning: File not found: {gene_set_path}")
            continue

        print(f"  Reading: {os.path.basename(gene_set_path)}")

        with open(gene_set_path, 'r') as f:
            lines = f.readlines()

        # First non-comment line is the gene set name
        gene_set_name = None
        genes = []

        for line in lines:
            line = line.strip()
            if not line:
                continue
            # Skip comment lines (starting with #)
            if line.startswith('#'):
                continue
            # First non-comment line is the gene set name
            if gene_set_name is None:
                gene_set_name = line
                continue
            # All other lines are genes
            genes.append(line)

        if gene_set_name and genes:
            gene_sets_dict[gene_set_name] = genes
            print(f"    Loaded: {gene_set_name} with {len(genes)} genes")
        else:
            print(f"    Warning: Could not parse gene set from {gene_set_path}")

    print(f"\nTotal gene sets loaded: {len(gene_sets_dict)}")
    for gene_set, genes in gene_sets_dict.items():
        print(f"  - {gene_set}: {len(genes)} genes")

    return gene_sets_dict


def calculate_gene_set_scores(adata, gene_sets_dict):
    """
    Calculate gene set scores using scanpy's score_genes function.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object containing the gene expression data
    gene_sets_dict : dict
        Dictionary with gene set names as keys and lists of genes as values

    Returns:
    --------
    AnnData
        Updated AnnData object with gene set scores in .obs
    """
    print("\nCalculating gene set scores...")

    # Make a copy to avoid modifying the original
    adata_copy = adata.copy()

    # Calculate scores for each gene set
    for gene_set_name, genes in gene_sets_dict.items():
        print(f"  Scoring gene set: {gene_set_name}")

        # Filter genes to only those present in the dataset
        genes_in_data = [gene for gene in genes if gene in adata_copy.var.index]
        genes_missing = [gene for gene in genes if gene not in adata_copy.var.index]

        if len(genes_missing) > 0:
            print(f"    Warning: {len(genes_missing)}/{len(genes)} genes not found in dataset")
            if len(genes_missing) <= 5:
                print(f"    Missing genes: {genes_missing}")

        if len(genes_in_data) == 0:
            print(f"    Error: No genes from {gene_set_name} found in dataset. Skipping.")
            continue

        print(f"    Using {len(genes_in_data)}/{len(genes)} genes")

        # Calculate gene set score using scanpy
        # Use the normalized, log-transformed data (default X)
        sc.tl.score_genes(adata_copy, gene_list=genes_in_data,
                         score_name=f"{gene_set_name}_score", use_raw=False)

    return adata_copy


def create_celltype_heatmap(adata, gene_sets_dict, cell_type_column,
                            condition_split=True, output_dir=OUTPUT_DIR):
    """
    Create publication-ready heatmap showing gene set scores across cell types.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names as keys and lists of genes as values
    cell_type_column : str
        Column name containing cell type annotations
    condition_split : bool
        If True, create separate rows for Control and Mutant conditions
    output_dir : str
        Directory to save the plots
    """
    print(f"\nCreating heatmap for cell type column: {cell_type_column}")

    # Check if column exists
    if cell_type_column not in adata.obs.columns:
        print(f"  Error: Column '{cell_type_column}' not found in adata.obs")
        return

    # Get cell types with sufficient cells
    cell_type_counts = adata.obs[cell_type_column].value_counts()
    valid_cell_types = cell_type_counts[cell_type_counts >= MIN_CELLS_PER_TYPE].index.tolist()

    print(f"  Including {len(valid_cell_types)} cell types with >={MIN_CELLS_PER_TYPE} cells")
    print(f"  Cell types: {', '.join(valid_cell_types)}")

    # Filter data
    adata_filtered = adata[adata.obs[cell_type_column].isin(valid_cell_types)].copy()

    if len(adata_filtered) == 0:
        print(f"  No data available after filtering")
        return

    # Prepare data for heatmap
    heatmap_data = []

    if condition_split and 'condition' in adata_filtered.obs.columns:
        # Create separate entries for Control and Mutant
        conditions = ['Control', 'Mutant']

        for cell_type in sorted(valid_cell_types):
            for condition in conditions:
                # Filter for this cell type and condition
                mask = (adata_filtered.obs[cell_type_column] == cell_type) & \
                       (adata_filtered.obs['condition'] == condition)
                cell_data = adata_filtered[mask]

                if len(cell_data) == 0:
                    continue

                # Calculate mean scores for each gene set
                row_data = {
                    'Cell_Type': f"{cell_type}_{condition}",
                    'Cell_Type_Base': cell_type,
                    'Condition': condition,
                    'N_Cells': len(cell_data)
                }

                for gene_set_name in gene_sets_dict.keys():
                    score_column = f"{gene_set_name}_score"
                    if score_column in cell_data.obs.columns:
                        row_data[gene_set_name] = cell_data.obs[score_column].mean()

                heatmap_data.append(row_data)
    else:
        # No condition split, just cell types
        for cell_type in sorted(valid_cell_types):
            cell_data = adata_filtered[adata_filtered.obs[cell_type_column] == cell_type]

            # Calculate mean scores for each gene set
            row_data = {
                'Cell_Type': cell_type,
                'N_Cells': len(cell_data)
            }

            for gene_set_name in gene_sets_dict.keys():
                score_column = f"{gene_set_name}_score"
                if score_column in cell_data.obs.columns:
                    row_data[gene_set_name] = cell_data.obs[score_column].mean()

            heatmap_data.append(row_data)

    # Convert to DataFrame
    heatmap_df = pd.DataFrame(heatmap_data)

    if heatmap_df.empty:
        print(f"  No data available for heatmap")
        return

    # Prepare data for plotting
    gene_set_columns = [col for col in heatmap_df.columns
                       if col not in ['Cell_Type', 'Cell_Type_Base', 'Condition', 'N_Cells']]

    if not gene_set_columns:
        print(f"  No gene set scores found")
        return

    # Set index for heatmap
    plot_data = heatmap_df.set_index('Cell_Type')[gene_set_columns]

    # Calculate appropriate figure size
    # Width: 0.4 inches per gene set + 3 inches for labels
    # Height: 0.35 inches per cell type + 1.5 inches for labels
    fig_width = min(0.4 * len(gene_set_columns) + 3, 10)  # Cap at 10 inches
    fig_height = min(0.35 * len(plot_data) + 1.5, 12)  # Cap at 12 inches

    # Create heatmap
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Use a diverging colormap centered at 0
    sns.heatmap(plot_data, ax=ax, cmap='RdBu_r', center=0,
               cbar_kws={'label': 'Mean Gene Set Score', 'shrink': 0.8},
               linewidths=0.5, linecolor='lightgray',
               annot=False, fmt='.2f')

    # Improve labels
    ax.set_xlabel('Gene Sets', fontweight='bold')
    ax.set_ylabel('Cell Types', fontweight='bold')

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    # Add title
    title_suffix = " (by Condition)" if condition_split else ""
    plt.title(f'Gene Set Enrichment Scores: {cell_type_column}{title_suffix}',
             fontweight='bold', pad=15)

    plt.tight_layout()

    # Save in multiple formats for publication
    base_filename = f'gene_set_heatmap_{cell_type_column}'
    if condition_split:
        base_filename += '_by_condition'

    # Save as PNG (high resolution)
    plt.savefig(os.path.join(output_dir, f'{base_filename}.png'),
               bbox_inches='tight', dpi=300)

    # Save as PDF (vector format for publication)
    plt.savefig(os.path.join(output_dir, f'{base_filename}.pdf'),
               bbox_inches='tight', dpi=300)

    plt.close()

    print(f"  Saved: {base_filename}.png and .pdf")

    # Also save the data as CSV
    # Add cell counts to the saved data
    save_df = heatmap_df.copy()
    save_df.to_csv(os.path.join(output_dir, f'{base_filename}_data.csv'), index=False)
    print(f"  Saved data: {base_filename}_data.csv")


def create_individual_geneset_heatmaps(adata, gene_sets_dict, cell_type_column,
                                       split_mode='condition', output_dir=OUTPUT_DIR):
    """
    Create separate heatmaps for each gene set with side-by-side condition/genotype comparisons.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names as keys and lists of genes as values
    cell_type_column : str
        Column name containing cell type annotations
    split_mode : str
        'condition' for Control/Mutant split (2 columns)
        'genotype_condition' for Mutant-Nestin, Control-Nestin, Mutant-Emx1, Control-Emx1 (4 columns)
    output_dir : str
        Directory to save the plots
    """
    print(f"\nCreating individual gene set heatmaps for: {cell_type_column}")
    print(f"  Split mode: {split_mode}")

    # Check if column exists
    if cell_type_column not in adata.obs.columns:
        print(f"  Error: Column '{cell_type_column}' not found")
        return

    # Get cell types with sufficient cells
    cell_type_counts = adata.obs[cell_type_column].value_counts()
    valid_cell_types = cell_type_counts[cell_type_counts >= MIN_CELLS_PER_TYPE].index.tolist()

    print(f"  Including {len(valid_cell_types)} cell types")

    # Filter data
    adata_filtered = adata[adata.obs[cell_type_column].isin(valid_cell_types)].copy()

    if len(adata_filtered) == 0:
        print(f"  No data available")
        return

    # Create a separate heatmap for each gene set
    for gene_set_name in gene_sets_dict.keys():
        print(f"\n  Processing gene set: {gene_set_name}")
        score_column = f"{gene_set_name}_score"

        if score_column not in adata_filtered.obs.columns:
            print(f"    Score column not found, skipping")
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
            # Four columns: Mutant-Nestin, Control-Nestin, Mutant-Emx1, Control-Emx1
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
                   cbar_kws={'label': 'Mean Gene Set Score', 'shrink': 0.8},
                   linewidths=0.5, linecolor='lightgray',
                   annot=False, fmt='.3f',
                   cbar=True)

        # Improve labels
        ax.set_xlabel('Condition', fontweight='bold', fontsize=9)
        ax.set_ylabel('Cell Type', fontweight='bold', fontsize=9)

        plt.xticks(rotation=45, ha='right', fontsize=8)
        plt.yticks(rotation=0, fontsize=8)

        # Add title
        title = f'{gene_set_name}\n{cell_type_column.replace("_", " ").title()}'
        plt.title(title, fontsize=10, fontweight='bold', pad=15)

        plt.tight_layout()

        # Save files
        safe_gene_set = gene_set_name.replace('/', '_').replace(' ', '_')
        base_filename = f'{safe_gene_set}_{cell_type_column}_{split_mode}'

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


def create_summary_report(adata, gene_sets_dict, output_dir):
    """
    Create a summary report of gene set scores.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names as keys and lists of genes as values
    output_dir : str
        Directory to save the report
    """
    print("\nCreating summary report...")

    summary_data = []

    for gene_set_name in gene_sets_dict.keys():
        score_column = f"{gene_set_name}_score"

        if score_column not in adata.obs.columns:
            continue

        scores = adata.obs[score_column]

        summary_stats = {
            'Gene_Set': gene_set_name,
            'N_Cells': len(scores),
            'Mean_Score': scores.mean(),
            'Std_Score': scores.std(),
            'Median_Score': scores.median(),
            'Min_Score': scores.min(),
            'Max_Score': scores.max(),
            'N_Genes_in_Set': len(gene_sets_dict[gene_set_name]),
            'N_Genes_Found': len([g for g in gene_sets_dict[gene_set_name]
                                 if g in adata.var.index])
        }
        summary_data.append(summary_stats)

    # Convert to DataFrame and save
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(output_dir, 'gene_set_scores_summary.csv'),
                     index=False)

    print(f"  Saved: gene_set_scores_summary.csv")

    # Create text report
    with open(os.path.join(output_dir, 'gene_set_scores_report.txt'), 'w') as f:
        f.write("Gene Set Scoring Summary Report - scRNA Data\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"Total gene sets analyzed: {len(gene_sets_dict)}\n")
        f.write(f"Total cells: {len(adata)}\n")
        f.write(f"Total genes in dataset: {adata.n_vars}\n\n")

        f.write("Gene Sets Overview:\n")
        f.write("-" * 40 + "\n")
        for gene_set_name, genes in gene_sets_dict.items():
            score_column = f"{gene_set_name}_score"
            if score_column in adata.obs.columns:
                genes_found = len([g for g in genes if g in adata.var.index])
                mean_score = adata.obs[score_column].mean()
                f.write(f"{gene_set_name}:\n")
                f.write(f"  Genes: {genes_found}/{len(genes)} found\n")
                f.write(f"  Mean score: {mean_score:.4f}\n")
            else:
                f.write(f"{gene_set_name}: SKIPPED\n")

        f.write(f"\nOutput directory: {output_dir}\n")

    print(f"  Saved: gene_set_scores_report.txt")


def main():
    """Main function to run the gene set scoring and heatmap visualization pipeline."""

    print("=" * 70)
    print("Gene Set Scoring and Heatmap Visualization for scRNA Data")
    print("=" * 70)

    # Load the scRNA-seq data
    if not os.path.exists(ADATA_PATH):
        print(f"Error: Data file not found at {ADATA_PATH}")
        sys.exit(1)

    print(f"\nLoading scRNA data from: {ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)
    print(f"Loaded data: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"Available obs columns: {list(adata.obs.columns)}")

    # Load gene sets from MSigDB .grp files
    print("\n" + "=" * 70)
    gene_sets_dict = load_gene_sets_from_grp(GENE_SETS_FILES)
    print("=" * 70)

    if not gene_sets_dict:
        print("No gene sets loaded. Exiting.")
        sys.exit(1)

    # Calculate gene set scores
    print("\n" + "=" * 70)
    adata_scored = calculate_gene_set_scores(adata, gene_sets_dict)
    print("=" * 70)

    # Create individual gene set heatmaps for each cell type column
    print("\n" + "=" * 70)
    print("Creating Individual Gene Set Heatmaps")
    print("=" * 70)

    for cell_type_col in CELL_TYPE_COLUMNS:
        if cell_type_col in adata_scored.obs.columns:
            # Create heatmaps with 2-column split (Control/Mutant)
            create_individual_geneset_heatmaps(
                adata_scored, gene_sets_dict, cell_type_col,
                split_mode='condition', output_dir=OUTPUT_DIR
            )

            # Create heatmaps with 4-column split (Control-Nestin, Mutant-Nestin, Control-Emx1, Mutant-Emx1)
            if 'genotype' in adata_scored.obs.columns:
                create_individual_geneset_heatmaps(
                    adata_scored, gene_sets_dict, cell_type_col,
                    split_mode='genotype_condition', output_dir=OUTPUT_DIR
                )
        else:
            print(f"\nWarning: Column '{cell_type_col}' not found in data")

    # Create summary report
    print("\n" + "=" * 70)
    create_summary_report(adata_scored, gene_sets_dict, OUTPUT_DIR)

    print("\n" + "=" * 70)
    print("Pipeline completed successfully!")
    print(f"Results saved to: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
