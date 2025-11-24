#!/usr/bin/env python3
"""
Script to generate marker gene expression dotplot for cell_type_L2_new annotations
Creates a comprehensive dotplot showing top biomarkers for each cell type from the annotation_final.h5ad

This script:
1. Loads annotation_final.h5ad
2. Loads biomarker genes from DEGs_cell_type_L2_newFC_0_25/biomarkers/sig_deg_lists
3. Creates a dotplot with top markers per cell type
"""

import os
import sys
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Configuration
ANNOTATION_H5AD_PATH = '/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/annotation_final.h5ad'
BIOMARKER_BASE_DIR = '/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/DEGs_cell_type_L2_newFC_0_25/biomarkers/sig_deg_lists'
OUTPUT_DIR = '/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/ADDITIONAL_PLOTS/dotplots/L2'
CELL_TYPE_COLUMN = 'cell_type_L2_new'
# Generate multiple plots with different numbers of top markers
TOP_N_GENES_LIST = [3, 5]  # Will create separate plots for each

def load_cell_type_markers(biomarker_dir, cell_type, top_n=5):
    """
    Load top marker genes for a specific cell type

    Args:
        biomarker_dir: Base directory containing biomarker files
        cell_type: Cell type name (matching folder name)
        top_n: Number of top genes to extract

    Returns:
        List of top marker gene names
    """
    # Handle cell types with special characters (spaces, parentheses)
    # Map from cell_type_L2_new values to folder names
    cell_type_folder_map = {
        'Lamp5 (DG)': 'Lamp5_DG',
        'VIP (DG)': 'VIP_DG',
        'ENT CTX': 'ENT_CTX',
        'Mature GC': 'Mature_GC',
        'Immature GC': 'Immature_GC',
        'Mossy cells': 'Mossy_cells'
    }

    # Get folder name (use mapping if exists, otherwise use cell_type as-is)
    folder_name = cell_type_folder_map.get(cell_type, cell_type)

    # Construct file path
    cell_type_dir = Path(biomarker_dir) / folder_name
    csv_file = cell_type_dir / f'Cluster_{cell_type}_vs_Rest_up_significant.csv'

    if not csv_file.exists():
        print(f"  ⚠️  Warning: Could not find {csv_file}")
        return []

    # Read CSV file
    try:
        df = pd.read_csv(csv_file)

        # Get top N genes by score (assuming first column is gene names)
        if 'names' in df.columns:
            top_genes = df.nsmallest(top_n, 'scores')['names'].tolist() if 'scores' in df.columns else df['names'].head(top_n).tolist()
        else:
            # If no 'names' column, assume first column is gene names
            top_genes = df.iloc[:top_n, 0].tolist()

        print(f"  ✓ Loaded {len(top_genes)} markers for {cell_type}")
        return top_genes

    except Exception as e:
        print(f"  ❌ Error loading markers for {cell_type}: {e}")
        return []


def create_cell_type_dotplot(adata, marker_genes_dict, output_dir, top_n_genes, cell_type_column='cell_type_L2_new'):
    """
    Create dotplot showing marker gene expression across all cell types

    Args:
        adata: AnnData object with cell type annotations
        marker_genes_dict: Dictionary mapping cell type to marker genes
        output_dir: Directory to save output plots
        top_n_genes: Number of top genes used (for filename)
        cell_type_column: Column name containing cell type annotations
    """

    print(f"\n{'='*80}")
    print(f"Creating dotplot for {cell_type_column} (top {top_n_genes} genes per cell type)")
    print(f"{'='*80}")

    print(f"  Total cells: {adata.n_obs}")
    print(f"  Total genes: {adata.n_vars}")

    # Get unique cell types in the data
    cell_types_in_data = sorted(adata.obs[cell_type_column].unique())
    print(f"  Cell types in data: {len(cell_types_in_data)}")

    # Collect all marker genes in order of cell types
    all_marker_genes = []
    cell_type_gene_count = {}

    for cell_type in cell_types_in_data:
        if cell_type in marker_genes_dict:
            genes = marker_genes_dict[cell_type]
            all_marker_genes.extend(genes)
            cell_type_gene_count[cell_type] = len(genes)
        else:
            print(f"  ⚠️  No markers found for: {cell_type}")
            cell_type_gene_count[cell_type] = 0

    # Remove duplicates while preserving order
    unique_marker_genes = []
    seen = set()
    for gene in all_marker_genes:
        if gene not in seen:
            unique_marker_genes.append(gene)
            seen.add(gene)

    # Filter genes that exist in the dataset
    available_genes = [g for g in unique_marker_genes if g in adata.var_names]
    missing_genes = [g for g in unique_marker_genes if g not in adata.var_names]

    print(f"\n  Gene availability:")
    print(f"    Total unique marker genes: {len(unique_marker_genes)}")
    print(f"    Available in dataset: {len(available_genes)}")
    print(f"    Missing from dataset: {len(missing_genes)}")

    if len(missing_genes) > 0 and len(missing_genes) <= 20:
        print(f"    Missing genes: {', '.join(missing_genes[:20])}")

    if len(available_genes) == 0:
        print(f"  ❌ Error: No marker genes found in dataset!")
        return

    # Create dotplot with appropriate sizing
    n_cell_types = len(cell_types_in_data)
    n_genes = len(available_genes)

    # Adjust figure size based on number of genes and cell types
    # For L2 with 34 cell types, we need more vertical space
    # More space per gene for better readability
    if n_genes <= 100:
        # Top 3 markers: ~100 genes
        fig_width = max(30, n_genes * 0.65)
        fig_height = max(16, n_cell_types * 0.6)
        gene_fontsize = 16
        celltype_fontsize = 16
        dot_max = 0.7
        dot_min = 0.0
    elif n_genes <= 170:
        # Top 5 markers: ~165 genes
        fig_width = max(40, n_genes * 0.60)
        fig_height = max(18, n_cell_types * 0.65)
        gene_fontsize = 14
        celltype_fontsize = 15
        dot_max = 0.7
        dot_min = 0.0
    else:
        # Top 10 markers: ~300 genes
        fig_width = max(50, n_genes * 0.50)
        fig_height = max(20, n_cell_types * 0.7)
        gene_fontsize = 12
        celltype_fontsize = 14
        dot_max = 0.7
        dot_min = 0.0

    print(f"\n  Creating dotplot...")
    print(f"    Figure size: {fig_width:.1f} x {fig_height:.1f}")
    print(f"    Gene label font size: {gene_fontsize}")
    print(f"    Cell type label font size: {celltype_fontsize}")

    try:
        # Set font sizes for better visibility - adjusted based on plot size
        plt.rcParams['font.size'] = gene_fontsize
        plt.rcParams['axes.labelsize'] = gene_fontsize + 4
        plt.rcParams['axes.titlesize'] = gene_fontsize + 6
        plt.rcParams['xtick.labelsize'] = gene_fontsize
        plt.rcParams['ytick.labelsize'] = celltype_fontsize
        plt.rcParams['legend.fontsize'] = gene_fontsize + 2
        plt.rcParams['legend.title_fontsize'] = gene_fontsize + 4

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

        # Prepare var_group_positions and labels to show gene groups by cell type
        var_group_positions = []
        var_group_labels = []
        start_idx = 0

        for cell_type in cell_types_in_data:
            if cell_type in marker_genes_dict:
                genes = marker_genes_dict[cell_type]
                # Count how many of these genes are actually in available_genes
                ct_genes_in_data = [g for g in genes if g in available_genes]
                if len(ct_genes_in_data) > 0:
                    # Find the indices of these genes in available_genes
                    indices = [available_genes.index(g) for g in ct_genes_in_data]
                    if indices:
                        var_group_positions.append((min(indices), max(indices)))
                        var_group_labels.append(cell_type)

        print(f"    Gene groupings: {len(var_group_positions)} cell types with markers")

        # Create the dotplot using scanpy's dotplot function (simpler, better layout handling)
        # Create figure with appropriate size
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        # Use scanpy's dotplot function which handles layout better
        sc.pl.dotplot(
            adata,
            var_names=available_genes,
            groupby=cell_type_column,
            standard_scale='var',
            dendrogram=True,
            var_group_positions=var_group_positions if len(var_group_positions) > 0 else None,
            var_group_labels=var_group_labels if len(var_group_labels) > 0 else None,
            var_group_rotation=45,
            dot_max=dot_max,
            dot_min=dot_min,
            smallest_dot=15,
            cmap='RdBu_r',
            ax=ax,
            show=False
        )

        # Get the actual figure object
        fig = plt.gcf()

        if fig is not None:
            title_fontsize = max(18, gene_fontsize + 6)
            # Add title
            fig.suptitle(f'Cell Type Markers (L2) - Top {top_n_genes} per cell type',
                        fontsize=title_fontsize, fontweight='bold')

            # Create dedicated folder for this top_n value
            top_n_output_dir = os.path.join(output_dir, f'top_{top_n_genes}_markers')
            os.makedirs(top_n_output_dir, exist_ok=True)

            # Let scanpy handle the layout automatically
            plt.tight_layout()

            # Save the plot as PNG
            output_path = os.path.join(top_n_output_dir, f'{cell_type_column}_markers_dotplot.png')
            fig.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.5)
            print(f"\n  ✓ Saved dotplot: {output_path}")

            # Save as PDF for vector graphics
            output_path_pdf = os.path.join(top_n_output_dir, f'{cell_type_column}_markers_dotplot.pdf')
            fig.savefig(output_path_pdf, bbox_inches='tight')
            print(f"  ✓ Saved dotplot (PDF): {output_path_pdf}")

            plt.close(fig)
        else:
            print(f"  ❌ Error: Figure object is None")

    except Exception as e:
        print(f"  ❌ Error creating dotplot: {e}")
        import traceback
        traceback.print_exc()
        plt.close()


def main():
    """Main execution function"""

    print("="*80)
    print("CELL TYPE L2 NEW MARKER GENE DOTPLOT GENERATION")
    print("="*80)

    # Load annotation data
    print(f"\n1. Loading annotation data from:")
    print(f"   {ANNOTATION_H5AD_PATH}")

    if not os.path.exists(ANNOTATION_H5AD_PATH):
        print(f"❌ Error: File not found: {ANNOTATION_H5AD_PATH}")
        return False

    adata = sc.read_h5ad(ANNOTATION_H5AD_PATH)

    print(f"   ✓ Loaded data shape: {adata.shape}")
    print(f"   ✓ Observations (cells): {adata.n_obs}")
    print(f"   ✓ Variables (genes): {adata.n_vars}")

    # Check if cell_type_L2_new exists
    if CELL_TYPE_COLUMN not in adata.obs.columns:
        print(f"❌ Error: Column '{CELL_TYPE_COLUMN}' not found in adata.obs")
        print(f"   Available columns: {list(adata.obs.columns)}")
        return False

    # Get unique cell types
    cell_types = sorted(adata.obs[CELL_TYPE_COLUMN].unique())
    print(f"\n2. Cell types in {CELL_TYPE_COLUMN}: {len(cell_types)}")
    for ct in cell_types:
        n_cells = (adata.obs[CELL_TYPE_COLUMN] == ct).sum()
        print(f"   - {ct}: {n_cells} cells")

    # Create dotplots with different numbers of top markers
    print(f"\n3. Creating dotplots with different numbers of top markers:")
    print(f"   Top N values: {TOP_N_GENES_LIST}")
    print(f"   Biomarker base directory: {BIOMARKER_BASE_DIR}")

    for top_n in TOP_N_GENES_LIST:
        print(f"\n{'='*80}")
        print(f"PROCESSING TOP {top_n} MARKERS")
        print(f"{'='*80}")

        # Load marker genes for each cell type
        marker_genes_dict = {}

        for cell_type in cell_types:
            markers = load_cell_type_markers(BIOMARKER_BASE_DIR, cell_type, top_n=top_n)
            if markers:
                marker_genes_dict[cell_type] = markers

        print(f"\n   ✓ Loaded markers for {len(marker_genes_dict)}/{len(cell_types)} cell types")

        if len(marker_genes_dict) == 0:
            print(f"   ⚠️  Warning: No marker genes loaded for top {top_n}!")
            continue

        # Create dotplot
        print(f"\n   Creating dotplot with top {top_n} markers...")
        create_cell_type_dotplot(
            adata=adata,
            marker_genes_dict=marker_genes_dict,
            output_dir=OUTPUT_DIR,
            top_n_genes=top_n,
            cell_type_column=CELL_TYPE_COLUMN
        )

    print("\n" + "="*80)
    print("COMPLETED")
    print("="*80)
    print(f"Output directory: {OUTPUT_DIR}")

    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
