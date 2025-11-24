#!/usr/bin/env python3
"""
Create Collages from Gene Set Enrichment Heatmaps

This script creates organized collages from individual gene set heatmaps,
grouping them by:
- Score type (positive, negative, aggregate)
- Cell type level (L1, L2)
- Split mode (condition, genotype_condition)

Input:
    - Individual heatmap PNG files from gene_set_enrichment_heatmaps_custom/

Output:
    - Multi-panel collage figures showing all gene sets together
    - Organized by score type and cell type for easy comparison
"""

import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import numpy as np
from PIL import Image
import re

# --- Configuration ---
WORKING_DIR = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data"
os.chdir(WORKING_DIR)

# Input directory containing individual heatmaps
INPUT_DIR = os.path.join(WORKING_DIR, "results_from_raw/gene_set_enrichment_heatmaps_custom")

# Output directory for collages
OUTPUT_DIR = os.path.join(WORKING_DIR, "results_from_raw/gene_set_enrichment_collages")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Configuration for collages
SCORE_TYPES = ['positive', 'negative', 'aggregate']
CELL_TYPE_LEVELS = ['cell_type_L1', 'cell_type_L2_new']
SPLIT_MODES = ['condition', 'genotype_condition']

# Publication-ready plot settings
plt.rcParams['font.size'] = 8
plt.rcParams['axes.labelsize'] = 9
plt.rcParams['axes.titlesize'] = 11
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 12
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def find_heatmap_files(input_dir, score_type, cell_type_level, split_mode):
    """
    Find all heatmap PNG files matching the specified criteria.

    Parameters:
    -----------
    input_dir : str
        Directory containing the heatmap files
    score_type : str
        Score type (positive, negative, aggregate)
    cell_type_level : str
        Cell type level (cell_type_L1, cell_type_L2_new)
    split_mode : str
        Split mode (condition, genotype_condition)

    Returns:
    --------
    list
        List of tuples: (gene_set_name, file_path)
    """
    # Pattern to match files
    pattern = f"*_{score_type}_{cell_type_level}_{split_mode}.png"

    # Find all matching files
    matching_files = []
    for file_path in Path(input_dir).glob(pattern):
        # Extract gene set name from filename
        filename = file_path.stem  # Remove .png extension

        # Remove the suffix to get gene set name
        suffix = f"_{score_type}_{cell_type_level}_{split_mode}"
        gene_set_name = filename.replace(suffix, '')

        # Clean up gene set name for display
        gene_set_display = gene_set_name.replace('_', ' ')

        matching_files.append((gene_set_display, str(file_path)))

    # Sort by gene set name
    matching_files.sort(key=lambda x: x[0])

    return matching_files


def create_collage(file_list, score_type, cell_type_level, split_mode,
                   output_dir, max_cols=3):
    """
    Create a collage from a list of heatmap files.

    Parameters:
    -----------
    file_list : list
        List of tuples: (gene_set_name, file_path)
    score_type : str
        Score type for labeling
    cell_type_level : str
        Cell type level for labeling
    split_mode : str
        Split mode for labeling
    output_dir : str
        Directory to save the collage
    max_cols : int
        Maximum number of columns in the collage
    """
    if not file_list:
        print(f"  No files found for {score_type}/{cell_type_level}/{split_mode}")
        return

    n_plots = len(file_list)
    print(f"  Creating collage with {n_plots} heatmaps...")

    # Calculate grid dimensions
    n_cols = min(max_cols, n_plots)
    n_rows = int(np.ceil(n_plots / n_cols))

    # Load images to determine aspect ratios
    images = []
    titles = []
    max_height = 0
    max_width = 0

    for gene_set_name, file_path in file_list:
        img = Image.open(file_path)
        images.append(img)
        titles.append(gene_set_name)
        max_height = max(max_height, img.height)
        max_width = max(max_width, img.width)

    # Calculate figure size
    # Base size on the number of plots and their dimensions
    subplot_width = 6  # inches per subplot
    subplot_height = subplot_width * (max_height / max_width) if max_width > 0 else 6

    fig_width = n_cols * subplot_width
    fig_height = n_rows * subplot_height + 1  # +1 for main title

    # Create figure
    fig = plt.figure(figsize=(fig_width, fig_height))

    # Add main title
    score_type_display = score_type.capitalize()
    cell_type_display = cell_type_level.replace('_', ' ').replace('cell type', 'Cell Type')
    split_mode_display = split_mode.replace('_', ' ').replace('condition', 'Condition').replace('genotype', 'Genotype')

    main_title = f'Gene Set Enrichment Heatmaps\n{score_type_display} Scores - {cell_type_display} - {split_mode_display}'
    fig.suptitle(main_title, fontsize=14, fontweight='bold', y=0.995)

    # Create grid
    gs = GridSpec(n_rows, n_cols, figure=fig,
                  left=0.05, right=0.95, top=0.96, bottom=0.02,
                  hspace=0.3, wspace=0.2)

    # Plot each heatmap
    for idx, (gene_set_name, file_path) in enumerate(file_list):
        row = idx // n_cols
        col = idx % n_cols

        ax = fig.add_subplot(gs[row, col])

        # Load and display image
        img = mpimg.imread(file_path)
        ax.imshow(img)
        ax.axis('off')

        # Add subplot title
        ax.set_title(gene_set_name, fontsize=10, fontweight='bold', pad=5)

    # Hide any empty subplots
    for idx in range(n_plots, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        ax = fig.add_subplot(gs[row, col])
        ax.axis('off')

    # Save the collage
    base_filename = f'collage_{score_type}_{cell_type_level}_{split_mode}'

    plt.savefig(os.path.join(output_dir, f'{base_filename}.png'),
               bbox_inches='tight', dpi=300, facecolor='white')
    plt.savefig(os.path.join(output_dir, f'{base_filename}.pdf'),
               bbox_inches='tight', dpi=300, facecolor='white')
    plt.close()

    print(f"  Saved: {base_filename}")


def create_combined_score_collage(input_dir, gene_set_names, cell_type_level,
                                   split_mode, output_dir):
    """
    Create a collage showing all three score types (positive, negative, aggregate)
    for each gene set side-by-side.

    Parameters:
    -----------
    input_dir : str
        Directory containing the heatmap files
    gene_set_names : list
        List of gene set names
    cell_type_level : str
        Cell type level
    split_mode : str
        Split mode
    output_dir : str
        Directory to save the collage
    """
    print(f"\n  Creating combined score type collage for {cell_type_level}/{split_mode}...")

    # Collect files for each gene set across all score types
    gene_set_files = {}

    for gene_set_name in gene_set_names:
        files = {}
        for score_type in SCORE_TYPES:
            # Find the file for this combination
            pattern = f"*_{score_type}_{cell_type_level}_{split_mode}.png"
            matching = list(Path(input_dir).glob(pattern))

            # Find the one matching this gene set
            for file_path in matching:
                filename = file_path.stem
                suffix = f"_{score_type}_{cell_type_level}_{split_mode}"
                extracted_name = filename.replace(suffix, '').replace('_', ' ')

                if extracted_name.lower() == gene_set_name.lower():
                    files[score_type] = str(file_path)
                    break

        if len(files) == 3:  # Only include if all three score types are present
            gene_set_files[gene_set_name] = files

    if not gene_set_files:
        print(f"    No complete gene sets found")
        return

    n_gene_sets = len(gene_set_files)
    print(f"    Found {n_gene_sets} gene sets with all three score types")

    # Create figure with 3 columns (one per score type) and n_gene_sets rows
    n_cols = 3
    n_rows = n_gene_sets

    # Calculate figure size
    subplot_width = 5
    subplot_height = 4

    fig_width = n_cols * subplot_width
    fig_height = n_rows * subplot_height + 1

    fig = plt.figure(figsize=(fig_width, fig_height))

    # Add main title
    cell_type_display = cell_type_level.replace('_', ' ').replace('cell type', 'Cell Type')
    split_mode_display = split_mode.replace('_', ' ').replace('condition', 'Condition').replace('genotype', 'Genotype')

    main_title = f'Gene Set Enrichment: All Score Types\n{cell_type_display} - {split_mode_display}'
    fig.suptitle(main_title, fontsize=14, fontweight='bold', y=0.995)

    # Create grid
    gs = GridSpec(n_rows, n_cols, figure=fig,
                  left=0.05, right=0.95, top=0.96, bottom=0.02,
                  hspace=0.35, wspace=0.15)

    # Add column headers
    for col_idx, score_type in enumerate(SCORE_TYPES):
        ax = fig.add_subplot(gs[0, col_idx])
        ax.text(0.5, 1.15, score_type.capitalize(),
               ha='center', va='bottom', transform=ax.transAxes,
               fontsize=12, fontweight='bold')

    # Plot heatmaps
    for row_idx, (gene_set_name, files) in enumerate(sorted(gene_set_files.items())):
        for col_idx, score_type in enumerate(SCORE_TYPES):
            ax = fig.add_subplot(gs[row_idx, col_idx])

            if score_type in files:
                img = mpimg.imread(files[score_type])
                ax.imshow(img)

            ax.axis('off')

            # Add row label (gene set name) on the left
            if col_idx == 0:
                ax.text(-0.05, 0.5, gene_set_name,
                       ha='right', va='center', rotation=0,
                       transform=ax.transAxes,
                       fontsize=9, fontweight='bold')

    # Save the collage
    base_filename = f'collage_all_scores_{cell_type_level}_{split_mode}'

    plt.savefig(os.path.join(output_dir, f'{base_filename}.png'),
               bbox_inches='tight', dpi=300, facecolor='white')
    plt.savefig(os.path.join(output_dir, f'{base_filename}.pdf'),
               bbox_inches='tight', dpi=300, facecolor='white')
    plt.close()

    print(f"    Saved: {base_filename}")


def main():
    """Main function to create all collages."""

    print("=" * 70)
    print("Creating Gene Set Enrichment Heatmap Collages")
    print("=" * 70)

    if not os.path.exists(INPUT_DIR):
        print(f"Error: Input directory not found: {INPUT_DIR}")
        sys.exit(1)

    # Check if there are any PNG files
    png_files = list(Path(INPUT_DIR).glob("*.png"))
    if not png_files:
        print(f"Error: No PNG files found in {INPUT_DIR}")
        sys.exit(1)

    print(f"Found {len(png_files)} heatmap files")

    # Create collages for each combination of parameters
    print("\n" + "=" * 70)
    print("Creating Individual Score Type Collages")
    print("=" * 70)

    for score_type in SCORE_TYPES:
        for cell_type_level in CELL_TYPE_LEVELS:
            for split_mode in SPLIT_MODES:
                print(f"\nProcessing: {score_type} / {cell_type_level} / {split_mode}")

                # Find matching files
                file_list = find_heatmap_files(INPUT_DIR, score_type,
                                              cell_type_level, split_mode)

                # Create collage
                create_collage(file_list, score_type, cell_type_level,
                             split_mode, OUTPUT_DIR, max_cols=3)

    # Create combined score type collages
    print("\n" + "=" * 70)
    print("Creating Combined Score Type Collages")
    print("=" * 70)

    # Get unique gene set names
    all_gene_sets = set()
    for file_path in Path(INPUT_DIR).glob("*_positive_*.png"):
        filename = file_path.stem
        # Extract gene set name
        for score_type in SCORE_TYPES:
            for cell_type_level in CELL_TYPE_LEVELS:
                for split_mode in SPLIT_MODES:
                    suffix = f"_{score_type}_{cell_type_level}_{split_mode}"
                    if suffix in filename:
                        gene_set = filename.replace(suffix, '').replace('_', ' ')
                        all_gene_sets.add(gene_set)

    for cell_type_level in CELL_TYPE_LEVELS:
        for split_mode in SPLIT_MODES:
            create_combined_score_collage(INPUT_DIR, list(all_gene_sets),
                                         cell_type_level, split_mode, OUTPUT_DIR)

    print("\n" + "=" * 70)
    print("Collage creation completed successfully!")
    print(f"Results saved to: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
