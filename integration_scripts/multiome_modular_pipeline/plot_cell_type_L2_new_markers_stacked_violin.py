#!/usr/bin/env python3
"""
Script to generate marker gene expression stacked violin plots for cell_type_L2_new annotations
Creates violin plots showing expression distributions across cell types

This script:
1. Loads annotation_final.h5ad
2. Loads biomarker genes from DEGs_cell_type_L2_newFC_0_25/biomarkers/sig_deg_lists
3. Creates stacked violin plots with top markers per cell type
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
OUTPUT_DIR = '/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/ADDITIONAL_PLOTS/stacked_violins/L2'
CELL_TYPE_COLUMN = 'cell_type_L2_new'
# Generate multiple plots with different numbers of top markers
TOP_N_GENES_LIST = [2]  # Will create separate plots for each

def load_cell_type_markers(biomarker_dir, cell_type, top_n=5):
    """Load top marker genes for a specific cell type"""
    cell_type_folder_map = {
        'Lamp5 (DG)': 'Lamp5_DG',
        'VIP (DG)': 'VIP_DG',
        'ENT CTX': 'ENT_CTX',
        'Mature GC': 'Mature_GC',
        'Immature GC': 'Immature_GC',
        'Mossy cells': 'Mossy_cells'
    }

    folder_name = cell_type_folder_map.get(cell_type, cell_type)
    cell_type_dir = Path(biomarker_dir) / folder_name
    csv_file = cell_type_dir / f'Cluster_{cell_type}_vs_Rest_up_significant.csv'

    if not csv_file.exists():
        print(f"  ⚠️  Warning: Could not find {csv_file}")
        return []

    try:
        df = pd.read_csv(csv_file)
        if 'names' in df.columns:
            top_genes = df.head(top_n)['names'].tolist()
        else:
            top_genes = df.iloc[:top_n, 0].tolist()

        print(f"  ✓ Loaded {len(top_genes)} markers for {cell_type}")
        return top_genes
    except Exception as e:
        print(f"  ❌ Error loading markers for {cell_type}: {e}")
        return []


def create_stacked_violin(adata, marker_genes_dict, output_dir, top_n_genes, cell_type_column='cell_type_L2_new'):
    """Create stacked violin plot showing expression distributions"""

    print(f"\n{'='*80}")
    print(f"Creating stacked violin plot for {cell_type_column} (top {top_n_genes} genes per cell type)")
    print(f"{'='*80}")

    cell_types_in_data = sorted(adata.obs[cell_type_column].unique())

    # Collect all marker genes
    all_marker_genes = []
    for cell_type in cell_types_in_data:
        if cell_type in marker_genes_dict:
            all_marker_genes.extend(marker_genes_dict[cell_type])

    # Remove duplicates while preserving order
    unique_marker_genes = []
    seen = set()
    for gene in all_marker_genes:
        if gene not in seen:
            unique_marker_genes.append(gene)
            seen.add(gene)

    available_genes = [g for g in unique_marker_genes if g in adata.var_names]
    print(f"  Total unique genes: {len(unique_marker_genes)}, Available: {len(available_genes)}")

    if len(available_genes) == 0:
        print(f"  ❌ Error: No marker genes found!")
        return

    # Figure sizing - stacked violins need more vertical space
    n_cell_types = len(cell_types_in_data)
    n_genes = len(available_genes)

    # Two orientations: standard (genes as columns) and swapped (genes as rows)
    for swap_axes in [False, True]:
        try:
            if swap_axes:
                # Genes on Y-axis, cell types on X-axis
                fig_width = max(14, n_cell_types * 0.8)
                fig_height = max(20, n_genes * 0.3)
                fontsize = 11
                suffix = "_swapped"
            else:
                # Genes on X-axis, cell types on Y-axis
                fig_width = max(20, n_genes * 0.40)
                fig_height = max(14, n_cell_types * 0.8)
                fontsize = 11
                suffix = ""

            plt.rcParams['font.size'] = fontsize
            plt.rcParams['axes.labelsize'] = fontsize + 2
            plt.rcParams['axes.titlesize'] = fontsize + 4

            os.makedirs(output_dir, exist_ok=True)

            # Prepare gene groupings
            var_group_positions = []
            var_group_labels = []

            for cell_type in cell_types_in_data:
                if cell_type in marker_genes_dict:
                    genes = marker_genes_dict[cell_type]
                    ct_genes_in_data = [g for g in genes if g in available_genes]
                    if len(ct_genes_in_data) > 0:
                        indices = [available_genes.index(g) for g in ct_genes_in_data]
                        if indices:
                            var_group_positions.append((min(indices), max(indices)))
                            var_group_labels.append(cell_type)

            # Create stacked violin plot using scanpy's convenience function
            # Store current figure count to identify the new figure
            figs_before = set(plt.get_fignums())

            sc.pl.stacked_violin(
                adata,
                var_names=available_genes,
                groupby=cell_type_column,
                dendrogram=False,  # Disable dendrogram to avoid matplotlib compatibility issues
                var_group_positions=var_group_positions if len(var_group_positions) > 0 else None,
                var_group_labels=var_group_labels if len(var_group_labels) > 0 else None,
                var_group_rotation=45 if not swap_axes else 0,
                swap_axes=swap_axes,
                figsize=(fig_width, fig_height),
                show=False
            )

            # Get the newly created figure
            figs_after = set(plt.get_fignums())
            new_figs = figs_after - figs_before
            if new_figs:
                fig = plt.figure(list(new_figs)[0])
            else:
                fig = plt.gcf()

            if fig is not None and hasattr(fig, 'get_axes') and len(fig.get_axes()) > 0:
                orientation_text = "(Swapped Axes)" if swap_axes else "(Standard)"
                fig.suptitle(f'Cell Type Markers (L2) - Stacked Violin {orientation_text}\nTop {top_n_genes} per cell type',
                            fontsize=fontsize + 6, fontweight='bold')

                # Save outputs
                top_n_output_dir = os.path.join(output_dir, f'top_{top_n_genes}_markers')
                os.makedirs(top_n_output_dir, exist_ok=True)

                plt.tight_layout()

                output_path = os.path.join(top_n_output_dir, f'{cell_type_column}_markers_stacked_violin{suffix}.png')
                fig.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.5)
                print(f"\n  ✓ Saved stacked violin{suffix}: {output_path}")

                output_path_pdf = os.path.join(top_n_output_dir, f'{cell_type_column}_markers_stacked_violin{suffix}.pdf')
                fig.savefig(output_path_pdf, bbox_inches='tight')
                print(f"  ✓ Saved stacked violin{suffix} (PDF): {output_path_pdf}")

                plt.close(fig)
            else:
                print(f"  ❌ Error: Figure object is None")

        except Exception as e:
            print(f"  ❌ Error creating stacked violin (swap={swap_axes}): {e}")
            import traceback
            traceback.print_exc()
            plt.close()


def main():
    """Main execution function"""

    print("="*80)
    print("CELL TYPE L2 MARKER GENE STACKED VIOLIN PLOT GENERATION")
    print("="*80)

    print(f"\n1. Loading annotation data from:")
    print(f"   {ANNOTATION_H5AD_PATH}")

    if not os.path.exists(ANNOTATION_H5AD_PATH):
        print(f"❌ Error: File not found: {ANNOTATION_H5AD_PATH}")
        return False

    adata = sc.read_h5ad(ANNOTATION_H5AD_PATH)
    print(f"   ✓ Loaded data: {adata.shape}")

    if CELL_TYPE_COLUMN not in adata.obs.columns:
        print(f"❌ Error: Column '{CELL_TYPE_COLUMN}' not found")
        return False

    cell_types = sorted(adata.obs[CELL_TYPE_COLUMN].unique())
    print(f"\n2. Cell types: {len(cell_types)}")

    for top_n in TOP_N_GENES_LIST:
        print(f"\n{'='*80}")
        print(f"PROCESSING TOP {top_n} MARKERS")
        print(f"{'='*80}")

        marker_genes_dict = {}
        for cell_type in cell_types:
            markers = load_cell_type_markers(BIOMARKER_BASE_DIR, cell_type, top_n=top_n)
            if markers:
                marker_genes_dict[cell_type] = markers

        print(f"\n   ✓ Loaded markers for {len(marker_genes_dict)}/{len(cell_types)} cell types")

        if len(marker_genes_dict) == 0:
            print(f"   ⚠️  Warning: No marker genes loaded!")
            continue

        create_stacked_violin(
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
