#!/usr/bin/env python
# Analysis script for Rosa26_SBP1 transgene expression in scRNA-seq data
# This script contains interactive code blocks that can be run with VS Code, PyCharm, or Jupyter

# %% Import libraries
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set plotting style
sc.settings.set_figure_params(dpi=100, frameon=False)
plt.rcParams['figure.figsize'] = (8, 6)

MIN_GENES_PER_CELL = 500
MIN_CELLS_PER_GENE = 3
MAX_MITO_PERCENT = 20

# SELECTED_SAMPLE = "Emx1_Mut"

# %% Load the Cell Ranger count data

adatas = {}

samples = ["Emx1_Ctrl", "Emx1_Mut"]
for sample in samples:
    print(f"Processing {sample}...")
    data_dir = os.path.join("counts_trans", f"cellranger_counts_R26_{sample}_adult", "outs", "filtered_feature_bc_matrix")
    adata = sc.read_10x_mtx(
        data_dir,
        var_names='gene_symbols',
        cache=False
    )
    adatas[sample] = adata


# %% Basic information about the dataset
for sample, adata in adatas.items():
    print(f"Dataset shape: {adata.shape}")
    print(f"Number of cells: {adata.n_obs}")
    print(f"Number of genes: {adata.n_vars}")
    OUTPUT_ADATA = os.path.join("counts_trans", f"{sample}.h5ad")
    print(f"Output AnnData path: {OUTPUT_ADATA}")
    
    # adata = adatas[SELECTED_SAMPLE]
    transgene_name = "Rosa26_SBP1"
    
    adata.obs["sample"] = sample

    print("\nPerforming QC and Filtering on adata...")

    # Calculate mitochondrial gene percentage
    adata.var_names = adata.var_names.astype(str)
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # --- Filtering ---
    print(f"Initial cell count: {adata.n_obs}")
    print(f"Initial gene count: {adata.n_vars}")

    # Calculate the 95th percentile for max_genes
    # Handle case where n_genes_by_counts might be all zeros or NaN
    if 'n_genes_by_counts' in adata.obs and adata.obs['n_genes_by_counts'].notna().any():
        max_genes_percentile = np.percentile(adata.obs['n_genes_by_counts'].dropna(), 95)
        print(f"Calculated 95th percentile for max genes: {max_genes_percentile:.0f}")
    else:
        print("Warning: 'n_genes_by_counts' not found or contains NaNs. Skipping max genes filter.")
        max_genes_percentile = np.inf # Set to infinity to effectively disable the filter

    # --- Apply Filters ---
    # 1. Filter cells based on min number of genes
    n_obs_before = adata.n_obs
    if 'n_genes_by_counts' in adata.obs:
        adata = adata[adata.obs.n_genes_by_counts > MIN_GENES_PER_CELL, :]
        print(f"Cells after min_genes filter (> {MIN_GENES_PER_CELL}): {adata.n_obs} (removed {n_obs_before - adata.n_obs})")
    else:
        print("Skipping min_genes filter: 'n_genes_by_counts' not found.")

    # 2. Filter cells based on max number of genes (percentile)
    n_obs_before = adata.n_obs
    if 'n_genes_by_counts' in adata.obs and np.isfinite(max_genes_percentile):
        adata = adata[adata.obs.n_genes_by_counts < max_genes_percentile, :]
        print(f"Cells after max_genes filter (< {max_genes_percentile:.0f}): {adata.n_obs} (removed {n_obs_before - adata.n_obs})")
    else:
        print("Skipping max_genes filter.")

    # 3. Filter cells based on mitochondrial content
    n_obs_before = adata.n_obs
    if 'pct_counts_mt' in adata.obs:
        adata = adata[adata.obs.pct_counts_mt < MAX_MITO_PERCENT, :]
        print(f"Cells after mito filter (< {MAX_MITO_PERCENT}%): {adata.n_obs} (removed {n_obs_before - adata.n_obs})")
    else:
        print("Skipping mito filter: 'pct_counts_mt' not found.")

    # 4. Filter genes based on minimum number of cells expressing them
    n_vars_before = adata.n_vars
    if adata.n_obs > 0:
        sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
        print(f"Genes after min_cells filter ({MIN_CELLS_PER_GENE}): {adata.n_vars} (removed {n_vars_before - adata.n_vars})")
    else:
        print("Skipping gene filtering as no cells remain after QC.")

    print(f"Final shape of filtered adata: {adata.shape}")
    print("QC and Filtering complete.")


    # Analyze transgene expression

    # Extract transgene expression values
    transgene_expr = adata[:, transgene_name].X.toarray().flatten()

    # Basic statistics
    positive_cells = np.sum(transgene_expr > 0)
    percent_positive = (positive_cells / adata.n_obs) * 100

    print(f"\nTransgene: {transgene_name}")
    print(f"Cells expressing transgene: {positive_cells} out of {adata.n_obs} ({percent_positive:.2f}%)")
    print(f"Mean expression (all cells): {np.mean(transgene_expr):.4f}")
    print(f"Mean expression (positive cells): {np.mean(transgene_expr[transgene_expr > 0]):.4f}")
    print(f"Max expression: {np.max(transgene_expr):.4f}")

    # Visualize transgene expression
    plt.figure(figsize=(10, 6))

    # Distribution of expression values (non-zero)
    plt.subplot(1, 2, 1)
    positive_expr = transgene_expr[transgene_expr > 0]
    if len(positive_expr) > 0:
        sns.histplot(positive_expr, kde=True)
        plt.title(f"{transgene_name} Expression\n(Positive Cells Only)")
        plt.xlabel("Expression Level")
        plt.ylabel("Number of Cells")
    else:
        plt.text(0.5, 0.5, "No positive cells", ha='center', va='center')

    # Pie chart of positive vs negative cells
    plt.subplot(1, 2, 2)
    plt.pie([positive_cells, adata.n_obs - positive_cells], 
            labels=[f'Positive\n{percent_positive:.1f}%', f'Negative\n{100-percent_positive:.1f}%'],
            autopct='%1.1f%%', colors=['#ff9999','#66b3ff'])
    plt.title(f"{transgene_name} Expression")

    plt.tight_layout()
    plt.show()

    # Process data for clustering and visualization
    # The data will be saved with raw counts.
    # Normalization and log transformation are removed from this step.
    # If needed for specific visualizations within this script,
    # consider applying them to a temporary copy of adata.

    # Add transgene expression as observation
    adata.obs['transgene_expr'] = transgene_expr
    adata.obs['transgene_positive'] = adata.obs['transgene_expr'] > 0


    # Save adata
    print(f"\nSaving updated AnnData object to {OUTPUT_ADATA}...")
    adata.write_h5ad(OUTPUT_ADATA)
    print("Updated AnnData object saved successfully.")


    # create additional normalized version
    # Pre-process the data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Save adata
    OUTPUT_ADATA = os.path.join("counts_trans", f"{sample}_normalized.h5ad")
    print(f"\nSaving updated AnnData object to {OUTPUT_ADATA}...")
    adata.write_h5ad(OUTPUT_ADATA)
    print("Updated AnnData object saved successfully.")