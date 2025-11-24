# %%
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import seaborn as sns
from scipy.sparse import csr_matrix
import random

# PROJECT_DIR = "/home/michal/Github/SRF_Linda_RNA"
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

# Set seeds for all random number generators
random_seed = 0
np.random.seed(random_seed)
random.seed(random_seed)

REMOVE_DOUBLETS = False

cluster_key = 'leiden_0.4'

# %%
# Set up directories
# Set up directories
if REMOVE_DOUBLETS:
    BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw", "doublets_removed")
else:
    BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")

INPUT_DIR = BASE_RESULTS_DIR
adata_path = os.path.join(INPUT_DIR, 'annotation_final.h5ad')

PARENT_OUTPUT_DIR = f"{INPUT_DIR}/clusters_similarity"
PLOT_OUTPUT_DIR = os.path.join(PARENT_OUTPUT_DIR)

# Create output directories if they don't exist
os.makedirs(PLOT_OUTPUT_DIR, exist_ok=True)

print(f"Input directory: {INPUT_DIR}")
print(f"Plot output directory: {PLOT_OUTPUT_DIR}")

# Configure scanpy settings
sc.settings.figdir = PLOT_OUTPUT_DIR
sc.settings.set_figure_params(dpi=150, facecolor='white') # Adjust DPI if needed

# %%
# Load merged data
print(f"\nLoading merged dataset from {adata_path}")
adata = sc.read_h5ad(adata_path)

# %%
# Calculate cluster similarity based on average gene expression
print(f"\nCalculating cluster similarity using '{cluster_key}'...")

# Calculate average expression per cluster
# Ensure adata.X is not sparse or convert it if necessary
# Using layer='log1p' if normalized data is stored there, otherwise use adata.X
# Assuming log1p normalized data is suitable for correlation
if 'log1p' in adata.uns:
    if isinstance(adata.X, csr_matrix):
        expression_matrix = adata.X.toarray() # Use raw counts if log1p not directly in X
    else:
        expression_matrix = adata.X # Use raw counts if log1p not directly in X
    print("Using adata.X for expression matrix (assuming log1p was applied).")
    # If you have log1p data stored in a layer, use:
    # expression_matrix = adata.layers['log1p'].toarray() if isinstance(adata.layers['log1p'], csr_matrix) else adata.layers['log1p']
    # print("Using adata.layers['log1p'] for expression matrix.")
else:
     if isinstance(adata.X, csr_matrix):
        expression_matrix = adata.X.toarray()
     else:
        expression_matrix = adata.X
     print("Using adata.X for expression matrix.")


# Create a DataFrame for easier manipulation
expression_df = pd.DataFrame(expression_matrix, index=adata.obs_names, columns=adata.var_names)   # type: ignore
expression_df[cluster_key] = adata.obs[cluster_key].astype('category').values # Ensure cluster key is categorical

# Calculate mean expression per cluster
cluster_means = expression_df.groupby(cluster_key, observed=True).mean() # Use observed=True for categorical groupby

# Calculate pairwise correlation between clusters (using Pearson correlation)
cluster_corr = cluster_means.T.corr() # Transpose because corr calculates column correlations

print("Cluster similarity calculation complete.")

# %%
# Visualize the similarity matrix as a heatmap
print(f"\nPlotting cluster similarity heatmap for '{cluster_key}'...")
# Increase figure size slightly and adjust annotation font size
sns.clustermap(cluster_corr, cmap='viridis', annot=True, fmt=".2f", linewidths=.5,
               figsize=(16, 16), annot_kws={"size": 6}) # Increased figsize
plt.suptitle(f'Cluster Similarity based on Mean Expression ({cluster_key})', y=1.03) # Adjust title position
heatmap_path = os.path.join(PLOT_OUTPUT_DIR, f'cluster_similarity_heatmap_{cluster_key}.png')
plt.savefig(heatmap_path, bbox_inches='tight')
print(f"Heatmap saved to {heatmap_path}")
# plt.show() # Keep commented out for script execution

# %%
# Calculate cluster similarity based on Highly Variable Genes (HVGs)
# This focuses the comparison on genes driving variability in the dataset.

print(f"\nCalculating cluster similarity using HVGs for '{cluster_key}'...")

if 'highly_variable' not in adata.var.columns:
    print("Warning: 'highly_variable' column not found in adata.var. Skipping HVG analysis.")
else:
    hvg_mask = adata.var['highly_variable']
    if hvg_mask.sum() == 0:
         print("Warning: No highly variable genes found. Skipping HVG analysis.")
    else:
        # Filter expression DataFrame for HVGs
        expression_df_hvg = expression_df[adata.var_names[hvg_mask]]
        expression_df_hvg[cluster_key] = adata.obs[cluster_key].astype('category').values # Add cluster labels back

        # Calculate mean expression per cluster using only HVGs
        cluster_means_hvg = expression_df_hvg.groupby(cluster_key, observed=True).mean()

        # Calculate pairwise correlation between clusters based on HVGs
        cluster_corr_hvg = cluster_means_hvg.T.corr()

        print("HVG-based cluster similarity calculation complete.")

        # Visualize the HVG similarity matrix as a heatmap
        print(f"\nPlotting HVG-based cluster similarity heatmap for '{cluster_key}'...")
        sns.clustermap(cluster_corr_hvg, cmap='viridis', annot=True, fmt=".2f", linewidths=.5,
                       figsize=(16, 16), annot_kws={"size": 6}) # Increased figsize
        plt.suptitle(f'Cluster Similarity based on HVG Mean Expression ({cluster_key})', y=1.03)
        hvg_heatmap_path = os.path.join(PLOT_OUTPUT_DIR, f'cluster_similarity_heatmap_hvg_{cluster_key}.png')
        plt.savefig(hvg_heatmap_path, bbox_inches='tight')
        print(f"HVG heatmap saved to {hvg_heatmap_path}")
        # plt.show()

# %%
# Alternative: Use scanpy's dendrogram function
# This calculates correlation and plots a dendrogram.
# Using use_rep='X_pca' calculates similarity based on the principal components,
# which often captures major biological variation robustly.
print(f"\nGenerating dendrogram for '{cluster_key}' using scanpy...")
sc.tl.dendrogram(adata, groupby=cluster_key, use_rep='X_pca') # Use PCA representation for robustness
dendro_path_key = f'_dendrogram_{cluster_key}.png'

# Temporarily adjust figure params for dendrogram clarity
# Get current settings from matplotlib's rcParams, which scanpy modifies
original_figsize = plt.rcParams['figure.figsize']
original_dpi = plt.rcParams['figure.dpi']
sc.settings.set_figure_params(figsize=(14.0, 6.0), dpi=150) # type: ignore

# Generate and save the plot
sc.pl.dendrogram(adata, groupby=cluster_key, show=False) # Use show=False, save happens automatically via settings

# Add title manually after plotting but before saving is finalized by show=False/implicit save
plt.title(f'Cluster Dendrogram ({cluster_key})', y=1.05) # Add a title

# Construct path and explicitly save with title
expected_dendro_path = os.path.join(PLOT_OUTPUT_DIR, f'dendrogram{dendro_path_key}')
plt.savefig(expected_dendro_path, bbox_inches='tight')
plt.close() # Close the figure to prevent display in interactive environments

# Restore original settings using set_figure_params
sc.settings.set_figure_params(figsize=original_figsize, dpi=original_dpi)
print(f"Dendrogram saved to {expected_dendro_path}")
