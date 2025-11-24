# %% [markdown]
# # 2. Subcluster Leiden Cluster 1
#
# This script loads the full dataset, subsets it to cells in Leiden cluster 1,
# and then re-runs dimensionality reduction and clustering to find sub-populations.

# %%
import scanpy as sc
import os
import sys
import anndata as ad
from pathlib import Path

# %%
# --- Configuration ---
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = os.path.join(PROJECT_DIR, "combine_data")
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

# Define Input Directory for Processed Data
INPUT_DIR = os.path.join(WORKING_DIR, "results_from_raw")
input_file = os.path.join(INPUT_DIR, 'merged_raw_processed.h5ad')

# Define Output Directory for this analysis
OUTPUT_DIR = os.path.join(WORKING_DIR, "results_subclustering")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define Cluster to Subcluster
CLUSTER_ID = '1'
LEIDEN_RESOLUTION_KEY = 'leiden_0.4' # The key from the initial clustering

# Define Subclustering Parameters
N_PCS = 30 # Number of PCs to use for subclustering
SUBCLUSTER_RESOLUTION = 0.2 # Resolution for sub-clustering

output_file_subcluster = os.path.join(OUTPUT_DIR, f'subcluster_{CLUSTER_ID}_processed.h5ad')

# %% [markdown]
# ## Load Full Dataset

# %%
print(f"Loading full dataset from {input_file}...")
if not os.path.exists(input_file):
    print(f"Error: Input file not found at {input_file}")
    sys.exit(1)

adata_full = ad.read_h5ad(input_file)
print("Full dataset loaded:")
print(adata_full)

# %% [markdown]
# ## Subset to Cluster

# %%
print(f"Subsetting data to cluster '{CLUSTER_ID}' from '{LEIDEN_RESOLUTION_KEY}'...")

# Ensure the leiden column is treated as categorical
if LEIDEN_RESOLUTION_KEY in adata_full.obs:
    adata_full.obs[LEIDEN_RESOLUTION_KEY] = adata_full.obs[LEIDEN_RESOLUTION_KEY].astype('category')
else:
    print(f"Error: Leiden key '{LEIDEN_RESOLUTION_KEY}' not found in adata.obs")
    sys.exit(1)

# Subset the anndata object
adata_sub = adata_full[adata_full.obs[LEIDEN_RESOLUTION_KEY] == CLUSTER_ID].copy()

print(f"Subset created with {adata_sub.n_obs} cells.")

# %% [markdown]
# ## Re-run Analysis on Subset

# %%
print("\nRe-running analysis on the subset...")

# It's often a good practice to re-identify highly variable genes on the subset
print("Finding highly variable genes on the subset...")
# Using 'cell_ranger' flavor as 'seurat_v3' can cause errors on subsets.
sc.pp.highly_variable_genes(adata_sub, n_top_genes=3000, flavor='cell_ranger')

# Re-run PCA
print(f"Running PCA using {adata_sub.n_vars} genes...")
sc.tl.pca(adata_sub, n_comps=50, use_highly_variable=True, svd_solver='arpack')

# Re-compute neighborhood graph
print(f"Computing neighborhood graph using {N_PCS} PCs...")
sc.pp.neighbors(adata_sub, n_neighbors=15, n_pcs=N_PCS)

# Re-run UMAP
print("Running UMAP...")
sc.tl.umap(adata_sub)

# Re-run Leiden for subclustering
print(f"Running Leiden clustering with resolution {SUBCLUSTER_RESOLUTION}...")
sc.tl.leiden(adata_sub, resolution=SUBCLUSTER_RESOLUTION, key_added=f'leiden_sub_{CLUSTER_ID}')

print("Sub-analysis complete.")

# %% [markdown]
# ## Visualization

# %%
print("\nGenerating visualizations for subclusters...")

# Create directory for plots if it doesn't exist
plot_dir = os.path.join(OUTPUT_DIR, "plots_subclustering")
os.makedirs(plot_dir, exist_ok=True)
sc.settings.figdir = plot_dir

# Plot UMAP colored by the new subclusters
sc.pl.umap(adata_sub, color=f'leiden_sub_{CLUSTER_ID}', legend_loc='on data',
           save=f"_umap_subcluster_{CLUSTER_ID}.png", show=True,
           title=f'Subclusters of Leiden {CLUSTER_ID}')

# Plot UMAP colored by original metadata to see if there's a pattern
for col in ['sample', 'condition', 'genotype']:
    if col in adata_sub.obs:
        sc.pl.umap(adata_sub, color=col,
                   save=f"_umap_subcluster_{CLUSTER_ID}_{col}.png", show=True,
                   title=f'Subcluster {CLUSTER_ID} - colored by {col}')

# %% [markdown]
# ## Save Subclustered Data

# %%
print(f"\nSaving subclustered dataset to {output_file_subcluster}")
try:
    adata_sub.write(Path(output_file_subcluster))
    print("Successfully saved the subclustered AnnData object.")
except Exception as e:
    print(f"Error saving AnnData object: {e}")

print("\nScript finished!")