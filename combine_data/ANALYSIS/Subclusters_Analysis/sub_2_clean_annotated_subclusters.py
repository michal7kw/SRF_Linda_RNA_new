# %% [markdown]
# # 13. Clean and Simplify Annotated Subclusters
#
# This script loads a fully annotated subcluster AnnData object, filters out
# sex-related genes, and simplifies the object by removing unnecessary data
# slots before saving the final cleaned version.

# %%
import scanpy as sc
import numpy as np
import pandas as pd
import os
import sys
import anndata as ad
import random

# %%
# --- Configuration: SET THE CLUSTER ID TO PROCESS ---
CLUSTER_ID = '0' # <-- CHANGE THIS to '1' or '16' to process other clusters

# --- Script Configuration ---
# Set seeds for reproducibility
random_seed = 0
np.random.seed(random_seed)
random.seed(random_seed)

# Define Project Directories
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = os.path.join(PROJECT_DIR, "combine_data")
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

# Define input/output directories and file paths
BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_subclustering")

input_path = os.path.join(BASE_RESULTS_DIR, f'subcluster_{CLUSTER_ID}_annotated.h5ad')
output_path = os.path.join(BASE_RESULTS_DIR, f'subcluster_{CLUSTER_ID}_cleaned.h5ad')

print(f"--- Cleaning Subcluster {CLUSTER_ID} ---")
print(f"Input file: {input_path}")
print(f"Output file: {output_path}")

# %% [markdown]
# ## Load Annotated Subcluster Data

# %%
print(f"\nLoading annotated dataset from {input_path}")
if not os.path.exists(input_path):
    raise FileNotFoundError(f"Input file not found: {input_path}. Please run the annotation scripts first.")
adata = sc.read_h5ad(input_path)
print("Dataset loaded.")

# %% [markdown]
# ## Filter Sex-Related Genes

# %%
genes_to_remove = ["Xist", "Tsix", "Eif2s3x", "Ddx3x", "Kdm5c", "Kdm6a", "Rps6ka3", "Uty", "Ddx3y", "Eif2s3y", "Kdm5d", "Rps4y1", "Jarid1d"]
print(f"\nFiltering out {len(genes_to_remove)} sex-related genes...")

genes_present = adata.var_names.intersection(genes_to_remove)
print(f"Found {len(genes_present)} sex-related genes in adata: {list(genes_present)}")

if len(genes_present) > 0:
    print(f"Original number of genes: {adata.n_vars}")
    adata = adata[:, ~adata.var_names.isin(genes_present)].copy()
    print(f"Number of genes after filtering: {adata.n_vars}")

    if adata.raw is not None:
        print(f"Original number of genes in raw: {adata.raw.n_vars}")
        raw_genes_present = adata.raw.var_names.intersection(genes_to_remove)
        if len(raw_genes_present) > 0:
            # Create a new AnnData object for the filtered raw data, which is the expected type for the .raw attribute
            adata_raw_filtered = adata.raw[:, ~adata.raw.var_names.isin(raw_genes_present)]
            adata.raw = ad.AnnData(X=adata_raw_filtered.X.copy(), var=adata_raw_filtered.var.copy())
            print(f"Number of genes in raw after filtering: {adata.raw.n_vars}")
else:
    print("No specified sex-related genes found to filter.")

# %% [markdown]
# ## Simplify AnnData Object

# %%
print("\nSimplifying the AnnData object...")

# Define keys to keep, including new subcluster and annotation columns
SUBCLUSTER_KEY = f'leiden_sub_{CLUSTER_ID}'
MANUAL_ANNOTATION_KEY = f'manual_annotation_subcluster_{CLUSTER_ID}'

obs_keep = [
    'sample', 'condition', 'genotype', 'n_genes_by_counts', 'total_counts',
    'pct_counts_mt', 'leiden_0.4', SUBCLUSTER_KEY,
    'DG_majority_voting', 'DG_conf_score', 'ISO_majority_voting', 'ISO_conf_score'
]
# The manual annotation key is not expected at this stage anymore.
# It will be added in the next script.

var_keep = [
    'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts',
    'total_counts', 'n_cells', 'highly_variable', 'highly_variable_rank',
    'means', 'variances', 'variances_norm'
]
uns_keep = [
    'leiden_0.4_colors', f'{SUBCLUSTER_KEY}_colors', 'DG_majority_voting_colors',
    'ISO_majority_voting_colors', 'condition_colors', 'genotype_colors', 'hvg',
    'leiden_0.4', 'neighbors', 'pca', 'sample_colors', 'umap'
]
obsm_keep = ['X_pca', 'X_umap']
varm_keep = ['PCs']
obsp_keep = ['connectivities', 'distances']

# Filter the lists to only include keys that actually exist in the object
obs_keep_existing = [k for k in obs_keep if k in adata.obs.columns]
var_keep_existing = [k for k in var_keep if k in adata.var.columns]
uns_keep_existing = [k for k in uns_keep if k in adata.uns]

# Create a new AnnData object with the selected slots
adata_simple = ad.AnnData(
    X=adata.X.copy(), # type: ignore
    obs=adata.obs[obs_keep_existing].copy(),
    var=adata.var[var_keep_existing].copy()
)

adata_simple.uns = {k: adata.uns[k] for k in uns_keep_existing}
adata_simple.obsm = {k: adata.obsm[k] for k in obsm_keep if k in adata.obsm}
adata_simple.varm = {k: adata.varm[k] for k in varm_keep if k in adata.varm}
adata_simple.obsp = {k: adata.obsp[k] for k in obsp_keep if k in adata.obsp}

if adata.raw is not None:
    print("Copying raw data...")
    # Re-wrap the raw data in a new AnnData object to satisfy the type requirement for the .raw attribute
    adata_simple.raw = ad.AnnData(X=adata.raw.X.copy(), var=adata.raw.var.copy()) # type: ignore

print("AnnData object simplified.")

# %% [markdown]
# ## Save Cleaned Data

# %%
print(f"\nSaving cleaned and simplified dataset to {output_path}")
adata_simple.write_h5ad(output_path, compression='gzip') # type: ignore
print("Script finished successfully.")
# %%
