# %%
import scanpy as sc
import numpy as np
import pandas as pd
import os
import sys
import anndata as ad
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

# REMOVE_DOUBLETS = True
REMOVE_DOUBLETS = False


# %%
# Set up directories based on the output of the previous script (1_annotate_merged_raw.py)
if REMOVE_DOUBLETS:
    BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw", "doublets_removed")
else:
    BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")

INPUT_DIR = BASE_RESULTS_DIR # Use output from script 1
OUTPUT_DIR = BASE_RESULTS_DIR

adata_path = os.path.join(INPUT_DIR, 'annotated.h5ad') # Input is the annotated file from script 1
output_path = os.path.join(OUTPUT_DIR, 'annotated_cleaned.h5ad')

os.makedirs(OUTPUT_DIR, exist_ok=True)

LEIDEN_KEY = 'leiden_0.4'
UMAP_KEY = 'X_umap'

# %%
# Load merged data from script 1
print(f"\nLoading annotated dataset from {adata_path}")
adata = sc.read_h5ad(adata_path)

# %%
# Filter out sex-related genes
genes_to_remove = ["Xist", "Tsix", "Eif2s3x", "Ddx3x", "Kdm5c", "Kdm6a", "Rps6ka3", "Uty", "Ddx3y", "Eif2s3y", "Kdm5d", "Rps4y1", "Jarid1d"]
print(f"\nFiltering out {len(genes_to_remove)} sex-related genes...")

# Check which genes are actually present in the data
genes_present = adata.var_names.intersection(genes_to_remove)
print(f"Found {len(genes_present)} sex-related genes in adata: {list(genes_present)}")

if len(genes_present) > 0:
    print(f"Original number of genes: {adata.n_vars}")
    adata = adata[:, ~adata.var_names.isin(genes_present)].copy()
    print(f"Number of genes after filtering: {adata.n_vars}")

    # Filter raw data if it exists
    if adata.raw is not None:
        # Check if raw.var is defined and has names
        if hasattr(adata.raw, 'var_names') and adata.raw.var_names is not None:
            print(f"Original number of genes in raw: {adata.raw.n_vars}")
            # Filter raw based on its own var_names list
            raw_genes_present_in_raw = adata.raw.var_names.intersection(genes_to_remove)
            print(f"Found {len(raw_genes_present_in_raw)} sex-related genes in adata.raw: {list(raw_genes_present_in_raw)}")
            if len(raw_genes_present_in_raw) > 0:
                 # Create a new AnnData object for the filtered raw data
                 adata_raw_filtered = adata.raw[:, ~adata.raw.var_names.isin(raw_genes_present_in_raw)].copy()
                 adata.raw = ad.AnnData(X=adata_raw_filtered.X, var=adata_raw_filtered.var, obs=adata.obs[[]]) # type: ignore
                 print(f"Number of genes in raw after filtering: {adata.raw.n_vars}")
            else:
                print("No specified sex-related genes found in adata.raw.var_names.")
        else:
             print("adata.raw.var_names not found or is None. Skipping raw data filtering.")
else:
    print("No specified sex-related genes found in adata.var_names. No filtering applied to main data.")

# %%
print(adata)

# %%
# Simplify the AnnData object
print("\nSimplifying the AnnData object...")

# Define columns/keys to keep
obs_keep = ['sample', 'condition', 'genotype', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', LEIDEN_KEY, 'ISO_majority_voting', 'ISO_conf_score', 'DG_majority_voting', 'DG_conf_score']
var_keep = ['mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'n_cells', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'mean', 'std']
uns_keep = ['DG_majority_voting_colors', 'ISO_majority_voting_colors', 'condition_colors', 'genotype_colors', 'hvg', LEIDEN_KEY, 'leiden_0.8_colors', 'log1p', 'neighbors', 'pca', 'sample_colors', 'umap']
obsm_keep = ['X_pca', 'X_umap']
varm_keep = ['PCs']
obsp_keep = ['connectivities', 'distances']

# Create a new AnnData object with selected slots
adata_simple = ad.AnnData(X=adata.X.copy(), # type: ignore
                          obs=adata.obs[obs_keep].copy(),
                          var=adata.var[var_keep].copy())

# Copy selected unstructured annotation
adata_simple.uns = {k: adata.uns[k] for k in uns_keep if k in adata.uns}

# Copy selected multidimensional annotation
adata_simple.obsm = {k: adata.obsm[k] for k in obsm_keep if k in adata.obsm}

# Copy selected variance annotation
adata_simple.varm = {k: adata.varm[k] for k in varm_keep if k in adata.varm}

# Copy selected pairwise annotation
adata_simple.obsp = {k: adata.obsp[k] for k in obsp_keep if k in adata.obsp}

# Copy raw data if it exists
if adata.raw is not None:
    print("Copying raw data...")
    # Ensure raw.var has the correct columns (usually just gene names/IDs)
    # It might be different from the main var dataframe
    # Create a minimal obs DataFrame with the index from adata_simple to set obs_names correctly
    raw_obs_minimal = pd.DataFrame(index=adata_simple.obs_names)
    adata_simple.raw = ad.AnnData(X=adata.raw.X.copy(), # type: ignore
                                  var=adata.raw.var.copy(),
                                  obs=raw_obs_minimal)


# Save the simplified AnnData object

print(f"\nSaving simplified annotated dataset to {output_path}")
adata_simple.write_h5ad(output_path, compression='gzip') # type: ignore