# %% [markdown]
# # 11. Annotate Subclusters with CellTypist
#
# This script loads a processed subcluster dataset, annotates it using
# CellTypist models, and saves the newly annotated data.

# %%
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import warnings
from celltypist import models, annotate
import random

# %%
# --- Configuration: SET THE CLUSTER ID TO PROCESS ---
CLUSTER_ID = '1' # <-- CHANGE THIS to '1' or '16' to process other clusters

# --- Script Configuration ---
# Set seeds for reproducibility
random_seed = 0
np.random.seed(random_seed)
random.seed(random_seed)
sc.settings.seed = random_seed # type: ignore

# Define Project Directories
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = os.path.join(PROJECT_DIR, "combine_data")
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

# Define models for annotation
MODELS = {
    "Mouse_Isocortex_Hippocampus": "Cell types from the adult mouse isocortex and hippocampal formation",
    "Mouse_Dentate_Gyrus": "Cell types from the dentate gyrus in perinatal, juvenile, and adult mice"
}

# Set up directories
BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_subclustering")
INPUT_H5AD_FILE = os.path.join(BASE_RESULTS_DIR, f'subcluster_{CLUSTER_ID}_processed.h5ad')
MODEL_DIR = os.path.join(WORKING_DIR, "models")
OUTPUT_H5AD_FILE = os.path.join(BASE_RESULTS_DIR, f'subcluster_{CLUSTER_ID}_annotated.h5ad')
PLOT_DIR = os.path.join(BASE_RESULTS_DIR, "plots_subclustering_annotated")
os.makedirs(PLOT_DIR, exist_ok=True)
sc.settings.figdir = PLOT_DIR

# CellTypist parameters
MAJORITY_VOTING = True
PROB_THRESHOLD = 0.5

# %% [markdown]
# ## Load Subclustered Data

# %%
print(f"--- Processing Subcluster {CLUSTER_ID} ---")
print(f"Loading dataset from {INPUT_H5AD_FILE}...")
if not os.path.exists(INPUT_H5AD_FILE):
    raise FileNotFoundError(f"Input file not found: {INPUT_H5AD_FILE}. Please run the subclustering script first.")

adata = sc.read_h5ad(INPUT_H5AD_FILE)
print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes.")

# %% [markdown]
# ## Prepare Data for CellTypist
# CellTypist requires log1p normalized data from raw counts.

# %%
adata_ct = adata.copy()

if adata_ct.raw is None:
    raise ValueError("AnnData object must contain raw counts in `.raw` for CellTypist preparation.")

print("Preparing data for CellTypist...")
# 1. Use raw counts for the current genes
adata_ct.X = adata_ct.raw[:, adata_ct.var_names].X.copy() # type: ignore
# 2. Normalize total counts to 10,000
sc.pp.normalize_total(adata_ct, target_sum=10000)
# 3. Log1p transform
sc.pp.log1p(adata_ct)
print("Data preparation complete.")

# %% [markdown]
# ## Run CellTypist Annotation

# %%
print(f"\nRunning CellTypist annotation with majority_voting={MAJORITY_VOTING}, prob_threshold={PROB_THRESHOLD}")

for model_name, description in MODELS.items():
    print(f"\nProcessing model: {model_name}")
    model_path = os.path.join(MODEL_DIR, f"{model_name}.pkl")
    if not os.path.exists(model_path):
        print(f"Model file not found at {model_path}. Attempting to download...")
        try:
            models.celltypist_models(force_update=False) # type: ignore
            model = models.Model.load(model_name)
            model.write(model_path)
            print(f"Model downloaded and saved to {model_path}")
        except Exception as e:
            raise FileNotFoundError(f"Model file not found and download failed: {e}")
    
    print("Running annotation...")
    predictions = annotate(adata_ct, model=model_path, majority_voting=MAJORITY_VOTING, p_thres=PROB_THRESHOLD)
    
    # Add annotations to the temporary adata_ct object
    prefix = 'ISO_' if 'Isocortex' in model_name else 'DG_'
    predictions.to_adata(adata_ct, prefix=prefix) # type: ignore
    
    # Define the new columns that were added
    majority_voting_col = f"{prefix}majority_voting"
    conf_score_col = f"{prefix}conf_score"
    
    # Copy the new annotation columns from adata_ct back to the original adata object
    if majority_voting_col in adata_ct.obs.columns:
        adata.obs[majority_voting_col] = adata_ct.obs[majority_voting_col]
        print(f"Added '{majority_voting_col}' to the main AnnData object.")
    if conf_score_col in adata_ct.obs.columns:
        adata.obs[conf_score_col] = adata_ct.obs[conf_score_col]
        print(f"Added '{conf_score_col}' to the main AnnData object.")

    # Visualize results on the original UMAP
    if majority_voting_col in adata.obs.columns:
        sc.pl.umap(adata, color=majority_voting_col,
                   legend_loc='on data', legend_fontsize='xx-small',
                   title=f"Subcluster {CLUSTER_ID} - {model_name}",
                   save=f"_subcluster_{CLUSTER_ID}_{model_name}_celltype.png", show=True)

# %% [markdown]
# ## Save Annotated Data

# %%
print(f"\nSaving annotated subcluster data to {OUTPUT_H5AD_FILE}")
try:
    adata.write(OUTPUT_H5AD_FILE)
    print("Successfully saved annotated AnnData object.")
except Exception as e:
    print(f"Error saving AnnData object: {e}")

print(f"\n--- Annotation for Subcluster {CLUSTER_ID} Finished ---")
# %%
