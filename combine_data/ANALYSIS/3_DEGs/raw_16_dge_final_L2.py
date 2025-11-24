# %% [markdown]
# # Environment Setup and Data Loading

# %%
import scanpy as sc
import numpy as np
import pandas as pd
import os
import sys
import anndata as ad
import random
import functions_degs
import importlib

PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = os.path.join(PROJECT_DIR, "combine_data")
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

importlib.reload(functions_degs)

# Set seeds for all random number generators
random_seed = 0
np.random.seed(random_seed)
random.seed(random_seed)

# %%
# Set up directories

# REMOVE_DOUBLETS = True
REMOVE_DOUBLETS = False

FIX_TRESHOLD = True
# FIX_TRESHOLD = False

if FIX_TRESHOLD:
    BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")
else:
    if REMOVE_DOUBLETS:
        BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw_percentile_threshold", "doublets_removed")
    else:
        BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw_percentile_threshold")

INPUT_DIR = BASE_RESULTS_DIR
adata_path = os.path.join(INPUT_DIR, 'annotation_final.h5ad')

PARENT_OUTPUT_DIR = os.path.join(INPUT_DIR, "DEGs_cell_type_L2")

# CUSTOM_ANALYSIS =  None
CUSTOM_ANALYSIS =  "FC_0_25"

DEG_BY = "cell_type_L2"

if CUSTOM_ANALYSIS is not None:
    PARENT_OUTPUT_DIR = PARENT_OUTPUT_DIR + CUSTOM_ANALYSIS

PLOT_OUTPUT_DIR = os.path.join(PARENT_OUTPUT_DIR, 'plots')
DGE_OUTPUT_DIR = os.path.join(PARENT_OUTPUT_DIR, 'dge_res')
BIOMARKER_OUTPUT_DIR = os.path.join(PARENT_OUTPUT_DIR, 'biomarkers')

# Create output directories if they don't exist
os.makedirs(PLOT_OUTPUT_DIR, exist_ok=True)
os.makedirs(DGE_OUTPUT_DIR, exist_ok=True)
os.makedirs(BIOMARKER_OUTPUT_DIR, exist_ok=True)

print(f"Input directory: {INPUT_DIR}")
print(f"Plot output directory (Overall/Genotype): {PLOT_OUTPUT_DIR}")
print(f"DGE output directory (Overall/Genotype): {DGE_OUTPUT_DIR}")
print(f"Biomarker output directory: {BIOMARKER_OUTPUT_DIR}")

# Configure scanpy settings
sc.settings.figdir = PLOT_OUTPUT_DIR
sc.settings.set_figure_params(dpi=150, facecolor='white')

# %%
# Load data
print(f"\nLoading dataset from {adata_path}")
adata = sc.read_h5ad(adata_path)


# %% [markdown]
# # Data Preprocessing

# %%
# Save current settings
original_max_rows = pd.get_option('display.max_rows')

# Set to None to display all rows
pd.set_option('display.max_rows', None)

# Check original clusters:
print(f"Original {DEG_BY} clusters:")
print(adata.obs[DEG_BY].value_counts())

# Reset to original settings
pd.set_option('display.max_rows', original_max_rows)

# %%
# Check adata structure
print("\nAnnData object summary:")
print(adata)
print("\nAvailable layers:", list(adata.layers.keys()))
print("Raw data available:", adata.raw is not None)
if adata.raw:
    print("Raw data shape:", adata.raw.X.shape)

# %%
# Create the 'for_DEGs' layer: Normalized and log1p transformed counts
# Assuming raw counts are in adata.raw.X
print("\nCreating 'for_DEGs' layer...")
if adata.raw is not None and adata.raw.X is not None:
    # Create a temporary AnnData with raw counts to perform normalization and log1p
    adata_for_dge = ad.AnnData(adata.raw.X.copy()) # type: ignore
    adata_for_dge.obs_names = adata.obs_names # type: ignore
    adata_for_dge.var_names = adata.raw.var_names # type: ignore

    print("Normalizing total counts (target_sum=1e4)...")
    sc.pp.normalize_total(adata_for_dge, target_sum=1e4)

    print("Applying log1p transformation...")
    sc.pp.log1p(adata_for_dge)

    # Ensure the gene order matches the main adata object
    adata_for_dge = adata_for_dge[:, adata.var_names].copy() # type: ignore

    print("Storing result in adata.layers['for_DEGs']...")
    adata.layers['for_DEGs'] = adata_for_dge.X.copy() # type: ignore
    print("'for_DEGs' layer created with shape:", adata.layers['for_DEGs'].shape)
else:
    print("Warning: adata.raw.X not found. Cannot create 'for_DEGs' layer from raw counts.")
    print("DGE analysis will proceed using adata.X, which might be scaled.")

# %%
print("\nMetadata check:")
print("condition:", list(adata.obs.condition.unique()))
print("genotype:", list(adata.obs.genotype.unique()))
print(f"{DEG_BY}:", list(adata.obs[DEG_BY].unique()))

# %% [markdown]
# # Differential Gene Expression Analysis

# %%
# Overall DGE (Mutant vs Control) using the 'for_DEGs' layer, grouped by the DEG_BY
print(f"\nRunning Overall DGE (grouped by {DEG_BY})...")
dge_results = functions_degs.run_overall_dge(
    adata,
    grouping_key=DEG_BY,
    dge_output_dir=DGE_OUTPUT_DIR,
    plot_output_dir=PLOT_OUTPUT_DIR,
    layer='for_DEGs' if 'for_DEGs' in adata.layers else None
)

# %%
# Genotype-Specific DGE (Mutant vs Control within each genotype) using the 'for_DEGs' layer, grouped by the DEG_BY
print(f"\nRunning Genotype-Specific DGE (grouped by {DEG_BY})...")
dge_by_genotype = functions_degs.run_genotype_specific_dge(
    adata,
    grouping_key=DEG_BY,
    dge_output_dir=DGE_OUTPUT_DIR,
    plot_output_dir=PLOT_OUTPUT_DIR,
    layer='for_DEGs' if 'for_DEGs' in adata.layers else None
)

# %%
# Genotype Comparison DGE (Nestin vs Emx1 within each condition) using the 'for_DEGs' layer, grouped by the DEG_BY
print(f"\nRunning Genotype Comparison DGE (grouped by {DEG_BY})...")
dge_genotype_within_condition = functions_degs.run_genotype_comparison_dge(
    adata,
    grouping_key=DEG_BY,
    dge_output_dir=DGE_OUTPUT_DIR,
    plot_output_dir=PLOT_OUTPUT_DIR,
    layer='for_DEGs' if 'for_DEGs' in adata.layers else None
)

# %%
# Cell Type Comparison DGE (Each cell type vs Rest) using the 'for_DEGs' layer, using the DEG_BY
print(f"\nRunning Cell Type Comparison DGE (Markers for {DEG_BY})...")
dge_markers = functions_degs.run_cluster_comparison_dge(
    adata,
    grouping_key=DEG_BY, 
    dge_output_dir=BIOMARKER_OUTPUT_DIR, 
    plot_output_dir=BIOMARKER_OUTPUT_DIR,
    layer='for_DEGs' if 'for_DEGs' in adata.layers else None,
    method='wilcoxon' # or 't-test'
)

# %%
print("\nAll analysis and output generation complete.")
print(f"Overall/Genotype DGE outputs saved in: {DGE_OUTPUT_DIR}")
print(f"Biomarker DGE outputs saved in: {BIOMARKER_OUTPUT_DIR}")
