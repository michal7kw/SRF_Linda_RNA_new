# %% [markdown]
# # 12. Add Manual Annotations to Subclusters
#
# This script loads a CellTypist-annotated subcluster AnnData object and
# merges it with manual annotations provided in a CSV file.

# %%
import scanpy as sc
import pandas as pd
import os
import sys
import warnings

# %%
# --- Configuration: SET THE CLUSTER ID TO PROCESS ---
CLUSTER_ID = '0' # <-- CHANGE THIS to '1' or '16' to process other clusters

# --- User Input: SET THE PATH TO YOUR MANUAL ANNOTATION FILE ---
# You must provide the CSV file from MapMyCells.
MANUAL_ANNOTATION_CSV_PATH = f'models/subcluster_{CLUSTER_ID}_mapmycells_annotations.csv' # <-- UPDATE THIS PATH

# --- Script Configuration ---
# Define Project Directories
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = os.path.join(PROJECT_DIR, "combine_data")
os.chdir(WORKING_DIR)

# Define input/output directories and file paths
BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_subclustering")
INPUT_H5AD_FILE = os.path.join(BASE_RESULTS_DIR, f'subcluster_{CLUSTER_ID}_cleaned.h5ad')
OUTPUT_H5AD_FILE = os.path.join(BASE_RESULTS_DIR, f'subcluster_{CLUSTER_ID}_final_annotated.h5ad')

print(f"--- Processing Subcluster {CLUSTER_ID} ---")
print(f"Input AnnData: {INPUT_H5AD_FILE}")
print(f"Manual Annotation CSV: {MANUAL_ANNOTATION_CSV_PATH}")
print(f"Output AnnData: {OUTPUT_H5AD_FILE}")

# %% [markdown]
# ## Load Annotated Subcluster Data

# %%
print(f"\nLoading AnnData object from {INPUT_H5AD_FILE}...")
if not os.path.exists(INPUT_H5AD_FILE):
    raise FileNotFoundError(f"Input AnnData file not found at {INPUT_H5AD_FILE}. Please run script 13 (cleaning) first.")
adata = sc.read_h5ad(INPUT_H5AD_FILE)
print("AnnData object loaded.")

# %% [markdown]
# ## Load and Merge Manual Annotations

# %%
print(f"\nLoading manual annotations from {MANUAL_ANNOTATION_CSV_PATH}...")
if not os.path.exists(MANUAL_ANNOTATION_CSV_PATH):
    warnings.warn(f"Warning: Manual annotation file not found at {MANUAL_ANNOTATION_CSV_PATH}. Skipping this step.")
    # If the file doesn't exist, we can just save the input file to the output location and finish.
    adata.write(OUTPUT_H5AD_FILE)
    print(f"\nNo manual annotation file found. Copied input to {OUTPUT_H5AD_FILE} and exiting.")
    sys.exit(0)

df_manual = pd.read_csv(MANUAL_ANNOTATION_CSV_PATH)
print(f"Loaded {len(df_manual)} annotations from manual CSV.")

# --- Data Validation ---
# --- Data Validation ---
required_cols = ['cell_id', 'subclass_name', 'supertype_name']
if not all(col in df_manual.columns for col in required_cols):
    raise ValueError(f"Manual annotation CSV must contain the following columns: {', '.join(required_cols)}")

# The 'cell_id' from MapMyCells should match the AnnData index.
df_manual = df_manual.set_index('cell_id')

# --- Merge Annotations ---
new_cols = {
    'mapmycells_subclass_name': 'subclass_name',
    'mapmycells_supertype_name': 'supertype_name'
}

for new_col, source_col in new_cols.items():
    if new_col in adata.obs.columns:
        warnings.warn(f"Column '{new_col}' already exists. It will be overwritten.")
    
    adata.obs[new_col] = adata.obs.index.map(df_manual[source_col])
    
    # --- Verification ---
    mapped_count = adata.obs[new_col].notna().sum()
    total_cells = adata.n_obs
    print(f"Mapped {mapped_count}/{total_cells} cells for '{new_col}'.")

    # Fill unmapped cells with a default value, e.g., 'Unassigned'
    adata.obs[new_col] = adata.obs[new_col].fillna('Unassigned').astype('category')
    print(f"Filled unmapped cells in '{new_col}' with 'Unassigned'.")

print("\nUpdated adata.obs head:")
print(adata.obs[list(new_cols.keys())].head())

# %% [markdown]
# ## Save Final Annotated Data

# %%
print(f"\nSaving final annotated AnnData object to {OUTPUT_H5AD_FILE}...")
adata.write(OUTPUT_H5AD_FILE, compression='gzip')
print("Final annotated AnnData object saved successfully.")

print(f"\n--- Manual Annotation for Subcluster {CLUSTER_ID} Finished ---")

# %%
