# %%
import scanpy as sc
import pandas as pd
import os
import warnings

# %%
# Define project directory and working directory
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = os.path.join(PROJECT_DIR, "combine_data")
os.chdir(WORKING_DIR)

# Define input/output directories and file paths
# Input AnnData file (output from raw_5_clean_data.py)

# REMOVE_DOUBLETS = True
REMOVE_DOUBLETS = False

# Set up directories
if REMOVE_DOUBLETS:
    BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw", "doublets_removed")
else:
    BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")

INPUT_DIR = BASE_RESULTS_DIR
ADATA_INPUT_PATH = os.path.join(INPUT_DIR, 'annotated_cleaned.h5ad')

# Input CSV annotation files
MODELS_DIR = os.path.join(WORKING_DIR, "models")
FIRST_LAYER_CSV_PATH = os.path.join(MODELS_DIR, 'first_layer_LB_cleaned.csv')
SECOND_LAYER_CSV_PATH = os.path.join(MODELS_DIR, 'second_layer_LB_cleaned.csv')

# Output AnnData file path
OUTPUT_DIR = INPUT_DIR
ADATA_OUTPUT_PATH = os.path.join(OUTPUT_DIR, 'mapmycells.h5ad')

print(f"Input AnnData: {ADATA_INPUT_PATH}")
print(f"First layer CSV: {FIRST_LAYER_CSV_PATH}")
print(f"Second layer CSV: {SECOND_LAYER_CSV_PATH}")
print(f"Output AnnData: {ADATA_OUTPUT_PATH}")

# %%
# Load the AnnData object
print(f"\nLoading AnnData object from {ADATA_INPUT_PATH}...")
if not os.path.exists(ADATA_INPUT_PATH):
    raise FileNotFoundError(f"Input AnnData file not found at {ADATA_INPUT_PATH}")
adata = sc.read_h5ad(ADATA_INPUT_PATH)
print("AnnData object loaded:")
print(adata)

# %%
# Load the first layer annotations
print(f"\nLoading first layer annotations from {FIRST_LAYER_CSV_PATH}...")
if not os.path.exists(FIRST_LAYER_CSV_PATH):
    raise FileNotFoundError(f"First layer CSV file not found at {FIRST_LAYER_CSV_PATH}")
df_first_layer = pd.read_csv(FIRST_LAYER_CSV_PATH)
print(f"Loaded {len(df_first_layer)} annotations from first layer CSV.")
print(df_first_layer.head())

# Check expected columns
if 'Barcode' not in df_first_layer.columns or 'first layer LB' not in df_first_layer.columns:
    raise ValueError("First layer CSV must contain 'Barcode' and 'first layer LB' columns.")

# Modify Barcode: Replace the last underscore with a hyphen
print("Modifying first layer barcodes (replacing last '_' with '-')...")
df_first_layer['Barcode'] = df_first_layer['Barcode'].str.replace(r'_([^_]*)$', r'-\1', regex=True)
print(df_first_layer.head()) # Show modified barcodes

# Set Barcode as index
df_first_layer = df_first_layer.set_index('Barcode')

# %%
# Load the second layer annotations
print(f"\nLoading second layer annotations from {SECOND_LAYER_CSV_PATH}...")
if not os.path.exists(SECOND_LAYER_CSV_PATH):
    raise FileNotFoundError(f"Second layer CSV file not found at {SECOND_LAYER_CSV_PATH}")
df_second_layer = pd.read_csv(SECOND_LAYER_CSV_PATH)
print(f"Loaded {len(df_second_layer)} annotations from second layer CSV.")
print(df_second_layer.head())

# Check expected columns
if 'Barcode' not in df_second_layer.columns or 'second layer LB' not in df_second_layer.columns:
    raise ValueError("Second layer CSV must contain 'Barcode' and 'second layer LB' columns.")

# Modify Barcode: Replace the last underscore with a hyphen
print("Modifying second layer barcodes (replacing last '_' with '-')...")
df_second_layer['Barcode'] = df_second_layer['Barcode'].str.replace(r'_([^_]*)$', r'-\1', regex=True)
print(df_second_layer.head()) # Show modified barcodes

# Set Barcode as index
df_second_layer = df_second_layer.set_index('Barcode')

# %%
# Add annotations to adata.obs
# Use .map() to align annotations based on the index (Barcode)
# This handles cases where barcodes might be in the CSV but not in adata, or vice-versa
print("\nAdding annotations to adata.obs...")

# Check if columns already exist and warn if overwriting
new_col_first = 'mapmycells_first_layer'
new_col_second = 'mapmycells_second_layer'

if new_col_first in adata.obs.columns:
    warnings.warn(f"Column '{new_col_first}' already exists in adata.obs. It will be overwritten.")
if new_col_second in adata.obs.columns:
    warnings.warn(f"Column '{new_col_second}' already exists in adata.obs. It will be overwritten.")

adata.obs[new_col_first] = adata.obs.index.map(df_first_layer['first layer LB'])
adata.obs[new_col_second] = adata.obs.index.map(df_second_layer['second layer LB'])

# Check how many annotations were successfully mapped
mapped_first = adata.obs[new_col_first].notna().sum()
mapped_second = adata.obs[new_col_second].notna().sum()
total_cells = adata.n_obs

print(f"Mapped {mapped_first}/{total_cells} cells for '{new_col_first}'.")
print(f"Mapped {mapped_second}/{total_cells} cells for '{new_col_second}'.")

# Optionally, fill NaN values if needed (e.g., with 'Unknown')
# adata.obs[new_col_first] = adata.obs[new_col_first].fillna('Unknown')
# adata.obs[new_col_second] = adata.obs[new_col_second].fillna('Unknown')

print("\nUpdated adata.obs head:")
print(adata.obs[[new_col_first, new_col_second]].head())

# %%
# Save the updated AnnData object
print(f"\nSaving updated AnnData object to {ADATA_OUTPUT_PATH}...")
# Ensure output directory exists
os.makedirs(os.path.dirname(ADATA_OUTPUT_PATH), exist_ok=True)
adata.write_h5ad(ADATA_OUTPUT_PATH, compression='gzip') # type: ignore
print("Updated AnnData object saved successfully.")

# %%
print("\nScript finished.")
# %%
