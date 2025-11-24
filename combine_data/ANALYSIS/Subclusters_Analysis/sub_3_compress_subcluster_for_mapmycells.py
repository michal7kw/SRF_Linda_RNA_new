# %% [markdown]
# # 14. Compress Subcluster for MapMyCells
#
# This script prepares a cleaned subcluster AnnData object for submission to
# the MapMyCells service. It extracts raw counts, ensures the matrix is in
# CSR format, and saves a minimal, compressed .h5ad file. It also handles
# splitting the file if it exceeds a size limit.

# %%
import anndata
import os
import sys
import scipy.sparse
import math
import logging

# %%
# --- Configuration ---
CLUSTER_ID = '0' # <-- CHANGE THIS to '1' or '16' to process other clusters
MAX_FILE_SIZE_GB = 2

# --- Script Configuration ---
# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define Project Directories
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = os.path.join(PROJECT_DIR, "combine_data")
os.chdir(WORKING_DIR)

# Define input/output directories and file paths
BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_subclustering")
INPUT_H5AD_PATH = os.path.join(BASE_RESULTS_DIR, f'subcluster_{CLUSTER_ID}_cleaned.h5ad')

# Define a dedicated output directory for MapMyCells files
MAPMYCELLS_OUTPUT_DIR = os.path.join(BASE_RESULTS_DIR, "for_mapmycells")
os.makedirs(MAPMYCELLS_OUTPUT_DIR, exist_ok=True)
OUTPUT_H5AD_PATH_BASE = os.path.join(MAPMYCELLS_OUTPUT_DIR, f'subcluster_{CLUSTER_ID}_for_mapmycells')

MAX_FILE_SIZE_BYTES = MAX_FILE_SIZE_GB * 1024 * 1024 * 1024

print(f"--- Compressing Subcluster {CLUSTER_ID} for MapMyCells ---")
print(f"Input file: {INPUT_H5AD_PATH}")
print(f"Output base: {OUTPUT_H5AD_PATH_BASE}")

# %% [markdown]
# ## Load Cleaned Subcluster Data

# %%
logging.info(f"Loading AnnData object from: {INPUT_H5AD_PATH}")
try:
    adata = anndata.read_h5ad(INPUT_H5AD_PATH)
    logging.info("AnnData object loaded successfully.")
except FileNotFoundError:
    logging.error(f"Error: Input file not found at {INPUT_H5AD_PATH}. Please run the cleaning script (13) first.")
    sys.exit(1)

# %% [markdown]
# ## Prepare Data for MapMyCells

# %%
# 1. Extract raw counts
if adata.raw is None:
    logging.error("Error: .raw attribute not found in the AnnData object. Cannot proceed.")
    sys.exit(1)

logging.info("Extracting raw counts from .raw attribute.")
count_matrix = adata.raw.X
obs_data = adata.obs
var_data = adata.raw.var # Use var from .raw to match the matrix

# 2. Ensure CSR format
logging.info("Ensuring count matrix is in CSR format.")
if not isinstance(count_matrix, scipy.sparse.csr_matrix):
    count_matrix = scipy.sparse.csr_matrix(count_matrix)
    logging.info("Converted count matrix to CSR format.")
else:
    logging.info("Count matrix is already in CSR format.")

# 3. Create a new, minimal AnnData object
logging.info("Creating a new minimal AnnData object for submission.")
ad_to_save = anndata.AnnData(
    X=count_matrix,
    obs=obs_data.copy(),
    var=var_data.copy()
)

# %% [markdown]
# ## Save and Split Compressed File

# %%
# Write the full data to a compressed file initially to check its size
output_path_single = f"{OUTPUT_H5AD_PATH_BASE}.h5ad"
logging.info(f"Writing initial compressed file to: {output_path_single}")
ad_to_save.write_h5ad(output_path_single, compression='gzip')

# Check file size
file_size_bytes = os.path.getsize(output_path_single)
logging.info(f"Compressed file size: {file_size_bytes / (1024**3):.2f} GB")

# Split if necessary
if file_size_bytes > MAX_FILE_SIZE_BYTES:
    logging.info(f"File size exceeds {MAX_FILE_SIZE_GB}GB limit. Splitting into multiple files.")
    num_files_needed = math.ceil(file_size_bytes / MAX_FILE_SIZE_BYTES)
    rows_per_subset = ad_to_save.shape[0] // num_files_needed
    logging.info(f"Splitting into {num_files_needed} files with ~{rows_per_subset} rows each.")

    for i in range(num_files_needed):
        start_idx = i * rows_per_subset
        end_idx = (i + 1) * rows_per_subset if i < num_files_needed - 1 else ad_to_save.shape[0]
        output_path_part = f"{OUTPUT_H5AD_PATH_BASE}_part_{i+1}_of_{num_files_needed}.h5ad"
        logging.info(f"Writing part {i+1}/{num_files_needed} to {output_path_part}")
        
        ad_chunk = anndata.AnnData(
            X=count_matrix[start_idx:end_idx, :],
            obs=obs_data.iloc[start_idx:end_idx].copy(),
            var=var_data.copy()
        )
        ad_chunk.write_h5ad(output_path_part, compression='gzip')

    logging.info(f"Successfully split data. Removing original large file: {output_path_single}")
    os.remove(output_path_single)
else:
    logging.info(f"File size is within the limit. No splitting needed.")
    logging.info(f"Final output file: {output_path_single}")

logging.info("Script finished successfully.")