import anndata
import os
import scipy.sparse
import math
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- Configuration ---
INPUT_H5AD_PATH = os.path.join('results_from_raw', 'final_annotation', 'merged_raw_final_annotated_simple.h5ad') # Input file from raw_5_finalize_annotation.py
OUTPUT_H5AD_PATH_BASE = 'merged_raw_final_annotated_simple_compressed' # Base name for output file(s)
MAX_FILE_SIZE_GB = 2
MAX_FILE_SIZE_BYTES = MAX_FILE_SIZE_GB * 1024 * 1024 * 1024 # 2 GB in bytes

# --- PART 1: DATA INPUT ---
logging.info(f"Loading AnnData object from: {INPUT_H5AD_PATH}")
try:
    adata = anndata.read_h5ad(INPUT_H5AD_PATH)
    logging.info("AnnData object loaded successfully.")
    logging.info(f"Original shape: {adata.shape}")
except FileNotFoundError:
    logging.error(f"Error: Input file not found at {INPUT_H5AD_PATH}")
    exit()
except Exception as e:
    logging.error(f"An error occurred while loading the AnnData object: {e}")
    exit()

# Extract necessary components (adjust if your .obs and .var index/columns are different)
count_matrix = adata.X
obs_data = adata.obs # Assuming index contains cell IDs
var_data = adata.var # Assuming index contains gene IDs (e.g., Ensembl IDs, gene symbols)

# --- PART 2: DATA OUTPUT ---

logging.info("Ensuring count matrix is in CSR format.")
# Convert count matrix to CSR sparse format if it isn't already
if not isinstance(count_matrix, scipy.sparse.csr_matrix):
    try:
        count_matrix = scipy.sparse.csr_matrix(count_matrix)
        logging.info("Converted count matrix to CSR format.")
    except Exception as e:
        logging.error(f"Failed to convert matrix to CSR format: {e}")
        exit()
else:
    logging.info("Count matrix is already in CSR format.")


# Create a new AnnData object with the potentially converted matrix
# Ensure obs and var indices match the matrix dimensions
if count_matrix.shape[0] != len(obs_data):
     logging.warning(f"Mismatch in number of observations (matrix rows: {count_matrix.shape[0]}, obs rows: {len(obs_data)}). Attempting to slice obs.")
     # Handle potential mismatch - this assumes obs needs slicing if matrix was sliced, adjust if needed
     # If the input adata wasn't modified, this shouldn't be necessary unless obs was detached
     # obs_data = obs_data.iloc[:count_matrix.shape[0]]

if count_matrix.shape[1] != len(var_data):
    logging.warning(f"Mismatch in number of variables (matrix columns: {count_matrix.shape[1]}, var rows: {len(var_data)}). Attempting to slice var.")
    # Handle potential mismatch - similar logic as above
    # var_data = var_data.iloc[:count_matrix.shape[1]]


# Create the AnnData object for saving
try:
    ad_to_save = anndata.AnnData(
        X=count_matrix,
        obs=obs_data.copy(), # Use copies to avoid modifying original adata's obs/var
        var=var_data.copy()
    )
    logging.info("Created new AnnData object for saving.")
except Exception as e:
    logging.error(f"Failed to create AnnData object for saving: {e}")
    exit()

# Define the initial output path
output_path_single = f"{OUTPUT_H5AD_PATH_BASE}.h5ad"

# Write the full data to a compressed file initially
logging.info(f"Writing compressed AnnData object to: {output_path_single}")
try:
    ad_to_save.write_h5ad(output_path_single, compression='gzip')
    logging.info("Compressed file written successfully.")
except Exception as e:
    logging.error(f"Failed to write compressed h5ad file: {e}")
    exit()

# Determine the file size
try:
    file_size_bytes = os.path.getsize(output_path_single)
    file_size_gb = file_size_bytes / (1024 * 1024 * 1024)
    logging.info(f"Compressed file size: {file_size_bytes} bytes ({file_size_gb:.2f} GB)")
except FileNotFoundError:
    logging.error(f"Error: Could not find the written file {output_path_single} to check its size.")
    exit()

# Check if file size exceeds the limit
if file_size_bytes > MAX_FILE_SIZE_BYTES:
    logging.info(f"File size exceeds {MAX_FILE_SIZE_GB}GB limit. Splitting into multiple files.")

    # Determine the number of files needed and rows per file
    num_files_needed = math.ceil(file_size_bytes / MAX_FILE_SIZE_BYTES)
    num_total_rows = ad_to_save.shape[0]
    rows_per_subset = num_total_rows // num_files_needed

    logging.info(f"Splitting into {num_files_needed} files with approximately {rows_per_subset} rows each.")
    if num_files_needed > 8:
        logging.warning("Splitting into more than 8 files. Consider using a code-based mapping approach if possible.")

    # Write each chunk to a separate file
    split_success = True
    for i in range(num_files_needed):
        start_idx = i * rows_per_subset
        # Ensure the last chunk includes any remaining rows
        end_idx = (i + 1) * rows_per_subset if i < num_files_needed - 1 else num_total_rows

        output_path_part = f"{OUTPUT_H5AD_PATH_BASE}_part_{i+1}_of_{num_files_needed}.h5ad"
        logging.info(f"Writing part {i+1}/{num_files_needed} (rows {start_idx} to {end_idx}) to {output_path_part}")

        try:
            # Create AnnData object for the chunk
            # Slicing CSR matrix preserves sparsity
            ad_chunk = anndata.AnnData(
                X=count_matrix[start_idx:end_idx, :],
                obs=obs_data.iloc[start_idx:end_idx].copy(), # Slice obs dataframe
                var=var_data.copy() # Var remains the same for all chunks
            )
            # Write the chunk with compression
            ad_chunk.write_h5ad(output_path_part, compression='gzip')
            logging.info(f"Successfully wrote {output_path_part}")
        except Exception as e:
            logging.error(f"Failed to write chunk {i+1} to {output_path_part}: {e}")
            split_success = False
            break # Stop splitting if one part fails

    if split_success:
        logging.info("Successfully split the data into multiple files.")
        # Optionally remove the original large file if splitting was successful
        logging.info(f"Removing the original large file: {output_path_single}")
        try:
            os.remove(output_path_single)
        except OSError as e:
            logging.warning(f"Could not remove original large file {output_path_single}: {e}")
    else:
        logging.error("File splitting failed. The original large file was kept.")

else:
    logging.info(f"File size is within the {MAX_FILE_SIZE_GB}GB limit. No splitting needed.")
    logging.info(f"Final output file: {output_path_single}")

logging.info("Script finished.")
