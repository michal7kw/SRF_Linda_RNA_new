#!/usr/bin/env Rscript
# %% [markdown]
# # 10. Create .cloupe File from Subclustered Seurat Object
#
# This script converts a subclustered Seurat object (.rds) into a
# Loupe file (.cloupe) for visualization in the Loupe Browser.

# %% [code]
# --- Configuration: SET THE CLUSTER ID TO PROCESS ---
CLUSTER_ID <- '0' # <-- CHANGE THIS to '1' or '16' to process other clusters
# The USE_FINAL_ANNOTATED_FILE flag is no longer needed.
# The script will now auto-detect the correct input file.

# --- Load Libraries ---
library(loupeR)
library(Seurat)

cat("Libraries loaded.\n")
cat("Processing for Cluster ID:", CLUSTER_ID, "\n")

# %% [code]
# --- Path Configuration ---
# Auto-detect environment (WSL vs. Windows)
is_wsl <- FALSE
try({
  if (Sys.info()["sysname"] == "Linux" && file.exists("/proc/version")) {
    proc_version <- readLines("/proc/version", n = 1)
    is_wsl <- grepl("microsoft", tolower(proc_version))
  }
}, silent = TRUE)

if (is_wsl) {
  project_dir <- "/mnt/d/Github/SRF_Linda_RNA"
  cat("Running in WSL environment, using path:", project_dir, "\n")
} else {
  project_dir <- "D:/Github/SRF_Linda_RNA"
  cat("Running in Windows environment, using path:", project_dir, "\n")
}

# Define input and output directories
seurat_input_dir <- file.path(project_dir, "combine_data", "results_subclustering")
output_dir <- seurat_input_dir

# Define input and output file paths based on CLUSTER_ID
# --- Auto-detect input file ---
# This makes the script more robust by checking for the final annotated file first,
# and falling back to the cleaned file if it's not found.
final_annotated_path <- file.path(seurat_input_dir, paste0("subcluster_", CLUSTER_ID, "_final_annotated.rds"))
cleaned_path <- file.path(seurat_input_dir, paste0("subcluster_", CLUSTER_ID, "_cleaned.rds"))

if (file.exists(final_annotated_path)) {
  seurat_input_file <- final_annotated_path
  output_name <- paste0("subcluster_", CLUSTER_ID, "_final_annotated_loupe")
  cat("Found FINAL annotated RDS file. Using it as input.\n")
} else if (file.exists(cleaned_path)) {
  seurat_input_file <- cleaned_path
  output_name <- paste0("subcluster_", CLUSTER_ID, "_cleaned_loupe")
  cat("Found CLEANED RDS file. Using it as input.\n")
} else {
  stop(paste("ERROR: Neither '..._final_annotated.rds' nor '..._cleaned.rds' file found for cluster", CLUSTER_ID))
}

cat("Configuration set.\n")
cat("Input Seurat RDS path:", seurat_input_file, "\n")
cat("Output Loupe name:", output_name, "\n")

# %% [code]
# --- Check Input File ---
cat("\nChecking for input Seurat file...\n")
if (!file.exists(seurat_input_file)) {
  stop(paste("ERROR: Input Seurat file not found at:", seurat_input_file))
}
cat("Input Seurat file confirmed at:", seurat_input_file, "\n")

# %% [code]
# --- Load Seurat Object ---
cat("\nLoading Seurat object from:", seurat_input_file, "\n")
seurat_obj <- NULL
tryCatch({
  seurat_obj <- readRDS(seurat_input_file)
  cat("Seurat object loaded successfully.\n")
  print(seurat_obj)
}, error = function(e) {
  stop(paste("ERROR loading Seurat object:", e$message))
})

# %% [code]
# --- Prepare Seurat Object for Loupe Conversion ---
if (is.null(seurat_obj)) {
    stop("Seurat object is NULL, cannot proceed.")
}

# Verify that the 'counts' slot exists and is not empty
default_assay <- DefaultAssay(seurat_obj)
cat("\nChecking for raw counts in default assay '", default_assay, "'...\n", sep="")
if (!"counts" %in% slotNames(seurat_obj@assays[[default_assay]])) {
    stop("ERROR: 'counts' slot not found. loupeR requires raw counts.")
}
counts_data <- GetAssayData(seurat_obj, assay = default_assay, layer = "counts")
if (nrow(counts_data) == 0 || ncol(counts_data) == 0) {
  stop("ERROR: 'counts' slot is present but empty.")
}
cat("'counts' slot verified.\n")

# Correct barcode formatting if necessary (from scanpy's -0, -1 suffixes)
cat("\nCorrecting cell barcodes for LoupeR compatibility if needed...\n")
original_barcodes <- colnames(seurat_obj)
# Regex to replace the last hyphen before final digits with an underscore
corrected_barcodes <- gsub("(.*)-(\\d+)$", "\\1_\\2", original_barcodes)
num_changed <- sum(original_barcodes != corrected_barcodes)

if (num_changed > 0) {
  cat(paste("Corrected", num_changed, "barcodes.\n"))
  seurat_obj <- RenameCells(seurat_obj, new.names = corrected_barcodes)
  cat("Cell names updated in the Seurat object.\n")
} else {
  cat("No barcode correction was necessary.\n")
}

# %% [code]
# --- Convert Seurat to Loupe ---
cat("\nStarting conversion to .cloupe file...\n")
cat("Output will be saved in:", output_dir, "\n")

tryCatch({
  create_loupe_from_seurat(
    obj = seurat_obj,
    output_dir = output_dir,
    output_name = output_name,
    force = TRUE # Overwrite existing file
  )
  cat("\nConversion successful! Loupe file created at:", file.path(output_dir, paste0(output_name, ".cloupe")), "\n")
}, error = function(e) {
  if (grepl("Error creating loupe file.*system command", e$message, ignore.case = TRUE)) {
      warning("Loupe conversion failed. This might be due to the loupeR executable not being found or the EULA not being accepted. Please run setup() interactively in R.", immediate. = TRUE)
  }
  stop(paste("ERROR during Loupe conversion:", e$message))
})

cat("\nScript finished for cluster", CLUSTER_ID, ".\n")