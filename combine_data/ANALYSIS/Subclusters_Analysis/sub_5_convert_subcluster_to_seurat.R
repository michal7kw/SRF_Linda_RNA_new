#!/usr/bin/env Rscript
# %% [markdown]
# # 9. Convert Subclustered H5AD to Seurat Object
#
# This script converts a specific subclustered AnnData object (.h5ad)
# into a Seurat object (.rds), preparing it for further analysis or
# conversion to other formats like .cloupe.

# %% [code]
# --- Configuration: SET THE CLUSTER ID TO PROCESS ---
CLUSTER_ID <- '0' # <-- CHANGE THIS to '1' or '16' to process other clusters
USE_FINAL_ANNOTATED_FILE <- TRUE # <-- SET TO FALSE if you skipped the manual annotation step (script 12)

# --- Load Libraries ---
library(schard)
library(Seurat)
library(Matrix)

cat("Libraries loaded.\n")

# --- Dynamic Keys ---
SUBCLUSTER_KEY <- paste0('leiden_sub_', CLUSTER_ID)
cat("Processing for Cluster ID:", CLUSTER_ID, "\n")
cat("Using subcluster metadata key:", SUBCLUSTER_KEY, "\n")
MANUAL_ANNOTATION_KEY <- paste0('manual_annotation_subcluster_', CLUSTER_ID)

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
input_dir <- file.path(project_dir, "combine_data", "results_subclustering")
output_dir <- input_dir

# Define input and output file paths based on CLUSTER_ID
# Define input and output file paths based on whether the final annotated file is used
if (USE_FINAL_ANNOTATED_FILE) {
  input_h5ad_path <- file.path(input_dir, paste0("subcluster_", CLUSTER_ID, "_final_annotated.h5ad"))
  output_seurat_path <- file.path(output_dir, paste0("subcluster_", CLUSTER_ID, "_final_annotated.rds"))
  cat("Using FINAL annotated H5AD file (includes manual annotations).\n")
} else {
  # If not using manual annotations, the input is the cleaned file from script 13
  input_h5ad_path <- file.path(input_dir, paste0("subcluster_", CLUSTER_ID, "_cleaned.h5ad"))
  output_seurat_path <- file.path(output_dir, paste0("subcluster_", CLUSTER_ID, "_cleaned.rds"))
  cat("Using CLEANED H5AD file (no manual annotations).\n")
}

cat("Configuration set.\n")
cat("Input H5AD path:", input_h5ad_path, "\n")
cat("Output Seurat path:", output_seurat_path, "\n")

# %% [code]
# --- Check Input File ---
cat("\nChecking for input file...\n")
if (!file.exists(input_h5ad_path)) {
  stop(paste("ERROR: Input file not found at:", input_h5ad_path))
}
cat("Input file confirmed at:", input_h5ad_path, "\n")

# %% [code]
# --- Load h5ad into Seurat ---
cat("\nLoading h5ad file using schard::h5ad2seurat...\n")
cat("Using use.raw = TRUE to load raw counts from the .raw slot.\n")

seurat_obj <- NULL
tryCatch({
  seurat_obj <- schard::h5ad2seurat(input_h5ad_path, use.raw = TRUE)
  cat("\nAnnData object loaded into Seurat object.\n")
}, error = function(e) {
  stop(paste("ERROR during h5ad loading:", e$message))
})

# %% [code]
# --- Verification ---
cat("\n--- Seurat Object Summary (Post-Loading) ---\n")
if (is.null(seurat_obj)) {
    stop("ERROR: Seurat object is NULL after loading.")
}
print(seurat_obj)

# Verify assay and counts
default_assay <- DefaultAssay(seurat_obj)
cat("\nDefault Assay:", default_assay, "\n")
if (!"counts" %in% slotNames(seurat_obj@assays[[default_assay]])) {
    stop(paste("ERROR: 'counts' slot not found in the '", default_assay, "' assay."))
}
counts_matrix <- GetAssayData(seurat_obj, assay = default_assay, layer = "counts")
if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
    stop(paste("ERROR: 'counts' slot in assay '", default_assay, "' is empty."))
}
cat("'counts' slot is populated (Dimensions: ", nrow(counts_matrix), "x", ncol(counts_matrix), ").\n")

# Verify metadata
cat("\nVerifying metadata columns...\n")
# Define required metadata columns, including the new annotations
required_meta <- c("sample", "condition", "genotype", "leiden_0.4", SUBCLUSTER_KEY,
                   "DG_majority_voting", "ISO_majority_voting")

# Add the manual annotation key only if we expect it to be there
if (USE_FINAL_ANNOTATED_FILE) {
  required_meta <- c(required_meta, MANUAL_ANNOTATION_KEY)
}

all_meta_cols <- names(seurat_obj@meta.data)
missing_meta <- required_meta[!required_meta %in% all_meta_cols]

if (length(missing_meta) > 0) {
    warning("Expected metadata columns MISSING: ", paste(missing_meta, collapse=", "))
} else {
    cat("All required metadata columns found.\n")
}

# Verify reductions
cat("\nAvailable reductions:", paste(names(seurat_obj@reductions), collapse=", "), "\n")
if (!"pca" %in% names(seurat_obj@reductions) || !"umap" %in% names(seurat_obj@reductions)) {
    warning("Expected reductions 'pca' or 'umap' may be missing.")
}

# %% [code]
# --- Save Seurat Object ---
cat("\nSaving Seurat object to:", output_seurat_path, "\n")
tryCatch({
    saveRDS(seurat_obj, file = output_seurat_path)
    cat("Seurat object saved successfully.\n")
}, error = function(e) {
    stop(paste("ERROR during saving Seurat object:", e$message))
})

cat("\nScript finished for cluster", CLUSTER_ID, ".\n")