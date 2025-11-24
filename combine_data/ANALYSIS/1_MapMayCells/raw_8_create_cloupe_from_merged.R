#!/usr/bin/env Rscript
# Script to convert a Seurat object (.rds) derived from merged *annotated* data
# to a Loupe file (.cloupe) using the loupeR package.

# --- Dependencies & Setup ---

# Load libraries
library(loupeR)
library(Seurat)
library(utils) # For read.delim


# --- Configuration ---
# Detect if running in WSL or Windows and set paths accordingly
is_wsl <- FALSE
try({
  if (Sys.info()["sysname"] == "Linux" && file.exists("/proc/version")) {
    proc_version <- readLines("/proc/version", n = 1)
    is_wsl <- grepl("microsoft", tolower(proc_version))
  }
}, silent = TRUE)

if (is_wsl) {
  # WSL path format
  project_dir <- "/mnt/d/Github/SRF_Linda_RNA"
  cat("Running in WSL environment, using path:", project_dir, "\n")
} else {
  # Windows path format
  project_dir <- "D:/Github/SRF_Linda_RNA"
  cat("Running in Windows environment, using path:", project_dir, "\n")
}

# Input Seurat object (created by raw_6_convert_merged_to_seurat.R from simplified annotated H5AD)
seurat_input_dir <- file.path(project_dir, "combine_data", "results_from_raw")
seurat_input_file <- file.path(seurat_input_dir, "mapmycells.rds")

output_dir <- seurat_input_dir
output_name <- "mapmycells_loupe"


# --- Check Input Files ---
cat("Checking for input Seurat file at:", seurat_input_file, "\n")
if (!file.exists(seurat_input_file)) {
  # Try alternative path format as fallback
  alt_path <- if (is_wsl) {
    gsub("^/mnt/d/", "D:/", seurat_input_file)
  } else {
    gsub("^D:/", "/mnt/d/", seurat_input_file)
  }
  cat("File not found. Trying alternative path:", alt_path, "\n")

  if (file.exists(alt_path)) {
    cat("File found at alternative path. Using:", alt_path, "\n")
    seurat_input_file <- alt_path
  } else {
    stop(paste("ERROR: Input Seurat file not found at either path:",
               seurat_input_file, "or", alt_path))
  }
}
cat("Input Seurat file confirmed at:", seurat_input_file, "\n")

# NOTE: We are NOT using an external feature file here.
# loupeR should use the rownames from the Seurat object's count matrix.
feature_ids <- NULL
cat("Using feature names (rownames) directly from the Seurat object.\n")


# --- Load Seurat Object ---
cat("\nLoading Seurat object from:", seurat_input_file, "\n")
tryCatch({
  seurat_obj <- readRDS(seurat_input_file)
  cat("Seurat object loaded successfully.\n")
  print(seurat_obj)

  # --- Verify Seurat Object Structure ---
  # Check 1: Ensure the default assay (usually RNA) has raw counts
  default_assay <- DefaultAssay(seurat_obj)
  cat(paste("Checking assay:", default_assay, "\n"))
  # Corrected check: Verify "counts" is a slot within the Assay object
  if (!default_assay %in% names(seurat_obj@assays)){
      stop(paste("ERROR: Default assay '", default_assay, "' not found in Seurat object.", sep=""))
  }
  if (!"counts" %in% slotNames(seurat_obj@assays[[default_assay]])) {
      stop(paste("ERROR: 'counts' slot not found in the default assay (", default_assay, "). loupeR requires raw counts."))
  }
  # Check if the counts slot is empty
  counts_data <- GetAssayData(seurat_obj, assay = default_assay, slot = "counts")
  if (nrow(counts_data) == 0 || ncol(counts_data) == 0) {
    stop(paste("ERROR: 'counts' slot in the default assay (", default_assay, ") is present but empty."))
  }
  cat("'counts' slot found and contains data in the default assay.\n")
  cat("Number of features in Seurat object:", nrow(counts_data), "\n")

  # --- Correct Barcode Formatting for LoupeR ---
  # Re-introducing correction: LoupeR requires specific barcode formats.
  # scanpy.concat adds suffixes like '-0', '-1' which need adjustment.
  # We will replace the hyphen preceding the final numeric suffix with an underscore.
  # Example: AAACAGCCAAGCTTTG-1-0  ->  AAACAGCCAAGCTTTG-1_0
  cat("\nCorrecting cell barcodes for LoupeR compatibility...\n")
  original_barcodes <- colnames(seurat_obj)

  # Regex: Capture everything up to the last hyphen (.*) and the final digits (\d+)
  # Replace the last hyphen with an underscore.
  corrected_barcodes <- gsub("(.*)-(\\d+)$", "\\1_\\2", original_barcodes)

  num_changed <- sum(original_barcodes != corrected_barcodes)
  if (num_changed > 0) {
    cat(paste("Corrected", num_changed, "barcodes by replacing final hyphen with underscore.\n"))
    # Rename cells in the Seurat object
    seurat_obj <- RenameCells(seurat_obj, new.names = corrected_barcodes)
    cat("Cell names updated in the Seurat object.\n")
    # Optional: Print an example
    example_index <- which(original_barcodes != corrected_barcodes)[1]
    if (!is.na(example_index)){
        old_example <- original_barcodes[example_index]
        new_example <- corrected_barcodes[example_index]
        cat(paste("  Example: '", old_example, "' -> '", new_example, "'\n", sep=""))
    }
  } else {
    # This case might occur if scanpy did not add suffixes (e.g., only one sample merged)
    cat("No barcodes matched the pattern for correction (e.g., name-digits$). Checking raw format compliance.\n")
    # Add a check for standard Cell Ranger format if no changes were made
    if (!all(grepl("^[ACGT]+-[0-9]+$", original_barcodes))) {
        warning("Barcodes did not require suffix correction, but some do not match standard Cell Ranger format (e.g., ACGTACGT-1). LoupeR might still fail.")
    }
  }

}, error = function(e) {
  stop(paste("ERROR loading or processing Seurat object:", e$message))
})

# --- Debugging: Inspect Feature Names ---
if (!is.null(seurat_obj)) {
    cat("\n--- Debugging Feature Names in Seurat Object Before Loupe Conversion ---\n")
    default_assay_debug <- DefaultAssay(seurat_obj)
    if (!is.null(default_assay_debug) && default_assay_debug %in% names(seurat_obj@assays)) {
        counts_matrix_debug <- GetAssayData(seurat_obj, assay = default_assay_debug, slot = "counts")
        all_feature_names <- rownames(counts_matrix_debug)
        cat(paste("Total features in counts matrix:", length(all_feature_names), "\n"))
        
        # Check for expected names
        expected_name_hyphen <- "Rosa26-SBP1"
        expected_name_underscore <- "Rosa26_SBP1"
        
        if (expected_name_hyphen %in% all_feature_names) {
            cat(paste("SUCCESS (Debug): Found '", expected_name_hyphen, "' in feature names.\n", sep=""))
        } else {
            cat(paste("WARNING (Debug): Did NOT find '", expected_name_hyphen, "' in feature names.\n", sep=""))
        }
        
        if (expected_name_underscore %in% all_feature_names) {
            cat(paste("SUCCESS (Debug): Found '", expected_name_underscore, "' in feature names.\n", sep=""))
        } else {
            cat(paste("WARNING (Debug): Did NOT find '", expected_name_underscore, "' in feature names.\n", sep=""))
        }
        
        # Print a few feature names that contain "Rosa" or "SBP"
        rosa_features <- all_feature_names[grepl("rosa", all_feature_names, ignore.case = TRUE)]
        sbp_features <- all_feature_names[grepl("sbp", all_feature_names, ignore.case = TRUE)]
        
        if(length(rosa_features) > 0) {
            cat("Features containing 'Rosa':", paste(rosa_features, collapse=", "), "\n")
        } else {
            cat("No features found containing 'Rosa'.\n")
        }
        if(length(sbp_features) > 0) {
            cat("Features containing 'SBP':", paste(sbp_features, collapse=", "), "\n")
        } else {
            cat("No features found containing 'SBP'.\n")
        }
        cat("Last 10 feature names:", paste(tail(all_feature_names, 10), collapse=", "), "\n")

        # Print summary stats for Rosa26-SBP1 if found
        if (expected_name_hyphen %in% all_feature_names) {
            rosa_expression_values <- counts_matrix_debug[expected_name_hyphen, ]
            cat(paste("Summary for '", expected_name_hyphen, "' expression in counts matrix:\n", sep=""))
            print(summary(rosa_expression_values))
            cat(paste("Number of cells with >0 expression for '", expected_name_hyphen, "': ", sum(rosa_expression_values > 0), "\n", sep=""))
        }

    } else {
        cat("WARNING (Debug): Default assay not found or invalid for debugging.\n")
    }
    cat("---------------------------------------------------------------------\n")
} else {
    cat("WARNING (Debug): Seurat object is NULL, skipping feature name debugging.\n")
}

# --- Convert Seurat to Loupe ---
cat("\nStarting conversion to .cloupe file...\n")
cat("Output directory:", output_dir, "\n")
cat("Output name:", output_name, "\n")

tryCatch({
  create_loupe_from_seurat(
    obj = seurat_obj,
    output_dir = output_dir,
    output_name = output_name,
    feature_ids = feature_ids, # Use NULL to rely on Seurat object rownames
    force = TRUE               # Overwrite existing file
  )
  cat("\nConversion successful! Output file:", output_dir, "\n")

}, error = function(e) {
  # Provide more context if the error is about the executable path
  if (grepl("Error creating loupe file.*system command", e$message, ignore.case = TRUE)) {
      warning("Loupe conversion failed. This might be due to the loupeR executable not being found or the EULA not being accepted. Please run setup() interactively in R.", immediate. = TRUE)
  }
  stop(paste("ERROR during Loupe conversion:", e$message))
})

cat("\nScript finished.\n")