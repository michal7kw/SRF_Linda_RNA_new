#!/usr/bin/env Rscript
# %% [markdown]
# Convert Merged *Annotated* H5AD to Seurat Object

# %% [code]
# Load libraries
library(schard)
library(Seurat)
library(Matrix)

cat("Libraries loaded.\n")

LEIDEN_KEY <- 'leiden_0.4' 


# %% [code]
# --- Configuration ---
# Set working directory (adjust if necessary, assumes running from project root)
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

# Define input and output directories relative to the project directory
input_dir <- file.path(project_dir, "combine_data", "results_from_raw")
output_dir <- input_dir

input_h5ad_path <- file.path(input_dir, "mapmycells.h5ad")
output_seurat_path <- file.path(output_dir, "mapmycells.rds")

cat("Configuration set.\n")
cat("Input H5AD path:", input_h5ad_path, "\n")
cat("Output Seurat path:", output_seurat_path, "\n")


# %% [code]
# --- Check Input File ---
cat("\nChecking for input file at:", input_h5ad_path, "\n")
if (!file.exists(input_h5ad_path)) {
  alt_path <- if (is_wsl) {
    gsub("^/mnt/d/", "D:/", input_h5ad_path) # Convert WSL to Windows
  } else {
    gsub("^D:/", "/mnt/d/", input_h5ad_path) # Convert Windows to WSL
  }
  cat("File not found. Trying alternative path:", alt_path, "\n")

  if (file.exists(alt_path)) {
    cat("File found at alternative path. Using:", alt_path, "\n")
    input_h5ad_path <- alt_path # Update the path to the one that works
  } else {
    stop(paste("ERROR: Input file not found at either path:",
               input_h5ad_path, "or", alt_path))
  }
}
cat("Input file confirmed at:", input_h5ad_path, "\n")

# %% [code]
# --- Load h5ad into Seurat using schard ---
# Load h5ad into Seurat using schard, loading the main .X matrix
cat("\nLoading h5ad file:", input_h5ad_path, "\n")
cat("Using schard::h5ad2seurat with use.raw = TRUE to load raw counts directly from .raw slot.\n")

# Initialize seurat_obj to NULL
seurat_obj <- NULL

# Force check file existence just before loading
cat("\nPerforming explicit file existence check just before loading...\n")
stopifnot(file.exists(input_h5ad_path))
cat("Explicit file check PASSED.\n")

tryCatch({
  # Load directly using use.raw = TRUE. This should populate the 'counts' slot
  # from adata.raw.X and potentially the 'data' slot from adata.X if present.
  # NOTE: If adata.X contains normalized/scaled data, it might overwrite the 'data' slot.
  # The 'counts' slot should reliably contain the raw counts from adata.raw.X.
  # Load using use.raw = FALSE. This should load adata.X into the 'data' slot.
  # The Python script aims to put the definitive raw counts (including transgene) in .raw.
  seurat_obj <- schard::h5ad2seurat(input_h5ad_path, use.raw = TRUE)

  cat("\nAnnData object loaded into Seurat object.\n")

  # --- Basic check after loading ---
  if (is.null(seurat_obj)) {
      stop("ERROR: Seurat object is NULL after loading attempt.")
  }

  default_assay <- DefaultAssay(seurat_obj)
  if (is.null(default_assay) || !default_assay %in% names(seurat_obj@assays)) {
      # Try to guess the assay name if default is missing (e.g., 'RNA')
      if ("RNA" %in% names(seurat_obj@assays)) {
          default_assay <- "RNA"
          DefaultAssay(seurat_obj) <- default_assay
          cat("Default assay was missing, set to 'RNA'.\n")
      } else {
          stop(paste("Default assay '", default_assay, "' not found or invalid after loading, and 'RNA' assay not found either.", sep=""))
      }
  }
  cat(paste("Working with assay:", default_assay, "\n"))

  # Verify that the 'counts' slot was populated (should come from .raw.X)
  assay_obj <- seurat_obj@assays[[default_assay]]
  if (!"counts" %in% slotNames(assay_obj)) {
      stop(paste("ERROR: 'counts' slot not found in the '", default_assay, "' assay after loading with use.raw = TRUE. Check input AnnData .raw slot.", sep=""))
  }
  counts_matrix <- GetAssayData(seurat_obj, assay = default_assay, layer = "counts")
  if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
      stop(paste("ERROR: 'counts' slot in assay '", default_assay, "' is empty after loading.", sep=""))
  }
  cat("'counts' slot populated from adata.raw.X (Dimensions: ", nrow(counts_matrix), "x", ncol(counts_matrix), ").\n")

  # Optional: Check if 'data' slot was also populated (from adata.X)
  # This might happen if .raw didn't exist and .X was used by schard as a fallback,
  # or if schard populates both.
  if ("data" %in% slotNames(assay_obj)) {
      data_matrix_check <- GetAssayData(seurat_obj, assay = default_assay, layer = "data")
      if (!is.null(data_matrix_check) && nrow(data_matrix_check) > 0 && ncol(data_matrix_check) > 0) {
          cat("'data' slot was also populated (likely from adata.X, Dimensions: ", nrow(data_matrix_check), "x", ncol(data_matrix_check), ").\n")
      } else {
          cat("'data' slot is present but empty or invalid.\n")
      }
  } else {
      cat("'data' slot was not populated.\n")
  }

}, error = function(e) {
  stop(paste("ERROR during h5ad loading:", e$message))
})

# %% [code]
# --- Verification ---
cat("\n--- Seurat Object Summary (Post-Loading) ---\n")
if (!is.null(seurat_obj)) {
    print(seurat_obj)

    # Check assays
    cat("\nAvailable Assays:", names(seurat_obj@assays), "\n")
    default_assay <- DefaultAssay(seurat_obj) # Re-confirm default assay

    if (!is.null(default_assay) && default_assay %in% names(seurat_obj@assays)) {
        cat(paste("\nVerifying", default_assay, "Assay:\n"))
        assay_obj_verify <- seurat_obj@assays[[default_assay]]

        # Check if counts slot exists and looks like raw counts (integers, max > 1)
        if (!"counts" %in% slotNames(assay_obj_verify)) {
             stop(paste("ERROR: Verification failed - 'counts' slot not found in the '", default_assay, "' assay.", sep=""))
        }
        counts_matrix_final <- GetAssayData(seurat_obj, assay = default_assay, layer = "counts")
        cat("'counts' slot dimensions:", dim(counts_matrix_final)[1], "x", dim(counts_matrix_final)[2], "\n")

        if (nrow(counts_matrix_final) > 0 && ncol(counts_matrix_final) > 0) {
            max_count <- max(counts_matrix_final)
            # Sample check for integers (more robust check)
            # Check a sample for non-integer values, as checking all can be slow
            sample_indices <- sample(length(counts_matrix_final), min(1000, length(counts_matrix_final)))
            is_integer_final <- all(counts_matrix_final[sample_indices] == floor(counts_matrix_final[sample_indices]))
            cat("Max value in 'counts' slot:", max_count, "\n")
            cat("Sample check for integers in 'counts' slot:", is_integer_final, "\n")
            if (!is_integer_final || max_count <= 1) {
                warning(paste("WARNING: Verification - 'counts' slot in assay '", default_assay, "' might not contain raw counts (non-integer or max <= 1). Please double-check.", sep=""))
                # Consider stopping if raw counts are strictly required:
                # stop(...)
            } else {
                 cat("'counts' slot appears to contain raw counts.\n")
            }
        } else {
            stop(paste("ERROR: Verification failed - 'counts' slot in assay '", default_assay, "' is empty.", sep=""))
        }

        # Check data slot presence
        if ("data" %in% slotNames(assay_obj_verify)) {
             data_matrix_final <- GetAssayData(seurat_obj, assay = default_assay, layer = "data")
             if (!is.null(data_matrix_final) && nrow(data_matrix_final) > 0 && ncol(data_matrix_final) > 0) {
                cat("'data' slot dimensions:", dim(data_matrix_final)[1], "x", dim(data_matrix_final)[2], "\n")
             } else {
                cat("'data' slot is present but empty or invalid.\n")
             }
        } else {
             cat("'data' slot not found in assay (might have been cleared or not loaded).\n")
        }
    } else {
        stop(paste("ERROR: Default assay '", default_assay, "' not found or invalid during verification.", sep=""))
    }

    # Check metadata - should include 'sample', 'condition', 'genotype', and CellTypist annotations
    cat("\nAvailable metadata columns (head):", paste(names(seurat_obj@meta.data), collapse=", "), "\n")
    # print(head(seurat_obj@meta.data))
    # Check for expected base columns, CellTypist columns, AND the alternative Leiden clustering
    # Updated for simplified AnnData object
    required_meta <- c("sample", "condition", "genotype",
                       "n_genes_by_counts", "total_counts", "pct_counts_mt", # Basic QC metrics
                       "DG_majority_voting", "ISO_majority_voting",
                       "DG_conf_score", "ISO_conf_score",
                       LEIDEN_KEY, "mapmycells_first_layer", "mapmycells_second_layer") # Clustering and annotation
                       # "Rosa26_SBP1_expression" is now a feature, not metadata
    all_meta_cols <- names(seurat_obj@meta.data)
    missing_meta <- required_meta[!required_meta %in% all_meta_cols]
    found_meta <- required_meta[required_meta %in% all_meta_cols]

    if (length(found_meta) > 0) {
        cat("Found expected metadata columns:", paste(found_meta, collapse=", "), "\n")
    }
    if (length(missing_meta) > 0) {
        warning("Expected metadata columns MISSING: ", paste(missing_meta, collapse=", "), "\nCheck the input H5AD file and the schard conversion.")
    }


    # Check reductions (like UMAP, PCA) - names might be prefixed with 'X' by schard
    cat("\nAvailable reductions:", paste(names(seurat_obj@reductions), collapse=", "), "\n")

    # Function to find reduction name flexibly
    find_reduction <- function(possible_names, reductions_list) {
      for (name in possible_names) {
        if (name %in% reductions_list) {
          return(name)
        }
      }
      return(NULL)
    }
}

# %% [code]

# --- Save Seurat Object ---
if (!is.null(seurat_obj)) {
    # --- Final Check for Transgene Feature ---
    cat("\nPerforming final check for 'Rosa26_SBP1' (expected as 'Rosa26-SBP1' in Seurat) feature before saving...\n")
    final_counts_matrix <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "counts")
    # R/Seurat often converts underscores to hyphens in feature names
    expected_seurat_feature_name <- "Rosa26-SBP1"
    if (expected_seurat_feature_name %in% rownames(final_counts_matrix)) {
        cat(paste("SUCCESS: Feature (expected as '", expected_seurat_feature_name, "' in Seurat, originally 'Rosa26_SBP1') found in the feature names (rownames) of the counts matrix.\n", sep=""))
    } else {
        warning(paste("WARNING: Feature 'Rosa26_SBP1' (expected as '", expected_seurat_feature_name, "' in Seurat) NOT found in the feature names (rownames) of the final counts matrix before saving. Check H5AD and conversion steps.", sep=""))
        # Optionally print some rownames to debug:
        cat("First 10 rownames:", paste(head(rownames(final_counts_matrix), 10), collapse=", "), "\n")
        cat("Last 10 rownames:", paste(tail(rownames(final_counts_matrix), 10), collapse=", "), "\n")
    }

    cat("\nSaving Seurat object to:", output_seurat_path, "\n")
    tryCatch({
        saveRDS(seurat_obj, file = output_seurat_path)
        cat("Seurat object saved successfully.\n")
    }, error = function(e) {
        stop(paste("ERROR during saving Seurat object:", e$message))
    })
} else {
    cat("\nSeurat object is NULL, skipping save.\n")
}