#!/usr/bin/env Rscript
# %% [markdown]
# # 9. Convert Subclustered H5AD to Seurat Object (Parallelized)
#
# This script converts specific subclustered AnnData objects (.h5ad)
# for multiple cluster IDs in parallel into Seurat objects (.rds).
#
# %% [code]
# --- Configuration ---
USE_FINAL_ANNOTATED_FILE <- FALSE # <-- SET TO FALSE if you skipped the manual annotation step (script 12)
CLUSTER_IDS_TO_PROCESS <- c('0', '1', '16') # <-- IDs to process in parallel

# --- Load Libraries ---
library(schard)
library(Seurat)
library(Matrix)
library(parallel)

cat("Libraries loaded.\n")

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

# --- Main Processing Function ---
process_cluster <- function(CLUSTER_ID, use_final, proj_dir) {
  cat("\n--- Starting processing for Cluster ID:", CLUSTER_ID, "---\n")

  # --- Dynamic Keys ---
  SUBCLUSTER_KEY <- paste0('leiden_sub_', CLUSTER_ID)
  cat("Using subcluster metadata key:", SUBCLUSTER_KEY, "for cluster", CLUSTER_ID, "\n")
  MANUAL_ANNOTATION_KEY <- paste0('manual_annotation_subcluster_', CLUSTER_ID)

  # --- Path Configuration for the current cluster ---
  input_dir <- file.path(proj_dir, "combine_data", "results_subclustering")
  output_dir <- input_dir

  if (use_final) {
    input_h5ad_path <- file.path(input_dir, paste0("subcluster_", CLUSTER_ID, "_final_annotated.h5ad"))
    output_seurat_path <- file.path(output_dir, paste0("subcluster_", CLUSTER_ID, "_final_annotated.rds"))
    cat("Using FINAL annotated H5AD file for cluster", CLUSTER_ID, ".\n")
  } else {
    input_h5ad_path <- file.path(input_dir, paste0("subcluster_", CLUSTER_ID, "_cleaned.h5ad"))
    output_seurat_path <- file.path(output_dir, paste0("subcluster_", CLUSTER_ID, "_cleaned.rds"))
    cat("Using CLEANED H5AD file for cluster", CLUSTER_ID, ".\n")
  }
  cat("Input H5AD path:", input_h5ad_path, "\n")
  cat("Output Seurat path:", output_seurat_path, "\n")

  # --- Check Input File ---
  if (!file.exists(input_h5ad_path)) {
    stop(paste("ERROR: Input file not found for cluster", CLUSTER_ID, "at:", input_h5ad_path))
  }
  cat("Input file confirmed for cluster", CLUSTER_ID, ".\n")

  # --- Load h5ad into Seurat ---
  cat("\nLoading h5ad file for cluster", CLUSTER_ID, "...\n")
  seurat_obj <- NULL
  tryCatch({
    seurat_obj <- schard::h5ad2seurat(input_h5ad_path, use.raw = TRUE)
    cat("AnnData object loaded into Seurat for cluster", CLUSTER_ID, ".\n")
  }, error = function(e) {
    stop(paste("ERROR during h5ad loading for cluster", CLUSTER_ID, ":", e$message))
  })

  # --- Verification ---
  cat("\n--- Verifying Seurat Object for Cluster", CLUSTER_ID, "---\n")
  if (is.null(seurat_obj)) {
      stop(paste("ERROR: Seurat object is NULL after loading for cluster", CLUSTER_ID))
  }
  print(seurat_obj)

  default_assay <- DefaultAssay(seurat_obj)
  cat("\nDefault Assay:", default_assay, "\n")
  if (!"counts" %in% slotNames(seurat_obj@assays[[default_assay]])) {
      stop(paste("ERROR: 'counts' slot not found in the '", default_assay, "' assay for cluster", CLUSTER_ID))
  }
  counts_matrix <- GetAssayData(seurat_obj, assay = default_assay, layer = "counts")
  if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
      stop(paste("ERROR: 'counts' slot in assay '", default_assay, "' is empty for cluster", CLUSTER_ID))
  }
  cat("'counts' slot is populated.\n")

  required_meta <- c("sample", "condition", "genotype", "leiden_0.4", SUBCLUSTER_KEY,
                     "DG_majority_voting", "ISO_majority_voting")
  if (use_final) {
    required_meta <- c(required_meta, MANUAL_ANNOTATION_KEY)
  }
  all_meta_cols <- names(seurat_obj@meta.data)
  missing_meta <- required_meta[!required_meta %in% all_meta_cols]
  if (length(missing_meta) > 0) {
      warning(
          "For cluster ", CLUSTER_ID, ", expected metadata columns MISSING: ",
          paste(missing_meta, collapse=", "),
          "\n  Available columns: ", paste(all_meta_cols, collapse=", ")
      )
  } else {
      cat("All required metadata columns found for cluster", CLUSTER_ID, ".\n")
  }

  cat("\nAvailable reductions:", paste(names(seurat_obj@reductions), collapse=", "), "\n")
  # Flexible check for dimensional reductions
  find_reduction <- function(possible_names, available_reductions) {
    for (name in possible_names) {
      if (name %in% available_reductions) return(TRUE)
    }
    return(FALSE)
  }
  
  pca_found <- find_reduction(c("pca", "X_pca", "Xpca_"), names(seurat_obj@reductions))
  umap_found <- find_reduction(c("umap", "X_umap", "Xumap_"), names(seurat_obj@reductions))

  if (!pca_found) {
    warning("For cluster ", CLUSTER_ID, ", PCA reduction not found (checked for pca, X_pca, Xpca_).")
  }
  if (!umap_found) {
    warning("For cluster ", CLUSTER_ID, ", UMAP reduction not found (checked for umap, X_umap, Xumap_).")
  }

  # --- Save Seurat Object ---
  cat("\nSaving Seurat object for cluster", CLUSTER_ID, "to:", output_seurat_path, "\n")
  tryCatch({
      saveRDS(seurat_obj, file = output_seurat_path)
      cat("Seurat object for cluster", CLUSTER_ID, "saved successfully.\n")
  }, error = function(e) {
      stop(paste("ERROR during saving Seurat object for cluster", CLUSTER_ID, ":", e$message))
  })

  cat("--- Finished processing for Cluster ID:", CLUSTER_ID, "---\n")
  return(paste("Success:", CLUSTER_ID))
}

# --- Parallel Execution ---
num_cores <- detectCores()
if (is.na(num_cores) || num_cores <= 1) {
    num_cores <- 1
} else {
    num_cores <- num_cores - 1 # Leave one core free
}

cat("\nStarting script execution. Using", num_cores, "core(s) for processing.\n")

if (num_cores > 1 && .Platform$OS.type != "windows") {
    # Use mclapply for Linux, macOS, and WSL
    cat("Running in parallel using mclapply...\n")
    results <- mclapply(CLUSTER_IDS_TO_PROCESS, function(cid) {
        tryCatch({
            process_cluster(cid, USE_FINAL_ANNOTATED_FILE, project_dir)
        }, error = function(e) {
            cat("ERROR processing cluster", cid, ":", e$message, "\n")
            return(paste("Failure:", cid, "Error:", e$message))
        })
    }, mc.cores = num_cores)
} else if (num_cores > 1 && .Platform$OS.type == "windows") {
    # Use parLapply for Windows
    cat("Running in parallel using parLapply on Windows...\n")
    cl <- makeCluster(num_cores)
    
    # Export necessary objects and functions to the cluster workers
    clusterExport(cl, varlist=c("process_cluster", "USE_FINAL_ANNOTATED_FILE", "project_dir"))
    
    # Load required libraries on each worker
    clusterEvalQ(cl, {
        library(schard)
        library(Seurat)
        library(Matrix)
    })
    
    results <- parLapply(cl, CLUSTER_IDS_TO_PROCESS, function(cid) {
        tryCatch({
            process_cluster(cid, USE_FINAL_ANNOTATED_FILE, project_dir)
        }, error = function(e) {
            return(paste("Failure:", cid, "Error:", e$message))
        })
    })
    
    stopCluster(cl)
} else {
    # Fallback to sequential execution
    cat("Running sequentially...\n")
    results <- lapply(CLUSTER_IDS_TO_PROCESS, function(cid) {
        tryCatch({
            process_cluster(cid, USE_FINAL_ANNOTATED_FILE, project_dir)
        }, error = function(e) {
            cat("ERROR processing cluster", cid, ":", e$message, "\n")
            return(paste("Failure:", cid, "Error:", e$message))
        })
    })
}

cat("\n--- All processing finished. --- \n")
cat("Results:\n")
print(results)
cat("\nScript complete.\n")