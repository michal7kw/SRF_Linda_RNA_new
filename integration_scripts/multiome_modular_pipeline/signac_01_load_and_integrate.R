#!/usr/bin/env Rscript

################################################################################
# Signac Multiome Analysis - Step 1: Data Loading and Integration
################################################################################
#
# This script loads scRNA-seq and scATAC-seq data and creates an integrated
# Seurat object following Signac best practices.
#
# Inputs:
#   - ATAC peak-barcode matrices:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_ATAC_data/chromap_final_output/{sample}_filtered_peak_bc_matrix_pseudobulk/
#       * matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz
#   - ATAC fragment files:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_ATAC_data/chromap_final_output/fragments/{sample}_fragments.tsv.gz
#   - RNA h5ad file:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/annotation_final.h5ad
#   - Cell Ranger ARC barcode whitelists:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/DATA/Data_integration/arc_barcode_lists/gex_737K-arc-v1.txt
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_ATAC_data/barcode_whitelists/atac_737K-arc-v1_rc.txt
#
# Outputs:
#   - Integrated Seurat object:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/integrated_seurat_raw.rds
#   - QC plots:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/qc/01_atac_qc_overview.pdf
#   - Integration summary:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/integration_summary.csv
#   - Translated fragment files (if bgzip/tabix available):
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/translated_fragments/{sample}_fragments_translated.tsv.gz
#   - Cache files (for faster reruns):
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/cache/integrated_seurat_object.rds
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/cache/consensus_peaks.rds
#
# Based on: https://stuartlab.org/signac/articles/pbmc_multiomic
#
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(EnsDb.Mmusculus.v79)  # Mouse gene annotations
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(GenomicRanges)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(reticulate)  # For reading h5ad files
  library(data.table)  # For fast fragment file processing
})

cat("================================================================================\n")
cat("  SIGNAC MULTIOME ANALYSIS - STEP 1: DATA LOADING AND INTEGRATION\n")
cat("================================================================================\n\n")

################################################################################
# Configure Python/Conda for reticulate
################################################################################

# Set conda path explicitly (helps reticulate find conda)
# Look for conda in common locations
possible_conda_paths <- c(
  "~/conda/bin/conda",
  "~/miniconda3/bin/conda",
  "~/anaconda3/bin/conda",
  "/beegfs/scratch/ric.sessa/kubacki.michal/conda/bin/conda",
  Sys.which("conda")
)

conda_found <- FALSE
for (conda_path in possible_conda_paths) {
  if (file.exists(conda_path)) {
    cat(sprintf("Found conda at: %s\n", conda_path))
    Sys.setenv(RETICULATE_CONDA = conda_path)
    use_condaenv("seurat_full2", conda = conda_path, required = FALSE)
    conda_found <- TRUE
    break
  }
}

if (!conda_found) {
  cat("Warning: Could not find conda installation\n")
  cat("Reticulate will try to auto-detect Python\n")
}

cat("\n")

################################################################################
# Configuration
################################################################################

# Directories
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top"
ATAC_DIR <- file.path(BASE_DIR, "SRF_ATAC_data/chromap_final_output")
RNA_DIR <- file.path(BASE_DIR, "SRF_Linda_RNA/combine_data/results_from_raw")
OUTPUT_DIR <- file.path(BASE_DIR, "SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results")

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "qc"), recursive = TRUE, showWarnings = FALSE)

# ALL samples are in single h5ad file
RNA_H5AD_FILE <- "annotation_final.h5ad"  # Contains ALL samples (Nestin + Emx1)

# ATAC samples mapping (all 4 samples)
ATAC_SAMPLES <- list(
  nestin_ctrl = "R26-Nestin-Ctrl-adult",
  nestin_mut = "R26-Nestin-Mut-adult",
  emx1_ctrl = "R26-Emx1-Ctrl-adult",
  emx1_mut = "R26-Emx1-Mut-adult"
)

# Use pseudobulk QC mode (more cells for better power)
USE_PSEUDOBULK <- TRUE

# QC thresholds (will be refined based on data distribution)
QC_THRESHOLDS <- list(
  min_atac_counts = 1000,
  max_atac_counts = 100000,
  min_rna_counts = 500,
  max_rna_counts = 50000,
  min_tss_enrichment = 2,
  max_nucleosome_signal = 2,
  max_blacklist_ratio = 0.05,
  min_pct_reads_in_peaks = 20
)

cat("Configuration:\n")
cat(sprintf("  Base directory: %s\n", BASE_DIR))
cat(sprintf("  RNA h5ad file: %s\n", RNA_H5AD_FILE))
cat(sprintf("  ATAC samples: %d\n", length(ATAC_SAMPLES)))
cat(sprintf("  Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("  Pseudobulk mode: %s\n", USE_PSEUDOBULK))
cat("\n")

################################################################################
# Function: Load ATAC Data
################################################################################

load_atac_data <- function(sample_name, use_pseudobulk = FALSE) {
  cat(sprintf("\n=== Loading ATAC data for %s ===\n", sample_name))

  # Determine matrix directory
  if (use_pseudobulk) {
    matrix_dir <- file.path(ATAC_DIR, paste0(sample_name, "_filtered_peak_bc_matrix_pseudobulk"))
    cat("  Using PSEUDOBULK QC mode (more cells for aggregation)\n")
  } else {
    matrix_dir <- file.path(ATAC_DIR, paste0(sample_name, "_filtered_peak_bc_matrix_qc"))
    cat("  Using STRINGENT QC mode\n")
  }

  if (!dir.exists(matrix_dir)) {
    stop(sprintf("Matrix directory not found: %s", matrix_dir))
  }

  # Load peak-barcode matrix
  cat("  Loading peak-barcode matrix...\n")
  counts <- ReadMtx(
    mtx = file.path(matrix_dir, "matrix.mtx.gz"),
    cells = file.path(matrix_dir, "barcodes.tsv.gz"),
    features = file.path(matrix_dir, "features.tsv.gz"),
    feature.column = 1  # Peak coordinates are in first column
  )

  # Fragment file path
  fragment_path <- file.path(ATAC_DIR, "fragments", paste0(sample_name, "_fragments.tsv.gz"))

  if (!file.exists(fragment_path)) {
    stop(sprintf("Fragment file not found: %s", fragment_path))
  }

  # Check for tabix index
  if (!file.exists(paste0(fragment_path, ".tbi"))) {
    stop(sprintf("Fragment file index not found: %s.tbi", fragment_path))
  }

  cat(sprintf("  Loaded matrix: %d peaks x %d cells\n", nrow(counts), ncol(counts)))
  cat(sprintf("  Fragment file: %s\n", fragment_path))

  return(list(
    counts = counts,
    fragment_path = fragment_path,
    sample = sample_name
  ))
}

################################################################################
# Function: Setup Python Environment for h5ad Reading
################################################################################

setup_python_env <- function() {
  cat("  Setting up Python environment for h5ad reading...\n")

  # Strategy: Try conda environments FIRST (where scanpy is installed)
  # Then fall back to system Python

  tryCatch({
    # Option 1: Try seurat_full2 conda environment (YOUR environment!)
    cat("  Trying conda environment: seurat_full2...\n")
    use_condaenv("seurat_full2", required = TRUE)
    sc <- import("scanpy", delay_load = TRUE)
    cat("  ✓ Using conda environment 'seurat_full2' with scanpy\n")
    return(TRUE)
  }, error = function(e1) {
    cat(sprintf("    Failed: %s\n", e1$message))
    tryCatch({
      # Option 2: Try base conda environment
      cat("  Trying conda environment: base...\n")
      use_condaenv("base", required = TRUE)
      sc <- import("scanpy", delay_load = TRUE)
      cat("  ✓ Using conda environment 'base' with scanpy\n")
      return(TRUE)
    }, error = function(e2) {
      cat(sprintf("    Failed: %s\n", e2$message))
      tryCatch({
        # Option 3: Auto-detect available conda environments
        cat("  Searching for conda environments with scanpy...\n")
        conda_envs <- conda_list()
        if (nrow(conda_envs) > 0) {
          for (i in 1:nrow(conda_envs)) {
            env_name <- conda_envs$name[i]
            cat(sprintf("    Trying: %s...\n", env_name))
            result <- tryCatch({
              use_condaenv(env_name, required = TRUE)
              sc <- import("scanpy", delay_load = TRUE)
              cat(sprintf("  ✓ Using conda environment '%s' with scanpy\n", env_name))
              return(TRUE)
            }, error = function(e) {
              NULL
            })
            if (!is.null(result) && result) {
              return(TRUE)
            }
          }
        }
        stop("No conda environment with scanpy found")
      }, error = function(e3) {
        cat(sprintf("    Failed: %s\n", e3$message))
        tryCatch({
          # Option 4: System Python (least likely)
          cat("  Trying system Python...\n")
          use_python("/usr/bin/python3", required = TRUE)
          sc <- import("scanpy", delay_load = TRUE)
          cat("  ✓ Using system Python with scanpy\n")
          return(TRUE)
        }, error = function(e4) {
          cat("\n  ✗ Could not find scanpy in any Python environment!\n\n")
          cat("  Environments tried:\n")
          cat("    1. conda: seurat_full2\n")
          cat("    2. conda: base\n")
          cat("    3. All available conda environments\n")
          cat("    4. System Python: /usr/bin/python3\n\n")
          cat("  Please install scanpy in your conda environment:\n")
          cat("    conda activate seurat_full2\n")
          cat("    pip install scanpy anndata\n\n")
          cat("  Then verify:\n")
          cat("    python -c 'import scanpy; print(scanpy.__version__)'\n")
          return(FALSE)
        })
      })
    })
  })
}

################################################################################
# Function: Load RNA Data from h5ad
################################################################################

load_rna_from_h5ad <- function(h5ad_path) {
  cat("\n=== Loading RNA data from h5ad ===\n")

  if (!file.exists(h5ad_path)) {
    stop(sprintf("RNA h5ad file not found: %s", h5ad_path))
  }

  cat(sprintf("  Loading: %s\n", h5ad_path))

  # Setup Python environment
  if (!setup_python_env()) {
    stop("Failed to setup Python environment. Please install scanpy:\n  pip install scanpy anndata")
  }

  # Use scanpy to read h5ad
  sc <- import("scanpy")
  sp <- import("scipy.sparse")
  adata <- sc$read_h5ad(h5ad_path)

  # Get cell names
  cell_names <- adata$obs_names$to_list()

  # CRITICAL DECISION: Use HVG-filtered gene set (BEST PRACTICE)
  # ────────────────────────────────────────────────────────────────────────
  # adata.raw.X = 33,685 genes (all genes) with raw counts
  # adata.X     = 26,870 genes (HVG-filtered) with log-normalized data
  #
  # STRATEGY: Use only 26,870 HVG genes for multiome integration
  #
  # RATIONALE:
  # 1. Non-HVG genes (6,815 excluded) are low-variance/housekeeping genes
  # 2. LinkPeaks relies on expression variance - non-HVGs give poor correlations
  # 3. Maintains consistency with scanpy pipeline (UMAP/clustering used HVGs)
  # 4. Standard practice in published multiome studies
  # 5. Better signal-to-noise ratio for peak-gene linkage
  #
  # IMPLEMENTATION:
  # - Extract raw counts for HVG genes from adata.raw.X (for RNA@counts)
  # - Use adata.X for log-normalized data (for RNA@data)
  # - Both slots will have exactly 26,870 genes
  # ────────────────────────────────────────────────────────────────────────

  if (!is.null(adata$raw)) {
    cat("  STRATEGY: Using HVG-filtered gene set (BEST PRACTICE)\n")
    cat("  ═══════════════════════════════════════════════════════\n")

    # Get HVG gene names (26,870 genes from adata.var)
    gene_names <- adata$var_names$to_list()
    cat(sprintf("  Target gene set: %d HVG genes\n", length(gene_names)))

    # Get all gene names from adata.raw (33,685 genes)
    all_genes <- adata$raw$var_names$to_list()
    cat(sprintf("  Available in adata.raw: %d genes\n", length(all_genes)))

    # Find which columns in adata.raw.X correspond to our HVGs
    gene_indices <- match(gene_names, all_genes)

    # Validate all HVGs are present in raw data
    if (any(is.na(gene_indices))) {
      n_missing <- sum(is.na(gene_indices))
      stop(sprintf("ERROR: %d HVG genes not found in adata.raw!", n_missing))
    }

    cat("  ✓ All HVG genes found in adata.raw\n")
    cat("  Extracting raw counts for HVG genes...\n")

    # Subset adata.raw.X to only HVG genes
    # Strategy: Pass both adata and gene names to Python, then subset

    # Assign to Python environment
    py$adata_obj <- adata
    py$hvg_genes <- gene_names

    # Use Python to subset adata.raw to HVG genes
    py_run_string("
# Find indices of HVG genes in adata_obj.raw.var_names
hvg_indices = [i for i, g in enumerate(adata_obj.raw.var_names) if g in hvg_genes]

# Extract the subset matrix (cells x genes)
# This works correctly with scipy sparse matrices
X_subset_hvg = adata_obj.raw.X[:, hvg_indices]
")

    # Get the subset matrix back to R
    X_data <- py$X_subset_hvg

    cat(sprintf("  ✓ Extracted raw counts for %d HVG genes\n", length(gene_names)))
    cat(sprintf("  Excluded %d non-HVG genes (low variance/housekeeping)\n",
                length(all_genes) - length(gene_names)))

  } else {
    cat("  WARNING: No adata.raw present, using adata.X\n")
    X_data <- adata$X
    gene_names <- adata$var_names$to_list()
  }

  cat(sprintf("\n  Final dimensions: %d cells × %d genes\n", length(cell_names), length(gene_names)))

  # Convert matrix with proper error handling
  cat("  Converting matrix to R format...\n")
  library(Matrix)

  counts <- tryCatch({
    # Try sparse conversion first (memory efficient)
    cat("    Attempting sparse matrix conversion...\n")
    X_csc <- sp$csc_matrix(X_data)

    data_vec <- as.numeric(X_csc$data)
    indices_vec <- as.integer(X_csc$indices) + 1L
    indptr_vec <- as.integer(X_csc$indptr)
    shape_tuple <- X_csc$shape
    n_cells <- as.integer(shape_tuple[1])
    n_genes <- as.integer(shape_tuple[2])

    cat(sprintf("    Sparse matrix: %d cells x %d genes (%d non-zero)\n",
                n_cells, n_genes, length(data_vec)))

    # Create sparse matrix (cells x genes)
    counts_cxg <- sparseMatrix(
      i = indices_vec,
      p = indptr_vec,
      x = data_vec,
      dims = c(n_cells, n_genes),
      index1 = TRUE
    )

    # Transpose to genes x cells and set names
    cat("    Transposing to genes x cells...\n")
    counts_gxc <- t(counts_cxg)
    rownames(counts_gxc) <- gene_names
    colnames(counts_gxc) <- cell_names

    cat("  ✓ Sparse matrix created (memory efficient)\n")
    counts_gxc

  }, error = function(e) {
    # Fall back to dense conversion
    cat(sprintf("    Sparse failed (%s), trying dense...\n", e$message))

    # Transpose FIRST, then set names
    X_mat <- as.matrix(X_data)
    cat(sprintf("    Matrix shape: %d x %d\n", nrow(X_mat), ncol(X_mat)))

    counts_t <- t(X_mat)
    cat(sprintf("    After transpose: %d x %d\n", nrow(counts_t), ncol(counts_t)))

    # Verify dimensions match
    if (nrow(counts_t) != length(gene_names)) {
      stop(sprintf("Row mismatch: matrix has %d rows but %d gene names",
                   nrow(counts_t), length(gene_names)))
    }
    if (ncol(counts_t) != length(cell_names)) {
      stop(sprintf("Col mismatch: matrix has %d cols but %d cell names",
                   ncol(counts_t), length(cell_names)))
    }

    # Set names
    rownames(counts_t) <- gene_names
    colnames(counts_t) <- cell_names

    cat("  ✓ Dense matrix created\n")
    counts_t
  })

  # Extract metadata
  metadata <- as.data.frame(adata$obs)
  rownames(metadata) <- cell_names

  cat(sprintf("  ✓ Loaded RNA counts: %d genes x %d cells\n", nrow(counts), ncol(counts)))
  cat(sprintf("  Metadata columns: %d\n", ncol(metadata)))

  # Check for cell type column
  celltype_cols <- grep("cell_type|celltype|annotation", colnames(metadata),
                        ignore.case = TRUE, value = TRUE)
  if (length(celltype_cols) > 0) {
    cat(sprintf("  Cell type columns found: %s\n", paste(celltype_cols, collapse = ", ")))
  }

  # Extract embeddings (UMAP, PCA, etc.) from adata.obsm
  cat("\n  Loading embeddings from adata.obsm...\n")
  embeddings <- list()

  # Get available embeddings - use Python method keys()
  obsm_keys <- tryCatch({
    # Call the keys() method on obsm and convert to list
    adata$obsm$keys()$to_list()
  }, error = function(e) {
    character(0)
  })

  if (length(obsm_keys) > 0) {
    cat(sprintf("    Available embeddings: %s\n", paste(obsm_keys, collapse = ", ")))

    for (key in obsm_keys) {
      tryCatch({
        # Access embedding data using the key
        emb_data <- as.matrix(adata$obsm[[key]])
        rownames(emb_data) <- cell_names

        # Convert key name (X_umap → umap, X_pca → pca)
        clean_key <- gsub("^X_", "", key)
        embeddings[[clean_key]] <- emb_data

        cat(sprintf("    ✓ Loaded %s: %d cells x %d dimensions\n",
                   key, nrow(emb_data), ncol(emb_data)))
      }, error = function(e) {
        cat(sprintf("    ⚠ Failed to load %s: %s\n", key, e$message))
      })
    }
  } else {
    cat("    No embeddings found in adata.obsm\n")
  }

  # Extract normalized/scaled data from adata.layers
  cat("\n  Loading additional data layers...\n")
  layers <- list()

  # Get available layers - use Python method keys()
  layer_keys <- tryCatch({
    # Call the keys() method on layers and convert to list
    adata$layers$keys()$to_list()
  }, error = function(e) {
    character(0)
  })

  if (length(layer_keys) > 0) {
    cat(sprintf("    Available layers: %s\n", paste(layer_keys, collapse = ", ")))

    # Commonly used layer names for normalized/scaled data
    important_layers <- c("normalized", "scaled", "log1p", "lognorm")

    for (key in layer_keys) {
      # Only load important layers to save memory
      if (any(grepl(paste(important_layers, collapse = "|"), key, ignore.case = TRUE))) {
        tryCatch({
          layer_data <- adata$layers[[key]]

          # Convert to sparse matrix
          if (sp$issparse(layer_data)) {
            layer_csr <- sp$csr_matrix(layer_data)
            indices_vec <- as.integer(layer_csr$indices)
            indptr_vec <- as.integer(layer_csr$indptr)
            data_vec <- as.numeric(layer_csr$data)

            layer_counts <- sparseMatrix(
              i = indices_vec + 1,  # R is 1-indexed
              p = indptr_vec,
              x = data_vec,
              dims = c(layer_csr$shape[1], layer_csr$shape[2]),
              index1 = TRUE,
              repr = "C"
            )
          } else {
            layer_counts <- as.matrix(layer_data)
          }

          # Transpose to genes x cells
          layer_counts <- t(layer_counts)
          rownames(layer_counts) <- gene_names
          colnames(layer_counts) <- cell_names

          layers[[key]] <- layer_counts

          cat(sprintf("    ✓ Loaded layer '%s': %d genes x %d cells\n",
                     key, nrow(layer_counts), ncol(layer_counts)))
        }, error = function(e) {
          cat(sprintf("    ⚠ Failed to load layer '%s': %s\n", key, e$message))
        })
      }
    }
  } else {
    cat("    No additional layers found\n")
  }

  # Extract adata.X (log-normalized HVG data) - MATCHES gene_names from above
  if (!is.null(adata$raw)) {
    cat("\n  Loading adata.X (log-normalized HVG data)...\n")
    cat("  ═══════════════════════════════════════════\n")
    tryCatch({
      X_data_norm <- adata$X

      # Get gene names - MUST match gene_names used for counts above
      X_gene_names <- adata$var_names$to_list()

      # Verify gene consistency
      if (!identical(X_gene_names, gene_names)) {
        warning("Gene name mismatch between counts and normalized data!")
        cat(sprintf("    Counts genes: %d\n", length(gene_names)))
        cat(sprintf("    Normalized genes: %d\n", length(X_gene_names)))
      } else {
        cat(sprintf("  ✓ Gene set consistency verified: %d HVG genes\n", length(X_gene_names)))
      }

      # Use same robust approach as for adata.raw.X subset
      X_counts <- tryCatch({
        # Try sparse conversion first
        if (sp$issparse(X_data_norm)) {
          X_csr <- sp$csr_matrix(X_data_norm)
          indices_vec <- as.integer(X_csr$indices)
          indptr_vec <- as.integer(X_csr$indptr)
          data_vec <- as.numeric(X_csr$data)

          X_sparse <- sparseMatrix(
            i = indices_vec + 1,
            p = indptr_vec,
            x = data_vec,
            dims = c(X_csr$shape[1], X_csr$shape[2]),
            index1 = TRUE,
            repr = "C"
          )

          # Transpose to genes x cells
          X_sparse_t <- t(X_sparse)
          X_sparse_t
        } else {
          # Dense matrix
          X_mat <- as.matrix(X_data_norm)
          t(X_mat)
        }
      }, error = function(e) {
        cat(sprintf("    Sparse failed (%s), trying dense...\n", e$message))
        X_mat <- as.matrix(X_data_norm)
        t(X_mat)
      })

      # Set dimnames
      rownames(X_counts) <- X_gene_names
      colnames(X_counts) <- cell_names

      # CRITICAL: This is log-normalized data - store it for RNA@data slot
      layers[["data"]] <- X_counts

      cat(sprintf("    ✓ Loaded adata.X (log-normalized): %d genes x %d cells\n",
                 nrow(X_counts), ncol(X_counts)))
      cat(sprintf("    Data range: [%.3f, %.3f] (should be log1p scale)\n",
                 min(X_counts@x), max(X_counts@x)))
    }, error = function(e) {
      cat(sprintf("    ⚠ Failed to load adata.X: %s\n", e$message))
      cat("    Will use NormalizeData on raw counts instead\n")
    })
  }

  cat("\n  ✓ RNA data loading complete\n")

  return(list(
    counts = counts,
    metadata = metadata,
    embeddings = embeddings,
    layers = layers
  ))
}

################################################################################
# Function: Load Cell Ranger ARC Barcode Translation
################################################################################

load_arc_barcode_translation <- function() {
  cat("\n=== Loading Cell Ranger ARC Barcode Translation ===\n")

  # Paths to Cell Ranger ARC whitelist files
  rna_whitelist_file <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/DATA/Data_integration/arc_barcode_lists/gex_737K-arc-v1.txt"
  atac_whitelist_file <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_ATAC_data/barcode_whitelists/atac_737K-arc-v1_rc.txt"

  if (!file.exists(rna_whitelist_file)) {
    stop(sprintf("RNA whitelist not found: %s", rna_whitelist_file))
  }

  if (!file.exists(atac_whitelist_file)) {
    stop(sprintf("ATAC whitelist not found: %s", atac_whitelist_file))
  }

  cat(sprintf("  Loading RNA whitelist: %s\n", basename(rna_whitelist_file)))
  cat(sprintf("  Loading ATAC whitelist: %s (reverse complement used by chromap)\n",
              basename(atac_whitelist_file)))

  cat("\n  ⚠️  CRITICAL ASSUMPTION:\n")
  cat("     ATAC data MUST have been processed with: atac_737K-arc-v1_rc.txt\n")
  cat("     If a different whitelist was used, barcode translation will be incorrect!\n")
  cat("     Verify with your chromap command or 10x pipeline settings.\n\n")

  # Load whitelist files
  rna_barcodes <- readLines(rna_whitelist_file)
  atac_barcodes <- readLines(atac_whitelist_file)

  if (length(rna_barcodes) != length(atac_barcodes)) {
    stop(sprintf("❌ Barcode lists have different lengths: RNA=%d, ATAC=%d\n",
                   length(rna_barcodes), length(atac_barcodes)),
         "   This indicates incompatible whitelist files!\n",
         "   Both must be from the same Cell Ranger ARC version.\n")
  }

  cat(sprintf("  ✓ Loaded %d RNA → %d ATAC barcode mappings\n",
              length(rna_barcodes), length(atac_barcodes)))

  # Create translation mappings (ATAC → RNA)
  atac_to_rna <- setNames(rna_barcodes, atac_barcodes)

  # Show examples
  cat("\n  Translation examples:\n")
  for (i in 1:min(5, length(rna_barcodes))) {
    cat(sprintf("    RNA: %s ↔ ATAC: %s\n", rna_barcodes[i], atac_barcodes[i]))
  }

  return(list(
    atac_to_rna = atac_to_rna,
    rna_whitelist = rna_barcodes,
    atac_whitelist = atac_barcodes
  ))
}

################################################################################
# Function: Translate ATAC Barcodes to RNA Format
################################################################################

translate_atac_barcodes <- function(atac_barcodes, barcode_mapping, sample_name) {
  cat(sprintf("\n  Translating ATAC barcodes for %s...\n", sample_name))

  # Sample suffix mapping (extracted from annotation_final.h5ad cell barcodes)
  # Each sample has a unique gem group ID in the format: BARCODE-1-X
  # These suffixes match the actual RNA cell barcodes in the h5ad file
  sample_suffix_map <- list(
    "R26-Emx1-Ctrl-adult" = "1-0",     # Emx1 Ctrl gem group
    "R26-Emx1-Mut-adult" = "1-1",      # Emx1 Mut gem group
    "R26-Nestin-Ctrl-adult" = "1-2",   # Nestin Ctrl gem group
    "R26-Nestin-Mut-adult" = "1-3"     # Nestin Mut gem group
  )

  # Lookup suffix for this sample
  suffix <- sample_suffix_map[[sample_name]]

  if (is.null(suffix)) {
    stop(sprintf("Unknown sample name: %s\nKnown samples: %s",
                 sample_name,
                 paste(names(sample_suffix_map), collapse = ", ")))
  }

  cat(sprintf("    Using gem group suffix: %s (from h5ad barcode analysis)\n", suffix))

  atac_to_rna <- barcode_mapping$atac_to_rna

  # VECTORIZED TRANSLATION (1000x+ faster than loop!)
  # Instead of looping through 70,000+ barcodes one-by-one,
  # use vectorized lookup which is orders of magnitude faster

  # Handle edge case: empty input
  if (length(atac_barcodes) == 0) {
    return(data.frame(
      original_atac_barcode = character(0),
      rna_full_barcode = character(0),
      translation_successful = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  # Vectorized lookup: returns NA for barcodes not in mapping
  rna_cores <- atac_to_rna[atac_barcodes]

  # Vectorized string concatenation: add suffix to all at once
  # paste0 is vectorized, so this processes all 70k barcodes in milliseconds
  rna_full_barcodes <- paste0(rna_cores, "-", suffix)

  # Vectorized logical test: which barcodes were successfully translated
  # (those that aren't NA after lookup)
  translation_successful <- !is.na(rna_cores)

  # Create results data frame (now takes <1 second instead of 2+ minutes)
  translation_results <- data.frame(
    original_atac_barcode = atac_barcodes,
    rna_full_barcode = rna_full_barcodes,
    translation_successful = translation_successful,
    stringsAsFactors = FALSE
  )

  # Set failed translations to NA (instead of "NA-1-2")
  translation_results$rna_full_barcode[!translation_successful] <- NA_character_

  # Calculate translation rate
  n_success <- sum(translation_successful)
  translation_rate <- n_success / length(atac_barcodes) * 100

  cat(sprintf("    Translation results:\n"))
  cat(sprintf("      Total ATAC barcodes: %d\n", length(atac_barcodes)))
  cat(sprintf("      Successfully translated: %d\n", n_success))
  cat(sprintf("      Translation rate: %.1f%%\n", translation_rate))

  return(translation_results)
}

################################################################################
# Function: Translate Fragment File Barcodes
################################################################################

translate_fragment_file <- function(fragment_file, atac_to_rna, sample_suffix,
                                     output_file, cells_to_keep = NULL) {
  cat(sprintf("\n  Translating fragment file: %s\n", basename(fragment_file)))
  cat(sprintf("    Sample suffix: %s\n", sample_suffix))

  if (!file.exists(fragment_file)) {
    stop(sprintf("Fragment file not found: %s", fragment_file))
  }

  # Check for required tools
  if (system("which bgzip", ignore.stdout = TRUE, ignore.stderr = TRUE) != 0) {
    stop("bgzip not found. Install with: conda install -c bioconda htslib")
  }

  if (system("which tabix", ignore.stdout = TRUE, ignore.stderr = TRUE) != 0) {
    stop("tabix not found. Install with: conda install -c bioconda htslib")
  }

  # Read fragment file using data.table for speed
  cat("    Reading fragment file...\n")

  # Fragment format: chr, start, end, barcode, count
  fragments <- data.table::fread(
    cmd = sprintf("zcat %s", fragment_file),
    header = FALSE,
    col.names = c("chr", "start", "end", "barcode", "count"),
    showProgress = FALSE
  )

  cat(sprintf("    Read %d fragments\n", nrow(fragments)))

  # Translate barcodes using vectorized lookup
  cat("    Translating barcodes...\n")

  # Get RNA core barcodes
  rna_cores <- atac_to_rna[fragments$barcode]

  # Create full RNA barcodes with suffix
  rna_full_barcodes <- paste0(rna_cores, "-", sample_suffix)

  # Track translation success
  translation_success <- !is.na(rna_cores)

  cat(sprintf("    Translated %d/%d fragments (%.1f%%)\n",
              sum(translation_success), nrow(fragments),
              100 * sum(translation_success) / nrow(fragments)))

  # Keep only successfully translated fragments
  fragments_translated <- fragments[translation_success, ]
  fragments_translated$barcode <- rna_full_barcodes[translation_success]

  # Filter to cells in Seurat object if provided
  if (!is.null(cells_to_keep)) {
    cat(sprintf("    Filtering to %d cells in analysis...\n", length(cells_to_keep)))

    fragments_translated <- fragments_translated[fragments_translated$barcode %in% cells_to_keep, ]

    cat(sprintf("    Retained %d fragments for %d cells\n",
                nrow(fragments_translated),
                length(unique(fragments_translated$barcode))))
  }

  # Sort by chromosome and position
  data.table::setorder(fragments_translated, chr, start)

  # Write output file
  cat(sprintf("    Writing translated fragments to: %s\n", output_file))

  data.table::fwrite(fragments_translated, file = output_file,
                    sep = "\t", col.names = FALSE, quote = FALSE)

  # Compress with bgzip
  cat("    Compressing with bgzip...\n")
  system(sprintf("bgzip -f %s", output_file))

  output_file_gz <- paste0(output_file, ".gz")

  # Create tabix index
  cat("    Creating tabix index...\n")
  system(sprintf("tabix -p bed %s", output_file_gz))

  cat(sprintf("  ✓ Translated fragment file created: %s\n", output_file_gz))

  return(output_file_gz)
}

################################################################################
# Function: Create Integrated Seurat Object with Barcode Translation
################################################################################

create_integrated_object <- function(rna_data, atac_data_list, cache_dir = "signac_results/cache") {
  cat("\n=== Creating integrated multiome object with barcode translation ===\n")

  # Create cache directory if it doesn't exist
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Check for cached integrated object
  cache_file <- file.path(cache_dir, "integrated_seurat_object.rds")
  if (file.exists(cache_file)) {
    cat("\n  ✓ Found cached integrated object, loading...\n")
    cat(sprintf("    Cache file: %s\n", cache_file))
    cat("    Set cache_dir=NULL or delete cache to recompute\n")
    seurat_obj <- readRDS(cache_file)
    cat(sprintf("    Loaded: %d cells, %d RNA features, %d ATAC features\n",
                ncol(seurat_obj),
                nrow(seurat_obj[["RNA"]]),
                nrow(seurat_obj[["ATAC"]])))
    return(seurat_obj)
  }

  # STEP 1: Load barcode translation mapping
  barcode_mapping <- load_arc_barcode_translation()

  # STEP 2: Create Seurat object from RNA data
  cat("\n  Creating Seurat object from RNA data...\n")
  cat(sprintf("    RNA assay: %d features x %d cells\n",
              nrow(rna_data$counts), ncol(rna_data$counts)))

  seurat_obj <- CreateSeuratObject(
    counts = rna_data$counts,
    assay = "RNA",
    meta.data = rna_data$metadata
  )

  # STEP 2b: Validate/create required metadata columns
  cat("\n  Validating metadata columns...\n")
  cat("  ═══════════════════════════════\n")

  # Check for required columns
  required_cols <- c("genotype", "condition", "cell_type_L2_new")
  existing_cols <- colnames(seurat_obj@meta.data)

  cat(sprintf("    Existing columns: %s\n", paste(head(existing_cols, 10), collapse = ", ")))
  if (length(existing_cols) > 10) {
    cat(sprintf("    ... and %d more\n", length(existing_cols) - 10))
  }

  # Validate/create genotype column
  if (!"genotype" %in% existing_cols) {
    cat("    ⚠ 'genotype' column missing - attempting to create from 'sample' column\n")

    if ("sample" %in% existing_cols) {
      seurat_obj$genotype <- ifelse(
        grepl("Nestin", seurat_obj$sample, ignore.case = TRUE), "Nestin",
        ifelse(grepl("Emx1", seurat_obj$sample, ignore.case = TRUE), "Emx1", NA)
      )
      cat(sprintf("    ✓ Created 'genotype' column from 'sample'\n"))
      cat(sprintf("      Genotypes: %s\n", paste(unique(seurat_obj$genotype), collapse = ", ")))
    } else {
      stop("ERROR: Cannot create 'genotype' column - 'sample' column not found!")
    }
  } else {
    cat(sprintf("    ✓ 'genotype' column exists: %s\n",
                paste(unique(seurat_obj$genotype), collapse = ", ")))
  }

  # Validate/create condition column
  if (!"condition" %in% existing_cols) {
    cat("    ⚠ 'condition' column missing - attempting to create from 'sample' column\n")

    if ("sample" %in% existing_cols) {
      seurat_obj$condition <- ifelse(
        grepl("Ctrl", seurat_obj$sample, ignore.case = TRUE), "Ctrl",
        ifelse(grepl("Mut", seurat_obj$sample, ignore.case = TRUE), "Mut", NA)
      )
      cat(sprintf("    ✓ Created 'condition' column from 'sample'\n"))
      cat(sprintf("      Conditions: %s\n", paste(unique(seurat_obj$condition), collapse = ", ")))
    } else {
      stop("ERROR: Cannot create 'condition' column - 'sample' column not found!")
    }
  } else {
    cat(sprintf("    ⚠ 'condition' column exists with values: %s\n",
                paste(unique(seurat_obj$condition), collapse = ", ")))
    cat("    Standardizing condition values to 'Ctrl' and 'Mut'...\n")

    # Standardize condition values (handles Control/Ctrl/control and Mutant/Mut/mutant)
    seurat_obj$condition <- ifelse(
      grepl("^Contr", seurat_obj$condition, ignore.case = TRUE), "Ctrl",
      ifelse(grepl("^Mut", seurat_obj$condition, ignore.case = TRUE), "Mut",
             seurat_obj$condition)
    )
    cat(sprintf("    ✓ Standardized 'condition' column: %s\n",
                paste(unique(seurat_obj$condition), collapse = ", ")))
  }

  # Validate cell type column
  if (!"cell_type_L2_new" %in% existing_cols) {
    cat("    ⚠ WARNING: 'cell_type_L2_new' column missing!\n")
    cat("      This will cause issues in downstream analysis.\n")
    cat("      Please ensure cell type annotations are present in the h5ad file.\n")
  } else {
    n_celltypes <- length(unique(seurat_obj$cell_type_L2_new))
    cat(sprintf("    ✓ 'cell_type_L2_new' column exists: %d cell types\n", n_celltypes))
  }

  # Verify no NAs in critical columns
  if (any(is.na(seurat_obj$genotype))) {
    n_na <- sum(is.na(seurat_obj$genotype))
    cat(sprintf("    ⚠ WARNING: %d cells have NA genotype!\n", n_na))
  }
  if (any(is.na(seurat_obj$condition))) {
    n_na <- sum(is.na(seurat_obj$condition))
    cat(sprintf("    ⚠ WARNING: %d cells have NA condition!\n", n_na))
  }

  cat("  ✓ Metadata validation complete\n")

  # Add pre-computed normalized data from scanpy (adata.X)
  if (!is.null(rna_data$layers) && length(rna_data$layers) > 0) {
    cat("\n  Adding scanpy-normalized RNA data...\n")
    cat("  ════════════════════════════════════\n")

    for (layer_name in names(rna_data$layers)) {
      layer_data <- rna_data$layers[[layer_name]]

      tryCatch({
        if (layer_name == "data") {
          # This is log-normalized data from scanpy (adata.X)
          # Both counts and data should have SAME gene set (26,870 HVGs)

          counts_genes <- rownames(seurat_obj[["RNA"]])
          norm_genes <- rownames(layer_data)

          # Verify perfect match (they should be identical since both use HVGs)
          if (identical(counts_genes, norm_genes)) {
            cat(sprintf("  ✓ Perfect gene match: %d genes\n", length(counts_genes)))

            # Add normalized data using Seurat v5 compatible method
            # Try SetAssayData first (works for both v4 and v5)
            seurat_obj <- SetAssayData(
              object = seurat_obj,
              slot = "data",
              new.data = layer_data,
              assay = "RNA"
            )

            cat(sprintf("  ✓ Added scanpy log-normalized data to RNA data slot\n"))

            # Validate data range (should be log1p scale)
            # Handle both sparse and dense matrices
            if (inherits(layer_data, "sparseMatrix")) {
              data_range <- range(layer_data@x)
            } else {
              data_range <- range(as.vector(layer_data))
            }
            cat(sprintf("  Data range: [%.3f, %.3f] (log1p scale)\n",
                       data_range[1], data_range[2]))

            if (data_range[2] > 20) {
              warning("Unexpected data range - may not be log-normalized!")
            }

          } else {
            # Fallback: use common genes (shouldn't happen with HVG approach)
            common_genes <- intersect(counts_genes, norm_genes)
            cat(sprintf("  ⚠ Gene mismatch - using %d common genes\n", length(common_genes)))
            seurat_obj <- SetAssayData(
              object = seurat_obj,
              slot = "data",
              new.data = layer_data[common_genes, colnames(seurat_obj)],
              assay = "RNA"
            )
          }

        } else {
          # Store other layers for reference
          cat(sprintf("  Loaded layer '%s' (%d genes)\n", layer_name, nrow(layer_data)))
        }
      }, error = function(e) {
        cat(sprintf("  ⚠ Failed to add layer '%s': %s\n", layer_name, e$message))
      })
    }
  }

  # Add embeddings (UMAP, PCA, etc.)
  if (!is.null(rna_data$embeddings) && length(rna_data$embeddings) > 0) {
    cat("\n  Adding RNA embeddings...\n")

    for (emb_name in names(rna_data$embeddings)) {
      emb_data <- rna_data$embeddings[[emb_name]]

      tryCatch({
        # Create Seurat dimensional reduction object
        colnames(emb_data) <- paste0(emb_name, "_", seq_len(ncol(emb_data)))

        reduction_obj <- CreateDimReducObject(
          embeddings = emb_data,
          key = paste0(emb_name, "_"),
          assay = "RNA"
        )

        seurat_obj[[emb_name]] <- reduction_obj

        cat(sprintf("    ✓ Added %s embedding (%d dimensions)\n", emb_name, ncol(emb_data)))
      }, error = function(e) {
        cat(sprintf("    ⚠ Failed to add %s embedding: %s\n", emb_name, e$message))
      })
    }
  }

  # CRITICAL: Verify RNA normalization exists (required for LinkPeaks)
  cat("\n  ═══════════════════════════════════════════════\n")
  cat("  FINAL VERIFICATION: RNA Normalization Status\n")
  cat("  ═══════════════════════════════════════════════\n")

  rna_data_slot <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  rna_counts_slot <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

  # Helper function to get range from sparse or dense matrix
  safe_range <- function(mat) {
    if (inherits(mat, "sparseMatrix")) {
      if (length(mat@x) == 0) {
        return(c(Inf, -Inf))  # Empty
      }
      return(range(mat@x))
    } else {
      # Dense matrix
      return(range(as.vector(mat)))
    }
  }

  # Check 1: Data slot is not empty
  if (inherits(rna_data_slot, "sparseMatrix")) {
    data_is_empty <- (length(rna_data_slot@x) == 0)
  } else {
    data_is_empty <- (length(rna_data_slot) == 0) || all(dim(rna_data_slot) == 0)
  }

  # Check 2: Data slot differs from counts (normalized vs raw)
  data_differs_from_counts <- !identical(rna_data_slot, rna_counts_slot)

  # Check 3: Data values are in log scale (max < 20, typical for log1p)
  data_range <- safe_range(rna_data_slot)
  data_is_log_scale <- (data_range[2] < 20) && (data_range[2] > 0)

  cat(sprintf("  RNA @counts: %d genes × %d cells\n",
             nrow(rna_counts_slot), ncol(rna_counts_slot)))
  cat(sprintf("  RNA @data:   %d genes × %d cells\n",
             nrow(rna_data_slot), ncol(rna_data_slot)))
  cat(sprintf("  Data range: [%.3f, %.3f]\n", data_range[1], data_range[2]))

  # Comprehensive check
  has_valid_normalization <- !data_is_empty && data_differs_from_counts && data_is_log_scale

  if (has_valid_normalization) {
    cat("\n  ✅ RNA normalization VERIFIED\n")
    cat("     Source: scanpy log-normalization (adata.X)\n")
    cat("     Status: Ready for LinkPeaks correlation analysis\n")

  } else {
    cat("\n  ⚠️  RNA normalization issue detected:\n")
    if (data_is_empty) {
      cat("     - Data slot is empty\n")
    }
    if (!data_differs_from_counts) {
      cat("     - Data slot identical to counts (not normalized)\n")
    }
    if (!data_is_log_scale) {
      cat("     - Data not in log scale (range unexpected)\n")
    }

    cat("\n  → Running Seurat NormalizeData as fallback...\n")

    seurat_obj <- NormalizeData(seurat_obj, assay = "RNA",
                                normalization.method = "LogNormalize",
                                scale.factor = 10000, verbose = FALSE)

    # Re-verify
    data_range_new <- range(GetAssayData(seurat_obj, assay = "RNA", slot = "data")@x)
    cat(sprintf("  ✓ Re-normalized: range [%.3f, %.3f]\n",
               data_range_new[1], data_range_new[2]))
  }

  # Get RNA cell barcodes
  rna_barcodes <- colnames(seurat_obj)
  cat(sprintf("\n  ✓ RNA Seurat object ready: %d cells\n", length(rna_barcodes)))

  # STEP 3: Process ATAC data with barcode translation
  cat("\n  Processing ATAC data with barcode translation...\n")

  translated_atac_list <- list()

  for (i in seq_along(atac_data_list)) {
    atac_info <- atac_data_list[[i]]
    sample_name <- atac_info$sample

    cat(sprintf("\n  [%d/%d] Processing %s\n", i, length(atac_data_list), sample_name))

    # Get original ATAC barcodes
    atac_barcodes <- colnames(atac_info$counts)
    cat(sprintf("    Original ATAC barcodes: %d\n", length(atac_barcodes)))

    # Translate ATAC barcodes to RNA format
    translation_df <- translate_atac_barcodes(atac_barcodes, barcode_mapping, sample_name)

    # Filter to successfully translated cells
    success_idx <- which(translation_df$translation_successful)

    if (length(success_idx) == 0) {
      warning(sprintf("No successful translations for %s - skipping", sample_name))
      next
    }

    # Subset ATAC matrix to translated cells
    translated_counts <- atac_info$counts[, success_idx, drop = FALSE]

    # Rename columns to RNA barcode format
    colnames(translated_counts) <- translation_df$rna_full_barcode[success_idx]

    cat(sprintf("    Translated ATAC matrix: %d peaks x %d cells\n",
                nrow(translated_counts), ncol(translated_counts)))

    translated_atac_list[[i]] <- list(
      sample = sample_name,
      counts = translated_counts,
      fragment_path = atac_info$fragment_path,
      translation_df = translation_df
    )
  }

  if (length(translated_atac_list) == 0) {
    stop("No ATAC samples were successfully translated")
  }

  # STEP 4: Find paired cells (cells present in both RNA and ATAC)
  cat("\n  Matching cells between RNA and translated ATAC...\n")

  # Combine all translated ATAC cell barcodes
  all_atac_barcodes <- unique(unlist(lapply(translated_atac_list, function(x) {
    colnames(x$counts)
  })))

  cat(sprintf("    Total translated ATAC barcodes: %d\n", length(all_atac_barcodes)))

  # Find overlap with RNA barcodes
  paired_barcodes <- intersect(rna_barcodes, all_atac_barcodes)

  cat(sprintf("    Paired cell barcodes: %d\n", length(paired_barcodes)))

  # Calculate pairing efficiency
  rna_pairing_efficiency <- length(paired_barcodes) / length(rna_barcodes) * 100
  atac_pairing_efficiency <- length(paired_barcodes) / length(all_atac_barcodes) * 100

  cat(sprintf("\n  Pairing Efficiency:\n"))
  cat(sprintf("    RNA perspective: %d/%d = %.1f%%\n",
              length(paired_barcodes), length(rna_barcodes), rna_pairing_efficiency))
  cat(sprintf("    ATAC perspective: %d/%d = %.1f%%\n",
              length(paired_barcodes), length(all_atac_barcodes), atac_pairing_efficiency))

  if (length(paired_barcodes) == 0) {
    stop("No paired cells found after barcode translation!")
  }

  # STEP 5: Subset RNA to paired cells
  cat("\n  Subsetting RNA data to paired cells...\n")
  seurat_obj <- subset(seurat_obj, cells = paired_barcodes)

  cat(sprintf("    RNA after subsetting: %d cells\n", ncol(seurat_obj)))

  # DIAGNOSTIC: Check sample distribution after subsetting
  cat("    Sample distribution after subsetting:\n")
  if ("sample" %in% colnames(seurat_obj@meta.data)) {
    sample_counts <- table(seurat_obj@meta.data$sample)
    for (samp in names(sample_counts)) {
      cat(sprintf("      %s: %d cells\n", samp, sample_counts[samp]))
    }
  }

  # Create sample name mappings FIRST (before renaming cells)
  # RNA sample names: "Nestin_Ctrl", "Nestin_Mut", "Emx1_Ctrl", "Emx1_Mut"
  # ATAC sample names: "R26-Nestin-Ctrl-adult", "R26-Nestin-Mut-adult", etc.
  atac_to_rna_sample <- list(
    "R26-Nestin-Ctrl-adult" = "Nestin_Ctrl",
    "R26-Nestin-Mut-adult" = "Nestin_Mut",
    "R26-Emx1-Ctrl-adult" = "Emx1_Ctrl",
    "R26-Emx1-Mut-adult" = "Emx1_Mut"
  )

  # Create reverse mapping: RNA sample name → ATAC prefix (for cell renaming)
  rna_to_atac_prefix <- list(
    "Nestin_Ctrl" = "Nestin-Ctrl-adult",
    "Nestin_Mut" = "Nestin-Mut-adult",
    "Emx1_Ctrl" = "Emx1-Ctrl-adult",
    "Emx1_Mut" = "Emx1-Mut-adult"
  )

  # CRITICAL: Save sample assignments BEFORE renaming cells
  # (after renaming, metadata rownames change, making lookups unreliable)
  original_cell_to_sample <- as.character(seurat_obj$sample)  # Convert factor to character!
  names(original_cell_to_sample) <- colnames(seurat_obj)

  # DIAGNOSTIC: Check sample assignment logic
  cat("\n  DIAGNOSTIC: Verifying sample assignments...\n")
  for (samp in unique(original_cell_to_sample)) {
    sample_cells <- names(original_cell_to_sample)[original_cell_to_sample == samp]
    cat(sprintf("    Sample '%s': %d cells, barcode examples: %s\n",
                samp, length(sample_cells),
                paste(head(sample_cells, 3), collapse = ", ")))
  }

  # Rename RNA cells to match ATAC naming scheme (add sample prefix for uniqueness)
  cat("\n  Renaming RNA cells to match ATAC sample prefixes...\n")

  # DIAGNOSTIC: Test the mapping for first cell of each sample
  cat("  DIAGNOSTIC: Testing rna_to_atac_prefix mapping...\n")
  for (test_sample in unique(original_cell_to_sample)) {
    test_prefix <- rna_to_atac_prefix[[test_sample]]
    cat(sprintf("    '%s' → '%s'\n", test_sample, test_prefix))
  }

  rna_cell_names_new <- sapply(colnames(seurat_obj), function(cell_bc) {
    # Find which sample this cell belongs to (from saved mapping)
    rna_sample_name <- original_cell_to_sample[cell_bc]

    # Get the corresponding ATAC prefix
    atac_prefix <- rna_to_atac_prefix[[rna_sample_name]]

    if (is.null(atac_prefix)) {
      stop(sprintf("Unknown RNA sample name: %s (cell: %s)", rna_sample_name, cell_bc))
    }

    paste0(atac_prefix, "_", cell_bc)
  })
  # Update Seurat object cell names
  seurat_obj <- RenameCells(seurat_obj, new.names = rna_cell_names_new)
  cat(sprintf("    ✓ Renamed %d RNA cells with sample prefixes\n", length(rna_cell_names_new)))

  # DIAGNOSTIC: Check renaming results
  cat("\n  DIAGNOSTIC: Checking renamed cells...\n")
  for (samp in unique(original_cell_to_sample)) {
    sample_original_cells <- names(original_cell_to_sample)[original_cell_to_sample == samp]
    sample_renamed_cells <- rna_cell_names_new[sample_original_cells]
    cat(sprintf("    Sample '%s': renamed examples: %s\n",
                samp, paste(head(sample_renamed_cells, 2), collapse = ", ")))
  }

  # Create mapping: RNA sample name -> list of paired RNA cells (with ATAC-style prefixes)
  # Use the original cell-to-sample mapping (before renaming)
  rna_cells_by_rna_sample <- split(rna_cell_names_new, original_cell_to_sample)

  cat(sprintf("\n    RNA cells per sample:\n"))
  for (rna_samp in names(rna_cells_by_rna_sample)) {
    cat(sprintf("      %s: %d cells\n", rna_samp, length(rna_cells_by_rna_sample[[rna_samp]])))
  }

  # STEP 6: Combine ATAC data for paired cells
  cat("\n  Combining ATAC data for paired cells...\n")

  # Check for cached consensus peaks
  consensus_cache_file <- file.path(cache_dir, "consensus_peaks.rds")
  if (file.exists(consensus_cache_file)) {
    cat("    ✓ Loading cached consensus peaks...\n")
    consensus_peaks <- readRDS(consensus_cache_file)
    cat(sprintf("      Loaded %d consensus peaks from cache\n", length(consensus_peaks)))
  } else {
    # Create consensus peak set using GenomicRanges::reduce()
    # This is the Signac-recommended approach for merging samples with different peaks
    # Reference: https://stuartlab.org/signac/articles/merging

    cat("    Converting peak names to genomic ranges...\n")
    all_peaks <- lapply(translated_atac_list, function(x) rownames(x$counts))

    # Convert peak names to GRanges objects
    peak_granges_list <- lapply(seq_along(all_peaks), function(i) {
    peaks <- all_peaks[[i]]
    sample_name <- names(translated_atac_list)[i]
    cat(sprintf("      Sample %s: %d peaks before merging\n", sample_name, length(peaks)))

    # Use Signac's StringToGRanges to parse chr:start-end format
    gr <- tryCatch({
      StringToGRanges(peaks, sep = c(":", "-"))
    }, error = function(e) {
      stop(sprintf("Failed to parse peak names for sample %s: %s\n", sample_name, e$message))
    })

    return(gr)
  })

  # Create unified consensus peak set using reduce()
  # reduce() merges all overlapping peaks into single unified regions
  cat("\n    Creating consensus peak set with reduce()...\n")

  # Combine all peaks from all samples
  all_peaks_combined <- do.call(c, peak_granges_list)
  cat(sprintf("      Total peaks before reduce: %d\n", length(all_peaks_combined)))

  # Merge overlapping peaks
  consensus_peaks <- reduce(all_peaks_combined)
  cat(sprintf("      Consensus peaks after reduce: %d\n", length(consensus_peaks)))

  # Optional: Filter by peak width (remove very small or very large peaks)
  peak_widths <- width(consensus_peaks)
  cat(sprintf("      Peak width range: %d - %d bp (median: %d bp)\n",
              min(peak_widths), max(peak_widths), median(peak_widths)))

  # Filter peaks: remove very small (<20bp) or very large (>10kb) peaks
  peaks_before_filter <- length(consensus_peaks)
  consensus_peaks <- consensus_peaks[peak_widths >= 20 & peak_widths <= 10000]
  cat(sprintf("      Peaks after width filtering: %d (removed %d peaks)\n",
              length(consensus_peaks), peaks_before_filter - length(consensus_peaks)))

  # Validate overlap coverage for each sample
  cat("\n    Validating consensus peaks cover original peaks...\n")
  overlap_stats <- lapply(seq_along(peak_granges_list), function(i) {
    sample_gr <- peak_granges_list[[i]]
    overlaps <- findOverlaps(sample_gr, consensus_peaks)
    n_covered <- length(unique(queryHits(overlaps)))
    pct_covered <- 100 * n_covered / length(sample_gr)

    sample_name <- names(translated_atac_list)[i]
    cat(sprintf("      Sample %s: %d/%d peaks covered (%.1f%%)\n",
                sample_name, n_covered, length(sample_gr), pct_covered))

    return(list(
      sample = sample_name,
      original_peaks = length(sample_gr),
      covered_peaks = n_covered,
      pct_covered = pct_covered
    ))
  })

    # Warn if coverage is poor
    min_coverage <- min(sapply(overlap_stats, function(x) x$pct_covered))
    if (min_coverage < 80) {
      warning(sprintf("⚠️  Poor consensus peak coverage: minimum %.1f%% of original peaks covered.\n",
                     min_coverage),
             "   This may indicate very different chromatin accessibility between samples.\n")
    }

    # Save consensus peaks to cache
    cat("\n    💾 Saving consensus peaks to cache...\n")
    saveRDS(consensus_peaks, file = consensus_cache_file)
    cat(sprintf("      ✓ Saved %d peaks to: %s\n", length(consensus_peaks), consensus_cache_file))
  }

  # CRITICAL: Validate that consensus peaks exist
  if (length(consensus_peaks) == 0) {
    stop("❌ No consensus peaks created!\n",
         "   This means peak calling produced completely different peaks for each sample.\n",
         "   Please check:\n",
         "   1. Peak calling parameters (too stringent?)\n",
         "   2. Sample quality (enough fragments?)\n",
         "   3. Whether samples are from the same genome/species\n")
  }

  if (length(consensus_peaks) < 10000) {
    warning(sprintf("⚠️  Fewer consensus peaks than expected (%d)! This may indicate:\n", length(consensus_peaks)),
           "   - Different peak calling parameters between samples\n",
           "   - Very different chromatin accessibility profiles\n",
           "   - Poor quality ATAC data\n",
           "   Expected: 40,000-100,000 peaks for typical ATAC-seq\n")
  }

  # MODIFIED APPROACH: Map existing chromap quantifications to consensus peaks
  # Instead of requantifying with FeatureMatrix (which requires matching fragment barcodes),
  # we use the existing chromap-quantified peak-cell matrices and map them to consensus peaks
  cat("\n    Mapping existing peak quantifications to consensus peak set...\n")
  cat("    Note: Using chromap-quantified data (fragment files have incompatible barcodes)\n")

  atac_counts_list <- list()

  for (i in seq_along(translated_atac_list)) {
    atac_info <- translated_atac_list[[i]]
    atac_sample_name <- atac_info$sample  # e.g., "R26-Nestin-Ctrl-adult"

    # Map ATAC sample name to RNA sample name
    rna_sample_name <- atac_to_rna_sample[[atac_sample_name]]

    if (is.null(rna_sample_name)) {
      cat(sprintf("      Sample %s: no mapping to RNA sample name, skipping...\n", atac_sample_name))
      next
    }

    # Get RNA cells for this specific sample (these are already prefixed with ATAC-style names)
    rna_cells_this_sample <- rna_cells_by_rna_sample[[rna_sample_name]]

    if (is.null(rna_cells_this_sample) || length(rna_cells_this_sample) == 0) {
      cat(sprintf("      Sample %s (RNA: %s): no RNA cells, skipping...\n",
                  atac_sample_name, rna_sample_name))
      next
    }

    # Get ATAC cells and create prefixed names using correct mapping
    atac_cells <- colnames(atac_info$counts)
    sample_prefix <- rna_to_atac_prefix[[rna_sample_name]]  # e.g., "Nestin-Ctrl-adult"

    if (is.null(sample_prefix)) {
      cat(sprintf("      Sample %s (RNA: %s): no prefix mapping, skipping...\n",
                  atac_sample_name, rna_sample_name))
      next
    }

    atac_cells_prefixed <- paste0(sample_prefix, "_", atac_cells)

    # DIAGNOSTIC: Show what we're trying to match
    cat(sprintf("      DEBUG: Sample prefix used: %s\n", sample_prefix))
    cat(sprintf("      DEBUG: ATAC cells: %d (e.g., %s)\n",
                length(atac_cells_prefixed),
                paste(head(atac_cells_prefixed, 2), collapse = ", ")))
    cat(sprintf("      DEBUG: RNA cells for sample '%s': %d (e.g., %s)\n",
                rna_sample_name,
                length(rna_cells_this_sample),
                paste(head(rna_cells_this_sample, 2), collapse = ", ")))

    # Only include ATAC cells that have matching RNA cells from this sample
    paired_atac_cells_unique <- intersect(atac_cells_prefixed, rna_cells_this_sample)

    if (length(paired_atac_cells_unique) == 0) {
      cat(sprintf("      Sample %s: no matching RNA-ATAC pairs, skipping...\n", atac_sample_name))
      next
    }

    cat(sprintf("      Sample %s: processing %d paired cells...\n",
                atac_sample_name, length(paired_atac_cells_unique)))

    # Extract the original unprefixed barcodes for indexing into the ATAC matrix
    paired_atac_cells <- sub(paste0("^", sample_prefix, "_"), "", paired_atac_cells_unique)

    # Get the existing peak-cell matrix from chromap
    existing_counts <- atac_info$counts
    existing_peaks_str <- rownames(existing_counts)

    # Convert existing peak names to GRanges
    # Peak format: chr10:100006770-100006970
    cat(sprintf("        Converting %d existing peaks to GRanges...\n", length(existing_peaks_str)))
    time_start_conversion <- Sys.time()
    existing_peaks_gr <- StringToGRanges(existing_peaks_str, sep = c(":", "-"))
    time_conversion <- as.numeric(Sys.time() - time_start_conversion)
    cat(sprintf("        ✓ Converted in %.2f sec\n", time_conversion))

    # Ensure GRanges are sorted for fast overlap detection
    if (is.unsorted(consensus_peaks)) {
      cat("        Sorting consensus peaks...\n")
      consensus_peaks <- sort(consensus_peaks)
    }
    if (is.unsorted(existing_peaks_gr)) {
      cat("        Sorting existing peaks...\n")
      existing_peaks_gr <- sort(existing_peaks_gr)
    }

    # Find overlaps between existing peaks and consensus peaks
    cat(sprintf("        Finding overlaps with %d consensus peaks...\n", length(consensus_peaks)))
    time_start_overlap <- Sys.time()
    overlaps <- findOverlaps(existing_peaks_gr, consensus_peaks, type = "any")
    time_overlap <- as.numeric(Sys.time() - time_start_overlap)
    cat(sprintf("        ✓ Found %d overlaps in %.2f sec\n", length(overlaps), time_overlap))

    # For each consensus peak, sum counts from all overlapping existing peaks
    # This handles cases where consensus peaks may merge multiple sample-specific peaks
    consensus_peak_names <- GRangesToString(consensus_peaks)

    # OPTIMIZED AGGREGATION USING TRIPLET FORMAT
    # Instead of row-by-row assignment (very slow for sparse matrices),
    # we build the matrix in one shot from triplets
    cat("        Aggregating counts using vectorized triplet approach...\n")
    time_start_agg <- Sys.time()

    if (length(overlaps) == 0) {
      warning("No overlaps found between existing and consensus peaks!")
      consensus_counts <- Matrix::sparseMatrix(
        i = integer(0),
        j = integer(0),
        x = numeric(0),
        dims = c(length(consensus_peaks), length(paired_atac_cells)),
        dimnames = list(consensus_peak_names, paired_atac_cells_unique)
      )
    } else {
      # Convert existing counts to COO (coordinate) format
      # This gives us triplets: (peak_idx, cell_idx, count)
      existing_coo <- summary(existing_counts[, paired_atac_cells, drop = FALSE])
      # Columns: i (peak index), j (cell index), x (count)

      # Create mapping: existing_peak_idx -> consensus_peak_idx
      peak_mapping <- data.frame(
        existing_idx = queryHits(overlaps),
        consensus_idx = subjectHits(overlaps)
      )

      # Map each non-zero entry to consensus peaks
      existing_coo$existing_idx <- existing_coo$i
      coo_mapped <- merge(existing_coo, peak_mapping,
                          by = "existing_idx",
                          all.x = FALSE)  # Only keep peaks with overlaps

      if (nrow(coo_mapped) == 0) {
        warning("No counts mapped to consensus peaks!")
        consensus_counts <- Matrix::sparseMatrix(
          i = integer(0),
          j = integer(0),
          x = numeric(0),
          dims = c(length(consensus_peaks), length(paired_atac_cells)),
          dimnames = list(consensus_peak_names, paired_atac_cells_unique)
        )
      } else {
        # Aggregate counts by (consensus_peak, cell)
        # Multiple existing peaks can map to same consensus peak
        aggregated <- aggregate(
          x ~ consensus_idx + j,
          data = coo_mapped,
          FUN = sum
        )

        # Create sparse matrix from aggregated triplets
        consensus_counts <- Matrix::sparseMatrix(
          i = aggregated$consensus_idx,
          j = aggregated$j,
          x = aggregated$x,
          dims = c(length(consensus_peaks), length(paired_atac_cells)),
          dimnames = list(consensus_peak_names, paired_atac_cells_unique)
        )
      }
    }

    time_agg <- as.numeric(Sys.time() - time_start_agg)
    cat(sprintf("        ✓ Aggregated in %.2f sec\n", time_agg))

    # Report mapping statistics
    n_consensus_with_counts <- sum(Matrix::rowSums(consensus_counts) > 0)
    pct_mapped <- 100 * n_consensus_with_counts / length(consensus_peaks)

    cat(sprintf("        Mapping results:\n"))
    cat(sprintf("          Existing peaks: %d\n", nrow(existing_counts)))
    cat(sprintf("          Consensus peaks: %d\n", length(consensus_peaks)))
    cat(sprintf("          Overlapping peak pairs: %d\n", length(overlaps)))
    cat(sprintf("          Consensus peaks with counts: %d (%.1f%%)\n",
                n_consensus_with_counts, pct_mapped))

    if (pct_mapped < 50) {
      warning(sprintf("Sample %s: Only %.1f%% of consensus peaks have counts after mapping!\n",
                     atac_sample_name, pct_mapped),
             "   This suggests poor overlap between sample peaks and consensus peaks.\n")
    }

    cat(sprintf("        Final matrix: %d peaks x %d cells\n",
                nrow(consensus_counts), ncol(consensus_counts)))

    atac_counts_list[[i]] <- consensus_counts
  }

  # Combine ATAC counts
  combined_atac <- do.call(cbind, atac_counts_list)

  cat(sprintf("    Combined ATAC: %d peaks x %d cells\n",
              nrow(combined_atac), ncol(combined_atac)))

  # CRITICAL VERIFICATION: Ensure exact match between RNA and ATAC cells
  cat("\n  Verifying RNA-ATAC cell matching...\n")
  rna_cells <- colnames(seurat_obj)
  atac_cells <- colnames(combined_atac)

  cat(sprintf("    RNA cells: %d\n", length(rna_cells)))
  cat(sprintf("    ATAC cells: %d\n", length(atac_cells)))

  if (length(rna_cells) != length(atac_cells)) {
    stop(sprintf("Cell count mismatch! RNA: %d cells, ATAC: %d cells\n",
                 length(rna_cells), length(atac_cells)),
         "This should not happen - check the pairing logic above.")
  }

  # Ensure cells are in the same order
  if (!all(rna_cells == atac_cells)) {
    cat("    Reordering ATAC cells to match RNA cell order...\n")
    combined_atac <- combined_atac[, rna_cells, drop = FALSE]
  }

  cat(sprintf("    ✓ Perfect match: %d cells in both RNA and ATAC\n", length(rna_cells)))

  # Track sample origin for each cell
  cat("\n  Creating sample metadata...\n")
  sample_metadata <- data.frame(
    cell_barcode = colnames(combined_atac),
    atac_sample = NA_character_,
    stringsAsFactors = FALSE
  )
  rownames(sample_metadata) <- sample_metadata$cell_barcode

  # Assign sample labels based on which matrix each cell came from
  for (i in seq_along(translated_atac_list)) {
    atac_info <- translated_atac_list[[i]]
    sample_name <- atac_info$sample

    if (!is.null(atac_counts_list[[i]])) {
      sample_cells <- colnames(atac_counts_list[[i]])
      sample_metadata[sample_cells, "atac_sample"] <- sample_name
    }
  }

  cat(sprintf("    ✓ Sample assignments:\n"))
  sample_table <- table(sample_metadata$atac_sample)
  for (samp in names(sample_table)) {
    cat(sprintf("      %s: %d cells\n", samp, sample_table[samp]))
  }


  # STEP 7: Create ChromatinAssay
  cat("\n================================================================================\n")
  cat("STEP 7: Creating ChromatinAssay\n")
  cat("================================================================================\n")

  # Get gene annotations
  # Using Ensembl v79 (corresponds to mm10/GRCm38 genome build)
  # IMPORTANT: This should match the genome annotation used for RNA-seq quantification
  # If your RNA data was processed with a different Ensembl version, update this accordingly
  cat("    Loading gene annotations: EnsDb.Mmusculus.v79 (Ensembl 79, mm10/GRCm38)\n")
  cat("    ⚠️  Verify this matches your RNA-seq reference annotation!\n")

  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
  genome(annotations) <- "mm10"

  cat(sprintf("    ✓ Loaded %d genes from annotation\n", length(annotations)))

  # STEP 7.5: Translate Fragment Files (if tools available)
  cat("\n  Checking for fragment translation capability...\n")

  has_bgzip <- system("which bgzip", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  has_tabix <- system("which tabix", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0

  translated_fragments <- NULL
  fragment_objects <- list()

  if (!has_bgzip || !has_tabix) {
    cat("    ⚠️  bgzip/tabix not found - fragment translation disabled\n")
    cat("       Fragment-based QC (TSS enrichment, nucleosome signal) will be unavailable\n")
    cat("       To enable: conda install -c bioconda htslib\n")
  } else {
    cat("    ✓ Tools available: bgzip, tabix\n")
    cat("    Translating fragment files for paired cells...\n")

    translated_frag_dir <- file.path(OUTPUT_DIR, "translated_fragments")
    dir.create(translated_frag_dir, recursive = TRUE, showWarnings = FALSE)

    translated_fragments <- list()

    # Extract atac_to_rna mapping for fragment translation
    atac_to_rna <- barcode_mapping$atac_to_rna

    for (atac_info in translated_atac_list) {
      sample_name <- atac_info$sample
      fragment_file <- atac_info$fragment_path

      # Determine sample suffix using explicit gem group mapping (from h5ad barcode analysis)
      sample_suffix_map <- list(
        "R26-Emx1-Ctrl-adult" = "1-0",     # Emx1 Ctrl gem group
        "R26-Emx1-Mut-adult" = "1-1",      # Emx1 Mut gem group
        "R26-Nestin-Ctrl-adult" = "1-2",   # Nestin Ctrl gem group
        "R26-Nestin-Mut-adult" = "1-3"     # Nestin Mut gem group
      )

      sample_suffix <- sample_suffix_map[[sample_name]]
      if (is.null(sample_suffix)) {
        cat(sprintf("        ⚠ Unknown sample name: %s\n", sample_name))
        cat(sprintf("        Known samples: %s\n", paste(names(sample_suffix_map), collapse = ", ")))
        next
      }

      cat(sprintf("        Using gem group suffix: %s\n", sample_suffix))

      cat(sprintf("      Processing %s...\n", sample_name))

      if (!file.exists(fragment_file)) {
        cat(sprintf("        ⚠ Fragment file not found: %s\n", fragment_file))
        next
      }

      # Get paired cells for this sample
      sample_cells <- paired_barcodes[sample_metadata$atac_sample == sample_name]

      if (length(sample_cells) == 0) {
        cat(sprintf("        ⚠ No paired cells for %s\n", sample_name))
        next
      }

      # Output file path
      output_file <- file.path(translated_frag_dir,
                               sprintf("%s_fragments_translated.tsv", sample_name))

      # Translate fragment file
      tryCatch({
        translated_file <- translate_fragment_file(
          fragment_file = fragment_file,
          atac_to_rna = atac_to_rna,
          sample_suffix = sample_suffix,
          output_file = output_file,
          cells_to_keep = sample_cells
        )

        translated_fragments[[sample_name]] <- translated_file

      }, error = function(e) {
        cat(sprintf("        ⚠ Translation failed: %s\n", e$message))
      })
    }

    if (length(translated_fragments) > 0) {
      cat(sprintf("    ✓ Successfully translated %d fragment files\n",
                  length(translated_fragments)))

      # Create Fragment objects for Signac
      for (sample_name in names(translated_fragments)) {
        frag_file <- translated_fragments[[sample_name]]
        sample_cells <- paired_barcodes[sample_metadata$atac_sample == sample_name]

        tryCatch({
          frag_obj <- CreateFragmentObject(
            path = frag_file,
            cells = sample_cells,
            validate.fragments = FALSE  # Skip validation for speed
          )
          fragment_objects[[sample_name]] <- frag_obj
          cat(sprintf("      ✓ Created Fragment object for %s\n", sample_name))
        }, error = function(e) {
          cat(sprintf("      ⚠ Failed to create Fragment object for %s: %s\n",
                      sample_name, e$message))
        })
      }
    } else {
      cat("    ⚠ No fragment files were successfully translated\n")
    }
  }

  # Create ChromatinAssay with translated fragment files (if available)
  if (length(fragment_objects) > 0) {
    cat("\n  Creating ChromatinAssay WITH translated fragment files...\n")
    cat("    This enables fragment-based QC and coverage analysis!\n")

    chrom_assay <- CreateChromatinAssay(
      counts = combined_atac,
      sep = c(":", "-"),
      annotation = annotations,
      genome = "mm10",
      fragments = fragment_objects
    )

    cat("    ✓ ChromatinAssay created with fragment file support\n")
    cat("      TSS enrichment: ENABLED\n")
    cat("      Nucleosome signal: ENABLED\n")
    cat("      Coverage plots: ENABLED\n")

  } else {
    cat("\n  Creating ChromatinAssay WITHOUT fragment files...\n")
    cat("    (Fragment-based QC will be unavailable)\n")

    chrom_assay <- CreateChromatinAssay(
      counts = combined_atac,
      sep = c(":", "-"),
      annotation = annotations,
      genome = "mm10"
    )

    cat("    ⚠ ChromatinAssay created without fragment support\n")
    cat("      TSS enrichment: DISABLED\n")
    cat("      Nucleosome signal: DISABLED\n")
    cat("      Coverage plots: Will need separate script\n")
  }

  # CRITICAL: Multi-layer ATAC normalization for proper gene-peak linkage
  cat("\n  Creating ATAC normalization layers...\n")
  cat("  Following Python pipeline approach (corrected_03_load_translate_atac.py:414-510)\n")

  # LAYER 1: Raw counts (already in chrom_assay@counts)
  cat("    [1/3] Raw counts already in ChromatinAssay@counts\n")

  # For now, don't apply any normalization to the ChromatinAssay data slot
  # We'll apply TF-IDF after adding to Seurat object (Signac best practice)
  cat("    Normalization will be applied after adding to Seurat object\n")

  cat("\n  ✅ ATAC data prepared:\n")
  cat("     • counts: Raw consensus peak counts (mapped from chromap quantifications)\n")
  cat("     • Normalization: Will be applied after adding to Seurat object\n")

  # STEP 8: Store barcode translation metadata
  cat("\n  Storing barcode translation metadata...\n")

  # Build complete ATAC → RNA barcode mapping for all paired cells
  atac_to_rna_full_map <- data.frame(
    rna_barcode = character(),
    original_atac_barcode = character(),
    stringsAsFactors = FALSE
  )

  for (atac_info in translated_atac_list) {
    translation_df <- atac_info$translation_df
    successful <- translation_df$translation_successful

    for (i in which(successful)) {
      atac_bc <- translation_df$original_atac_barcode[i]
      rna_bc <- translation_df$rna_full_barcode[i]
      if (rna_bc %in% paired_barcodes) {
        atac_to_rna_full_map <- rbind(atac_to_rna_full_map,
                                       data.frame(rna_barcode = rna_bc,
                                                 original_atac_barcode = atac_bc,
                                                 stringsAsFactors = FALSE))
      }
    }
  }

  # Add barcode translation to metadata for each cell
  seurat_obj$original_atac_barcode <- atac_to_rna_full_map$original_atac_barcode[
    match(colnames(seurat_obj), atac_to_rna_full_map$rna_barcode)
  ]

  # Add sample metadata
  seurat_obj$atac_sample <- sample_metadata$atac_sample[
    match(colnames(seurat_obj), sample_metadata$cell_barcode)
  ]

  cat(sprintf("    ✓ Barcode translation stored in metadata for %d cells\n",
              sum(!is.na(seurat_obj$original_atac_barcode))))
  cat(sprintf("    ✓ Sample assignments stored in metadata\n"))

  # STEP 9: Add ATAC assay to Seurat object
  cat("\n  Adding ATAC assay to Seurat object...\n")
  seurat_obj[["ATAC"]] <- chrom_assay
  DefaultAssay(seurat_obj) <- "ATAC"

  # STEP 10: Apply TF-IDF normalization (Signac best practice: after merging)
  cat("\n  Applying TF-IDF normalization to merged ATAC data...\n")
  cat("    💡 Running AFTER merging ensures same IDF weights across all samples\n")
  cat("    💡 LinkPeaks will use raw counts (peak.slot='counts') by default\n")

  seurat_obj <- RunTFIDF(seurat_obj, assay = "ATAC")

  # Verify TF-IDF was applied
  tfidf_data <- GetAssayData(seurat_obj, assay = "ATAC", slot = "data")
  if (inherits(tfidf_data, "sparseMatrix")) {
    tfidf_range <- range(tfidf_data@x)
  } else {
    tfidf_range <- range(as.vector(tfidf_data))
  }
  cat(sprintf("    ✓ TF-IDF applied to data slot: range [%.6f, %.6f]\n", tfidf_range[1], tfidf_range[2]))

  # NOTE about Fragment files
  if (length(fragment_objects) > 0) {
    cat("\n  ✅ FRAGMENT FILES STATUS:\n")
    cat("      Fragment files successfully translated and linked!\n")
    cat("      Fragment-dependent QC metrics (NucleosomeSignal, TSS enrichment): AVAILABLE\n")
    cat("      Coverage plots and footprinting: AVAILABLE\n")
    cat("      Peak-based analysis (LinkPeaks, differential accessibility): AVAILABLE\n")
  } else {
    cat("\n  ⚠️  FRAGMENT FILES STATUS:\n")
    cat("      Fragment files not available (bgzip/tabix not found or translation failed)\n")
    cat("      Fragment-dependent QC metrics (NucleosomeSignal, TSS enrichment): UNAVAILABLE\n")
    cat("      Peak-based analysis (LinkPeaks, differential accessibility): Still works fine using quantified counts\n")
    cat("      To enable fragments: Install htslib (conda install -c bioconda htslib) and rerun\n")
  }

  cat("\n  ✅ Integrated multiome object created successfully!\n")
  cat(sprintf("    Paired cells: %d\n", ncol(seurat_obj)))
  cat(sprintf("    RNA features: %d\n", nrow(seurat_obj[["RNA"]])))
  cat(sprintf("    ATAC features: %d peaks (consensus across samples)\n", nrow(seurat_obj[["ATAC"]])))
  cat(sprintf("    Pairing efficiency: %.1f%%\n", rna_pairing_efficiency))

  cat("\n  📊 ATAC Data Layers:\n")
  cat("     • counts: Raw consensus peak counts (mapped from chromap quantifications)\n")
  cat("     • data: TF-IDF normalized (via RunTFIDF) for dimensionality reduction\n")
  cat("     • LinkPeaks: Uses raw counts slot by default (standard Signac workflow)\n")

  # Save integrated object to cache for faster reruns
  if (!is.null(cache_dir)) {
    cat("\n  💾 Saving integrated object to cache...\n")
    cache_file <- file.path(cache_dir, "integrated_seurat_object.rds")
    saveRDS(seurat_obj, file = cache_file, compress = "gzip")
    cat(sprintf("    ✓ Saved to: %s\n", cache_file))
    cat(sprintf("    File size: %.1f MB\n", file.size(cache_file) / 1e6))
    cat("    Next run will load from cache (delete to recompute)\n")
  }

  return(seurat_obj)
}

################################################################################
# Function: Compute ATAC QC Metrics
################################################################################

compute_atac_qc <- function(seurat_obj) {
  cat("\n=== Computing ATAC QC metrics ===\n")

  # Check if ATAC assay exists
  if (!"ATAC" %in% names(seurat_obj@assays)) {
    stop("ERROR: ATAC assay not found in Seurat object!\n",
         "   Please check ChromatinAssay creation step.\n")
  }

  # Check if ATAC assay is valid (not NULL)
  if (is.null(seurat_obj[["ATAC"]])) {
    stop("ERROR: ATAC assay is NULL!\n",
         "   ChromatinAssay creation failed.\n")
  }

  DefaultAssay(seurat_obj) <- "ATAC"

  # Check if fragment files are attached
  fragments_list <- Fragments(seurat_obj)
  has_fragments <- !is.null(fragments_list) && length(fragments_list) > 0

  if (!has_fragments) {
    cat("  ⚠️  NO fragment files detected - skipping fragment-dependent QC metrics\n")
    cat("     (Nucleosome signal and TSS enrichment require fragment files)\n")
    cat("     To enable: Install htslib and rerun Step 1\n")
  } else {
    cat(sprintf("  ✓ Fragment files detected (%d samples) - computing fragment-based QC\n",
                length(fragments_list)))

    # Nucleosome signal
    cat("  Computing nucleosome signal...\n")
    tryCatch({
      seurat_obj <- NucleosomeSignal(seurat_obj)
      cat("    ✓ Nucleosome signal computed\n")

      # Report stats
      ns_median <- median(seurat_obj$nucleosome_signal, na.rm = TRUE)
      ns_q25 <- quantile(seurat_obj$nucleosome_signal, 0.25, na.rm = TRUE)
      ns_q75 <- quantile(seurat_obj$nucleosome_signal, 0.75, na.rm = TRUE)
      cat(sprintf("      Median: %.2f (IQR: %.2f - %.2f)\n", ns_median, ns_q25, ns_q75))
      cat(sprintf("      Cells with NS < 2: %d (%.1f%%)\n",
                  sum(seurat_obj$nucleosome_signal < 2, na.rm = TRUE),
                  100 * mean(seurat_obj$nucleosome_signal < 2, na.rm = TRUE)))

    }, error = function(e) {
      cat(sprintf("    ⚠ Nucleosome signal failed: %s\n", e$message))
    })

    # TSS enrichment
    cat("  Computing TSS enrichment...\n")
    tryCatch({
      seurat_obj <- TSSEnrichment(seurat_obj, fast = TRUE)
      cat("    ✓ TSS enrichment computed\n")

      # Report stats
      tss_median <- median(seurat_obj$TSS.enrichment, na.rm = TRUE)
      tss_q25 <- quantile(seurat_obj$TSS.enrichment, 0.25, na.rm = TRUE)
      tss_q75 <- quantile(seurat_obj$TSS.enrichment, 0.75, na.rm = TRUE)
      cat(sprintf("      Median: %.2f (IQR: %.2f - %.2f)\n", tss_median, tss_q25, tss_q75))
      cat(sprintf("      Cells with TSS > 2: %d (%.1f%%)\n",
                  sum(seurat_obj$TSS.enrichment > 2, na.rm = TRUE),
                  100 * mean(seurat_obj$TSS.enrichment > 2, na.rm = TRUE)))

    }, error = function(e) {
      cat(sprintf("    ⚠ TSS enrichment failed: %s\n", e$message))
    })
  }

  # Blacklist ratio (should work without fragments)
  cat("  Computing blacklist ratio...\n")
  tryCatch({
    seurat_obj$blacklist_ratio <- FractionCountsInRegion(
      object = seurat_obj,
      assay = "ATAC",
      regions = blacklist_mm10
    )
    cat("    ✓ Blacklist ratio computed\n")
  }, error = function(e) {
    cat(sprintf("    ⚠ Blacklist ratio failed: %s\n", e$message))
    cat("    → Skipping\n")
  })

  # Basic QC metrics that don't require fragments
  cat("  Computing basic peak counts...\n")
  if (!"nFeature_ATAC" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$nFeature_ATAC <- colSums(GetAssayData(seurat_obj, assay = "ATAC", slot = "counts") > 0)
  }
  if (!"nCount_ATAC" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$nCount_ATAC <- colSums(GetAssayData(seurat_obj, assay = "ATAC", slot = "counts"))
  }
  cat("    ✓ Basic peak counts computed\n")

  cat("\n  ✓ ATAC QC metrics completed (with fallbacks for fragment-dependent metrics)\n")

  return(seurat_obj)
}

################################################################################
# Function: Create QC Plots
################################################################################

create_qc_plots <- function(seurat_obj, output_dir) {
  cat("\n=== Creating QC plots ===\n")

  # Check for NA/Inf values in QC metrics
  qc_metrics <- c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal', 'blacklist_ratio')
  for (metric in qc_metrics) {
    if (metric %in% colnames(seurat_obj@meta.data)) {
      vals <- seurat_obj@meta.data[[metric]]
      n_na <- sum(is.na(vals))
      n_inf <- sum(is.infinite(vals))
      if (n_na > 0 || n_inf > 0) {
        cat(sprintf("  ⚠ %s has %d NA and %d Inf values - will filter for plots\n",
                   metric, n_na, n_inf))
      }
    }
  }

  # Check which QC metrics are available
  has_tss <- "TSS.enrichment" %in% colnames(seurat_obj@meta.data) &&
             sum(!is.na(seurat_obj$TSS.enrichment)) > 0
  has_nucleosome <- "nucleosome_signal" %in% colnames(seurat_obj@meta.data) &&
                    sum(!is.na(seurat_obj$nucleosome_signal)) > 0
  has_blacklist <- "blacklist_ratio" %in% colnames(seurat_obj@meta.data) &&
                   sum(!is.na(seurat_obj$blacklist_ratio)) > 0

  # Check if fragment files are attached
  fragments_list <- Fragments(seurat_obj)
  has_fragments <- !is.null(fragments_list) && length(fragments_list) > 0

  if (!has_tss && !has_nucleosome) {
    cat("  ⚠️  Fragment-dependent QC metrics not available - skipping density plots\n")
    cat("     (This is expected when using chromap-quantified data without fragment files)\n")
    p1 <- ggplot() +
      ggtitle("TSS/Nucleosome metrics unavailable") +
      annotate("text", x = 0.5, y = 0.5,
               label = "Fragment-dependent QC metrics\nnot computed (no fragment files)",
               size = 5) +
      theme_void()
    p2 <- p1
  } else {
    # Filter cells with valid QC metrics for density plots
    valid_cells <- !is.na(seurat_obj$nCount_ATAC) &
                   !is.infinite(seurat_obj$nCount_ATAC)

    if (has_tss) {
      valid_cells <- valid_cells &
                     !is.na(seurat_obj$TSS.enrichment) &
                     !is.infinite(seurat_obj$TSS.enrichment)
    }

    if (has_nucleosome) {
      valid_cells <- valid_cells &
                     !is.na(seurat_obj$nucleosome_signal) &
                     !is.infinite(seurat_obj$nucleosome_signal)
    }

    n_valid <- sum(valid_cells)
    n_total <- ncol(seurat_obj)
    cat(sprintf("  Using %d/%d cells with valid QC metrics for density plots\n", n_valid, n_total))

    if (n_valid < 10) {
      cat("  ⚠ Too few valid cells for density plots - skipping\n")
      p1 <- p2 <- ggplot() + ggtitle("Insufficient valid data")
    } else {
      # Extract valid data directly from metadata
      meta_valid <- seurat_obj@meta.data[valid_cells, ]

      # Create density plots from metadata
      tryCatch({
        if (has_tss) {
          p1 <- ggplot(meta_valid, aes(x = nCount_ATAC, y = TSS.enrichment)) +
            geom_point(alpha = 0.3, size = 0.5) +
            geom_density_2d(color = "blue") +
            scale_x_log10() +
            theme_classic() +
            ggtitle("ATAC Counts vs TSS Enrichment") +
            labs(x = "nCount_ATAC (log10)", y = "TSS Enrichment")
        } else {
          p1 <- ggplot() + ggtitle("TSS enrichment not available") + theme_void()
        }

        if (has_nucleosome && has_tss) {
          p2 <- ggplot(meta_valid, aes(x = nucleosome_signal, y = TSS.enrichment)) +
            geom_point(alpha = 0.3, size = 0.5) +
            geom_density_2d(color = "blue") +
            theme_classic() +
            ggtitle("Nucleosome Signal vs TSS Enrichment") +
            labs(x = "Nucleosome Signal", y = "TSS Enrichment")
        } else {
          p2 <- ggplot() + ggtitle("Nucleosome signal not available") + theme_void()
        }
      }, error = function(e) {
        cat(sprintf("  ⚠ Failed to create density plots: %s\n", e$message))
        p1 <<- ggplot() + ggtitle("Plot generation failed")
        p2 <<- ggplot() + ggtitle("Plot generation failed")
      })
    }
  }

  # Violin plots - only plot available metrics
  vln_features <- c("nCount_RNA", "nCount_ATAC")
  if (has_tss) vln_features <- c(vln_features, "TSS.enrichment")
  if (has_nucleosome) vln_features <- c(vln_features, "nucleosome_signal")
  if (has_blacklist) vln_features <- c(vln_features, "blacklist_ratio")

  p3 <- VlnPlot(
    object = seurat_obj,
    features = vln_features,
    ncol = length(vln_features),
    pt.size = 0
  )

  # Fragment length histogram
  if (has_fragments) {
    cat("  Creating fragment length histogram...\n")
    tryCatch({
      # Try with group.by first
      if ("atac_sample" %in% colnames(seurat_obj@meta.data) &&
          length(unique(seurat_obj$atac_sample)) > 0) {
        p4 <- FragmentHistogram(object = seurat_obj, group.by = "atac_sample")
      } else {
        # Fallback: no grouping
        p4 <- FragmentHistogram(object = seurat_obj)
      }
      cat("    ✓ Fragment histogram created\n")
    }, error = function(e) {
      cat(sprintf("    ⚠ Fragment histogram failed: %s\n", e$message))
      p4 <<- ggplot() +
        ggtitle("Fragment Histogram") +
        annotate("text", x = 0.5, y = 0.5,
                label = paste0("Fragment histogram generation failed:\n", e$message),
                size = 4, hjust = 0.5) +
        theme_void()
    })
  } else {
    cat("  Fragment histogram skipped (no fragments available)\n")
    p4 <- ggplot() +
      ggtitle("Fragment Histogram Unavailable") +
      annotate("text", x = 0.5, y = 0.5,
              label = "Fragment files not available.\nInstall htslib and rerun to enable fragment-based analysis.",
              size = 5, hjust = 0.5) +
      theme_void()
  }

  # Save plots
  pdf(file.path(output_dir, "qc/01_atac_qc_overview.pdf"), width = 20, height = 12)
  tryCatch({
    print((p1 | p2) / p3)
  }, error = function(e) {
    cat(sprintf("  ⚠ Failed to print QC plots: %s\n", e$message))
    plot.new()
    text(0.5, 0.5, paste0("QC plot generation failed:\n", e$message), cex = 1.5)
  })

  tryCatch({
    print(p4)
  }, error = function(e) {
    cat(sprintf("  ⚠ Failed to print fragment histogram: %s\n", e$message))
    plot.new()
    text(0.5, 0.5, paste0("Fragment histogram failed:\n", e$message), cex = 1.5)
  })
  dev.off()

  cat(sprintf("  ✓ QC plots saved to: %s/qc/\n", output_dir))

  # Return updated object
  return(seurat_obj)
}

################################################################################
# Function: Filter Cells by QC
################################################################################

filter_cells_by_qc <- function(seurat_obj, thresholds) {
  cat("\n=== Filtering cells by QC metrics ===\n")

  n_cells_before <- ncol(seurat_obj)
  cat(sprintf("  Cells before filtering: %d\n", n_cells_before))

  # Show actual data distribution
  cat("\n  Actual QC metric distributions:\n")
  cat(sprintf("    nCount_ATAC: [%g, %g] (median: %g)\n",
              min(seurat_obj$nCount_ATAC, na.rm = TRUE),
              max(seurat_obj$nCount_ATAC, na.rm = TRUE),
              median(seurat_obj$nCount_ATAC, na.rm = TRUE)))
  cat(sprintf("    nCount_RNA: [%g, %g] (median: %g)\n",
              min(seurat_obj$nCount_RNA, na.rm = TRUE),
              max(seurat_obj$nCount_RNA, na.rm = TRUE),
              median(seurat_obj$nCount_RNA, na.rm = TRUE)))

  # Check if fragment-dependent metrics are available
  has_tss <- "TSS.enrichment" %in% colnames(seurat_obj@meta.data)
  has_nucleosome <- "nucleosome_signal" %in% colnames(seurat_obj@meta.data)
  has_blacklist <- "blacklist_ratio" %in% colnames(seurat_obj@meta.data)

  if (has_tss) {
    cat(sprintf("    TSS.enrichment: [%g, %g] (median: %g)\n",
                min(seurat_obj$TSS.enrichment, na.rm = TRUE),
                max(seurat_obj$TSS.enrichment, na.rm = TRUE),
                median(seurat_obj$TSS.enrichment, na.rm = TRUE)))
  } else {
    cat("    TSS.enrichment: NOT AVAILABLE (fragment files incompatible)\n")
  }

  if (has_nucleosome) {
    cat(sprintf("    nucleosome_signal: [%g, %g] (median: %g)\n",
                min(seurat_obj$nucleosome_signal[is.finite(seurat_obj$nucleosome_signal)], na.rm = TRUE),
                max(seurat_obj$nucleosome_signal[is.finite(seurat_obj$nucleosome_signal)], na.rm = TRUE),
                median(seurat_obj$nucleosome_signal, na.rm = TRUE)))
  } else {
    cat("    nucleosome_signal: NOT AVAILABLE (fragment files incompatible)\n")
  }

  if (has_blacklist) {
    cat(sprintf("    blacklist_ratio: [%g, %g] (median: %g, %d NA)\n",
                min(seurat_obj$blacklist_ratio, na.rm = TRUE),
                max(seurat_obj$blacklist_ratio, na.rm = TRUE),
                median(seurat_obj$blacklist_ratio, na.rm = TRUE),
                sum(is.na(seurat_obj$blacklist_ratio))))
  } else {
    cat("    blacklist_ratio: NOT AVAILABLE\n")
  }

  # Adapt thresholds based on data
  cat("\n  ⚠️  ADAPTIVE THRESHOLDS (based on actual data):\n")
  # Use 5th percentile for min, 95th percentile for max
  adapted_min_atac <- max(quantile(seurat_obj$nCount_ATAC, 0.05, na.rm = TRUE), 100)
  adapted_max_atac <- quantile(seurat_obj$nCount_ATAC, 0.95, na.rm = TRUE)

  if (has_tss) {
    tss_median <- median(seurat_obj$TSS.enrichment, na.rm = TRUE)
    # Only use TSS filter if we have real data (median > 0)
    if (tss_median > 0) {
      adapted_min_tss <- max(quantile(seurat_obj$TSS.enrichment, 0.05, na.rm = TRUE), 1)
    } else {
      adapted_min_tss <- NULL  # Skip TSS filter if fragments were empty
      cat("    ⚠ TSS enrichment is 0 for all cells (empty fragments?) - skipping TSS filter\n")
    }
  } else {
    adapted_min_tss <- NULL
  }

  cat(sprintf("    ATAC counts: > %g (was %d)\n", adapted_min_atac, thresholds$min_atac_counts))
  cat(sprintf("    ATAC counts: < %g (was %d)\n", adapted_max_atac, thresholds$max_atac_counts))
  cat(sprintf("    RNA counts: > %d\n", thresholds$min_rna_counts))
  cat(sprintf("    RNA counts: < %d\n", thresholds$max_rna_counts))

  if (has_tss) {
    cat(sprintf("    TSS enrichment: > %g (was %g)\n", adapted_min_tss, thresholds$min_tss_enrichment))
  } else {
    cat("    TSS enrichment: SKIPPED (not available)\n")
  }

  if (has_nucleosome) {
    cat(sprintf("    Nucleosome signal: < %g\n", thresholds$max_nucleosome_signal))
  } else {
    cat("    Nucleosome signal: SKIPPED (not available)\n")
  }

  if (has_blacklist) {
    cat(sprintf("    Blacklist ratio: < %g (or NA)\n", thresholds$max_blacklist_ratio))
  } else {
    cat("    Blacklist ratio: SKIPPED (not available)\n")
  }

  cat("\n  Applying filters...\n")

  # Build filter step-by-step by subsetting cells
  # Start with all cells
  cells_keep <- rep(TRUE, ncol(seurat_obj))

  # Filter by ATAC counts
  cells_keep <- cells_keep &
                seurat_obj$nCount_ATAC > adapted_min_atac &
                seurat_obj$nCount_ATAC < adapted_max_atac
  cat(sprintf("    After ATAC count filter: %d cells\n", sum(cells_keep)))

  # Filter by RNA counts
  cells_keep <- cells_keep &
                seurat_obj$nCount_RNA > thresholds$min_rna_counts &
                seurat_obj$nCount_RNA < thresholds$max_rna_counts
  cat(sprintf("    After RNA count filter: %d cells\n", sum(cells_keep)))

  # Add TSS filter if available
  if (has_tss && !is.null(adapted_min_tss)) {
    cells_keep <- cells_keep & seurat_obj$TSS.enrichment > adapted_min_tss
    cat(sprintf("    After TSS enrichment filter: %d cells\n", sum(cells_keep)))
  }

  # Add nucleosome filter if available
  if (has_nucleosome) {
    cells_keep <- cells_keep & seurat_obj$nucleosome_signal < thresholds$max_nucleosome_signal
    cat(sprintf("    After nucleosome signal filter: %d cells\n", sum(cells_keep)))
  }

  # Add blacklist filter if available
  if (has_blacklist) {
    cells_keep <- cells_keep &
                  (is.na(seurat_obj$blacklist_ratio) |
                   seurat_obj$blacklist_ratio < thresholds$max_blacklist_ratio)
    cat(sprintf("    After blacklist ratio filter: %d cells\n", sum(cells_keep)))
  }

  # Get cell names to keep
  cells_to_keep <- colnames(seurat_obj)[cells_keep]

  # Subset object
  seurat_obj <- subset(seurat_obj, cells = cells_to_keep)

  n_cells_after <- ncol(seurat_obj)

  cat(sprintf("\n  Cells after filtering: %d\n", n_cells_after))
  cat(sprintf("  Cells removed: %d (%.1f%%)\n",
              n_cells_before - n_cells_after,
              100 * (n_cells_before - n_cells_after) / n_cells_before))

  return(seurat_obj)
}

################################################################################
# Main Execution
################################################################################

main <- function() {
  cat("\n================================================================================\n")
  cat("STEP 1: Loading ATAC Data (All 4 Samples)\n")
  cat("================================================================================\n")

  # Load ATAC data for all 4 samples
  atac_data_list <- list()

  for (sample_key in names(ATAC_SAMPLES)) {
    sample_name <- ATAC_SAMPLES[[sample_key]]
    cat(sprintf("\n--- Loading ATAC for %s ---\n", sample_name))
    atac_data <- load_atac_data(sample_name, USE_PSEUDOBULK)
    atac_data_list[[sample_key]] <- atac_data
  }

  cat(sprintf("\n✓ Loaded ATAC data for %d samples\n", length(atac_data_list)))

  cat("\n================================================================================\n")
  cat("STEP 2: Loading RNA Data (All Samples from Single h5ad)\n")
  cat("================================================================================\n")

  # Load single h5ad file containing ALL samples (Nestin + Emx1)
  rna_h5ad_path <- file.path(RNA_DIR, RNA_H5AD_FILE)
  cat(sprintf("  Loading: %s\n", rna_h5ad_path))

  rna_data <- load_rna_from_h5ad(rna_h5ad_path)

  if (is.null(rna_data)) {
    stop("Failed to load RNA data from h5ad file")
  }

  cat(sprintf("\n✓ Loaded RNA data: %d cells x %d genes\n",
              ncol(rna_data$counts), nrow(rna_data$counts)))

  # Check for expected metadata columns
  if ("genotype" %in% colnames(rna_data$metadata)) {
    genotypes <- unique(rna_data$metadata$genotype)
    cat(sprintf("  Genotypes: %s\n", paste(genotypes, collapse = ", ")))
  }

  if ("condition" %in% colnames(rna_data$metadata)) {
    conditions <- unique(rna_data$metadata$condition)
    cat(sprintf("  Conditions: %s\n", paste(conditions, collapse = ", ")))
  }

  if ("sample" %in% colnames(rna_data$metadata)) {
    samples <- table(rna_data$metadata$sample)
    cat("  Sample distribution:\n")
    for (sample_name in names(samples)) {
      cat(sprintf("    %s: %d cells\n", sample_name, samples[sample_name]))
    }
  }

  cat("\n================================================================================\n")
  cat("STEP 3: Creating Integrated Object\n")
  cat("================================================================================\n")

  # Create integrated Seurat object
  seurat_obj <- create_integrated_object(rna_data, atac_data_list)

  cat("\n================================================================================\n")
  cat("STEP 3B: Validation Checks\n")
  cat("================================================================================\n")

  # Validate fragment files for all 4 samples
  cat("\n--- Fragment File Validation ---\n")
  translated_frag_dir <- file.path(OUTPUT_DIR, "translated_fragments")

  fragment_check_passed <- TRUE
  for (sample_key in names(ATAC_SAMPLES)) {
    sample_name <- ATAC_SAMPLES[[sample_key]]
    expected_frag_file <- file.path(translated_frag_dir,
                                    sprintf("%s_fragments_translated.tsv.gz", sample_name))

    if (file.exists(expected_frag_file)) {
      file_size_mb <- file.size(expected_frag_file) / 1e6
      cat(sprintf("  ✓ %s: %.2f MB\n", sample_name, file_size_mb))
    } else {
      cat(sprintf("  ✗ %s: MISSING\n", sample_name))
      fragment_check_passed <- FALSE
    }
  }

  if (!fragment_check_passed) {
    cat("\n  ⚠ WARNING: Some fragment files are missing!\n")
    cat("    BigWig track generation (step 6) may fail for missing samples.\n")
    cat("    Fragment-dependent QC metrics may also be unavailable.\n")
  } else {
    cat("\n  ✓ All 4 fragment files found\n")
  }

  # Validate barcode overlap between RNA and ATAC
  cat("\n--- RNA-ATAC Barcode Overlap Validation ---\n")

  rna_barcodes <- colnames(seurat_obj[["RNA"]])
  atac_barcodes <- colnames(seurat_obj[["ATAC"]])
  overlap_barcodes <- intersect(rna_barcodes, atac_barcodes)

  overlap_pct <- 100 * length(overlap_barcodes) / length(rna_barcodes)

  cat(sprintf("  RNA cells: %d\n", length(rna_barcodes)))
  cat(sprintf("  ATAC cells: %d\n", length(atac_barcodes)))
  cat(sprintf("  Overlap: %d (%.1f%% of RNA cells)\n",
              length(overlap_barcodes), overlap_pct))

  if (overlap_pct < 70) {
    cat("\n  ⚠ WARNING: Low barcode overlap (<70%)!\n")
    cat("    This suggests barcode translation may have failed.\n")
    cat("    Expected overlap: >80% for successful integration.\n")
    cat("\n    Troubleshooting:\n")
    cat("    1. Check that ATAC fragment files have correct format\n")
    cat("    2. Verify barcode translation function is working\n")
    cat("    3. Check for suffix mismatches (e.g., -1 vs -2)\n")
  } else if (overlap_pct < 90) {
    cat("\n  ⚠ Note: Moderate overlap (70-90%)\n")
    cat("    Some cells may be RNA-only or ATAC-only.\n")
    cat("    This is acceptable but verify it's expected for your data.\n")
  } else {
    cat("\n  ✓ Excellent barcode overlap (>90%)\n")
  }

  # Validate genotype distribution
  cat("\n--- Genotype Distribution Validation ---\n")

  if ("genotype" %in% colnames(seurat_obj@meta.data) &&
      "condition" %in% colnames(seurat_obj@meta.data)) {

    geno_cond_table <- table(seurat_obj$genotype, seurat_obj$condition)
    cat("\n  Genotype × Condition distribution:\n")
    print(geno_cond_table)

    # Check for balanced distribution
    total_cells <- sum(geno_cond_table)
    expected_per_group <- total_cells / 4

    imbalance_detected <- FALSE
    for (geno in rownames(geno_cond_table)) {
      for (cond in colnames(geno_cond_table)) {
        n_cells <- geno_cond_table[geno, cond]
        deviation_pct <- abs(n_cells - expected_per_group) / expected_per_group * 100

        if (deviation_pct > 50) {
          cat(sprintf("\n  ⚠ %s-%s has %.0f%% deviation from expected\n",
                     geno, cond, deviation_pct))
          imbalance_detected <- TRUE
        }
      }
    }

    if (!imbalance_detected) {
      cat("\n  ✓ Groups are reasonably balanced\n")
    }
  } else {
    cat("  ⚠ Cannot validate - genotype or condition column missing\n")
  }

  cat("\n  ✓ Validation checks complete\n")

  cat("\n================================================================================\n")
  cat("STEP 4: Computing QC Metrics\n")
  cat("================================================================================\n")

  # Compute ATAC QC metrics
  seurat_obj <- compute_atac_qc(seurat_obj)

  # Create QC plots
  seurat_obj <- create_qc_plots(seurat_obj, OUTPUT_DIR)

  cat("\n================================================================================\n")
  cat("STEP 5: Quality Control Filtering\n")
  cat("================================================================================\n")

  # Filter cells by QC
  seurat_obj <- filter_cells_by_qc(seurat_obj, QC_THRESHOLDS)

  cat("\n================================================================================\n")
  cat("STEP 6: Saving Results\n")
  cat("================================================================================\n")

  # Save Seurat object
  output_file <- file.path(OUTPUT_DIR, "integrated_seurat_raw.rds")
  cat(sprintf("  Saving integrated object to: %s\n", output_file))
  saveRDS(seurat_obj, file = output_file)

  # Save summary
  summary_df <- data.frame(
    assay = c("RNA", "ATAC"),
    n_features = c(nrow(seurat_obj[["RNA"]]), nrow(seurat_obj[["ATAC"]])),
    n_cells_total = c(ncol(seurat_obj), length(Cells(seurat_obj[["ATAC"]])))
  )

  write.csv(summary_df, file.path(OUTPUT_DIR, "integration_summary.csv"),
            row.names = FALSE)

  cat("\n✓ Step 1 complete!\n")
  cat(sprintf("  Output saved to: %s\n", OUTPUT_DIR))
  cat("\nNext step: Run signac_02_process_and_reduce.R\n\n")
}

# Run main function
if (!interactive()) {
  main()
}
