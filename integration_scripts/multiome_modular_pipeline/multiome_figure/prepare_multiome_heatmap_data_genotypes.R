#!/usr/bin/env Rscript

################################################################################
# Prepare Multiome Heatmap Data - Genotype-Stratified Version
################################################################################
#
# This script integrates DEG, DA, and peak-gene linkage results to identify
# high-confidence chromatin-driven DEGs for BOTH Nestin and Emx1 genotypes.
#
# Strategy:
#   1. Load Seurat object with RNA and ATAC data
#   2. For each GENOTYPE (Nestin, Emx1):
#      - For each cell type, identify "chromatin-driven DEGs":
#        * DEG with significant expression change
#        * Linked to DA peak with CONCORDANT change (same direction)
#        * Meaningful correlation (|r| > threshold)
#   3. Select TOP cell types by abundance (automatic detection)
#   4. Extract normalized RNA expression and ATAC gene activity
#   5. Calculate per-cell-type, per-condition averages
#   6. Save structured data for visualization (SEPARATE for each genotype)
#
# File Detection Strategy:
#   - First tries genotype-specific files: {CellType}_{Genotype}_DEGs.csv
#   - Falls back to combined files if genotype-specific not available:
#     {CellType}_DEG.csv (contains data from both genotypes)
#   - NOTE: Using combined files means results for both genotypes will be
#     identical for that cell type (no genotype stratification)
#
# Output: Separate .rds files for Nestin and Emx1
#
#
# Input Files:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/integrated_seurat_processed.rds
#     Integrated Seurat object with processed RNA and ATAC data (L2 annotations)
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/celltype_results/DEG/{CellType}_{Genotype}_DEGs.csv
#     Differential expression results for each cell type and genotype
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/celltype_results/DA/{CellType}_{Genotype}_DA_peaks.csv
#     Differential accessibility results for each cell type and genotype
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/celltype_results/peak_gene_links/{CellType}_{Genotype}_peak_gene_links.csv
#     Peak-gene linkage results with correlation scores
#
# Output Files:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/multiome_figure/multiome_heatmap_data_{Genotype}.rds
#     Processed data matrices for each genotype containing RNA and ATAC z-scores
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/multiome_figure/output_genotypes/chromatin_driven_degs_{Genotype}.csv
#     Summary of selected chromatin-driven DEGs for each genotype
#
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(tidyverse)
  library(Matrix)
})

cat("================================================================================\n")
cat("  PREPARE MULTIOME HEATMAP DATA: GENOTYPE-STRATIFIED ANALYSIS\n")
cat("================================================================================\n\n")

################################################################################
# Configuration
################################################################################

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
RESULTS_DIR <- file.path(BASE_DIR, "signac_results/celltype_results")
FIGURE_DIR <- file.path(BASE_DIR, "multiome_figure")
OUTPUT_DIR <- file.path(FIGURE_DIR, "output_genotypes")

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", OUTPUT_DIR))
}

# Analysis parameters
CELLTYPE_COL <- "cell_type_L2_new"
CONDITION_COL <- "condition"
GENOTYPE_COL <- "genotype"

# Selection criteria for chromatin-driven DEGs (CONCORDANT ONLY)
MIN_GENE_LOG2FC <- 0.5     # Minimum gene expression change
MIN_PEAK_LOG2FC <- 0.5     # Minimum peak accessibility change
MIN_CORRELATION <- 0.2     # Minimum |peak-gene correlation|
MAX_PADJ <- 0.05           # Maximum adjusted p-value

# Number of top cell types to include (by abundance)
TOP_N_CELLTYPES <- 8

# Genes to select per cell type
GENES_PER_CELLTYPE <- 15

cat("Configuration:\n")
cat(sprintf("  Min gene log2FC: %.2f\n", MIN_GENE_LOG2FC))
cat(sprintf("  Min peak log2FC: %.2f\n", MIN_PEAK_LOG2FC))
cat(sprintf("  Min |correlation|: %.2f\n", MIN_CORRELATION))
cat(sprintf("  Max adjusted p-value: %.2e\n", MAX_PADJ))
cat(sprintf("  Top N cell types: %d\n", TOP_N_CELLTYPES))
cat(sprintf("  Genes per cell type: %d\n", GENES_PER_CELLTYPE))
cat(sprintf("  Filter mode: CONCORDANT ONLY\n"))
cat("\n")

################################################################################
# Load Seurat Object
################################################################################

cat("=== Loading Seurat object ===\n")
seurat_file <- file.path(BASE_DIR, "signac_results/integrated_seurat_processed.rds")

if (!file.exists(seurat_file)) {
  stop(sprintf("Seurat object not found: %s\nPlease run signac_02_process_and_reduce.R first.", seurat_file))
}

seurat_obj <- readRDS(seurat_file)
cat(sprintf("  Loaded: %d cells, %d RNA features, %d ATAC features\n",
            ncol(seurat_obj),
            nrow(seurat_obj[["RNA"]]),
            nrow(seurat_obj[["ATAC"]])))

# Check required columns
required_cols <- c(CELLTYPE_COL, CONDITION_COL, GENOTYPE_COL)
missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required metadata columns: %s", paste(missing_cols, collapse = ", ")))
}

# Get genotypes
genotypes <- unique(seurat_obj@meta.data[[GENOTYPE_COL]])
genotypes <- genotypes[!is.na(genotypes)]
cat(sprintf("  Found %d genotypes: %s\n", length(genotypes), paste(genotypes, collapse = ", ")))

# Compute gene activity from ATAC BEFORE subsetting (to preserve ChromatinAssay structure)
if (!"GeneActivity" %in% Assays(seurat_obj)) {
  cat("\n=== Computing gene activity from ATAC (full object) ===\n")
  cat("  This preserves the ChromatinAssay structure...\n")
  gene.activities <- GeneActivity(seurat_obj, assay = "ATAC")
  seurat_obj[["GeneActivity"]] <- CreateAssayObject(counts = gene.activities)
  seurat_obj <- NormalizeData(
    object = seurat_obj,
    assay = "GeneActivity",
    normalization.method = "LogNormalize",
    scale.factor = median(seurat_obj$nCount_ATAC)
  )
  cat("  ✓ Gene activity computed and normalized\n")
} else {
  cat("\n=== Gene activity assay already exists ===\n")
}

################################################################################
# Determine Top Cell Types by Abundance
################################################################################

cat("\n=== Determining top cell types by abundance ===\n")

# Count cells per cell type (across all genotypes and conditions)
celltype_counts <- table(seurat_obj@meta.data[[CELLTYPE_COL]])
celltype_counts <- sort(celltype_counts, decreasing = TRUE)

cat(sprintf("  Total unique cell types: %d\n", length(celltype_counts)))
cat(sprintf("  Selecting top %d most abundant cell types:\n", TOP_N_CELLTYPES))

top_celltypes <- names(celltype_counts)[1:min(TOP_N_CELLTYPES, length(celltype_counts))]

for (i in 1:length(top_celltypes)) {
  ct <- top_celltypes[i]
  count <- celltype_counts[ct]
  cat(sprintf("    %d. %s: %d cells\n", i, ct, count))
}

cat("\n")

################################################################################
# Function: Process Single Genotype
################################################################################

process_genotype <- function(geno) {

  cat("\n================================================================================\n")
  cat(sprintf("PROCESSING GENOTYPE: %s\n", geno))
  cat("================================================================================\n\n")

  # Initialize storage
  chromatin_driven_degs_all <- data.frame()
  file_type_used <- list()  # Track which file type was used for each cell type

  # Process each top cell type
  for (ct in top_celltypes) {

    cat(sprintf("--- Cell type: %s (%s) ---\n", ct, geno))

    # Construct file paths for this genotype + cell type
    analysis_id <- sprintf("%s_%s", ct, geno)
    deg_file <- file.path(RESULTS_DIR, "DEG", sprintf("%s_DEGs.csv", analysis_id))
    da_file <- file.path(RESULTS_DIR, "DA", sprintf("%s_DA_peaks.csv", analysis_id))
    links_file <- file.path(RESULTS_DIR, "peak_gene_links", sprintf("%s_peak_gene_links.csv", analysis_id))

    # Check if genotype-specific files exist, otherwise try combined files
    files_found <- c(deg = file.exists(deg_file),
                     da = file.exists(da_file),
                     links = file.exists(links_file))

    if (!all(files_found)) {
      # Try combined (legacy) files as fallback
      cat(sprintf("  ⚠ Genotype-specific files not found, trying combined files...\n"))
      deg_file_combined <- file.path(RESULTS_DIR, "DEG", sprintf("%s_DEG.csv", ct))
      da_file_combined <- file.path(RESULTS_DIR, "DA", sprintf("%s_DA.csv", ct))
      links_file_combined <- file.path(RESULTS_DIR, "peak_gene_links", sprintf("%s_peak_gene_links.csv", ct))

      files_found_combined <- c(deg = file.exists(deg_file_combined),
                                da = file.exists(da_file_combined),
                                links = file.exists(links_file_combined))

      if (all(files_found_combined)) {
        deg_file <- deg_file_combined
        da_file <- da_file_combined
        links_file <- links_file_combined
        cat(sprintf("  ✓ Using combined files (legacy format)\n"))
        cat(sprintf("  ⚠ WARNING: Combined files contain data from BOTH genotypes!\n"))
        file_type_used[[ct]] <- "combined"
      } else {
        cat(sprintf("  ⚠ Skipping (neither genotype-specific nor combined files available)\n"))
        cat(sprintf("      Missing genotype-specific: %s\n",
                    paste(names(files_found)[!files_found], collapse = ", ")))
        cat(sprintf("      Missing combined: %s\n",
                    paste(names(files_found_combined)[!files_found_combined], collapse = ", ")))
        next
      }
    } else {
      cat(sprintf("  ✓ Using genotype-specific files\n"))
      file_type_used[[ct]] <- "genotype-specific"
    }

    # Load results
    deg_res <- read.csv(deg_file, stringsAsFactors = FALSE)
    da_res <- read.csv(da_file, stringsAsFactors = FALSE)
    links_res <- read.csv(links_file, stringsAsFactors = FALSE)

    cat(sprintf("  Loaded: %d DEGs, %d DA peaks, %d links\n",
                nrow(deg_res), nrow(da_res), nrow(links_res)))

    # Filter DEGs by significance and log2FC
    deg_sig <- deg_res %>%
      filter(!is.na(p_val_adj) & p_val_adj < MAX_PADJ & abs(avg_log2FC) >= MIN_GENE_LOG2FC)

    if (nrow(deg_sig) == 0) {
      cat("  No significant DEGs pass filters\n")
      next
    }

    # Filter DA peaks
    da_sig <- da_res %>%
      filter(!is.na(p_val_adj) & p_val_adj < MAX_PADJ & abs(avg_log2FC) >= MIN_PEAK_LOG2FC)

    # Filter links by correlation
    links_filt <- links_res %>%
      filter(!is.na(score) & abs(score) >= MIN_CORRELATION)

    if (nrow(links_filt) == 0) {
      cat("  No peak-gene links pass correlation filter\n")
      next
    }

    # Join DEG with links
    # Links should have columns: gene, peak, score (correlation)
    deg_with_links <- deg_sig %>%
      inner_join(links_filt, by = "gene")

    if (nrow(deg_with_links) == 0) {
      cat("  No DEGs have linked peaks\n")
      next
    }

    # Join with DA peaks
    # Need to match peak names
    # DA peaks have column "peak", links have column "peak"
    deg_with_da <- deg_with_links %>%
      inner_join(da_sig %>% select(peak, da_log2FC = avg_log2FC, da_padj = p_val_adj),
                 by = "peak")

    if (nrow(deg_with_da) == 0) {
      cat("  No linked DEGs have DA peaks\n")
      next
    }

    # Filter for CONCORDANT changes only
    # Concordant: gene and peak change in same direction
    deg_with_da <- deg_with_da %>%
      mutate(
        gene_direction = sign(avg_log2FC),
        peak_direction = sign(da_log2FC),
        concordant = gene_direction == peak_direction
      ) %>%
      filter(concordant == TRUE)

    if (nrow(deg_with_da) == 0) {
      cat("  No concordant gene-peak pairs\n")
      next
    }

    cat(sprintf("  Found %d concordant chromatin-driven DEGs\n", nrow(deg_with_da)))

    # For each gene, select the peak with strongest correlation
    deg_best <- deg_with_da %>%
      group_by(gene) %>%
      slice_max(order_by = abs(score), n = 1, with_ties = FALSE) %>%
      ungroup()

    # Add cell type and genotype
    deg_best$celltype <- ct
    deg_best$genotype <- geno

    # Append to master list
    chromatin_driven_degs_all <- rbind(chromatin_driven_degs_all, deg_best)
  }

  # Select top genes per cell type
  cat("\n=== Selecting top genes per cell type ===\n")

  selected_genes_all <- data.frame()

  for (ct in top_celltypes) {
    ct_degs <- chromatin_driven_degs_all %>%
      filter(celltype == ct)

    if (nrow(ct_degs) == 0) {
      cat(sprintf("  %s: No chromatin-driven DEGs\n", ct))
      next
    }

    # Select top by absolute gene log2FC
    ct_top <- ct_degs %>%
      arrange(desc(abs(avg_log2FC))) %>%
      head(GENES_PER_CELLTYPE)

    cat(sprintf("  %s: Selected %d genes (from %d candidates)\n",
                ct, nrow(ct_top), nrow(ct_degs)))

    selected_genes_all <- rbind(selected_genes_all, ct_top)
  }

  if (nrow(selected_genes_all) == 0) {
    warning(sprintf("No chromatin-driven DEGs found for genotype %s", geno))
    return(NULL)
  }

  cat(sprintf("\nTotal genes selected for %s: %d\n", geno, nrow(selected_genes_all)))

  ################################################################################
  # Extract RNA and ATAC Signals
  ################################################################################

  cat("\n=== Extracting RNA and ATAC signals ===\n")

  # Subset Seurat object to this genotype
  cells_geno <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[GENOTYPE_COL]] == geno, ])
  seurat_geno <- subset(seurat_obj, cells = cells_geno)

  cat(sprintf("  Cells for %s: %d\n", geno, ncol(seurat_geno)))

  # Get unique genes
  genes_to_plot <- unique(selected_genes_all$gene)
  cat(sprintf("  Unique genes to extract: %d\n", length(genes_to_plot)))

  # Extract RNA expression (log-normalized)
  DefaultAssay(seurat_geno) <- "RNA"

  # Check which genes are available
  genes_available <- intersect(genes_to_plot, rownames(seurat_geno[["RNA"]]))
  genes_missing <- setdiff(genes_to_plot, genes_available)

  if (length(genes_missing) > 0) {
    cat(sprintf("  ⚠ Warning: %d genes not found in RNA assay\n", length(genes_missing)))
  }

  rna_expr <- GetAssayData(seurat_geno, assay = "RNA", layer = "data")[genes_available, , drop = FALSE]

  # Extract ATAC gene activity (already computed on full object)
  DefaultAssay(seurat_geno) <- "GeneActivity"
  genes_in_atac <- intersect(genes_available, rownames(seurat_geno[["GeneActivity"]]))
  atac_activity <- GetAssayData(seurat_geno, assay = "GeneActivity", layer = "data")[genes_in_atac, , drop = FALSE]

  cat(sprintf("  Extracted RNA: %d genes x %d cells\n", nrow(rna_expr), ncol(rna_expr)))
  cat(sprintf("  Extracted ATAC: %d genes x %d cells\n", nrow(atac_activity), ncol(atac_activity)))

  ################################################################################
  # Calculate Per-Celltype, Per-Condition Averages
  ################################################################################

  cat("\n=== Calculating averages per cell type and condition ===\n")

  # Create matrices for heatmap
  rna_matrix <- matrix(NA, nrow = length(genes_in_atac), ncol = length(top_celltypes) * 2,
                       dimnames = list(genes_in_atac, paste(rep(top_celltypes, each = 2),
                                                             rep(c("Ctrl", "Mut"), length(top_celltypes)),
                                                             sep = "_")))

  atac_matrix <- rna_matrix  # Same structure

  for (ct in top_celltypes) {
    for (cond in c("Ctrl", "Mut")) {

      # Get cells for this combination
      cells_subset <- rownames(seurat_geno@meta.data[
        seurat_geno@meta.data[[CELLTYPE_COL]] == ct &
          seurat_geno@meta.data[[CONDITION_COL]] == cond, ])

      if (length(cells_subset) < 5) {
        cat(sprintf("  ⚠ %s_%s: Too few cells (%d), skipping\n", ct, cond, length(cells_subset)))
        next
      }

      col_name <- paste(ct, cond, sep = "_")

      # Calculate averages
      rna_avg <- rowMeans(rna_expr[genes_in_atac, cells_subset, drop = FALSE])
      atac_avg <- rowMeans(atac_activity[genes_in_atac, cells_subset, drop = FALSE])

      rna_matrix[, col_name] <- rna_avg
      atac_matrix[, col_name] <- atac_avg
    }
  }

  # Remove rows/columns with all NA
  rna_matrix <- rna_matrix[rowSums(!is.na(rna_matrix)) > 0, colSums(!is.na(rna_matrix)) > 0, drop = FALSE]
  atac_matrix <- atac_matrix[rowSums(!is.na(atac_matrix)) > 0, colSums(!is.na(atac_matrix)) > 0, drop = FALSE]

  cat(sprintf("  Final matrix dimensions: %d genes x %d conditions\n",
              nrow(rna_matrix), ncol(rna_matrix)))

  ################################################################################
  # Z-score Normalization (within each gene)
  ################################################################################

  cat("\n=== Z-score normalizing ===\n")

  # Z-score normalize each row (gene)
  rna_matrix_z <- t(scale(t(rna_matrix)))
  atac_matrix_z <- t(scale(t(atac_matrix)))

  # Replace NaN with 0 (genes with no variance)
  rna_matrix_z[is.nan(rna_matrix_z)] <- 0
  atac_matrix_z[is.nan(atac_matrix_z)] <- 0

  ################################################################################
  # Package Results
  ################################################################################

  results <- list(
    genotype = geno,
    rna_matrix_raw = rna_matrix,
    atac_matrix_raw = atac_matrix,
    rna_matrix_zscore = rna_matrix_z,
    atac_matrix_zscore = atac_matrix_z,
    selected_genes = selected_genes_all,
    top_celltypes = top_celltypes,
    file_type_used = file_type_used,
    config = list(
      min_gene_log2fc = MIN_GENE_LOG2FC,
      min_peak_log2fc = MIN_PEAK_LOG2FC,
      min_correlation = MIN_CORRELATION,
      max_padj = MAX_PADJ,
      genes_per_celltype = GENES_PER_CELLTYPE,
      filter_mode = "concordant"
    )
  )

  return(results)
}

################################################################################
# Process Both Genotypes
################################################################################

all_results <- list()

for (geno in genotypes) {
  geno_results <- process_genotype(geno)

  if (!is.null(geno_results)) {
    all_results[[geno]] <- geno_results

    # Save individual RDS file
    output_file <- file.path(FIGURE_DIR, sprintf("multiome_heatmap_data_%s.rds", geno))
    saveRDS(geno_results, output_file)
    cat(sprintf("\n✓ Saved %s results to: %s\n", geno, output_file))

    # Save summary CSV
    summary_file <- file.path(OUTPUT_DIR, sprintf("chromatin_driven_degs_%s.csv", geno))
    write.csv(geno_results$selected_genes, summary_file, row.names = FALSE)
    cat(sprintf("✓ Saved %s gene list to: %s\n", geno, summary_file))
  }
}

################################################################################
# Summary
################################################################################

cat("\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n\n")

for (geno in names(all_results)) {
  res <- all_results[[geno]]
  cat(sprintf("%s:\n", geno))
  cat(sprintf("  Cell types analyzed: %d\n", length(res$top_celltypes)))
  cat(sprintf("  Total genes selected: %d\n", nrow(res$selected_genes)))
  cat(sprintf("  Matrix dimensions: %d genes x %d conditions\n",
              nrow(res$rna_matrix_zscore), ncol(res$rna_matrix_zscore)))

  # Report file types used
  n_genotype_specific <- sum(unlist(res$file_type_used) == "genotype-specific")
  n_combined <- sum(unlist(res$file_type_used) == "combined")

  if (length(res$file_type_used) > 0) {
    cat(sprintf("  File types used:\n"))
    cat(sprintf("    - Genotype-specific: %d cell types\n", n_genotype_specific))
    cat(sprintf("    - Combined (legacy): %d cell types\n", n_combined))

    if (n_combined > 0) {
      combined_cts <- names(res$file_type_used)[unlist(res$file_type_used) == "combined"]
      cat(sprintf("    ⚠ Cell types using combined files: %s\n",
                  paste(combined_cts, collapse = ", ")))
      cat(sprintf("      (These show identical results for both genotypes)\n"))
    }
  }
  cat("\n")
}

cat("✓ Data preparation complete!\n")
cat("\nNext step: Run plot_multiome_heatmap_genotypes.R to create visualizations\n\n")
