#!/usr/bin/env Rscript

################################################################################
# Plot Multiome Heatmaps - Genotype-Stratified Version
################################################################################
#
# Creates separate heatmaps for Nestin and Emx1 genotypes showing:
# - RNA expression changes (DEGs)
# - ATAC gene activity changes (chromatin accessibility)
# - Side-by-side comparison for concordant chromatin-driven genes
#
#
# Input Files:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/multiome_figure/multiome_heatmap_data_{Genotype}.rds
#     Processed data for each genotype containing RNA and ATAC matrices with z-scores
#     (Genotype = Nestin or Emx1)
#
# Output Files:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/multiome_figure/output_genotypes/multiome_heatmap_{Genotype}.pdf
#     Combined RNA+ATAC heatmaps for each genotype (L2 cell types)
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/multiome_figure/output_genotypes/multiome_heatmap_{Genotype}.png
#     PNG versions for presentations
#
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(tidyverse)
})

cat("================================================================================\n")
cat("  PLOT MULTIOME HEATMAPS: GENOTYPE-STRATIFIED VISUALIZATION\n")
cat("================================================================================\n\n")

################################################################################
# Configuration
################################################################################

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
FIGURE_DIR <- file.path(BASE_DIR, "multiome_figure")
OUTPUT_DIR <- file.path(FIGURE_DIR, "output_genotypes")

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Plotting parameters
HEATMAP_WIDTH_CM <- 18
HEATMAP_HEIGHT_CM <- 25
FONT_SIZE_ROW <- 8
FONT_SIZE_COL <- 10
FONT_SIZE_TITLE <- 12

# Color scales
color_scale_rna <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
color_scale_atac <- colorRamp2(c(-2, 0, 2), c("darkgreen", "white", "orange"))

################################################################################
# Function: Create Heatmap for Single Genotype
################################################################################

create_genotype_heatmap <- function(geno) {

  cat(sprintf("\n=== Creating heatmap for %s ===\n", geno))

  # Load data
  data_file <- file.path(FIGURE_DIR, sprintf("multiome_heatmap_data_%s.rds", geno))

  if (!file.exists(data_file)) {
    warning(sprintf("Data file not found: %s\nSkipping %s", data_file, geno))
    return(NULL)
  }

  data <- readRDS(data_file)

  rna_matrix <- data$rna_matrix_zscore
  atac_matrix <- data$atac_matrix_zscore
  selected_genes <- data$selected_genes
  top_celltypes <- data$top_celltypes

  cat(sprintf("  Loaded data: %d genes x %d conditions\n",
              nrow(rna_matrix), ncol(rna_matrix)))

  # Ensure matrices have same genes and columns
  common_genes <- intersect(rownames(rna_matrix), rownames(atac_matrix))
  common_cols <- intersect(colnames(rna_matrix), colnames(atac_matrix))

  if (length(common_genes) == 0 || length(common_cols) == 0) {
    warning(sprintf("No common genes/columns for %s", geno))
    return(NULL)
  }

  rna_matrix <- rna_matrix[common_genes, common_cols]
  atac_matrix <- atac_matrix[common_genes, common_cols]

  # Create row annotations (cell type)
  gene_celltype <- selected_genes$celltype[match(common_genes, selected_genes$gene)]
  row_annotation <- rowAnnotation(
    CellType = gene_celltype,
    col = list(CellType = setNames(
      rainbow(length(top_celltypes)),
      top_celltypes
    )),
    show_legend = TRUE,
    annotation_name_gp = gpar(fontsize = FONT_SIZE_COL)
  )

  # Create column split (cell type)
  col_celltype <- sapply(strsplit(colnames(rna_matrix), "_"), function(x) {
    paste(x[-length(x)], collapse = "_")
  })

  # Column annotations (Ctrl vs Mut)
  col_condition <- sapply(strsplit(colnames(rna_matrix), "_"), function(x) x[length(x)])
  col_annotation <- HeatmapAnnotation(
    Condition = col_condition,
    col = list(Condition = c("Ctrl" = "lightblue", "Mut" = "pink")),
    show_legend = TRUE,
    annotation_name_gp = gpar(fontsize = FONT_SIZE_COL)
  )

  # RNA heatmap
  ht_rna <- Heatmap(
    rna_matrix,
    name = "RNA\nZ-score",
    col = color_scale_rna,

    # Row settings
    cluster_rows = TRUE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = FONT_SIZE_ROW),
    row_split = gene_celltype,
    row_title_gp = gpar(fontsize = FONT_SIZE_COL),

    # Column settings
    cluster_columns = FALSE,
    column_split = col_celltype,
    column_title_gp = gpar(fontsize = FONT_SIZE_COL),
    column_names_gp = gpar(fontsize = FONT_SIZE_COL - 2),
    column_names_rot = 45,

    # Annotations
    left_annotation = row_annotation,
    top_annotation = col_annotation,

    # Title
    column_title = sprintf("%s - RNA Expression (Chromatin-Driven DEGs)", geno),
    column_title_gp = gpar(fontsize = FONT_SIZE_TITLE, fontface = "bold"),

    # Heatmap body
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = FONT_SIZE_COL),
      labels_gp = gpar(fontsize = FONT_SIZE_COL - 2)
    )
  )

  # ATAC heatmap
  ht_atac <- Heatmap(
    atac_matrix,
    name = "ATAC\nZ-score",
    col = color_scale_atac,

    # Row settings
    cluster_rows = TRUE,
    show_row_names = FALSE,  # Don't duplicate gene names
    row_split = gene_celltype,
    row_title_gp = gpar(fontsize = FONT_SIZE_COL),

    # Column settings
    cluster_columns = FALSE,
    column_split = col_celltype,
    column_title_gp = gpar(fontsize = FONT_SIZE_COL),
    column_names_gp = gpar(fontsize = FONT_SIZE_COL - 2),
    column_names_rot = 45,

    # Title
    column_title = sprintf("%s - ATAC Gene Activity", geno),
    column_title_gp = gpar(fontsize = FONT_SIZE_TITLE, fontface = "bold"),

    # Heatmap body
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = FONT_SIZE_COL),
      labels_gp = gpar(fontsize = FONT_SIZE_COL - 2)
    )
  )

  # Combine heatmaps
  ht_list <- ht_rna + ht_atac

  # Save to PDF
  output_file <- file.path(OUTPUT_DIR, sprintf("multiome_heatmap_%s.pdf", geno))

  pdf(output_file, width = HEATMAP_WIDTH_CM / 2.54, height = HEATMAP_HEIGHT_CM / 2.54)
  draw(ht_list,
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom",
       merge_legend = TRUE)
  dev.off()

  cat(sprintf("  ✓ Saved heatmap to: %s\n", output_file))

  # Also save as PNG for presentations
  output_png <- file.path(OUTPUT_DIR, sprintf("multiome_heatmap_%s.png", geno))
  png(output_png, width = HEATMAP_WIDTH_CM * 100, height = HEATMAP_HEIGHT_CM * 100, res = 300)
  draw(ht_list,
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom",
       merge_legend = TRUE)
  dev.off()

  cat(sprintf("  ✓ Saved PNG to: %s\n", output_png))

  return(ht_list)
}

################################################################################
# Main Execution
################################################################################

# Detect available genotypes
genotype_files <- list.files(FIGURE_DIR, pattern = "^multiome_heatmap_data_.*\\.rds$", full.names = FALSE)
genotypes <- gsub("multiome_heatmap_data_(.*)\\.rds", "\\1", genotype_files)

if (length(genotypes) == 0) {
  stop("No data files found. Please run prepare_multiome_heatmap_data_genotypes.R first.")
}

cat(sprintf("Found data for %d genotype(s): %s\n\n", length(genotypes), paste(genotypes, collapse = ", ")))

# Create heatmaps for each genotype
results <- list()
for (geno in genotypes) {
  results[[geno]] <- create_genotype_heatmap(geno)
}

################################################################################
# Summary
################################################################################

cat("\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n\n")

for (geno in genotypes) {
  if (!is.null(results[[geno]])) {
    cat(sprintf("✓ %s: Heatmap created successfully\n", geno))
  } else {
    cat(sprintf("✗ %s: Failed to create heatmap\n", geno))
  }
}

cat(sprintf("\nOutput directory: %s\n", OUTPUT_DIR))
cat("\nFiles created:\n")
cat("  - multiome_heatmap_<genotype>.pdf\n")
cat("  - multiome_heatmap_<genotype>.png\n")
cat("\n✓ Visualization complete!\n\n")
