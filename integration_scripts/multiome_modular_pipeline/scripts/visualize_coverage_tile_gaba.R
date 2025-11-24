#!/usr/bin/env Rscript

################################################################################
# Coverage and Tile Plot Visualization for GABA Neurons (Combined Controls)
#
# This script creates clean visualizations without peaks:
# - CoveragePlot (peaks=FALSE): Clean coverage tracks + gene annotation + RNA
# - TilePlot: Per-cell fragment patterns (top 100 cells per group)
#
# Purpose:
# - Exploration: See clean coverage patterns
# - Publication: Generate cleaner figures without peak clutter
# - QC: Assess single-cell heterogeneity with TilePlot
# - Comparison: Compare patterns between conditions
#
# Input:
#   - genes_inter.txt (gene list for GABA neurons)
#   - integrated_seurat_processed.rds
#
# Output:
#   - Coverage plots without peaks (cleaner view)
#   - Tile plots showing individual cell fragment patterns
#
# INPUT FILES:
#   - Integrated Seurat object with processed data:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/integrated_seurat_processed.rds
#   - Gene list file (GABA-specific genes):
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/genes_inter.txt
#
# OUTPUT FILES:
#   - Coverage plots without peaks (clean view):
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/visualizations/GABA/coverage_no_peaks/<gene>_coverage_clean.pdf
#   - Tile plots for each group:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/visualizations/GABA/tile_plots/<gene>_tile_<group>.pdf
#   - Summary CSV with processing status:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/visualizations/GABA/visualization_summary.csv
#
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(EnsDb.Mmusculus.v79)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# Configuration
RESULTS_DIR <- "signac_results_L1"
SEURAT_FILE <- file.path(RESULTS_DIR, "integrated_seurat_processed.rds")
GENE_LIST_FILE <- "../genes_inter.txt"
OUTPUT_DIR <- file.path(RESULTS_DIR, "visualizations", "GABA")

# Cell type and condition columns
CELLTYPE_COL <- "cell_type_L1"
CELL_TYPE <- "GABA"
CONDITION_COL <- "condition"
GENOTYPE_COL <- "genotype"

# Plotting parameters
EXTEND_UPSTREAM <- 50000      # 50kb upstream
EXTEND_DOWNSTREAM <- 50000    # 50kb downstream
TILE_CELLS <- 100            # Top N cells for TilePlot

# Groups to visualize
GROUPS <- c("Combined_Ctrl", "Nestin_Mut", "Emx1_Mut")

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "coverage_no_peaks"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tile_plots"), recursive = TRUE, showWarnings = FALSE)

cat("================================================================================\n")
cat("COVERAGE AND TILE PLOT VISUALIZATION: GABA NEURONS (COMBINED CONTROLS)\n")
cat("================================================================================\n\n")
cat("Visualization types:\n")
cat("  1. CoveragePlot (peaks=FALSE): Clean coverage + gene annotation + RNA\n")
cat("  2. TilePlot: Per-cell fragment patterns (top 100 cells per group)\n\n")

################################################################################
# Load Data
################################################################################

cat("Loading data...\n")
seurat_obj <- readRDS(SEURAT_FILE)
cat(sprintf("  ✓ Loaded Seurat object with %d cells\n", ncol(seurat_obj)))

cat("\nLoading gene list...\n")
genes_of_interest <- unique(read.table(GENE_LIST_FILE, header = FALSE, stringsAsFactors = FALSE)$V1)
cat(sprintf("  ✓ Loaded %d genes\n", length(genes_of_interest)))

# Filter out mitochondrial genes
mt_genes <- genes_of_interest[grepl("^mt-", genes_of_interest, ignore.case = TRUE)]
if (length(mt_genes) > 0) {
  cat(sprintf("\n  ⚠ Filtering %d mitochondrial genes (incompatible with coverage plots):\n", length(mt_genes)))
  cat(sprintf("    %s\n", paste(mt_genes, collapse = ", ")))
  genes_of_interest <- genes_of_interest[!grepl("^mt-", genes_of_interest, ignore.case = TRUE)]
  cat(sprintf("  ✓ %d genes remaining\n", length(genes_of_interest)))
}

################################################################################
# Subset to GABA Neurons
################################################################################

cat("\nSubsetting to GABA neurons...\n")
gaba_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[CELLTYPE_COL]] == CELL_TYPE, ])
seurat_gaba <- subset(seurat_obj, cells = gaba_cells)
cat(sprintf("  ✓ Selected %d GABA neurons\n", ncol(seurat_gaba)))

# Create combined grouping
seurat_gaba@meta.data$combined_group <- ifelse(
  seurat_gaba@meta.data[[CONDITION_COL]] == "Ctrl",
  "Combined_Ctrl",
  paste(seurat_gaba@meta.data[[GENOTYPE_COL]], "Mut", sep = "_")
)

group_dist <- table(seurat_gaba@meta.data$combined_group)
cat("\n  Group distribution:\n")
print(group_dist)

################################################################################
# Visualize Each Gene
################################################################################

cat("\n================================================================================\n")
cat("CREATING VISUALIZATIONS\n")
cat("================================================================================\n\n")

# Summary tracking
summary_stats <- data.frame(
  gene = character(),
  status = character(),
  coverage_plot_created = logical(),
  tile_plots_created = integer(),
  stringsAsFactors = FALSE
)

for (gene in genes_of_interest) {
  cat(sprintf("\n--- Gene: %s ---\n", gene))

  # Validate gene exists in ATAC annotation
  gene_in_atac <- tryCatch({
    coords <- LookupGeneCoords(seurat_gaba, gene = gene, assay = "ATAC")
    !is.null(coords) && !any(is.na(coords)) && length(coords) > 0
  }, error = function(e) FALSE)

  if (!gene_in_atac) {
    cat("  ⚠ Skipping - gene not found in ATAC annotation\n")
    summary_stats <- rbind(summary_stats, data.frame(
      gene = gene,
      status = "missing_annotation",
      coverage_plot_created = FALSE,
      tile_plots_created = 0
    ))
    next
  }

  # Check if gene exists in RNA assay
  if (!gene %in% rownames(seurat_gaba[["RNA"]])) {
    cat("  ⚠ Skipping - gene not in RNA assay\n")
    summary_stats <- rbind(summary_stats, data.frame(
      gene = gene,
      status = "missing_rna",
      coverage_plot_created = FALSE,
      tile_plots_created = 0
    ))
    next
  }

  # Create CoveragePlot (peaks=FALSE)
  cat("  Creating CoveragePlot (peaks=FALSE)...\n")
  coverage_success <- FALSE

  tryCatch({
    DefaultAssay(seurat_gaba) <- "ATAC"

    p_coverage <- CoveragePlot(
      object = seurat_gaba,
      region = gene,
      features = gene,              # Show RNA expression
      expression.assay = "RNA",
      assay = "ATAC",
      group.by = "combined_group",
      extend.upstream = EXTEND_UPSTREAM,
      extend.downstream = EXTEND_DOWNSTREAM,
      annotation = TRUE,            # Show gene annotation
      peaks = FALSE,                # NO peaks - clean view
      links = FALSE,                # NO links
      window = 100
    )

    # Add title
    plot_title <- sprintf("GABA Neurons - %s\nClean Coverage View (no peaks)", gene)

    # Save
    pdf_file <- file.path(OUTPUT_DIR, "coverage_no_peaks",
                         sprintf("%s_coverage_clean.pdf", gene))
    ggsave(pdf_file, plot = p_coverage + ggtitle(plot_title),
           width = 16, height = 10, device = "pdf")

    cat(sprintf("    ✓ Saved coverage plot: %s\n", basename(pdf_file)))
    coverage_success <- TRUE

  }, error = function(e) {
    cat(sprintf("    ✗ Error creating coverage plot: %s\n", e$message))
  })

  # Create TilePlots for each group
  cat("  Creating TilePlots (per-cell fragments)...\n")
  tile_count <- 0

  # Set active identity to combined_group for TilePlot filtering
  Idents(seurat_gaba) <- "combined_group"

  for (group in GROUPS) {
    cat(sprintf("    Group: %s...", group))

    tryCatch({
      DefaultAssay(seurat_gaba) <- "ATAC"

      p_tile <- TilePlot(
        object = seurat_gaba,
        region = gene,
        idents = group,
        extend.upstream = EXTEND_UPSTREAM,
        extend.downstream = EXTEND_DOWNSTREAM
      )

      # Add title
      plot_title <- sprintf("GABA - %s - %s\nPer-Cell Fragment Patterns (top %d cells)",
                           gene, group, TILE_CELLS)

      # Save
      pdf_file <- file.path(OUTPUT_DIR, "tile_plots",
                           sprintf("%s_tile_%s.pdf", gene, group))
      ggsave(pdf_file, plot = p_tile + ggtitle(plot_title),
             width = 14, height = 10, device = "pdf")

      cat(" ✓\n")
      tile_count <- tile_count + 1

    }, error = function(e) {
      cat(sprintf(" ✗ Error: %s\n", e$message))
    })
  }

  cat(sprintf("    Created %d/%d tile plots\n", tile_count, length(GROUPS)))

  # Update summary
  summary_stats <- rbind(summary_stats, data.frame(
    gene = gene,
    status = ifelse(coverage_success && tile_count > 0, "success",
                   ifelse(coverage_success, "coverage_only",
                         ifelse(tile_count > 0, "tile_only", "failed"))),
    coverage_plot_created = coverage_success,
    tile_plots_created = tile_count
  ))
}

################################################################################
# Save Summary
################################################################################

cat("\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n\n")

write.csv(summary_stats,
          file.path(OUTPUT_DIR, "visualization_summary.csv"),
          row.names = FALSE)

cat("Summary statistics:\n")
print(summary_stats)

# Status breakdown
cat("\n\nStatus breakdown:\n")
status_counts <- table(summary_stats$status)
print(status_counts)

cat(sprintf("\n✓ Visualization complete!\n"))
cat(sprintf("  Total genes processed: %d\n", nrow(summary_stats)))
cat(sprintf("  Coverage plots created: %d\n", sum(summary_stats$coverage_plot_created)))
cat(sprintf("  Tile plots created: %d\n", sum(summary_stats$tile_plots_created)))
cat(sprintf("  Fully successful: %d\n", sum(summary_stats$status == "success")))

cat(sprintf("\nResults saved to: %s\n", OUTPUT_DIR))
cat("\nOutput structure:\n")
cat(sprintf("  - Coverage (no peaks): %s/coverage_no_peaks/\n", OUTPUT_DIR))
cat(sprintf("  - Tile plots: %s/tile_plots/\n", OUTPUT_DIR))
cat(sprintf("  - Summary: %s/visualization_summary.csv\n", OUTPUT_DIR))
