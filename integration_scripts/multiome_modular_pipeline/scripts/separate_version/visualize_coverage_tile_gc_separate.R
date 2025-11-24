#!/usr/bin/env Rscript

################################################################################
# Coverage and Tile Plot Visualization for Granule Cells (Separate Genotypes)
#
# This version uses 4 separate groups: Nestin-Ctrl, Nestin-Mut, Emx1-Ctrl, Emx1-Mut
# Cell types: Immature GC + Mature GC
#
# INPUT FILES:
#   - Integrated Seurat object with processed data:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/integrated_seurat_processed.rds
#   - Gene list file (GC-specific genes):
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/genes_gc.txt
#
# OUTPUT FILES:
#   - Coverage plots without peaks (clean view):
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/visualizations_separate/GC/coverage_no_peaks/<gene>_coverage_clean_separate.pdf
#   - Tile plots for each genotype group:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/visualizations_separate/GC/tile_plots/<gene>_tile_<group>.pdf
#   - Summary CSV with processing status:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/visualizations_separate/GC/visualization_summary_separate.csv
#
################################################################################

suppressPackageStartupMessages({
  library(Signac); library(Seurat); library(GenomicRanges)
  library(BSgenome.Mmusculus.UCSC.mm10); library(EnsDb.Mmusculus.v79)
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

# Configuration
RESULTS_DIR <- "signac_results_L1"
SEURAT_FILE <- file.path(RESULTS_DIR, "integrated_seurat_processed.rds")
GENE_LIST_FILE <- "../genes_gc.txt"
OUTPUT_DIR <- file.path(RESULTS_DIR, "visualizations_separate", "GC")

CELLTYPE_COL <- "cell_type_L2_new"
CELL_TYPES <- c("Immature GC", "Mature GC")
CONDITION_COL <- "condition"; GENOTYPE_COL <- "genotype"

EXTEND_UPSTREAM <- 50000; EXTEND_DOWNSTREAM <- 50000; TILE_CELLS <- 100

# 4 separate groups
GROUPS <- c("Nestin-Ctrl", "Nestin-Mut", "Emx1-Ctrl", "Emx1-Mut")

# Create directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "coverage_no_peaks"), recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "tile_plots"), recursive = TRUE)

cat("================================================================================\n")
cat("COVERAGE AND TILE VISUALIZATION: GRANULE CELLS (SEPARATE GENOTYPES)\n")
cat("================================================================================\n\n")

# Load data
seurat_obj <- readRDS(SEURAT_FILE)
genes_of_interest <- unique(read.table(GENE_LIST_FILE, header = FALSE)$V1)

# Filter mt-genes
mt_genes <- genes_of_interest[grepl("^mt-", genes_of_interest, ignore.case = TRUE)]
if (length(mt_genes) > 0) {
  cat(sprintf("⚠ Filtering %d mitochondrial genes\n", length(mt_genes)))
  genes_of_interest <- genes_of_interest[!grepl("^mt-", genes_of_interest, ignore.case = TRUE)]
}

cat(sprintf("✓ Loaded %d genes\n", length(genes_of_interest)))

# Subset to Granule Cells
gc_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[CELLTYPE_COL]] %in% CELL_TYPES, ])
seurat_gc <- subset(seurat_obj, cells = gc_cells)

# Create separate genotype grouping
seurat_gc@meta.data$genotype_condition <- paste(
  seurat_gc@meta.data[[GENOTYPE_COL]],
  seurat_gc@meta.data[[CONDITION_COL]], sep = "-"
)

cat(sprintf("✓ Selected %d granule cells\n", ncol(seurat_gc)))
print(table(seurat_gc@meta.data$genotype_condition))

# Summary tracking
summary_stats <- data.frame(
  gene = character(), status = character(),
  coverage_plot_created = logical(), tile_plots_created = integer(),
  stringsAsFactors = FALSE
)

cat("\n=== CREATING VISUALIZATIONS ===\n")

for (gene in genes_of_interest) {
  cat(sprintf("\n-- %s --\n", gene))

  # Validate
  gene_in_atac <- tryCatch({
    coords <- LookupGeneCoords(seurat_gc, gene = gene, assay = "ATAC")
    !is.null(coords) && !any(is.na(coords)) && length(coords) > 0
  }, error = function(e) FALSE)

  if (!gene_in_atac) {
    cat("  ⚠ Missing ATAC annotation\n")
    summary_stats <- rbind(summary_stats, data.frame(
      gene = gene, status = "missing_annotation",
      coverage_plot_created = FALSE, tile_plots_created = 0))
    next
  }

  if (!gene %in% rownames(seurat_gc[["RNA"]])) {
    cat("  ⚠ Missing RNA\n")
    summary_stats <- rbind(summary_stats, data.frame(
      gene = gene, status = "missing_rna",
      coverage_plot_created = FALSE, tile_plots_created = 0))
    next
  }

  # Coverage plot
  cat("  CoveragePlot...")
  coverage_success <- FALSE
  tryCatch({
    DefaultAssay(seurat_gc) <- "ATAC"
    p_coverage <- CoveragePlot(
      object = seurat_gc, region = gene, features = gene,
      expression.assay = "RNA", assay = "ATAC",
      group.by = "genotype_condition",
      extend.upstream = EXTEND_UPSTREAM,
      extend.downstream = EXTEND_DOWNSTREAM,
      annotation = TRUE, peaks = FALSE, links = FALSE, window = 100
    )
    ggsave(file.path(OUTPUT_DIR, "coverage_no_peaks", sprintf("%s_coverage_clean_separate.pdf", gene)),
           plot = p_coverage + ggtitle(sprintf("Granule Cells - %s (Separate Genotypes)", gene)),
           width = 16, height = 12, device = "pdf")
    cat(" ✓\n")
    coverage_success <- TRUE
  }, error = function(e) cat(sprintf(" ✗ %s\n", e$message)))

  # Tile plots for each group
  cat("  TilePlots: ")
  tile_count <- 0

  # Set active identity to genotype_condition for TilePlot filtering
  Idents(seurat_gc) <- "genotype_condition"

  for (group in GROUPS) {
    tryCatch({
      DefaultAssay(seurat_gc) <- "ATAC"
      p_tile <- TilePlot(object = seurat_gc, region = gene, idents = group,
                        extend.upstream = EXTEND_UPSTREAM,
                        extend.downstream = EXTEND_DOWNSTREAM)
      ggsave(file.path(OUTPUT_DIR, "tile_plots", sprintf("%s_tile_%s.pdf", gene, group)),
             plot = p_tile + ggtitle(sprintf("GC - %s - %s", gene, group)),
             width = 14, height = 10, device = "pdf")
      cat(".")
      tile_count <- tile_count + 1
    }, error = function(e) cat("x"))
  }
  cat(sprintf(" %d/%d ✓\n", tile_count, length(GROUPS)))

  summary_stats <- rbind(summary_stats, data.frame(
    gene = gene,
    status = ifelse(coverage_success && tile_count > 0, "success",
                   ifelse(coverage_success, "coverage_only",
                         ifelse(tile_count > 0, "tile_only", "failed"))),
    coverage_plot_created = coverage_success,
    tile_plots_created = tile_count))
}

# Save summary
write.csv(summary_stats, file.path(OUTPUT_DIR, "visualization_summary_separate.csv"), row.names = FALSE)

cat("\n✓ Complete! Results in:", OUTPUT_DIR, "\n")
cat(sprintf("  Coverage plots: %d\n", sum(summary_stats$coverage_plot_created)))
cat(sprintf("  Tile plots: %d\n", sum(summary_stats$tile_plots_created)))
cat(sprintf("  Fully successful: %d\n", sum(summary_stats$status == "success")))
