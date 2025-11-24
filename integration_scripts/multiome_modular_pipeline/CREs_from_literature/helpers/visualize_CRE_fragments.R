#!/usr/bin/env Rscript

################################################################################
# Visualize ATAC Fragment Coverage at Literature CRE Locations
#
# Creates coverage plots directly from fragment files, not requiring peaks
# Shows raw ATAC signal to see if there's accessibility even without peak calls
#
# INPUT FILES:
# - ../signac_results_L1/integrated_seurat_processed.rds: Seurat object with ATAC-seq fragment data
# - output/unique_CREs_for_genes_of_interest_with_header.bed: BED file with literature CRE coordinates
#
# OUTPUT FILES:
# - output/CRE_fragment_coverage/*_separate_genotypes.pdf: Individual CRE plots with 4 genotype groups
# - output/CRE_fragment_coverage/*_combined_controls.pdf: Individual CRE plots with combined controls
# - output/CRE_fragment_coverage/*_zoomed.pdf: Zoomed views (±10kb) around each CRE
# - output/CRE_fragment_coverage/all_CREs_summary.pdf: Multi-panel summary of all CREs
#
################################################################################

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(EnsDb.Mmusculus.v79)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

# Configuration
RESULTS_DIR <- "../signac_results_L1"
SEURAT_FILE <- file.path(RESULTS_DIR, "integrated_seurat_processed.rds")
LITERATURE_CRE_FILE <- "output/unique_CREs_for_genes_of_interest_with_header.bed"
OUTPUT_DIR <- "output/CRE_fragment_coverage"

CELL_TYPE <- "GABA"
CELLTYPE_COL <- "cell_type_L1"
CONDITION_COL <- "condition"
GENOTYPE_COL <- "genotype"

# Plotting parameters
EXTEND_UPSTREAM <- 50000      # 50kb
EXTEND_DOWNSTREAM <- 50000    # 50kb

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("================================================================================\n")
cat("VISUALIZE ATAC FRAGMENT COVERAGE AT LITERATURE CRE LOCATIONS\n")
cat("================================================================================\n\n")

# ============================================================================
# STEP 1: Load Seurat object
# ============================================================================
cat("STEP 1: Loading Seurat object...\n")
cat(strrep("-", 80), "\n")

seurat_obj <- readRDS(SEURAT_FILE)
cat(sprintf("✓ Loaded Seurat object (%d cells)\n", ncol(seurat_obj)))

# Filter for GABA cells
seurat_obj <- subset(seurat_obj, subset = !!sym(CELLTYPE_COL) == CELL_TYPE)
cat(sprintf("✓ Filtered to %d GABA cells\n", ncol(seurat_obj)))

# Create grouping variables
seurat_obj$plot_group <- paste0(
  seurat_obj@meta.data[[GENOTYPE_COL]], "-",
  seurat_obj@meta.data[[CONDITION_COL]]
)

seurat_obj$plot_group_combined <- ifelse(
  seurat_obj@meta.data[[CONDITION_COL]] == "Ctrl",
  "Combined-Ctrl",
  seurat_obj$plot_group
)

cat("\nGroup distribution:\n")
table_separate <- table(seurat_obj$plot_group)
cat(sprintf("  Nestin-Ctrl: %d\n", table_separate["Nestin-Ctrl"]))
cat(sprintf("  Nestin-Mut: %d\n", table_separate["Nestin-Mut"]))
cat(sprintf("  Emx1-Ctrl: %d\n", table_separate["Emx1-Ctrl"]))
cat(sprintf("  Emx1-Mut: %d\n", table_separate["Emx1-Mut"]))

table_combined <- table(seurat_obj$plot_group_combined)
cat("\nCombined groups:\n")
cat(sprintf("  Combined-Ctrl: %d\n", table_combined["Combined-Ctrl"]))
cat(sprintf("  Nestin-Mut: %d\n", table_combined["Nestin-Mut"]))
cat(sprintf("  Emx1-Mut: %d\n", table_combined["Emx1-Mut"]))

# ============================================================================
# STEP 2: Load literature CREs
# ============================================================================
cat("\n")
cat(strrep("=", 80), "\n")
cat("STEP 2: Loading literature CREs...\n")
cat(strrep("-", 80), "\n")

lit_cres <- read.table(LITERATURE_CRE_FILE, header=TRUE, sep="\t", stringsAsFactors=FALSE)
cat(sprintf("✓ Loaded %d literature CREs\n", nrow(lit_cres)))

# Convert to GRanges
cre_ranges <- GRanges(
  seqnames = lit_cres$chr,
  ranges = IRanges(start = lit_cres$start, end = lit_cres$end),
  cre_id = lit_cres$cCRE1,
  genes = lit_cres$Gene,
  n_genes = lit_cres$n_genes
)

cat("\nCREs to visualize:\n")
for (i in 1:length(cre_ranges)) {
  cat(sprintf("  [%d] %s: %s:%d-%d (genes: %s)\n",
              i,
              cre_ranges$cre_id[i],
              as.character(seqnames(cre_ranges)[i]),
              start(cre_ranges)[i],
              end(cre_ranges)[i],
              cre_ranges$genes[i]))
}

# ============================================================================
# STEP 3: Get annotation
# ============================================================================
cat("\n")
cat(strrep("=", 80), "\n")
cat("STEP 3: Loading gene annotation...\n")
cat(strrep("-", 80), "\n")

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
annotation <- keepStandardChromosomes(annotation, pruning.mode = "coarse")
cat("✓ Loaded gene annotation\n")

# ============================================================================
# STEP 4: Create coverage plots for each CRE
# ============================================================================
cat("\n")
cat(strrep("=", 80), "\n")
cat("STEP 4: Creating coverage plots...\n")
cat(strrep("-", 80), "\n")

# Check if fragment files are available
fragment_path <- Fragments(seurat_obj[["ATAC"]])
if (length(fragment_path) == 0 || is.null(fragment_path)) {
  cat("⚠️  Warning: No fragment files found in Seurat object\n")
  cat("   Will attempt to use CoveragePlot with peaks/bins\n")
  use_fragments <- FALSE
} else {
  # Handle fragment path extraction safely
  tryCatch({
    if (is.list(fragment_path)) {
      fragment_paths <- sapply(fragment_path, function(x) {
        if (is.character(x)) return(basename(x))
        if (is.list(x) && !is.null(x$path)) return(basename(x$path))
        return("unknown")
      })
      cat(sprintf("✓ Fragment files available: %s\n",
                  paste(fragment_paths, collapse=", ")))
    } else if (is.character(fragment_path)) {
      cat(sprintf("✓ Fragment files available: %s\n", basename(fragment_path)))
    } else {
      cat("✓ Fragment files available (complex structure)\n")
    }
    use_fragments <- TRUE
  }, error = function(e) {
    cat("✓ Fragment files present (path extraction failed)\n")
    use_fragments <- TRUE
  })
}

# Process each CRE
for (i in 1:length(cre_ranges)) {
  cre <- cre_ranges[i]
  cre_id <- cre$cre_id
  cre_genes <- cre$genes
  cre_chr <- as.character(seqnames(cre))
  cre_start <- start(cre)
  cre_end <- end(cre)

  cat(sprintf("\n[%d/%d] Processing %s (genes: %s)\n",
              i, length(cre_ranges), cre_id, cre_genes))
  cat(sprintf("  Location: %s:%s-%s\n", cre_chr,
              format(cre_start, big.mark=","),
              format(cre_end, big.mark=",")))

  # Define extended region
  region_start <- max(1, cre_start - EXTEND_UPSTREAM)
  region_end <- cre_end + EXTEND_DOWNSTREAM

  # Create region string (Signac format: chr:start-end)
  region_string <- sprintf("%s:%d-%d", cre_chr, region_start, region_end)

  cat(sprintf("  Extended region: %s (%.1f kb)\n",
              region_string,
              (region_end - region_start) / 1000))

  # ============================================================================
  # Plot 1: Separate genotypes (4 groups)
  # ============================================================================
  cat("  Creating plot 1: Separate genotypes (4 groups)...\n")

  tryCatch({
    p1 <- CoveragePlot(
      object = seurat_obj,
      region = region_string,
      annotation = annotation,
      group.by = "plot_group",
      extend.upstream = 0,
      extend.downstream = 0,
      window = 200,
      ncol = 1
    )

    # Add CRE location highlighting
    p1 <- p1 +
      annotate("rect",
               xmin = cre_start, xmax = cre_end,
               ymin = -Inf, ymax = Inf,
               alpha = 0.15, fill = "red") +
      annotate("segment",
               x = cre_start, xend = cre_end,
               y = Inf, yend = Inf,
               color = "red", linewidth = 2) +
      annotate("text",
               x = (cre_start + cre_end) / 2,
               y = Inf,
               label = sprintf("Literature CRE\n%s", cre_id),
               vjust = 1.5, size = 3, color = "red", fontface = "bold") +
      ggtitle(sprintf("%s - ATAC Coverage (Separate Genotypes)\nGenes: %s | Location: %s",
                     cre_id, cre_genes, region_string)) +
      theme(
        plot.title = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 9, face = "bold")
      )

    # Save plot 1
    output_file_1 <- file.path(OUTPUT_DIR, sprintf("%s_separate_genotypes.pdf", cre_id))
    ggsave(output_file_1, p1, width = 14, height = 10)
    cat(sprintf("  ✓ Saved: %s\n", basename(output_file_1)))

  }, error = function(e) {
    cat(sprintf("  ✗ Error creating separate genotypes plot: %s\n", e$message))
  })

  # ============================================================================
  # Plot 2: Combined controls (3 groups)
  # ============================================================================
  cat("  Creating plot 2: Combined controls (3 groups)...\n")

  tryCatch({
    p2 <- CoveragePlot(
      object = seurat_obj,
      region = region_string,
      annotation = annotation,
      group.by = "plot_group_combined",
      extend.upstream = 0,
      extend.downstream = 0,
      window = 200,
      ncol = 1
    )

    # Add CRE location highlighting
    p2 <- p2 +
      annotate("rect",
               xmin = cre_start, xmax = cre_end,
               ymin = -Inf, ymax = Inf,
               alpha = 0.15, fill = "red") +
      annotate("segment",
               x = cre_start, xend = cre_end,
               y = Inf, yend = Inf,
               color = "red", linewidth = 2) +
      annotate("text",
               x = (cre_start + cre_end) / 2,
               y = Inf,
               label = sprintf("Literature CRE\n%s", cre_id),
               vjust = 1.5, size = 3, color = "red", fontface = "bold") +
      ggtitle(sprintf("%s - ATAC Coverage (Combined Controls)\nGenes: %s | Location: %s",
                     cre_id, cre_genes, region_string)) +
      theme(
        plot.title = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 9, face = "bold")
      )

    # Save plot 2
    output_file_2 <- file.path(OUTPUT_DIR, sprintf("%s_combined_controls.pdf", cre_id))
    ggsave(output_file_2, p2, width = 14, height = 8)
    cat(sprintf("  ✓ Saved: %s\n", basename(output_file_2)))

  }, error = function(e) {
    cat(sprintf("  ✗ Error creating combined controls plot: %s\n", e$message))
  })

  # ============================================================================
  # Plot 3: Zoomed on CRE only (±10kb)
  # ============================================================================
  cat("  Creating plot 3: Zoomed on CRE (±10kb)...\n")

  tryCatch({
    zoom_start <- max(1, cre_start - 10000)
    zoom_end <- cre_end + 10000
    zoom_string <- sprintf("%s:%d-%d", cre_chr, zoom_start, zoom_end)

    p3 <- CoveragePlot(
      object = seurat_obj,
      region = zoom_string,
      annotation = annotation,
      group.by = "plot_group_combined",
      extend.upstream = 0,
      extend.downstream = 0,
      window = 50,
      ncol = 1
    )

    p3 <- p3 +
      annotate("rect",
               xmin = cre_start, xmax = cre_end,
               ymin = -Inf, ymax = Inf,
               alpha = 0.15, fill = "red") +
      annotate("segment",
               x = cre_start, xend = cre_end,
               y = Inf, yend = Inf,
               color = "red", linewidth = 3) +
      annotate("text",
               x = (cre_start + cre_end) / 2,
               y = Inf,
               label = sprintf("CRE: %s", cre_id),
               vjust = 1.5, size = 4, color = "red", fontface = "bold") +
      ggtitle(sprintf("%s - Zoomed View (±10kb)\nGenes: %s",
                     cre_id, cre_genes)) +
      theme(
        plot.title = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 9, face = "bold")
      )

    output_file_3 <- file.path(OUTPUT_DIR, sprintf("%s_zoomed.pdf", cre_id))
    ggsave(output_file_3, p3, width = 14, height = 8)
    cat(sprintf("  ✓ Saved: %s\n", basename(output_file_3)))

  }, error = function(e) {
    cat(sprintf("  ✗ Error creating zoomed plot: %s\n", e$message))
  })
}

# ============================================================================
# STEP 5: Create multi-panel summary
# ============================================================================
cat("\n")
cat(strrep("=", 80), "\n")
cat("STEP 5: Creating multi-panel summary...\n")
cat(strrep("-", 80), "\n")

# Create a summary showing all CREs in one figure (zoomed view)
cat("Creating compact multi-panel summary...\n")

plot_list <- list()
success_count <- 0

for (i in 1:min(length(cre_ranges), 12)) {  # Limit to 12 for readability
  cre <- cre_ranges[i]
  cre_id <- cre$cre_id
  cre_chr <- as.character(seqnames(cre))
  cre_start <- start(cre)
  cre_end <- end(cre)
  cre_genes <- cre$genes

  tryCatch({
    zoom_start <- max(1, cre_start - 25000)
    zoom_end <- cre_end + 25000
    zoom_string <- sprintf("%s:%d-%d", cre_chr, zoom_start, zoom_end)

    p <- CoveragePlot(
      object = seurat_obj,
      region = zoom_string,
      annotation = annotation,
      group.by = "plot_group_combined",
      extend.upstream = 0,
      extend.downstream = 0,
      window = 200,
      ncol = 1
    )

    p <- p +
      annotate("rect",
               xmin = cre_start, xmax = cre_end,
               ymin = -Inf, ymax = Inf,
               alpha = 0.2, fill = "red") +
      ggtitle(sprintf("%s (%s)", cre_id, cre_genes)) +
      theme(
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        legend.position = "none",
        strip.text = element_text(size = 6)
      )

    plot_list[[length(plot_list) + 1]] <- p
    success_count <- success_count + 1

  }, error = function(e) {
    cat(sprintf("  Skipping %s for summary (error: %s)\n", cre_id, e$message))
  })
}

if (length(plot_list) > 0) {
  # Arrange in grid
  n_plots <- length(plot_list)
  n_cols <- min(3, n_plots)
  n_rows <- ceiling(n_plots / n_cols)

  summary_plot <- wrap_plots(plot_list, ncol = n_cols)

  output_file_summary <- file.path(OUTPUT_DIR, "all_CREs_summary.pdf")
  ggsave(output_file_summary, summary_plot, width = 18, height = 6 * n_rows)
  cat(sprintf("✓ Saved multi-panel summary: %s (%d CREs)\n",
              basename(output_file_summary), success_count))
} else {
  cat("⚠️  No plots created for summary\n")
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================
cat("\n")
cat(strrep("=", 80), "\n")
cat("VISUALIZATION COMPLETE!\n")
cat(strrep("=", 80), "\n")

cat(sprintf("\nProcessed %d literature CREs\n", length(cre_ranges)))
cat(sprintf("\nOutput directory: %s\n", OUTPUT_DIR))

cat("\nFiles created per CRE:\n")
cat("  - *_separate_genotypes.pdf (4 groups: Nestin-Ctrl, Nestin-Mut, Emx1-Ctrl, Emx1-Mut)\n")
cat("  - *_combined_controls.pdf (3 groups: Combined-Ctrl, Nestin-Mut, Emx1-Mut)\n")
cat("  - *_zoomed.pdf (±10kb view around CRE)\n")
cat("\nSummary file:\n")
cat("  - all_CREs_summary.pdf (all CREs in one figure)\n")

cat("\nInterpretation guide:\n")
cat("  - Red shaded box = Literature CRE location\n")
cat("  - Coverage tracks = ATAC-seq signal (fragment pileup)\n")
cat("  - Higher coverage = More accessible chromatin\n")
cat("  - Compare Ctrl vs Mut to see accessibility changes\n")
cat("  - Even without peak calls, you can see if signal is present\n")

cat("\n")
cat(strrep("=", 80), "\n")
