#!/usr/bin/env Rscript

################################################################################
# Visualization of Peak-Gene Linkages for GABA Neurons (Separate Genotypes)
# LEGACY VERSION
#
# ⚠ COMPATIBILITY NOTE: This script expects gene-specific linkage files
# organized differently from the new array job outputs:
#   - Expected: peak_gene_linkage_analysis_separate/GABA/by_genotype/{gene}_{genotype}_peak_gene_links.csv
#   - Array jobs produce: celltype_results/peak_gene_links/GABA_{Genotype}_peak_gene_links.csv
#
# RECOMMENDATION: Use analyze_peak_gene_linkage_gaba_separate_UPDATED.R instead,
# which is compatible with array job outputs.
#
# This script creates visualizations for peak-gene linkages with separate
# analysis per genotype (no combined controls).
#
# Generates:
# 1. Coverage plots showing ATAC signal for all 4 groups (Nestin-Ctrl, Nestin-Mut, Emx1-Ctrl, Emx1-Mut)
# 2. Per-cell fragment abundance plots by genotype and condition
# 3. Correlation scatter plots for peak-gene pairs (by genotype)
#
# Input files:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/integrated_seurat_processed.rds
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/genes_inter.txt
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/peak_gene_linkage_analysis_separate/GABA/by_genotype/{gene}_{genotype}_peak_gene_links.csv
#
# Output:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/peak_gene_linkage_analysis_separate/GABA/visualizations/coverage_{gene}_gaba_separate.pdf
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/peak_gene_linkage_analysis_separate/GABA/visualizations/fragment_abundance_gaba_separate.pdf
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/peak_gene_linkage_analysis_separate/GABA/visualizations/correlation_scatter_{gene}_{genotype}_gaba_separate.pdf
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
  library(gridExtra)
})

# Configuration
RESULTS_DIR <- "signac_results_L1"
SEURAT_FILE <- file.path(RESULTS_DIR, "integrated_seurat_processed.rds")
GENE_LIST_FILE <- "../genes_inter.txt"
LINKAGE_DIR <- file.path(RESULTS_DIR, "peak_gene_linkage_analysis_separate", "GABA")
OUTPUT_DIR <- file.path(LINKAGE_DIR, "visualizations")

# Cell type and condition columns
CELLTYPE_COL <- "cell_type_L1"
CELL_TYPE <- "GABA"
CONDITION_COL <- "condition"
GENOTYPE_COL <- "genotype"

# Plotting parameters
EXTEND_UPSTREAM <- 100000    # Extend 100kb upstream
EXTEND_DOWNSTREAM <- 100000  # Extend 100kb downstream

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("================================================================================\n")
cat("VISUALIZATION: PEAK-GENE LINKAGES - GABA NEURONS (SEPARATE GENOTYPES)\n")
cat("================================================================================\n\n")

################################################################################
# Load Data
################################################################################

cat("Loading integrated Seurat object...\n")
if (!file.exists(SEURAT_FILE)) {
  stop(sprintf("Error: Seurat file not found: %s\n", SEURAT_FILE))
}
seurat_obj <- readRDS(SEURAT_FILE)
cat(sprintf("  ✓ Loaded Seurat object with %d cells\n", ncol(seurat_obj)))

cat("\nLoading gene list...\n")
genes_of_interest <- read.table(GENE_LIST_FILE, header = FALSE, stringsAsFactors = FALSE)$V1
genes_of_interest <- unique(genes_of_interest)
cat(sprintf("  ✓ Loaded %d genes\n", length(genes_of_interest)))

cat("\nSubsetting to GABA neurons...\n")
gaba_cells <- seurat_obj@meta.data %>%
  filter(.data[[CELLTYPE_COL]] == CELL_TYPE) %>%
  rownames()
seurat_gaba <- subset(seurat_obj, cells = gaba_cells)
cat(sprintf("  ✓ Selected %d GABA neurons\n", ncol(seurat_gaba)))

# Add genotype-condition grouping (4 groups, no combined controls)
seurat_gaba@meta.data$genotype_condition <- paste(
  seurat_gaba@meta.data[[GENOTYPE_COL]],
  seurat_gaba@meta.data[[CONDITION_COL]],
  sep = "-"
)

################################################################################
# Per-Cell Fragment Abundance Plots
################################################################################

cat("\n================================================================================\n")
cat("CREATING FRAGMENT ABUNDANCE PLOTS\n")
cat("================================================================================\n\n")

# Get fragment counts per cell
if ("nCount_ATAC" %in% colnames(seurat_gaba@meta.data)) {
  cat("Creating per-cell fragment abundance plots...\n")

  # Prepare data for plotting
  plot_data <- seurat_gaba@meta.data %>%
    select(nCount_ATAC, genotype_condition, all_of(GENOTYPE_COL), all_of(CONDITION_COL)) %>%
    mutate(
      log10_fragments = log10(nCount_ATAC + 1),
      genotype = .data[[GENOTYPE_COL]],
      condition = .data[[CONDITION_COL]]
    )

  # Violin plot by genotype-condition (all 4 groups)
  p1 <- ggplot(plot_data, aes(x = genotype_condition, y = log10_fragments, fill = genotype_condition)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c(
      "Nestin-Ctrl" = "#4CAF50",
      "Nestin-Mut" = "#F44336",
      "Emx1-Ctrl" = "#2196F3",
      "Emx1-Mut" = "#FF9800"
    )) +
    labs(
      title = "ATAC Fragment Counts per Cell - GABA Neurons",
      subtitle = "Separate Groups (No Combined Controls)",
      x = "Group",
      y = "log10(Fragment Count + 1)",
      fill = "Group"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  # Violin plot by genotype and condition (faceted)
  p2 <- ggplot(plot_data, aes(x = condition, y = log10_fragments, fill = condition)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
    facet_wrap(~ genotype, nrow = 1) +
    scale_fill_manual(values = c("Ctrl" = "#4CAF50", "Mut" = "#F44336")) +
    labs(
      title = "ATAC Fragment Counts by Genotype",
      x = "Condition",
      y = "log10(Fragment Count + 1)",
      fill = "Condition"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "lightgray")
    )

  # Save combined plot
  pdf(file.path(OUTPUT_DIR, "fragment_abundance_gaba_separate.pdf"), width = 14, height = 6)
  print(p1 + p2 + plot_layout(widths = c(1, 1)))
  dev.off()

  cat("  ✓ Saved fragment abundance plots\n")

  # Print summary statistics
  cat("\n  Fragment count summary by group:\n")
  summary_by_group <- plot_data %>%
    group_by(genotype_condition) %>%
    summarise(
      n_cells = n(),
      median_fragments = median(nCount_ATAC),
      mean_fragments = mean(nCount_ATAC),
      sd_fragments = sd(nCount_ATAC),
      .groups = "drop"
    )
  print(summary_by_group)
} else {
  cat("  ⚠ Warning: nCount_ATAC not found in metadata\n")
}

################################################################################
# Coverage Plots for Genes with Linkages
################################################################################

cat("\n================================================================================\n")
cat("CREATING COVERAGE PLOTS\n")
cat("================================================================================\n\n")

for (gene in genes_of_interest) {
  # Check if any genotype has linkage results for this gene
  has_links <- FALSE
  for (geno in c("Nestin", "Emx1")) {
    linkage_file <- file.path(LINKAGE_DIR, "by_genotype", sprintf("%s_%s_peak_gene_links.csv", gene, geno))
    if (file.exists(linkage_file)) {
      links <- read.csv(linkage_file)
      if (nrow(links) > 0) {
        has_links <- TRUE
        break
      }
    }
  }

  if (!has_links) {
    cat(sprintf("  Skipping %s (no significant links in any genotype)\n", gene))
    next
  }

  cat(sprintf("\n--- Creating coverage plot for %s ---\n", gene))

  # Create coverage plot by genotype-condition (all 4 groups)
  tryCatch({
    # Set default assay
    DefaultAssay(seurat_gaba) <- "ATAC"

    # Create coverage plot
    p_cov <- CoveragePlot(
      object = seurat_gaba,
      region = gene,
      features = gene,
      expression.assay = "RNA",
      assay = "ATAC",
      group.by = "genotype_condition",
      extend.upstream = EXTEND_UPSTREAM,
      extend.downstream = EXTEND_DOWNSTREAM,
      annotation = TRUE,
      peaks = TRUE,
      window = 100
    )

    # Save plot
    pdf_file <- file.path(OUTPUT_DIR, sprintf("coverage_%s_gaba_separate.pdf", gene))
    ggsave(pdf_file, plot = p_cov, width = 16, height = 12, device = "pdf")

    cat(sprintf("  ✓ Saved coverage plot: %s\n", basename(pdf_file)))
  }, error = function(e) {
    cat(sprintf("  ✗ Error creating coverage plot for %s: %s\n", gene, e$message))
  })
}

################################################################################
# Correlation Scatter Plots (Separate by Genotype)
################################################################################

cat("\n================================================================================\n")
cat("CREATING CORRELATION SCATTER PLOTS\n")
cat("================================================================================\n\n")

for (gene in genes_of_interest) {
  # Create separate scatter plots for each genotype
  for (geno in c("Nestin", "Emx1")) {
    linkage_file <- file.path(LINKAGE_DIR, "by_genotype", sprintf("%s_%s_peak_gene_links.csv", gene, geno))

    if (!file.exists(linkage_file)) {
      next
    }

    links <- read.csv(linkage_file)

    if (nrow(links) == 0) {
      next
    }

    cat(sprintf("\n--- Creating scatter plots for %s (%s) ---\n", gene, geno))

    # Subset to this genotype
    cells_geno <- colnames(seurat_gaba)[seurat_gaba@meta.data[[GENOTYPE_COL]] == geno]
    seurat_geno <- subset(seurat_gaba, cells = cells_geno)

    # Create scatter plot for top 6 peaks (by correlation strength)
    top_links <- links %>%
      arrange(desc(abs(correlation))) %>%
      head(6)

    plots_list <- list()

    for (i in 1:min(6, nrow(top_links))) {
      peak <- top_links$peak[i]
      cor_val <- top_links$correlation[i]
      pval <- top_links$padj[i]

      # Get peak accessibility and gene expression
      peak_acc <- GetAssayData(seurat_geno, assay = "ATAC", layer = "data")[peak, ]
      gene_expr <- GetAssayData(seurat_geno, assay = "RNA", layer = "data")[gene, ]

      # Create data frame
      scatter_data <- data.frame(
        peak_accessibility = as.numeric(peak_acc),
        gene_expression = as.numeric(gene_expr),
        condition = seurat_geno@meta.data[[CONDITION_COL]]
      )

      # Create scatter plot
      p <- ggplot(scatter_data, aes(x = peak_accessibility, y = gene_expression, color = condition)) +
        geom_point(alpha = 0.5, size = 1.5) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 0.5) +
        scale_color_manual(values = c("Ctrl" = "#4CAF50", "Mut" = "#F44336")) +
        labs(
          title = sprintf("Peak %d (%s)", i, geno),
          subtitle = sprintf("r = %.3f, FDR = %.2e", cor_val, pval),
          x = "Peak Accessibility",
          y = sprintf("%s Expression", gene),
          color = "Condition"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8),
          legend.position = if (i == 1) "right" else "none"
        )

      plots_list[[i]] <- p
    }

    # Arrange plots
    if (length(plots_list) > 0) {
      pdf_file <- file.path(OUTPUT_DIR, sprintf("correlation_scatter_%s_%s_gaba_separate.pdf", gene, geno))
      pdf(pdf_file, width = 14, height = 10)
      grid.arrange(grobs = plots_list, ncol = 3)
      dev.off()

      cat(sprintf("  ✓ Saved scatter plots: %s\n", basename(pdf_file)))
    }
  }
}

cat("\n================================================================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("================================================================================\n\n")
cat(sprintf("All plots saved to: %s\n", OUTPUT_DIR))
