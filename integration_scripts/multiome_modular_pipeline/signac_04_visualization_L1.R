#!/usr/bin/env Rscript

################################################################################
# Signac Multiome Analysis - Step 4: Advanced Visualization (L1 VERSION)
################################################################################
#
# This script creates publication-quality visualizations for L1 cell types:
#   1. Volcano plots for DEG and DA (9 broad categories)
#   2. Coverage plots for key genes showing ATAC signal
#   3. Peak-gene linkage visualizations
#   4. Integration plots with existing DEG results
#
# L1 VERSION: Visualizations for 9 broad cell categories
# Outputs to: signac_results_L1/plots/
#
# Inputs:
#   - Integrated Seurat object from Step 2:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/integrated_seurat_processed.rds
#   - Summary statistics from Step 3:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/summary_statistics.csv
#   - Differential expression results (DEGs) from Step 3:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/DEG/{celltype}_{genotype}_DEGs.csv
#   - Differential accessibility results (DA peaks) from Step 3:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/DA/{celltype}_{genotype}_DA_peaks.csv
#   - Peak-gene linkage results from Step 3:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/peak_gene_links/{celltype}_{genotype}_peak_gene_links.csv
#
# Outputs:
#   - DEG volcano plots:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/plots/03_DEG_volcano_plots.pdf
#   - DA volcano plots:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/plots/04_DA_volcano_plots.pdf
#   - Peak-gene link summary plots:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/plots/06_peak_gene_link_summary.pdf
#   - Summary barplot:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/plots/07_summary_barplot.pdf
#   - Peak-gene link summary table:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/plots/peak_gene_link_summary.csv
#
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

cat("================================================================================\n")
cat("  SIGNAC MULTIOME ANALYSIS - STEP 4: ADVANCED VISUALIZATION\n")
cat("================================================================================\n\n")

################################################################################
# Configuration
################################################################################

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
INPUT_DIR <- file.path(BASE_DIR, "signac_results_L1")
RESULTS_DIR <- file.path(INPUT_DIR, "celltype_results")
OUTPUT_DIR <- file.path(INPUT_DIR, "plots")

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Cell types for detailed visualization (top by cell count)
TOP_N_CELLTYPES <- 5

CELLTYPE_COL <- "cell_type_L1"
CONDITION_COL <- "atac_sample"

################################################################################
# Load Data
################################################################################

cat("\n=== Loading data ===\n")

# Load Seurat object
seurat_obj <- readRDS(file.path(INPUT_DIR, "integrated_seurat_processed.rds"))
cat(sprintf("  ✓ Loaded Seurat object: %d cells\n", ncol(seurat_obj)))

# Load summary statistics
summary_stats <- read.csv(file.path(RESULTS_DIR, "summary_statistics.csv"))
cat(sprintf("  ✓ Loaded summary statistics: %d cell types\n", nrow(summary_stats)))

################################################################################
# Function: Create Volcano Plot
################################################################################

create_volcano_plot <- function(results_df, title, logfc_col = "avg_log2FC",
                                 pval_col = "p_val_adj", logfc_thresh = 0.25,
                                 pval_thresh = 0.05, max_labels = 10) {

  # Add significance column
  results_df$significant <- results_df[[pval_col]] < pval_thresh &
    abs(results_df[[logfc_col]]) > logfc_thresh

  results_df$direction <- "NS"
  results_df$direction[results_df$significant & results_df[[logfc_col]] > 0] <- "Up"
  results_df$direction[results_df$significant & results_df[[logfc_col]] < 0] <- "Down"

  # -log10 transform p-values
  results_df$neg_log10_pval <- -log10(results_df[[pval_col]] + 1e-300)

  # Label top genes
  if ("gene" %in% colnames(results_df)) {
    results_df <- results_df %>%
      arrange(!!sym(pval_col)) %>%
      mutate(label = ifelse(row_number() <= max_labels & significant, gene, ""))
  } else {
    results_df$label <- ""
  }

  # Create plot
  p <- ggplot(results_df, aes(x = !!sym(logfc_col), y = neg_log10_pval,
                               color = direction, label = label)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "grey")) +
    geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed", color = "black") +
    labs(
      title = title,
      x = "Log2 Fold Change (Mut vs Ctrl)",
      y = "-Log10(Adjusted P-value)",
      color = "Direction"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  # Add labels if present
  if (any(results_df$label != "")) {
    p <- p + ggrepel::geom_text_repel(max.overlaps = 20, size = 3)
  }

  return(p)
}


################################################################################
# Volcano Plots for Top Cell Types
################################################################################

cat("\n================================================================================\n")
cat("STEP 1: Creating Volcano Plots\n")
cat("================================================================================\n")

# Get top cell types by total cells (considering genotype-stratified results)
summary_stats_filtered <- summary_stats %>%
  filter(!is.na(n_degs)) %>%
  mutate(total_cells = n_cells_ctrl + n_cells_mut) %>%
  arrange(desc(total_cells))

# For genotype-stratified results, we have celltype + genotype
# Let's get top cell-type-genotype combinations
top_analyses <- head(summary_stats_filtered, TOP_N_CELLTYPES * 2)  # Get more to show both genotypes

# Check if genotype column exists (new format)
has_genotype <- "genotype" %in% colnames(summary_stats_filtered)

if (has_genotype) {
  cat(sprintf("  Creating volcano plots for top %d cell type × genotype combinations:\n",
              nrow(top_analyses)))
  for (i in 1:nrow(top_analyses)) {
    cat(sprintf("    %s (%s): %d cells\n",
                top_analyses$celltype[i],
                top_analyses$genotype[i],
                top_analyses$total_cells[i]))
  }
} else {
  # Old format without genotype
  top_celltypes <- head(summary_stats_filtered$celltype, TOP_N_CELLTYPES)
  cat(sprintf("  Creating volcano plots for top %d cell types:\n", TOP_N_CELLTYPES))
  cat(sprintf("    %s\n", paste(top_celltypes, collapse = ", ")))
}

# DEG volcano plots
deg_volcano_plots <- list()

if (has_genotype) {
  # New format with genotype stratification
  for (i in 1:nrow(top_analyses)) {
    ct <- top_analyses$celltype[i]
    geno <- top_analyses$genotype[i]
    analysis_id <- sprintf("%s_%s", ct, geno)

    deg_file <- file.path(RESULTS_DIR, "DEG", sprintf("%s_DEGs.csv", analysis_id))

    if (file.exists(deg_file)) {
      deg_res <- read.csv(deg_file)

      if (nrow(deg_res) > 0) {
        p <- create_volcano_plot(
          deg_res,
          title = sprintf("%s (%s) - DEG", ct, geno),
          max_labels = 15
        )
        deg_volcano_plots[[analysis_id]] <- p
      }
    } else {
      cat(sprintf("  ⚠ File not found: %s\n", deg_file))
    }
  }
} else {
  # Old format without genotype
  for (ct in top_celltypes) {
    deg_file <- file.path(RESULTS_DIR, "DEG", sprintf("%s_DEGs.csv", ct))

    if (file.exists(deg_file)) {
      deg_res <- read.csv(deg_file)

      if (nrow(deg_res) > 0) {
        p <- create_volcano_plot(
          deg_res,
          title = sprintf("%s - Differential Gene Expression", ct),
          max_labels = 15
        )
        deg_volcano_plots[[ct]] <- p
      }
    }
  }
}

# Save DEG volcano plots
if (length(deg_volcano_plots) > 0) {
  pdf(file.path(OUTPUT_DIR, "03_DEG_volcano_plots.pdf"), width = 20, height = 12)
  if (length(deg_volcano_plots) >= 4) {
    print(wrap_plots(deg_volcano_plots[1:4], ncol = 2))
    if (length(deg_volcano_plots) > 4) {
      print(wrap_plots(deg_volcano_plots[5:length(deg_volcano_plots)], ncol = 2))
    }
  } else {
    print(wrap_plots(deg_volcano_plots, ncol = 2))
  }
  dev.off()

  cat(sprintf("  ✓ DEG volcano plots saved: %s\n",
              file.path(OUTPUT_DIR, "03_DEG_volcano_plots.pdf")))
}

# DA volcano plots
da_volcano_plots <- list()

if (has_genotype) {
  # New format with genotype stratification
  for (i in 1:nrow(top_analyses)) {
    ct <- top_analyses$celltype[i]
    geno <- top_analyses$genotype[i]
    analysis_id <- sprintf("%s_%s", ct, geno)

    da_file <- file.path(RESULTS_DIR, "DA", sprintf("%s_DA_peaks.csv", analysis_id))

    if (file.exists(da_file)) {
      da_res <- read.csv(da_file)

      if (nrow(da_res) > 0) {
        p <- create_volcano_plot(
          da_res,
          title = sprintf("%s (%s) - DA", ct, geno),
          max_labels = 0  # Don't label peaks
        )
        da_volcano_plots[[analysis_id]] <- p
      }
    }
  }
} else {
  # Old format without genotype
  for (ct in top_celltypes) {
    da_file <- file.path(RESULTS_DIR, "DA", sprintf("%s_DA_peaks.csv", ct))

    if (file.exists(da_file)) {
      da_res <- read.csv(da_file)

      if (nrow(da_res) > 0) {
        p <- create_volcano_plot(
          da_res,
          title = sprintf("%s - Differential Accessibility", ct),
          max_labels = 0  # Don't label peaks
        )
        da_volcano_plots[[ct]] <- p
      }
    }
  }
}

# Save DA volcano plots
if (length(da_volcano_plots) > 0) {
  pdf(file.path(OUTPUT_DIR, "04_DA_volcano_plots.pdf"), width = 20, height = 12)
  if (length(da_volcano_plots) >= 4) {
    print(wrap_plots(da_volcano_plots[1:4], ncol = 2))
    if (length(da_volcano_plots) > 4) {
      print(wrap_plots(da_volcano_plots[5:length(da_volcano_plots)], ncol = 2))
    }
  } else {
    print(wrap_plots(da_volcano_plots, ncol = 2))
  }
  dev.off()

  cat(sprintf("  ✓ DA volcano plots saved: %s\n",
              file.path(OUTPUT_DIR, "04_DA_volcano_plots.pdf")))
}

################################################################################
# Note: Coverage Plots Moved to Separate Script
################################################################################
# Coverage plots require fragment files with translated barcodes.
# See signac_05_coverage_plots.R for proper coverage plot generation
# using barcode-translated fragment files.

################################################################################
# Peak-Gene Link Summary Plots
################################################################################

cat("\n================================================================================\n")
cat("STEP 3: Summarizing Peak-Gene Links\n")
cat("================================================================================\n")

# Compile all peak-gene links
all_links <- list()

if (has_genotype) {
  # New format with genotype stratification
  for (i in 1:nrow(top_analyses)) {
    ct <- top_analyses$celltype[i]
    geno <- top_analyses$genotype[i]
    analysis_id <- sprintf("%s_%s", ct, geno)

    links_file <- file.path(RESULTS_DIR, "peak_gene_links", sprintf("%s_peak_gene_links.csv", analysis_id))

    if (file.exists(links_file)) {
      links <- read.csv(links_file)
      if (nrow(links) > 0) {
        links$celltype <- ct
        links$genotype <- geno
        links$analysis_id <- analysis_id
        all_links[[analysis_id]] <- links
      }
    }
  }
} else {
  # Old format without genotype
  for (ct in top_celltypes) {
    links_file <- file.path(RESULTS_DIR, "peak_gene_links", sprintf("%s_peak_gene_links.csv", ct))

    if (file.exists(links_file)) {
      links <- read.csv(links_file)
      if (nrow(links) > 0) {
        all_links[[ct]] <- links
      }
    }
  }
}

if (length(all_links) > 0) {
  # Combine all links
  if (has_genotype) {
    # For genotype-stratified data, we already added celltype/genotype/analysis_id columns
    combined_links <- bind_rows(all_links)
  } else {
    # For old format, bind with .id to create celltype column
    combined_links <- bind_rows(all_links, .id = "celltype")
  }

  # Extract peak start coordinates from peak column (format: "chr-start-end")
  combined_links <- combined_links %>%
    mutate(
      peak_start = as.numeric(sub(".*-(\\d+)-\\d+$", "\\1", peak)),
      gene_start = start
    )

  # Summary statistics
  if (has_genotype) {
    link_summary <- combined_links %>%
      group_by(analysis_id, celltype, genotype) %>%
      summarise(
        n_links = n(),
        n_unique_genes = n_distinct(gene),
        n_unique_peaks = n_distinct(peak),
        mean_correlation = mean(score, na.rm = TRUE),
        median_distance = median(abs(gene_start - peak_start), na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    link_summary <- combined_links %>%
      group_by(celltype) %>%
      summarise(
        n_links = n(),
        n_unique_genes = n_distinct(gene),
        n_unique_peaks = n_distinct(peak),
        mean_correlation = mean(score, na.rm = TRUE),
        median_distance = median(abs(gene_start - peak_start), na.rm = TRUE)
      )
  }

  # Plot: Number of links per cell type (or cell type x genotype)
  if (has_genotype) {
    p1 <- ggplot(link_summary, aes(x = reorder(analysis_id, n_links), y = n_links, fill = genotype)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(
        title = "Peak-Gene Links per Cell Type × Genotype",
        x = "Cell Type × Genotype",
        y = "Number of Links"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  } else {
    p1 <- ggplot(link_summary, aes(x = reorder(celltype, n_links), y = n_links)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(
        title = "Peak-Gene Links per Cell Type",
        x = "Cell Type",
        y = "Number of Links"
      ) +
      theme_bw()
  }

  # Plot: Distribution of correlations
  if (has_genotype) {
    p2 <- ggplot(combined_links, aes(x = score, fill = genotype)) +
      geom_density(alpha = 0.5) +
      labs(
        title = "Distribution of Peak-Gene Correlations",
        x = "Correlation Score",
        y = "Density"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  } else {
    p2 <- ggplot(combined_links, aes(x = score, fill = celltype)) +
      geom_density(alpha = 0.5) +
      labs(
        title = "Distribution of Peak-Gene Correlations",
        x = "Correlation Score",
        y = "Density"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  }

  # Plot: Distribution of distances
  combined_links$distance <- abs(combined_links$gene_start - combined_links$peak_start)

  if (has_genotype) {
    p3 <- ggplot(combined_links, aes(x = distance / 1000, fill = genotype)) +
      geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
      labs(
        title = "Distribution of Peak-Gene Distances",
        x = "Distance (kb)",
        y = "Count"
      ) +
      scale_x_continuous(limits = c(0, 200)) +
      theme_bw() +
      theme(legend.position = "bottom")
  } else {
    p3 <- ggplot(combined_links, aes(x = distance / 1000, fill = celltype)) +
      geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
      labs(
        title = "Distribution of Peak-Gene Distances",
        x = "Distance (kb)",
        y = "Count"
      ) +
      scale_x_continuous(limits = c(0, 200)) +
      theme_bw() +
      theme(legend.position = "bottom")
  }

  # Save plots
  pdf(file.path(OUTPUT_DIR, "06_peak_gene_link_summary.pdf"), width = 14, height = 10)
  print(p1)
  print(p2)
  print(p3)
  print(p1 / (p2 | p3))
  dev.off()

  # Save summary table
  write.csv(link_summary, file.path(OUTPUT_DIR, "peak_gene_link_summary.csv"),
            row.names = FALSE)

  cat(sprintf("  ✓ Peak-gene link plots saved: %s\n",
              file.path(OUTPUT_DIR, "06_peak_gene_link_summary.pdf")))
}

################################################################################
# Summary Barplot
################################################################################

cat("\n================================================================================\n")
cat("STEP 4: Creating Summary Barplot\n")
cat("================================================================================\n")

# Reshape summary stats for plotting
if (has_genotype) {
  # For genotype-stratified data, create analysis_id for x-axis
  summary_long <- summary_stats %>%
    filter(!is.na(n_degs)) %>%
    mutate(analysis_id = sprintf("%s_%s", celltype, genotype)) %>%
    select(analysis_id, celltype, genotype, n_degs_up, n_degs_down, n_da_more, n_da_less) %>%
    pivot_longer(cols = c(n_degs_up, n_degs_down, n_da_more, n_da_less),
                 names_to = "metric", values_to = "count") %>%
    mutate(
      category = case_when(
        grepl("deg", metric) ~ "Gene Expression",
        grepl("da", metric) ~ "Accessibility"
      ),
      direction = case_when(
        grepl("up|more", metric) ~ "Up/More",
        grepl("down|less", metric) ~ "Down/Less"
      )
    )

  # Order by total changes
  analysis_order <- summary_stats %>%
    filter(!is.na(n_degs)) %>%
    mutate(
      analysis_id = sprintf("%s_%s", celltype, genotype),
      total_changes = n_degs + n_da_peaks
    ) %>%
    arrange(desc(total_changes)) %>%
    pull(analysis_id)

  summary_long$analysis_id <- factor(summary_long$analysis_id, levels = analysis_order)

  # Create plot with genotype faceting
  p <- ggplot(summary_long, aes(x = analysis_id, y = count, fill = interaction(category, direction))) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(
      values = c(
        "Gene Expression.Up/More" = "#E41A1C",
        "Gene Expression.Down/Less" = "#377EB8",
        "Accessibility.Up/More" = "#FF7F00",
        "Accessibility.Down/Less" = "#4DAF4A"
      ),
      labels = c(
        "Gene Expression.Up/More" = "DEG Up",
        "Gene Expression.Down/Less" = "DEG Down",
        "Accessibility.Up/More" = "DA More Accessible",
        "Accessibility.Down/Less" = "DA Less Accessible"
      )
    ) +
    labs(
      title = "Summary of Differential Changes by Cell Type × Genotype",
      x = "Cell Type × Genotype",
      y = "Number of Changes",
      fill = "Category"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
} else {
  # Old format without genotype
  summary_long <- summary_stats %>%
    filter(!is.na(n_degs)) %>%
    select(celltype, n_degs_up, n_degs_down, n_da_more, n_da_less) %>%
    pivot_longer(cols = -celltype, names_to = "metric", values_to = "count") %>%
    mutate(
      category = case_when(
        grepl("deg", metric) ~ "Gene Expression",
        grepl("da", metric) ~ "Accessibility"
      ),
      direction = case_when(
        grepl("up|more", metric) ~ "Up/More",
        grepl("down|less", metric) ~ "Down/Less"
      )
    )

  # Order by total changes
  celltype_order <- summary_stats %>%
    filter(!is.na(n_degs)) %>%
    mutate(total_changes = n_degs + n_da_peaks) %>%
    arrange(desc(total_changes)) %>%
    pull(celltype)

  summary_long$celltype <- factor(summary_long$celltype, levels = celltype_order)

  # Create plot
  p <- ggplot(summary_long, aes(x = celltype, y = count, fill = interaction(category, direction))) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(
      values = c(
        "Gene Expression.Up/More" = "#E41A1C",
        "Gene Expression.Down/Less" = "#377EB8",
        "Accessibility.Up/More" = "#FF7F00",
        "Accessibility.Down/Less" = "#4DAF4A"
      ),
      labels = c(
        "Gene Expression.Up/More" = "DEG Up",
        "Gene Expression.Down/Less" = "DEG Down",
        "Accessibility.Up/More" = "DA More Accessible",
        "Accessibility.Down/Less" = "DA Less Accessible"
      )
    ) +
    labs(
      title = "Summary of Differential Changes by Cell Type",
      x = "Cell Type",
      y = "Number of Changes",
      fill = "Category"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

pdf(file.path(OUTPUT_DIR, "07_summary_barplot.pdf"), width = 14, height = 8)
print(p)
dev.off()

cat(sprintf("  ✓ Summary barplot saved: %s\n",
            file.path(OUTPUT_DIR, "07_summary_barplot.pdf")))

################################################################################
# Final Summary
################################################################################

cat("\n================================================================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("================================================================================\n\n")

cat("Output files created:\n")
cat(sprintf("  - %s/03_DEG_volcano_plots.pdf\n", OUTPUT_DIR))
cat(sprintf("  - %s/04_DA_volcano_plots.pdf\n", OUTPUT_DIR))
cat(sprintf("  - %s/06_peak_gene_link_summary.pdf\n", OUTPUT_DIR))
cat(sprintf("  - %s/07_summary_barplot.pdf\n", OUTPUT_DIR))
cat(sprintf("  - %s/peak_gene_link_summary.csv\n", OUTPUT_DIR))
cat("\nNote: For coverage plots with fragment files, run signac_05_coverage_plots.R\n")

cat("\n✓ Step 4 complete!\n")
cat(sprintf("  All plots saved to: %s\n", OUTPUT_DIR))
cat("\nAnalysis pipeline complete! 🎉\n\n")
