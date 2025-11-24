#!/usr/bin/env Rscript

################################################################################
# POST-PROCESSING FILTERS FOR MULTIOME RESULTS (L1 VERSION)
################################################################################
#
# This script applies recommended filters to address three critical issues:
#
# 1. EXTREME LOG2FC VALUES in DA results (|log2FC| > 10)
#    - Caused by binary ON/OFF peaks (0% in one group)
#    - Solution: Filter peaks with min.pct threshold
#
# 2. PERMISSIVE CORRELATION THRESHOLD in peak-gene links (r > 0.05)
#    - Current threshold explains only 0.25% of variance
#    - Solution: Filter to r > 0.2 (explains 4% of variance)
#
# 3. NO FDR CORRECTION on peak-gene correlations
#    - Millions of tests → many false positives
#    - Solution: Apply Benjamini-Hochberg FDR correction
#
# L1 VERSION: Filters results from signac_results_L1/ (9 broad categories)
#
# Usage:
#   Rscript POST_PROCESSING_FILTERS_L1.R
#
# Output:
#   - Creates signac_results_L1/celltype_results/filtered/ subdirectory
#   - Generates comparison plots (before/after filtering)
#   - Produces summary statistics
#
# Inputs:
#   - Differential accessibility (DA) results from Step 3:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/DA/{celltype}_{genotype}_DA_peaks.csv
#   - Peak-gene linkage results from Step 3:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/peak_gene_links/{celltype}_{genotype}_peak_gene_links.csv
#
# Outputs:
#   - Filtered DA results:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/filtered/{celltype}_{genotype}_DA_peaks.csv
#   - Filtered peak-gene links:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/filtered/{celltype}_{genotype}_peak_gene_links.csv
#   - DA filtering summary:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/filtered/DA_filtering_summary.csv
#   - Links filtering summary:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/filtered/links_filtering_summary.csv
#   - DA filtering comparison plots:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/filtered/DA_filtering_comparison.pdf
#   - Correlation distribution plots:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/filtered/correlation_distribution.pdf
#
################################################################################

library(dplyr)
library(ggplot2)
library(patchwork)

# Configuration
BASE_DIR <- "signac_results_L1/celltype_results"
OUTPUT_DIR <- "signac_results_L1/celltype_results/filtered"

# Filter thresholds
MIN_PCT_THRESHOLD <- 0.05     # Peaks must be accessible in ≥5% of cells (at least one group)
MIN_PCT_BOTH_GROUPS <- 0.01   # OR ≥1% in both groups (less stringent for binary peaks)
MAX_LOG2FC_DA <- 10           # Maximum biologically plausible log2FC for chromatin (increased from 8)
MIN_CORRELATION <- 0.2        # Minimum correlation for peak-gene links
FDR_THRESHOLD <- 0.05         # FDR q-value threshold

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

################################################################################
# FUNCTION: Filter DA Results
################################################################################

filter_da_results <- function(da_file, celltype) {
  cat(sprintf("\n=== Filtering DA results: %s ===\n", celltype))

  # Read DA results
  da <- read.csv(da_file, stringsAsFactors = FALSE)
  n_original <- nrow(da)

  cat(sprintf("  Original DA peaks: %d\n", n_original))

  # Diagnostic: Check for extreme values
  n_extreme <- sum(abs(da$avg_log2FC) > MAX_LOG2FC_DA, na.rm = TRUE)
  if (n_extreme > 0) {
    cat(sprintf("  ⚠ Found %d peaks with |log2FC| > %d (%.1f%%)\n",
                n_extreme, MAX_LOG2FC_DA, 100 * n_extreme / n_original))

    # Show examples of extreme peaks
    extreme_peaks <- da[abs(da$avg_log2FC) > MAX_LOG2FC_DA, ]
    extreme_peaks <- extreme_peaks[order(-abs(extreme_peaks$avg_log2FC)), ]
    cat("    Example extreme peaks:\n")
    for (i in 1:min(3, nrow(extreme_peaks))) {
      cat(sprintf("      %s: log2FC=%.1f, pct.1=%.1f%%, pct.2=%.1f%%\n",
                  extreme_peaks$peak[i],
                  extreme_peaks$avg_log2FC[i],
                  extreme_peaks$pct.1[i] * 100,
                  extreme_peaks$pct.2[i] * 100))
    }
  }

  # Apply filters (relaxed for scATAC-seq)
  da_filtered <- da %>%
    filter(
      # Remove extreme log2FC (likely technical artifacts)
      abs(avg_log2FC) <= MAX_LOG2FC_DA,

      # At least one group must have ≥5% cells with accessible peak
      # OR both groups have ≥1% (allows condition-specific peaks)
      (pct.1 >= MIN_PCT_THRESHOLD | pct.2 >= MIN_PCT_THRESHOLD) &
      (pct.1 >= MIN_PCT_BOTH_GROUPS & pct.2 >= MIN_PCT_BOTH_GROUPS),

      # Keep significance threshold
      p_val_adj < 0.05
    )

  n_filtered <- nrow(da_filtered)
  pct_retained <- 100 * n_filtered / n_original

  cat(sprintf("  Filtered DA peaks: %d (%.1f%% retained)\n",
              n_filtered, pct_retained))
  cat(sprintf("  Removed: %d peaks (%.1f%%)\n",
              n_original - n_filtered, 100 - pct_retained))

  # Breakdown of what was removed
  n_extreme_removed <- sum(abs(da$avg_log2FC) > MAX_LOG2FC_DA, na.rm = TRUE)
  n_low_pct <- sum(da$pct.1 < MIN_PCT_THRESHOLD | da$pct.2 < MIN_PCT_THRESHOLD, na.rm = TRUE)

  cat(sprintf("    - Extreme log2FC: %d peaks\n", n_extreme_removed))
  cat(sprintf("    - Low percentage: %d peaks\n", n_low_pct))

  # Save filtered results
  output_file <- file.path(OUTPUT_DIR, basename(da_file))
  write.csv(da_filtered, output_file, row.names = FALSE, quote = TRUE)
  cat(sprintf("  ✓ Saved: %s\n", output_file))

  # Return summary stats
  return(list(
    celltype = celltype,
    original = n_original,
    filtered = n_filtered,
    pct_retained = pct_retained,
    n_extreme = n_extreme_removed,
    n_low_pct = n_low_pct
  ))
}

################################################################################
# FUNCTION: Filter Peak-Gene Links
################################################################################

filter_peak_gene_links <- function(link_file, celltype) {
  cat(sprintf("\n=== Filtering peak-gene links: %s ===\n", celltype))

  # Read links
  links <- read.csv(link_file, stringsAsFactors = FALSE)
  n_original <- nrow(links)

  cat(sprintf("  Original links: %d\n", n_original))

  # Apply FDR correction
  links$qvalue <- p.adjust(links$pvalue, method = "BH")

  # Check current correlation distribution
  cat(sprintf("  Correlation range: %.3f to %.3f\n",
              min(links$score, na.rm = TRUE),
              max(links$score, na.rm = TRUE)))
  cat(sprintf("  Median |correlation|: %.3f\n",
              median(abs(links$score), na.rm = TRUE)))

  # Count links below threshold
  n_weak <- sum(abs(links$score) < MIN_CORRELATION, na.rm = TRUE)
  cat(sprintf("  ⚠ Links with |r| < %.2f: %d (%.1f%%)\n",
              MIN_CORRELATION, n_weak, 100 * n_weak / n_original))

  # Apply filters
  links_filtered <- links %>%
    filter(
      # Stricter correlation threshold
      abs(score) >= MIN_CORRELATION,

      # FDR-corrected p-value
      qvalue < FDR_THRESHOLD
    )

  n_filtered <- nrow(links_filtered)
  pct_retained <- 100 * n_filtered / n_original

  cat(sprintf("  Filtered links: %d (%.1f%% retained)\n",
              n_filtered, pct_retained))
  cat(sprintf("  Removed: %d links (%.1f%%)\n",
              n_original - n_filtered, 100 - pct_retained))

  # Breakdown
  n_weak_removed <- sum(abs(links$score) < MIN_CORRELATION, na.rm = TRUE)
  n_fdr_fail <- sum(links$qvalue >= FDR_THRESHOLD, na.rm = TRUE)

  cat(sprintf("    - Weak correlation (|r| < %.2f): %d links\n",
              MIN_CORRELATION, n_weak_removed))
  cat(sprintf("    - Failed FDR (q ≥ %.2f): %d links\n",
              FDR_THRESHOLD, n_fdr_fail))

  # Save filtered results
  output_file <- file.path(OUTPUT_DIR, basename(link_file))
  write.csv(links_filtered, output_file, row.names = FALSE, quote = TRUE)
  cat(sprintf("  ✓ Saved: %s\n", output_file))

  # Return summary stats
  return(list(
    celltype = celltype,
    original = n_original,
    filtered = n_filtered,
    pct_retained = pct_retained,
    n_weak = n_weak_removed,
    n_fdr = n_fdr_fail
  ))
}

################################################################################
# MAIN EXECUTION
################################################################################

cat("\n================================================================================\n")
cat("POST-PROCESSING FILTERS FOR MULTIOME RESULTS\n")
cat("================================================================================\n")
cat("\nConfiguration:\n")
cat(sprintf("  Min percentage (DA): %.1f%% (at least one group) OR %.1f%% (both groups)\n",
            MIN_PCT_THRESHOLD * 100, MIN_PCT_BOTH_GROUPS * 100))
cat(sprintf("  Max |log2FC| (DA): %d\n", MAX_LOG2FC_DA))
cat(sprintf("  Min |correlation| (links): %.2f\n", MIN_CORRELATION))
cat(sprintf("  FDR threshold: %.2f\n", FDR_THRESHOLD))
cat("\nNote: Filters relaxed for scATAC-seq to retain condition-specific peaks\n")

# Process DA results
cat("\n################################################################################\n")
cat("FILTERING DIFFERENTIAL ACCESSIBILITY RESULTS\n")
cat("################################################################################\n")

da_files <- list.files(file.path(BASE_DIR, "DA"), pattern = "_(Nestin|Emx1)_DA_peaks\\.csv$", full.names = TRUE)
da_stats <- list()

for (da_file in da_files) {
  # Extract celltype and genotype from filename
  # Format: {celltype}_{genotype}_DA_peaks.csv
  analysis_id <- gsub("_DA_peaks\\.csv$", "", basename(da_file))
  da_stats[[analysis_id]] <- filter_da_results(da_file, analysis_id)
}

# Process peak-gene links
cat("\n################################################################################\n")
cat("FILTERING PEAK-GENE LINKAGE RESULTS\n")
cat("################################################################################\n")

link_files <- list.files(file.path(BASE_DIR, "peak_gene_links"),
                         pattern = "_(Nestin|Emx1)_peak_gene_links\\.csv$", full.names = TRUE)
link_stats <- list()

for (link_file in link_files) {
  # Extract celltype and genotype from filename
  # Format: {celltype}_{genotype}_peak_gene_links.csv
  analysis_id <- gsub("_peak_gene_links\\.csv$", "", basename(link_file))
  link_stats[[analysis_id]] <- filter_peak_gene_links(link_file, analysis_id)
}

################################################################################
# SUMMARY STATISTICS
################################################################################

cat("\n================================================================================\n")
cat("SUMMARY STATISTICS\n")
cat("================================================================================\n")

# DA summary
cat("\n=== Differential Accessibility ===\n\n")

# Convert list of stats to proper data frame
da_summary_df <- data.frame(
  celltype = sapply(da_stats, function(x) x$celltype),
  original = sapply(da_stats, function(x) x$original),
  filtered = sapply(da_stats, function(x) x$filtered),
  pct_retained = sapply(da_stats, function(x) x$pct_retained),
  n_extreme = sapply(da_stats, function(x) x$n_extreme),
  n_low_pct = sapply(da_stats, function(x) x$n_low_pct),
  stringsAsFactors = FALSE
)

cat(sprintf("%-15s %10s %10s %10s\n", "Cell Type", "Original", "Filtered", "% Retained"))
cat(strrep("-", 50), "\n")
for (i in 1:nrow(da_summary_df)) {
  cat(sprintf("%-15s %10d %10d %10.1f%%\n",
              da_summary_df$celltype[i],
              da_summary_df$original[i],
              da_summary_df$filtered[i],
              da_summary_df$pct_retained[i]))
}
cat(strrep("-", 50), "\n")
cat(sprintf("%-15s %10d %10d %10.1f%%\n",
            "TOTAL",
            sum(da_summary_df$original),
            sum(da_summary_df$filtered),
            100 * sum(da_summary_df$filtered) / sum(da_summary_df$original)))

# Links summary
cat("\n=== Peak-Gene Links ===\n\n")

# Convert list of stats to proper data frame
link_summary_df <- data.frame(
  celltype = sapply(link_stats, function(x) x$celltype),
  original = sapply(link_stats, function(x) x$original),
  filtered = sapply(link_stats, function(x) x$filtered),
  pct_retained = sapply(link_stats, function(x) x$pct_retained),
  n_weak = sapply(link_stats, function(x) x$n_weak),
  n_fdr = sapply(link_stats, function(x) x$n_fdr),
  stringsAsFactors = FALSE
)

cat(sprintf("%-15s %10s %10s %10s\n", "Cell Type", "Original", "Filtered", "% Retained"))
cat(strrep("-", 50), "\n")
for (i in 1:nrow(link_summary_df)) {
  cat(sprintf("%-15s %10d %10d %10.1f%%\n",
              link_summary_df$celltype[i],
              link_summary_df$original[i],
              link_summary_df$filtered[i],
              link_summary_df$pct_retained[i]))
}
cat(strrep("-", 50), "\n")
cat(sprintf("%-15s %10d %10d %10.1f%%\n",
            "TOTAL",
            sum(link_summary_df$original),
            sum(link_summary_df$filtered),
            100 * sum(link_summary_df$filtered) / sum(link_summary_df$original)))

# Save summary tables
write.csv(da_summary_df, file.path(OUTPUT_DIR, "DA_filtering_summary.csv"),
          row.names = FALSE)
write.csv(link_summary_df, file.path(OUTPUT_DIR, "links_filtering_summary.csv"),
          row.names = FALSE)

################################################################################
# VISUALIZATION: Before/After Comparison
################################################################################

cat("\n================================================================================\n")
cat("GENERATING COMPARISON PLOTS\n")
cat("================================================================================\n")

# Dynamically select an example cell type with sufficient data
# Prioritize cell types with both DA peaks and peak-gene links
example_analysis_id <- NULL

# First, try to find a cell type with both DA peaks and peak-gene links
if (length(link_stats) > 0) {
  # Get cell types with links that have reasonable retention
  link_celltypes <- names(link_stats)
  link_retention <- sapply(link_stats, function(x) x$pct_retained)

  # Find cell types with >5% retention and >20 original links
  good_link_celltypes <- link_celltypes[
    link_retention > 5 &
    sapply(link_stats, function(x) x$original) > 20
  ]

  if (length(good_link_celltypes) > 0) {
    # Pick the one with most original links
    link_counts <- sapply(link_stats[good_link_celltypes], function(x) x$original)
    example_analysis_id <- good_link_celltypes[which.max(link_counts)]
    cat(sprintf("\n  Selected '%s' as example (has %d peak-gene links)\n",
                example_analysis_id, max(link_counts)))
  }
}

# If no good link examples, fall back to DA peaks only
if (is.null(example_analysis_id) && length(da_stats) > 0) {
  # Get cell types with reasonable DA retention
  da_celltypes <- names(da_stats)
  da_retention <- sapply(da_stats, function(x) x$pct_retained)

  # Find cell types with >5% retention and >100 original peaks
  good_da_celltypes <- da_celltypes[
    da_retention > 5 &
    sapply(da_stats, function(x) x$original) > 100
  ]

  if (length(good_da_celltypes) > 0) {
    # Pick the one with most filtered peaks
    da_filtered_counts <- sapply(da_stats[good_da_celltypes], function(x) x$filtered)
    example_analysis_id <- good_da_celltypes[which.max(da_filtered_counts)]
    cat(sprintf("\n  Selected '%s' as example (has %d filtered DA peaks, but no links)\n",
                example_analysis_id, max(da_filtered_counts)))
  }
}

# Proceed with visualization if we found a suitable example
da_file_example <- if (!is.null(example_analysis_id)) {
  file.path(BASE_DIR, "DA", paste0(example_analysis_id, "_DA_peaks.csv"))
} else {
  cat("\n  ⚠ Warning: No suitable cell types found for visualization\n")
  cat("    Skipping comparison plots.\n")
  ""
}

if (!is.null(example_analysis_id) && file.exists(da_file_example)) {
  da_orig <- read.csv(da_file_example)
  da_filt <- read.csv(file.path(OUTPUT_DIR, paste0(example_analysis_id, "_DA_peaks.csv")))

  # Add filter status to original data
  da_orig$filtered_out <- !(da_orig$peak %in% da_filt$peak)

  # Create comparison plot
  p1 <- ggplot(da_orig, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = filtered_out), alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("TRUE" = "gray70", "FALSE" = "red3"),
                       labels = c("Retained", "Removed"),
                       name = "Status") +
    geom_vline(xintercept = c(-MAX_LOG2FC_DA, MAX_LOG2FC_DA),
               linetype = "dashed", color = "blue") +
    labs(title = paste0(example_analysis_id, " - Before Filtering"),
         subtitle = sprintf("%d peaks total", nrow(da_orig)),
         x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)") +
    theme_minimal() +
    theme(legend.position = "bottom")

  p2 <- ggplot(da_filt, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = direction), alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("More_Accessible" = "red3",
                                   "Less_Accessible" = "blue3")) +
    labs(title = paste0(example_analysis_id, " - After Filtering"),
         subtitle = sprintf("%d peaks retained (%.1f%%)",
                           nrow(da_filt),
                           100 * nrow(da_filt) / nrow(da_orig)),
         x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)") +
    theme_minimal() +
    theme(legend.position = "bottom")

  p_combined <- p1 + p2

  ggsave(file.path(OUTPUT_DIR, "DA_filtering_comparison.pdf"),
         p_combined, width = 12, height = 6)
  cat("  ✓ Saved: filtered/DA_filtering_comparison.pdf\n")
}

# Correlation distribution comparison
link_file_example <- if (!is.null(example_analysis_id)) {
  file.path(BASE_DIR, "peak_gene_links",
            paste0(example_analysis_id, "_peak_gene_links.csv"))
} else {
  ""
}

if (!is.null(example_analysis_id) && file.exists(link_file_example)) {
  links_orig <- read.csv(link_file_example)
  links_filt_path <- file.path(OUTPUT_DIR,
                                paste0(example_analysis_id, "_peak_gene_links.csv"))

  if (file.exists(links_filt_path)) {
    links_filt <- read.csv(links_filt_path)

    # Add FDR if not present in original
    if (!"qvalue" %in% colnames(links_orig)) {
      links_orig$qvalue <- p.adjust(links_orig$pvalue, method = "BH")
    }

    # Add filter status
    links_orig$filtered_out <- !(paste(links_orig$peak, links_orig$gene) %in%
                                  paste(links_filt$peak, links_filt$gene))

    p3 <- ggplot(links_orig, aes(x = abs(score))) +
      geom_histogram(aes(fill = filtered_out), bins = 50, alpha = 0.7) +
      geom_vline(xintercept = MIN_CORRELATION, linetype = "dashed", color = "red") +
      scale_fill_manual(values = c("TRUE" = "gray70", "FALSE" = "blue3"),
                        labels = c("Retained", "Removed"),
                        name = "Status") +
      labs(title = paste0(example_analysis_id, " - Correlation Distribution"),
           subtitle = sprintf("Original: %d links, Filtered: %d links (%.1f%%)",
                             nrow(links_orig), nrow(links_filt),
                             100 * nrow(links_filt) / nrow(links_orig)),
           x = "Absolute Correlation (|r|)", y = "Count") +
      theme_minimal() +
      theme(legend.position = "bottom")

    ggsave(file.path(OUTPUT_DIR, "correlation_distribution.pdf"),
           p3, width = 8, height = 6)
    cat("  ✓ Saved: filtered/correlation_distribution.pdf\n")
  } else {
    cat(sprintf("\n  ⚠ Warning: No filtered links found for '%s'\n", example_analysis_id))
    cat("    Skipping correlation distribution plot.\n")
  }
} else {
  if (length(link_stats) == 0) {
    cat("\n  ⚠ Warning: No peak-gene links were processed\n")
    cat("    Skipping correlation distribution plot.\n")
    cat("    This may indicate that Step 3 did not generate peak-gene links.\n")
  }
}

cat("\n================================================================================\n")
cat("FILTERING COMPLETE\n")
cat("================================================================================\n")
cat("\nFiltered results saved to:\n")
cat(sprintf("  %s/\n", OUTPUT_DIR))
cat("\nKey files:\n")
cat("  - DA_filtering_summary.csv\n")
cat("  - links_filtering_summary.csv\n")
cat("  - DA_filtering_comparison.pdf\n")
cat("  - correlation_distribution.pdf\n")
cat("\nNext steps:\n")
cat("  1. Review summary statistics above\n")
cat("  2. Examine comparison plots\n")
cat("  3. Use filtered/ files for downstream analysis\n")
cat("  4. Repeat INTERPRETING_MULTIOME_RESULTS.md analyses with filtered data\n")
cat("\n================================================================================\n")
