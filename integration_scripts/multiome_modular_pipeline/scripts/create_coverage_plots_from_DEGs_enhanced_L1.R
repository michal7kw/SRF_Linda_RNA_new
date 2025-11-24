#!/usr/bin/env Rscript

################################################################################
# Enhanced Coverage Plots for Significant DEGs (Combined Controls)
# REFACTORED VERSION - Uses Pre-computed Links
#
# This script creates comprehensive visualizations for genes identified as
# significant in differential expression analysis.
#
# KEY CHANGE: Loads pre-computed LinkPeaks from array job outputs instead of manual cor.test()
#
# PREREQUISITES: Run peak-gene linkage array jobs first:
#   sbatch run_peak_gene_linkage_high_priority_L1.sh  OR
#   sbatch run_peak_gene_linkage_all_L1.sh
#
# Features:
# 1. Per-cell fragment abundance plots by cell type and group
# 2. Coverage plots showing ATAC signal + RNA expression
# 3. Pre-computed peak-gene linkage from Step 3
# 4. Genomic links visualization (arcs connecting peaks to genes)
#
# Input:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/DEG/{celltype}_{genotype}_DEGs.csv
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/peak_gene_links/{celltype}_{genotype}_peak_gene_links.csv
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/integrated_seurat_processed.rds
#
# Output:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/DEG_coverage_plots_enhanced/fragment_abundance/{celltype}_fragment_abundance_combined.pdf
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/DEG_coverage_plots_enhanced/{celltype}_{genotype}/{gene}_coverage_combined.pdf
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/DEG_coverage_plots_enhanced/peak_gene_links/{celltype}_{genotype}_{gene}_links.csv
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
  library(ggforce)  # Required for genomic link visualization
})

# Configuration
RESULTS_DIR <- "signac_results_L1"
SEURAT_FILE <- file.path(RESULTS_DIR, "integrated_seurat_processed.rds")
DEG_DIR <- file.path(RESULTS_DIR, "celltype_results", "DEG")
LINKS_DIR <- file.path(RESULTS_DIR, "celltype_results", "peak_gene_links")
OUTPUT_DIR <- file.path(RESULTS_DIR, "DEG_coverage_plots_enhanced")

# Analysis parameters
PVAL_THRESHOLD <- 0.05       # FDR threshold
LOG2FC_THRESHOLD <- 0.5      # Minimum log2 fold change
TOP_N_GENES <- 20            # Number of top genes to plot per comparison

# Cell types to analyze (L1 annotation)
CELL_TYPES <- c("GABA", "Excitatory", "Oligo", "Astrocytes", "Immune")
GENOTYPES <- c("Nestin", "Emx1")

# Plotting parameters
EXTEND_UPSTREAM <- 100000
EXTEND_DOWNSTREAM <- 100000

# Metadata columns
CELLTYPE_COL <- "cell_type_L1"
CONDITION_COL <- "condition"
GENOTYPE_COL <- "genotype"

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "fragment_abundance"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "peak_gene_links"), recursive = TRUE, showWarnings = FALSE)

cat("================================================================================\n")
cat("ENHANCED DEG COVERAGE PLOTS (Combined Controls - Using Pre-computed Links)\n")
cat("================================================================================\n\n")
cat("KEY CHANGE: Loads pre-computed LinkPeaks from array jobs (NO manual cor.test())\n\n")

################################################################################
# Load Data
################################################################################

cat("Loading integrated Seurat object...\n")
if (!file.exists(SEURAT_FILE)) {
  stop(sprintf("Error: Seurat file not found: %s\n", SEURAT_FILE))
}
seurat_obj <- readRDS(SEURAT_FILE)
cat(sprintf("  ✓ Loaded Seurat object with %d cells\n", ncol(seurat_obj)))

# Cells with condition="Ctrl" → combined_group="Combined_Ctrl"
# Cells with genotype="Nestin" and condition≠"Ctrl" → combined_group="Nestin_Mut"
# Cells with genotype="Emx1" and condition≠"Ctrl" → combined_group="Emx1_Mut"

# Add combined grouping
seurat_obj@meta.data$combined_group <- ifelse(
  seurat_obj@meta.data[[CONDITION_COL]] == "Ctrl",
  "Combined_Ctrl",
  paste(seurat_obj@meta.data[[GENOTYPE_COL]], "Mut", sep = "_")
)

################################################################################
# Load Pre-computed Links from Array Jobs
################################################################################

cat("\nLoading pre-computed LinkPeaks from array job outputs...\n")

if (!dir.exists(LINKS_DIR)) {
  cat("✗ ERROR: Pre-computed links directory not found!\n")
  cat(sprintf("  Expected: %s\n", LINKS_DIR))
  cat("\nYou must run peak-gene linkage array jobs first to generate LinkPeaks.\n")
  cat("  Option 1 (High Priority): sbatch run_peak_gene_linkage_high_priority_L1.sh\n")
  cat("  Option 2 (All Cell Types): sbatch run_peak_gene_linkage_all_L1.sh\n\n")
  stop("Pre-computed links not found. Run array jobs first.")
}

# Load all link files for the cell types we're analyzing
all_links <- list()
for (celltype in CELL_TYPES) {
  for (genotype in GENOTYPES) {
    link_file <- file.path(LINKS_DIR, sprintf("%s_%s_peak_gene_links.csv", celltype, genotype))

    if (file.exists(link_file)) {
      links_df <- read.csv(link_file)
      links_df$celltype <- celltype
      links_df$genotype <- genotype
      all_links[[paste(celltype, genotype, sep = "_")]] <- links_df
      cat(sprintf("  ✓ Loaded %s %s: %d links\n", celltype, genotype, nrow(links_df)))
    } else {
      cat(sprintf("  ⚠ Missing: %s %s\n", celltype, genotype))
    }
  }
}

if (length(all_links) == 0) {
  stop("No pre-computed links found for any cell type. Run array jobs first.")
}

# Combine all links
all_links_df <- do.call(rbind, all_links)
cat(sprintf("\n✓ Total pre-computed links loaded: %d\n", nrow(all_links_df)))
cat(sprintf("  Unique genes with links: %d\n", length(unique(all_links_df$gene))))

################################################################################
# Create Fragment Abundance Plots Per Cell Type
################################################################################

cat("\n================================================================================\n")
cat("CREATING FRAGMENT ABUNDANCE PLOTS\n")
cat("================================================================================\n\n")

if ("nCount_ATAC" %in% colnames(seurat_obj@meta.data)) {
  for (celltype in CELL_TYPES) {
    cat(sprintf("Creating fragment abundance plot for %s...\n", celltype))

    cells_ct <- seurat_obj@meta.data %>%
      filter(.data[[CELLTYPE_COL]] == celltype) %>%
      rownames()

    if (length(cells_ct) == 0) {
      cat(sprintf("  ⚠ No cells found for %s\n", celltype))
      next
    }

    plot_data <- seurat_obj@meta.data[cells_ct, ] %>%
      select(nCount_ATAC, combined_group, all_of(GENOTYPE_COL), all_of(CONDITION_COL)) %>%
      mutate(
        log10_fragments = log10(nCount_ATAC + 1),
        genotype = .data[[GENOTYPE_COL]],
        condition = .data[[CONDITION_COL]]
      )

    # Violin plot by combined group
    p1 <- ggplot(plot_data, aes(x = combined_group, y = log10_fragments, fill = combined_group)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
      scale_fill_manual(values = c(
        "Combined_Ctrl" = "#4CAF50",
        "Nestin_Mut" = "#F44336",
        "Emx1_Mut" = "#FF9800"
      )) +
      labs(
        title = sprintf("ATAC Fragment Counts - %s", celltype),
        subtitle = "Combined Controls (Nestin-Ctrl + Emx1-Ctrl merged)",
        x = "Group",
        y = "log10(Fragment Count + 1)",
        fill = "Group"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )

    # Violin plot by genotype (faceted)
    p2 <- ggplot(plot_data, aes(x = condition, y = log10_fragments, fill = condition)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
      facet_wrap(~ genotype, nrow = 1) +
      scale_fill_manual(values = c("Ctrl" = "#4CAF50", "Mut" = "#F44336")) +
      labs(
        title = "By Genotype",
        x = "Condition",
        y = "log10(Fragment Count + 1)",
        fill = "Condition"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(face = "bold")
      )

    # Save
    pdf_file <- file.path(OUTPUT_DIR, "fragment_abundance",
                         sprintf("%s_fragment_abundance_combined.pdf", celltype))
    ggsave(pdf_file, plot = p1 + p2, width = 16, height = 6, device = "pdf")

    cat(sprintf("  ✓ Saved: %s\n", basename(pdf_file)))
  }
} else {
  cat("  ⚠ Warning: nCount_ATAC not found in metadata\n")
}

################################################################################
# Process Each Cell Type and Genotype
################################################################################

cat("\n================================================================================\n")
cat("PROCESSING DEG RESULTS WITH PRE-COMPUTED PEAK-GENE LINKAGE\n")
cat("================================================================================\n\n")

for (celltype in CELL_TYPES) {
  cat(sprintf("\n=== Cell Type: %s ===\n", celltype))

  for (genotype in GENOTYPES) {
    deg_file <- file.path(DEG_DIR, sprintf("%s_%s_DEGs.csv", celltype, genotype))

    if (!file.exists(deg_file)) {
      cat(sprintf("  ⚠ Skipping %s (DEG file not found)\n", genotype))
      next
    }

    cat(sprintf("\n--- %s Genotype ---\n", genotype))

    # Read DEG results
    deg_results <- read.csv(deg_file)

    # Filter for significant genes
    sig_degs <- deg_results %>%
      filter(p_val_adj < PVAL_THRESHOLD & abs(avg_log2FC) > LOG2FC_THRESHOLD) %>%
      arrange(p_val_adj)

    n_sig <- nrow(sig_degs)
    cat(sprintf("  Found %d significant DEGs\n", n_sig))

    if (n_sig == 0) {
      cat("  ⚠ No significant DEGs, skipping\n")
      next
    }

    # Select top N genes
    top_genes <- head(sig_degs, TOP_N_GENES)
    n_plot <- nrow(top_genes)

    # Subset to cell type
    cells_ct <- seurat_obj@meta.data %>%
      filter(.data[[CELLTYPE_COL]] == celltype) %>%
      rownames()

    if (length(cells_ct) == 0) {
      cat(sprintf("  ⚠ No cells found for cell type %s\n", celltype))
      next
    }

    seurat_subset <- subset(seurat_obj, cells = cells_ct)
    cat(sprintf("  Subset: %d cells\n", ncol(seurat_subset)))

    # Filter pre-computed links to this cell type and genotype
    links_ct_geno <- all_links_df %>%
      filter(celltype == !!celltype & genotype == !!genotype)

    cat(sprintf("  Pre-computed links available: %d (%d unique genes)\n",
                nrow(links_ct_geno),
                length(unique(links_ct_geno$gene))))

    # Create output directory
    output_subdir <- file.path(OUTPUT_DIR, sprintf("%s_%s", celltype, genotype))
    dir.create(output_subdir, recursive = TRUE, showWarnings = FALSE)

    # Create coverage plots with pre-computed links
    for (i in 1:n_plot) {
      gene <- top_genes$gene[i]
      log2fc <- top_genes$avg_log2FC[i]
      pval_adj <- top_genes$p_val_adj[i]
      direction <- top_genes$direction[i]

      cat(sprintf("    [%d/%d] %s (log2FC=%.2f, FDR=%.2e)\n",
                  i, n_plot, gene, log2fc, pval_adj))

      # Skip mitochondrial genes
      if (grepl("^mt-", gene, ignore.case = TRUE)) {
        cat(sprintf("      ⚠ Skipping mitochondrial gene\n"))
        next
      }

      # Check if gene exists in ATAC annotation
      gene_in_atac <- tryCatch({
        coords <- LookupGeneCoords(seurat_subset, gene = gene, assay = "ATAC")
        !is.null(coords) && !any(is.na(coords)) && length(coords) > 0
      }, error = function(e) FALSE)

      if (!gene_in_atac) {
        cat(sprintf("      ⚠ Skipping - gene not found in ATAC annotation\n"))
        next
      }

      # Get pre-computed links for this gene
      gene_links <- links_ct_geno %>% filter(gene == !!gene)

      if (nrow(gene_links) > 0) {
        cat(sprintf("      Found %d pre-computed links (from Step 3)\n", nrow(gene_links)))

        # Save filtered links for this DEG
        write.csv(gene_links,
                  file.path(OUTPUT_DIR, "peak_gene_links",
                           sprintf("%s_%s_%s_links.csv", celltype, genotype, gene)),
                  row.names = FALSE)
      } else {
        cat(sprintf("      ⓘ No pre-computed links found for this gene\n"))
      }

      # Create coverage plot
      tryCatch({
        DefaultAssay(seurat_subset) <- "ATAC"

        # Add genomic links if available
        if (nrow(gene_links) > 0) {
          # Parse peak coordinates
          if ("seqnames" %in% colnames(gene_links)) {
            peak_gr <- GRanges(
              seqnames = gene_links$seqnames,
              ranges = IRanges(
                start = as.numeric(gene_links$start),
                end = as.numeric(gene_links$end)
              ),
              score = gene_links$score
            )
          } else if ("peak" %in% colnames(gene_links)) {
            peak_coords <- do.call(rbind, strsplit(gene_links$peak, "-"))
            peak_gr <- GRanges(
              seqnames = peak_coords[, 1],
              ranges = IRanges(
                start = as.numeric(peak_coords[, 2]),
                end = as.numeric(peak_coords[, 3])
              ),
              score = gene_links$score
            )
          } else {
            cat("      ⚠ Cannot parse peak coordinates\n")
            peak_gr <- NULL
          }

          if (!is.null(peak_gr)) {
            # Get gene coordinates
            gene_coords <- LookupGeneCoords(seurat_subset, gene = gene, assay = "ATAC")
            gene_start <- start(gene_coords)

            # Create Link objects
            link_list <- lapply(seq_along(peak_gr), function(i) {
              peak_center <- (start(peak_gr[i]) + end(peak_gr[i])) / 2
              GRanges(
                seqnames = as.character(seqnames(peak_gr[i])),
                ranges = IRanges(
                  start = min(peak_center, gene_start),
                  end = max(peak_center, gene_start)
                ),
                score = peak_gr[i]$score
              )
            })

            all_link_gr <- do.call(c, link_list)
            Links(seurat_subset[["ATAC"]]) <- all_link_gr
            cat(sprintf("      Added %d genomic links to ATAC assay\n", length(all_link_gr)))
          }
        }

        # Create coverage plot
        p_cov <- CoveragePlot(
          object = seurat_subset,
          region = gene,
          features = gene,
          expression.assay = "RNA",
          assay = "ATAC",
          group.by = "combined_group",
          extend.upstream = EXTEND_UPSTREAM,
          extend.downstream = EXTEND_DOWNSTREAM,
          annotation = TRUE,
          peaks = TRUE,
          links = if(nrow(gene_links) > 0) TRUE else FALSE,
          window = 100
        )

        # Add title
        plot_title <- sprintf("%s - %s %s\nlog2FC = %.2f, FDR = %.2e (%s)\n%d pre-computed links (Step 3)",
                             celltype, genotype, gene, log2fc, pval_adj, direction,
                             nrow(gene_links))

        # Save
        pdf_file <- file.path(output_subdir,
                             sprintf("%02d_%s_coverage_combined.pdf", i, gene))
        ggsave(pdf_file, plot = p_cov + ggtitle(plot_title),
               width = 16, height = 12, device = "pdf")

        cat(sprintf("      ✓ Saved coverage plot\n"))

      }, error = function(e) {
        cat(sprintf("      ✗ Error for %s: %s\n", gene, e$message))
      })
    }

    cat(sprintf("  ✓ Completed %s %s\n", celltype, genotype))
  }
}

################################################################################
# Summary
################################################################################

cat("\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n\n")

# Count files
n_fragment_plots <- length(list.files(file.path(OUTPUT_DIR, "fragment_abundance"), pattern = "\\.pdf$"))
n_coverage_plots <- length(list.files(OUTPUT_DIR, pattern = "coverage.*\\.pdf$", recursive = TRUE))
n_link_files <- length(list.files(file.path(OUTPUT_DIR, "peak_gene_links"), pattern = "\\.csv$"))

cat(sprintf("Fragment abundance plots: %d\n", n_fragment_plots))
cat(sprintf("Coverage plots: %d\n", n_coverage_plots))
cat(sprintf("Peak-gene link files (filtered to DEGs): %d\n", n_link_files))
cat(sprintf("\n✓ Analysis complete!\n"))
cat(sprintf("  Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("\nNOTE: All peak-gene links were pre-computed by Step 3.\n"))
cat(sprintf("      No manual correlation computation was performed.\n"))
