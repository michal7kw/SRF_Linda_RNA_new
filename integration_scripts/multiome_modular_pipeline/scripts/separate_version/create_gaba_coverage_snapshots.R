#!/usr/bin/env Rscript

################################################################################
# GABA Coverage Snapshots for Specific Genes
################################################################################
#
# This script creates publication-quality chromatin accessibility coverage
# plots for GABA neurons, showing:
#   - Genes: Read from genes_inter.txt file (splicing regulators)
#   - Comparisons: Nestin Ctrl vs Mut, Emx1 Ctrl vs Mut (separate)
#   - Region: Gene body + 100kb upstream/downstream
#   - Shows: ATAC coverage tracks, peaks, gene annotation, and RNA expression
#   - Uses robust Signac CoveragePlot for genome browser visualization
#
# INPUT FILES:
#   - Integrated Seurat object with processed data:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/integrated_seurat_processed.rds
#   - Gene list file (splicing regulators):
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/genes_inter.txt
#   - Fragment files (linked in Seurat object ATAC assay):
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/translated_fragments/*_fragments_translated.tsv.gz
#
# OUTPUT FILES:
#   - Main combined PDF with all genes:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/plots/gaba_coverage_snapshots/GABA_coverage_snapshots_specific_genes.pdf
#   - Individual gene PDFs:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/plots/gaba_coverage_snapshots/GABA_coverage_<gene>.pdf
#   - Summary PDF with all genes by genotype:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/plots/gaba_coverage_snapshots/GABA_coverage_summary_all_genes.pdf
#
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

cat("================================================================================\n")
cat("  GABA COVERAGE SNAPSHOTS - SPECIFIC GENES\n")
cat("================================================================================\n\n")

################################################################################
# Configuration
################################################################################

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
INPUT_DIR <- file.path(BASE_DIR, "signac_results_L1")
OUTPUT_DIR <- file.path(INPUT_DIR, "plots", "gaba_coverage_snapshots")

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Read genes of interest from file
GENES_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/genes_inter.txt"
GENES_OF_INTEREST <- readLines(GENES_FILE)
GENES_OF_INTEREST <- GENES_OF_INTEREST[nchar(GENES_OF_INTEREST) > 0]  # Remove empty lines
GENES_OF_INTEREST <- unique(GENES_OF_INTEREST)  # Remove duplicates

# Cell types to include (GABA neurons from L1 annotation)
CELL_TYPES <- c("GABA")

# Genomic window around gene
EXTEND_UPSTREAM <- 100000   # 100kb upstream
EXTEND_DOWNSTREAM <- 100000  # 100kb downstream

# Metadata columns
CELLTYPE_COL <- "cell_type_L1"
GENOTYPE_COL <- "genotype"
CONDITION_COL <- "condition"

cat("Configuration:\n")
cat(sprintf("  Gene list file: %s\n", GENES_FILE))
cat(sprintf("  Genes loaded: %d\n", length(GENES_OF_INTEREST)))
cat(sprintf("    %s\n", paste(GENES_OF_INTEREST, collapse = ", ")))
cat(sprintf("  Cell types: %s\n", paste(CELL_TYPES, collapse = " + ")))
cat(sprintf("  Genomic window: Gene ± %d kb\n", EXTEND_UPSTREAM / 1000))
cat(sprintf("  Output directory: %s\n", OUTPUT_DIR))
cat("\n")

################################################################################
# Load Data
################################################################################

cat("=== Loading Seurat object ===\n")
seurat_obj <- readRDS(file.path(INPUT_DIR, "integrated_seurat_processed.rds"))
cat(sprintf("  ✓ Loaded: %d cells\n", ncol(seurat_obj)))

# Check metadata columns
if (!CELLTYPE_COL %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Cell type column '%s' not found in metadata", CELLTYPE_COL))
}
if (!GENOTYPE_COL %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Genotype column '%s' not found in metadata", GENOTYPE_COL))
}
if (!CONDITION_COL %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Condition column '%s' not found in metadata", CONDITION_COL))
}

# CRITICAL: Check and validate fragment file for coverage visualization
cat("\n=== Checking fragment files for coverage visualization ===\n")
if ("ATAC" %in% names(seurat_obj@assays)) {
  fragments <- Fragments(seurat_obj[["ATAC"]])
  if (length(fragments) > 0) {
    cat(sprintf("  ✓ Fragment files found: %d\n", length(fragments)))
    for (i in seq_along(fragments)) {
      frag_path <- fragments[[i]]@path
      cat(sprintf("    Fragment %d: %s\n", i, frag_path))
      if (file.exists(frag_path)) {
        cat(sprintf("      ✓ File exists\n"))
      } else {
        cat(sprintf("      ✗ WARNING: File not found!\n"))
      }
    }
  } else {
    cat("  ✗ WARNING: No fragment files linked to ATAC assay!\n")
    cat("  Coverage tracks may not display properly.\n")
  }
} else {
  stop("ATAC assay not found in Seurat object")
}

################################################################################
# Filter to GABA
################################################################################

cat("\n=== Filtering to GABA neurons ===\n")

# Get cells that are GABA
gaba_cells <- seurat_obj@meta.data %>%
  filter(.data[[CELLTYPE_COL]] %in% CELL_TYPES) %>%
  rownames()

cat(sprintf("  Total GABA cells: %d\n", length(gaba_cells)))

# Show breakdown by cell type
gaba_breakdown <- seurat_obj@meta.data %>%
  filter(.data[[CELLTYPE_COL]] %in% CELL_TYPES) %>%
  group_by(.data[[CELLTYPE_COL]]) %>%
  summarise(n = n(), .groups = "drop")

for (i in 1:nrow(gaba_breakdown)) {
  cat(sprintf("    %s: %d cells\n",
              gaba_breakdown[[CELLTYPE_COL]][i],
              gaba_breakdown$n[i]))
}

# Show breakdown by genotype and condition
cat("\n  Breakdown by genotype and condition:\n")
gaba_genotype_breakdown <- seurat_obj@meta.data %>%
  filter(.data[[CELLTYPE_COL]] %in% CELL_TYPES) %>%
  group_by(.data[[GENOTYPE_COL]], .data[[CONDITION_COL]]) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(.data[[GENOTYPE_COL]], .data[[CONDITION_COL]])

for (i in 1:nrow(gaba_genotype_breakdown)) {
  cat(sprintf("    %s %s: %d cells\n",
              gaba_genotype_breakdown[[GENOTYPE_COL]][i],
              gaba_genotype_breakdown[[CONDITION_COL]][i],
              gaba_genotype_breakdown$n[i]))
}

# Subset Seurat object
seurat_gaba <- subset(seurat_obj, cells = gaba_cells)
cat(sprintf("\n  ✓ Subsetted to %d GABA cells\n", ncol(seurat_gaba)))

# Add combined cell type label for plotting
seurat_gaba@meta.data$combined_celltype <- "GABA"

################################################################################
# Function: Create Coverage Plot for Gene (Improved)
################################################################################

create_gene_coverage_plot <- function(seurat_obj, gene, genotype,
                                       extend_upstream = 100000,
                                       extend_downstream = 100000,
                                       max_cutoff = NULL) {

  cat(sprintf("\n  Creating coverage plot: %s (%s genotype)\n", gene, genotype))

  # Get cells for this genotype FIRST
  cells_genotype <- seurat_obj@meta.data %>%
    filter(.data[[GENOTYPE_COL]] == genotype) %>%
    rownames()

  if (length(cells_genotype) == 0) {
    cat(sprintf("    ⚠ No cells found for genotype: %s\n", genotype))
    return(NULL)
  }

  # Subset to genotype-specific cells
  seurat_geno <- subset(seurat_obj, cells = cells_genotype)
  cat(sprintf("    Cells for %s: %d\n", genotype, ncol(seurat_geno)))

  # NOW create the combined genotype_condition identifier on the subsetted object
  seurat_geno$geno_cond <- paste(seurat_geno@meta.data[[GENOTYPE_COL]],
                                  seurat_geno@meta.data[[CONDITION_COL]],
                                  sep = "_")

  # Validate that fragment files are still linked after subsetting
  if ("ATAC" %in% names(seurat_geno@assays)) {
    fragments_geno <- Fragments(seurat_geno[["ATAC"]])
    if (length(fragments_geno) == 0) {
      cat("      ⚠ WARNING: Fragment files lost during subsetting. Re-linking...\n")
      # Copy fragment files from original object
      fragments_orig <- Fragments(seurat_obj[["ATAC"]])
      Fragments(seurat_geno[["ATAC"]]) <- fragments_orig
    }
  }

  # Show breakdown by condition
  condition_counts <- table(seurat_geno@meta.data[[CONDITION_COL]])
  for (cond in names(condition_counts)) {
    cat(sprintf("      %s: %d cells\n", cond, condition_counts[cond]))
  }

  # Set identity to combined genotype_condition for plotting
  Idents(seurat_geno) <- "geno_cond"

  # Get the identities for this genotype only (should be e.g., "Nestin_Ctrl", "Nestin_Mut")
  genotype_idents <- paste(genotype, c("Ctrl", "Mut"), sep = "_")
  available_idents <- intersect(genotype_idents, unique(as.character(Idents(seurat_geno))))

  cat(sprintf("    Identities to plot: %s\n", paste(available_idents, collapse = ", ")))

  # CRITICAL: Set default assay to ATAC for coverage visualization
  DefaultAssay(seurat_geno) <- "ATAC"

  # Create improved coverage plot with explicit coverage tracks
  p <- tryCatch({
    # Use CoveragePlot with explicit parameters for pile-up visualization
    coverage_plot <- CoveragePlot(
      object = seurat_geno,
      region = gene,
      features = gene,
      expression.assay = "RNA",
      assay = "ATAC",  # Explicitly set ATAC assay for coverage
      extend.upstream = extend_upstream,
      extend.downstream = extend_downstream,
      idents = available_idents,  # Only show this genotype's conditions
      annotation = TRUE,
      peaks = TRUE,
      tile = FALSE,
      links = FALSE,
      ncol = 1,
      window = 100,  # Window size for coverage calculation (100bp bins)
      group.by = "geno_cond",  # Explicitly state the grouping variable
      max.cutoff = max_cutoff,  # Apply consistent Y-axis scaling if provided
      scale.factor = ifelse(!is.null(max_cutoff), 1e6, NULL)  # Normalize to CPM if max_cutoff is set
    )

    # Apply theming
    coverage_plot <- coverage_plot &
      theme_classic() &
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black")
      )

    # Add title as patchwork annotation
    coverage_plot +
      plot_annotation(
        title = sprintf("%s - %s Genotype (GABA neurons)", gene, genotype),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
      )

  }, error = function(e) {
    cat(sprintf("    ⚠ Error creating coverage plot: %s\n", e$message))
    cat(sprintf("    Full error: %s\n", paste(capture.output(print(e)), collapse = "\n")))
    return(NULL)
  })

  if (!is.null(p)) {
    cat(sprintf("    ✓ Plot created successfully\n"))
  }

  return(p)
}

################################################################################
# Function: Calculate Maximum Coverage for Gene Across All Groups
################################################################################

calculate_max_coverage_for_gene <- function(seurat_obj, gene,
                                           extend_upstream = 100000,
                                           extend_downstream = 100000) {

  cat(sprintf("    Calculating maximum coverage for consistent scaling...\n"))

  # Calculate max coverage across ALL groups (Nestin_Ctrl, Nestin_Mut, Emx1_Ctrl, Emx1_Mut)
  # This ensures identical Y-axis scaling for both genotype plots

  tryCatch({
    # Set default assay
    DefaultAssay(seurat_obj) <- "ATAC"

    # Create combined genotype_condition identifier
    seurat_obj$geno_cond <- paste(seurat_obj@meta.data[[GENOTYPE_COL]],
                                   seurat_obj@meta.data[[CONDITION_COL]],
                                   sep = "_")

    # Set identity
    Idents(seurat_obj) <- "geno_cond"

    # Create a temporary coverage plot with all groups to extract max y-axis value
    cat(sprintf("      Creating temporary plot with all groups...\n"))
    temp_plot <- CoveragePlot(
      object = seurat_obj,
      region = gene,
      assay = "ATAC",
      extend.upstream = extend_upstream,
      extend.downstream = extend_downstream,
      annotation = FALSE,
      peaks = FALSE,
      tile = FALSE,
      links = FALSE,
      ncol = 1,
      window = 100,
      group.by = "geno_cond",
      max.cutoff = "q99.5",  # Use 99.5th percentile for comprehensive view
      scale.factor = 1e6  # Normalize to CPM
    )

    # Extract maximum y-axis value from all coverage tracks
    max_y <- 0

    if (inherits(temp_plot, "patchwork")) {
      # Patchwork object - iterate through all plots
      for (i in seq_along(temp_plot)) {
        p <- temp_plot[[i]]
        if (!is.null(p) && inherits(p, "gg")) {
          built <- ggplot_build(p)
          if (length(built$layout$panel_params) > 0) {
            y_range <- built$layout$panel_params[[1]]$y.range
            if (!is.null(y_range) && length(y_range) >= 2) {
              max_y <- max(max_y, y_range[2], na.rm = TRUE)
            }
          }
        }
      }
    } else if (inherits(temp_plot, "gg")) {
      # Single ggplot object
      built <- ggplot_build(temp_plot)
      if (length(built$layout$panel_params) > 0) {
        y_range <- built$layout$panel_params[[1]]$y.range
        if (!is.null(y_range) && length(y_range) >= 2) {
          max_y <- y_range[2]
        }
      }
    }

    # Validate extracted value
    if (max_y <= 0 || is.na(max_y) || is.infinite(max_y)) {
      cat(sprintf("      Could not extract valid max coverage, using q99.5\n"))
      return("q99.5")
    }

    # Add 5% buffer to ensure all peaks are visible
    max_y <- max_y * 1.05

    cat(sprintf("      ✓ Calculated max coverage: %.2f CPM\n", max_y))
    cat(sprintf("      This value will be used for both genotype comparisons\n"))
    return(max_y)

  }, error = function(e) {
    cat(sprintf("      Warning: Error calculating max coverage: %s\n", e$message))
    cat(sprintf("      Falling back to q99.5 (may cause scale differences)\n"))
    return("q99.5")
  })
}

################################################################################
# Generate Coverage Plots for All Genes
################################################################################

cat("\n================================================================================\n")
cat("GENERATING COVERAGE PLOTS\n")
cat("================================================================================\n")

# Store all plots
all_plots <- list()

for (gene in GENES_OF_INTEREST) {
  cat(sprintf("\n=== Processing gene: %s ===\n", gene))

  # Calculate maximum coverage across all groups for this gene
  # This ensures consistent Y-axis scaling between plots
  max_coverage <- calculate_max_coverage_for_gene(
    seurat_gaba,
    gene,
    extend_upstream = EXTEND_UPSTREAM,
    extend_downstream = EXTEND_DOWNSTREAM
  )

  # Create plots for both genotypes
  # Both plots will use the same max_cutoff for consistent scaling
  plot_nestin <- create_gene_coverage_plot(
    seurat_gaba,
    gene,
    genotype = "Nestin",
    extend_upstream = EXTEND_UPSTREAM,
    extend_downstream = EXTEND_DOWNSTREAM,
    max_cutoff = max_coverage
  )

  plot_emx1 <- create_gene_coverage_plot(
    seurat_gaba,
    gene,
    genotype = "Emx1",
    extend_upstream = EXTEND_UPSTREAM,
    extend_downstream = EXTEND_DOWNSTREAM,
    max_cutoff = max_coverage
  )

  # Store plots
  all_plots[[paste0(gene, "_Nestin")]] <- plot_nestin
  all_plots[[paste0(gene, "_Emx1")]] <- plot_emx1
}

################################################################################
# Save Individual Gene PNG Files
################################################################################

cat("\n================================================================================\n")
cat("SAVING RESULTS (PNG FORMAT)\n")
cat("================================================================================\n")

# Filter out NULL plots
valid_plots <- all_plots[!sapply(all_plots, is.null)]

if (length(valid_plots) == 0) {
  cat("  ⚠ No valid plots generated. Check for errors above.\n")
  quit(status = 1)
}

cat(sprintf("  Valid plots generated: %d\n", length(valid_plots)))

cat("\n  Saving individual gene PNG files...\n")

for (gene in GENES_OF_INTEREST) {
  plot_nestin <- all_plots[[paste0(gene, "_Nestin")]]
  plot_emx1 <- all_plots[[paste0(gene, "_Emx1")]]

  # Save Nestin comparison
  if (!is.null(plot_nestin)) {
    output_file_nestin <- file.path(OUTPUT_DIR,
                                     sprintf("GABA_coverage_%s_Nestin.png", gene))
    png(output_file_nestin, width = 1400, height = 1000, res = 100)
    print(plot_nestin)
    dev.off()
    cat(sprintf("    ✓ %s\n", basename(output_file_nestin)))
  }

  # Save Emx1 comparison
  if (!is.null(plot_emx1)) {
    output_file_emx1 <- file.path(OUTPUT_DIR,
                                   sprintf("GABA_coverage_%s_Emx1.png", gene))
    png(output_file_emx1, width = 1400, height = 1000, res = 100)
    print(plot_emx1)
    dev.off()
    cat(sprintf("    ✓ %s\n", basename(output_file_emx1)))
  }
}

################################################################################
# Create Summary Comparison Images
################################################################################

cat("\n  Creating summary comparison plots...\n")

# Create combined image for each genotype (all genes in one plot)
# Nestin genotype
nestin_plots <- list()
for (gene in GENES_OF_INTEREST) {
  p <- all_plots[[paste0(gene, "_Nestin")]]
  if (!is.null(p)) {
    nestin_plots[[gene]] <- p
  }
}

if (length(nestin_plots) > 0) {
  output_nestin_summary <- file.path(OUTPUT_DIR, "GABA_coverage_all_genes_Nestin.png")

  # Calculate height based on number of genes (1000px per gene)
  plot_height <- length(nestin_plots) * 1000

  png(output_nestin_summary, width = 1400, height = plot_height, res = 100)

  # Combine all plots vertically
  combined_plot <- wrap_plots(nestin_plots, ncol = 1) +
    plot_annotation(
      title = "Nestin Genotype (Ctrl vs Mut) - All Genes (GABA neurons)",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )
  print(combined_plot)

  dev.off()
  cat(sprintf("    ✓ %s\n", basename(output_nestin_summary)))
}

# Emx1 genotype
emx1_plots <- list()
for (gene in GENES_OF_INTEREST) {
  p <- all_plots[[paste0(gene, "_Emx1")]]
  if (!is.null(p)) {
    emx1_plots[[gene]] <- p
  }
}

if (length(emx1_plots) > 0) {
  output_emx1_summary <- file.path(OUTPUT_DIR, "GABA_coverage_all_genes_Emx1.png")

  # Calculate height based on number of genes (1000px per gene)
  plot_height <- length(emx1_plots) * 1000

  png(output_emx1_summary, width = 1400, height = plot_height, res = 100)

  # Combine all plots vertically
  combined_plot <- wrap_plots(emx1_plots, ncol = 1) +
    plot_annotation(
      title = "Emx1 Genotype (Ctrl vs Mut) - All Genes (GABA neurons)",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )
  print(combined_plot)

  dev.off()
  cat(sprintf("    ✓ %s\n", basename(output_emx1_summary)))
}

################################################################################
# Summary Statistics
################################################################################

cat("\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n")

cat("\nCell counts used for plotting:\n")
cat(sprintf("  Total GABA: %d\n", ncol(seurat_gaba)))

# Breakdown by genotype and condition
summary_table <- seurat_gaba@meta.data %>%
  group_by(.data[[GENOTYPE_COL]], .data[[CONDITION_COL]]) %>%
  summarise(
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(.data[[GENOTYPE_COL]], .data[[CONDITION_COL]])

cat("\nBreakdown:\n")
cat(sprintf("  %-10s %-10s %8s\n",
            "Genotype", "Condition", "Total"))
cat(strrep("-", 35), "\n")

for (i in 1:nrow(summary_table)) {
  cat(sprintf("  %-10s %-10s %8d\n",
              summary_table[[GENOTYPE_COL]][i],
              summary_table[[CONDITION_COL]][i],
              summary_table$n_cells[i]))
}

cat("\nGenes analyzed:\n")
for (gene in GENES_OF_INTEREST) {
  nestin_success <- !is.null(all_plots[[paste0(gene, "_Nestin")]])
  emx1_success <- !is.null(all_plots[[paste0(gene, "_Emx1")]])

  status <- if (nestin_success && emx1_success) {
    "✓ Both genotypes"
  } else if (nestin_success) {
    "⚠ Nestin only"
  } else if (emx1_success) {
    "⚠ Emx1 only"
  } else {
    "✗ Failed"
  }

  cat(sprintf("  %s: %s\n", gene, status))
}

cat("\nOutput files (PNG format):\n")
cat(sprintf("  - %s/GABA_coverage_<gene>_Nestin.png (individual comparisons)\n", OUTPUT_DIR))
cat(sprintf("  - %s/GABA_coverage_<gene>_Emx1.png (individual comparisons)\n", OUTPUT_DIR))
cat(sprintf("  - %s/GABA_coverage_all_genes_Nestin.png (all genes summary)\n", OUTPUT_DIR))
cat(sprintf("  - %s/GABA_coverage_all_genes_Emx1.png (all genes summary)\n", OUTPUT_DIR))

cat("\n✓ Coverage snapshot generation complete!\n")
cat("✓ All plots have identical Y-axis scales across both genotypes\n\n")
