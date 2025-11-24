#!/usr/bin/env Rscript

################################################################################
# GC Coverage Snapshots - Combined Controls Version
################################################################################
#
# This script creates publication-quality chromatin accessibility coverage
# plots for GC (Immature + Mature), showing:
#   - Genes: Read from genes_gc.txt file
#   - Comparisons: Combined Controls (Nestin-Ctrl + Emx1-Ctrl) vs each mutant
#     * Combined Controls vs Nestin-Mut
#     * Combined Controls vs Emx1-Mut
#   - Region: Gene body + 100kb upstream/downstream
#   - Shows: ATAC coverage tracks, peaks, gene annotation, and RNA expression
#   - Uses robust Signac CoveragePlot for genome browser visualization
#
# INPUT FILES:
#   - Integrated Seurat object with processed data:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/integrated_seurat_processed.rds
#   - Gene list file (GC-specific genes):
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/genes_gc.txt
#   - Fragment files (linked in Seurat object ATAC assay):
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/translated_fragments/*_fragments_translated.tsv.gz
#
# OUTPUT FILES:
#   - Main combined PDF with all genes:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/plots/gc_coverage_snapshots_combined_ctrl/GC_coverage_snapshots_combined_controls.pdf
#   - Individual gene PDFs:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/plots/gc_coverage_snapshots_combined_ctrl/GC_coverage_<gene>_combined_ctrl.pdf
#   - Summary PDF with all genes by comparison:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/plots/gc_coverage_snapshots_combined_ctrl/GC_coverage_summary_all_genes_combined_ctrl.pdf
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
cat("  GC COVERAGE SNAPSHOTS - COMBINED CONTROLS VERSION\n")
cat("================================================================================\n\n")

################################################################################
# Configuration
################################################################################

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
INPUT_DIR <- file.path(BASE_DIR, "signac_results")
OUTPUT_DIR <- file.path(INPUT_DIR, "plots", "gc_coverage_snapshots_combined_ctrl")

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Read genes of interest from file
GENES_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/genes_gc.txt"
GENES_OF_INTEREST <- readLines(GENES_FILE)
GENES_OF_INTEREST <- GENES_OF_INTEREST[nchar(GENES_OF_INTEREST) > 0]  # Remove empty lines
GENES_OF_INTEREST <- unique(GENES_OF_INTEREST)  # Remove duplicates

# Cell types to include (Granule Cells)
CELL_TYPES <- c("Immature GC", "Mature GC")

# Genomic window around gene
EXTEND_UPSTREAM <- 100000   # 100kb upstream
EXTEND_DOWNSTREAM <- 100000  # 100kb downstream

# Metadata columns
CELLTYPE_COL <- "cell_type_L2_new"
GENOTYPE_COL <- "genotype"
CONDITION_COL <- "condition"

cat("Configuration:\n")
cat(sprintf("  Gene list file: %s\n", GENES_FILE))
cat(sprintf("  Genes loaded: %d\n", length(GENES_OF_INTEREST)))
cat(sprintf("    %s\n", paste(GENES_OF_INTEREST, collapse = ", ")))
cat(sprintf("  Cell types: %s\n", paste(CELL_TYPES, collapse = " + ")))
cat(sprintf("  Genomic window: Gene ± %d kb\n", EXTEND_UPSTREAM / 1000))
cat(sprintf("  Comparison strategy: Combined Controls (Nestin+Emx1) vs each mutant\n"))
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
# Filter to GC
################################################################################

cat("\n=== Filtering to GC ===\n")

# Get cells that are Immature GC or Mature GC
gc_cells <- seurat_obj@meta.data %>%
  filter(.data[[CELLTYPE_COL]] %in% CELL_TYPES) %>%
  rownames()

cat(sprintf("  Total GC cells: %d\n", length(gc_cells)))

# Show breakdown by cell type
gc_breakdown <- seurat_obj@meta.data %>%
  filter(.data[[CELLTYPE_COL]] %in% CELL_TYPES) %>%
  group_by(.data[[CELLTYPE_COL]]) %>%
  summarise(n = n(), .groups = "drop")

for (i in 1:nrow(gc_breakdown)) {
  cat(sprintf("    %s: %d cells\n",
              gc_breakdown[[CELLTYPE_COL]][i],
              gc_breakdown$n[i]))
}

# Show breakdown by genotype and condition
cat("\n  Breakdown by genotype and condition:\n")
gc_genotype_breakdown <- seurat_obj@meta.data %>%
  filter(.data[[CELLTYPE_COL]] %in% CELL_TYPES) %>%
  group_by(.data[[GENOTYPE_COL]], .data[[CONDITION_COL]]) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(.data[[GENOTYPE_COL]], .data[[CONDITION_COL]])

for (i in 1:nrow(gc_genotype_breakdown)) {
  cat(sprintf("    %s %s: %d cells\n",
              gc_genotype_breakdown[[GENOTYPE_COL]][i],
              gc_genotype_breakdown[[CONDITION_COL]][i],
              gc_genotype_breakdown$n[i]))
}

# Subset Seurat object
seurat_gc <- subset(seurat_obj, cells = gc_cells)
cat(sprintf("\n  ✓ Subsetted to %d GC cells\n", ncol(seurat_gc)))

# Add combined cell type label for plotting
seurat_gc@meta.data$combined_celltype <- "GC (Immature + Mature)"

################################################################################
# Create Combined Grouping Variable
################################################################################

cat("\n=== Creating combined control grouping ===\n")

# Create a new grouping variable:
# - "Combined_Ctrl" for both Nestin-Ctrl and Emx1-Ctrl
# - "Nestin_Mut" for Nestin-Mut
# - "Emx1_Mut" for Emx1-Mut

seurat_gc@meta.data$combined_group <- ifelse(
  seurat_gc@meta.data[[CONDITION_COL]] == "Ctrl",
  "Combined_Ctrl",  # Combine all controls (Nestin + Emx1)
  paste(seurat_gc@meta.data[[GENOTYPE_COL]], "Mut", sep = "_")  # Keep mutants separate
)

# Show breakdown by combined group
combined_breakdown <- seurat_gc@meta.data %>%
  group_by(combined_group) %>%
  summarise(
    n = n(),
    nestin = sum(.data[[GENOTYPE_COL]] == "Nestin"),
    emx1 = sum(.data[[GENOTYPE_COL]] == "Emx1"),
    immature = sum(.data[[CELLTYPE_COL]] == "Immature GC"),
    mature = sum(.data[[CELLTYPE_COL]] == "Mature GC"),
    .groups = "drop"
  ) %>%
  arrange(combined_group)

cat("\n  Combined grouping breakdown:\n")
for (i in 1:nrow(combined_breakdown)) {
  cat(sprintf("    %s: %d cells (Nestin: %d, Emx1: %d, Immature: %d, Mature: %d)\n",
              combined_breakdown$combined_group[i],
              combined_breakdown$n[i],
              combined_breakdown$nestin[i],
              combined_breakdown$emx1[i],
              combined_breakdown$immature[i],
              combined_breakdown$mature[i]))
}

################################################################################
# Function: Create Coverage Plot for Gene (Combined Controls)
################################################################################

create_gene_coverage_plot_combined <- function(seurat_obj, gene, mutant_genotype,
                                               extend_upstream = 100000,
                                               extend_downstream = 100000,
                                               max_cutoff = NULL) {

  cat(sprintf("\n  Creating coverage plot: %s (Combined Controls vs %s Mut)\n", gene, mutant_genotype))

  # Filter to cells we want to compare:
  # - All controls (Combined_Ctrl)
  # - Specific mutant (e.g., Nestin_Mut or Emx1_Mut)
  mutant_label <- paste(mutant_genotype, "Mut", sep = "_")

  cells_to_include <- seurat_obj@meta.data %>%
    filter(combined_group %in% c("Combined_Ctrl", mutant_label)) %>%
    rownames()

  if (length(cells_to_include) == 0) {
    cat(sprintf("    ⚠ No cells found for comparison\n"))
    return(NULL)
  }

  # Subset to selected cells
  seurat_subset <- subset(seurat_obj, cells = cells_to_include)
  cat(sprintf("    Cells in comparison: %d\n", ncol(seurat_subset)))

  # Validate that fragment files are still linked after subsetting
  if ("ATAC" %in% names(seurat_subset@assays)) {
    fragments_subset <- Fragments(seurat_subset[["ATAC"]])
    if (length(fragments_subset) == 0) {
      cat("      ⚠ WARNING: Fragment files lost during subsetting. Re-linking...\n")
      # Copy fragment files from original object
      fragments_orig <- Fragments(seurat_obj[["ATAC"]])
      Fragments(seurat_subset[["ATAC"]]) <- fragments_orig
    }
  }

  # Show breakdown by group
  group_counts <- table(seurat_subset@meta.data$combined_group)
  for (grp in names(group_counts)) {
    cat(sprintf("      %s: %d cells\n", grp, group_counts[grp]))
  }

  # Set identity to combined group for plotting
  Idents(seurat_subset) <- "combined_group"

  # Get the identities to plot (should be "Combined_Ctrl" and mutant_label)
  available_idents <- c("Combined_Ctrl", mutant_label)
  available_idents <- intersect(available_idents, unique(as.character(Idents(seurat_subset))))

  cat(sprintf("    Identities to plot: %s\n", paste(available_idents, collapse = ", ")))

  # CRITICAL: Set default assay to ATAC for coverage visualization
  DefaultAssay(seurat_subset) <- "ATAC"

  # Create improved coverage plot with explicit coverage tracks
  p <- tryCatch({
    # Use CoveragePlot with explicit parameters for pile-up visualization
    coverage_plot <- CoveragePlot(
      object = seurat_subset,
      region = gene,
      features = gene,
      expression.assay = "RNA",
      assay = "ATAC",  # Explicitly set ATAC assay for coverage
      extend.upstream = extend_upstream,
      extend.downstream = extend_downstream,
      idents = available_idents,  # Combined controls vs specific mutant
      annotation = TRUE,
      peaks = TRUE,
      tile = FALSE,
      links = FALSE,
      ncol = 1,
      window = 100,  # Window size for coverage calculation (100bp bins)
      group.by = "combined_group",  # Explicitly state the grouping variable
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
        title = sprintf("%s - Combined Controls vs %s Mut (GC: Immature + Mature)", gene, mutant_genotype),
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

  # Calculate max coverage across ALL groups (Combined_Ctrl, Nestin_Mut, Emx1_Mut)
  # This ensures identical Y-axis scaling for both comparison plots

  tryCatch({
    # Set default assay
    DefaultAssay(seurat_obj) <- "ATAC"

    # Set identity to combined_group
    Idents(seurat_obj) <- "combined_group"

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
      group.by = "combined_group",
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
    cat(sprintf("      This value will be used for both comparisons\n"))
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
    seurat_gc,
    gene,
    extend_upstream = EXTEND_UPSTREAM,
    extend_downstream = EXTEND_DOWNSTREAM
  )

  # Create plots comparing combined controls against each mutant
  # Both plots will use the same max_cutoff for consistent scaling
  plot_nestin_mut <- create_gene_coverage_plot_combined(
    seurat_gc,
    gene,
    mutant_genotype = "Nestin",
    extend_upstream = EXTEND_UPSTREAM,
    extend_downstream = EXTEND_DOWNSTREAM,
    max_cutoff = max_coverage
  )

  plot_emx1_mut <- create_gene_coverage_plot_combined(
    seurat_gc,
    gene,
    mutant_genotype = "Emx1",
    extend_upstream = EXTEND_UPSTREAM,
    extend_downstream = EXTEND_DOWNSTREAM,
    max_cutoff = max_coverage
  )

  # Store plots
  all_plots[[paste0(gene, "_vs_Nestin_Mut")]] <- plot_nestin_mut
  all_plots[[paste0(gene, "_vs_Emx1_Mut")]] <- plot_emx1_mut
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
  plot_nestin <- all_plots[[paste0(gene, "_vs_Nestin_Mut")]]
  plot_emx1 <- all_plots[[paste0(gene, "_vs_Emx1_Mut")]]

  # Save Nestin comparison
  if (!is.null(plot_nestin)) {
    output_file_nestin <- file.path(OUTPUT_DIR,
                                     sprintf("GC_coverage_%s_vs_Nestin_Mut.png", gene))
    png(output_file_nestin, width = 1400, height = 1000, res = 100)
    print(plot_nestin)
    dev.off()
    cat(sprintf("    ✓ %s\n", basename(output_file_nestin)))
  }

  # Save Emx1 comparison
  if (!is.null(plot_emx1)) {
    output_file_emx1 <- file.path(OUTPUT_DIR,
                                   sprintf("GC_coverage_%s_vs_Emx1_Mut.png", gene))
    png(output_file_emx1, width = 1400, height = 1000, res = 100)
    print(plot_emx1)
    dev.off()
    cat(sprintf("    ✓ %s\n", basename(output_file_emx1)))
  }
}

################################################################################
# Create Summary Comparison Images
################################################################################

# cat("\n  Creating summary comparison plots...\n")

# # Create combined image for each comparison (all genes in one plot)
# # Nestin comparison
# nestin_plots <- list()
# for (gene in GENES_OF_INTEREST) {
#   p <- all_plots[[paste0(gene, "_vs_Nestin_Mut")]]
#   if (!is.null(p)) {
#     nestin_plots[[gene]] <- p
#   }
# }

# if (length(nestin_plots) > 0) {
#   output_nestin_summary <- file.path(OUTPUT_DIR, "GC_coverage_all_genes_vs_Nestin_Mut.png")

#   # Calculate height based on number of genes (1000px per gene)
#   plot_height <- length(nestin_plots) * 1000

#   png(output_nestin_summary, width = 1400, height = plot_height, res = 100)

#   # Combine all plots vertically
#   combined_plot <- wrap_plots(nestin_plots, ncol = 1) +
#     plot_annotation(
#       title = "Combined Controls vs Nestin Mut - All Genes (GC: Immature + Mature)",
#       theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
#     )
#   print(combined_plot)

#   dev.off()
#   cat(sprintf("    ✓ %s\n", basename(output_nestin_summary)))
# }

# # Emx1 comparison
# emx1_plots <- list()
# for (gene in GENES_OF_INTEREST) {
#   p <- all_plots[[paste0(gene, "_vs_Emx1_Mut")]]
#   if (!is.null(p)) {
#     emx1_plots[[gene]] <- p
#   }
# }

# if (length(emx1_plots) > 0) {
#   output_emx1_summary <- file.path(OUTPUT_DIR, "GC_coverage_all_genes_vs_Emx1_Mut.png")

#   # Calculate height based on number of genes (1000px per gene)
#   plot_height <- length(emx1_plots) * 1000

#   png(output_emx1_summary, width = 1400, height = plot_height, res = 100)

#   # Combine all plots vertically
#   combined_plot <- wrap_plots(emx1_plots, ncol = 1) +
#     plot_annotation(
#       title = "Combined Controls vs Emx1 Mut - All Genes (GC: Immature + Mature)",
#       theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
#     )
#   print(combined_plot)

#   dev.off()
#   cat(sprintf("    ✓ %s\n", basename(output_emx1_summary)))
# }

################################################################################
# Summary Statistics
################################################################################

cat("\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n")

cat("\nCell counts used for plotting:\n")
cat(sprintf("  Total GC: %d\n", ncol(seurat_gc)))

# Breakdown by combined group
summary_table <- seurat_gc@meta.data %>%
  group_by(combined_group) %>%
  summarise(
    n_cells = n(),
    nestin = sum(.data[[GENOTYPE_COL]] == "Nestin"),
    emx1 = sum(.data[[GENOTYPE_COL]] == "Emx1"),
    immature = sum(.data[[CELLTYPE_COL]] == "Immature GC"),
    mature = sum(.data[[CELLTYPE_COL]] == "Mature GC"),
    .groups = "drop"
  ) %>%
  arrange(combined_group)

cat("\nBreakdown by combined group:\n")
cat(sprintf("  %-18s %8s %10s %10s %10s %10s\n",
            "Group", "Total", "Nestin", "Emx1", "Immature", "Mature"))
cat(strrep("-", 75), "\n")

for (i in 1:nrow(summary_table)) {
  cat(sprintf("  %-18s %8d %10d %10d %10d %10d\n",
              summary_table$combined_group[i],
              summary_table$n_cells[i],
              summary_table$nestin[i],
              summary_table$emx1[i],
              summary_table$immature[i],
              summary_table$mature[i]))
}

cat("\nGenes analyzed:\n")
for (gene in GENES_OF_INTEREST) {
  nestin_success <- !is.null(all_plots[[paste0(gene, "_vs_Nestin_Mut")]])
  emx1_success <- !is.null(all_plots[[paste0(gene, "_vs_Emx1_Mut")]])

  status <- if (nestin_success && emx1_success) {
    "✓ Both comparisons"
  } else if (nestin_success) {
    "⚠ Nestin comparison only"
  } else if (emx1_success) {
    "⚠ Emx1 comparison only"
  } else {
    "✗ Failed"
  }

  cat(sprintf("  %s: %s\n", gene, status))
}

cat("\nOutput files (PNG format):\n")
cat(sprintf("  - %s/GC_coverage_<gene>_vs_Nestin_Mut.png (individual comparisons)\n", OUTPUT_DIR))
cat(sprintf("  - %s/GC_coverage_<gene>_vs_Emx1_Mut.png (individual comparisons)\n", OUTPUT_DIR))
cat(sprintf("  - %s/GC_coverage_all_genes_vs_Nestin_Mut.png (all genes summary)\n", OUTPUT_DIR))
cat(sprintf("  - %s/GC_coverage_all_genes_vs_Emx1_Mut.png (all genes summary)\n", OUTPUT_DIR))

cat("\n✓ Coverage snapshot generation complete!\n")
cat("✓ All plots have identical Y-axis scales for Combined Controls\n\n")
