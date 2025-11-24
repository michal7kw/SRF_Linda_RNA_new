#!/usr/bin/env Rscript

################################################################################
# Peak-Gene Linkage Analysis for Granule Cells (UPDATED VERSION - Using Pre-computed Links)
#
# This script analyzes peak-gene linkages for genes in genes_gc.txt,
# using PRE-COMPUTED LinkPeaks results from array job outputs.
#
# KEY CHANGE: Loads pre-computed LinkPeaks from array jobs instead of manual cor.test()
#
# USES L2 PIPELINE: This script analyzes L2 cell types (Immature GC, Mature GC).
# Granule cell subtypes only exist at L2 annotation level.
#
# PREREQUISITES: Run peak-gene linkage array jobs first:
#   sbatch run_peak_gene_linkage_high_priority.sh  OR
#   sbatch run_peak_gene_linkage_all.sh
#
# Note: For broad excitatory neuron analysis (L1 level), use
#       create_coverage_plots_from_DEGs_enhanced_L1.R with "Excitatory" cell type.
#
# Input files:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/integrated_seurat_processed.rds
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/genes_gc.txt
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/celltype_results/peak_gene_links/Immature_GC_*_peak_gene_links.csv
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/celltype_results/peak_gene_links/Mature_GC_*_peak_gene_links.csv
#
# Output:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/peak_gene_linkage_analysis_UPDATED/GC/peak_gene_links/{gene}_links.csv
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/peak_gene_linkage_analysis_UPDATED/GC/coverage_plots/{gene}_coverage.pdf
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/peak_gene_linkage_analysis_UPDATED/GC/analysis_summary.csv
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
  library(ggforce)
})

# Configuration
RESULTS_DIR <- "signac_results"  # L2 pipeline - GC subtypes only exist at L2 level
SEURAT_FILE <- file.path(RESULTS_DIR, "integrated_seurat_processed.rds")
GENE_LIST_FILE <- "../genes_gc.txt"
OUTPUT_DIR <- file.path(RESULTS_DIR, "peak_gene_linkage_analysis_UPDATED", "GC")

# L2 peak-gene links directory (where GC-specific links are stored)
STEP3_LINKS_DIR <- file.path(RESULTS_DIR, "celltype_results", "peak_gene_links")

# Cell type columns
CELLTYPE_COL <- "cell_type_L2_new"
GC_CELL_TYPES <- c("Immature GC", "Mature GC", "GC")  # Try all possible GC names

CONDITION_COL <- "condition"
GENOTYPE_COL <- "genotype"

# Plotting parameters
EXTEND_UPSTREAM <- 50000
EXTEND_DOWNSTREAM <- 50000

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "peak_gene_links"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "coverage_plots"), recursive = TRUE, showWarnings = FALSE)

cat("================================================================================\n")
cat("GRANULE CELLS: Combined Controls (Using Pre-computed Links from Step 3)\n")
cat("================================================================================\n\n")

################################################################################
# Load Data
################################################################################

cat("Loading data...\n")
seurat_obj <- readRDS(SEURAT_FILE)
genes_of_interest <- read.table(GENE_LIST_FILE, header = FALSE)$V1 %>% unique()

# Filter mt-genes
mt_genes <- grepl("^mt-", genes_of_interest, ignore.case = TRUE)
if (any(mt_genes)) {
  cat(sprintf("  Filtering %d mt-genes\n", sum(mt_genes)))
  genes_of_interest <- genes_of_interest[!mt_genes]
}

cat(sprintf("  ✓ %d genes to analyze\n", length(genes_of_interest)))

################################################################################
# Load Pre-computed Links
################################################################################

cat("\nLoading pre-computed links from Step 3...\n")

# Determine which results directory to use
if (dir.exists(STEP3_LINKS_DIR_L2)) {
  links_dir <- STEP3_LINKS_DIR_L2
  cat("  Using L2 array job results (signac_results/)\n")
} else if (dir.exists(STEP3_LINKS_DIR_L1)) {
  links_dir <- STEP3_LINKS_DIR_L1
  cat("  Using L1 array job results (signac_results_L1/)\n")
} else {
  cat("\n✗ ERROR: Peak-gene linkage results directory not found!\n\n")
  cat("You must run peak-gene linkage array jobs first:\n")
  cat("  Option 1 (High Priority): sbatch run_peak_gene_linkage_high_priority_L1.sh\n")
  cat("  Option 2 (All Cell Types): sbatch run_peak_gene_linkage_all_L1.sh\n\n")
  stop("Peak-gene linkage results not found. Run array jobs first.")
}

# Find all GC-related link files
all_files <- list.files(links_dir, pattern = "_peak_gene_links\\.csv$", full.names = TRUE)

# Filter to GC-related files (Immature GC, Mature GC, or just GC)
gc_files <- all_files[grepl("GC", basename(all_files), ignore.case = FALSE)]

if (length(gc_files) == 0) {
  cat("\n✗ ERROR: No GC-related LinkPeaks files found!\n")
  cat(sprintf("  Searched in: %s\n", links_dir))
  cat("  Available files:\n")
  for (f in basename(all_files)[1:min(10, length(all_files))]) {
    cat(sprintf("    - %s\n", f))
  }
  cat("\nThis script expects L2 cell type files (Immature GC, Mature GC).\n")
  cat("Current array jobs generate L1 cell types only.\n\n")
  cat("ALTERNATIVES:\n")
  cat("  1. For L1 analysis, use: create_coverage_plots_from_DEGs_enhanced.R\n")
  cat("     with Excitatory neurons (likely includes granule cells)\n")
  cat("  2. Run L2-specific peak-gene linkage (not currently implemented)\n\n")
  stop("No GC LinkPeaks files found")
}

cat("  Found GC link files:\n")
for (f in basename(gc_files)) {
  cat(sprintf("    - %s\n", f))
}

# Load and combine all GC link files
all_links_list <- list()
for (file in gc_files) {
  file_base <- tools::file_path_sans_ext(basename(file))
  df <- read.csv(file)

  # Ensure genotype column exists with capitalized names
  # Array job outputs already include genotype column with proper capitalization
  if (!"genotype" %in% colnames(df)) {
    # If missing, extract from filename (e.g., "Immature_GC_Nestin_peak_gene_links" -> "Nestin")
    if (grepl("Nestin", file_base)) {
      df$genotype <- "Nestin"
    } else if (grepl("Emx1", file_base)) {
      df$genotype <- "Emx1"
    } else {
      df$genotype <- "unknown"
    }
  }

  all_links_list[[file_base]] <- df
  cat(sprintf("    Loaded %d links from %s\n", nrow(df), basename(file)))
}

all_links <- do.call(rbind, all_links_list)
cat(sprintf("\n✓ Total links loaded: %d\n", nrow(all_links)))

# Filter to genes of interest
links_filtered <- all_links %>% filter(gene %in% genes_of_interest)

cat(sprintf("✓ Filtered to genes of interest: %d links for %d genes\n",
            nrow(links_filtered),
            length(unique(links_filtered$gene))))

if (nrow(links_filtered) == 0) {
  stop("No links found for genes of interest")
}

################################################################################
# Subset to Granule Cells
################################################################################

cat("\nSubsetting to Granule Cells...\n")

# Try to find GC cells using all possible GC cell type names
gc_cells <- c()
for (ct_name in GC_CELL_TYPES) {
  cells <- seurat_obj@meta.data %>%
    filter(.data[[CELLTYPE_COL]] == ct_name) %>%
    rownames()
  if (length(cells) > 0) {
    cat(sprintf("  Found %d cells with cell type '%s'\n", length(cells), ct_name))
    gc_cells <- c(gc_cells, cells)
  }
}

gc_cells <- unique(gc_cells)

if (length(gc_cells) == 0) {
  stop(sprintf("Error: No GC cells found. Tried cell types: %s\n",
               paste(GC_CELL_TYPES, collapse = ", ")))
}

seurat_gc <- subset(seurat_obj, cells = gc_cells)
cat(sprintf("✓ Selected %d Granule Cells\n", ncol(seurat_gc)))

# Create combined grouping
seurat_gc@meta.data$combined_group <- ifelse(
  seurat_gc@meta.data[[CONDITION_COL]] == "Ctrl",
  "Combined_Ctrl",
  paste0(seurat_gc@meta.data[[GENOTYPE_COL]], "_Mut")
)

cat("\nGroup distribution:\n")
print(table(seurat_gc@meta.data$combined_group))

Idents(seurat_gc) <- seurat_gc@meta.data$combined_group

################################################################################
# Create Visualizations
################################################################################

cat("\n================================================================================\n")
cat("CREATING VISUALIZATIONS\n")
cat("================================================================================\n\n")

summary_stats <- data.frame(
  gene = character(),
  status = character(),
  n_links = integer(),
  stringsAsFactors = FALSE
)

genes_to_process <- unique(links_filtered$gene)

for (gene in genes_to_process) {
  cat(sprintf("--- %s ---\n", gene))

  # Validate gene
  gene_in_atac <- tryCatch({
    coords <- LookupGeneCoords(seurat_gc, gene = gene, assay = "ATAC")
    !is.null(coords) && length(coords) > 0
  }, error = function(e) FALSE)

  if (!gene_in_atac) {
    cat("  ⚠ Gene not in ATAC annotation\n")
    summary_stats <- rbind(summary_stats, data.frame(
      gene = gene, status = "missing_annotation", n_links = 0
    ))
    next
  }

  # Get links
  gene_links <- links_filtered %>% filter(gene == !!gene)
  cat(sprintf("  %d links\n", nrow(gene_links)))

  # Save CSV
  write.csv(gene_links,
            file.path(OUTPUT_DIR, "peak_gene_links", sprintf("%s_links.csv", gene)),
            row.names = FALSE)

  # Create coverage plot
  plot_success <- tryCatch({
    DefaultAssay(seurat_gc) <- "ATAC"

    # Parse peaks
    if ("seqnames" %in% colnames(gene_links)) {
      peak_gr <- GRanges(
        seqnames = gene_links$seqnames,
        ranges = IRanges(start = as.numeric(gene_links$start),
                        end = as.numeric(gene_links$end)),
        score = gene_links$score
      )
    } else if ("peak" %in% colnames(gene_links)) {
      coords <- do.call(rbind, strsplit(as.character(gene_links$peak), "-"))
      peak_gr <- GRanges(
        seqnames = coords[, 1],
        ranges = IRanges(start = as.numeric(coords[, 2]),
                        end = as.numeric(coords[, 3])),
        score = gene_links$score
      )
    } else {
      cat("  ⚠ Cannot parse peak coordinates\n")
      return(FALSE)
    }

    # Create links
    gene_coords <- LookupGeneCoords(seurat_gc, gene = gene, assay = "ATAC")
    gene_start <- start(gene_coords)

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
    Links(seurat_gc[["ATAC"]]) <- all_link_gr

    # Plot
    p_cov <- CoveragePlot(
      object = seurat_gc,
      region = gene,
      features = gene,
      expression.assay = "RNA",
      assay = "ATAC",
      group.by = "combined_group",
      extend.upstream = EXTEND_UPSTREAM,
      extend.downstream = EXTEND_DOWNSTREAM,
      annotation = TRUE,
      peaks = TRUE,
      links = TRUE,
      window = 100
    )

    plot_file <- file.path(OUTPUT_DIR, "coverage_plots",
                          sprintf("%s_coverage.pdf", gene))
    ggsave(plot_file, plot = p_cov, width = 14, height = 10, device = "pdf")
    cat(sprintf("  ✓ Saved: %s\n", basename(plot_file)))

    TRUE
  }, error = function(e) {
    cat(sprintf("  ✗ Plot error: %s\n", e$message))
    FALSE
  })

  summary_stats <- rbind(summary_stats, data.frame(
    gene = gene,
    status = ifelse(plot_success, "success", "failed"),
    n_links = nrow(gene_links)
  ))
}

################################################################################
# Save Summary
################################################################################

cat("\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n\n")

summary_file <- file.path(OUTPUT_DIR, "analysis_summary.csv")
write.csv(summary_stats, summary_file, row.names = FALSE)

cat(sprintf("Total genes: %d\n", nrow(summary_stats)))
cat(sprintf("Successful: %d\n", sum(summary_stats$status == "success")))
cat(sprintf("Failed: %d\n", sum(summary_stats$status == "failed")))

cat(sprintf("\n✓ Results saved to: %s\n", OUTPUT_DIR))
