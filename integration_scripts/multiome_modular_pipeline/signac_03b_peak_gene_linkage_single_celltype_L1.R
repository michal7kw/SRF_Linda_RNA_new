#!/usr/bin/env Rscript

################################################################################
# Signac Multiome Analysis - Step 3b: Peak-Gene Linkage for Single Cell Type (L1)
################################################################################
#
# This script performs peak-gene linkage for a SINGLE cell type × genotype
# combination, allowing for longer runtime with more computational resources.
#
# Usage:
#   Rscript signac_03b_peak_gene_linkage_single_celltype_L1.R <CELLTYPE> <GENOTYPE>
#
# Example:
#   Rscript signac_03b_peak_gene_linkage_single_celltype_L1.R Excitatory Nestin
#
# This addresses the peak count limitation by dedicating resources to one
# cell type at a time, removing the MAX_PEAKS_FOR_LINKAGE threshold.
#
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(dplyr)
})

cat("================================================================================\n")
cat("  PEAK-GENE LINKAGE - SINGLE CELL TYPE (L1)\n")
cat("================================================================================\n\n")

################################################################################
# Parse Command Line Arguments
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("ERROR: Incorrect number of arguments\n\n")
  cat("Usage: Rscript signac_03b_peak_gene_linkage_single_celltype_L1.R <CELLTYPE> <GENOTYPE>\n")
  cat("\nExample:\n")
  cat("  Rscript signac_03b_peak_gene_linkage_single_celltype_L1.R Excitatory Nestin\n\n")
  cat("Available cell types: Unknown, Excitatory, Immune, Oligo, Astrocytes, GABA, Vascular, Progenitors, Ependymal\n")
  cat("Available genotypes: Emx1, Nestin\n\n")
  quit(status = 1)
}

CELLTYPE_TARGET <- args[1]
GENOTYPE_TARGET <- args[2]

cat(sprintf("Target cell type: %s\n", CELLTYPE_TARGET))
cat(sprintf("Target genotype: %s\n", GENOTYPE_TARGET))
cat("\n")

################################################################################
# Configuration
################################################################################

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
INPUT_DIR <- file.path(BASE_DIR, "signac_results_L1")
RESULTS_DIR <- file.path(INPUT_DIR, "celltype_results")

# Create output directory
dir.create(file.path(RESULTS_DIR, "peak_gene_links"), recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
CELLTYPE_COL <- "cell_type_L1"
CONDITION_COL <- "condition"
GENOTYPE_COL <- "genotype"
MIN_CELLS_PER_GROUP <- 30

# Peak-gene linkage parameters (OPTIMIZED for comprehensive analysis)
PEAK_GENE_DISTANCE <- 200000         # 100kb upstream/downstream
MIN_CORRELATION <- 0.05               # Minimum correlation threshold
PVALUE_CUTOFF <- 0.05                # P-value cutoff
MIN_CELLS_FOR_LINKAGE <- 10          # Minimum cells expressing both peak and gene

# NO MAX_PEAKS_FOR_LINKAGE - we'll process however many peaks exist

cat("Configuration:\n")
cat(sprintf("  Cell type column: %s\n", CELLTYPE_COL))
cat(sprintf("  Target: %s (%s)\n", CELLTYPE_TARGET, GENOTYPE_TARGET))
cat(sprintf("  Peak-gene distance: %d bp\n", PEAK_GENE_DISTANCE))
cat(sprintf("  Min correlation: %.2f\n", MIN_CORRELATION))
cat(sprintf("  P-value cutoff: %.2f\n", PVALUE_CUTOFF))
cat(sprintf("  No peak count limit (dedicated processing)\n"))
cat("\n")

################################################################################
# Load Data
################################################################################

cat("\n=== Loading processed Seurat object ===\n")
input_file <- file.path(INPUT_DIR, "integrated_seurat_processed.rds")

if (!file.exists(input_file)) {
  stop(sprintf("Input file not found: %s\nPlease run signac_02_process_and_reduce.R first.", input_file))
}

seurat_obj <- readRDS(input_file)
cat(sprintf("  Loaded object: %d cells\n", ncol(seurat_obj)))

# Verify target cell type and genotype exist
if (!CELLTYPE_TARGET %in% seurat_obj@meta.data[[CELLTYPE_COL]]) {
  stop(sprintf("Cell type '%s' not found in metadata.", CELLTYPE_TARGET))
}

if (!GENOTYPE_TARGET %in% seurat_obj@meta.data[[GENOTYPE_COL]]) {
  stop(sprintf("Genotype '%s' not found in metadata.", GENOTYPE_TARGET))
}

################################################################################
# Subset to Target Cell Type × Genotype
################################################################################

cat(sprintf("\n=== Subsetting to %s (%s) ===\n", CELLTYPE_TARGET, GENOTYPE_TARGET))

# Filter to target genotype and cell type
cells_target <- seurat_obj@meta.data[[CELLTYPE_COL]] == CELLTYPE_TARGET &
                seurat_obj@meta.data[[GENOTYPE_COL]] == GENOTYPE_TARGET

n_cells_total <- sum(cells_target)
n_cells_ctrl <- sum(cells_target & seurat_obj@meta.data[[CONDITION_COL]] == "Ctrl")
n_cells_mut <- sum(cells_target & seurat_obj@meta.data[[CONDITION_COL]] == "Mut")

cat(sprintf("  Total cells: %d (Ctrl: %d, Mut: %d)\n", n_cells_total, n_cells_ctrl, n_cells_mut))

# Check minimum cells
if (n_cells_total < MIN_CELLS_PER_GROUP * 2) {
  stop(sprintf("Insufficient cells for analysis (need at least %d total, have %d)",
               MIN_CELLS_PER_GROUP * 2, n_cells_total))
}

# Subset Seurat object
cells_to_keep <- rownames(seurat_obj@meta.data[cells_target, ])
seurat_subset <- subset(seurat_obj, cells = cells_to_keep)

cat(sprintf("  Subset created: %d cells\n", ncol(seurat_subset)))

################################################################################
# Filter Invalid Peaks (Critical for mm10)
################################################################################

cat("\n=== Validating peak coordinates ===\n")

# Set default assay to ATAC
DefaultAssay(seurat_subset) <- "ATAC"

# Get mm10 chromosome sizes
genome <- BSgenome.Mmusculus.UCSC.mm10
genome_info <- seqlengths(genome)

# Get peak coordinates using granges
peaks <- granges(seurat_subset[["ATAC"]])
n_peaks_before <- length(peaks)

# Standard chromosomes (including both chrM and chrMT to handle naming variations)
standard_chroms <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM", "chrMT")

cat(sprintf("  Using standard chromosomes (chr1-19, chrX, chrY, chrM/chrMT)\n"))

# Convert peaks to data frame for filtering
peaks_df <- as.data.frame(peaks)

# Normalize chrMT to chrM for consistency with genome
peaks_df$seqnames <- as.character(peaks_df$seqnames)
peaks_df$seqnames[peaks_df$seqnames == "chrMT"] <- "chrM"

# Filter peaks
valid_peaks <- peaks_df %>%
  filter(
    seqnames %in% c(paste0("chr", 1:19), "chrX", "chrY", "chrM"),  # Standard chromosomes only (after normalization)
    start > 0,                       # Positive coordinates
    end > start,                     # Valid range
    end <= genome_info[as.character(seqnames)]  # Within chromosome bounds
  )

valid_peak_names <- paste(valid_peaks$seqnames, valid_peaks$start, valid_peaks$end, sep = "-")

# Subset ATAC assay to valid peaks
seurat_subset[["ATAC"]] <- subset(seurat_subset[["ATAC"]], features = valid_peak_names)

n_peaks_after <- nrow(seurat_subset[["ATAC"]])

cat(sprintf("  Peaks before filtering: %s\n", format(n_peaks_before, big.mark = ",")))
cat(sprintf("  Peaks after filtering:  %s\n", format(n_peaks_after, big.mark = ",")))
cat(sprintf("  Removed: %s invalid peaks\n", format(n_peaks_before - n_peaks_after, big.mark = ",")))

if (n_peaks_after == 0) {
  cat("\nERROR: No valid peaks remaining after filtering\n")
  quit(status = 1)
}

################################################################################
# Compute GC Content (Important for Normalization)
################################################################################

cat("\n================================================================================\n")
cat("COMPUTING PEAK GC CONTENT\n")
cat("================================================================================\n\n")

cat("Computing GC content for peaks...\n")
cat("  This improves peak-gene linkage by accounting for GC bias\n\n")

tryCatch({
  # Suppress seqlevel warnings from genome mismatch
  suppressWarnings({
    seurat_subset <- RegionStats(seurat_subset, genome = genome)
  })
  cat("  ✓ GC content computed successfully\n")
}, error = function(e) {
  cat(sprintf("  ⚠ Warning: GC content computation failed: %s\n", e$message))
  cat("  Continuing without GC content (may affect linkage quality)\n")
})

################################################################################
# Peak-Gene Linkage
################################################################################

cat("\n================================================================================\n")
cat("PEAK-GENE LINKAGE\n")
cat("================================================================================\n\n")

n_peaks <- nrow(seurat_subset[["ATAC"]])
n_genes <- nrow(seurat_subset[["RNA"]])

cat(sprintf("Computing peak-gene links for %s × %s...\n", CELLTYPE_TARGET, GENOTYPE_TARGET))
cat(sprintf("  Peaks: %s\n", format(n_peaks, big.mark = ",")))
cat(sprintf("  Genes: %s\n", format(n_genes, big.mark = ",")))
cat(sprintf("  Peak-gene distance: %d bp (±%d kb)\n",
            PEAK_GENE_DISTANCE, PEAK_GENE_DISTANCE / 1000))
cat(sprintf("  Min correlation: %.2f\n", MIN_CORRELATION))
cat(sprintf("  This may take 1-3 hours depending on cell count and peak number...\n"))
cat("\n")

start_time <- Sys.time()

tryCatch({
  seurat_subset <- LinkPeaks(
    object = seurat_subset,
    peak.assay = "ATAC",
    expression.assay = "RNA",
    genes.use = NULL,  # All genes (explicit)
    distance = PEAK_GENE_DISTANCE,
    min.cells = MIN_CELLS_FOR_LINKAGE,
    pvalue_cutoff = PVALUE_CUTOFF,
    score_cutoff = MIN_CORRELATION,
    gene.coords = NULL,  # Uses default gene coordinates
    gene.id = FALSE,     # Uses gene symbols instead of IDs
    verbose = TRUE
  )

  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

  cat(sprintf("\n  ✓ LinkPeaks completed in %.1f minutes\n", elapsed_time))

}, error = function(e) {
  cat(sprintf("ERROR in LinkPeaks: %s\n", e$message))
  quit(status = 1)
})

################################################################################
# Extract and Save Links
################################################################################

cat("\n=== Extracting peak-gene links ===\n")

links <- Links(seurat_subset)

if (length(links) == 0) {
  cat("  ⚠ No significant peak-gene links found\n")
  cat(sprintf("  This may indicate:\n"))
  cat(sprintf("    - Weak correlation between accessibility and expression\n"))
  cat(sprintf("    - Insufficient cells for robust correlation\n"))
  cat(sprintf("    - Correlation threshold (%.2f) may be too stringent\n", MIN_CORRELATION))

  # Save empty results file
  analysis_id <- sprintf("%s_%s", CELLTYPE_TARGET, GENOTYPE_TARGET)
  output_file <- file.path(RESULTS_DIR, "peak_gene_links", sprintf("%s_peak_gene_links.csv", analysis_id))

  empty_df <- data.frame(
    celltype = character(),
    genotype = character(),
    gene = character(),
    peak = character(),
    score = numeric(),
    pvalue = numeric(),
    stringsAsFactors = FALSE
  )

  write.csv(empty_df, output_file, row.names = FALSE)
  cat(sprintf("\n  ✓ Empty results saved to: %s\n", output_file))

} else {
  cat(sprintf("  ✓ Found %s peak-gene links\n", format(length(links), big.mark = ",")))

  # Convert to data frame
  links_df <- as.data.frame(links)

  # Add metadata
  links_df$celltype <- CELLTYPE_TARGET
  links_df$genotype <- GENOTYPE_TARGET

  # Save results
  analysis_id <- sprintf("%s_%s", CELLTYPE_TARGET, GENOTYPE_TARGET)
  output_file <- file.path(RESULTS_DIR, "peak_gene_links", sprintf("%s_peak_gene_links.csv", analysis_id))

  write.csv(links_df, output_file, row.names = FALSE)

  cat(sprintf("\n  ✓ Results saved to: %s\n", output_file))

  # Print comprehensive summary statistics
  cat("\n=== Summary Statistics ===\n")
  cat(sprintf("  Total links: %s\n", format(nrow(links_df), big.mark = ",")))
  cat(sprintf("  Unique genes: %s\n", format(n_distinct(links_df$gene), big.mark = ",")))
  cat(sprintf("  Unique peaks: %s\n", format(n_distinct(links_df$peak), big.mark = ",")))

  if ("score" %in% colnames(links_df)) {
    cat(sprintf("\nCorrelation statistics:\n"))
    cat(sprintf("  Mean:            %.3f\n", mean(links_df$score, na.rm = TRUE)))
    cat(sprintf("  Median:          %.3f\n", median(links_df$score, na.rm = TRUE)))
    cat(sprintf("  Range:           %.3f to %.3f\n",
                min(links_df$score, na.rm = TRUE),
                max(links_df$score, na.rm = TRUE)))
  }

  if ("pvalue" %in% colnames(links_df)) {
    cat(sprintf("\nP-value range:     %.2e to %.2e\n",
                min(links_df$pvalue, na.rm = TRUE),
                max(links_df$pvalue, na.rm = TRUE)))
  }

  # Calculate peak-gene distances
  if ("start" %in% colnames(links_df) && "peak" %in% colnames(links_df)) {
    cat(sprintf("\nPeak-gene distance statistics:\n"))

    # Extract peak start positions from peak names (format: chr-start-end)
    peak_starts <- as.numeric(sub(".*-(\\d+)-\\d+$", "\\1", links_df$peak))
    distances <- abs(links_df$start - peak_starts)

    cat(sprintf("  Mean distance:   %.1f kb\n", mean(distances, na.rm = TRUE) / 1000))
    cat(sprintf("  Median distance: %.1f kb\n", median(distances, na.rm = TRUE) / 1000))
    cat(sprintf("  Max distance:    %.1f kb\n", max(distances, na.rm = TRUE) / 1000))

    # Distance distribution
    within_10kb <- sum(distances <= 10000, na.rm = TRUE)
    within_50kb <- sum(distances <= 50000, na.rm = TRUE)
    within_100kb <- sum(distances <= 100000, na.rm = TRUE)

    cat(sprintf("\nDistance distribution:\n"))
    cat(sprintf("  Within 10 kb:    %s (%.1f%%)\n",
                format(within_10kb, big.mark = ","),
                100 * within_10kb / nrow(links_df)))
    cat(sprintf("  Within 50 kb:    %s (%.1f%%)\n",
                format(within_50kb, big.mark = ","),
                100 * within_50kb / nrow(links_df)))
    cat(sprintf("  Within 100 kb:   %s (%.1f%%)\n",
                format(within_100kb, big.mark = ","),
                100 * within_100kb / nrow(links_df)))
  }
}

################################################################################
# Complete
################################################################################

cat("\n================================================================================\n")
cat("PEAK-GENE LINKAGE COMPLETE\n")
cat("================================================================================\n\n")

cat(sprintf("Cell type:         %s\n", CELLTYPE_TARGET))
cat(sprintf("Genotype:          %s\n", GENOTYPE_TARGET))
cat(sprintf("Total cells:       %s (Ctrl: %s, Mut: %s)\n",
            format(n_cells_total, big.mark = ","),
            format(n_cells_ctrl, big.mark = ","),
            format(n_cells_mut, big.mark = ",")))
cat(sprintf("Total peaks:       %s\n", format(n_peaks_after, big.mark = ",")))

if (exists("links_df") && nrow(links_df) > 0) {
  cat(sprintf("Peak-gene links:   %s\n", format(nrow(links_df), big.mark = ",")))
  cat(sprintf("Unique genes:      %s\n", format(n_distinct(links_df$gene), big.mark = ",")))
} else {
  cat(sprintf("Peak-gene links:   0\n"))
}

cat(sprintf("\nOutput file:       %s\n", basename(output_file)))
cat(sprintf("File size:         %.2f MB\n", file.size(output_file) / 1024^2))

cat("\n✓ Analysis complete!\n")
