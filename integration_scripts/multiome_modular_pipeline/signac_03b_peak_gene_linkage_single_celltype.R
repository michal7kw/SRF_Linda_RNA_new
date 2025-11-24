#!/usr/bin/env Rscript

################################################################################
# Signac Multiome Analysis - Step 3b: Peak-Gene Linkage for Single Cell Type
# (L2 Cell Type Version)
################################################################################
#
# This script performs peak-gene linkage analysis for a SINGLE cell type × genotype
# combination, allowing for parallel execution via SLURM array jobs.
#
# WHY THIS SCRIPT EXISTS:
# The main Step 3 (signac_03_celltype_analysis.R) skips peak-gene linkage
# when a cell type has >30,000 peaks (due to runtime concerns). This script
# removes that limitation and processes individual cell types with more resources.
#
# USAGE:
#   Rscript signac_03b_peak_gene_linkage_single_celltype.R <CELLTYPE> <GENOTYPE>
#
# ARGUMENTS:
#   CELLTYPE: L2 cell type name (e.g., "Astrocytes", "Mature GC", "CA1")
#   GENOTYPE: Genotype name (e.g., "Nestin" or "Emx1")
#
# EXAMPLES:
#   Rscript signac_03b_peak_gene_linkage_single_celltype.R "Astrocytes" "Nestin"
#   Rscript signac_03b_peak_gene_linkage_single_celltype.R "Mature GC" "Emx1"
#
# INPUTS:
#   - signac_results/integrated_seurat_processed.rds
#
# OUTPUTS:
#   - signac_results/celltype_results/peak_gene_links/{CELLTYPE}_{GENOTYPE}_peak_gene_links.csv
#
################################################################################

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("\n")
  cat("ERROR: Incorrect number of arguments\n\n")
  cat("USAGE:\n")
  cat("  Rscript signac_03b_peak_gene_linkage_single_celltype.R <CELLTYPE> <GENOTYPE>\n\n")
  cat("ARGUMENTS:\n")
  cat("  CELLTYPE: L2 cell type name (e.g., 'Astrocytes', 'Mature GC', 'CA1')\n")
  cat("  GENOTYPE: Genotype name ('Nestin' or 'Emx1')\n\n")
  cat("EXAMPLES:\n")
  cat("  Rscript signac_03b_peak_gene_linkage_single_celltype.R 'Astrocytes' 'Nestin'\n")
  cat("  Rscript signac_03b_peak_gene_linkage_single_celltype.R 'Mature GC' 'Emx1'\n\n")
  quit(status = 1)
}

CELLTYPE_TARGET <- args[1]
GENOTYPE_TARGET <- args[2]

# Load required libraries
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(GenomicRanges)
  library(dplyr)
})

cat("================================================================================\n")
cat("  PEAK-GENE LINKAGE: SINGLE CELL TYPE (L2 VERSION)\n")
cat("================================================================================\n\n")
cat(sprintf("Target cell type: %s\n", CELLTYPE_TARGET))
cat(sprintf("Target genotype:  %s\n", GENOTYPE_TARGET))
cat("\n")

################################################################################
# Configuration
################################################################################

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
INPUT_DIR <- file.path(BASE_DIR, "signac_results")
OUTPUT_DIR <- INPUT_DIR
RESULTS_DIR <- file.path(OUTPUT_DIR, "celltype_results")

# Ensure output directories exist
dir.create(file.path(RESULTS_DIR, "peak_gene_links"), recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
CELLTYPE_COL <- "cell_type_L2_new"  # L2 cell type annotation column
CONDITION_COL <- "condition"        # Ctrl vs Mut
GENOTYPE_COL <- "genotype"          # Nestin vs Emx1
MIN_CELLS_PER_GROUP <- 30           # Minimum cells per condition

# Peak-gene linkage parameters (OPTIMIZED for comprehensive analysis)
PEAK_GENE_DISTANCE <- 200000         # 100kb upstream/downstream
MIN_CORRELATION <- 0.05               # Minimum correlation threshold
PVALUE_CUTOFF <- 0.05                # P-value cutoff
MIN_CELLS_FOR_LINKAGE <- 10          # Minimum cells expressing both peak and gene
# NO MAX_PEAKS_FOR_LINKAGE - we'll process however many peaks exist

cat("Configuration:\n")
cat(sprintf("  Cell type column: %s\n", CELLTYPE_COL))
cat(sprintf("  Peak-gene distance: %d bp\n", PEAK_GENE_DISTANCE))
cat(sprintf("  Min correlation: %.2f\n", MIN_CORRELATION))
cat(sprintf("  Min cells: %d\n", MIN_CELLS_FOR_LINKAGE))
cat(sprintf("  P-value cutoff: %.3f\n", PVALUE_CUTOFF))
cat("\n")

################################################################################
# Load Data
################################################################################

cat("Loading Seurat object...\n")
SEURAT_FILE <- file.path(INPUT_DIR, "integrated_seurat_processed.rds")

if (!file.exists(SEURAT_FILE)) {
  cat(sprintf("ERROR: Seurat file not found: %s\n", SEURAT_FILE))
  quit(status = 1)
}

seurat_obj <- readRDS(SEURAT_FILE)
cat(sprintf("  ✓ Loaded: %d cells, %d genes (RNA), %d peaks (ATAC)\n",
            ncol(seurat_obj),
            nrow(seurat_obj[["RNA"]]),
            nrow(seurat_obj[["ATAC"]])))

################################################################################
# Add Metadata (if missing)
################################################################################

# Create condition column if missing
if (!CONDITION_COL %in% colnames(seurat_obj@meta.data)) {
  cat("Creating condition column...\n")
  seurat_obj@meta.data[[CONDITION_COL]] <- ifelse(
    grepl("Ctrl", seurat_obj@meta.data$orig.ident, ignore.case = TRUE),
    "Ctrl",
    "Mut"
  )
}

# Create genotype column if missing
if (!GENOTYPE_COL %in% colnames(seurat_obj@meta.data)) {
  cat("Creating genotype column...\n")
  seurat_obj@meta.data[[GENOTYPE_COL]] <- ifelse(
    grepl("Nestin", seurat_obj@meta.data$orig.ident, ignore.case = TRUE),
    "Nestin",
    "Emx1"
  )
}

################################################################################
# Filter to Target Cell Type and Genotype
################################################################################

cat(sprintf("\nFiltering to %s × %s...\n", CELLTYPE_TARGET, GENOTYPE_TARGET))

# Get cells matching both cell type and genotype
target_cells <- seurat_obj@meta.data %>%
  filter(.data[[CELLTYPE_COL]] == CELLTYPE_TARGET,
         .data[[GENOTYPE_COL]] == GENOTYPE_TARGET) %>%
  rownames()

if (length(target_cells) == 0) {
  cat(sprintf("ERROR: No cells found for %s × %s\n", CELLTYPE_TARGET, GENOTYPE_TARGET))
  quit(status = 1)
}

# Subset Seurat object
seurat_ct <- subset(seurat_obj, cells = target_cells)

cat(sprintf("  ✓ Filtered to %d cells\n", ncol(seurat_ct)))

# Check sample distribution
cat("\nSample distribution:\n")
sample_dist <- table(seurat_ct@meta.data[[CONDITION_COL]])
print(sample_dist)

# Validate minimum cells
n_ctrl <- sum(seurat_ct@meta.data[[CONDITION_COL]] == "Ctrl")
n_mut <- sum(seurat_ct@meta.data[[CONDITION_COL]] == "Mut")

if (n_ctrl < MIN_CELLS_PER_GROUP || n_mut < MIN_CELLS_PER_GROUP) {
  cat(sprintf("\nWARNING: Insufficient cells (Ctrl: %d, Mut: %d, minimum: %d)\n",
              n_ctrl, n_mut, MIN_CELLS_PER_GROUP))
  cat("Skipping linkage analysis.\n")
  quit(status = 0)
}

################################################################################
# Filter Invalid Peaks (Critical for mm10)
################################################################################

cat("\nValidating peak coordinates...\n")

# Get mm10 chromosome sizes
genome <- BSgenome.Mmusculus.UCSC.mm10
genome_info <- seqlengths(genome)

# Get peak coordinates
peaks <- granges(seurat_ct[["ATAC"]])
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
seurat_ct[["ATAC"]] <- subset(seurat_ct[["ATAC"]], features = valid_peak_names)

n_peaks_after <- nrow(seurat_ct[["ATAC"]])

cat(sprintf("  Peaks before filtering: %d\n", n_peaks_before))
cat(sprintf("  Peaks after filtering:  %d\n", n_peaks_after))
cat(sprintf("  Removed: %d invalid peaks\n", n_peaks_before - n_peaks_after))

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
    seurat_ct <- RegionStats(seurat_ct, genome = genome)
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

cat(sprintf("Computing peak-gene links for %s × %s...\n", CELLTYPE_TARGET, GENOTYPE_TARGET))
cat(sprintf("  Peaks: %s\n", format(n_peaks_after, big.mark = ",")))
cat(sprintf("  Genes: %s\n", format(nrow(seurat_ct[["RNA"]]), big.mark = ",")))
cat(sprintf("  Peak-gene distance: %d bp (±%d kb)\n",
            PEAK_GENE_DISTANCE, PEAK_GENE_DISTANCE / 1000))
cat(sprintf("  Min correlation: %.2f\n", MIN_CORRELATION))
cat(sprintf("  This may take 1-3 hours depending on cell count and peak number...\n"))
cat("\n")

# Set default assay to ATAC
DefaultAssay(seurat_ct) <- "ATAC"

# Perform peak-gene linkage
start_time <- Sys.time()

seurat_ct <- tryCatch({
  LinkPeaks(
    object = seurat_ct,
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
}, error = function(e) {
  cat(sprintf("ERROR in LinkPeaks: %s\n", e$message))
  quit(status = 1)
})

end_time <- Sys.time()
elapsed_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
cat(sprintf("\n  ✓ LinkPeaks completed in %.1f minutes\n", elapsed_time))

# Extract links from the Seurat object
cat("\nExtracting peak-gene links...\n")
links <- Links(seurat_ct)

if (is.null(links) || length(links) == 0) {
  cat("\n⚠ No significant peak-gene links found\n")
  cat("This could mean:\n")
  cat("  1. No strong correlations between peaks and genes in this cell type\n")
  cat("  2. Insufficient cells for reliable correlation\n")
  cat("  3. Technical issues with peak/gene accessibility\n")

  # Create empty output file
  empty_df <- data.frame(
    seqnames = character(),
    start = integer(),
    end = integer(),
    width = integer(),
    strand = character(),
    gene = character(),
    score = numeric(),
    pvalue = numeric(),
    celltype = character(),
    genotype = character(),
    stringsAsFactors = FALSE
  )

  output_file <- file.path(
    RESULTS_DIR,
    "peak_gene_links",
    sprintf("%s_%s_peak_gene_links.csv", CELLTYPE_TARGET, GENOTYPE_TARGET)
  )
  write.csv(empty_df, output_file, row.names = FALSE)
  cat(sprintf("\n✓ Saved empty results: %s\n", basename(output_file)))
  quit(status = 0)
}

# Convert Links to data frame
cat("\nProcessing links...\n")
links_df <- as.data.frame(links)
links_df$celltype <- CELLTYPE_TARGET
links_df$genotype <- GENOTYPE_TARGET

cat(sprintf("  ✓ Found %d significant peak-gene links\n", nrow(links_df)))
cat(sprintf("  Unique genes: %d\n", length(unique(links_df$gene))))

# Save results
output_file <- file.path(
  RESULTS_DIR,
  "peak_gene_links",
  sprintf("%s_%s_peak_gene_links.csv", CELLTYPE_TARGET, GENOTYPE_TARGET)
)

write.csv(links_df, output_file, row.names = FALSE)

cat(sprintf("\n✓ Saved: %s\n", basename(output_file)))
cat(sprintf("  File size: %.2f MB\n", file.size(output_file) / 1024^2))

################################################################################
# Summary Statistics
################################################################################

cat("\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n\n")

cat(sprintf("Cell type:         %s\n", CELLTYPE_TARGET))
cat(sprintf("Genotype:          %s\n", GENOTYPE_TARGET))
cat(sprintf("Total cells:       %s (Ctrl: %s, Mut: %s)\n",
            format(ncol(seurat_ct), big.mark = ","),
            format(n_ctrl, big.mark = ","),
            format(n_mut, big.mark = ",")))
cat(sprintf("Total peaks:       %s\n", format(n_peaks_after, big.mark = ",")))
cat(sprintf("Peak-gene links:   %s\n", format(nrow(links_df), big.mark = ",")))
cat(sprintf("Unique genes:      %s\n", format(length(unique(links_df$gene)), big.mark = ",")))
cat(sprintf("Unique peaks:      %s\n", format(length(unique(rownames(links_df))), big.mark = ",")))

if (nrow(links_df) > 0) {
  cat(sprintf("\nCorrelation statistics:\n"))
  cat(sprintf("  Mean:            %.3f\n", mean(links_df$score, na.rm = TRUE)))
  cat(sprintf("  Median:          %.3f\n", median(links_df$score, na.rm = TRUE)))
  cat(sprintf("  Range:           %.3f to %.3f\n",
              min(links_df$score), max(links_df$score)))

  cat(sprintf("\nP-value range:     %.2e to %.2e\n",
              min(links_df$pvalue), max(links_df$pvalue)))

  # Calculate peak-gene distances
  if ("start" %in% colnames(links_df)) {
    cat(sprintf("\nPeak-gene distance statistics:\n"))

    tryCatch({
      # Check if we have peak names in rownames or a peak column
      if (length(rownames(links_df)) > 0 && grepl("-", rownames(links_df)[1])) {
        # Extract peak start positions from rownames (format: chr-start-end)
        peak_starts <- as.numeric(sub(".*-(\\d+)-\\d+$", "\\1", rownames(links_df)))
      } else if ("peak" %in% colnames(links_df)) {
        # Extract from peak column if available
        peak_starts <- as.numeric(sub(".*-(\\d+)-\\d+$", "\\1", links_df$peak))
      } else {
        stop("Could not find peak coordinates in expected format")
      }

      # Gene start positions are in links_df$start
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

    }, error = function(e) {
      cat(sprintf("  Warning: Could not compute distance statistics: %s\n", e$message))
      cat(sprintf("           Peak coordinate columns may not be in expected format\n"))
    })
  }
}

cat("\n✓ Analysis complete!\n")
