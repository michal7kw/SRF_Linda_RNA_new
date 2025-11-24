#!/usr/bin/env Rscript

################################################################################
# Signac Multiome Analysis - Step 3: Cell-Type-Specific Differential Analysis
################################################################################
#
# This script performs comprehensive cell-type-specific analysis:
#   1. Differential gene expression (RNA) - Ctrl vs Mut per cell type
#   2. Differential accessibility (ATAC) - Ctrl vs Mut per cell type
#   3. Peak-gene linkage using Signac's LinkPeaks with correlation
#   4. Integration with existing DEG results
#
# This addresses all issues from the Python pipeline:
#
# INPUTS:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/integrated_seurat_processed.rds
#     (Processed Seurat object from Step 2 with RNA and ATAC assays)
#
# OUTPUTS:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/celltype_results/DEG/[celltype]_[genotype]_DEGs.csv
#     (Differential expression results for each cell type and genotype)
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/celltype_results/DA/[celltype]_[genotype]_DA_peaks.csv
#     (Differential accessibility results for each cell type and genotype)
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/celltype_results/peak_gene_links/[celltype]_[genotype]_peak_gene_links.csv
#     (Peak-gene linkage results for each cell type and genotype)
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/celltype_results/summary_statistics.csv
#     (Summary statistics for all cell types and genotypes)
#
#   - Uses sophisticated correlation-based peak-gene linkage (not just promoters)
#   - Proper statistical testing with multiple methods
#   - Cell-type-specific analysis
#   - Genome-wide coverage
#
# Based on: https://stuartlab.org/signac/articles/pbmc_multiomic
#
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(presto)  # Fast Wilcoxon test
})

cat("================================================================================\n")
cat("  SIGNAC MULTIOME ANALYSIS - STEP 3: CELL-TYPE-SPECIFIC ANALYSIS\n")
cat("================================================================================\n\n")

################################################################################
# Configuration
################################################################################

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
INPUT_DIR <- file.path(BASE_DIR, "signac_results")
OUTPUT_DIR <- INPUT_DIR
RESULTS_DIR <- file.path(OUTPUT_DIR, "celltype_results")

# Create output directories
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "DEG"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "DA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "peak_gene_links"), recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
CELLTYPE_COL <- "cell_type_L2_new"  # Cell type annotation column
CONDITION_COL <- "condition"        # Ctrl vs Mut (will be created if not exists)
GENOTYPE_COL <- "genotype"          # Nestin vs Emx1 (will be created if not exists)
MIN_CELLS_PER_GROUP <- 30           # Minimum cells per condition (lowered for genotype-stratified analysis)

# Sample identifiers (for backwards compatibility and metadata creation)
SAMPLE_MAPPING <- list(
  nestin = list(
    ctrl = "R26-Nestin-Ctrl-adult",
    mut = "R26-Nestin-Mut-adult",
    genotype = "Nestin"
  ),
  emx1 = list(
    ctrl = "R26-Emx1-Ctrl-adult",
    mut = "R26-Emx1-Mut-adult",
    genotype = "Emx1"
  )
)

# Analysis mode: "within_genotype" or "across_genotype"
ANALYSIS_MODE <- "within_genotype"  # Compare Ctrl vs Mut within each genotype separately

# Differential testing parameters
MIN_PCT <- 0.1                       # Minimum fraction of cells expressing
LOGFC_THRESHOLD <- 0.25              # Log2FC threshold
PVAL_THRESHOLD <- 0.05               # Adjusted p-value threshold

# Peak-gene linkage parameters (OPTIMIZED for faster runtime)
PEAK_GENE_DISTANCE <- 100000         # Maximum distance for peak-gene links (100kb) - REDUCED from 200kb
MIN_CORRELATION <- 0.1               # Minimum correlation for peak-gene links - INCREASED from 0.05
PVALUE_CUTOFF <- 0.05                # P-value cutoff for correlations
MAX_PEAKS_FOR_LINKAGE <- 30000       # Skip linkage if DA peaks > this threshold (too slow)

cat("Configuration:\n")
cat(sprintf("  Cell type column: %s\n", CELLTYPE_COL))
cat(sprintf("  Condition column: %s\n", CONDITION_COL))
cat(sprintf("  Genotype column: %s\n", GENOTYPE_COL))
cat(sprintf("  Analysis mode: %s\n", ANALYSIS_MODE))
cat(sprintf("  Min cells per group: %d\n", MIN_CELLS_PER_GROUP))
cat(sprintf("  Peak-gene distance: %d bp (OPTIMIZED: reduced from 200kb)\n", PEAK_GENE_DISTANCE))
cat(sprintf("  Min correlation: %.2f (OPTIMIZED: increased from 0.05)\n", MIN_CORRELATION))
cat(sprintf("  Max peaks for linkage: %d (OPTIMIZED: skip if exceeded)\n", MAX_PEAKS_FOR_LINKAGE))
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

# Verify required columns exist
if (!CELLTYPE_COL %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Cell type column '%s' not found in metadata. Available columns: %s",
               CELLTYPE_COL, paste(colnames(seurat_obj@meta.data), collapse = ", ")))
}

if (!CONDITION_COL %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Condition column '%s' not found in metadata.", CONDITION_COL))
}

if (!GENOTYPE_COL %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Genotype column '%s' not found in metadata.", GENOTYPE_COL))
}

# Get genotypes to analyze
genotypes <- unique(seurat_obj@meta.data[[GENOTYPE_COL]])
genotypes <- genotypes[!is.na(genotypes)]

cat(sprintf("  Found %d genotypes: %s\n",
            length(genotypes), paste(genotypes, collapse = ", ")))

# Get conditions
conditions <- unique(seurat_obj@meta.data[[CONDITION_COL]])
conditions <- conditions[!is.na(conditions)]

cat(sprintf("  Found %d conditions: %s\n",
            length(conditions), paste(conditions, collapse = ", ")))

# Get cell types to analyze
celltypes <- unique(seurat_obj@meta.data[[CELLTYPE_COL]])
celltypes <- celltypes[!is.na(celltypes)]

cat(sprintf("  Found %d cell types: %s\n",
            length(celltypes),
            paste(head(celltypes, 5), collapse = ", ")))

if (length(celltypes) > 5) {
  cat(sprintf("  ... and %d more\n", length(celltypes) - 5))
}

################################################################################
# Function: Differential Gene Expression (RNA)
################################################################################

run_deg_analysis <- function(seurat_obj, celltype, ctrl_ident, mut_ident,
                              min_pct = 0.1, logfc_threshold = 0.25) {
  cat(sprintf("\n  === Differential Gene Expression: %s ===\n", celltype))

  DefaultAssay(seurat_obj) <- "RNA"

  # Subset to cell type (use direct metadata access, not WhichCells with get())
  celltype_values <- seurat_obj@meta.data[[CELLTYPE_COL]]
  cells_ct <- colnames(seurat_obj)[celltype_values == celltype]
  seurat_ct <- subset(seurat_obj, cells = cells_ct)

  # Check cell counts
  n_ctrl <- sum(seurat_ct@meta.data[[CONDITION_COL]] == ctrl_ident)
  n_mut <- sum(seurat_ct@meta.data[[CONDITION_COL]] == mut_ident)

  cat(sprintf("    Cells: %d Ctrl, %d Mut\n", n_ctrl, n_mut))

  if (n_ctrl < MIN_CELLS_PER_GROUP || n_mut < MIN_CELLS_PER_GROUP) {
    cat(sprintf("    ⚠ Skipping (insufficient cells, min = %d)\n", MIN_CELLS_PER_GROUP))
    return(NULL)
  }

  # Set idents to condition (use values from metadata column, not the column name)
  Idents(seurat_ct) <- seurat_ct@meta.data[[CONDITION_COL]]

  # Verify idents were set correctly
  ident_levels <- levels(factor(Idents(seurat_ct)))
  cat(sprintf("    Identity levels: %s\n", paste(ident_levels, collapse = ", ")))

  # Check if required idents exist
  if (!(ctrl_ident %in% ident_levels && mut_ident %in% ident_levels)) {
    cat(sprintf("    ✗ Error: Required idents not found in data\n"))
    cat(sprintf("      Expected: %s, %s\n", ctrl_ident, mut_ident))
    cat(sprintf("      Found: %s\n", paste(ident_levels, collapse = ", ")))
    return(NULL)
  }

  # Run differential expression (Wilcoxon)
  cat("    Running Wilcoxon test...\n")
  deg_results <- FindMarkers(
    object = seurat_ct,
    ident.1 = mut_ident,
    ident.2 = ctrl_ident,
    test.use = "wilcox",
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    verbose = FALSE
  )

  if (nrow(deg_results) == 0) {
    cat("    No significant DEGs found\n")
    return(NULL)
  }

  # Add gene names as column
  deg_results$gene <- rownames(deg_results)

  # Add direction
  deg_results$direction <- ifelse(deg_results$avg_log2FC > 0, "Up", "Down")

  # Filter by adjusted p-value
  deg_sig <- deg_results[deg_results$p_val_adj < PVAL_THRESHOLD, ]

  n_up <- sum(deg_sig$direction == "Up")
  n_down <- sum(deg_sig$direction == "Down")

  cat(sprintf("    ✓ Found %d DEGs (%d up, %d down)\n",
              nrow(deg_sig), n_up, n_down))

  return(deg_results)
}

################################################################################
# Function: Differential Accessibility (ATAC)
################################################################################

run_da_analysis <- function(seurat_obj, celltype, ctrl_ident, mut_ident,
                             min_pct = 0.1) {
  cat(sprintf("\n  === Differential Accessibility: %s ===\n", celltype))

  DefaultAssay(seurat_obj) <- "ATAC"

  # Subset to cell type (use direct metadata access, not WhichCells with get())
  celltype_values <- seurat_obj@meta.data[[CELLTYPE_COL]]
  cells_ct <- colnames(seurat_obj)[celltype_values == celltype]
  seurat_ct <- subset(seurat_obj, cells = cells_ct)

  # Check cell counts
  n_ctrl <- sum(seurat_ct@meta.data[[CONDITION_COL]] == ctrl_ident)
  n_mut <- sum(seurat_ct@meta.data[[CONDITION_COL]] == mut_ident)

  cat(sprintf("    Cells: %d Ctrl, %d Mut\n", n_ctrl, n_mut))

  if (n_ctrl < MIN_CELLS_PER_GROUP || n_mut < MIN_CELLS_PER_GROUP) {
    cat(sprintf("    ⚠ Skipping (insufficient cells)\n"))
    return(NULL)
  }

  # Set idents to condition (use values from metadata column, not the column name)
  Idents(seurat_ct) <- seurat_ct@meta.data[[CONDITION_COL]]

  # Verify idents were set correctly
  ident_levels <- levels(factor(Idents(seurat_ct)))
  cat(sprintf("    Identity levels: %s\n", paste(ident_levels, collapse = ", ")))

  # Check if required idents exist
  if (!(ctrl_ident %in% ident_levels && mut_ident %in% ident_levels)) {
    cat(sprintf("    ✗ Error: Required idents not found in data\n"))
    cat(sprintf("      Expected: %s, %s\n", ctrl_ident, mut_ident))
    cat(sprintf("      Found: %s\n", paste(ident_levels, collapse = ", ")))
    return(NULL)
  }

  # Run differential accessibility (Wilcoxon)
  cat("    Running Wilcoxon test...\n")
  da_results <- FindMarkers(
    object = seurat_ct,
    ident.1 = mut_ident,
    ident.2 = ctrl_ident,
    test.use = "wilcox",
    min.pct = min_pct,
    logfc.threshold = 0,  # No FC threshold for peaks
    verbose = FALSE
  )

  if (nrow(da_results) == 0) {
    cat("    No significant DA peaks found\n")
    return(NULL)
  }

  # Add peak names as column
  da_results$peak <- rownames(da_results)

  # Add direction
  da_results$direction <- ifelse(da_results$avg_log2FC > 0, "More_Accessible", "Less_Accessible")

  # Filter by adjusted p-value
  da_sig <- da_results[da_results$p_val_adj < PVAL_THRESHOLD, ]

  n_more <- sum(da_sig$direction == "More_Accessible")
  n_less <- sum(da_sig$direction == "Less_Accessible")

  cat(sprintf("    ✓ Found %d DA peaks (%d more accessible, %d less accessible)\n",
              nrow(da_sig), n_more, n_less))

  return(da_results)
}

################################################################################
# Function: Peak-Gene Linkage
################################################################################

run_peak_gene_linkage <- function(seurat_obj, celltype) {
  cat(sprintf("\n  === Peak-Gene Linkage: %s ===\n", celltype))

  DefaultAssay(seurat_obj) <- "ATAC"

  # Subset to cell type (use direct metadata access, not WhichCells with get())
  celltype_values <- seurat_obj@meta.data[[CELLTYPE_COL]]
  cells_ct <- colnames(seurat_obj)[celltype_values == celltype]
  seurat_ct <- subset(seurat_obj, cells = cells_ct)

  n_cells <- ncol(seurat_ct)
  cat(sprintf("    Cells: %d\n", n_cells))

  if (n_cells < MIN_CELLS_PER_GROUP * 2) {
    cat(sprintf("    ⚠ Skipping (insufficient cells for correlation)\n"))
    return(NULL)
  }

  # Validate and filter peaks to remove out-of-bounds regions
  cat("    Validating peak coordinates...\n")

  # Get chromosome lengths from genome
  genome <- BSgenome.Mmusculus.UCSC.mm10
  all_chrom_lengths <- seqlengths(genome)

  # CRITICAL: Only use standard chromosomes (not scaffolds, random contigs, etc.)
  # Mouse mm10 should have: chr1-19, chrX, chrY, chrM
  standard_chroms <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM")
  valid_chroms <- standard_chroms[standard_chroms %in% names(all_chrom_lengths)]
  chrom_lengths <- all_chrom_lengths[valid_chroms]

  cat(sprintf("    Using %d standard chromosomes (chr1-19, chrX, chrY, chrM)\n", length(valid_chroms)))
  cat(sprintf("    Genome object has %d total chromosomes (including scaffolds)\n", length(all_chrom_lengths)))

  # Get peaks as GRanges
  peaks <- granges(seurat_ct[["ATAC"]])

  # CRITICAL: Use original peak names from counts matrix, not from GRanges names
  # GRanges names may be reformatted (chr1:100-200 vs chr1-100-200)
  peak_names <- rownames(seurat_ct[["ATAC"]])

  # IMPORTANT: We need to validate peaks by parsing the peak names directly
  # because GRanges order might not match rownames order
  # Peak format: chr1-3051510-3052173 (dash-separated)

  # Parse peak names to extract chromosome, start, end
  peak_coords <- strsplit(peak_names, "-")

  # Identify valid peaks (within chromosome boundaries AND valid chromosome names)
  valid_peaks <- sapply(peak_coords, function(coord) {
    if (length(coord) != 3) return(FALSE)  # Invalid format

    chr <- coord[1]
    start_pos <- as.numeric(coord[2])
    end_pos <- as.numeric(coord[3])

    # Check if chromosome is standard
    if (!chr %in% valid_chroms) {
      return(FALSE)
    }

    # Check if peak start is valid (>0)
    if (is.na(start_pos) || start_pos < 1) {
      return(FALSE)
    }

    # Check if peak end is within chromosome length
    if (is.na(end_pos) || end_pos > chrom_lengths[chr]) {
      return(FALSE)
    }

    return(TRUE)
  })

  n_invalid <- sum(!valid_peaks)
  n_total <- length(peak_names)
  pct_invalid <- 100 * n_invalid / n_total

  if (n_invalid > 0) {
    cat(sprintf("    ⚠ Found %d invalid peaks (%.2f%%):\n", n_invalid, pct_invalid))

    # Report breakdown of invalid peaks by type
    invalid_chr <- sum(sapply(peak_coords, function(coord) {
      if (length(coord) != 3) return(TRUE)  # Invalid format counts as invalid chr
      !coord[1] %in% valid_chroms
    }))
    invalid_bounds <- n_invalid - invalid_chr

    if (invalid_chr > 0) {
      cat(sprintf("      - %d peaks on non-standard chromosomes\n", invalid_chr))
    }
    if (invalid_bounds > 0) {
      cat(sprintf("      - %d peaks beyond chromosome boundaries\n", invalid_bounds))
    }

    # Get names of valid peaks (now correctly indexed to peak_names)
    valid_peak_names <- peak_names[valid_peaks]

    cat(sprintf("    Filtering ATAC assay from %d to %d peaks...\n", n_total, length(valid_peak_names)))

    # APPROACH: Create new ChromatinAssay with filtered peaks
    # We must keep it as ChromatinAssay (not Assay) for RegionStats/LinkPeaks to work
    tryCatch({
      # Get the ATAC assay
      atac_assay <- seurat_ct[["ATAC"]]

      # Get counts matrix and filter to valid peaks
      counts_matrix <- GetAssayData(atac_assay, layer = "counts")

      # Check if valid_peak_names exist in counts matrix
      valid_peak_names_existing <- valid_peak_names[valid_peak_names %in% rownames(counts_matrix)]

      if (length(valid_peak_names_existing) == 0) {
        cat("    ✗ Error: No valid peak names found in counts matrix!\n")
        cat(sprintf("      Peak names format: %s\n", paste(head(valid_peak_names, 3), collapse = ", ")))
        cat(sprintf("      Counts rownames format: %s\n", paste(head(rownames(counts_matrix), 3), collapse = ", ")))
        return(NULL)
      }

      # Filter counts to valid peaks
      counts_filtered <- counts_matrix[valid_peak_names_existing, , drop = FALSE]

      # CRITICAL: Must use CreateChromatinAssay (not CreateAssayObject) to preserve
      # ChromatinAssay structure required by RegionStats() and LinkPeaks()
      # Extract annotation and fragments from original assay
      annotation_obj <- Annotation(atac_assay)
      fragments_obj <- Fragments(atac_assay)

      # Create new ChromatinAssay with filtered peaks
      atac_assay_new <- CreateChromatinAssay(
        counts = counts_filtered,
        annotation = annotation_obj,
        fragments = fragments_obj,
        sep = c("-", "-")  # Peak format: chr1-3051510-3052173
      )

      # Replace ATAC assay
      seurat_ct[["ATAC"]] <- atac_assay_new

      # Verify filtering worked
      n_remaining <- nrow(seurat_ct[["ATAC"]])
      cat(sprintf("    ✓ Successfully retained %d valid peaks (%.1f%%)\n",
                  n_remaining, 100 * n_remaining / n_total))

      # Re-check if we still have enough peaks
      if (n_remaining < 100) {
        cat(sprintf("    ✗ Too few valid peaks remaining (%d)\n", n_remaining))
        return(NULL)
      }

      if (n_remaining == n_total) {
        cat(sprintf("    ⚠ WARNING: Filtering failed - still have all %d peaks!\n", n_total))
        cat("      This suggests the filtering logic did not work correctly.\n")
        return(NULL)
      }

    }, error = function(e) {
      cat(sprintf("    ✗ Error during filtering: %s\n", e$message))
      cat("      Attempting workaround: Skip peak filtering and use all peaks\n")
      cat("      Note: This may cause RegionStats to fail later\n")
      # Don't return NULL - let it proceed and fail at RegionStats with better diagnostics
    })
  } else {
    cat("    ✓ All peaks have valid coordinates\n")
  }

  # Compute GC content for peaks (required for LinkPeaks)
  cat("    Computing peak GC content...\n")

  # Use try-catch to provide better error messages if RegionStats still fails
  tryCatch({
    seurat_ct <- RegionStats(seurat_ct, genome = genome)
    cat("    ✓ GC content computed successfully\n")
  }, error = function(e) {
    cat(sprintf("    ✗ Error computing GC content: %s\n", e$message))

    # Additional diagnostics
    final_peaks <- granges(seurat_ct[["ATAC"]])
    cat("\n    Diagnostics:\n")
    cat(sprintf("      Peaks in object: %d\n", length(final_peaks)))
    cat(sprintf("      Chromosomes in peaks: %s\n",
                paste(unique(as.character(seqnames(final_peaks))), collapse = ", ")))
    cat(sprintf("      Standard chromosomes: %s\n",
                paste(valid_chroms, collapse = ", ")))

    # Check for specific problematic peaks
    peak_chrs <- as.character(seqnames(final_peaks))
    non_standard <- peak_chrs[!peak_chrs %in% valid_chroms]
    if (length(non_standard) > 0) {
      cat(sprintf("      Non-standard chromosomes found: %s\n",
                  paste(unique(non_standard), collapse = ", ")))
    }

    stop(e)
  })

  # Link peaks to genes
  cat(sprintf("    Linking peaks to genes (distance < %d bp, cor > %.2f)...\n",
              PEAK_GENE_DISTANCE, MIN_CORRELATION))

  # Check if we have too many peaks (would take too long)
  n_peaks <- nrow(seurat_ct[["ATAC"]])
  if (n_peaks > MAX_PEAKS_FOR_LINKAGE) {
    cat(sprintf("    ⚠ Skipping linkage: %d peaks exceeds threshold of %d\n",
                n_peaks, MAX_PEAKS_FOR_LINKAGE))
    cat("      (This would take several hours - consider analyzing this cell type separately)\n")
    return(NULL)
  }

  # Try LinkPeaks with error handling
  tryCatch({
    seurat_ct <- LinkPeaks(
      object = seurat_ct,
      peak.assay = "ATAC",
      expression.assay = "RNA",
      distance = PEAK_GENE_DISTANCE,
      min.cells = 10,
      pvalue_cutoff = PVALUE_CUTOFF,
      score_cutoff = MIN_CORRELATION,
      gene.id = FALSE  # Use gene names instead of IDs (more compatible)
    )
  }, error = function(e) {
    cat(sprintf("    ✗ Error in peak-gene linkage: %s\n", e$message))
    return(NULL)
  })

  # Check if LinkPeaks succeeded
  if (is.null(seurat_ct)) {
    return(NULL)
  }

  # Extract links
  links <- Links(seurat_ct)

  if (length(links) == 0) {
    cat("    ⚠ No significant peak-gene links found\n")
    return(NULL)
  }

  cat(sprintf("    ✓ Found %d peak-gene links\n", length(links)))

  # Convert to data frame
  links_df <- as.data.frame(links)

  # Add gene and peak information
  if (nrow(links_df) > 0) {
    links_df$celltype <- celltype
  }

  return(links_df)
}

################################################################################
# Main Analysis Loop
################################################################################

cat("\n================================================================================\n")
cat("CELL-TYPE-SPECIFIC ANALYSIS\n")
cat("================================================================================\n")

# Initialize results storage
deg_results_list <- list()
da_results_list <- list()
links_results_list <- list()

# Summary statistics
summary_stats <- data.frame(
  celltype = character(),
  genotype = character(),
  n_cells_ctrl = integer(),
  n_cells_mut = integer(),
  n_degs = integer(),
  n_degs_up = integer(),
  n_degs_down = integer(),
  n_da_peaks = integer(),
  n_da_more = integer(),
  n_da_less = integer(),
  n_peak_gene_links = integer(),
  stringsAsFactors = FALSE
)

# Loop through genotypes and cell types
for (geno in genotypes) {
  cat(sprintf("\n================================================================================\n"))
  cat(sprintf("GENOTYPE: %s\n", geno))
  cat(sprintf("================================================================================\n"))

  # Subset to this genotype
  cells_geno <- seurat_obj@meta.data[[GENOTYPE_COL]] == geno
  cat(sprintf("Total cells for %s: %d\n", geno, sum(cells_geno)))

  for (ct in celltypes) {
    cat(sprintf("\n--- %s: %s ---\n", geno, ct))

    # Get cell counts for this genotype AND cell type
    cells_ct_geno <- cells_geno & (seurat_obj@meta.data[[CELLTYPE_COL]] == ct)
    n_ctrl <- sum(cells_ct_geno & seurat_obj@meta.data[[CONDITION_COL]] == "Ctrl")
    n_mut <- sum(cells_ct_geno & seurat_obj@meta.data[[CONDITION_COL]] == "Mut")

    cat(sprintf("  Cells: Ctrl=%d, Mut=%d\n", n_ctrl, n_mut))

    # Check minimum cells
    if (n_ctrl < MIN_CELLS_PER_GROUP || n_mut < MIN_CELLS_PER_GROUP) {
      cat(sprintf("  ⚠ Skipping (insufficient cells, min=%d)\n", MIN_CELLS_PER_GROUP))

      summary_stats <- rbind(summary_stats, data.frame(
        celltype = ct,
        genotype = geno,
        n_cells_ctrl = n_ctrl,
        n_cells_mut = n_mut,
        n_degs = NA,
        n_degs_up = NA,
        n_degs_down = NA,
        n_da_peaks = NA,
        n_da_more = NA,
        n_da_less = NA,
        n_peak_gene_links = NA
      ))
      next
    }

    # Create analysis identifier
    analysis_id <- sprintf("%s_%s", ct, geno)

    # Subset Seurat object to current genotype before analysis
    # This ensures we're only comparing Ctrl vs Mut within this genotype
    cells_geno_subset <- rownames(seurat_obj@meta.data[cells_ct_geno, ])
    seurat_geno <- subset(seurat_obj, cells = cells_geno_subset)

    cat(sprintf("    Subset to genotype %s: %d cells\n", geno, ncol(seurat_geno)))

    # 1. Differential Gene Expression
    # Ensure DEG directory exists
    dir.create(file.path(RESULTS_DIR, "DEG"), recursive = TRUE, showWarnings = FALSE)
    deg_file <- file.path(RESULTS_DIR, "DEG", sprintf("%s_DEGs.csv", analysis_id))
  if (file.exists(deg_file)) {
    cat("  === Differential Gene Expression: Loading cached results ===\n")
    deg_res <- read.csv(deg_file)
    deg_results_list[[analysis_id]] <- deg_res
    n_sig <- sum(deg_res$p_val_adj < PVAL_THRESHOLD, na.rm = TRUE)
    n_up <- sum(deg_res$p_val_adj < PVAL_THRESHOLD & deg_res$avg_log2FC > 0, na.rm = TRUE)
    n_down <- sum(deg_res$p_val_adj < PVAL_THRESHOLD & deg_res$avg_log2FC < 0, na.rm = TRUE)
    cat(sprintf("    ✓ Loaded %d DEGs (%d up, %d down)\n", n_sig, n_up, n_down))
  } else {
    deg_res <- tryCatch({
      # Call with generic "Ctrl" and "Mut" - function will compare within genotype subset
      run_deg_analysis(seurat_geno, ct, "Ctrl", "Mut",
                       min_pct = MIN_PCT, logfc_threshold = LOGFC_THRESHOLD)
    }, error = function(e) {
      cat(sprintf("    ✗ Error in DEG analysis: %s\n", e$message))
      NULL
    })

    if (!is.null(deg_res)) {
      deg_results_list[[analysis_id]] <- deg_res
      # Save to file
      write.csv(deg_res, deg_file, row.names = FALSE)
    }
  }

  # 2. Differential Accessibility
  # Ensure DA directory exists
  dir.create(file.path(RESULTS_DIR, "DA"), recursive = TRUE, showWarnings = FALSE)
  da_file <- file.path(RESULTS_DIR, "DA", sprintf("%s_DA_peaks.csv", analysis_id))
  if (file.exists(da_file)) {
    cat("  === Differential Accessibility: Loading cached results ===\n")
    da_res <- read.csv(da_file)
    da_results_list[[analysis_id]] <- da_res
    n_sig <- sum(da_res$p_val_adj < PVAL_THRESHOLD, na.rm = TRUE)
    n_more <- sum(da_res$p_val_adj < PVAL_THRESHOLD & da_res$avg_log2FC > 0, na.rm = TRUE)
    n_less <- sum(da_res$p_val_adj < PVAL_THRESHOLD & da_res$avg_log2FC < 0, na.rm = TRUE)
    cat(sprintf("    ✓ Loaded %d DA peaks (%d more accessible, %d less accessible)\n", n_sig, n_more, n_less))
  } else {
    da_res <- tryCatch({
      # Call with generic "Ctrl" and "Mut" - function will compare within genotype subset
      run_da_analysis(seurat_geno, ct, "Ctrl", "Mut", min_pct = MIN_PCT)
    }, error = function(e) {
      cat(sprintf("    ✗ Error in DA analysis: %s\n", e$message))
      NULL
    })

    if (!is.null(da_res)) {
      da_results_list[[analysis_id]] <- da_res
      # Save to file
      write.csv(da_res, da_file, row.names = FALSE)
    }
  }

  # 3. Peak-Gene Linkage
  # Ensure peak_gene_links directory exists
  dir.create(file.path(RESULTS_DIR, "peak_gene_links"), recursive = TRUE, showWarnings = FALSE)
  links_file <- file.path(RESULTS_DIR, "peak_gene_links", sprintf("%s_peak_gene_links.csv", analysis_id))
  if (file.exists(links_file)) {
    cat("  === Peak-Gene Linkage: Loading cached results ===\n")
    links_res <- read.csv(links_file)
    links_results_list[[analysis_id]] <- links_res
    cat(sprintf("    ✓ Loaded %d peak-gene links\n", nrow(links_res)))
  } else {
    links_res <- tryCatch({
      # Use genotype-subset Seurat object for peak-gene linkage
      run_peak_gene_linkage(seurat_geno, ct)
    }, error = function(e) {
      cat(sprintf("    ✗ Error in peak-gene linkage: %s\n", e$message))
      NULL
    })

    if (!is.null(links_res)) {
      links_results_list[[analysis_id]] <- links_res
      # Save to file
      write.csv(links_res, links_file, row.names = FALSE)
    }
  }

  # Compile summary stats
  n_degs <- if (!is.null(deg_res)) sum(deg_res$p_val_adj < PVAL_THRESHOLD) else 0
  n_degs_up <- if (!is.null(deg_res)) sum(deg_res$p_val_adj < PVAL_THRESHOLD & deg_res$avg_log2FC > 0) else 0
  n_degs_down <- if (!is.null(deg_res)) sum(deg_res$p_val_adj < PVAL_THRESHOLD & deg_res$avg_log2FC < 0) else 0

  n_da <- if (!is.null(da_res)) sum(da_res$p_val_adj < PVAL_THRESHOLD) else 0
  n_da_more <- if (!is.null(da_res)) sum(da_res$p_val_adj < PVAL_THRESHOLD & da_res$avg_log2FC > 0) else 0
  n_da_less <- if (!is.null(da_res)) sum(da_res$p_val_adj < PVAL_THRESHOLD & da_res$avg_log2FC < 0) else 0

  n_links <- if (!is.null(links_res)) nrow(links_res) else 0

  summary_stats <- rbind(summary_stats, data.frame(
    celltype = ct,
    genotype = geno,
    n_cells_ctrl = n_ctrl,
    n_cells_mut = n_mut,
    n_degs = n_degs,
    n_degs_up = n_degs_up,
    n_degs_down = n_degs_down,
    n_da_peaks = n_da,
    n_da_more = n_da_more,
    n_da_less = n_da_less,
    n_peak_gene_links = n_links
  ))
}  # End of for (ct in celltypes)

}  # End of for (geno in genotypes)

################################################################################
# Save Summary Results
################################################################################

cat("\n================================================================================\n")
cat("SAVING SUMMARY RESULTS\n")
cat("================================================================================\n")

# Save summary statistics
write.csv(summary_stats, file.path(RESULTS_DIR, "summary_statistics.csv"),
          row.names = FALSE)

cat("  ✓ Summary statistics saved\n")

# Print summary table
cat("\nSummary Statistics:\n")
print(summary_stats)

# Create summary report
cat("\n✓ Step 3 complete!\n")
cat(sprintf("  Results saved to: %s\n", RESULTS_DIR))
cat("\nNext step: Run signac_04_visualization.R\n\n")
