#!/usr/bin/env Rscript

################################################################################
# Signac Multiome Analysis - Step 6: Create BigWig Tracks for IGV (L1 VERSION)
################################################################################
#
# This script creates pseudobulk bigwig files from ATAC-seq data for
# visualization in genome browsers like IGV.
#
# L1 VERSION: Uses broad cell lineage annotations (9 categories: Excitatory,
# GABA, Oligo, Astrocytes, Immune, Vascular, Progenitors, Unknown, Ependymal)
#
# Features:
#   - Pseudobulk aggregation by L1 cell lineage and condition
#   - Normalized coverage tracks (CPM, RPKM)
#   - Separate tracks for Ctrl vs Mut comparison
#   - Uses translated fragment files with correct cell barcodes
#   - Higher cell counts per track → better signal-to-noise
#
# Output:
#   - BigWig files organized by L1 lineage and condition
#   - IGV session file for easy loading
#
# Inputs:
#   - Integrated Seurat object from Step 2:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/integrated_seurat_processed.rds
#   - Translated fragment files from Step 1:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/translated_fragments/{sample}_fragments_translated.tsv.gz
#   - Summary statistics from Step 3:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/summary_statistics.csv
#   - Chromosome sizes file:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/references/mm10.chrom.sizes
#
# Outputs:
#   - Per-sample BigWig tracks (all cells):
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/bigwig_tracks_L1/by_sample/{Genotype}-{Condition}_all_cells.bw
#   - Per-celltype BigWig tracks:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/bigwig_tracks_L1/by_celltype/{CellType}_{Genotype}-{Condition}.bw
#   - IGV session file:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/bigwig_tracks_L1/igv_session.xml
#   - README with usage instructions:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/bigwig_tracks_L1/README.md
#
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(rtracklayer)  # For bigwig export
  library(data.table)
  library(dplyr)
  library(parallel)
})

cat("================================================================================\n")
cat("  SIGNAC MULTIOME ANALYSIS - STEP 6: CREATE BIGWIG TRACKS (L1 VERSION)\n")
cat("================================================================================\n\n")

################################################################################
# Configuration
################################################################################

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
INPUT_DIR <- file.path(BASE_DIR, "signac_results_L1")
TRANSLATED_FRAG_DIR <- file.path(BASE_DIR, "signac_results", "translated_fragments")
OUTPUT_DIR <- file.path(INPUT_DIR, "bigwig_tracks_L1")

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "by_celltype"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "by_sample"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "combined"), recursive = TRUE, showWarnings = FALSE)

# Cell type and condition columns
CELLTYPE_COL <- "cell_type_L1"  # L1 version: broad cell lineages
CONDITION_COL <- "atac_sample"

# Sample identifiers (all 4 samples: Nestin + Emx1)
SAMPLES <- c(
  "R26-Nestin-Ctrl-adult",
  "R26-Nestin-Mut-adult",
  "R26-Emx1-Ctrl-adult",
  "R26-Emx1-Mut-adult"
)

# Genome information
GENOME <- "mm10"
CHR_SIZES_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/references/mm10.chrom.sizes"

# Processing parameters
BINSIZE <- 50  # Base pairs for binning (50bp is standard for ATAC-seq)
N_CORES <- 32  # Parallel processing cores

# Normalization method
NORM_METHOD <- "CPM"  # CPM (counts per million) or RPKM

cat("Configuration:\n")
cat(sprintf("  Input directory: %s\n", INPUT_DIR))
cat(sprintf("  Fragment directory: %s\n", TRANSLATED_FRAG_DIR))
cat(sprintf("  Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("  Bin size: %d bp\n", BINSIZE))
cat(sprintf("  Normalization: %s\n", NORM_METHOD))
cat(sprintf("  Cores: %d\n", N_CORES))
cat("\n")

################################################################################
# Function: Get Chromosome Sizes
################################################################################

get_chromosome_sizes <- function(genome = "mm10", sizes_file = NULL) {
  cat("  Loading chromosome sizes...\n")

  if (!is.null(sizes_file) && file.exists(sizes_file)) {
    # Read from file
    chrom_sizes <- read.table(sizes_file, header = FALSE,
                              col.names = c("chrom", "size"))
    cat(sprintf("    ✓ Loaded from file: %d chromosomes\n", nrow(chrom_sizes)))
  } else {
    # Get from BSgenome
    library(BSgenome.Mmusculus.UCSC.mm10)
    genome_obj <- BSgenome.Mmusculus.UCSC.mm10

    # Get standard chromosomes
    standard_chroms <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM")
    all_lengths <- seqlengths(genome_obj)

    chrom_sizes <- data.frame(
      chrom = names(all_lengths[standard_chroms]),
      size = as.numeric(all_lengths[standard_chroms])
    )

    cat(sprintf("    ✓ Loaded from BSgenome: %d chromosomes\n", nrow(chrom_sizes)))
  }

  return(chrom_sizes)
}

################################################################################
# Function: Create Genome Bins
################################################################################

create_genome_bins <- function(chrom_sizes, binsize = 50) {
  cat(sprintf("  Creating genome bins (%d bp)...\n", binsize))

  bins_list <- list()

  for (i in 1:nrow(chrom_sizes)) {
    chr <- chrom_sizes$chrom[i]
    chr_len <- chrom_sizes$size[i]

    # Create bins for this chromosome
    starts <- seq(1, chr_len, by = binsize)
    ends <- pmin(starts + binsize - 1, chr_len)

    bins_list[[chr]] <- GRanges(
      seqnames = chr,
      ranges = IRanges(start = starts, end = ends)
    )
  }

  # Combine all chromosomes
  # Use GRangesList and unlist for safer combining
  genome_bins <- unlist(GRangesList(bins_list))

  # Set seqlengths (required for BigWig export)
  seqlengths_vec <- setNames(chrom_sizes$size, chrom_sizes$chrom)
  seqlengths(genome_bins) <- seqlengths_vec

  cat(sprintf("    ✓ Created %d bins across %d chromosomes\n",
              length(genome_bins), length(bins_list)))

  return(genome_bins)
}

################################################################################
# Function: Load Fragments for Cell Subset
################################################################################

load_fragments_for_cells <- function(fragment_file, cell_barcodes) {
  cat(sprintf("    Loading fragments from: %s\n", basename(fragment_file)))
  cat(sprintf("    Target cells: %d\n", length(cell_barcodes)))

  if (!file.exists(fragment_file)) {
    cat(sprintf("      ✗ Fragment file not found: %s\n", fragment_file))
    return(NULL)
  }

  # Read all fragments (memory permitting - for large files, use filtering)
  # Fragment format: chr, start, end, barcode, count
  fragments <- fread(
    cmd = sprintf("zcat %s", fragment_file),
    header = FALSE,
    col.names = c("chr", "start", "end", "barcode", "count"),
    showProgress = FALSE
  )

  cat(sprintf("      Total fragments read: %d\n", nrow(fragments)))

  # Filter to cells of interest
  fragments_filtered <- fragments[fragments$barcode %in% cell_barcodes, ]

  cat(sprintf("      Fragments for target cells: %d (%.1f%%)\n",
              nrow(fragments_filtered),
              100 * nrow(fragments_filtered) / nrow(fragments)))

  if (nrow(fragments_filtered) == 0) {
    cat("      ⚠ No fragments found for target cells\n")
    return(NULL)
  }

  # Convert to GRanges
  frags_gr <- GRanges(
    seqnames = fragments_filtered$chr,
    ranges = IRanges(
      start = fragments_filtered$start,
      end = fragments_filtered$end
    ),
    count = fragments_filtered$count
  )

  return(frags_gr)
}

################################################################################
# Function: Compute Coverage in Bins
################################################################################

compute_binned_coverage <- function(fragments_gr, genome_bins, norm_method = "CPM") {
  cat("      Computing coverage...\n")

  if (is.null(fragments_gr) || length(fragments_gr) == 0) {
    cat("        ⚠ No fragments to process\n")
    return(NULL)
  }

  # Count overlaps between fragments and bins
  overlaps <- findOverlaps(fragments_gr, genome_bins)

  # Sum counts per bin (handling fragment counts if present)
  if ("count" %in% names(mcols(fragments_gr))) {
    counts_per_bin <- aggregate(
      mcols(fragments_gr)$count[queryHits(overlaps)],
      by = list(bin = subjectHits(overlaps)),
      FUN = sum
    )
  } else {
    # Just count number of fragments per bin
    counts_per_bin <- as.data.frame(table(subjectHits(overlaps)))
    colnames(counts_per_bin) <- c("bin", "x")
    counts_per_bin$bin <- as.integer(as.character(counts_per_bin$bin))
  }

  # Initialize coverage vector
  coverage_vec <- rep(0, length(genome_bins))
  coverage_vec[counts_per_bin$bin] <- counts_per_bin$x

  # Normalize
  if (norm_method == "CPM") {
    # Counts per million reads
    total_counts <- sum(coverage_vec)
    coverage_vec <- coverage_vec * 1e6 / total_counts
    cat(sprintf("        Normalized: CPM (total counts: %.2e)\n", total_counts))
  } else if (norm_method == "RPKM") {
    # Reads per kilobase per million reads
    total_counts <- sum(coverage_vec)
    bin_kb <- width(genome_bins) / 1000
    coverage_vec <- coverage_vec / bin_kb * 1e6 / total_counts
    cat(sprintf("        Normalized: RPKM (total counts: %.2e)\n", total_counts))
  }

  # Add coverage to bins
  mcols(genome_bins)$score <- coverage_vec

  cat(sprintf("        Bins with coverage > 0: %d (%.1f%%)\n",
              sum(coverage_vec > 0),
              100 * sum(coverage_vec > 0) / length(coverage_vec)))

  return(genome_bins)
}

################################################################################
# Function: Export BigWig
################################################################################

export_bigwig <- function(bins_with_coverage, output_file, description = "") {
  cat(sprintf("      Exporting to: %s\n", basename(output_file)))

  if (is.null(bins_with_coverage)) {
    cat("        ⚠ No coverage data to export\n")
    return(FALSE)
  }

  # Remove bins with zero coverage (reduces file size)
  bins_nonzero <- bins_with_coverage[mcols(bins_with_coverage)$score > 0]

  cat(sprintf("        Non-zero bins: %d\n", length(bins_nonzero)))

  if (length(bins_nonzero) == 0) {
    cat("        ⚠ No non-zero coverage bins\n")
    return(FALSE)
  }

  # Export as bigwig
  tryCatch({
    export.bw(bins_nonzero, con = output_file)
    cat(sprintf("        ✓ Exported: %s (%.2f MB)\n",
                basename(output_file),
                file.size(output_file) / 1e6))
    return(TRUE)
  }, error = function(e) {
    cat(sprintf("        ✗ Export failed: %s\n", e$message))
    return(FALSE)
  })
}

################################################################################
# Load Data
################################################################################

cat("\n================================================================================\n")
cat("STEP 1: Loading Data\n")
cat("================================================================================\n\n")

# Load Seurat object
cat("  Loading Seurat object...\n")
seurat_obj <- readRDS(file.path(INPUT_DIR, "integrated_seurat_processed.rds"))
cat(sprintf("    ✓ Loaded: %d cells\n", ncol(seurat_obj)))

# Get metadata
metadata <- seurat_obj@meta.data

# Check for required columns
if (!CELLTYPE_COL %in% colnames(metadata)) {
  stop(sprintf("Cell type column '%s' not found in metadata", CELLTYPE_COL))
}

if (!CONDITION_COL %in% colnames(metadata)) {
  stop(sprintf("Condition column '%s' not found in metadata", CONDITION_COL))
}

# Get cell types
celltypes <- unique(metadata[[CELLTYPE_COL]])
celltypes <- celltypes[!is.na(celltypes)]

cat(sprintf("  Cell types: %d\n", length(celltypes)))
cat(sprintf("    %s\n", paste(head(celltypes, 5), collapse = ", ")))
if (length(celltypes) > 5) {
  cat(sprintf("    ... and %d more\n", length(celltypes) - 5))
}

# Load chromosome sizes
chrom_sizes <- get_chromosome_sizes(genome = GENOME, sizes_file = CHR_SIZES_FILE)

# Create genome bins
genome_bins <- create_genome_bins(chrom_sizes, binsize = BINSIZE)

################################################################################
# Identify Fragment Files
################################################################################

cat("\n================================================================================\n")
cat("STEP 2: Identifying Fragment Files\n")
cat("================================================================================\n\n")

# Look for translated fragment files
translated_files <- list.files(
  TRANSLATED_FRAG_DIR,
  pattern = "_fragments_translated\\.tsv\\.gz$",
  full.names = TRUE
)

if (length(translated_files) == 0) {
  stop(sprintf("No translated fragment files found in: %s\n\nPlease run Step 1 (signac_01_load_and_integrate.R) first to create translated fragment files.", TRANSLATED_FRAG_DIR))
}

cat(sprintf("  Found %d translated fragment files:\n", length(translated_files)))
for (f in translated_files) {
  cat(sprintf("    - %s\n", basename(f)))
}

# Map fragment files to samples
fragment_map <- list()
for (sample in SAMPLES) {
  # Expected fragment file name
  frag_file <- file.path(TRANSLATED_FRAG_DIR,
                         sprintf("%s_fragments_translated.tsv.gz", sample))

  if (file.exists(frag_file)) {
    fragment_map[[sample]] <- frag_file
  } else {
    cat(sprintf("  ⚠ Warning: Fragment file not found for %s\n", sample))
  }
}

cat("\n  Fragment file mapping:\n")
for (sample_name in names(fragment_map)) {
  cat(sprintf("    %s -> %s\n", sample_name, basename(fragment_map[[sample_name]])))
}

################################################################################
# Create BigWig Tracks
################################################################################

cat("\n================================================================================\n")
cat("STEP 3: Creating BigWig Tracks\n")
cat("================================================================================\n\n")

# Summary tracking
bigwig_files <- list()
n_success <- 0
n_failed <- 0

################################################################################
# 3A: Per Sample (All Cells)
################################################################################

cat("\n--- Per Sample Tracks (All Cells) ---\n\n")

for (sample_name in names(fragment_map)) {
  cat(sprintf("=== %s ===\n", sample_name))

  # Get all cells for this sample
  cells_sample <- rownames(metadata[metadata[[CONDITION_COL]] == sample_name, ])

  cat(sprintf("  Total cells: %d\n", length(cells_sample)))

  if (length(cells_sample) == 0) {
    cat("  ⚠ No cells for this sample, skipping\n\n")
    n_failed <- n_failed + 1
    next
  }

  # Strip sample prefix from barcodes to match fragment file format
  # Seurat barcodes: "Sample_BARCODE-1-2"
  # Fragment barcodes: "BARCODE-1-2"
  cells_sample_stripped <- sub("^[^_]+_", "", cells_sample)

  # Load fragments
  fragment_file <- fragment_map[[sample_name]]
  frags <- load_fragments_for_cells(fragment_file, cells_sample_stripped)

  if (is.null(frags)) {
    cat("  ⚠ Failed to load fragments, skipping\n\n")
    n_failed <- n_failed + 1
    next
  }

  # Compute coverage
  bins_cov <- compute_binned_coverage(frags, genome_bins, norm_method = NORM_METHOD)

  # Export
  # Extract genotype and condition from sample name
  # Format: R26-{Genotype}-{Condition}-adult
  sample_parts <- strsplit(sample_name, "-")[[1]]
  genotype <- sample_parts[2]  # "Nestin" or "Emx1"
  condition <- sample_parts[3]  # "Ctrl" or "Mut"
  sample_label <- sprintf("%s-%s", genotype, condition)

  output_file <- file.path(OUTPUT_DIR, "by_sample",
                           sprintf("%s_all_cells.bw", sample_label))

  success <- export_bigwig(bins_cov, output_file,
                          description = sprintf("%s - All Cells", sample_name))

  if (success) {
    bigwig_files[[sprintf("%s_all", sample_label)]] <- output_file
    n_success <- n_success + 1
  } else {
    n_failed <- n_failed + 1
  }

  cat("\n")
}

################################################################################
# 3B: Per Cell Type per Sample
################################################################################

cat("\n--- Per Cell Type Tracks ---\n\n")

# Load summary statistics to prioritize cell types with enough cells
summary_stats <- tryCatch({
  read.csv(file.path(INPUT_DIR, "celltype_results", "summary_statistics.csv"))
}, error = function(e) {
  NULL
})

# Determine which cell types to process
if (!is.null(summary_stats)) {
  # Filter to cell types with enough cells in both conditions
  celltypes_to_process <- summary_stats %>%
    filter(!is.na(n_degs) & n_cells_ctrl >= 50 & n_cells_mut >= 50) %>%
    arrange(desc(n_cells_ctrl + n_cells_mut)) %>%
    pull(celltype)

  cat(sprintf("  Processing %d cell types with ≥50 cells per condition\n\n",
              length(celltypes_to_process)))
} else {
  # Process all cell types
  celltypes_to_process <- celltypes
  cat(sprintf("  Processing all %d cell types\n\n", length(celltypes_to_process)))
}

for (ct in celltypes_to_process) {
  cat(sprintf("=== %s ===\n", ct))

  # Get cells for this cell type
  cells_ct <- rownames(metadata[metadata[[CELLTYPE_COL]] == ct &
                                 !is.na(metadata[[CELLTYPE_COL]]), ])

  if (length(cells_ct) == 0) {
    cat("  ⚠ No cells for this cell type, skipping\n\n")
    next
  }

  # Process each sample
  for (sample_name in names(fragment_map)) {
    # Extract genotype and condition from sample name
    # Format: R26-{Genotype}-{Condition}-adult
    sample_parts <- strsplit(sample_name, "-")[[1]]
    genotype <- sample_parts[2]  # "Nestin" or "Emx1"
    condition <- sample_parts[3]  # "Ctrl" or "Mut"
    sample_label <- sprintf("%s-%s", genotype, condition)

    cat(sprintf("  --- %s ---\n", sample_label))

    # Get cells for this cell type AND sample
    cells_ct_sample <- intersect(
      cells_ct,
      rownames(metadata[metadata[[CONDITION_COL]] == sample_name, ])
    )

    cat(sprintf("    Cells: %d\n", length(cells_ct_sample)))

    if (length(cells_ct_sample) < 10) {
      cat("    ⚠ Too few cells (<10), skipping\n")
      n_failed <- n_failed + 1
      next
    }

    # Strip sample prefix from barcodes to match fragment file format
    # Seurat barcodes: "Sample_BARCODE-1-2"
    # Fragment barcodes: "BARCODE-1-2"
    cells_ct_sample_stripped <- sub("^[^_]+_", "", cells_ct_sample)

    # Load fragments
    fragment_file <- fragment_map[[sample_name]]
    frags <- load_fragments_for_cells(fragment_file, cells_ct_sample_stripped)

    if (is.null(frags)) {
      cat("    ⚠ Failed to load fragments, skipping\n")
      n_failed <- n_failed + 1
      next
    }

    # Compute coverage
    bins_cov <- compute_binned_coverage(frags, genome_bins, norm_method = NORM_METHOD)

    # Export
    ct_safe <- gsub(" ", "_", ct)
    ct_safe <- gsub("/", "_", ct_safe)
    output_file <- file.path(OUTPUT_DIR, "by_celltype",
                            sprintf("%s_%s.bw", ct_safe, sample_label))

    success <- export_bigwig(bins_cov, output_file,
                            description = sprintf("%s - %s", ct, sample_label))

    if (success) {
      bigwig_files[[sprintf("%s_%s", ct_safe, sample_label)]] <- output_file
      n_success <- n_success + 1
    } else {
      n_failed <- n_failed + 1
    }
  }

  cat("\n")
}

################################################################################
# Create IGV Session File
################################################################################

cat("\n================================================================================\n")
cat("STEP 4: Creating IGV Session File\n")
cat("================================================================================\n\n")

if (length(bigwig_files) > 0) {
  cat("  Generating IGV session XML...\n")

  # Create session XML
  session_xml <- c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    '<Session genome="mm10" version="8">',
    '  <Resources>'
  )

  # Add each bigwig file as a resource
  for (name in names(bigwig_files)) {
    file_path <- bigwig_files[[name]]
    session_xml <- c(session_xml,
                    sprintf('    <Resource path="%s"/>', file_path))
  }

  session_xml <- c(session_xml,
                  '  </Resources>',
                  '  <Panel height="800" name="DataPanel" width="1600">',
                  '    <Track altColor="0,0,178" autoScale="true" color="0,0,178" displayMode="COLLAPSED" ',
                  '           featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" ',
                  '           sortable="false" visible="true"/>')

  # Add track definitions
  for (name in names(bigwig_files)) {
    file_path <- bigwig_files[[name]]

    # Determine color based on genotype and condition
    # Nestin: blue/green, Emx1: red/orange
    if (grepl("Nestin.*Ctrl", name)) {
      color <- "0,150,0"    # Green for Nestin-Ctrl
    } else if (grepl("Nestin.*Mut", name)) {
      color <- "0,100,200"  # Blue for Nestin-Mut
    } else if (grepl("Emx1.*Ctrl", name)) {
      color <- "255,140,0"  # Orange for Emx1-Ctrl
    } else if (grepl("Emx1.*Mut", name)) {
      color <- "200,0,0"    # Red for Emx1-Mut
    } else {
      color <- "128,128,128"  # Gray for other
    }

    session_xml <- c(session_xml,
                    sprintf('    <Track altColor="%s" autoScale="true" color="%s" displayMode="COLLAPSED" ',
                            color, color),
                    sprintf('           featureVisibilityWindow="-1" fontSize="10" height="50" id="%s" name="%s" ',
                            name, name),
                    sprintf('           normalize="false" renderer="BAR_CHART" sortable="true" visible="true" windowFunction="mean"/>',
                            name))
  }

  session_xml <- c(session_xml,
                  '  </Panel>',
                  '  <PanelLayout dividerFractions="0.8"/>',
                  '  <HiddenAttributes>',
                  '    <Attribute name="DATA FILE"/>',
                  '    <Attribute name="DATA TYPE"/>',
                  '  </HiddenAttributes>',
                  '</Session>')

  # Write session file
  session_file <- file.path(OUTPUT_DIR, "igv_session.xml")
  writeLines(session_xml, session_file)

  cat(sprintf("  ✓ Session file created: %s\n", session_file))

  # Create README with loading instructions
  readme_text <- sprintf("
# BigWig Tracks for IGV Visualization

## Overview

This directory contains pseudobulk ATAC-seq coverage tracks in BigWig format.

- **Normalization**: %s
- **Bin size**: %d bp
- **Genome**: %s
- **Total tracks**: %d

## Directory Structure

- `by_sample/` - Tracks aggregated across all cells per sample
- `by_celltype/` - Tracks per cell type per sample (Ctrl vs Mut)
- `combined/` - Additional combined tracks (if generated)

## Loading in IGV

### Option 1: Load Session File (Recommended)
1. Open IGV Desktop
2. File → Load from File → Select `igv_session.xml`
3. All tracks will be loaded automatically

### Option 2: Manual Loading
1. Open IGV Desktop
2. Tracks → Load from File
3. Select individual .bw files

## Track Naming Convention

- Per-sample format: `{Genotype}-{Condition}_all_cells.bw`
  - Example: `Nestin-Ctrl_all_cells.bw`, `Emx1-Mut_all_cells.bw`
- Per-cell-type format: `{CellType}_{Genotype}-{Condition}.bw`
  - Example: `Microglia_Nestin-Ctrl.bw`, `Microglia_Emx1-Mut.bw`

## Color Coding in IGV

- **Nestin-Ctrl**: Green (0,150,0)
- **Nestin-Mut**: Blue (0,100,200)
- **Emx1-Ctrl**: Orange (255,140,0)
- **Emx1-Mut**: Red (200,0,0)

## Visualization Tips

1. **Compare Ctrl vs Mut**: Load both tracks for same cell type, align vertically
2. **Adjust scale**: Right-click track → Set Data Range
3. **Gene of interest**: Use search box to navigate to specific genes
4. **Export images**: File → Save Image

## Files Generated

", NORM_METHOD, BINSIZE, GENOME, length(bigwig_files))

  # Add file list
  readme_text <- paste0(readme_text, "\n### Generated BigWig Files:\n\n")
  for (name in sort(names(bigwig_files))) {
    file_path <- bigwig_files[[name]]
    file_size <- file.size(file_path) / 1e6
    readme_text <- paste0(readme_text,
                         sprintf("- `%s` (%.2f MB)\n",
                                basename(file_path), file_size))
  }

  readme_text <- paste0(readme_text, sprintf("
## Analysis Date

Generated: %s

## Questions?

For issues with track visualization or questions about the analysis,
refer to the main pipeline documentation.
", Sys.Date()))

  writeLines(readme_text, file.path(OUTPUT_DIR, "README.md"))

  cat(sprintf("  ✓ README created: %s\n", file.path(OUTPUT_DIR, "README.md")))
}

################################################################################
# Final Summary
################################################################################

cat("\n================================================================================\n")
cat("BIGWIG TRACK GENERATION COMPLETE\n")
cat("================================================================================\n\n")

cat("Summary:\n")
cat(sprintf("  Successfully created: %d tracks\n", n_success))
cat(sprintf("  Failed: %d tracks\n", n_failed))
cat(sprintf("  Total tracks: %d\n", length(bigwig_files)))
cat("\n")

cat("Output directory structure:\n")
cat(sprintf("  %s/\n", OUTPUT_DIR))
cat("    ├── by_sample/          # Per-sample tracks\n")
cat("    ├── by_celltype/        # Cell-type-specific tracks\n")
cat("    ├── igv_session.xml     # IGV session file\n")
cat("    └── README.md           # Usage instructions\n")
cat("\n")

if (length(bigwig_files) > 0) {
  cat("Next steps:\n")
  cat("  1. Open IGV Desktop (https://software.broadinstitute.org/software/igv/)\n")
  cat("  2. Load the session file:\n")
  cat(sprintf("     File → Load from File → %s\n", file.path(OUTPUT_DIR, "igv_session.xml")))
  cat("  3. Navigate to genes of interest (e.g., Srf, Fos, Jun)\n")
  cat("  4. Compare Ctrl vs Mut tracks side-by-side\n")
  cat("\n")

  cat("Example genes to explore:\n")
  cat("  - Srf (chr17:46386621-46397775)\n")
  cat("  - Fos (chr12:85473795-85476830)\n")
  cat("  - Jun (chr4:95057867-95060858)\n")
  cat("\n")
}

cat("✓ Step 6 complete!\n")
cat(sprintf("  All BigWig files saved to: %s\n", OUTPUT_DIR))
cat("\n")
