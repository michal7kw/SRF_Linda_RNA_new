#!/usr/bin/env Rscript

################################################################################
# Peak-Gene Linkage Analysis for GABA Neurons (UPDATED VERSION - Using Pre-computed Links)
#
# This script analyzes peak-gene linkages for specific genes of interest in
# GABA neurons, using PRE-COMPUTED LinkPeaks results from array job scripts.
#
# KEY CHANGE: Instead of manually computing correlations with cor.test(),
# this script loads the LinkPeaks results already computed by array jobs,
# ensuring consistency and avoiding redundant computation.
#
# PREREQUISITES: Run peak-gene linkage array jobs first:
#   sbatch run_peak_gene_linkage_high_priority_L1.sh  OR
#   sbatch run_peak_gene_linkage_all_L1.sh            OR
#   sbatch run_signac_peak_gene_linkage_single_L1.sh GABA Nestin
#   sbatch run_signac_peak_gene_linkage_single_L1.sh GABA Emx1
#
# Updates in this version:
# - Loads pre-computed LinkPeaks from array job outputs
# - Filters to genes of interest from genes_inter.txt
# - Creates coverage plots with genomic link visualization
# - Creates fragment abundance plots
# - Ensures consistency with array job parameters
#
# Analysis approach:
# 1. Load integrated Seurat object with RNA + ATAC modalities
# 2. Load pre-computed LinkPeaks CSV files from Step 3 for GABA neurons
# 3. Filter to genes in genes_inter.txt
# 4. Merge links from Nestin and Emx1 genotypes
# 5. For each gene with significant links:
#    a. Filter mitochondrial genes and validate ATAC annotation
#    b. Add links to Seurat object
#    c. Create coverage plots with linked peaks highlighted
#    d. Create per-cell fragment abundance plots
#
# Input files:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/integrated_seurat_processed.rds
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/genes_inter.txt
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/peak_gene_links/GABA_Nestin_peak_gene_links.csv
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/celltype_results/peak_gene_links/GABA_Emx1_peak_gene_links.csv
#
# Output:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/peak_gene_linkage_analysis_UPDATED/GABA/peak_gene_links/{gene}_peak_gene_links.csv
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/peak_gene_linkage_analysis_UPDATED/GABA/coverage_plots/{gene}_coverage_with_links.pdf
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/peak_gene_linkage_analysis_UPDATED/GABA/fragment_abundance/GABA_fragment_abundance.pdf
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/peak_gene_linkage_analysis_UPDATED/GABA/analysis_summary.csv
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
  library(ggforce)  # Required for genomic link visualization
})

# Configuration
RESULTS_DIR <- "signac_results_L1"
SEURAT_FILE <- file.path(RESULTS_DIR, "integrated_seurat_processed.rds")
GENE_LIST_FILE <- "../genes_inter.txt"
OUTPUT_DIR <- file.path(RESULTS_DIR, "peak_gene_linkage_analysis_UPDATED", "GABA")

# Step 3 pre-computed links directory
STEP3_LINKS_DIR <- file.path(RESULTS_DIR, "celltype_results", "peak_gene_links")

# Cell type and condition columns
CELLTYPE_COL <- "cell_type_L1"
CELL_TYPE <- "GABA"
CONDITION_COL <- "condition"
GENOTYPE_COL <- "genotype"

# Analysis parameters
MIN_CORRELATION <- 0.05       # Minimum correlation for linkage (matches Step 3 threshold)
PVALUE_CUTOFF <- 0.05        # P-value cutoff for significance

# Plotting parameters
EXTEND_UPSTREAM <- 50000      # 50kb upstream for coverage plots
EXTEND_DOWNSTREAM <- 50000    # 50kb downstream for coverage plots

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "peak_gene_links"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "coverage_plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "fragment_abundance"), recursive = TRUE, showWarnings = FALSE)

cat("================================================================================\n")
cat("PEAK-GENE LINKAGE ANALYSIS: GABA NEURONS (Using Pre-computed Links from Step 3)\n")
cat("================================================================================\n\n")
cat("Key changes in this version:\n")
cat("  - Loads pre-computed LinkPeaks from Step 3 (NO manual cor.test())\n")
cat("  - Ensures consistency with Step 3 parameters\n")
cat("  - Faster execution (no redundant correlation computation)\n")
cat("  - Filters to genes of interest from gene list\n\n")

################################################################################
# Load Data
################################################################################

cat("Loading integrated Seurat object...\n")
if (!file.exists(SEURAT_FILE)) {
  stop(sprintf("Error: Seurat file not found: %s\n", SEURAT_FILE))
}
seurat_obj <- readRDS(SEURAT_FILE)
cat(sprintf("  ✓ Loaded Seurat object with %d cells\n", ncol(seurat_obj)))

cat("\nLoading gene list...\n")
if (!file.exists(GENE_LIST_FILE)) {
  stop(sprintf("Error: Gene list file not found: %s\n", GENE_LIST_FILE))
}
genes_of_interest <- read.table(GENE_LIST_FILE, header = FALSE, stringsAsFactors = FALSE)$V1
genes_of_interest <- unique(genes_of_interest)
cat(sprintf("  ✓ Loaded %d genes of interest\n", length(genes_of_interest)))
cat(sprintf("  Genes: %s\n", paste(head(genes_of_interest, 10), collapse = ", ")))
if (length(genes_of_interest) > 10) {
  cat(sprintf("  ... and %d more\n", length(genes_of_interest) - 10))
}

# Filter out mitochondrial genes
mt_genes <- genes_of_interest[grepl("^mt-", genes_of_interest, ignore.case = TRUE)]
if (length(mt_genes) > 0) {
  cat(sprintf("\n  ⚠ Filtering %d mitochondrial genes (incompatible with coverage plots):\n", length(mt_genes)))
  cat(sprintf("    %s\n", paste(mt_genes, collapse = ", ")))
  genes_of_interest <- genes_of_interest[!grepl("^mt-", genes_of_interest, ignore.case = TRUE)]
  cat(sprintf("  ✓ %d genes remaining after filtering\n", length(genes_of_interest)))
}

################################################################################
# Load Pre-computed LinkPeaks from Step 3
################################################################################

cat("\n================================================================================\n")
cat("LOADING PRE-COMPUTED LINKPEAKS FROM STEP 3\n")
cat("================================================================================\n\n")

# Define expected link files
link_files <- c(
  nestin = file.path(STEP3_LINKS_DIR, "GABA_Nestin_peak_gene_links.csv"),
  emx1 = file.path(STEP3_LINKS_DIR, "GABA_Emx1_peak_gene_links.csv")
)

# Check if files exist
missing_files <- !file.exists(link_files)
if (any(missing_files)) {
  cat("✗ ERROR: Pre-computed LinkPeaks files not found!\n\n")
  cat("Missing files:\n")
  for (file in link_files[missing_files]) {
    cat(sprintf("  - %s\n", file))
  }
  cat("\nThese files should be generated by peak-gene linkage array job scripts.\n")
  cat("Please run array jobs first before running this script.\n\n")
  cat("To run array jobs:\n")
  cat("  Option 1 (High Priority): sbatch run_peak_gene_linkage_high_priority_L1.sh\n")
  cat("  Option 2 (All Cell Types): sbatch run_peak_gene_linkage_all_L1.sh\n")
  cat("  Option 3 (Individual): sbatch run_signac_peak_gene_linkage_single_L1.sh GABA Nestin\n")
  cat("                          sbatch run_signac_peak_gene_linkage_single_L1.sh GABA Emx1\n\n")
  stop("Pre-computed LinkPeaks not found. Run array jobs first.")
}

# Load link files
all_links_list <- list()
for (geno_name in names(link_files)) {
  cat(sprintf("Loading %s links...\n", geno_name))
  links_df <- read.csv(link_files[[geno_name]])

  # Check if gene column exists
  if (!"gene" %in% colnames(links_df)) {
    cat(sprintf("  ⚠ Warning: 'gene' column not found in %s\n", geno_name))
    cat(sprintf("  Available columns: %s\n", paste(colnames(links_df), collapse = ", ")))
    next
  }

  # Ensure genotype column exists with capitalized names
  # Array job outputs already include genotype column with proper capitalization
  if (!"genotype" %in% colnames(links_df)) {
    # If missing, add it with capitalized name
    links_df$genotype <- tools::toTitleCase(geno_name)
  }
  all_links_list[[geno_name]] <- links_df

  cat(sprintf("  ✓ Loaded %d links\n", nrow(links_df)))
  cat(sprintf("    Unique genes: %d\n", length(unique(links_df$gene))))
}

# Merge all links
all_links <- do.call(rbind, all_links_list)
cat(sprintf("\n✓ Total links loaded: %d\n", nrow(all_links)))
cat(sprintf("  Unique genes with links: %d\n", length(unique(all_links$gene))))

# Filter to genes of interest
links_filtered <- all_links %>%
  filter(gene %in% genes_of_interest)

cat(sprintf("\n✓ Filtered to genes of interest: %d links\n", nrow(links_filtered)))
cat(sprintf("  Genes with links: %d / %d\n",
            length(unique(links_filtered$gene)),
            length(genes_of_interest)))

# Check if any genes have links
if (nrow(links_filtered) == 0) {
  cat("\n⚠ WARNING: No peak-gene links found for any genes in the gene list!\n")
  cat("This could mean:\n")
  cat("  1. Step 3 correlation threshold was too stringent (try lowering MIN_CORRELATION in Step 3)\n")
  cat("  2. These genes are not accessible in GABA neurons\n")
  cat("  3. Gene names don't match (check capitalization)\n\n")
  stop("No links found for genes of interest")
}

# Show summary by gene
genes_with_links <- links_filtered %>%
  group_by(gene) %>%
  summarise(
    n_links = n(),
    genotypes = paste(unique(genotype), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_links))

cat("\nGenes with most links:\n")
print(head(genes_with_links, 10))

################################################################################
# Subset to GABA Neurons
################################################################################

cat("\n================================================================================\n")
cat("SUBSETTING TO GABA NEURONS\n")
cat("================================================================================\n\n")

# Get GABA cells
gaba_cells <- seurat_obj@meta.data %>%
  filter(.data[[CELLTYPE_COL]] == CELL_TYPE) %>%
  rownames()

if (length(gaba_cells) == 0) {
  stop(sprintf("Error: No cells found for cell type '%s'\n", CELL_TYPE))
}

seurat_gaba <- subset(seurat_obj, cells = gaba_cells)
cat(sprintf("✓ Selected %d GABA neurons\n", ncol(seurat_gaba)))

# Check sample distribution
sample_dist <- table(seurat_gaba@meta.data[[GENOTYPE_COL]],
                     seurat_gaba@meta.data[[CONDITION_COL]])
cat("\nSample distribution:\n")
print(sample_dist)

################################################################################
# Create Combined Grouping Variable
################################################################################

cat("\n================================================================================\n")
cat("CREATING COMBINED GROUPING VARIABLE\n")
cat("================================================================================\n\n")

# Combined_Ctrl: merge Nestin-Ctrl + Emx1-Ctrl
# Nestin_Mut: Nestin-Mut only
# Emx1_Mut: Emx1-Mut only
seurat_gaba@meta.data$combined_group <- ifelse(
  seurat_gaba@meta.data[[CONDITION_COL]] == "Ctrl",
  "Combined_Ctrl",
  paste0(seurat_gaba@meta.data[[GENOTYPE_COL]], "_Mut")
)

# Check distribution
group_dist <- table(seurat_gaba@meta.data$combined_group)
cat("Combined group distribution:\n")
print(group_dist)

# Set identity to combined_group for plotting
Idents(seurat_gaba) <- seurat_gaba@meta.data$combined_group

################################################################################
# Fragment Abundance Analysis
################################################################################

cat("\n================================================================================\n")
cat("CREATING FRAGMENT ABUNDANCE PLOTS\n")
cat("================================================================================\n\n")

if ("nCount_ATAC" %in% colnames(seurat_gaba@meta.data)) {
  plot_data <- seurat_gaba@meta.data %>%
    select(all_of(c(GENOTYPE_COL, CONDITION_COL, "combined_group", "nCount_ATAC"))) %>%
    mutate(log10_fragments = log10(nCount_ATAC + 1))

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
      title = "ATAC Fragment Abundance in GABA Neurons",
      subtitle = "By Combined Group",
      x = "Group",
      y = "log10(Fragment Count + 1)",
      fill = "Group"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Violin plot by genotype
  p2 <- ggplot(plot_data, aes(x = !!sym(CONDITION_COL), y = log10_fragments, fill = !!sym(CONDITION_COL))) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
    facet_wrap(as.formula(paste("~", GENOTYPE_COL)), nrow = 1) +
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
  pdf_file <- file.path(OUTPUT_DIR, "fragment_abundance", "GABA_fragment_abundance.pdf")
  ggsave(pdf_file, plot = p1 + p2, width = 14, height = 6, device = "pdf")
  cat(sprintf("✓ Saved: %s\n", basename(pdf_file)))
} else {
  cat("⚠ Warning: nCount_ATAC not found in metadata\n")
}

################################################################################
# Create Visualizations for Each Gene
################################################################################

cat("\n================================================================================\n")
cat("CREATING GENE-SPECIFIC VISUALIZATIONS\n")
cat("================================================================================\n\n")

# Initialize summary stats
summary_stats <- data.frame(
  gene = character(),
  status = character(),
  n_links_nestin = integer(),
  n_links_emx1 = integer(),
  n_links_total = integer(),
  max_correlation = numeric(),
  stringsAsFactors = FALSE
)

# Get unique genes to process
genes_to_process <- unique(links_filtered$gene)
cat(sprintf("Processing %d genes with links...\n\n", length(genes_to_process)))

for (gene in genes_to_process) {
  cat(sprintf("--- Analyzing gene: %s ---\n", gene))

  # Check if gene exists in ATAC annotation
  gene_in_atac <- tryCatch({
    coords <- LookupGeneCoords(seurat_gaba, gene = gene, assay = "ATAC")
    !is.null(coords) && !any(is.na(coords)) && length(coords) > 0
  }, error = function(e) FALSE)

  if (!gene_in_atac) {
    cat(sprintf("  ⚠ Skipping - gene not found in ATAC annotation\n"))
    summary_stats <- rbind(summary_stats, data.frame(
      gene = gene,
      status = "missing_annotation",
      n_links_nestin = 0,
      n_links_emx1 = 0,
      n_links_total = 0,
      max_correlation = NA
    ))
    next
  }

  # Get links for this gene
  gene_links <- links_filtered %>% filter(gene == !!gene)

  n_links_nestin <- sum(gene_links$genotype == "Nestin")
  n_links_emx1 <- sum(gene_links$genotype == "Emx1")

  cat(sprintf("  Found %d links (Nestin: %d, Emx1: %d)\n",
              nrow(gene_links), n_links_nestin, n_links_emx1))

  # Get max correlation
  max_cor <- max(abs(gene_links$score), na.rm = TRUE)
  cat(sprintf("  Max |correlation|: %.3f\n", max_cor))

  # Save per-gene CSV
  gene_csv <- file.path(OUTPUT_DIR, "peak_gene_links", sprintf("%s_peak_gene_links.csv", gene))
  write.csv(gene_links, gene_csv, row.names = FALSE)

  # Create coverage plot with genomic links
  cat("  Creating coverage plot with genomic links...\n")
  plot_success <- tryCatch({
    DefaultAssay(seurat_gaba) <- "ATAC"

    # Parse peak coordinates from gene_links
    # Expected format: seqnames, start, end columns OR peak column like "chr1-1000-2000"
    if ("seqnames" %in% colnames(gene_links) && "start" %in% colnames(gene_links) && "end" %in% colnames(gene_links)) {
      # GRanges format already available
      peak_gr <- GRanges(
        seqnames = gene_links$seqnames,
        ranges = IRanges(
          start = as.numeric(gene_links$start),
          end = as.numeric(gene_links$end)
        ),
        score = gene_links$score
      )
    } else if ("peak" %in% colnames(gene_links)) {
      # Parse peak string format
      peak_coords <- do.call(rbind, strsplit(as.character(gene_links$peak), "-"))
      peak_gr <- GRanges(
        seqnames = peak_coords[, 1],
        ranges = IRanges(
          start = as.numeric(peak_coords[, 2]),
          end = as.numeric(peak_coords[, 3])
        ),
        score = gene_links$score
      )
    } else {
      cat("    ⚠ Cannot parse peak coordinates from link data\n")
      cat(sprintf("    Available columns: %s\n", paste(colnames(gene_links), collapse = ", ")))
      return(FALSE)
    }

    # Get gene coordinates
    gene_coords <- LookupGeneCoords(seurat_gaba, gene = gene, assay = "ATAC")

    # Create Link objects for each peak-gene pair
    link_list <- lapply(seq_along(peak_gr), function(i) {
      # Create a GRanges spanning from peak to gene TSS
      peak_start <- start(peak_gr[i])
      peak_end <- end(peak_gr[i])
      gene_start <- start(gene_coords)

      # Link connects peak center to gene TSS
      link_center <- (peak_start + peak_end) / 2

      link_gr <- GRanges(
        seqnames = as.character(seqnames(peak_gr[i])),
        ranges = IRanges(
          start = min(link_center, gene_start),
          end = max(link_center, gene_start)
        ),
        peak1 = peak_start,
        peak2 = gene_start,
        score = peak_gr[i]$score
      )
      return(link_gr)
    })

    # Combine into GRanges list
    all_link_gr <- do.call(c, link_list)

    # Add links to ATAC assay
    Links(seurat_gaba[["ATAC"]]) <- all_link_gr
    cat(sprintf("    Added %d genomic links to ATAC assay\n", length(all_link_gr)))

    # Create coverage plot with links
    p_cov <- CoveragePlot(
      object = seurat_gaba,
      region = gene,
      features = gene,
      expression.assay = "RNA",
      assay = "ATAC",
      group.by = "combined_group",
      extend.upstream = EXTEND_UPSTREAM,
      extend.downstream = EXTEND_DOWNSTREAM,
      annotation = TRUE,
      peaks = TRUE,
      links = TRUE,  # Show genomic links as arcs
      window = 100
    )

    # Save plot
    plot_file <- file.path(OUTPUT_DIR, "coverage_plots", sprintf("%s_coverage_with_links.pdf", gene))
    ggsave(plot_file, plot = p_cov, width = 14, height = 10, device = "pdf")
    cat(sprintf("    ✓ Saved: %s\n", basename(plot_file)))

    TRUE
  }, error = function(e) {
    cat(sprintf("    ✗ Error creating coverage plot: %s\n", e$message))
    FALSE
  })

  # Update summary stats
  summary_stats <- rbind(summary_stats, data.frame(
    gene = gene,
    status = ifelse(plot_success, "success", "plot_failed"),
    n_links_nestin = n_links_nestin,
    n_links_emx1 = n_links_emx1,
    n_links_total = nrow(gene_links),
    max_correlation = max_cor
  ))
}

################################################################################
# Save Summary
################################################################################

cat("\n================================================================================\n")
cat("SAVING SUMMARY\n")
cat("================================================================================\n\n")

summary_file <- file.path(OUTPUT_DIR, "analysis_summary.csv")
write.csv(summary_stats, summary_file, row.names = FALSE)
cat(sprintf("✓ Saved summary: %s\n", summary_file))

# Print summary statistics
cat("\nAnalysis Summary:\n")
cat(sprintf("  Total genes processed: %d\n", nrow(summary_stats)))
cat(sprintf("  Successful plots: %d\n", sum(summary_stats$status == "success")))
cat(sprintf("  Failed plots: %d\n", sum(summary_stats$status == "plot_failed")))
cat(sprintf("  Missing annotation: %d\n", sum(summary_stats$status == "missing_annotation")))

cat("\n================================================================================\n")
cat("✓ ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")
cat(sprintf("Results saved to: %s\n", OUTPUT_DIR))
