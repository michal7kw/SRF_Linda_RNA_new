#!/usr/bin/env Rscript

################################################################################
# Create Coverage Plots with Literature CREs
#
# Purpose:
# - Visualize GABA DEGs with their literature-validated CREs
# - Show ATAC coverage tracks (Combined Ctrl vs Nestin Mut vs Emx1 Mut)
# - Display genomic links from CREs to gene TSS
# - Include gene expression and annotations
#
# Prerequisites:
# - Run analyze_GABA_DEGs_with_literature_CREs.py first
# - Run extract_bigwig_signals_and_visualize.py first
#
# Input:
# - CREs_from_literature/output/GABA_DEG_analysis/gene_summary_{genotype}.tsv
# - CREs_from_literature/output/GABA_DEG_analysis/DEG_CRE_links_{genotype}.tsv
# - signac_results_L1/integrated_seurat_processed.rds
#
# Output:
# - CREs_from_literature/output/GABA_DEG_analysis/coverage_plots_literature_CREs/{gene}_{genotype}_coverage.pdf
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
})

# Configuration
BASE_DIR <- "CREs_from_literature"
INPUT_DIR <- file.path(BASE_DIR, "output", "GABA_DEG_analysis")
OUTPUT_DIR <- file.path(INPUT_DIR, "coverage_plots_literature_CREs")
SEURAT_FILE <- "signac_results_L1/integrated_seurat_processed.rds"

# Analysis parameters
GENOTYPES <- c("Nestin", "Emx1")
CELL_TYPE <- "GABA"
CELLTYPE_COL <- "cell_type_L1"
CONDITION_COL <- "condition"
GENOTYPE_COL <- "genotype"

# Plotting parameters
EXTEND_UPSTREAM <- 50000      # 50kb upstream
EXTEND_DOWNSTREAM <- 50000    # 50kb downstream
TOP_N_GENES <- 10             # Number of top genes to plot per genotype

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("================================================================================\n")
cat("CREATE COVERAGE PLOTS WITH LITERATURE CREs\n")
cat("================================================================================\n\n")

################################################################################
# Load Data
################################################################################

cat("Loading Seurat object...\n")
if (!file.exists(SEURAT_FILE)) {
  stop(sprintf("Error: Seurat file not found: %s\n", SEURAT_FILE))
}
seurat_obj <- readRDS(SEURAT_FILE)
cat(sprintf("  ✓ Loaded Seurat object with %d cells\n", ncol(seurat_obj)))

# Subset to GABA cells
cat(sprintf("\nSubsetting to %s cells...\n", CELL_TYPE))
if (!CELLTYPE_COL %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Error: Cell type column '%s' not found in metadata\n", CELLTYPE_COL))
}

seurat_gaba <- subset(seurat_obj, subset = !!sym(CELLTYPE_COL) == CELL_TYPE)
cat(sprintf("  ✓ %d %s cells\n", ncol(seurat_gaba), CELL_TYPE))

# Check required metadata columns
required_cols <- c(CONDITION_COL, GENOTYPE_COL)
missing_cols <- setdiff(required_cols, colnames(seurat_gaba@meta.data))
if (length(missing_cols) > 0) {
  stop(sprintf("Error: Required metadata columns missing: %s\n",
               paste(missing_cols, collapse = ", ")))
}

# Create combined grouping variable
# Combined_Ctrl: merge Nestin-Ctrl + Emx1-Ctrl
# Nestin_Mut: Nestin-Mut only
# Emx1_Mut: Emx1-Mut only
seurat_gaba@meta.data$combined_group <- ifelse(
  seurat_gaba@meta.data[[CONDITION_COL]] == "Ctrl",
  "Combined_Ctrl",
  paste0(seurat_gaba@meta.data[[GENOTYPE_COL]], "_Mut")
)

cat("\nCombined group distribution:\n")
print(table(seurat_gaba@meta.data$combined_group))

# Set identity
Idents(seurat_gaba) <- seurat_gaba@meta.data$combined_group

################################################################################
# Load Literature CRE Analysis Results
################################################################################

cat("\n================================================================================\n")
cat("LOADING LITERATURE CRE ANALYSIS RESULTS\n")
cat("================================================================================\n\n")

deg_cre_data <- list()
gene_summary_data <- list()

for (genotype in GENOTYPES) {
  cat(sprintf("%s:\n", genotype))

  # Load DEG-CRE links
  deg_cre_file <- file.path(INPUT_DIR, sprintf("DEG_CRE_links_%s.tsv", genotype))
  if (!file.exists(deg_cre_file)) {
    cat(sprintf("  ✗ File not found: %s\n", basename(deg_cre_file)))
    next
  }

  deg_cre <- read.delim(deg_cre_file)
  deg_cre_data[[genotype]] <- deg_cre
  cat(sprintf("  ✓ Loaded %d gene-CRE pairs\n", nrow(deg_cre)))

  # Load gene summary
  summary_file <- file.path(INPUT_DIR, sprintf("gene_summary_%s.tsv", genotype))
  if (file.exists(summary_file)) {
    gene_summary <- read.delim(summary_file)
    gene_summary_data[[genotype]] <- gene_summary
    cat(sprintf("  ✓ Loaded summary for %d genes\n", nrow(gene_summary)))
  }
}

if (length(deg_cre_data) == 0) {
  stop("Error: No DEG-CRE data loaded. Run Python analysis scripts first!")
}

################################################################################
# Select Top Genes for Visualization
################################################################################

cat("\n================================================================================\n")
cat("SELECTING TOP GENES FOR VISUALIZATION\n")
cat("================================================================================\n\n")

genes_to_plot <- list()

for (genotype in GENOTYPES) {
  if (!genotype %in% names(gene_summary_data)) next

  cat(sprintf("%s:\n", genotype))

  gene_summary <- gene_summary_data[[genotype]]

  # Sort by significance (padj) and select top N
  top_genes <- gene_summary %>%
    arrange(gene_padj) %>%
    filter(n_cres > 0) %>%  # Must have at least 1 CRE
    head(TOP_N_GENES)

  genes_to_plot[[genotype]] <- top_genes$gene

  cat(sprintf("  Selected %d genes:\n", length(genes_to_plot[[genotype]])))
  for (i in seq_along(genes_to_plot[[genotype]])) {
    gene <- genes_to_plot[[genotype]][i]
    gene_info <- top_genes[i, ]
    cat(sprintf("    %d. %s: log2FC=%.3f, padj=%.2e, %d CREs\n",
                i, gene, gene_info$gene_log2fc, gene_info$gene_padj, gene_info$n_cres))
  }
}

################################################################################
# Create Coverage Plots
################################################################################

cat("\n================================================================================\n")
cat("CREATING COVERAGE PLOTS\n")
cat("================================================================================\n\n")

# Initialize summary
plot_summary <- data.frame(
  genotype = character(),
  gene = character(),
  n_cres = integer(),
  plot_status = character(),
  stringsAsFactors = FALSE
)

for (genotype in GENOTYPES) {
  if (!genotype %in% names(genes_to_plot)) next

  cat(sprintf("\n%s:\n", genotype))
  cat(sprintf("%s\n", paste(rep("-", nchar(genotype) + 1), collapse = "")))

  genes <- genes_to_plot[[genotype]]
  deg_cre <- deg_cre_data[[genotype]]

  for (gene in genes) {
    cat(sprintf("\nProcessing: %s\n", gene))

    # Check if gene exists in ATAC annotation
    gene_in_atac <- tryCatch({
      coords <- LookupGeneCoords(seurat_gaba, gene = gene, assay = "ATAC")
      !is.null(coords) && !any(is.na(coords)) && length(coords) > 0
    }, error = function(e) FALSE)

    if (!gene_in_atac) {
      cat(sprintf("  ✗ Gene not found in ATAC annotation\n"))
      plot_summary <- rbind(plot_summary, data.frame(
        genotype = genotype,
        gene = gene,
        n_cres = 0,
        plot_status = "missing_annotation"
      ))
      next
    }

    # Get CREs for this gene
    gene_cres <- deg_cre %>% filter(gene == !!gene)

    if (nrow(gene_cres) == 0) {
      cat(sprintf("  ✗ No CREs found\n"))
      plot_summary <- rbind(plot_summary, data.frame(
        genotype = genotype,
        gene = gene,
        n_cres = 0,
        plot_status = "no_cres"
      ))
      next
    }

    cat(sprintf("  Found %d CREs\n", nrow(gene_cres)))

    # Create coverage plot with genomic links
    plot_success <- tryCatch({
      DefaultAssay(seurat_gaba) <- "ATAC"

      # Parse CRE coordinates and create GRanges
      cre_gr <- GRanges(
        seqnames = gene_cres$cre_chr,
        ranges = IRanges(
          start = as.numeric(gene_cres$cre_start),
          end = as.numeric(gene_cres$cre_end)
        ),
        cre_id = gene_cres$cre_id,
        pcc = gene_cres$cre_pcc
      )

      # Get gene coordinates
      gene_coords <- LookupGeneCoords(seurat_gaba, gene = gene, assay = "ATAC")

      # Create Link objects (arcs from CREs to gene TSS)
      link_list <- lapply(seq_along(cre_gr), function(i) {
        cre_start <- start(cre_gr[i])
        cre_end <- end(cre_gr[i])
        gene_start <- start(gene_coords)

        # Link from CRE center to gene TSS
        cre_center <- (cre_start + cre_end) / 2

        link_gr <- GRanges(
          seqnames = as.character(seqnames(cre_gr[i])),
          ranges = IRanges(
            start = min(cre_center, gene_start),
            end = max(cre_center, gene_start)
          ),
          peak1 = cre_start,
          peak2 = gene_start,
          score = abs(cre_gr[i]$pcc),  # Use correlation as score
          cre_id = cre_gr[i]$cre_id
        )
        return(link_gr)
      })

      all_link_gr <- do.call(c, link_list)

      # Add links to ATAC assay
      Links(seurat_gaba[["ATAC"]]) <- all_link_gr
      cat(sprintf("  Added %d genomic links\n", length(all_link_gr)))

      # Create coverage plot
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

      # Add title with gene info
      gene_info <- gene_summary_data[[genotype]] %>% filter(gene == !!gene)
      title_text <- sprintf(
        "%s (%s) - log2FC=%.3f, padj=%.2e\n%d Literature CREs (Hippocampal/GABAergic)",
        gene, genotype,
        gene_info$gene_log2fc[1],
        gene_info$gene_padj[1],
        nrow(gene_cres)
      )

      p_final <- p_cov +
        plot_annotation(
          title = title_text,
          theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
        )

      # Save plot
      plot_file <- file.path(OUTPUT_DIR, sprintf("%s_%s_coverage.pdf", gene, genotype))
      ggsave(plot_file, plot = p_final, width = 16, height = 10, device = "pdf")
      cat(sprintf("  ✓ Saved: %s\n", basename(plot_file)))

      TRUE
    }, error = function(e) {
      cat(sprintf("  ✗ Error: %s\n", e$message))
      FALSE
    })

    # Update summary
    plot_summary <- rbind(plot_summary, data.frame(
      genotype = genotype,
      gene = gene,
      n_cres = nrow(gene_cres),
      plot_status = ifelse(plot_success, "success", "error")
    ))
  }
}

################################################################################
# Save Summary
################################################################################

cat("\n================================================================================\n")
cat("SAVING SUMMARY\n")
cat("================================================================================\n\n")

summary_file <- file.path(OUTPUT_DIR, "coverage_plots_summary.tsv")
write.table(plot_summary, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("✓ Saved summary: %s\n", basename(summary_file)))

# Print summary statistics
cat("\nSummary Statistics:\n")
print(table(plot_summary$plot_status))

cat("\nSuccessful plots by genotype:\n")
success_summary <- plot_summary %>%
  filter(plot_status == "success") %>%
  group_by(genotype) %>%
  summarise(n_genes = n(), .groups = "drop")
print(success_summary)

################################################################################
# Final Summary
################################################################################

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("================================================================================\n\n")

n_success <- sum(plot_summary$plot_status == "success")
n_total <- nrow(plot_summary)

cat(sprintf("Coverage plots created: %d / %d\n", n_success, n_total))
cat(sprintf("\nOutput directory: %s\n", OUTPUT_DIR))
cat("  - {gene}_{genotype}_coverage.pdf - Individual coverage plots\n")
cat("  - coverage_plots_summary.tsv - Summary table\n\n")

if (n_success > 0) {
  cat("Plot interpretation:\n")
  cat("  - Top track: Gene expression (RNA) in Combined_Ctrl, Nestin_Mut, Emx1_Mut\n")
  cat("  - Middle track: ATAC coverage in three groups\n")
  cat("  - Arcs: Genomic links from literature CREs to gene TSS\n")
  cat("    (Arc height = correlation strength)\n")
  cat("  - Bottom: Gene annotation\n\n")
  cat("Look for:\n")
  cat("  - Decreased ATAC signal at CREs in Mut samples\n")
  cat("  - Correspondence between CRE accessibility and gene expression\n")
  cat("  - Multiple CREs showing consistent patterns\n")
}

cat("\nCompleted!\n")
cat("================================================================================\n")
