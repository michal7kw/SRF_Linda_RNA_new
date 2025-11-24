#!/usr/bin/env Rscript

#
# Input Files:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results_L1/integrated_seurat_processed.rds
#     Integrated Seurat object with processed RNA and ATAC data
#
# Output Files:
#   - /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/test_links_plot.pdf
#     Test CoveragePlot with links to verify correct format
#
# Test script to verify correct Link format for Signac CoveragePlot

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

cat("Loading Seurat object...\n")
seurat_obj <- readRDS("signac_results_L1/integrated_seurat_processed.rds")

# Test with a known gene that has links in our computed results
test_gene <- "Eef1a1"  # From the log, this had 1 link
cat(sprintf("Testing with gene: %s\n\n", test_gene))

# Get gene coordinates
gene_coords <- LookupGeneCoords(seurat_obj, gene = test_gene, assay = "ATAC")
cat("Gene coordinates:\n")
print(gene_coords)

# Create a test link manually
# From the log: Eef1a1 had 1 peak-gene link
# Let's create a fake link to see if it displays

# Create a peak region (fake example)
test_peak <- GRanges(
  seqnames = as.character(seqnames(gene_coords)),
  ranges = IRanges(
    start = start(gene_coords) - 50000,
    end = start(gene_coords) - 49000
  )
)

cat("\nTest peak:\n")
print(test_peak)

# Now create the link - METHOD 1: Simple GRanges with metadata
cat("\n=== Testing METHOD 1: GRanges with metadata ===\n")
link_gr1 <- test_peak
mcols(link_gr1)$gene <- test_gene
mcols(link_gr1)$score <- 0.25
mcols(link_gr1)$peak <- sprintf("%s-%d-%d",
                                seqnames(test_peak),
                                start(test_peak),
                                end(test_peak))

cat("Link object (Method 1):\n")
print(link_gr1)

# Add to ATAC assay
cat("\nAdding link to ATAC assay...\n")
Links(seurat_obj[["ATAC"]]) <- link_gr1

# Check if it was added
cat("Links in ATAC assay after adding:\n")
added_links <- Links(seurat_obj[["ATAC"]])
print(added_links)
cat("Number of links:", length(added_links), "\n")

# Try to plot with links
cat("\nCreating CoveragePlot with links=TRUE...\n")
DefaultAssay(seurat_obj) <- "ATAC"

pdf("test_links_plot.pdf", width = 12, height = 8)
tryCatch({
  p <- CoveragePlot(
    object = seurat_obj,
    region = test_gene,
    extend.upstream = 100000,
    extend.downstream = 100000,
    annotation = TRUE,
    peaks = TRUE,
    links = TRUE
  )
  print(p)
  cat("✓ Plot created successfully\n")
}, error = function(e) {
  cat("✗ Error creating plot:", e$message, "\n")
})
dev.off()

cat("\nTest complete. Check test_links_plot.pdf for results.\n")
