#!/usr/bin/env Rscript

################################################################################
# Signac Multiome Analysis - Step 2: Processing and Dimensionality Reduction
################################################################################
#
# This script processes RNA and ATAC data, performs dimensionality reduction,
# and creates joint UMAP visualizations.
#
# Inputs:
#   - Integrated Seurat object:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/integrated_seurat_raw.rds
#
# Outputs:
#   - Processed Seurat object:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/integrated_seurat_processed.rds
#   - UMAP overview plots:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/02_umap_overview.pdf
#   - RNA processing QC plots:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/qc/02_rna_processing.pdf
#   - ATAC processing QC plots:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/qc/02_atac_processing.pdf
#   - Processing summary:
#     /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/signac_results/02_processing_summary.txt
#
# Based on: https://stuartlab.org/signac/articles/pbmc_multiomic
#
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

cat("================================================================================\n")
cat("  SIGNAC MULTIOME ANALYSIS - STEP 2: PROCESSING & DIMENSIONALITY REDUCTION\n")
cat("================================================================================\n\n")

################################################################################
# Configuration
################################################################################

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline"
INPUT_DIR <- file.path(BASE_DIR, "signac_results")
OUTPUT_DIR <- INPUT_DIR

# Processing parameters
RNA_N_FEATURES <- 3000  # Number of variable features for RNA
ATAC_MIN_CUTOFF <- 5     # Minimum cutoff for ATAC feature selection
N_PCS <- 50              # Number of PCs for RNA
N_LSIS <- 40             # Number of LSI components for ATAC
UMAP_DIMS_RNA <- 1:50    # Dimensions to use for RNA UMAP
UMAP_DIMS_ATAC <- 2:40   # Dimensions to use for ATAC UMAP (skip first)

################################################################################
# Load Data
################################################################################

cat("\n=== Loading integrated Seurat object ===\n")
input_file <- file.path(INPUT_DIR, "integrated_seurat_raw.rds")

if (!file.exists(input_file)) {
  stop(sprintf("Input file not found: %s\nPlease run signac_01_load_and_integrate.R first.", input_file))
}

seurat_obj <- readRDS(input_file)

cat(sprintf("  Loaded object: %d RNA features, %d ATAC features, %d cells\n",
            nrow(seurat_obj[["RNA"]]),
            nrow(seurat_obj[["ATAC"]]),
            ncol(seurat_obj)))

################################################################################
# RNA Processing
################################################################################

cat("\n================================================================================\n")
cat("STEP 1: RNA Data Processing\n")
cat("================================================================================\n")

DefaultAssay(seurat_obj) <- "RNA"

# Normalize and find variable features
cat("  Normalizing RNA data...\n")
seurat_obj <- NormalizeData(seurat_obj)

cat(sprintf("  Finding top %d variable features...\n", RNA_N_FEATURES))
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = RNA_N_FEATURES)

# Scale data
cat("  Scaling RNA data...\n")
seurat_obj <- ScaleData(seurat_obj)

# PCA
cat(sprintf("  Running PCA (computing %d PCs)...\n", N_PCS))
seurat_obj <- RunPCA(seurat_obj, npcs = N_PCS, verbose = FALSE)

# Visualize PCA
pca_plot <- DimPlot(seurat_obj, reduction = "pca", group.by = "atac_sample") +
  ggtitle("RNA PCA by Sample")

# Check if cell type annotation exists
celltype_col <- NULL
for (col in c("cell_type_L2_new", "cell_type", "celltype", "predicted.id")) {
  if (col %in% colnames(seurat_obj@meta.data)) {
    celltype_col <- col
    cat(sprintf("  Found cell type column: %s\n", celltype_col))
    break
  }
}

if (!is.null(celltype_col)) {
  pca_plot_celltype <- DimPlot(seurat_obj, reduction = "pca", group.by = celltype_col) +
    ggtitle("RNA PCA by Cell Type")
}

# Elbow plot
elbow_plot <- ElbowPlot(seurat_obj, ndims = N_PCS) +
  ggtitle("RNA PCA Elbow Plot")

# UMAP on RNA
cat("  Computing RNA UMAP...\n")
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = UMAP_DIMS_RNA,
                      reduction.name = "umap.rna", reduction.key = "rnaUMAP_")

cat("  ✓ RNA processing complete\n")

################################################################################
# ATAC Processing
################################################################################

cat("\n================================================================================\n")
cat("STEP 2: ATAC Data Processing\n")
cat("================================================================================\n")

DefaultAssay(seurat_obj) <- "ATAC"

# TF-IDF normalization
cat("  Running TF-IDF normalization...\n")
seurat_obj <- RunTFIDF(seurat_obj)

# Feature selection
cat(sprintf("  Finding top features (min.cutoff = %d)...\n", ATAC_MIN_CUTOFF))
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = ATAC_MIN_CUTOFF)

# LSI (SVD on TF-IDF matrix)
cat(sprintf("  Running LSI (computing %d components)...\n", N_LSIS))
seurat_obj <- RunSVD(seurat_obj, n = N_LSIS)

# Check correlation with sequencing depth
cat("  Checking LSI component correlation with sequencing depth...\n")
depth_cor_plot <- DepthCor(seurat_obj) +
  ggtitle("LSI Components vs Sequencing Depth")

# UMAP on ATAC (excluding first LSI component if correlated with depth)
cat("  Computing ATAC UMAP (excluding first LSI component)...\n")
seurat_obj <- RunUMAP(seurat_obj, reduction = "lsi", dims = UMAP_DIMS_ATAC,
                      reduction.name = "umap.atac", reduction.key = "atacUMAP_")

cat("  ✓ ATAC processing complete\n")

################################################################################
# Joint UMAP Visualization
################################################################################

cat("\n================================================================================\n")
cat("STEP 3: Joint Multimodal UMAP\n")
cat("================================================================================\n")

# Build joint neighbor graph
cat("  Building weighted nearest neighbor graph...\n")
seurat_obj <- FindMultiModalNeighbors(
  object = seurat_obj,
  reduction.list = list("pca", "lsi"),
  dims.list = list(UMAP_DIMS_RNA, UMAP_DIMS_ATAC),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# Joint UMAP
cat("  Computing joint UMAP...\n")
seurat_obj <- RunUMAP(
  object = seurat_obj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  verbose = TRUE
)

# Clustering on joint graph
cat("  Computing clusters on joint graph...\n")
seurat_obj <- FindClusters(
  object = seurat_obj,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.8,
  verbose = FALSE
)

cat("  ✓ Joint UMAP complete\n")

################################################################################
# Create Visualization Plots
################################################################################

cat("\n================================================================================\n")
cat("STEP 4: Creating Visualization Plots\n")
cat("================================================================================\n")

# RNA UMAP plots
p1 <- DimPlot(seurat_obj, reduction = "umap.rna", group.by = "atac_sample") +
  ggtitle("RNA UMAP - by Sample") + NoLegend()

p2 <- DimPlot(seurat_obj, reduction = "umap.rna", group.by = "seurat_clusters",
              label = TRUE, repel = TRUE) +
  ggtitle("RNA UMAP - by Cluster") + NoLegend()

# ATAC UMAP plots
p3 <- DimPlot(seurat_obj, reduction = "umap.atac", group.by = "atac_sample") +
  ggtitle("ATAC UMAP - by Sample") + NoLegend()

p4 <- DimPlot(seurat_obj, reduction = "umap.atac", group.by = "seurat_clusters",
              label = TRUE, repel = TRUE) +
  ggtitle("ATAC UMAP - by Cluster") + NoLegend()

# Joint WNN UMAP plots
p5 <- DimPlot(seurat_obj, reduction = "wnn.umap", group.by = "atac_sample") +
  ggtitle("Joint WNN UMAP - by Sample")

p6 <- DimPlot(seurat_obj, reduction = "wnn.umap", group.by = "seurat_clusters",
              label = TRUE, repel = TRUE) +
  ggtitle("Joint WNN UMAP - by Cluster") + NoLegend()

# By cell type if available
if (!is.null(celltype_col)) {
  p7 <- DimPlot(seurat_obj, reduction = "umap.rna", group.by = celltype_col,
                label = TRUE, repel = TRUE) +
    ggtitle("RNA UMAP - by Cell Type") + NoLegend()

  p8 <- DimPlot(seurat_obj, reduction = "umap.atac", group.by = celltype_col,
                label = TRUE, repel = TRUE) +
    ggtitle("ATAC UMAP - by Cell Type") + NoLegend()

  p9 <- DimPlot(seurat_obj, reduction = "wnn.umap", group.by = celltype_col,
                label = TRUE, repel = TRUE) +
    ggtitle("Joint WNN UMAP - by Cell Type") + NoLegend()
}

# RNA modality weight
p10 <- FeaturePlot(seurat_obj, features = "RNA.weight", reduction = "wnn.umap") +
  ggtitle("RNA Modality Weight")

# Save plots
cat("  Saving plots...\n")

# QC plots
pdf(file.path(OUTPUT_DIR, "qc/02_rna_processing.pdf"), width = 16, height = 10)
print(pca_plot | elbow_plot)
if (!is.null(celltype_col)) {
  print(pca_plot_celltype)
}
dev.off()

pdf(file.path(OUTPUT_DIR, "qc/02_atac_processing.pdf"), width = 12, height = 6)
print(depth_cor_plot)
dev.off()

# UMAP overview
pdf(file.path(OUTPUT_DIR, "02_umap_overview.pdf"), width = 20, height = 12)
if (!is.null(celltype_col)) {
  print((p1 | p2 | p7) / (p3 | p4 | p8) / (p5 | p6 | p9))
} else {
  print((p1 | p2) / (p3 | p4) / (p5 | p6))
}
print(p10)
dev.off()

cat(sprintf("  ✓ Plots saved to: %s\n", OUTPUT_DIR))

################################################################################
# Save Results
################################################################################

cat("\n================================================================================\n")
cat("STEP 5: Saving Results\n")
cat("================================================================================\n")

# Save processed object
output_file <- file.path(OUTPUT_DIR, "integrated_seurat_processed.rds")
cat(sprintf("  Saving processed object to: %s\n", output_file))
saveRDS(seurat_obj, file = output_file)

# Create summary report
summary_text <- sprintf("
================================================================================
SIGNAC MULTIOME ANALYSIS - PROCESSING SUMMARY
================================================================================

Dataset: Nestin Ctrl vs Mut

Cells after QC: %d

RNA Processing:
  - Features: %d
  - Variable features: %d
  - PCs computed: %d
  - UMAP dimensions: %s

ATAC Processing:
  - Features: %d
  - Top features selected (min.cutoff = %d)
  - LSI components: %d
  - UMAP dimensions: %s (excluding component 1)

Joint Analysis:
  - Weighted nearest neighbor graph: ✓
  - Joint UMAP: ✓
  - Clusters identified: %d

Cell Type Annotation:
  - Cell type column: %s
  - Unique cell types: %d

Output Files:
  - integrated_seurat_processed.rds
  - 02_umap_overview.pdf
  - qc/02_rna_processing.pdf
  - qc/02_atac_processing.pdf

================================================================================
",
  ncol(seurat_obj),
  nrow(seurat_obj[["RNA"]]),
  length(VariableFeatures(seurat_obj)),
  N_PCS,
  paste(range(UMAP_DIMS_RNA), collapse = ":"),
  nrow(seurat_obj[["ATAC"]]),
  ATAC_MIN_CUTOFF,
  N_LSIS,
  paste(range(UMAP_DIMS_ATAC), collapse = ":"),
  length(unique(seurat_obj$seurat_clusters)),
  ifelse(is.null(celltype_col), "None", celltype_col),
  ifelse(is.null(celltype_col), 0, length(unique(seurat_obj@meta.data[[celltype_col]])))
)

cat(summary_text)

# Save summary
writeLines(summary_text, file.path(OUTPUT_DIR, "02_processing_summary.txt"))

cat("\n✓ Step 2 complete!\n")
cat(sprintf("  Output saved to: %s\n", OUTPUT_DIR))
cat("\nNext step: Run signac_03_celltype_analysis.R\n\n")
