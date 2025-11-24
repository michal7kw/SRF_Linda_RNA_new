# 0_Merge_CellTypists_Clean

This folder contains scripts for processing raw scRNA-seq data from Cell Ranger output, including data merging, quality control, cell type annotation using CellTypist, and data cleaning.

## Scripts Overview

### individual_sample_process.py
**Purpose**: Processes individual samples separately before merging
**Input**: 
- Raw Cell Ranger count data from `cellranger_final_count_data/[SAMPLE_NAME]/outs/filtered_feature_bc_matrix/`
- Samples: "Emx1_Ctrl", "Emx1_Mut", "Nestin_Ctrl", "Nestin_Mut"

**Output**:
- Individual processed `.h5ad` files in each sample directory
- QC plots in `plots/` subdirectories
- Sample processing summaries as `.csv` files
- Combined summary report in `individual_sample_results/all_samples_processing_summary.csv`

**Key Features**:
- Quality control filtering (min genes per cell: 500, max mito %: 20)
- Normalization and scaling
- PCA and UMAP dimensionality reduction
- Leiden clustering (resolution: 0.4)
- Comprehensive QC visualizations before/after filtering

### raw_0_merge_raw_and_process.py
**Purpose**: Merges multiple samples and performs joint processing
**Input**: 
- Raw Cell Ranger count data for all samples
- Same input directory structure as individual processing

**Output**:
- `results_from_raw/merged_raw_processed.h5ad` - Main processed dataset
- `results_from_raw/plots/` - UMAP visualizations and QC plots
- `results_from_raw/post_filtering_sample_summary.csv` - Sample summary statistics

**Key Features**:
- Concatenates datasets with 'outer' join to keep all genes
- Quality control with percentile-based filtering
- Highly variable genes identification (top 3000)
- Comprehensive QC visualizations with before/after comparisons
- Sample complexity analysis

### raw_1_annotate_merged_raw.py
**Purpose**: Annotates cell types using CellTypist models
**Input**:
- `results_from_raw/merged_raw_processed.h5ad` - Output from previous script
- CellTypist models: "Mouse_Isocortex_Hippocampus" and "Mouse_Dentate_Gyrus"

**Output**:
- `results_from_raw/annotated.h5ad` - Annotated dataset
- Cell type UMAP plots in `results_from_raw/plots/`
- Cell type summary CSV files
- Model comparison tables and heatmaps

**Key Features**:
- CellTypist annotation with majority voting and probability thresholding
- Dual model annotation (Isocortex/Hippocampus and Dentate Gyrus)
- Confidence scoring for annotations
- Cross-model comparison analysis
- Integration with existing Leiden clustering

### raw_5_clean_data.py
**Purpose**: Final data cleaning and simplification
**Input**:
- `results_from_raw/annotated.h5ad` - Annotated dataset from previous script

**Output**:
- `results_from_raw/annotated_cleaned.h5ad` - Final cleaned dataset

**Key Features**:
- Removes sex-related genes (Xist, Tsix, Ddx3x, etc.)
- Simplifies AnnData object by keeping only relevant columns/keys
- Preserves essential data while reducing file size
- Maintains raw data slot for downstream analysis

### alternative_umap_embedding/add_alternative_umap_embedding.py
**Purpose**: Adds alternative UMAP embeddings from external data sources
**Input**:
- Alternative data file with UMAP coordinates and cluster annotations
- Current processed dataset

**Output**:
- Updated dataset with alternative UMAP coordinates stored in `.obsm['X_umap_alt']`
- Alternative cluster annotations with `_alt` suffix

**Key Features**:
- Barcode matching between datasets with suffix stripping
- Handles missing data by filtering cells without matches
- Preserves original data while adding alternative embeddings
- Comprehensive logging of alignment process

## Processing Pipeline

1. **Individual Sample Processing** (`individual_sample_process.py`)
   - Process each sample separately for QC and initial analysis

2. **Joint Processing** (`raw_0_merge_raw_and_process.py`)
   - Merge all samples into single dataset
   - Apply joint QC and normalization

3. **Cell Type Annotation** (`raw_1_annotate_merged_raw.py`)
   - Annotate cells using CellTypist models
   - Generate confidence scores and comparisons

4. **Data Cleaning** (`raw_5_clean_data.py`)
   - Remove unwanted genes and simplify structure
   - Prepare final dataset for downstream analysis

5. **Alternative Embeddings** (`alternative_umap_embedding/add_alternative_umap_embedding.py`)
   - Optional: Add alternative UMAP coordinates from external sources

## Key Parameters

- **QC Thresholds**: Min genes per cell: 500, Max mito %: 20
- **Processing**: Top 3000 HVGs, 50 PCs, Leiden resolution 0.4
- **CellTypist**: Majority voting enabled, probability threshold 0.5
- **Samples**: Emx1_Ctrl, Emx1_Mut, Nestin_Ctrl, Nestin_Mut

## Dependencies

- scanpy
- pandas
- numpy
- matplotlib
- seaborn
- celltypist
- anndata
- pathlib
- tqdm

## Output Structure

```
results_from_raw/
├── merged_raw_processed.h5ad
├── annotated.h5ad
├── annotated_cleaned.h5ad
├── plots/
│   ├── umap_*.png
│   ├── qc_before_after_comparison.png
│   └── sample_complexity_histograms.png
├── post_filtering_sample_summary.csv
├── *cell_type_summary.csv
└── model_comparison_raw.csv

individual_sample_results/
└── all_samples_processing_summary.csv

cellranger_final_count_data/
├── Emx1_Ctrl/
│   ├── Emx1_Ctrl_processed.h5ad
│   ├── plots/
│   └── Emx1_Ctrl_processing_summary.csv
└── [similar structure for other samples]