# 2_Annotation

This folder contains scripts for finalizing cell type annotations and analyzing cluster similarities in the processed scRNA-seq dataset. These scripts integrate multiple annotation sources and provide tools for validating clustering results.

## Scripts Overview

### raw_9_finalize_annotation.py
**Purpose**: Finalizes cell type annotations by integrating multiple annotation sources
**Input**: 
- `results_from_raw/mapmycells.h5ad` - Dataset with MapMyCells annotations
- `models/Final_First_Layer.csv` - Final first layer cell type annotations
- `models/Final_Second_Layer.csv` - Final second layer cell type annotations  
- `models/Final_Second_Layer_new.csv` - Updated second layer cell type annotations

**Output**:
- `results_from_raw/annotation_final.h5ad` - Dataset with finalized annotations

**Key Features**:
- Loads final cell type annotations from multiple CSV files
- Handles barcode format conversion between CSV (underscore) and AnnData (hyphen) formats
- Adds three annotation layers to adata.obs:
  - `cell_type_L1` - First layer cell type annotations
  - `cell_type_L2` - Second layer cell type annotations  
  - `cell_type_L2_new` - Updated second layer annotations
- Fills missing annotations with 'Unknown'
- Generates UMAP visualizations for each annotation layer
- Includes utility function for highlighting specific cell types on plots

### raw_10_clusters_similiarty.py
**Purpose**: Analyzes cluster similarities using multiple approaches to validate clustering results
**Input**: 
- `results_from_raw/annotation_final.h5ad` - Finalized annotated dataset

**Output**:
- `results_from_raw/clusters_similarity/cluster_similarity_heatmap_leiden_0.4.png` - Similarity heatmap based on all genes
- `results_from_raw/clusters_similarity/cluster_similarity_heatmap_hvg_leiden_0.4.png` - Similarity heatmap based on highly variable genes
- `results_from_raw/clusters_similarity/dendrogram_dendrogram_leiden_0.4.png` - Cluster dendrogram

**Key Features**:
- **Cluster Similarity Analysis**: Calculates pairwise correlations between clusters based on:
  - Average expression of all genes
  - Average expression of highly variable genes (HVGs) only
  - Principal component representation (PCA space)
- **Visualization Tools**:
  - Clustered heatmaps with correlation coefficients
  - Hierarchical dendrograms showing cluster relationships
- **Data Handling**:
  - Supports both sparse and dense expression matrices
  - Handles categorical cluster labels properly
  - Configurable cluster key (default: 'leiden_0.4')

## Processing Pipeline

1. **Finalize Annotations** (`raw_9_finalize_annotation.py`)
   - Load dataset with existing annotations
   - Integrate final manual annotations from CSV files
   - Handle barcode format conversions
   - Generate visualization plots
   - Save finalized dataset

2. **Analyze Cluster Similarities** (`raw_10_clusters_similiarty.py`)
   - Load finalized annotated dataset
   - Calculate cluster similarities using multiple methods
   - Generate comprehensive similarity visualizations
   - Validate clustering robustness

## Key Parameters

- **Cluster Key**: 'leiden_0.4' (configurable)
- **Annotation Layers**: Three hierarchical levels of cell type annotations
- **Similarity Metrics**: Pearson correlation of mean expression
- **Visualization**: High-resolution plots (150 DPI) with customizable sizing

## Data Integration

The scripts integrate multiple annotation sources:
- **MapMyCells**: Automated cell type predictions
- **Manual Annotations**: Expert-curated cell type assignments
- **Hierarchical Structure**: Multiple layers of cell type granularity
- **Quality Control**: Missing data handling and format standardization

## Output Structure

```
results_from_raw/
├── annotation_final.h5ad                    # Finalized annotated dataset
└── clusters_similarity/
    ├── cluster_similarity_heatmap_leiden_0.4.png    # All-gene similarity
    ├── cluster_similarity_heatmap_hvg_leiden_0.4.png # HVG similarity
    └── dendrogram_dendrogram_leiden_0.4.png          # Cluster dendrogram
```

## Dependencies

### Python Libraries
- scanpy
- pandas
- numpy
- matplotlib
- seaborn
- scipy
- pathlib

## Usage Notes

1. **Sequence Dependency**: Scripts must be run in order as cluster similarity analysis requires the finalized annotations
2. **Barcode Format**: Handles automatic conversion between different barcode formats used by various tools
3. **Memory Considerations**: Similarity calculations can be memory-intensive for large datasets
4. **Visualization Settings**: Plots are automatically saved with high resolution for publication quality

## Analysis Applications

### Annotation Validation
- Compare automated vs manual annotations
- Identify discrepancies between annotation methods
- Validate hierarchical cell type classifications

### Cluster Quality Assessment
- Assess cluster separation and similarity
- Identify potentially merged or split clusters
- Validate clustering parameter choices

### Biological Interpretation
- Understand relationships between cell types
- Identify transitional cell states
- Guide downstream differential expression analysis

## Customization Options

- **Alternative Clustering**: Change `cluster_key` to analyze different clustering resolutions
- **Similarity Metrics**: Modify correlation method or use distance-based metrics
- **Gene Subsets**: Customize HVG selection or use specific gene sets
- **Visualization**: Adjust figure sizes, color schemes, and formatting