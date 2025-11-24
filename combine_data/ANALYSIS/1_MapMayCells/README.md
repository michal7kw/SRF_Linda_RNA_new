# 1_MapMayCells

This folder contains scripts for adding MapMyCells annotations to the processed scRNA-seq data and converting the data to different formats (Seurat and Loupe) for downstream analysis and visualization.

## Scripts Overview

### raw_6_add_mapmycells_annotations.py
**Purpose**: Adds MapMyCells annotations to the processed AnnData object
**Input**: 
- `results_from_raw/annotated_cleaned.h5ad` - Cleaned annotated dataset from previous step
- `models/first_layer_LB_cleaned.csv` - First layer MapMyCells annotations
- `models/second_layer_LB_cleaned.csv` - Second layer MapMyCells annotations

**Output**:
- `results_from_raw/mapmycells.h5ad` - Dataset with MapMyCells annotations added

**Key Features**:
- Loads two-layer MapMyCells annotations from CSV files
- Modifies barcode format (replaces last underscore with hyphen) for compatibility
- Maps annotations to cells based on barcode matching
- Adds 'mapmycells_first_layer' and 'mapmycells_second_layer' columns to adata.obs
- Handles missing annotations gracefully (leaves as NaN)
- Provides detailed mapping statistics

### raw_7_convert_merged_to_seurat.R
**Purpose**: Converts AnnData object to Seurat format for R-based analysis
**Input**: 
- `results_from_raw/mapmycells.h5ad` - AnnData object with MapMyCells annotations

**Output**:
- `results_from_raw/mapmycells.rds` - Seurat object

**Key Features**:
- Uses schard package for AnnData to Seurat conversion
- Preserves raw counts from .raw slot in Seurat 'counts' slot
- Handles both WSL and Windows path environments
- Comprehensive verification of loaded data structure
- Checks for presence of expected metadata columns:
  - Sample information (sample, condition, genotype)
  - QC metrics (n_genes_by_counts, total_counts, pct_counts_mt)
  - CellTypist annotations (DG_majority_voting, ISO_majority_voting)
  - MapMyCells annotations (mapmycells_first_layer, mapmycells_second_layer)
  - Clustering information (leiden_0.4)
- Verifies presence of transgene feature (Rosa26-SBP1)

### raw_8_create_cloupe_from_merged.R
**Purpose**: Converts Seurat object to Loupe (.cloupe) format for visualization in 10x Genomics Loupe Browser
**Input**: 
- `results_from_raw/mapmycells.rds` - Seurat object from previous script

**Output**:
- `results_from_raw/mapmycells_loupe.cloupe` - Loupe file for 10x Genomics viewer

**Key Features**:
- Uses loupeR package for Seurat to Loupe conversion
- Corrects barcode format for LoupeR compatibility:
  - Changes final hyphen to underscore (e.g., AAACAGCCAAGCTTTG-1-0 → AAACAGCCAAGCTTTG-1_0)
- Handles Cell Ranger barcode format requirements
- Comprehensive debugging of feature names
- Verifies presence of transgene feature in final object
- Supports both WSL and Windows environments

## Processing Pipeline

1. **Add MapMyCells Annotations** (`raw_6_add_mapmycells_annotations.py`)
   - Load cleaned AnnData object
   - Add two-layer annotations from CSV files
   - Save annotated dataset

2. **Convert to Seurat** (`raw_7_convert_merged_to_seurat.R`)
   - Convert AnnData to Seurat format using schard
   - Preserve raw counts and metadata
   - Verify data integrity

3. **Create Loupe File** (`raw_8_create_cloupe_from_merged.R`)
   - Convert Seurat to Loupe format
   - Fix barcode formatting for compatibility
   - Generate visualization-ready file

## Key Parameters

- **Barcode Format**: Converts between formats for different tools
- **Annotation Layers**: First and second layer MapMyCells annotations
- **Environment Support**: Handles both WSL and Windows paths
- **Data Preservation**: Maintains raw counts and all metadata

## Dependencies

### Python Dependencies
- scanpy
- pandas
- os
- warnings

### R Dependencies
- schard (for AnnData to Seurat conversion)
- Seurat
- loupeR (for Loupe file creation)
- Matrix
- utils

## Input File Structure

```
models/
├── first_layer_LB_cleaned.csv
│   ├── Barcode (cell barcodes)
│   └── first layer LB (annotations)
└── second_layer_LB_cleaned.csv
    ├── Barcode (cell barcodes)
    └── second layer LB (annotations)
```

## Output Structure

```
results_from_raw/
├── mapmycells.h5ad          # AnnData with MapMyCells annotations
├── mapmycells.rds           # Seurat object
└── mapmycells_loupe.cloupe  # Loupe file for visualization
```

## Data Integration

The scripts integrate multiple annotation sources:
- **CellTypist**: Automated cell type annotations
- **MapMyCells**: Manual/expert annotations with two hierarchical layers
- **Metadata**: Sample, condition, genotype information
- **QC Metrics**: Quality control measurements

This comprehensive annotation approach enables:
- Cross-validation of annotation methods
- Hierarchical cell type classification
- Flexible analysis in both Python (AnnData) and R (Seurat) environments
- Interactive visualization in 10x Genomics Loupe Browser

## Usage Notes

1. Scripts must be run in sequence as each depends on the output of the previous
2. Barcode format modifications are essential for compatibility between tools
3. The conversion processes preserve all original data while adding new annotations
4. Environment detection (WSL/Windows) ensures path compatibility across systems