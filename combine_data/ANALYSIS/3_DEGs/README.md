# 3_DEGs

This folder contains scripts for performing differential gene expression (DGE) analysis across multiple cell type annotations and experimental conditions.

## Scripts Overview

### Core DGE Analysis Scripts

#### raw_16_dge_clusters.py
**Purpose**: Performs DGE analysis using Leiden clustering (leiden_0.4) as grouping
**Input**: 
- `results_from_raw/annotation_final.h5ad` - Finalized annotated dataset

**Output**:
- `results_from_raw/DEGs_leiden_04FC_0_25/` - Complete DGE results directory structure

**Key Features**:
- Uses 'leiden_0.4' clustering for group comparisons
- Creates normalized 'for_DEGs' layer from raw counts
- Four types of DGE analyses:
  1. Overall DGE (Mutant vs Control)
  2. Genotype-specific DGE (within each genotype)
  3. Genotype comparison DGE (Nestin vs Emx1 within conditions)
  4. Cell type marker identification (each cluster vs rest)

#### raw_16_dge_final_L1.py
**Purpose**: Performs DGE analysis using first layer cell type annotations
**Input**: 
- `results_from_raw/annotation_final.h5ad` - Finalized annotated dataset

**Output**:
- `results_from_raw/DEGs_cell_type_L1FC_0_25/` - Complete DGE results directory structure

**Key Features**:
- Uses 'cell_type_L1' annotations for group comparisons
- Same four types of DGE analyses as cluster-based approach
- Provides biologically meaningful cell type comparisons

#### raw_16_dge_final_L2.py
**Purpose**: Performs DGE analysis using second layer cell type annotations
**Input**: 
- `results_from_raw/annotation_final.h5ad` - Finalized annotated dataset

**Output**:
- `results_from_raw/DEGs_cell_type_L2FC_0_25/` - Complete DGE results directory structure

**Key Features**:
- Uses 'cell_type_L2' annotations for more granular cell type comparisons
- Enables analysis of subpopulations within broader cell types

#### raw_16_dge_final_L2_new.py
**Purpose**: Performs DGE analysis using updated second layer cell type annotations
**Input**: 
- `results_from_raw/annotation_final.h5ad` - Finalized annotated dataset

**Output**:
- `results_from_raw/DEGs_cell_type_L2_newFC_0_25/` - Complete DGE results directory structure

**Key Features**:
- Uses 'cell_type_L2_new' annotations for refined cell type classifications
- Incorporates latest manual curation and annotation updates

#### raw_16_dge_MapMyCells_1.py
**Purpose**: Performs DGE analysis using MapMyCells first layer annotations
**Input**: 
- `results_from_raw/annotation_final.h5ad` - Finalized annotated dataset

**Output**:
- `results_from_raw/DEGs_mapmycells_L2FC_0_25/` - Complete DGE results directory structure

**Key Features**:
- Uses 'mapmycells_first_layer' annotations for expert-curated cell type comparisons
- Provides alternative annotation perspective for validation

#### raw_16_dge_MapMyCells_2.py
**Purpose**: Performs DGE analysis using MapMyCells second layer annotations
**Input**: 
- `results_from_raw/annotation_final.h5ad` - Finalized annotated dataset

**Output**:
- `results_from_raw/DEGs_mapmycells_L2_2FC_0_25/` - Complete DGE results directory structure

**Key Features**:
- Uses 'mapmycells_second_layer' annotations for detailed hierarchical analysis
- Enables comparison across multiple annotation systems

#### raw_17_dge_merged_clusters.py
**Purpose**: Performs DGE analysis on merged clusters or specialized groupings
**Input**: 
- `results_from_raw/annotation_final.h5ad` - Finalized annotated dataset

**Output**:
- Additional DGE results for custom cluster groupings

**Key Features**:
- Handles custom cluster definitions or merged populations
- Supports specialized biological hypotheses

### Visualization Scripts

#### create_alternative_volcano_plots.py
**Purpose**: Creates custom-styled volcano plots from existing DGE results
**Input**: 
- DGE results CSV files from any of the core analysis scripts

**Output**:
- Alternative volcano plots in `plots_alt/` subdirectories
- Custom color scheme: upregulated (orange #FDAC5D), downregulated (purple #7C3294)
- Gene count annotations in plot regions

**Key Features**:
- Custom styling with publication-quality appearance
- Data point counts in six plot regions
- No gene labels for cleaner visualization
- Processes multiple comparison types:
  - Both genomes condition comparison
  - Genotype-specific condition comparison
  - Condition-specific genotype comparison

## Analysis Types

### 1. Overall DGE (Mutant vs Control)
- Compares all mutant samples against all control samples
- Grouped by cell type/cluster annotations
- Identifies overall treatment effects within each cell type

### 2. Genotype-Specific DGE
- Compares mutant vs control within each genotype (Emx1, Nestin)
- Reveals genotype-specific treatment responses
- Accounts for genotype-background interactions

### 3. Genotype Comparison DGE
- Compares Nestin vs Emx1 within each condition (Control, Mutant)
- Identifies baseline differences between genotypes
- Reveals condition-dependent genotype effects

### 4. Cell Type Marker Identification
- Identifies marker genes for each cell type vs all others
- Uses Wilcoxon rank-sum test (configurable to t-test)
- Provides cell type-specific gene signatures

## Key Parameters

- **Fold Change Threshold**: 0.25 (log2 scale)
- **P-value Threshold**: 0.05 (adjusted)
- **Normalization**: Total counts to 10,000, log1p transformation
- **Statistical Tests**: Wilcoxon rank-sum (default), t-test (optional)
- **Multiple Testing Correction**: Benjamini-Hochberg (FDR)

## Output Structure

Each DGE analysis creates a standardized directory structure:

```
results_from_raw/DEGs_[annotation_type]FC_0_25/
├── dge_res/
│   ├── both_geno_cond_comp/
│   │   ├── dge_both_genomes_[celltype]_mut_vs_ctrl.csv
│   │   └── dge_both_genomes_[celltype]_mut_vs_ctrl.xlsx
│   ├── geno_spec_cond_comp/
│   │   ├── dge_[genotype]_[celltype]_mut_vs_ctrl.csv
│   │   └── dge_[genotype]_[celltype]_mut_vs_ctrl.xlsx
│   ├── cond_spec_geno_comp/
│   │   ├── dge_[condition]_cond_[celltype]_Nestin_vs_Emx1.csv
│   │   └── dge_[condition]_cond_[celltype]_Nestin_vs_Emx1.xlsx
│   └── [comparison]_comparison_skipped_[annotation].csv
├── plots/
│   ├── both_geno_cond_comp/
│   ├── geno_spec_cond_comp/
│   └── cond_spec_geno_comp/
├── plots_alt/  # Alternative volcano plots
└── biomarkers/
    └── [celltype]_vs_rest.csv
```

## Dependencies

### Python Libraries
- scanpy
- pandas
- numpy
- matplotlib
- seaborn
- anndata
- functions_degs (custom module)

## Data Processing Pipeline

1. **Data Loading**: Load finalized annotated dataset
2. **Layer Creation**: Create 'for_DEGs' layer with normalized counts
3. **Quality Checks**: Validate metadata and data structure
4. **DGE Analysis**: Run four types of comparisons
5. **Visualization**: Generate standard and alternative plots
6. **Export**: Save results in CSV and Excel formats

## Applications

### Biological Discovery
- Identify treatment-responsive genes in specific cell types
- Discover cell type-specific markers
- Understand genotype-dependent effects

### Validation and Quality Control
- Cross-validate results across annotation systems
- Assess consistency of biological findings
- Identify robust transcriptional changes

### Hypothesis Generation
- Generate candidate genes for functional validation
- Identify pathways for further investigation
- Discover novel cell type-specific responses

## Customization Options

- **Annotation Systems**: Easy to add new cell type annotations
- **Comparison Types**: Framework supports additional contrast types
- **Statistical Parameters**: Configurable thresholds and tests
- **Visualization**: Customizable styling and color schemes
- **Output Formats**: Support for additional file formats

## Integration with Downstream Analysis

The DGE results serve as input for:
- Gene set enrichment analysis (GSEA)
- Pathway analysis
- Network analysis
- Transcription factor motif analysis
- Integration with other omics data types