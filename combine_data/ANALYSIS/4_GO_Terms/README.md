# 4_GO_Terms

This folder contains scripts for performing Gene Ontology (GO) enrichment analysis on differential gene expression results. The analysis uses the Enrichr API through gseapy to identify biological processes, cellular components, and molecular functions associated with differentially expressed genes.

## Scripts Overview

### Core GO Enrichment Scripts

#### raw_17_go_enrichment_final_L1.py
**Purpose**: Performs GO enrichment analysis on first layer cell type annotations
**Input**: 
- `results_from_raw/annotation_final.h5ad` - Finalized annotated dataset
- `results_from_raw/DEGs_cell_type_L1FC_0_25/biomarkers/sig_deg_lists/*_significant.csv` - Significant DEG lists

**Output**:
- `results_from_raw/DEGs_cell_type_L1FC_0_25/biomarkers/go_enrichment/[cell_type]/` - GO enrichment results

**Key Features**:
- Processes up- and down-regulated genes separately
- Uses 2023 GO gene sets (Biological Process, Cellular Component, Molecular Function)
- Background gene universe from complete dataset
- Custom dotplot visualizations with gene counts
- Excel output with filtered significant terms (Adj P < 0.05)

#### raw_17_go_enrichment_final_L1_combined.py
**Purpose**: Performs combined GO enrichment analysis (up + down genes together) for first layer annotations
**Input**: 
- Same inputs as individual L1 analysis

**Output**:
- `results_from_raw/DEGs_cell_type_L1FC_0_25/biomarkers/go_enrichment_combined/[cell_type]/` - Combined results

**Key Features**:
- Combines up- and down-regulated genes for each comparison
- Provides overall pathway analysis regardless of regulation direction
- Useful for identifying general biological processes affected
- Standardized output format matching individual analyses

#### raw_17_go_enrichment_final_L2.py
**Purpose**: Performs GO enrichment analysis on second layer cell type annotations
**Input**: 
- `results_from_raw/DEGs_cell_type_L2FC_0_25/biomarkers/sig_deg_lists/*_significant.csv` - Significant DEG lists

**Output**:
- `results_from_raw/DEGs_cell_type_L2FC_0_25/biomarkers/go_enrichment/[cell_type]/` - GO enrichment results

**Key Features**:
- More granular analysis using second layer cell type classifications
- Enables identification of specific pathways in subpopulations
- Same processing pipeline as L1 analysis

#### raw_17_go_enrichment_final_L2_combined.py
**Purpose**: Performs combined GO enrichment analysis for second layer annotations
**Input**: 
- Same inputs as individual L2 analysis

**Output**:
- `results_from_raw/DEGs_cell_type_L2FC_0_25/biomarkers/go_enrichment_combined/[cell_type]/` - Combined results

**Key Features**:
- Combined analysis for second layer cell types
- Provides broader pathway perspective for subpopulations

#### raw_17_go_enrichment_final_L2_new.py
**Purpose**: Performs GO enrichment analysis on updated second layer cell type annotations
**Input**: 
- `results_from_raw/DEGs_cell_type_L2_newFC_0_25/biomarkers/sig_deg_lists/*_significant.csv` - Updated DEG lists

**Output**:
- `results_from_raw/DEGs_cell_type_L2_newFC_0_25/biomarkers/go_enrichment/[cell_type]/` - GO enrichment results

**Key Features**:
- Uses latest refined cell type annotations
- Incorporates most recent manual curation
- Updated biological interpretations

#### raw_17_go_enrichment_MapMyCells_1.py
**Purpose**: Performs GO enrichment analysis using MapMyCells first layer annotations
**Input**: 
- `results_from_raw/DEGs_mapmycells_L2FC_0_25/biomarkers/sig_deg_lists/*_significant.csv` - MapMyCells DEG lists

**Output**:
- `results_from_raw/DEGs_mapmycells_L2FC_0_25/biomarkers/go_enrichment/[cell_type]/` - GO enrichment results

**Key Features**:
- Expert-curated cell type annotations
- Alternative annotation perspective for validation
- Cross-validation with automated annotations

#### raw_17_go_enrichment_MapMyCells_2.py
**Purpose**: Performs GO enrichment analysis using MapMyCells second layer annotations
**Input**: 
- `results_from_raw/DEGs_mapmycells_L2_2FC_0_25/biomarkers/sig_deg_lists/*_significant.csv` - MapMyCells L2 DEG lists

**Output**:
- `results_from_raw/DEGs_mapmycells_L2_2FC_0_25/biomarkers/go_enrichment/[cell_type]/` - GO enrichment results

**Key Features**:
- Detailed hierarchical analysis with MapMyCells annotations
- Most granular cell type resolution
- Comprehensive pathway mapping

## Analysis Pipeline

### 1. Input Processing
- Loads significant DEG lists from biomarker analysis
- Extracts cell type and comparison information from file paths
- Validates gene list format and content

### 2. Enrichment Analysis
- Uses gseapy.enrichr with Enrichr API
- Background gene universe from complete dataset
- Three GO categories: Biological Process, Cellular Component, Molecular Function
- Organism: Mouse (Mus musculus)

### 3. Result Filtering
- Filters by adjusted p-value < 0.05
- Calculates gene counts for each enriched term
- Removes redundant or non-significant terms

### 4. Visualization
- Custom dotplot generation with matplotlib
- Color coding by significance (-log10 adjusted p-value)
- Gene count annotations on plots
- Clean term names (GO IDs removed, underscores replaced)

### 5. Output Generation
- Excel files with complete filtered results
- High-resolution PNG plots (150 DPI)
- Organized directory structure by cell type

## Key Parameters

- **Organism**: Mouse (Mus musculus)
- **Gene Sets**: GO_Biological_Process_2023, GO_Cellular_Component_2023, GO_Molecular_Function_2023
- **P-value Threshold**: 0.05 (adjusted)
- **Top Terms for Plotting**: 15 most significant
- **Background**: All genes in the dataset
- **Enrichr Cutoff**: 0.1 (pre-filtering)

## Output Structure

```
results_from_raw/DEGs_[annotation_type]FC_0_25/biomarkers/
├── go_enrichment/
│   └── [cell_type]/
│       ├── [comparison]_up_go_enrichment.xlsx
│       ├── [comparison]_up_go_enrichment_dotplot.png
│       ├── [comparison]_down_go_enrichment.xlsx
│       └── [comparison]_down_go_enrichment_dotplot.png
└── go_enrichment_combined/
    └── [cell_type]/
        ├── [comparison]_combined_go_enrichment.xlsx
        └── [comparison]_combined_go_enrichment_dotplot.png
```

## Visualization Features

### Dotplot Characteristics
- **X-axis**: GeneRatio (normalized gene count)
- **Y-axis**: Clean GO term names
- **Color**: Significance (-log10 adjusted p-value) using red colormap
- **Size**: Constant dot size with gene count annotations
- **Layout**: Dynamic height based on number of terms
- **Grid**: Light grid on x-axis for better readability

### Plot Customization
- GO IDs removed for cleaner display
- Underscores replaced with spaces in term names
- Gene count annotations next to each point
- Informative titles with cell type and comparison information
- Professional publication-quality styling

## Dependencies

### Python Libraries
- scanpy
- pandas
- numpy
- matplotlib
- pathlib
- gseapy
- sys
- os

## Applications

### Biological Interpretation
- Identify affected biological processes in specific cell types
- Understand cellular component changes
- Discover molecular function alterations
- Compare pathway activation across conditions

### Validation and Quality Control
- Cross-validate findings across annotation systems
- Assess biological coherence of DEG results
- Identify potential artifacts or false positives

### Hypothesis Generation
- Generate mechanistic hypotheses for observed phenotypes
- Identify candidate pathways for functional validation
- Discover novel cell type-specific functions

## Integration with Downstream Analysis

GO enrichment results serve as input for:
- Pathway analysis tools (Reactome, KEGG)
- Network analysis and visualization
- Transcription factor binding analysis
- Drug target identification
- Integration with proteomics or other omics data

## Customization Options

- **Gene Sets**: Easy to add additional GO versions or other pathway databases
- **Thresholds**: Configurable p-value and gene count filters
- **Visualization**: Customizable colors, layouts, and plot types
- **Organisms**: Support for other species with appropriate Enrichr libraries
- **Background**: Alternative background gene sets if needed

## Error Handling and Robustness

- Comprehensive error handling for missing files or malformed data
- Validation of gene list format and content
- Graceful handling of API failures or network issues
- Detailed logging for debugging and troubleshooting
- Automatic skipping of problematic files with informative messages