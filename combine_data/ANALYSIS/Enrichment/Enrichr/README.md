# Enrichr Pathway Analysis

## Overview

This folder contains scripts for performing Enrichr-based pathway enrichment analysis on differentially expressed genes from single-cell RNA-seq data. The pipeline identifies upregulated and downregulated genes between conditions and performs pathway analysis using multiple gene set libraries.

## Directory Structure

```
Enrichr/
├── enrichr_conditions.py              # Main Enrichr analysis pipeline
├── enrichr_conditions_filter.py       # Focused analysis with keyword filtering
├── enrichr_conditions.sh              # Shell script for execution
└── enrichr_conditions_filter.sh       # Shell script for filtered analysis
```

## Scripts

### enrichr_conditions.py

**Purpose**: Main pipeline for performing Enrichr pathway analysis on differentially expressed genes between conditions.

**Key Features**:
- Performs differential expression analysis for specific genotypes and clusters
- Separates upregulated and downregulated genes
- Runs Enrichr analysis across multiple gene set libraries
- Supports cluster merging options
- Parallel processing of gene sets
- Generates comprehensive summary reports

**Key Parameters**:
- `CLUSTER_OF_INTEREST = 'Mature GC'`: Target cell cluster
- `CLUSTER_COLUMN = 'cell_type_L2_new'`: Annotation column
- `GENE_SETS = ['KEGG_2019_Mouse', 'GO_Biological_Process_2023', 'Reactome_2022', 'MSigDB_Hallmark_2020']`: Gene set libraries
- `PADJ_THRESHOLD = 0.05`: Significance threshold
- `LOGFC_THRESHOLD = 0.25`: Log fold change threshold
- `MIN_CELLS = 50`: Minimum cells per gene

**Input Files**:
- `results_from_raw/annotation_final.h5ad`: Annotated single-cell dataset

**Output Files**:
- `results_from_raw/enrichr_between_conditions/{genotype}/{cluster}_Mutant_vs_Control/`: Analysis results
  - `Enrichr_Summary_All_Genesets.xlsx/.csv`: Combined results across all gene sets
  - Individual gene set results (stored in memory, combined in summary)

**Command Line Arguments**:
- `--genotype`: Filter by genotype ('Emx1', 'Nestin', 'both')
- `--merge_clusters`: Merge Immature GC and Mature GC clusters ('True', 'False')
- `--cluster_name`: Custom cluster name (overrides merge_clusters)
- `--cluster_column`: Custom cluster column (e.g., 'cell_type_L1')

**Usage**:
```bash
# Standard analysis for Emx1 genotype
python enrichr_conditions.py --genotype Emx1 --merge_clusters True

# Analysis for specific cluster
python enrichr_conditions.py --genotype Nestin --cluster_name GABA --cluster_column cell_type_L1

# Analysis for both genotypes
python enrichr_conditions.py --genotype both --merge_clusters False
```

**Analysis Workflow**:

1. **Data Loading and Filtering**:
   - Loads annotated AnnData object
   - Applies genotype filtering if specified
   - Handles cluster merging or custom cluster selection

2. **Differential Expression Analysis**:
   - Subsets data to target cluster
   - Filters genes by minimum cell count
   - Prepares data from raw counts (normalization + log1p)
   - Runs DEG analysis (Mutant vs Control)
   - Extracts up/down regulated genes based on thresholds

3. **Enrichr Analysis**:
   - Separates upregulated and downregulated gene lists
   - Runs Enrichr for each gene set library in parallel
   - Filters results by adjusted p-value threshold
   - Adds gene set library and direction information

4. **Result Compilation**:
   - Combines results across all gene sets
   - Sorts by direction and p-value
   - Saves comprehensive summary files

### enrichr_conditions_filter.py

**Purpose**: Performs focused analysis on Enrichr results by filtering for specific biological themes using keyword matching.

**Key Features**:
- Filters Enrichr results for specific biological themes
- Pre-configured focus areas (neurodegeneration, RNA splicing)
- Generates focused summary tables and visualizations
- Publication-ready plots with significance indicators

**Focus Areas**:
- **neuro**: Neurodegeneration & Cell Death
  - Keywords: 'neurodegeneration', 'apoptosis', 'cell death', 'autophagy', etc.
- **splicing**: RNA Splicing
  - Keywords: 'splicing', 'spliceosome', 'alternative splicing', 'intron', etc.

**Key Parameters**:
- `FOCUS_CONFIG`: Dictionary defining focus areas and keywords
- Search patterns: Case-insensitive keyword matching
- Visualization: Top 20 terms with significance indicators

**Input Files**:
- `results_from_raw/enrichr_between_conditions/{genotype}/{cluster}_Mutant_vs_Control/Enrichr_Summary_All_Genesets.xlsx`: Summary from main analysis

**Output Files**:
- `Enrichr_Focus_{file_suffix}.xlsx/.csv`: Filtered results for focus area
- `Enrichr_Focus_{file_suffix}.png`: Visualization of top hits

**Command Line Arguments**:
- `--genotype`: Genotype used in original analysis
- `--merge_clusters`: Cluster merging option
- `--cluster_name`: Custom cluster name
- `--focus`: Focus area ('neuro', 'splicing')

**Usage**:
```bash
# Filter for neurodegeneration-related pathways
python enrichr_conditions_filter.py --genotype Emx1 --merge_clusters True --focus neuro

# Filter for RNA splicing pathways
python enrichr_conditions_filter.py --genotype Nestin --merge_clusters False --focus splicing

# Analysis with custom cluster
python enrichr_conditions_filter.py --genotype both --cluster_name GABA --focus neuro
```

**Analysis Workflow**:

1. **Result Loading**:
   - Locates and loads Enrichr summary from main analysis
   - Validates file existence and format

2. **Keyword Filtering**:
   - Applies case-insensitive keyword matching
   - Searches all pathway terms for relevant keywords
   - Creates focused result subset

3. **Result Processing**:
   - Sorts filtered results by direction and p-value
   - Saves focused summary tables
   - Displays top hits in console

4. **Visualization**:
   - Creates bar plot of top 20 significant terms
   - Color-codes by direction (up/down regulated)
   - Includes significance threshold line
   - Publication-ready formatting

### Shell Scripts

#### enrichr_conditions.sh
Shell wrapper for the main Enrichr analysis pipeline. Provides convenient execution with predefined parameters.

#### enrichr_conditions_filter.sh
Shell wrapper for the filtered analysis. Simplifies execution of focused pathway analysis.

## Gene Set Libraries

The analysis uses the following Enrichr gene set libraries:

1. **KEGG_2019_Mouse**: KEGG pathways for mouse
2. **GO_Biological_Process_2023**: Gene Ontology Biological Process terms
3. **Reactome_2022**: Reactome pathway database
4. **MSigDB_Hallmark_2020**: Hallmark gene sets from MSigDB

## Output Interpretation

### Result Columns
- **Term**: Pathway/gene set name
- **Adjusted P-value**: Benjamini-Hochberg corrected p-value
- **Direction**: 'Upregulated_in_Mutant' or 'Upregulated_in_Control'
- **Gene_Set_Library**: Source gene set database
- **Overlap**: Number of genes from input list in pathway
- **Odds Ratio**: Enrichment strength
- **Combined Score**: Composite enrichment score

### Visualization Features
- **Color Coding**: Red for upregulated in mutant, blue for upregulated in control
- **Significance Line**: Green dashed line at p=0.05 threshold
- **Term Ordering**: Sorted by -log10(p-value) for easy interpretation
- **Publication Quality**: High-resolution PNG output

## Integration with Pipeline

### Dependencies
- **2_Annotation**: Final annotated dataset
- **3_DEGs**: Differential expression framework

### Outputs Used By
- **6_Pathways_Analysis**: Pathway-specific analyses
- **Publication figures**: Enrichment visualizations
- **Biological interpretation**: Pathway-level insights

## Customization Options

### Adding New Focus Areas
```python
# Add to FOCUS_CONFIG in enrichr_conditions_filter.py
'custom_focus': {
    'keywords': ['keyword1', 'keyword2', 'keyword3'],
    'file_suffix': 'Custom_Theme',
    'title_name': 'Custom Biological Theme'
}
```

### Modifying Gene Sets
```python
# Update GENE_SETS in enrichr_conditions.py
GENE_SETS = ['KEGG_2019_Mouse', 'GO_Biological_Process_2023', 
             'Reactome_2022', 'MSigDB_Hallmark_2020', 'Custom_Library']
```

### Adjusting Thresholds
```python
# Modify significance thresholds
PADJ_THRESHOLD = 0.01  # More stringent
LOGFC_THRESHOLD = 0.5   # Higher fold change requirement
```

## Troubleshooting

### Common Issues
1. **Missing Input File**: Ensure `enrichr_conditions.py` runs first
2. **No Significant Results**: Check thresholds and sample sizes
3. **Internet Connection**: Enrichr requires internet access
4. **Memory Issues**: Reduce gene sets or process sequentially

### Error Messages
- "No cells for genotype": Check genotype spelling and data availability
- "No significant enrichment": Consider relaxing thresholds
- "Enrichr summary file not found": Run main analysis first

## Performance Considerations

- **Parallel Processing**: Gene sets processed simultaneously
- **Memory Usage**: Moderate, depends on cluster size
- **Runtime**: 5-15 minutes per analysis
- **Internet Dependency**: Required for Enrichr API access

## Biological Applications

This Enrichr analysis is particularly useful for:
- **Pathway Discovery**: Identifying affected biological processes
- **Mechanism Elucidation**: Understanding disease-related pathways
- **Therapeutic Target Identification**: Finding druggable pathways
- **Comparative Analysis**: Comparing pathway activation across conditions
- **Validation**: Confirming expected biological changes

The focused analysis options enable targeted investigation of specific biological themes like neurodegeneration or RNA processing, making it ideal for hypothesis-driven research.