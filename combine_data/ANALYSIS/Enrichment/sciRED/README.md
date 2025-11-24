# sciRED (Single-cell RNA-seq Ensemble Decomposition) Directory

This directory contains scripts and tools for performing sciRED analysis, a machine learning approach for discovering interpretable biological factors from single-cell RNA-seq data. sciRED combines factor analysis with ensemble classification to identify factors that explain biological covariates.

## Overview

The sciRED analysis pipeline implements a four-step approach:

1. **Factor Discovery**: Unsupervised factor analysis using GLM regression and PCA
2. **Factor-Covariate Association**: Ensemble classification to link factors with biological covariates
3. **Interpretability Scoring**: Quantification of factor interpretability using FIST metrics
4. **Biological Interpretation**: Pathway enrichment analysis of factor loadings

## Main Analysis Scripts

### Core Analysis Scripts

#### [`conditions_genotypes_factors.py`](conditions_genotypes_factors.py)
**Purpose**: Comprehensive sciRED analysis for both conditions and genotypes with full biological interpretation

**Input Files**:
- `combine_data/results_from_raw/annotation_final.h5ad` - Main scRNA-seq data with annotations

**Key Parameters**:
- `NUM_COMPONENTS`: Number of factors to extract (default: 30)
- `NUM_GENES`: Number of highly variable genes to use (default: 2000)
- `INTERESTING_FACTOR_IDS`: Factors for detailed analysis (default: [0, 2, 10, 12])
- Target cell type: `cell_type_L2_new == "Mature GC"`

**Output Files**:
- `enrichment_results_summary/` - Comprehensive enrichment results
  - `summary_all_factors.csv` - Summary statistics for all factors
  - `{factor}_enrichr_{direction}_{gene_set}.csv` - Enrichr results
  - `{factor}_gsea_{gene_set}.csv` - GSEA results
- `loadings_with_proteins/` - Factor loadings with protein annotations
  - `factor_{n}_loadings_with_proteins.csv` - Gene loadings with protein names
- `gsea_results/` - Detailed GSEA analysis results
- `enrichr/` - Detailed Enrichr analysis results
- `genes_for_string.txt` - Genes formatted for STRING database analysis

**Key Features**:
- Complete sciRED pipeline with all four steps
- Directional FCAT analysis preserving effect direction
- Both GSEA prerank and Enrichr enrichment analysis
- Protein name annotation using mygene API
- STRING database integration for protein interaction analysis
- Comprehensive visualization of factor-covariate relationships

#### [`conditions_factors.py`](conditions_factors.py)
**Purpose**: Focused sciRED analysis for condition effects only

**Input Files**:
- `combine_data/results_from_raw/annotation_final.h5ad` - Main scRNA-seq data

**Key Parameters**:
- `NUM_COMPONENTS`: Number of factors to extract (default: 30)
- `NUM_GENES`: Number of highly variable genes (default: 2000)
- Target: Emx1 genotype, Mature GC cells
- Focus: Condition effects (Control vs. Mutant)

**Output Files**:
- Factor-covariate association tables (FCAT)
- Interpretability scores (FIST table)
- Factor loading plots and visualizations
- UMAP projections colored by factors and conditions

**Key Features**:
- Streamlined analysis focusing on condition effects
- Emx1 genotype-specific analysis
- Factor interpretability quantification
- Directionality analysis for binary conditions

#### [`conditions_celltypes_factors.py`](conditions_celltypes_factors.py)
**Purpose**: Extended analysis including cell type covariates alongside conditions

**Key Features**:
- Multi-covariate analysis (conditions + cell types)
- Enhanced factor discovery with multiple biological variables
- Extended FCAT analysis across multiple covariate types

## Analysis Pipeline Steps

### Step 1: Factor Discovery
- **GLM Regression**: Removes technical covariates (library size)
- **PCA Analysis**: Extracts principal components from residuals
- **Varimax Rotation**: Improves factor interpretability
- **Key Functions**: [`glm.poissonGLM()`](conditions_genotypes_factors.py:98), [`rot.varimax()`](conditions_genotypes_factors.py:172)

### Step 2: Factor-Covariate Association
- **Ensemble Classification**: Uses multiple ML classifiers (Logistic Regression, Decision Tree, Random Forest, XGBoost)
- **FCAT Calculation**: Factor-Covariate Association Table
- **Directional Analysis**: Preserves effect direction for binary covariates
- **Key Functions**: [`efca.FCAT()`](conditions_genotypes_factors.py:229), [`FCAT_directional()`](conditions_genotypes_factors.py:735)

### Step 3: Interpretability Scoring
- **FIST Metrics**: Factor Interpretability Score Table
  - **Bimodality**: Separability of factor scores
  - **Specificity**: Simpson diversity index
  - **Effect Size**: Factor variance
  - **Homogeneity**: Average scaled variance
- **Key Functions**: [`met.FIST()`](conditions_genotypes_factors.py:500)

### Step 4: Biological Interpretation
- **Gene Loading Analysis**: Extracts top/bottom genes per factor
- **Enrichment Analysis**: Both GSEA prerank and Enrichr
- **Protein Annotation**: Maps genes to protein names
- **STRING Integration**: Formats genes for protein interaction analysis

## Directory Structure

```
sciRED/
├── conditions_genotypes_factors.py      # Comprehensive analysis
├── conditions_factors.py               # Condition-focused analysis
├── conditions_celltypes_factors.py     # Multi-covariate analysis
├── enrichment_results_summary/         # Summary results
│   ├── summary_all_factors.csv
│   ├── Factor_*_enrichr_*.csv
│   └── Factor_*_gsea_*.csv
├── loadings_with_proteins/             # Protein-annotated loadings
│   └── factor_*_loadings_with_proteins.csv
├── gsea_results/                       # Detailed GSEA results
│   └── Factor_*/{gene_set}/
├── enrichr/                           # Detailed Enrichr results
│   └── Factor_*/{direction}_{gene_set}/
├── genes_for_string.txt               # STRING database input
└── Factor_*_Enrichr_summary.pdf       # Visual summaries
```

## Key Configuration Parameters

### Factor Analysis
- `NUM_COMPONENTS`: Number of factors to extract (30)
- `NUM_GENES`: Number of HVGs to use (2000)
- `NUM_COMP_TO_VIS`: Number of components to visualize (3)

### Enrichment Analysis
- **Gene Sets**: GO Biological Process, GO Cellular Component, GO Molecular Function, KEGG, Reactome, WikiPathways
- **Top Genes**: 150 genes for Enrichr analysis
- **GSEA Parameters**: min_size=10, max_size=500, permutation_num=100

### Classification
- **Models**: Logistic Regression, Decision Tree, Random Forest, XGBoost, KNeighbors
- **Scaling**: Standard scaling of importance scores
- **Mean**: Arithmetic mean of classifier scores

## Output Interpretation

### Factor Loadings
- **Positive Loadings**: Genes upregulated in one condition
- **Negative Loadings**: Genes upregulated in the other condition
- **Magnitude**: Strength of association with the factor

### FCAT Scores
- **High Values**: Strong factor-covariate association
- **Direction**: Indicates which condition the factor represents
- **Threshold**: Determined by Otsu's method for significance

### FIST Scores
- **Bimodality**: Separation of factor scores (higher = better)
- **Specificity**: How specific the factor is to covariates
- **Effect Size**: Variance explained by the factor
- **Homogeneity**: Consistency within covariate groups

## Usage Examples

### Run comprehensive analysis:
```python
# Execute the main analysis script
python conditions_genotypes_factors.py
```

### Extract genes for STRING analysis:
```python
# Get top 20 genes from each factor
top_genes = extract_genes_for_string(
    combined_loading_tables, 
    method='top_n', 
    top_n=20, 
    unique_only=True
)
```

### Analyze specific factor:
```python
# Focus on Factor 1
INTERESTING_FACTOR_ID = 0
factor_loadings = factor_loading[:, INTERESTING_FACTOR_ID]
gene_loadings = pd.Series(factor_loadings, index=genes)
sorted_loadings = gene_loadings.sort_values(ascending=False)
```

## Dependencies

### Python Packages
- **Core**: numpy, pandas, matplotlib, seaborn
- **Machine Learning**: sklearn, xgboost
- **Statistical**: statsmodels
- **sciRED**: sciRED library (ensembleFCA, glm, rotations, metrics)
- **Enrichment**: gseapy, mygene
- **Data**: anndata

### External Tools
- **STRING Database**: For protein interaction analysis
- **Enrichr API**: For pathway enrichment
- **MSigDB**: For gene set databases

## Key Features

### Advanced Analysis
- **Directional FCAT**: Preserves effect direction for binary covariates
- **Multi-covariate**: Handles multiple biological variables simultaneously
- **Protein Integration**: Maps genes to proteins for network analysis
- **Parallel Processing**: Accelerates enrichment analysis

### Visualization
- **Factor Plots**: PCA, UMAP, loading scatter plots
- **FCAT Heatmaps**: Factor-covariate association visualization
- **FIST Tables**: Interpretability metric visualization
- **Enrichment Summaries**: Comprehensive pathway analysis plots

### Data Integration
- **STRING Database**: Protein interaction network analysis
- **Multiple Enrichment**: Both GSEA and Enrichr methods
- **Cross-platform**: Compatible with various annotation formats

## Notes

- All analyses focus on Mature GC cells from the dentate gyrus
- Factor number (30) balances interpretability and comprehensiveness
- Directional analysis requires careful interpretation of effect signs
- Protein annotations depend on mygene API availability
- Results are saved in multiple formats for downstream analysis