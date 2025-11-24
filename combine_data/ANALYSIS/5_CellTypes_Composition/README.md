# 5_CellTypes_Composition

This folder contains scripts for analyzing and visualizing cell type composition across different experimental conditions and genotypes. The analysis calculates proportions and generates comprehensive visualizations and summaries for multiple annotation levels.

## Scripts Overview

### raw_19_cell_types_composition_with_annotations.py
**Purpose**: Analyzes cell type composition across experimental groups using multiple annotation systems
**Input**: 
- `results_from_raw/annotation_final.h5ad` - Finalized annotated dataset
- Fallback: `results_from_raw/final_annotation/merged_raw_final_annotated_simple.h5ad`

**Output**:
- `results_from_raw/Cell_types_composition/` - Complete composition analysis results

**Key Features**:
- Processes multiple annotation levels (L1, L2, L2_new)
- Calculates cell type proportions by experimental group
- Generates stacked bar plots with cell counts
- Creates comprehensive summary reports in multiple formats
- Automatic group creation (condition_genotype combinations)

## Analysis Pipeline

### 1. Data Loading and Preparation
- Loads finalized annotated dataset
- Creates combined group identifiers (condition_genotype)
- Validates annotation columns and categories
- Handles missing data gracefully

### 2. Proportion Calculation
For each annotation level:
- Groups cells by experimental condition and genotype
- Counts cells in each category within groups
- Calculates proportions relative to group totals
- Handles division by zero and missing categories

### 3. Visualization Generation
- Stacked bar plots showing cell type proportions
- Cell count annotations on plot segments
- Dynamic color mapping based on category count
- Professional publication-quality styling
- Automatic legend positioning and formatting

### 4. Summary Reports
Creates three types of output for each annotation level:
- **CSV Tables**: Raw proportions and counts
- **Text Summaries**: Human-readable composition reports
- **Markdown Tables**: Formatted tables for documentation

## Key Functions

### calculate_and_plot_proportions()
Core function that performs complete composition analysis for a given annotation key.

**Parameters**:
- `adata`: AnnData object with annotations
- `annotation_key`: Column name for cell type annotations
- `group_key`: Column name for experimental groups
- `output_dir`: Directory for saving results
- `plot_output_dir`: Directory for saving plots

**Processing Steps**:
1. Validates annotation key exists in data
2. Converts annotations to categorical format
3. Calculates counts and proportions
4. Saves results in multiple formats
5. Generates stacked bar visualization
6. Creates comprehensive summaries

### df_to_markdown()
Utility function that converts pandas DataFrames to Markdown table format with proper formatting for float values.

## Output Structure

```
results_from_raw/Cell_types_composition/
├── cell_type_L1_proportions_by_group.csv
├── cell_type_L1_composition_summary_by_group.txt
├── cell_type_L1_composition_table_by_group.md
├── cell_type_L1_proportions_stacked_bar.png
├── cell_type_L2_proportions_by_group.csv
├── cell_type_L2_composition_summary_by_group.txt
├── cell_type_L2_composition_table_by_group.md
├── cell_type_L2_proportions_stacked_bar.png
├── cell_type_L2_new_proportions_by_group.csv
├── cell_type_L2_new_composition_summary_by_group.txt
├── cell_type_L2_new_composition_table_by_group.md
└── cell_type_L2_new_proportions_stacked_bar.png
```

## Visualization Features

### Stacked Bar Plots
- **X-axis**: Experimental groups (Condition_Genotype combinations)
- **Y-axis**: Cell type proportions (0-1)
- **Colors**: Dynamic colormap based on number of categories
  - tab20 for ≤20 categories
  - gist_rainbow for >20 categories
- **Annotations**: Cell counts displayed on each segment
- **Legend**: Positioned below plot with automatic column arrangement

### Plot Customization
- Figure size: 12×7 inches
- High resolution: 300 DPI
- Automatic margin adjustment for legend
- White background for publication quality
- Rotated x-axis labels (45 degrees) for readability

## Annotation Levels Processed

### cell_type_L1
First layer cell type annotations providing broad categorization of major cell populations.

### cell_type_L2
Second layer cell type annotations offering more detailed subpopulation classification.

### cell_type_L2_new
Updated second layer annotations incorporating latest manual curation and refinements.

## Experimental Groups

The script automatically creates combined group identifiers:
- `Control_Emx1`
- `Mutant_Emx1`
- `Control_Nestin`
- `Mutant_Nestin`

## Output Formats

### CSV Files
- Raw proportion values with 4 decimal precision
- Cell counts for each category
- Machine-readable format for downstream analysis

### Text Summaries
- Human-readable composition reports
- Dominant cell type identification
- Detailed breakdown with percentages and counts
- Narrative format for easy interpretation

### Markdown Tables
- Formatted tables for documentation
- Both count and proportion tables
- Compatible with GitHub, documentation systems
- Easy integration into reports

## Dependencies

### Python Libraries
- scanpy
- pandas
- numpy
- matplotlib
- os
- sys
- random

### Custom Modules
- `functions` - Custom utility functions for the project

## Applications

### Biological Interpretation
- Identify cell type composition changes across conditions
- Detect shifts in cellular populations
- Understand genotype-specific effects on cell populations

### Quality Control
- Validate sample processing consistency
- Identify potential batch effects
- Assess experimental reproducibility

### Hypothesis Generation
- Generate hypotheses about cellular responses
- Identify target cell populations for further study
- Guide functional validation experiments

## Customization Options

### Annotation Keys
- Easy to add new annotation levels
- Supports any categorical annotation in adata.obs
- Automatic handling of missing categories

### Grouping Variables
- Configurable grouping keys
- Support for multiple experimental factors
- Flexible group creation logic

### Visualization Parameters
- Adjustable figure size and resolution
- Customizable color schemes
- Configurable text formatting and positioning

### Output Formats
- Extensible to additional file formats
- Customizable summary content
- Flexible naming conventions

## Integration with Downstream Analysis

Composition analysis results serve as input for:
- Statistical testing of proportion differences
- Integration with differential expression results
- Cross-validation with other annotation systems
- Meta-analysis across multiple experiments