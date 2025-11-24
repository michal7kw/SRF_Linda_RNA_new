# Splicing Genes CRE Analysis Pipeline

This pipeline analyzes CREs (cis-regulatory elements) associated with splicing-related genes from Reactome pathways, using literature-based CRE-gene correlations from ENCODE Table 16.

## Overview

The workflow identifies CREs linked to splicing machinery genes and visualizes ATAC-seq signal patterns at these regulatory elements. The analysis focuses on understanding how SRF mutations affect chromatin accessibility at splicing gene regulatory regions.

## Pipeline Steps

### Step 1: Extract Splicing Gene CREs (`1_extract_splicing_gene_CREs.py`)
- **Input**: Splicing genes list (1,138 genes) + ENCODE Table 16 (literature CRE-gene correlations)
- **Output**: 
  - `splicing_genes_CREs_all_celltypes.tsv` - All CREs linked to splicing genes
  - `splicing_genes_CREs_GABA.tsv` - CREs in GABA/hippocampal cell types
  - `splicing_genes_CREs_GABA_specific.tsv` - GABA-exclusive CREs
- **Filters**: FDR < 0.05, |PCC| > 0.1
- **Key Feature**: Uses published literature correlations (not Signac peak-gene links)

### Step 2: Convert to BED Format (`2_convert_splicing_CREs_to_bed.py`)
- **Input**: TSV files from Step 1
- **Output**: 
  - `splicing_genes_CREs_all.bed` - All CREs in BED6 format
  - `splicing_genes_CREs_GABA.bed` - GABA CREs in BED6 format
- **Purpose**: Prepare CRE coordinates for deepTools analysis

### Step 3: Create ATAC Signal Profiles (`3_create_splicing_profiles.sh`)
- **Input**: BED files + GABA BigWig files
- **Output**: 
  - Metaprofiles: `metaprofile_nestin_ctrl_vs_mut.png`, `metaprofile_emx1_ctrl_vs_mut.png`
  - Individual CRE plots: `individual_Nestin_<gene>_<cre_id>.png`, `individual_Emx1_<gene>_<cre_id>.png`
  - Signal matrices: `matrix_GABA_nestin.gz`, `matrix_GABA_emx1.gz`
- **Performance**: Default skips individual plots (3-5 min), full mode creates all plots
- **Usage**: `sbatch 3_create_splicing_profiles.sh` (fast) or `SKIP_INDIVIDUAL=0 sbatch 3_create_splicing_profiles.sh` (full)

### Step 4: Create Traditional Heatmaps (`4_create_heatmaps_splicing_genes.sh`)
- **Input**: Same as Step 3
- **Output**: 
  - Heatmaps: `heatmap_all_celltypes.png`, `heatmap_GABA.png`
  - Metaprofiles: `metaprofile_all_celltypes.png`, `metaprofile_GABA.png`
  - Genotype-specific: `heatmap_GABA_nestin.png`, `heatmap_GABA_emx1.png`
- **Method**: Uses deepTools plotHeatmap/plotProfile commands

### Step 5: Visualize BigWig Signal (`5_visualize_bigwig_signal.sh`)
- **Input**: Direct BigWig file reading
- **Output**: 
  - Metaprofiles: `metaprofile_nestin_ctrl_vs_mut.png`, `metaprofile_emx1_ctrl_vs_mut.png`
  - Individual plots: `individual_Nestin_<gene>_<cre_id>.png`, `individual_Emx1_<gene>_<cre_id>.png`
  - Statistics: `summary_nestin.tsv`, `summary_emx1.tsv`
- **Method**: Direct BigWig reading with pyBigWig (alternative to deepTools)

### Step 6: Create Custom Comparisons (`6_create_custom_comparisons.sh`)
- **Input**: ALL CREs with custom genotype comparisons
- **Output**: 
  - Custom metaprofiles: `metaprofile_nestin_ctrl_vs_nestin_mut.png`, `metaprofile_nestin_ctrl_vs_emx1_mut.png`, `metaprofile_nestin_mut_vs_emx1_mut.png`
  - Individual plots: `individual_comparison1_<gene>_<cre_id>.png`, `individual_comparison2_<gene>_<cre_id>.png`, `individual_comparison3_<gene>_<cre_id>.png`
  - Statistics: `comparison_statistics.txt`
- **Comparisons**: 
  1. Nestin-Ctrl vs Nestin-Mut (within-genotype mutation effect)
  2. Nestin-Ctrl vs Emx1-Mut (cross-genotype comparison)
  3. Nestin-Mut vs Emx1-Mut (mutant-to-mutant comparison)

## Generated Plots

### Metaprofiles (Ctrl vs Mut Comparisons)
![Nestin Ctrl vs Mut](output/splicing_genes_analysis/heatmaps_deeptools/profiles/metaprofile_nestin_ctrl_vs_mut.png)
*ATAC signal at splicing gene CREs comparing Nestin control vs mutant samples*

![Emx1 Ctrl vs Mut](output/splicing_genes_analysis/heatmaps_deeptools/profiles/metaprofile_emx1_ctrl_vs_mut.png)
*ATAC signal at splicing gene CREs comparing Emx1 control vs mutant samples*

### Heatmaps (Traditional deepTools)
![All Cell Types Heatmap](output/splicing_genes_analysis/heatmaps_deeptools/heatmap_all_celltypes.png)
*ATAC signal heatmap for all CREs across all cell types*

![GABA Cell Types Heatmap](output/splicing_genes_analysis/heatmaps_deeptools/heatmap_GABA.png)
*ATAC signal heatmap for CREs in GABA/hippocampal cell types*

![GABA Nestin Heatmap](output/splicing_genes_analysis/heatmaps_deeptools/heatmap_GABA_nestin.png)
*ATAC signal heatmap for GABA CREs in Nestin samples*

![GABA Emx1 Heatmap](output/splicing_genes_analysis/heatmaps_deeptools/heatmap_GABA_emx1.png)
*ATAC signal heatmap for GABA CREs in Emx1 samples*

### Metaprofiles (Traditional deepTools)
![All Cell Types Metaprofile](output/splicing_genes_analysis/heatmaps_deeptools/metaprofile_all_celltypes.png)
*Average ATAC signal profile for all CREs across all cell types*

![GABA Cell Types Metaprofile](output/splicing_genes_analysis/heatmaps_deeptools/metaprofile_GABA.png)
*Average ATAC signal profile for CREs in GABA/hippocampal cell types*

![GABA Nestin Metaprofile](output/splicing_genes_analysis/heatmaps_deeptools/metaprofile_GABA_nestin.png)
*Average ATAC signal profile for GABA CREs in Nestin samples*

![GABA Emx1 Metaprofile](output/splicing_genes_analysis/heatmaps_deeptools/metaprofile_GABA_emx1.png)
*Average ATAC signal profile for GABA CREs in Emx1 samples*

### BigWig Direct Visualizations
![Nestin BigWig Metaprofile](output/splicing_genes_analysis/bigwig_profiles/metaprofile_nestin_ctrl_vs_mut.png)
*Direct BigWig signal visualization for Nestin samples*

![Emx1 BigWig Metaprofile](output/splicing_genes_analysis/bigwig_profiles/metaprofile_emx1_ctrl_vs_mut.png)
*Direct BigWig signal visualization for Emx1 samples*

### Custom Comparison Plots
![Nestin Ctrl vs Nestin Mut](output/splicing_genes_analysis/custom_comparisons/profiles/metaprofile_nestin_ctrl_vs_nestin_mut.png)
*Within-genotype comparison: effect of mutation in Nestin background*

![Nestin Ctrl vs Emx1 Mut](output/splicing_genes_analysis/custom_comparisons/profiles/metaprofile_nestin_ctrl_vs_emx1_mut.png)
*Cross-genotype comparison: wild-type Nestin vs mutant Emx1*

![Nestin Mut vs Emx1 Mut](output/splicing_genes_analysis/custom_comparisons/profiles/metaprofile_nestin_mut_vs_emx1_mut.png)
*Mutant-to-mutant comparison: genotype effect under mutation*

### Individual CRE Plots (Optional)
Individual CRE plots are generated for each CRE showing detailed signal profiles. These are optional and can be skipped for faster analysis.

## Key Features

### Data Source
- **Splicing Genes**: 1,138 genes from Reactome pathways
- **CRE-Gene Links**: ENCODE Table 16 (3,256,804 published correlations)
- **Cell Type Classification**: Based on SubType column from Table 16
- **Signal Data**: GABA BigWig files (Nestin-Ctrl, Nestin-Mut, Emx1-Ctrl, Emx1-Mut)

### Analysis Focus
- **Biological Question**: How do SRF mutations affect chromatin accessibility at splicing gene regulatory elements?
- **CRE Types**: All CREs linked to splicing genes (not cell-type restricted)
- **Signal Analysis**: GABA neuron ATAC-seq signal only
- **Comparisons**: Within-genotype and cross-genotype mutation effects

### Performance Options
- **Fast Mode** (default): Metaprofiles only (3-5 minutes)
- **Full Mode**: All individual plots (30-50 minutes)
- **Parallel Processing**: 8x faster with `PARALLEL_JOBS=8`
- **DPI Control**: `INDIVIDUAL_DPI=100` for faster individual plots

## Usage

### Quick Analysis (Recommended)
```bash
# Extract CREs and create metaprofiles only
sbatch 0_RUN_SPLICING_GENES_ANALYSIS.sh
```

### Full Analysis with Individual Plots
```bash
# Complete analysis with all individual plots
SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 0_RUN_SPLICING_GENES_ANALYSIS.sh
```

### Custom Comparisons Only
```bash
# Run custom cross-genotype comparisons
sbatch 6_create_custom_comparisons.sh
```

## Biological Context

### Splicing Machinery Genes
The analysis targets genes involved in:
- mRNA splicing and processing
- Spliceosome assembly
- Alternative splicing regulation
- Gene expression fine-tuning

### Expected Patterns
- **Normal Splicing**: High accessibility at splicing gene CREs in control samples
- **SRF Mutation Effects**: Altered accessibility patterns may indicate:
  - Dysregulated splicing factor expression
  - Compensatory responses to mutations
  - Genotype-specific regulatory rewiring

## Technical Details

### Parameters
- **Window Size**: ±2000 bp around CRE center
- **Bin Size**: 50 bp per bin
- **Reference Point**: CRE center
- **Sorting**: By mean signal (descending)
- **Processors**: 16 (parallel processing)

### File Formats
- **BED6**: chr, start, end, cre_id, score, strand
- **TSV**: Tab-separated CRE-gene linkage data
- **BigWig**: Binary ATAC-seq signal tracks
- **Matrix**: deepTools compressed matrices (.gz/.tab)

## Output Structure

```
output/splicing_genes_analysis/
├── splicing_genes_CREs_all_celltypes.tsv
├── splicing_genes_CREs_GABA.tsv
├── splicing_genes_CREs_GABA_specific.tsv
├── splicing_genes_CREs_all.bed
├── splicing_genes_CREs_GABA.bed
├── heatmaps_deeptools/
│   ├── matrix_*.gz (reusable matrices)
│   ├── heatmap_*.png
│   ├── metaprofile_*.png
│   └── README.txt
├── profiles/
│   ├── metaprofile_*.png
│   ├── individual_*.png (optional)
│   └── README.txt
├── bigwig_profiles/
│   ├── metaprofile_*.png
│   ├── individual_*.png (optional)
│   ├── summary_*.tsv
│   └── README.txt
└── custom_comparisons/
    ├── matrix_*.gz
    ├── profiles/
    │   ├── metaprofile_*.png
    │   ├── individual_comparison*_*.png (optional)
    │   └── README.txt
    ├── comparison_statistics.txt
    └── README.txt
```

## Notes

- **Emx1-Ctrl Excluded**: Failed sample quality, excluded from analysis
- **Literature Data**: All CREs from published ENCODE correlations (Table 16)
- **GABA Signal Only**: Analysis uses GABA BigWig files for all CRE types
- **Reusable Matrices**: deepTools matrices can be reused for additional analyses