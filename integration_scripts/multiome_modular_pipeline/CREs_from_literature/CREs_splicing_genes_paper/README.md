# Splicing Genes CRE Analysis Pipeline

This pipeline analyzes CREs (cis-regulatory elements) associated with splicing-related genes from Reactome pathways, using literature-based CRE-gene correlations from Table 16.

## Overview

The workflow identifies CREs linked to splicing machinery genes and visualizes ATAC-seq signal patterns at these regulatory elements. The analysis focuses on understanding how SRF mutations affect chromatin accessibility at splicing gene regulatory regions.

## Pipeline Steps

### Step 1: Extract Splicing Gene CREs (`1_extract_splicing_gene_CREs.py`)
- **Input**: Splicing genes list (1,138 genes) + Table 16 (literature CRE-gene correlations)
- **Output**:
  - `splicing_genes_CREs_all_celltypes.tsv` - All CREs linked to splicing genes (94,542 links, 832 unique CREs)
  - `splicing_genes_CREs_GABA.tsv` - CREs in GABA/hippocampal cell types (22,442 links, 772 unique CREs)
  - `splicing_genes_CREs_GABA_specific.tsv` - GABA-exclusive CREs (13 links, 6 unique CREs)
- **Filters**: FDR < 0.05, |PCC| > 0.2
- **Key Feature**: Uses published literature correlations (not Signac peak-gene links)

### Step 2: Convert to BED Format (`2_convert_splicing_CREs_to_bed.py`)
- **Input**: TSV files from Step 1
- **Output**:
  - `splicing_genes_CREs_all.bed` - All CREs in BED6 format (832 CREs)
  - `splicing_genes_CREs_GABA.bed` - GABA CREs in BED6 format (772 CREs)
- **Purpose**: Prepare CRE coordinates for deepTools analysis

### Step 3: Create ATAC Signal Profiles (`3_create_splicing_profiles.sh`)
- **Input**: BED files + GABA BigWig files
- **Output** (in `output/heatmaps_deeptools/profiles_minSig{X}_minFC{Y}/`):
  - Metaprofiles: `metaprofile_nestin_ctrl_vs_mut.png`, `metaprofile_emx1_ctrl_vs_mut.png`
  - Individual CRE plots: `individual_Nestin_<gene>_<cre_id>.png`, `individual_Emx1_<gene>_<cre_id>.png`
  - Signal matrices: `matrix_GABA_nestin.gz`, `matrix_GABA_emx1.gz`
- **Performance**: Default skips individual plots (3-5 min), full mode creates all plots
- **Usage**: `sbatch 3_create_splicing_profiles.sh` (fast) or `SKIP_INDIVIDUAL=0 sbatch 3_create_splicing_profiles.sh` (full)

### Step 4: Create Traditional Heatmaps (`4_create_heatmaps_splicing_genes.sh`)
- **Input**: Same as Step 3
- **Output** (in `output/heatmaps_deeptools/`):
  - Heatmaps: `heatmap_all_celltypes.png`, `heatmap_GABA.png`
  - Metaprofiles: `metaprofile_all_celltypes.png`, `metaprofile_GABA.png`
  - Genotype-specific: `heatmap_GABA_nestin.png`, `heatmap_GABA_emx1.png`
- **Method**: Uses deepTools plotHeatmap/plotProfile commands

### Step 5: Visualize BigWig Signal (`5_visualize_bigwig_signal.sh`)
- **Input**: Direct BigWig file reading
- **Output** (in `output/bigwig_profiles_minSig{X}_minFC{Y}/`):
  - Metaprofiles: `metaprofile_nestin_ctrl_vs_mut.png`, `metaprofile_emx1_ctrl_vs_mut.png`
  - Individual plots: `individual_Nestin_<gene>_<cre_id>.png`, `individual_Emx1_<gene>_<cre_id>.png`
  - Statistics: `summary_nestin.tsv`, `summary_emx1.tsv`
- **Method**: Direct BigWig reading with pyBigWig (alternative to deepTools)

### Step 6: Create Custom Comparisons (`6_create_custom_comparisons.sh`)
- **Input**: ALL CREs with custom genotype comparisons
- **Output** (in `output/custom_comparisons_minSig{X}_minFC{Y}/profiles/`):
  - Custom metaprofiles: `metaprofile_nestin_ctrl_vs_nestin_mut.png`, `metaprofile_nestin_ctrl_vs_emx1_mut.png`, `metaprofile_nestin_mut_vs_emx1_mut.png`
  - Individual plots: `individual_comparison1_<gene>_<cre_id>.png`, `individual_comparison2_<gene>_<cre_id>.png`, `individual_comparison3_<gene>_<cre_id>.png`
  - Statistics: `comparison_statistics.txt`
- **Comparisons**:
  1. Nestin-Ctrl vs Nestin-Mut (within-genotype mutation effect)
  2. Nestin-Ctrl vs Emx1-Mut (cross-genotype comparison)
  3. Nestin-Mut vs Emx1-Mut (mutant-to-mutant comparison)

## Analysis Results

### CRE Statistics (from actual analysis)
| Category | CRE-Gene Links | Unique CREs | Unique Genes | Coverage |
|----------|---------------|-------------|--------------|----------|
| All Cell Types | 94,542 | 832 | 819 | 72.0% |
| GABA Cell Types | 22,442 | 772 | 768 | 67.5% |
| GABA-Specific | 13 | 6 | 6 | 0.5% |

## Data Sources

### Primary CRE Data Files

| File | Used? | Purpose |
|------|-------|---------|
| **mm10-cCREs.bed** | **NO** | Not used in this pipeline |
| **table_16.txt** | **YES** | Primary CRE source - literature CRE-gene correlations |

**table_16.txt** (Literature CRE-Gene Correlations):
- Path: `../data/table_16.txt`
- Content: 3,256,804 published CRE-gene correlations from ENCODE consortium
- Columns: Coordinate1 (CRE coordinates), cCRE1, Gene, SubType, PCC, FDR, etc.
- Usage: **Primary CRE source** - CRE coordinates extracted directly from Coordinate1 column
- Statistical filters: FDR < 0.05, |PCC| > 0.2

### CRE-Gene Linkage Method (Literature-Based)

This pipeline uses **ONLY table_16.txt**:
1. Filter Table 16 for splicing genes with significant correlations
2. Extract CRE coordinates directly from Coordinate1 column (chr_start_end format)
3. Cell type classification based on SubType column
4. Result: CREs with **literature-validated** gene links (no ENCODE type annotations)

### Other Data Sources
- **Signal Data**: GABA BigWig files (Nestin-Ctrl, Nestin-Mut, Emx1-Mut)
- **Note**: Emx1-Ctrl is not used (Nestin-Ctrl serves as control for Emx1 comparisons)

## Genes Analyzed

**Gene Set**: Splicing-related genes from Reactome pathways
- **Source**: `/beegfs/.../CREs_splicing_genes_paper/extracted_genes_final.csv`
- **Count**: 1,138 genes
- **Pathways**: mRNA splicing, spliceosome assembly, RNA processing

## CRE Filtering Criteria

| Filter | Value | Applied? |
|--------|-------|----------|
| **FDR threshold** | < 0.05 | YES - Benjamini-Hochberg correction |
| **PCC threshold** | > 0.2 | YES - Pearson correlation coefficient |
| **Distance window** | N/A | NO - uses Table 16 correlations |
| **Cell type filter** | All + GABA subset | YES - creates separate GABA output |
| **ENCODE type filter** | N/A | NO - no ENCODE cCREs used |
| **ENCODE intersect** | N/A | NO - uses Table 16 coordinates directly |

**Cell Type Keywords** (for GABA subset):
- Hippocampal: CA1, CA2, CA3, DG, DGNBL, GRC
- GABAergic: LAMP5, LAMP, VIP, SST, PV, PVGA, SSTGA, VIPGA, LAMGA, GABA, INH

**Filtering Pipeline**:
1. Load splicing genes list (1,138 genes)
2. Filter Table 16 for splicing genes
3. Apply statistical filters (FDR < 0.05, |PCC| > 0.2)
4. Extract CRE coordinates from Coordinate1 column
5. Create GABA subset by filtering SubType column
6. Result: 832 CREs (all), 772 GABA CREs, 6 GABA-specific CREs

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

### Filter Thresholds
```bash
# Run with specific signal/fold-change thresholds
MIN_SIGNAL=2.0 MIN_FC=2.0 sbatch 3_create_splicing_profiles.sh
```

## Output Structure

```
output/
├── splicing_genes_CREs_all_celltypes.tsv    # All CRE-gene links (94,542 rows)
├── splicing_genes_CREs_GABA.tsv             # GABA CRE-gene links (22,442 rows)
├── splicing_genes_CREs_GABA_specific.tsv    # GABA-specific links (13 rows)
├── splicing_genes_CREs_all.bed              # All CREs (832 regions)
├── splicing_genes_CREs_GABA.bed             # GABA CREs (772 regions)
├── SUMMARY_splicing_genes_CREs.txt          # Summary statistics
├── heatmaps_deeptools/
│   ├── matrix_*.gz (reusable matrices)
│   ├── heatmap_*.png
│   ├── metaprofile_*.png
│   ├── profiles_minSig{X}_minFC{Y}/
│   │   ├── metaprofile_*.png
│   │   └── individual_*.png (optional)
│   └── README.txt
├── bigwig_profiles_minSig{X}_minFC{Y}/
│   ├── metaprofile_*.png
│   ├── individual_*.png (optional)
│   ├── summary_*.tsv
│   └── README.txt
└── custom_comparisons_minSig{X}_minFC{Y}/
    ├── matrix_*.gz
    ├── profiles/
    │   ├── metaprofile_*.png
    │   └── individual_comparison*_*.png (optional)
    ├── comparison_statistics.txt
    └── README.txt
```

## Notes

- **Emx1-Ctrl Handling**: Emx1-Ctrl sample not used; Nestin-Ctrl serves as control for Emx1 comparisons
- **Literature Data**: All CREs from published correlations Table 16
- **GABA Signal Only**: Analysis uses GABA BigWig files for all CRE types
- **Reusable Matrices**: deepTools matrices can be reused for additional analyses
- **Parameterized Output**: Output directories include filter thresholds (e.g., `minSig2.0_minFC2.0`)
