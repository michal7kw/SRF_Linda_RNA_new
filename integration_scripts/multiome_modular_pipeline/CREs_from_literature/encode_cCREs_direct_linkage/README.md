# ENCODE cCREs Analysis Pipeline

This pipeline analyzes ENCODE candidate cis-regulatory elements (cCREs) associated with splicing-related genes, comparing ATAC-seq signal patterns across different CRE types (enhancers, CTCF sites, promoters, etc.).

## Overview

The workflow identifies ENCODE cCREs linked to splicing machinery genes and visualizes ATAC-seq signal patterns at these regulatory elements. The analysis focuses on understanding how SRF mutations affect chromatin accessibility at different classes of regulatory elements.

**Key Features:**
- Uses high-confidence ENCODE cCREs from mm10-cCREs.bed
- Analyzes ALL CRE types (dELS, pELS, PLS, CTCF-only, DNase-H3K4me3, etc.)
- Creates separate visualizations for each CRE type
- Focuses on GABA neurons (same as splicing_genes pipeline)
- Links CREs to splicing genes using literature correlations (Table 16)
- **Optimized**: Uses bedtools intersect for fast overlap detection

## Results Summary

Based on the latest run (2025-11-25):

| Metric | Value |
|--------|-------|
| Input splicing genes | 1,138 |
| Table 16 links (after filtering) | 94,542 |
| **Unique ENCODE cCREs found** | **2,780** |
| **Unique splicing genes represented** | **819** |
| Total CRE-gene associations | 348,126 |
| GABA-specific cCREs | 2,605 |
| GABA-specific genes | 763 |

### CRE Type Distribution

| CRE Type | All Cell Types | GABA Only |
|----------|---------------|-----------|
| dELS (distal enhancers) | 153,609 | 32,742 |
| dELS,CTCF-bound | 125,431 | 26,575 |
| pELS (proximal enhancers) | 30,546 | 6,383 |
| pELS,CTCF-bound | 25,472 | 5,323 |
| PLS (promoter-like) | 6,643 | 1,296 |
| PLS,CTCF-bound | 5,690 | 1,235 |
| DNase-H3K4me3 | 513 | 97 |
| DNase-H3K4me3,CTCF-bound | 203 | 27 |
| CTCF-only,CTCF-bound | 19 | 0 |

### Top Genes by CRE Count

| Gene | Associated cCREs |
|------|-----------------|
| Atxn1 | 8,830 |
| Celf4 | 6,936 |
| Ptbp1 | 4,014 |
| Cdk9 | 3,416 |
| Tsr1 | 2,856 |
| Ago2 | 2,778 |
| Atxn7l3 | 2,400 |
| Hnrnpul2 | 2,358 |
| Khdrbs3 | 2,240 |
| Srsf2 | 2,196 |

## Pipeline Steps

### Step 1: Extract ENCODE cCREs (`1_extract_encode_cCREs.py`)

**Method** (Optimized):
1. Load splicing genes list (1,138 genes)
2. Load Table 16 in chunks (500K rows at a time) - filters immediately
3. Use **bedtools intersect** to find overlapping ENCODE cCREs (fast!)
4. Classify by cell type and CRE type

**Input**:
- Splicing genes list: `/beegfs/.../splicing_genes/extracted_genes_final.csv`
- ENCODE cCREs: `../data/mm10-cCREs.bed` (~340,000 CREs)
- Table 16: `../data/table_16.txt` (3.2M correlations, 567 MB)

**Output**:
- `encode_cCREs_all_celltypes.tsv` - All ENCODE cCREs linked to splicing genes (348,126 links)
- `encode_cCREs_GABA.tsv` - ENCODE cCREs in GABA/hippocampal cell types (73,678 links)
- `encode_cCREs_by_type.tsv` - All cCREs with type annotations
- `SUMMARY_encode_cCREs.txt` - Summary statistics

### Step 2: Convert to BED Format (`2_convert_encode_cCREs_to_bed.py`)

**Input**: TSV files from Step 1

**Output**:
- `encode_cCREs_all.bed` - All CREs (2,780 unique)
- `encode_cCREs_GABA.bed` - GABA CREs (2,605 unique)
- Type-specific BED files:
  - `encode_cCREs_dELS.bed` (1,257 CREs)
  - `encode_cCREs_dELS_CTCF_bound.bed` (941 CREs)
  - `encode_cCREs_pELS.bed` (281 CREs)
  - `encode_cCREs_pELS_CTCF_bound.bed` (198 CREs)
  - `encode_cCREs_PLS.bed` (41 CREs)
  - `encode_cCREs_PLS_CTCF_bound.bed` (51 CREs)
  - `encode_cCREs_DNase_H3K4me3.bed` (8 CREs)
  - `encode_cCREs_DNase_H3K4me3_CTCF_bound.bed` (2 CREs)
  - `encode_cCREs_CTCF_only_CTCF_bound.bed` (1 CRE)

### Step 4: Create Heatmaps and Metaprofiles (`4_create_heatmaps_encode_cCREs.sh`)

**Input**: BED files + GABA BigWig files

**Output**:
- Heatmaps and metaprofiles for ALL CREs
- Heatmaps and metaprofiles for GABA CREs
- Heatmaps and metaprofiles for each CRE type
- Genotype-specific analyses (Nestin, Emx1)

### Step 5: Custom Comparison Visualizations (`5_visualize_custom_comparisons.py`)

**Input**: deepTools matrix files from Step 4

**Output**:
- `custom_comparisons/overview_all_conditions_*.png` - Combined overview of all 4 conditions
- `custom_comparisons/profiles/metaprofile_nestin_ctrl_vs_mut_*.png` - Nestin Ctrl vs Mut
- `custom_comparisons/profiles/metaprofile_emx1_ctrl_vs_mut_*.png` - Emx1 Ctrl vs Mut
- `custom_comparisons/profiles/metaprofile_nestin_vs_emx1_ctrl_*.png` - Baseline genotype comparison
- `custom_comparisons/profiles/metaprofile_nestin_vs_emx1_mut_*.png` - Mutant genotype comparison
- `custom_comparisons/comparison_statistics_*.txt` - Statistical summaries

**Features**:
- Publication-quality metaprofiles (DPI 300)
- Difference plots showing signal changes
- Statistical tests (paired t-test)
- Optional individual CRE plots

## Usage

### Quick Analysis (Recommended)
```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/encode_cCREs

# Run complete pipeline with SLURM job dependencies
./0_RUN_ENCODE_CCRES_ANALYSIS.sh
```

This automatically runs:
1. Extract ENCODE cCREs (~1 min)
2. Convert to BED format (~2-3 min)
4. Create heatmaps for all types (~30-60 min)
5. Create custom comparison visualizations (~15-30 min)

**Total runtime: ~60-120 minutes**

### Individual Steps

If you need to re-run specific steps:

```bash
# Step 1: Extract cCREs (requires bedtools)
sbatch 1_extract_encode_cCREs.sh

# Step 2: Convert to BED
sbatch 2_convert_encode_cCREs_to_bed.sh

# Step 4: Create heatmaps (requires deepTools)
sbatch 4_create_heatmaps_encode_cCREs.sh

# Step 5: Create custom comparison visualizations
sbatch 5_visualize_custom_comparisons.sh              # Fast mode (metaprofiles only)
sbatch 5_visualize_custom_comparisons.sh --full       # With individual CRE plots
sbatch 5_visualize_custom_comparisons.sh --all        # Use all_celltypes matrix
```

## CRE Type Definitions

The ENCODE cCRE classifications used in this analysis:

| Type | Full Name | Description |
|------|-----------|-------------|
| **dELS** | Distal Enhancer-Like Signature | Enhancers >2kb from TSS |
| **pELS** | Proximal Enhancer-Like Signature | Enhancers within 2kb of TSS |
| **PLS** | Promoter-Like Signature | Active promoter regions |
| **CTCF-only** | CTCF-only | CTCF binding without enhancer activity |
| **DNase-H3K4me3** | DNase with H3K4me3 | Promoter-like with active histone marks |
| **CTCF-bound** | CTCF-bound suffix | Additional CTCF binding at the element |

## Key Differences from splicing_genes Pipeline

| Feature | splicing_genes | encode_cCREs |
|---------|----------------|--------------|
| **CRE Source** | Table 16 CRE coordinates directly | ENCODE mm10-cCREs.bed |
| **CRE Types** | Mixed, not classified | Classified by type (9 categories) |
| **Unique CREs** | ~5-10 | **2,780** |
| **Type-Specific Analysis** | No | Yes (separate analysis per type) |
| **Overlap Method** | Naive O(n×m) loops | **bedtools intersect** (fast) |

## Data Sources

**Splicing Genes:**
- Path: `/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/splicing_genes/extracted_genes_final.csv`
- Content: 1,138 genes from Reactome pathways (mRNA splicing, spliceosome assembly)

**ENCODE cCREs:**
- Path: `../data/mm10-cCREs.bed`
- Content: Mouse (mm10) candidate cis-regulatory elements
- Format: BED6 with type annotation in column 6

**CRE-Gene Links:**
- Path: `../data/table_16.txt`
- Content: 3,256,804 published CRE-gene correlations
- Statistical filters applied: FDR < 0.05, |PCC| > 0.2

**ATAC-seq Signal:**
- Path: `../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw`
- Files: GABA_Nestin-Ctrl.bw, GABA_Nestin-Mut.bw, GABA_Emx1-Ctrl.bw, GABA_Emx1-Mut.bw

## Technical Details

### SLURM Resources

| Step | Time | CPUs | Memory |
|------|------|------|--------|
| Step 1 | 2:00:00 | 16 | 64GB |
| Step 2 | 0:30:00 | 2 | 16GB |
| Step 4 | 4:00:00 | 16 | 64GB |
| Step 5 | 2:00:00 | 8 | 32GB |

### deepTools Parameters
- **Window Size**: ±2000 bp around CRE center
- **Bin Size**: 50 bp per bin
- **Reference Point**: CRE center
- **Sorting**: By mean signal (descending)
- **DPI**: 300 (publication quality)

### Statistical Filters
- **FDR Threshold**: < 0.05 (Benjamini-Hochberg correction)
- **PCC Threshold**: > 0.2 (Pearson correlation coefficient)

## Output Structure

```
encode_cCREs/
├── 0_RUN_ENCODE_CCRES_ANALYSIS.sh    # Master script
├── 1_extract_encode_cCREs.py         # Step 1: Extract cCREs (optimized)
├── 1_extract_encode_cCREs.sh         # SLURM wrapper
├── 2_convert_encode_cCREs_to_bed.py  # Step 2: Convert to BED
├── 2_convert_encode_cCREs_to_bed.sh  # SLURM wrapper
├── 4_create_heatmaps_encode_cCREs.sh # Step 4: Heatmaps
├── 5_visualize_custom_comparisons.py # Step 5: Custom comparisons
├── 5_visualize_custom_comparisons.sh # SLURM wrapper
├── README.md                          # This file
├── logs/                              # SLURM logs
│   ├── 1_extract_encode_cCREs.log
│   ├── 2_convert_encode_cCREs_to_bed.log
│   ├── 4_create_heatmaps_encode_cCREs.log
│   └── 5_visualize_custom_comparisons.log
└── output/
    ├── encode_cCREs_all_celltypes.tsv
    ├── encode_cCREs_GABA.tsv
    ├── encode_cCREs_by_type.tsv
    ├── SUMMARY_encode_cCREs.txt
    ├── encode_cCREs_all.bed
    ├── encode_cCREs_GABA.bed
    ├── encode_cCREs_{type}.bed        # One per CRE type
    ├── heatmaps_deeptools/
    │   ├── matrix_*.gz                 # Reusable matrices
    │   ├── matrix_*.tab                # Tab-delimited matrices
    │   ├── heatmap_*.png               # Heatmaps
    │   ├── metaprofile_*.png           # Metaprofiles
    │   └── computeMatrix_*.log         # deepTools logs
    └── custom_comparisons/             # Step 5 output
        ├── overview_all_conditions_*.png   # Combined overview
        ├── comparison_statistics_*.txt     # Statistical summaries
        └── profiles/
            ├── metaprofile_nestin_ctrl_vs_mut_*.png
            ├── metaprofile_emx1_ctrl_vs_mut_*.png
            ├── metaprofile_nestin_vs_emx1_ctrl_*.png
            ├── metaprofile_nestin_vs_emx1_mut_*.png
            └── individual_*_*.png          # Optional individual CRE plots
```

## Biological Interpretation

### CRE Type Insights

**dELS (Distal Enhancers) - 80% of all links:**
- Most abundant CRE type in splicing gene regulation
- Cell-type-specific regulation expected
- Key targets for understanding SRF mutation effects
- CTCF-bound variants suggest insulator/enhancer dual function

**pELS (Proximal Enhancers) - 16% of all links:**
- Close to promoters, important for basal expression
- More stable than dELS
- Direct regulation of nearby splicing genes

**PLS (Promoter-Like Signatures) - 4% of all links:**
- Mark active transcription start sites
- Direct correlation with splicing gene expression
- Changes indicate transcriptional effects

**CTCF Sites:**
- Organize chromatin 3D structure
- Should be relatively stable across conditions
- Changes indicate chromatin reorganization

### Expected Patterns

- **Normal Function**: High accessibility at enhancers (dELS/pELS) in control samples
- **SRF Mutation Effects**:
  - Altered enhancer accessibility → dysregulated splicing factor expression
  - CTCF site stability → maintained chromatin architecture
  - Promoter changes → direct transcriptional effects
- **Cell-Type Specificity**: GABA-specific enhancers vs ubiquitous promoters