# ENCODE cCREs Analysis Pipeline

This pipeline analyzes ENCODE candidate cis-regulatory elements (cCREs) associated with splicing-related genes, comparing ATAC-seq signal patterns across different CRE types (enhancers, CTCF sites, promoters, etc.).

## Overview

The workflow identifies ENCODE cCREs linked to splicing machinery genes and visualizes ATAC-seq signal patterns at these regulatory elements. The analysis focuses on understanding how SRF mutations affect chromatin accessibility at different classes of regulatory elements.

**Key Features:**
- Uses high-confidence ENCODE cCREs from mm10-cCREs.bed
- Analyzes ALL CRE types (dELS, pELS, PLS, CTCF-only, DNase-H3K4me3, etc.)
- Creates separate visualizations for each CRE type
- Focuses on GABA neurons (same as splicing_genes pipeline)
- **NEW**: Links CREs to splicing genes using **Genomic Proximity** (BioMart coordinates + bedtools window)
- **Optimized**: Uses bedtools window for fast overlap detection

## Results Summary

Based on the latest run (BioMart + Proximity):

| Metric | Value |
|--------|-------|
| Input splicing genes | 1,138 |
| Genes with mm10 coordinates | ~960 |
| **Unique ENCODE cCREs found** | **~275,000** |
| **Unique splicing genes represented** | **~960** |
| Total CRE-gene associations | ~275,000 |

### CRE Type Distribution

| CRE Type | Count |
|----------|-------|
| dELS (distal enhancers) | ~120,000 |
| dELS,CTCF-bound | ~90,000 |
| pELS (proximal enhancers) | ~25,000 |
| pELS,CTCF-bound | ~20,000 |
| PLS (promoter-like) | ~5,000 |
| PLS,CTCF-bound | ~4,000 |
| DNase-H3K4me3 | ~400 |
| DNase-H3K4me3,CTCF-bound | ~150 |
| CTCF-only,CTCF-bound | ~15 |

## Pipeline Steps

### Step 1: Extract ENCODE cCREs (`1_extract_encode_cCREs.py`)

**Method** (BioMart + Proximity):
1. Load splicing genes list (1,138 genes)
2. **Fetch mm10 coordinates** for these genes from Ensembl BioMart (using `urllib`)
3. Use **bedtools window** (±500kb) to find ENCODE cCREs near these genes
4. Classify by cell type and CRE type

**Input**:
- Splicing genes list: `/beegfs/.../splicing_genes/extracted_genes_final.csv`
- ENCODE cCREs: `../data/mm10-cCREs.bed` (~340,000 CREs)
- **NO LONGER USES**: `../data/table_16.txt`

**Output**:
- `encode_cCREs_all_celltypes.tsv` - All ENCODE cCREs linked to splicing genes
- `encode_cCREs_GABA.tsv` - Same as above (kept for compatibility)
- `encode_cCREs_by_type.tsv` - All cCREs with type annotations
- `SUMMARY_encode_cCREs.txt` - Summary statistics

**Runtime**: ~1-2 minutes (depends on BioMart response time)

### Step 2: Convert to BED Format (`2_convert_encode_cCREs_to_bed.py`)

**Input**: TSV files from Step 1

**Output**:
- `encode_cCREs_all.bed` - All CREs
- `encode_cCREs_GABA.bed` - GABA CREs
- Type-specific BED files (dELS, pELS, PLS, etc.)

**Runtime**: ~2-3 minutes

### Step 4: Create Heatmaps and Metaprofiles (`4_create_heatmaps_encode_cCREs.sh`)

**Input**: BED files + GABA BigWig files

**Output**:
- Heatmaps and metaprofiles for ALL CREs
- Heatmaps and metaprofiles for GABA CREs
- Heatmaps and metaprofiles for each CRE type
- Genotype-specific analyses (Nestin, Emx1)

**Method**: Uses deepTools (computeMatrix + plotHeatmap/plotProfile)

**Runtime**: ~30-60 minutes

### Step 5: Custom Comparison Visualizations (`5_visualize_custom_comparisons.py`)

**Input**: deepTools matrix files from Step 4

**IMPORTANT**: Emx1-Ctrl is a **FAILED SAMPLE** - Nestin-Ctrl is used as the reference control for both Nestin-Mut and Emx1-Mut.

**Comparisons** (3 total):
1. **Nestin-Ctrl vs Nestin-Mut** - Within-genotype mutation effect
2. **Nestin-Ctrl vs Emx1-Mut** - Cross-genotype mutation effect (using Nestin-Ctrl as reference)
3. **Nestin-Mut vs Emx1-Mut** - Mutant genotype comparison

**Output**:
- `custom_comparisons/overview_all_conditions_*.png` - Combined overview of all 3 conditions
- `custom_comparisons/profiles/metaprofile_nestin_ctrl_vs_mut_*.png` - Nestin Ctrl vs Mut
- `custom_comparisons/profiles/metaprofile_nestin_ctrl_vs_emx1_mut_*.png` - Nestin-Ctrl vs Emx1-Mut
- `custom_comparisons/profiles/metaprofile_nestin_vs_emx1_mut_*.png` - Mutant genotype comparison
- `custom_comparisons/comparison_statistics_*.txt` - Statistical summaries

**Features**:
- Publication-quality metaprofiles (DPI 300)
- Difference plots showing signal changes
- Statistical tests (paired t-test)
- Optional individual CRE plots

**Method**: Python (matplotlib, seaborn, scipy)

**Runtime**: ~15-30 minutes

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
| **Unique CREs** | ~5-10 | **~275,000** |
| **Type-Specific Analysis** | No | Yes (separate analysis per type) |
| **Overlap Method** | Naive O(n×m) loops | **bedtools window** (fast) |
| **Step 1 Runtime** | ~30-45 min (timeout risk) | **~1-2 minutes** |

## Data Sources

**Splicing Genes:**
- Path: `/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/splicing_genes/extracted_genes_final.csv`
- Content: 1,138 genes from Reactome pathways (mRNA splicing, spliceosome assembly)

**ENCODE cCREs:**
- Path: `../data/mm10-cCREs.bed`
- Content: Mouse (mm10) candidate cis-regulatory elements
- Format: BED6 with type annotation in column 6

**Gene Coordinates:**
- Source: **Ensembl BioMart (mm10/GRCm38)**
- Method: Fetched dynamically via `urllib`

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
├── 1_extract_encode_cCREs.py         # Step 1: Extract cCREs (BioMart + Proximity)
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
        ├── overview_all_conditions_*.png   # Combined overview (3 conditions)
        ├── comparison_statistics_*.txt     # Statistical summaries
        └── profiles/
            ├── metaprofile_nestin_ctrl_vs_mut_*.png      # Nestin Ctrl vs Mut
            ├── metaprofile_nestin_ctrl_vs_emx1_mut_*.png # Nestin-Ctrl vs Emx1-Mut
            ├── metaprofile_nestin_vs_emx1_mut_*.png      # Mutant comparison
            └── individual_*_*.png                        # Optional individual CRE plots
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

## Troubleshooting

### Common Issues

**"bedtools not found"**
- The script automatically tries multiple conda environments
- If still failing, manually activate an environment with bedtools

**"No overlapping CREs found"**
- Check chromosome naming (chr1 vs 1)
- Verify input files exist and are not empty
- Try relaxing thresholds (PCC_THRESHOLD, FDR_THRESHOLD)

**"computeMatrix failed"**
- Check BigWig files exist and are accessible
- Verify BED files are properly formatted (no headers)
- Check memory allocation (may need >64GB for large sets)

**"Job timed out"**
- Step 1 is now optimized (~1 min), increase time only if needed
- Step 4 may need more time for many CRE types

## Dependencies

- **Python 3**: pandas, numpy
- **bedtools**: v2.31+ (for fast overlap detection)
- **deepTools**: v3.5+ (for visualization)
- **Conda environments**: `sc-chromatin2` or `rna_seq_analysis_deep`

## Notes

- **Emx1-Ctrl Sample**: This is a **FAILED SAMPLE** and is excluded from all comparisons. Nestin-Ctrl is used as the reference control for both Nestin-Mut and Emx1-Mut.
- **Runtime**: Total pipeline ~45-90 minutes (much faster than original)
- **Memory**: Peak ~64GB during Table 16 processing
- **Storage**: Output ~100-500MB depending on PNG generation

---

**Last Updated**: 2025-11-25
**Pipeline Version**: 1.2 (Added custom comparison visualizations)
**Compatible With**: splicing_genes pipeline v1.0
