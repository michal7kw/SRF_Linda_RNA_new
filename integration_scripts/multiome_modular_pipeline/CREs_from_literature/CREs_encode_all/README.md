# ENCODE cCREs Analysis Pipeline

This pipeline analyzes ENCODE candidate cis-regulatory elements (cCREs) associated with splicing-related genes, comparing ATAC-seq signal patterns across different CRE types (enhancers, CTCF sites, promoters, etc.).

## Overview

The workflow identifies ENCODE cCREs linked to splicing machinery genes and visualizes ATAC-seq signal patterns at these regulatory elements. 
The analysis focuses on understanding how SRF mutations affect chromatin accessibility at different classes of regulatory elements.

**Key Features:**
- Uses high-confidence ENCODE cCREs from **mm10-cCREs.bed**
- Analyzes ALL CRE types (dELS, pELS, PLS, CTCF-only, DNase-H3K4me3, etc.)
- Creates separate visualizations for each **CRE type**
- Focuses on **CREs_splicing_genes_paper**
- Links CREs to splicing genes using **Genomic Proximity** (BioMart coordinates + bedtools window)

> **CRITICAL**: Emx1-Ctrl is a **FAILED SAMPLE** and is excluded from all analyses.
> Only 3 conditions are analyzed: **Nestin-Ctrl** (reference control), **Nestin-Mut**, and **Emx1-Mut**.


## Results Summary

| Metric | Value |
|--------|-------|
| Input splicing genes | 1,138 |
| Genes with mm10 coordinates | ~960 |
| Unique ENCODE cCREs found | ~275,000 |
| Unique splicing genes represented | ~960 |
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

> **Note**: Step 3 is intentionally skipped in the numbering for consistency with related pipelines.

### Step 1: Extract ENCODE cCREs (`1_extract_encode_cCREs.py`)

**Method** (BioMart + Proximity):
1. Load splicing genes list (1,138 genes)
2. Fetch mm10 coordinates for these genes from Ensembl BioMart (using `urllib`)
3. Use bedtools window (±500kb) to find ENCODE cCREs near these genes
4. Classify by CRE type

**Input**:
- Splicing genes list: `/beegfs/.../CREs_splicing_genes_paper/extracted_genes_final.csv`
- ENCODE cCREs: `../data/mm10-cCREs.bed` (~340,000 CREs)

**Output** (in `output/`):
- `encode_cCREs_linked.tsv` - All ENCODE cCREs linked to splicing genes
- `encode_cCREs_by_type.tsv` - All cCREs with type annotations
- `SUMMARY_encode_cCREs.txt` - Summary statistics

### Step 2: Convert to BED Format (`2_convert_encode_cCREs_to_bed.py`)

**Input**: TSV files from Step 1

**Output** (in `output/`):
- `encode_cCREs_all.bed` - All CREs (BED6 format)
- `encode_cCREs_{type}.bed` - Type-specific BED files (dELS, pELS, PLS, CTCF_only, etc.)

### Step 4: Create Heatmaps and Metaprofiles (`4_create_heatmaps_encode_cCREs.sh`)

**Input**: BED files + GABA BigWig files

**Output** (in `output/heatmaps_deeptools/`):
- `heatmap_all_celltypes.png` / `metaprofile_all_celltypes.png` - All CREs
- `heatmap_{type}.png` / `metaprofile_{type}.png` - Per CRE type
- `heatmap_GABA_nestin.png` / `metaprofile_GABA_nestin.png` - Nestin genotype
- `heatmap_GABA_emx1.png` / `metaprofile_GABA_emx1.png` - Emx1 genotype
- `matrix_*.gz` / `matrix_*.tab` - Reusable deepTools matrices

**Parameters**: Window ±2kb, bin size 50bp, 16 processors

**Method**: Uses deepTools (computeMatrix + plotHeatmap/plotProfile)

### Step 5a: Custom Comparison Matrices (`5_create_custom_comparisons.sh`)

**Input**: BED files from Step 2 + BigWig files

> **WARNING**: Emx1-Ctrl is a **FAILED SAMPLE** and must be excluded from all analyses.
> Only Nestin-Ctrl should be used as the reference control for both Nestin-Mut and Emx1-Mut.

**Valid Comparisons** (3 total):
1. **Nestin-Ctrl vs Nestin-Mut** - Within-genotype mutation effect
2. **Nestin-Ctrl vs Emx1-Mut** - Cross-genotype mutation effect (using Nestin-Ctrl as reference)
3. **Nestin-Mut vs Emx1-Mut** - Mutant genotype comparison

**Output** (in `output/custom_comparisons/`):
- `matrix_nestin_ctrl_vs_mut.gz/.tab` - Nestin comparison matrix
- `matrix_nestin_ctrl_vs_emx1_mut.gz/.tab` - Cross-genotype comparison (Nestin-Ctrl as reference)
- `matrix_nestin_mut_vs_emx1_mut.gz/.tab` - Mutant genotype comparison
- `matrix_all_conditions.gz` - Combined 3-sample matrix (excludes Emx1-Ctrl)
- `profiles/overview_all_conditions.png` - Overview metaprofile (3 conditions)

**Method**: deepTools (computeMatrix + plotProfile)

### Step 5b: Custom Comparison Visualizations (`5_visualize_custom_comparisons.py`)

**Input**: deepTools matrix files from Step 5a (3 conditions only, excludes Emx1-Ctrl)

**Output** (in `output/custom_comparisons_minSig{X}_minFC{Y}/`):
- `overview_all_conditions_*.png` - Combined overview of all 3 conditions
- `profiles/metaprofile_nestin_ctrl_vs_mut_*.png` - Nestin Ctrl vs Mut
- `profiles/metaprofile_nestin_ctrl_vs_emx1_mut_*.png` - Nestin-Ctrl vs Emx1-Mut
- `profiles/metaprofile_nestin_vs_emx1_mut_*.png` - Mutant genotype comparison
- `profiles/individual_*.png` - Individual CRE plots (if --full mode)
- `comparison_statistics_*.txt` - Statistical summaries

**Filtering Parameters** (for individual plots):
- `MIN_SIGNAL` (default: 1.0) - Minimum max signal required
- `MIN_FC` (default: 1.5) - Minimum fold change required

**Method**: Python (matplotlib, seaborn, scipy)

## Usage

### Quick Analysis (Recommended)
```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_encode_all

# Run complete pipeline with SLURM job dependencies
./0_RUN_ENCODE_CCRES_ANALYSIS.sh
```

**Job Chain:** Step 1 (Extract) → Step 2 (Convert) → Step 4 (Heatmaps) → Step 5a (Custom matrices) → Step 5b (Comparisons)

**Expected Runtime:** ~4-6 hours total

### Individual Steps

If you need to re-run specific steps:

```bash
# Step 1: Extract cCREs (requires bedtools)
sbatch 1_extract_encode_cCREs.sh

# Step 2: Convert to BED
sbatch 2_convert_encode_cCREs_to_bed.sh

# Step 4: Create heatmaps (requires deepTools)
sbatch 4_create_heatmaps_encode_cCREs.sh

# Step 5a: Create custom comparison matrices
sbatch 5_create_custom_comparisons.sh

# Step 5b: Create custom comparison visualizations
sbatch 5_visualize_custom_comparisons.sh              # Fast mode (metaprofiles only)
sbatch 5_visualize_custom_comparisons.sh --full       # With individual CRE plots
sbatch 5_visualize_custom_comparisons.sh --all        # Use all_celltypes matrix

# Environment variables for filtering (optional)
MIN_SIGNAL=2.0 MIN_FC=2.0 sbatch 5_visualize_custom_comparisons.sh --full
```

## Output Directory Structure

```
output/
├── encode_cCREs_linked.tsv           # Step 1: All CRE-gene links
├── encode_cCREs_by_type.tsv          # Step 1: CREs with type annotations
├── SUMMARY_encode_cCREs.txt          # Step 1: Summary statistics
├── encode_cCREs_all.bed              # Step 2: All CREs (BED format)
├── encode_cCREs_{type}.bed           # Step 2: Type-specific BED files
├── heatmaps_deeptools/               # Step 4: Heatmaps and matrices
│   ├── heatmap_*.png
│   ├── metaprofile_*.png
│   └── matrix_*.gz/.tab
├── custom_comparisons/               # Step 5a: Custom comparison matrices (3 conditions)
│   ├── matrix_nestin_ctrl_vs_mut.gz/.tab
│   ├── matrix_nestin_ctrl_vs_emx1_mut.gz/.tab
│   ├── matrix_nestin_mut_vs_emx1_mut.gz/.tab
│   ├── matrix_all_conditions.gz      # 3 samples (excludes Emx1-Ctrl)
│   └── profiles/overview_all_conditions.png
└── custom_comparisons_minSig*_minFC*/  # Step 5b: Visualization outputs
    ├── overview_all_conditions_*.png   # 3 conditions only
    ├── comparison_statistics_*.txt
    └── profiles/
        ├── metaprofile_*.png
        └── individual_*.png          # (if --full mode)
```

> **Note**: Emx1-Ctrl is excluded from all analyses (failed sample). Only 3 conditions are analyzed: Nestin-Ctrl, Nestin-Mut, Emx1-Mut.

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

## Genes Analyzed

**Gene Set**: Splicing-related genes from Reactome pathways
- **Source**: `/beegfs/.../CREs_splicing_genes_paper/extracted_genes_final.csv`
- **Count**: 1,138 genes
- **Pathways**: mRNA splicing, spliceosome assembly, RNA processing

## CRE Filtering Criteria

| Filter | Value | Applied? |
|--------|-------|----------|
| **FDR threshold** | N/A | NO - proximity-based, no correlation data |
| **PCC threshold** | N/A | NO - proximity-based, no correlation data |
| **Distance window** | +/- 500kb | YES - from gene TSS |
| **Cell type filter** | N/A | NO - all ENCODE cCREs near genes |
| **ENCODE type filter** | All types | NO - includes dELS, pELS, PLS, CTCF, etc. |

**Filtering Pipeline**:
1. Fetch gene coordinates from BioMart (mm10)
2. Use bedtools window (+/-500kb) to find nearby ENCODE cCREs
3. No statistical filtering - purely distance-based
4. Result: ~275,000 CRE-gene associations

## Key Differences from CREs_splicing_genes_paper Pipeline

| Feature | CREs_splicing_genes_paper | encode_cCREs |
|---------|----------------|--------------|
| **CRE Source** | **Table 16** CRE coordinates directly | ENCODE **mm10-cCREs.bed** |
| **CRE Types** | Mixed, not classified | Classified by type (9 categories) |
| **Type-Specific Analysis** | No | Yes (separate analysis per type) |

## Data Sources

### Primary CRE Data Files

| File | Used? | Purpose |
|------|-------|---------|
| **mm10-cCREs.bed** | **YES** | Primary CRE source - ENCODE cCREs with type annotations |
| **table_16.txt** | **NO** | Not used in this pipeline |

**mm10-cCREs.bed** (ENCODE cCREs):
- Path: `../data/mm10-cCREs.bed`
- Content: ~340,000 mouse (mm10) candidate cis-regulatory elements from ENCODE
- Format: BED6 with CRE type annotation in column 6 (dELS, pELS, PLS, CTCF-only, etc.)
- Usage: CREs are linked to genes by **genomic proximity** (not literature correlations)

**table_16.txt** (Literature CRE-Gene Correlations):
- **NOT USED** in this pipeline
- For literature-based correlations, see `CREs_encode_paper_intersection/` or `CREs_splicing_genes_paper/`

### CRE-Gene Linkage Method

This pipeline uses **genomic proximity** to link CREs to genes:
1. Gene coordinates fetched from **Ensembl BioMart (mm10/GRCm38)**
2. **bedtools window** (+/-500kb) identifies nearby ENCODE cCREs
3. No statistical filtering (FDR/PCC) - purely distance-based

### Other Input Files

**Splicing Genes:**
- Path: `/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/splicing_genes/extracted_genes_final.csv`
- Content: 1,138 genes from Reactome pathways (mRNA splicing, spliceosome assembly)

**ATAC-seq Signal:**
- Path: `../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw`
- Files used: GABA_Nestin-Ctrl.bw, GABA_Nestin-Mut.bw, GABA_Emx1-Mut.bw
- **Excluded**: GABA_Emx1-Ctrl.bw (failed sample - do not use)