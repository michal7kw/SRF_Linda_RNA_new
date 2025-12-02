# Splicing Gene ENCODE cCRE Analysis Pipeline

This pipeline analyzes ATAC-seq signal at ENCODE cCREs (candidate cis-regulatory elements) linked to splicing-related genes by genomic proximity.

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
- For literature-based correlations, see `CREs_splicing_genes_paper/`

### CRE-Gene Linkage Method (Proximity-Based)

This pipeline uses **genomic proximity** to link CREs to genes:
1. Gene coordinates fetched from **Ensembl BioMart (mm10/GRCm38)**
2. **bedtools window** (+/-500kb) identifies nearby ENCODE cCREs
3. No statistical filtering (FDR/PCC) - purely distance-based
4. Result: ENCODE cCREs with **type annotations** but no correlation evidence

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
1. Load splicing genes list (1,138 genes)
2. Fetch gene coordinates from BioMart (mm10)
3. Use bedtools window (+/-500kb) to find nearby ENCODE cCREs
4. No statistical filtering - purely distance-based
5. Result: ENCODE cCREs with type annotations but no correlation evidence

## Comparison with `CREs_splicing_genes_paper` Pipeline

| Feature | `CREs_splicing_genes_paper` | `CREs_splicing_genes_encode` (this pipeline) |
|---------|------------------|----------------------------------------|
| CRE Source | Table 16 (literature) | ENCODE cCREs |
| Linkage | PCC correlation | Genomic proximity |
| Filtering | FDR < 0.05, \|PCC\| > 0.2 | Distance-based (+/- 500kb) |
| CRE Types | Literature-defined | ENCODE annotations (pELS, dELS, etc.) |

## Pipeline Steps

1. **1_extract_encode_cCREs**: Extract ENCODE cCREs near splicing genes using BioMart coordinates
2. **2_convert_to_bed**: Convert TSV files to BED format for deepTools
3. **3_create_profiles**: Create ATAC signal profiles with Ctrl vs Mut comparisons
4. **4_create_heatmaps**: Create heatmaps and metaprofiles with deepTools
5. **5_visualize_bigwig_signal**: Direct BigWig visualization for detailed analysis
6. **6_create_custom_comparisons**: Custom cross-genotype comparisons

## Quick Start

```bash
# Run the complete pipeline
./0_RUN_CREs_splicing_genes_encode_ANALYSIS.sh

# Or run individual steps
sbatch 1_extract_encode_cCREs.sh
sbatch 2_convert_to_bed.sh
sbatch 3_create_profiles.sh
sbatch 4_create_heatmaps.sh
sbatch 5_visualize_bigwig_signal.sh
sbatch 6_create_custom_comparisons.sh
```

## Output Files

### Step 1 Output (output/)
- `CREs_splicing_genes_encode_all.tsv` - All cCRE-gene associations
- `CREs_splicing_genes_encode_GABA.tsv` - Same data for GABA analysis
- `CREs_splicing_genes_encode_by_type.tsv` - Associations with CRE type info
- `SUMMARY_CREs_splicing_genes_encode.txt` - Summary statistics

### Step 2 Output (output/)
- `CREs_splicing_genes_encode_all.bed` - BED6 format for all CREs
- `CREs_splicing_genes_encode_GABA.bed` - BED6 format for GABA analysis
- `CREs_splicing_genes_encode_{type}.bed` - Type-specific BED files

### Step 3-4 Output (output/heatmaps_deeptools/)
- Heatmaps: `heatmap_*.png`
- Metaprofiles: `metaprofile_*.png`
- Signal matrices: `matrix_*.gz`, `matrix_*.tab`

### Step 5 Output (output/bigwig_profiles_minSig{X}_minFC{Y}/)
- Metaprofiles with Ctrl vs Mut comparisons:
  - `metaprofile_nestin_ctrl_vs_mut.png`
  - `metaprofile_emx1_ctrl_vs_mut.png`
- Individual CRE plots (optional):
  - `individual_Nestin_<Gene>_<CRE_ID>.png`
  - `individual_Emx1_<Gene>_<CRE_ID>.png`
- Summary tables:
  - `summary_nestin.tsv` - Statistics for all Nestin CREs
  - `summary_emx1.tsv` - Statistics for all Emx1 CREs

### Step 6 Output (output/custom_comparisons/)
- Cross-genotype comparison plots:
  - `profiles/metaprofile_nestin_ctrl_vs_nestin_mut.png`
  - `profiles/metaprofile_nestin_ctrl_vs_emx1_mut.png`
  - `profiles/metaprofile_nestin_mut_vs_emx1_mut.png`
- Signal matrices for each comparison:
  - `matrix_nestin_ctrl_vs_mut.gz` / `.tab`
  - `matrix_nestin_ctrl_vs_emx1_mut.gz` / `.tab`
  - `matrix_nestin_mut_vs_emx1_mut.gz` / `.tab`
- Note: Emx1-Ctrl is excluded (failed sample)

## Performance Options

Steps 3, 5, and 6 support performance tuning:

```bash
# Skip individual plots (fast mode - default)
sbatch 3_create_profiles.sh

# Create individual plots
SKIP_INDIVIDUAL=0 sbatch 3_create_profiles.sh

# Parallel processing for individual plots
SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 3_create_profiles.sh

# Custom DPI for individual plots
SKIP_INDIVIDUAL=0 INDIVIDUAL_DPI=100 sbatch 3_create_profiles.sh
```

## Dependencies

- Python 3.x with: pandas, numpy, matplotlib, seaborn, scipy, pyBigWig
- deepTools (computeMatrix, plotHeatmap, plotProfile)
- bedtools
- Conda environment: `rna_seq_analysis_deep`

## Input Requirements

- Splicing genes list: `/beegfs/.../CREs_splicing_genes_paper/extracted_genes_final.csv`
- ENCODE cCREs: `../data/mm10-cCREs.bed`
- BigWig files: `../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw`

## Known Issues

**CRITICAL: Directory Path Mismatch**

All scripts currently hardcode paths to `CREs_splicing_genes_encode`, but this directory is named `CREs_splicing_genes_encode`. Before running, either:
1. Rename this directory to `CREs_splicing_genes_encode`, OR
2. Update the `cd` commands in all scripts to use the correct path

Scripts affected:
- `0_RUN_CREs_splicing_genes_encode_ANALYSIS.sh`
- `1_extract_encode_cCREs.py` / `.sh`
- `2_convert_to_bed.py` / `.sh`
- `3_create_profiles.py` / `.sh`
- `4_create_heatmaps.sh`
- `5_visualize_bigwig_signal.py` / `.sh`
- `6_visualize_custom_comparisons.py`
- `6_create_custom_comparisons.sh`
