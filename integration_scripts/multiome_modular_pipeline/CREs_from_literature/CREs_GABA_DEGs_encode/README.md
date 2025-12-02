# GABA DEG ENCODE cCRE Analysis Pipeline

This pipeline analyzes ENCODE candidate cis-regulatory elements (cCREs) associated with GABA differentially expressed genes (DEGs), comparing ATAC-seq signal patterns across different conditions.

## Overview

The workflow identifies ENCODE cCREs linked to GABA-specific up-regulated and down-regulated genes using literature CRE-gene correlations (Table 16), then visualizes ATAC-seq signal patterns at these regulatory elements.

**Key Features:**
- Uses high-confidence ENCODE cCREs from mm10-cCREs.bed
- Analyzes UP-regulated and DOWN-regulated DEGs separately
- Creates separate visualizations for each DEG set
- Focuses on GABA neurons (GABA BigWig files only)
- Links CREs to DEGs using literature correlations (Table 16)
- **Optimized**: Uses bedtools intersect for fast overlap detection

## Comparisons

**IMPORTANT**: Emx1-Ctrl is a **FAILED SAMPLE** and is excluded from all analyses. Nestin-Ctrl is used as the reference control for both Nestin-Mut and Emx1-Mut.

We analyze 3 comparisons for each DEG set (up/down):

| Comparison | Description |
|------------|-------------|
| Nestin-Ctrl vs Nestin-Mut | Within-genotype mutation effect (Nestin) |
| Nestin-Ctrl vs Emx1-Mut | Cross-genotype mutation effect (Nestin-Ctrl as reference) |
| Nestin-Mut vs Emx1-Mut | Genotype comparison (Mutants only) |

## Pipeline Steps

### Step 1: Extract ENCODE cCREs (`1_extract_encode_cCREs_for_DEGs.py`)

**Method** (Optimized):
1. Load GABA DEG lists (up-regulated: ~5,851 genes, down-regulated: ~797 genes)
2. Query Table 16 for CRE-gene links with statistical filters
3. Filter for GABA/hippocampal cell types
4. Use **bedtools intersect** to find overlapping ENCODE cCREs
5. Create separate output files for up/down DEGs

**Input**:
- GABA DEG lists: `Cluster_GABA_vs_Rest_up_significant.csv`, `Cluster_GABA_vs_Rest_down_significant.csv`
- ENCODE cCREs: `../data/mm10-cCREs.bed`
- Table 16: `../data/table_16.txt`

**Output** (in `output/`):
- `encode_cCREs_up_DEGs.tsv` - ENCODE cCREs linked to up-regulated DEGs
- `encode_cCREs_up_DEGs.bed` - BED format for deepTools
- `encode_cCREs_down_DEGs.tsv` - ENCODE cCREs linked to down-regulated DEGs
- `encode_cCREs_down_DEGs.bed` - BED format for deepTools
- `SUMMARY_GABA_DEGs_encode_cCREs.txt` - Summary statistics

**Conda environment**: `sc-chromatin2` (pandas, bedtools)

**Runtime**: ~2-5 minutes

### Step 2: Compute Signal Matrices (`2_compute_signal_matrices.sh`)

**Input**: BED files from Step 1 + GABA BigWig files

**Output** (in `output/heatmaps_deeptools/`):
- `matrix_up_DEGs.gz/.tab` - Matrix for up-DEG ENCODE cCREs (all 3 conditions)
- `matrix_down_DEGs.gz/.tab` - Matrix for down-DEG ENCODE cCREs (all 3 conditions)
- `matrix_*_nestin.gz/.tab` - Nestin-only matrices (Ctrl vs Mut)
- `heatmap_*.png` - Heatmaps for each DEG set
- `metaprofile_*.png` - Metaprofiles for each DEG set

**Conda environment**: `rna_seq_analysis_deep` (deepTools)

### Step 3: Custom Comparison Visualizations (`3_visualize_deg_comparisons.py`)

**Input**: deepTools matrix files from Step 2

**Output** (in `custom_comparisons_minSig{X}_minFC{Y}/`):
- `overview_all_conditions_*.png` - Combined overview
- `profiles/metaprofile_*_up_DEGs.png` - Up-DEG comparisons
- `profiles/metaprofile_*_down_DEGs.png` - Down-DEG comparisons
- `profiles/metaprofile_up_vs_down_*.png` - Up vs Down comparison
- `profiles/individual_*.png` - Individual CRE plots (filtered by signal and FC)
- `comparison_statistics.txt` - Statistical summaries

**Features**:
- Publication-quality metaprofiles (DPI 300)
- Difference plots showing signal changes
- Statistical tests (paired t-test)
- Optional individual CRE plots with filtering:
  - Minimum signal threshold (default 1.0)
  - Minimum fold change threshold (default 1.5)
  - Maximum plots per comparison (default 50)

**Conda environment**: `sc-chromatin2` (matplotlib, seaborn, scipy)

## Usage

### Quick Analysis (Recommended)
```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_GABA_DEGs_encode

# Run complete pipeline with SLURM job dependencies
./0_RUN_GABA_DEG_ENCODE_ANALYSIS.sh

# Or with options
./0_RUN_GABA_DEG_ENCODE_ANALYSIS.sh --skip-individual  # Fast mode
./0_RUN_GABA_DEG_ENCODE_ANALYSIS.sh --parallel 8       # Parallel individual plots
./0_RUN_GABA_DEG_ENCODE_ANALYSIS.sh --dry-run          # Show what would run
```

**Total runtime: ~60-120 minutes**

### Step 3 Visualization Options

The visualization script supports additional filtering options:

| Option | Default | Description |
|--------|---------|-------------|
| `--skip-individual` | False | Skip individual CRE plots (fast mode) |
| `--parallel N` | 1 | Number of parallel processes |
| `--individual-dpi` | 150 | DPI for individual plots (metaprofiles always 300) |
| `--max-individual` | 50 | Maximum individual plots per comparison |
| `--min-signal` | 1.0 | Minimum max signal required for individual plots |
| `--min-fc` | 1.5 | Minimum fold change required for individual plots |

### Individual Steps

If you need to re-run specific steps:

```bash
# Step 1: Extract ENCODE cCREs (requires bedtools)
sbatch 1_extract_encode_cCREs_for_DEGs.sh

# Step 2: Compute signal matrices (requires deepTools)
sbatch 2_compute_signal_matrices.sh

# Step 3: Create visualizations
sbatch 3_visualize_deg_comparisons.sh
```

## Data Sources

### Primary CRE Data Files

| File | Used? | Purpose |
|------|-------|---------|
| **mm10-cCREs.bed** | **YES** | ENCODE cCREs - provides CRE type annotations |
| **table_16.txt** | **YES** | Literature CRE-gene correlations - provides DEG-CRE linkages |

**mm10-cCREs.bed** (ENCODE cCREs):
- Path: `../data/mm10-cCREs.bed`
- Content: ~340,000 mouse (mm10) candidate cis-regulatory elements from ENCODE
- Format: BED6 with CRE type annotation in column 6 (dELS, pELS, PLS, CTCF-only, etc.)
- Usage: **bedtools intersect** identifies ENCODE cCREs that overlap with Table 16 regions

**table_16.txt** (Literature CRE-Gene Correlations):
- Path: `../data/table_16.txt`
- Content: 3,256,804 published CRE-gene correlations from ENCODE consortium
- Columns: Coordinate1, cCRE1, Gene, SubType, PCC, FDR, etc.
- Usage: Links GABA DEGs to CRE regions with correlation evidence
- Statistical filters: FDR < 0.05, |PCC| > 0.2
- Cell type filter: GABA/hippocampal SubTypes only

### CRE-Gene Linkage Method (Hybrid Approach)

This pipeline uses **BOTH** data sources:
1. Filter Table 16 for GABA DEGs (up/down) with significant correlations
2. Apply GABA/hippocampal cell type filter on SubType column
3. Create BED file from Table 16 genomic coordinates
4. Use **bedtools intersect** to find ENCODE cCREs overlapping these regions
5. Result: ENCODE cCREs linked to DEGs with **type annotations** and **literature evidence**

### Other Input Files

**GABA DEGs:**
- Up-regulated: `DEGs_cell_type_L1FC_0_25/biomarkers/sig_deg_lists/GABA/Cluster_GABA_vs_Rest_up_significant.csv`
- Down-regulated: `DEGs_cell_type_L1FC_0_25/biomarkers/sig_deg_lists/GABA/Cluster_GABA_vs_Rest_down_significant.csv`

**ATAC-seq Signal:**
- Path: `../../signac_results_L1/bigwig_tracks_L1/by_celltype/`
- Files: `GABA_Nestin-Ctrl.bw`, `GABA_Nestin-Mut.bw`, `GABA_Emx1-Mut.bw`
- Note: Emx1-Ctrl excluded due to quality issues

## Genes Analyzed

**Gene Set**: GABA-specific differentially expressed genes (DEGs)
- **Up-regulated DEGs**: ~5,851 genes enriched in GABA neurons (GABA markers)
- **Down-regulated DEGs**: ~797 genes depleted in GABA neurons
- **Source**: DEG analysis comparing GABA vs Rest cell types
- **Analysis**: Separate processing for up/down DEGs

## CRE Filtering Criteria

| Filter | Value | Applied? |
|--------|-------|----------|
| **FDR threshold** | < 0.05 | YES - Benjamini-Hochberg correction |
| **PCC threshold** | > 0.2 | YES - Pearson correlation coefficient |
| **Distance window** | N/A | NO - uses Table 16 correlations |
| **Cell type filter** | GABA/Hippocampal | YES - SubType keywords |
| **ENCODE type filter** | All types | NO - includes dELS, pELS, PLS, CTCF, etc. |
| **ENCODE intersect** | Required overlap | YES - CREs must overlap ENCODE cCREs |

**Cell Type Keywords** (SubType filter):
- Hippocampal: Hippocampus, CA1, CA2, CA3, CA4, DG, DGNBL, GRC
- GABAergic: GABA, GABAergic, Interneuron, PV, SST, VIP, LAMP5, LAMP, PVGA, SSTGA, VIPGA, LAMGA, INH

**Filtering Pipeline**:
1. Filter Table 16 for GABA DEGs (up or down)
2. Apply statistical filters (FDR < 0.05, |PCC| > 0.2)
3. Apply GABA/hippocampal cell type filter on SubType
4. Create BED file from Table 16 coordinates
5. bedtools intersect with mm10-cCREs.bed
6. Result: ENCODE cCREs linked to DEGs with type annotations

## Technical Details

### SLURM Resources

| Step | Time | CPUs | Memory |
|------|------|------|--------|
| Step 1 | 2:00:00 | 8 | 64GB |
| Step 2 | 4:00:00 | 16 | 64GB |
| Step 3 | 2:00:00 | 8 | 32GB |

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
CREs_GABA_DEGs_encode/
├── 0_RUN_GABA_DEG_ENCODE_ANALYSIS.sh  # Master script
├── 1_extract_encode_cCREs_for_DEGs.py # Step 1: Extract ENCODE cCREs
├── 1_extract_encode_cCREs_for_DEGs.sh # Step 1: SLURM wrapper
├── 2_compute_signal_matrices.sh       # Step 2: Compute matrices
├── 3_visualize_deg_comparisons.py     # Step 3: Visualizations
├── 3_visualize_deg_comparisons.sh     # Step 3: SLURM wrapper
├── README.md                          # This file
├── logs/                              # SLURM logs
│   ├── 1_extract_encode_cCREs_for_DEGs.log/.err
│   ├── 2_compute_signal_matrices.log/.err
│   └── 3_visualize_deg_comparisons.log/.err
└── output/
    ├── encode_cCREs_up_DEGs.tsv       # Up-DEG links
    ├── encode_cCREs_up_DEGs.bed       # Up-DEG BED
    ├── encode_cCREs_down_DEGs.tsv     # Down-DEG links
    ├── encode_cCREs_down_DEGs.bed     # Down-DEG BED
    ├── SUMMARY_GABA_DEGs_encode_cCREs.txt
    ├── heatmaps_deeptools/
    │   ├── matrix_up_DEGs.gz/.tab     # Up-DEG matrix
    │   ├── matrix_down_DEGs.gz/.tab   # Down-DEG matrix
    │   ├── matrix_*_nestin.gz/.tab    # Nestin-only matrices
    │   ├── heatmap_*.png              # Heatmaps
    │   └── metaprofile_*.png          # Metaprofiles
    └── custom_comparisons_minSig{X}_minFC{Y}/  # X=min_signal, Y=min_fc
        ├── overview_all_conditions_up_DEGs.png
        ├── overview_all_conditions_down_DEGs.png
        ├── comparison_statistics.txt
        └── profiles/
            ├── metaprofile_nestin_ctrl_vs_mut_up_DEGs.png
            ├── metaprofile_nestin_ctrl_vs_mut_down_DEGs.png
            ├── metaprofile_nestin_ctrl_vs_emx1_mut_up_DEGs.png
            ├── metaprofile_nestin_ctrl_vs_emx1_mut_down_DEGs.png
            ├── metaprofile_nestin_mut_vs_emx1_mut_up_DEGs.png
            ├── metaprofile_nestin_mut_vs_emx1_mut_down_DEGs.png
            ├── metaprofile_up_vs_down_Nestin_Ctrl.png
            ├── metaprofile_up_vs_down_Nestin_Mut.png
            ├── metaprofile_up_vs_down_Emx1_Mut.png
            └── individual_*.png       # Individual CRE plots (if not skipped)
```

**Note**: The custom_comparisons directory name includes the min_signal and min_fc parameters used (e.g., `custom_comparisons_minSig1.0_minFC1.5/`).

## Related Pipelines

| Pipeline | Description |
|----------|-------------|
| `CREs_GABA_DEGs_paper/` | Table 16 CREs directly (not ENCODE filtered) |
| `CREs_encode_GABA_specific/` | GABA-specific ENCODE cCREs analysis |
| `CREs_splicing_genes_encode/` | ENCODE cCREs for splicing genes |
| `CREs_splicing_genes_paper/` | Table 16 CREs for splicing genes |
| `CREs_paper_specific/` | Cell-type specific Table 16 CREs |
| `CREs_paper_exclusive/` | Exclusive Table 16 CREs analysis |