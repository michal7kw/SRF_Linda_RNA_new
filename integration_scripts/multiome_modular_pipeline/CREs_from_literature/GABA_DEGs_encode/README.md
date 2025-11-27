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

Since Emx1-Ctrl sample has quality issues, we analyze 3 comparisons for each DEG set:

| Comparison | Description |
|------------|-------------|
| Nestin-Ctrl vs Nestin-Mut | Within-genotype mutation effect (Nestin) |
| Nestin-Ctrl vs Emx1-Mut | Cross-genotype mutation effect |
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

**Output**:
- `encode_cCREs_up_DEGs.tsv` - ENCODE cCREs linked to up-regulated DEGs
- `encode_cCREs_up_DEGs.bed` - BED format for deepTools
- `encode_cCREs_down_DEGs.tsv` - ENCODE cCREs linked to down-regulated DEGs
- `encode_cCREs_down_DEGs.bed` - BED format for deepTools
- `SUMMARY_GABA_DEGs_encode_cCREs.txt` - Summary statistics

**Runtime**: ~2-5 minutes

### Step 2: Compute Signal Matrices (`2_compute_signal_matrices.sh`)

**Input**: BED files from Step 1 + GABA BigWig files

**Output**:
- `matrix_up_DEGs.gz/.tab` - Matrix for up-DEG ENCODE cCREs
- `matrix_down_DEGs.gz/.tab` - Matrix for down-DEG ENCODE cCREs
- `heatmap_*.png` - Heatmaps for each DEG set
- `metaprofile_*.png` - Metaprofiles for each DEG set

**Method**: Uses deepTools (computeMatrix + plotHeatmap/plotProfile)

**Runtime**: ~30-60 minutes

### Step 3: Custom Comparison Visualizations (`3_visualize_deg_comparisons.py`)

**Input**: deepTools matrix files from Step 2

**Output**:
- `custom_comparisons/overview_all_conditions_*.png` - Combined overview
- `custom_comparisons/profiles/metaprofile_*_up_DEGs.png` - Up-DEG comparisons
- `custom_comparisons/profiles/metaprofile_*_down_DEGs.png` - Down-DEG comparisons
- `custom_comparisons/profiles/metaprofile_up_vs_down_*.png` - Up vs Down comparison
- `custom_comparisons/comparison_statistics.txt` - Statistical summaries

**Features**:
- Publication-quality metaprofiles (DPI 300)
- Difference plots showing signal changes
- Statistical tests (paired t-test)
- Optional individual CRE plots

**Runtime**: ~15-30 minutes

## Usage

### Quick Analysis (Recommended)
```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/GABA_DEGs_encode

# Run complete pipeline with SLURM job dependencies
./0_RUN_GABA_DEG_ENCODE_ANALYSIS.sh

# Or with options
./0_RUN_GABA_DEG_ENCODE_ANALYSIS.sh --skip-individual  # Fast mode
./0_RUN_GABA_DEG_ENCODE_ANALYSIS.sh --parallel 8       # Parallel individual plots
./0_RUN_GABA_DEG_ENCODE_ANALYSIS.sh --dry-run          # Show what would run
```

**Total runtime: ~60-120 minutes**

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

**GABA DEGs:**
- Up-regulated: `DEGs_cell_type_L1FC_0_25/biomarkers/sig_deg_lists/GABA/Cluster_GABA_vs_Rest_up_significant.csv`
- Down-regulated: `DEGs_cell_type_L1FC_0_25/biomarkers/sig_deg_lists/GABA/Cluster_GABA_vs_Rest_down_significant.csv`

**ENCODE cCREs:**
- Path: `../data/mm10-cCREs.bed`
- Content: Mouse (mm10) candidate cis-regulatory elements
- Types: dELS, pELS, PLS, CTCF-only, DNase-H3K4me3

**CRE-Gene Links:**
- Path: `../data/table_16.txt`
- Filters: FDR < 0.05, |PCC| > 0.2
- Cell types: GABA/hippocampal only (CA1, CA2, CA3, DG, LAMP5, VIP, SST, PV, etc.)

**ATAC-seq Signal:**
- Path: `../../signac_results_L1/bigwig_tracks_L1/by_celltype/`
- Files: `GABA_Nestin-Ctrl.bw`, `GABA_Nestin-Mut.bw`, `GABA_Emx1-Mut.bw`
- Note: Emx1-Ctrl excluded due to quality issues

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
GABA_DEGs_encode/
├── 0_RUN_GABA_DEG_ENCODE_ANALYSIS.sh  # Master script
├── 1_extract_encode_cCREs_for_DEGs.py # Step 1: Extract ENCODE cCREs
├── 2_compute_signal_matrices.sh       # Step 2: Compute matrices
├── 3_visualize_deg_comparisons.py     # Step 3: Visualizations
├── README.md                          # This file
├── logs/                              # SLURM logs
│   ├── 1_extract_encode_cCREs_for_DEGs.log
│   ├── 2_compute_signal_matrices.log
│   └── 3_visualize_deg_comparisons.log
└── output/
    ├── encode_cCREs_up_DEGs.tsv       # Up-DEG links
    ├── encode_cCREs_up_DEGs.bed       # Up-DEG BED
    ├── encode_cCREs_down_DEGs.tsv     # Down-DEG links
    ├── encode_cCREs_down_DEGs.bed     # Down-DEG BED
    ├── SUMMARY_GABA_DEGs_encode_cCREs.txt
    ├── heatmaps_deeptools/
    │   ├── matrix_up_DEGs.gz/.tab     # Up-DEG matrix
    │   ├── matrix_down_DEGs.gz/.tab   # Down-DEG matrix
    │   ├── heatmap_*.png              # Heatmaps
    │   └── metaprofile_*.png          # Metaprofiles
    └── custom_comparisons/
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
            └── metaprofile_up_vs_down_Emx1_Mut.png
```

## Biological Interpretation

### CRE Type Insights

**dELS (Distal Enhancers):**
- Most abundant CRE type
- Cell-type-specific regulation
- Key targets for SRF mutation effects

**pELS (Proximal Enhancers):**
- Close to promoters
- Important for basal expression
- Direct regulation of nearby genes

**PLS (Promoter-Like Signatures):**
- Mark active transcription start sites
- Direct correlation with gene expression

### Expected Patterns

- **Up-DEGs**: Higher ATAC signal at associated enhancers in GABA neurons
- **Down-DEGs**: Lower ATAC signal at associated enhancers in GABA neurons
- **Mutation Effects**:
  - Ctrl → Mut changes indicate direct SRF mutation effects
  - Nestin vs Emx1 differences indicate genotype-specific responses

## Related Pipelines

| Pipeline | Description |
|----------|-------------|
| `GABA_DEGs/` | Table 16 CREs directly (not ENCODE filtered) |
| `encode_cCREs/` | ENCODE cCREs for splicing genes |
| `splicing_genes/` | Original splicing gene analysis |