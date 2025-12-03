# Differential Splicing CRE Analysis (All Cell Types)

This pipeline identifies **differential Cis-Regulatory Elements (CREs)** associated with splicing genes, without restricting the search to specific cell types (e.g., GABAergic).

## Goal
To find regulatory elements that show significant changes in chromatin accessibility between conditions (Nestin-Ctrl, Nestin-Mut, Emx1-Mut), which might explain splicing defects.

## Key Features
1.  **No Cell Type Filter**: Uses all CREs linked to splicing genes from the source paper (Table 16).
2.  **Directional**: Separates **Enhancers** (PCC > 0) and **Silencers** (PCC < 0).
3.  **Differential Analysis**: Performs pairwise comparisons with strict thresholds (`min_fc=3.0`).

## Pipeline Steps

### 1. Extraction (`1_extract_differential_candidates.py`)
*   Extracts CREs from Table 16.
*   Filters by FDR < 0.05.
*   Separates Enhancers (PCC > 0.2) and Silencers (PCC < -0.2).
*   **Output**: `output/enhancer_CREs_all_celltypes.bed`, `output/silencer_CREs_all_celltypes.bed`

### 2. Matrix Computation (`2_compute_signal_matrices.sh`)
*   Uses `deepTools` to compute ATAC-seq signal matrices for the extracted CREs.
*   Input: BigWig files (Nestin-Ctrl, Nestin-Mut, Emx1-Mut).
*   **Output**: `output/matrices/matrix_*.gz`

### 3. Analysis (`3_analyze_differential_cres.py`)
*   Performs pairwise comparisons:
    *   Nestin-Ctrl vs Nestin-Mut
    *   Nestin-Ctrl vs Emx1-Mut
    *   Nestin-Mut vs Emx1-Mut
*   Identifies differential CREs (FC >= 3.0, Min Signal >= 2.0).
*   **Output**:
    *   `output/differential_lists/`: TSV files of differential CREs.
    *   `output/visualizations/`: Volcano plots, metaprofiles, and individual CRE plots.

## How to Run

```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_differential_all_celltypes

# Run the entire pipeline
sbatch 0_RUN_ALL.sh
```
