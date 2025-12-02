# GABA DEG CRE Analysis Pipeline

This pipeline visualizes ATAC-seq signal at CREs (cis-regulatory elements) linked to GABA-specific differentially expressed genes (DEGs).

## Overview

The pipeline creates publication-quality metaprofile plots showing ATAC signal around CREs linked to:
- **Up-regulated DEGs**: Genes enriched in GABA neurons (GABA markers)
- **Down-regulated DEGs**: Genes depleted in GABA neurons

## Quick Start

Run the complete pipeline with a single command:

```bash
sbatch 0_RUN_GABA_DEG_ANALYSIS.sh
```

## Pipeline Steps

| Script | Description | Runtime |
|--------|-------------|---------|
| `1_extract_gaba_deg_CREs.py` | Extract CREs linked to DEGs from Table 16 | ~5 min |
| `2_compute_signal_matrices.sh` | Compute deepTools matrices | ~30 min |
| `3_visualize_deg_comparisons.py` | Create visualization plots | ~10 min |

## Data Sources

### Primary CRE Data Files

| File | Used? | Purpose |
|------|-------|---------|
| **mm10-cCREs.bed** | **NO** | Not used in this pipeline |
| **table_16.txt** | **YES** | Primary CRE source - literature CRE-gene correlations |

**mm10-cCREs.bed** (ENCODE cCREs):
- **NOT USED** in this pipeline
- For ENCODE cCRE type annotations, see `CREs_GABA_DEGs_encode/`

**table_16.txt** (Literature CRE-Gene Correlations):
- Path: `../data/table_16.txt`
- Content: 3,256,804 published CRE-gene correlations from ENCODE consortium
- Columns: Coordinate1 (CRE coordinates), cCRE1, Gene, SubType, PCC, FDR, etc.
- Usage: **Primary CRE source** - CRE coordinates extracted directly from Coordinate1 column
- Statistical filters: FDR < 0.05, |PCC| > 0.2
- Cell type filter: GABA/Hippocampal SubTypes only

### CRE-Gene Linkage Method (Literature-Based)

This pipeline uses **ONLY table_16.txt**:
1. Filter Table 16 for GABA DEGs (up/down) with significant correlations
2. Apply GABA/hippocampal cell type filter on SubType column
3. Extract CRE coordinates directly from Coordinate1 column (chr_start_end format)
4. Result: CREs with **literature-validated** gene links (no ENCODE type annotations)

### Other Input Files

**GABA DEG Lists**:
- `Cluster_GABA_vs_Rest_up_significant.csv` (~5,851 genes)
- `Cluster_GABA_vs_Rest_down_significant.csv` (~797 genes)

**BigWig Files**: GABA ATAC signal from Signac pipeline
- `GABA_Nestin-Ctrl.bw`
- `GABA_Nestin-Mut.bw`
- `GABA_Emx1-Mut.bw`

## Genes Analyzed

**Gene Set**: GABA-specific differentially expressed genes (DEGs)
- **Up-regulated DEGs**: 5,851 genes enriched in GABA neurons (GABA markers)
  - 4,089 genes (69.9%) have CRE links → 4,024 unique CREs, 140,087 total links
- **Down-regulated DEGs**: 797 genes depleted in GABA neurons
  - 605 genes (75.9%) have CRE links → 620 unique CREs, 21,317 total links
- **Source**: DEG analysis comparing GABA vs Rest cell types (`DEGs_cell_type_L1FC_0_25`)
- **Analysis**: Separate processing for up/down DEGs

## CRE Filtering Criteria

| Filter | Value | Applied? |
|--------|-------|----------|
| **FDR threshold** | < 0.05 | YES - Benjamini-Hochberg correction |
| **PCC threshold** | > 0.2 | YES - Pearson correlation coefficient |
| **Distance window** | N/A | NO - uses Table 16 correlations |
| **Cell type filter** | GABA/Hippocampal | YES - SubType keywords |
| **ENCODE type filter** | N/A | NO - no ENCODE cCREs used |
| **ENCODE intersect** | N/A | NO - uses Table 16 coordinates directly |

**Cell Type Keywords** (SubType filter):
- Hippocampal: CA1, CA2, CA3, DG, DGNBL, GRC
- GABAergic: LAMP5, LAMP, VIP, SST, PV, PVGA, SSTGA, VIPGA, LAMGA, GABA, INH, interneuron

**Filtering Pipeline**:
1. Filter Table 16 for GABA DEGs (up or down)
2. Apply statistical filters (FDR < 0.05, |PCC| > 0.2)
3. Apply GABA/hippocampal cell type filter on SubType
4. Extract CRE coordinates from Coordinate1 column
5. Result: Literature-validated CREs (no ENCODE type annotations)

## Comparisons

For each DEG set, three comparisons are created:

1. **Nestin-Ctrl vs Nestin-Mut**: Within-genotype mutation effect
2. **Nestin-Ctrl vs Emx1-Mut**: Cross-genotype comparison
3. **Nestin-Mut vs Emx1-Mut**: Mutant-to-mutant genotype effect

> **Note**: Emx1-Ctrl is excluded (failed sample)

## Output Files

```
output/
  GABA_DEGs_up_CREs.tsv        # CRE links for up-regulated DEGs (140,087 links)
  GABA_DEGs_up_CREs.bed        # BED format (4,024 unique CREs)
  GABA_DEGs_down_CREs.tsv      # CRE links for down-regulated DEGs (21,317 links)
  GABA_DEGs_down_CREs.bed      # BED format (620 unique CREs)

  # Profile directories include filter parameters in name
  profiles_up_minSig2.0_minFC2.0/    # Up-regulated DEG visualizations
    metaprofile_nestin_ctrl_vs_mut.png
    metaprofile_nestin_ctrl_vs_emx1_mut.png
    metaprofile_nestin_mut_vs_emx1_mut.png
    individual_up_*.png              # Individual CRE plots (if enabled)

  profiles_down_minSig2.0_minFC2.0/  # Down-regulated DEG visualizations
    metaprofile_nestin_ctrl_vs_mut.png
    metaprofile_nestin_ctrl_vs_emx1_mut.png
    metaprofile_nestin_mut_vs_emx1_mut.png
    individual_down_*.png            # Individual CRE plots (if enabled)

  matrix_*.gz / matrix_*.tab   # deepTools signal matrices
  computeMatrix_*.log          # deepTools log files
  comparison_statistics.txt    # Statistical summaries
  SUMMARY_GABA_DEGs_CREs.txt   # Analysis report
```

## Performance Options

### Fast Mode (Default)
Creates metaprofiles only (skips individual CRE plots):
```bash
sbatch 0_RUN_GABA_DEG_ANALYSIS.sh
```

### Full Mode
Creates both metaprofiles and individual CRE plots:
```bash
SKIP_INDIVIDUAL=0 sbatch 0_RUN_GABA_DEG_ANALYSIS.sh
```

### Full Mode with Parallelization
Faster individual plot generation:
```bash
SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 0_RUN_GABA_DEG_ANALYSIS.sh
```

### Custom Filtering
Adjust individual plot filtering thresholds:
```bash
SKIP_INDIVIDUAL=0 MIN_SIGNAL=2.0 MIN_FC=2.0 sbatch 0_RUN_GABA_DEG_ANALYSIS.sh
```

**Filter parameters:**
- `MIN_SIGNAL`: Minimum max ATAC signal to include CRE (script default: 1.0)
- `MIN_FC`: Minimum fold change to include CRE (script default: 1.5)

> **Note**: The current output uses `minSig2.0_minFC2.0` because stricter filters were applied during generation.

## Running Individual Steps

```bash
# Step 1: Extract CREs
python 1_extract_gaba_deg_CREs.py

# Step 2: Compute matrices
sbatch 2_compute_signal_matrices.sh

# Step 3: Create visualizations
# Fast mode (metaprofiles only):
python 3_visualize_deg_comparisons.py --skip-individual

# Full mode with default filters:
python 3_visualize_deg_comparisons.py --parallel 8

# Full mode with custom filters:
python 3_visualize_deg_comparisons.py --parallel 8 --min-signal 2.0 --min-fc 2.0
```

## Visualization Style

The plots follow the same style as the splicing genes analysis:
- **Top panel**: Mean ATAC signal with SEM error bands
- **Bottom panel**: Difference plot (Condition2 - Condition1)
- **Colors**: Blue (Ctrl), Red (Nestin-Mut), Orange (Emx1-Mut)

## Related Analyses

This pipeline is modeled after:
- `../CREs_splicing_genes_paper/` - CRE analysis for splicing genes
