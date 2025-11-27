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

## Input Data

1. **GABA DEG Lists**:
   - `Cluster_GABA_vs_Rest_up_significant.csv` (~5,851 genes)
   - `Cluster_GABA_vs_Rest_down_significant.csv` (~797 genes)

2. **CRE-Gene Correlations**: Table 16 (literature data)
   - Filters: FDR < 0.05, |PCC| > 0.2
   - Cell types: GABA/Hippocampal subtypes

3. **BigWig Files**: GABA ATAC signal from Signac pipeline
   - `GABA_Nestin-Ctrl.bw`
   - `GABA_Nestin-Mut.bw`
   - `GABA_Emx1-Mut.bw`

## Comparisons

For each DEG set, three comparisons are created:

1. **Nestin-Ctrl vs Nestin-Mut**: Within-genotype mutation effect
2. **Nestin-Ctrl vs Emx1-Mut**: Cross-genotype comparison
3. **Nestin-Mut vs Emx1-Mut**: Mutant-to-mutant genotype effect

> **Note**: Emx1-Ctrl is excluded (failed sample)

## Output Files

```
output/
  GABA_DEGs_up_CREs.tsv        # CRE links for up-regulated DEGs
  GABA_DEGs_up_CREs.bed        # BED format
  GABA_DEGs_down_CREs.tsv      # CRE links for down-regulated DEGs
  GABA_DEGs_down_CREs.bed      # BED format

  profiles_up/                  # Up-regulated DEG visualizations
    metaprofile_nestin_ctrl_vs_mut.png
    metaprofile_nestin_ctrl_vs_emx1_mut.png
    metaprofile_nestin_mut_vs_emx1_mut.png

  profiles_down/                # Down-regulated DEG visualizations
    metaprofile_nestin_ctrl_vs_mut.png
    metaprofile_nestin_ctrl_vs_emx1_mut.png
    metaprofile_nestin_mut_vs_emx1_mut.png

  comparison_statistics.txt     # Statistical summaries
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

## Running Individual Steps

```bash
# Step 1: Extract CREs
python 1_extract_gaba_deg_CREs.py

# Step 2: Compute matrices
sbatch 2_compute_signal_matrices.sh

# Step 3: Create visualizations
python 3_visualize_deg_comparisons.py --skip-individual
# OR for full analysis:
python 3_visualize_deg_comparisons.py --parallel 8
```

## Visualization Style

The plots follow the same style as the splicing genes analysis:
- **Top panel**: Mean ATAC signal with SEM error bands
- **Bottom panel**: Difference plot (Condition2 - Condition1)
- **Colors**: Blue (Ctrl), Red (Nestin-Mut), Orange (Emx1-Mut)

## Dependencies

- Python: numpy, pandas, matplotlib, seaborn, scipy
- deepTools: computeMatrix
- Conda environment: `rna_seq_analysis_deep`

## Related Analyses

This pipeline is modeled after:
- `../splicing_genes/` - CRE analysis for splicing genes
