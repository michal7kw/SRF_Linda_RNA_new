# Splicing Gene CRE Analysis Pipeline - WITH DIRECTIONALITY

## PCC Directionality

**This pipeline separates CREs by correlation direction**:

| CRE Type | PCC Sign | Biological Meaning | Loss of Accessibility Implies |
|----------|----------|-------------------|------------------------------|
| **Enhancer-like** | PCC > 0.2 | Accessibility ↑ = Expression ↑ | Reduced gene expression |
| **Silencer-like** | PCC < -0.2 | Accessibility ↑ = Expression ↓ | Increased gene expression |

## Question

Loss of chromatin accessibility at splicing gene regulatory elements specifically in **Nestin-MUT GABA neurons**? 

compared to:
1. **Nestin-Ctrl**
2. **Emx1-Mut**

```
┌─────────────────────────────────────────────────────────────────────────────────────┐
│  Step 1: Extract CREs with PCC Directionality                                       │
│  ├── Load ~1,141 splicing genes (Reactome/GO pathways)                              │
│  ├── Query Table 16 (3.2M CRE-gene correlations)                                    │
│  ├── SEPARATE by PCC sign:                                                          │
│  │   ├── Enhancer-like: PCC > 0.2 (positive correlation)                            │
│  │   └── Silencer-like: PCC < -0.2 (negative correlation)                           │
│  ├── Filter for GABA/hippocampal cell types AND all cell types                      │
│  └── Output: enhancer_CREs_GABA.bed, silencer_CREs_GABA.bed, *_all.bed              │
├─────────────────────────────────────────────────────────────────────────────────────┤
│  Step 2a: Compute Signal Matrices (GABA-specific)                                   │
│  ├── Load GABA BigWig files (Nestin-Ctrl, Nestin-Mut, Emx1-Mut)                     │
│  ├── Compute deepTools matrices for enhancers and silencers                         │
│  ├── Create heatmaps and metaprofiles                                               │
│  └── Generate pairwise comparisons                                                  │
├─────────────────────────────────────────────────────────────────────────────────────┤
│  Step 2b: Compute Signal Matrices (All Cell Types)                                  │
│  ├── Same process but using all cell type CREs                                      │
│  └── Output: matrix_{enhancer,silencer}_all.gz                                      │
├─────────────────────────────────────────────────────────────────────────────────────┤
│  Step 3: Visualize and Identify Nestin-Specific Loss                                │
│  ├── Create publication-quality comparison plots                                    │
│  ├── Identify CREs with Nestin-specific loss pattern                                │
│  ├── Generate three-way comparison scatter plots                                    │
│  ├── Create individual CRE plots for top hits                                       │
│  └── Output: nestin_specific_loss_enhancer.tsv, nestin_specific_loss_silencer.tsv   │
├─────────────────────────────────────────────────────────────────────────────────────┤
│  Step 4: Differential CRE Analysis                                                  │
│  ├── Identify differential CREs using min_signal and min_fc thresholds             │
│  ├── Three comparisons: Ctrl vs Mut, Ctrl vs Emx1-Mut, Nestin-Mut vs Emx1-Mut      │
│  ├── Separate analysis for enhancers and silencers                                  │
│  ├── Runs for BOTH all cell types AND GABA-specific CREs                           │
│  ├── Create volcano plots, MA plots, and metaprofiles                              │
│  └── Output: differential_CREs_minSig{X}_minFC{Y}/ organized directory             │
└─────────────────────────────────────────────────────────────────────────────────────┘
```

## Quick Start

```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_directional

# Run complete pipeline
./0_RUN_SPLICING_DIRECTIONAL_PIPELINE.sh

# Or with options
./0_RUN_SPLICING_DIRECTIONAL_PIPELINE.sh --skip-individual  # Faster (no individual plots)
./0_RUN_SPLICING_DIRECTIONAL_PIPELINE.sh --dry-run          # Show what would run
./0_RUN_SPLICING_DIRECTIONAL_PIPELINE.sh --step 3           # Run only step 3
```

## Data Sources

### Primary CRE Data Files

| File | Used? | Purpose |
|------|-------|---------|
| **table_16.txt** | **YES** | CRE-gene correlations with PCC values |
| **mm10-cCREs.bed** | Optional | ENCODE CRE type annotations |

**table_16.txt** (Literature CRE-Gene Correlations):
- Path: `../data/table_16.txt`
- Content: 3,256,804 published CRE-gene correlations from ENCODE consortium
- Columns: Coordinate1, cCRE1, Gene, SubType, **PCC**, FDR, etc.
- **Critical**: PCC sign determines enhancer vs silencer classification

### BigWig Files (ATAC-seq Signal)

| File | Size | Status |
|------|------|--------|
| GABA_Nestin-Ctrl.bw | 12.9 MB | **Reference Control** |
| GABA_Nestin-Mut.bw | 25.4 MB | Primary comparison |
| GABA_Emx1-Mut.bw | 28.0 MB | Genotype comparison |
| GABA_Emx1-Ctrl.bw | 2.7 MB | **EXCLUDED (failed sample)** |

### Splicing Genes

- **Source**: `/beegfs/.../splicing_genes/extracted_genes_final.csv`
- **Count**: ~1,141 unique genes
- **Pathways**: mRNA splicing, spliceosome assembly, RNA processing (Reactome, GO)

## CRE Filtering Criteria

| Filter | Value | Applied? |
|--------|-------|----------|
| **FDR threshold** | < 0.05 | YES |
| **PCC threshold (enhancers)** | > 0.2 | YES (positive only) |
| **PCC threshold (silencers)** | < -0.2 | YES (negative only) |
| **Cell type filter** | GABA/Hippocampal | YES |

**Cell Type Keywords** (SubType filter):
- Hippocampal: CA1, CA2, CA3, CA4, DG, DGNBL, GRC
- GABAergic: GABA, PV, SST, VIP, LAMP5, LAMP, PVGA, SSTGA, VIPGA, LAMGA, INH, CGE, MGE, LGE

## Nestin-Specific Loss Detection

Parameters for identifying CREs with Nestin-specific loss:

| Parameter | Default | Description |
|-----------|---------|-------------|
| MIN_CTRL_SIGNAL | 0.5 | Minimum signal in Nestin-Ctrl (must be detectable) |
| MIN_FOLD_CHANGE | 1.5 | Minimum fold change (Ctrl/Mut) for "loss" |
| MAX_EMX1_LOSS_RATIO | 0.7 | Emx1-Mut must retain ≥70% of Ctrl signal |

**Detection Logic**:
```
Nestin-specific loss = (
    Nestin-Ctrl ≥ 0.5  AND           # Baseline detectable
    Nestin-Ctrl / Nestin-Mut ≥ 1.5 AND  # Loss in Nestin-Mut
    Emx1-Mut / Nestin-Ctrl ≥ 0.7     # Preserved in Emx1-Mut
)
```

## Output Structure

```
CREs_splicing_directional/
├── 0_RUN_SPLICING_DIRECTIONAL_PIPELINE.sh  # Master script
├── 1_extract_splicing_CREs_directional.py   # Step 1
├── 1_extract_splicing_CREs_directional.sh   # SLURM wrapper
├── 2_compute_signal_matrices.sh             # Step 2
├── 3_visualize_directional_comparisons.py   # Step 3
├── 3_visualize_directional_comparisons.sh   # SLURM wrapper
├── 4_differential_CRE_analysis.py           # Step 4: Differential CRE analysis
├── 4_differential_CRE_analysis.sh           # SLURM wrapper for Step 4
├── README.md                                # This file
├── logs/                                    # SLURM logs
└── output/
    ├── enhancer_CREs_GABA.bed               # Unique enhancer CREs for deepTools (PCC > 0.2)
    ├── enhancer_CREs_GABA.tsv               # Full enhancer data with all columns
    ├── enhancer_CREs_GABA_full.bed          # Full BED with gene info
    ├── enhancer_CREs_GABA_gene_links.tsv    # CRE-gene linkage details
    ├── enhancer_CREs_all.bed                # All cell types (not just GABA)
    ├── enhancer_CREs_all_celltypes.tsv      # All cell type data
    ├── enhancer_CREs_all_full.bed           # All cell types full BED
    ├── enhancer_CREs_all_gene_links.tsv     # All cell type gene links
    ├── silencer_CREs_GABA.bed               # Unique silencer CREs for deepTools (PCC < -0.2)
    ├── silencer_CREs_GABA.tsv               # Full silencer data with all columns
    ├── silencer_CREs_GABA_full.bed          # Full BED with gene info
    ├── silencer_CREs_GABA_gene_links.tsv    # CRE-gene linkage details
    ├── silencer_CREs_all.bed                # All cell types (not just GABA)
    ├── silencer_CREs_all_celltypes.tsv      # All cell type data
    ├── silencer_CREs_all_full.bed           # All cell types full BED
    ├── silencer_CREs_all_gene_links.tsv     # All cell type gene links
    ├── SUMMARY_splicing_CREs_directional.txt
    │
    ├── matrices/                            # deepTools matrices
    │   ├── matrix_enhancer_GABA.gz          # Main enhancer matrix (all 3 conditions)
    │   ├── matrix_enhancer_GABA.tab         # Tab-separated version
    │   ├── matrix_enhancer_nestin_comparison.gz   # Nestin-Ctrl vs Nestin-Mut
    │   ├── matrix_enhancer_cross_comparison.gz    # Nestin-Ctrl vs Emx1-Mut
    │   ├── matrix_enhancer_mutant_comparison.gz   # Nestin-Mut vs Emx1-Mut
    │   ├── matrix_silencer_GABA.gz          # Main silencer matrix (all 3 conditions)
    │   ├── matrix_silencer_GABA.tab         # Tab-separated version
    │   ├── matrix_silencer_nestin_comparison.gz   # Pairwise comparisons
    │   ├── matrix_silencer_cross_comparison.gz
    │   ├── matrix_silencer_mutant_comparison.gz
    │   └── matrix_combined_enhancer_silencer.gz   # Combined enhancer + silencer
    │
    ├── heatmaps_deeptools/                  # deepTools visualizations
    │   ├── heatmap_enhancer_all_conditions.png
    │   ├── heatmap_silencer_all_conditions.png
    │   ├── heatmap_combined_enhancer_vs_silencer.png
    │   ├── metaprofile_enhancer_all_conditions.png
    │   ├── metaprofile_enhancer_nestin_ctrl_vs_mut.png
    │   ├── metaprofile_enhancer_nestin_ctrl_vs_emx1_mut.png
    │   ├── metaprofile_enhancer_nestin_mut_vs_emx1_mut.png
    │   ├── metaprofile_silencer_all_conditions.png
    │   ├── metaprofile_silencer_nestin_ctrl_vs_mut.png
    │   ├── metaprofile_silencer_nestin_ctrl_vs_emx1_mut.png
    │   ├── metaprofile_silencer_nestin_mut_vs_emx1_mut.png
    │   └── metaprofile_combined_enhancer_vs_silencer.png
    │
    ├── visualizations/                      # Custom comparison plots
    │   ├── metaprofile_enhancer_with_difference.png
    │   ├── metaprofile_silencer_with_difference.png
    │   ├── scatter_enhancer_three_way_comparison.png
    │   ├── scatter_silencer_three_way_comparison.png
    │   └── summary_barplot_all_conditions.png
    │
    ├── nestin_specific_loss/                # KEY RESULTS (Step 3)
    │   ├── nestin_specific_loss_enhancer.tsv
    │   ├── nestin_specific_loss_silencer.tsv
    │   ├── SUMMARY_nestin_specific_loss.txt
    │   ├── individual_enhancer_minSig2.0_minFC2.0/
    │   │   └── cre_001_GeneName_fc2.5.png
    │   └── individual_silencer_minSig2.0_minFC2.0/
    │       └── cre_001_GeneName_fc2.5.png
    │
    └── differential_CREs_minSig0.05_minFC1.5/  # DIFFERENTIAL ANALYSIS (Step 4)
        ├── SUMMARY_differential_CRE_analysis.txt
        │
        ├── results/                          # TSV files with statistics
        │   ├── differential_CREs_enhancer_nestin_ctrl_vs_mut.tsv
        │   ├── differential_CREs_enhancer_nestin_ctrl_vs_mut_significant.tsv
        │   ├── differential_CREs_enhancer_nestin_ctrl_vs_emx1_mut.tsv
        │   ├── differential_CREs_enhancer_nestin_mut_vs_emx1_mut.tsv
        │   ├── differential_CREs_silencer_*.tsv
        │   └── ...
        │
        ├── plots/                            # Visualization plots
        │   ├── summary_enhancer_differential.png  # Stacked bar overview
        │   ├── summary_silencer_differential.png
        │   ├── enhancer/                     # Per CRE type
        │   │   ├── nestin_ctrl_vs_mut/       # Per comparison
        │   │   │   ├── volcano_enhancer_nestin_ctrl_vs_mut.png
        │   │   │   ├── ma_plot_enhancer_nestin_ctrl_vs_mut.png
        │   │   │   ├── metaprofile_up_enhancer_nestin_ctrl_vs_mut.png
        │   │   │   └── metaprofile_down_enhancer_nestin_ctrl_vs_mut.png
        │   │   ├── nestin_ctrl_vs_emx1_mut/
        │   │   └── nestin_mut_vs_emx1_mut/
        │   └── silencer/
        │       └── ...
        │
        ├── bed_files/                        # BED files for IGV/deepTools
        │   ├── enhancer/
        │   │   ├── dCREs_enhancer_nestin_ctrl_vs_mut_significant.bed
        │   │   ├── dCREs_enhancer_nestin_ctrl_vs_mut_up.bed
        │   │   ├── dCREs_enhancer_nestin_ctrl_vs_mut_down.bed
        │   │   └── ...
        │   └── silencer/
        │       └── ...
        │
        └── individual_plots/                 # Individual CRE signal profiles
            ├── enhancer/
            │   ├── nestin_ctrl_vs_mut/
            │   │   └── dCRE_001_up_GeneName_fc3.5.png
            │   └── ...
            └── silencer/
                └── ...
```

## Key Output Files

### 1. Nestin-Specific Loss Tables

**`output/nestin_specific_loss/nestin_specific_loss_enhancer.tsv`**

Contains CREs where:
- Nestin-Ctrl has HIGH signal
- Nestin-Mut has LOW signal (≥1.5x reduction)
- Emx1-Mut retains ≥70% of Ctrl signal

Columns:
- `chrom`, `start`, `end`: CRE coordinates
- `genes`: Linked splicing gene(s)
- `Nestin-Ctrl`, `Nestin-Mut`, `Emx1-Mut`: Mean signals
- `fc_nestin`: Fold change (Ctrl/Mut)
- `emx1_retention`: Emx1-Mut / Nestin-Ctrl ratio

### 2. Visualization Plots

**Key plots to examine:**

| Plot | Location | What it shows |
|------|----------|---------------|
| `metaprofile_enhancer_with_difference.png` | visualizations/ | Mean signal with Mut-Ctrl difference panel |
| `metaprofile_silencer_with_difference.png` | visualizations/ | Mean signal with Mut-Ctrl difference panel |
| `scatter_enhancer_three_way_comparison.png` | visualizations/ | All three conditions as scatter plots |
| `scatter_silencer_three_way_comparison.png` | visualizations/ | All three conditions as scatter plots |
| `summary_barplot_all_conditions.png` | visualizations/ | Mean signals with significance (t-test) |
| `heatmap_combined_enhancer_vs_silencer.png` | heatmaps_deeptools/ | Side-by-side enhancer vs silencer comparison |
| `metaprofile_*_nestin_ctrl_vs_mut.png` | heatmaps_deeptools/ | Within-genotype mutation effect |
| `metaprofile_*_nestin_ctrl_vs_emx1_mut.png` | heatmaps_deeptools/ | Cross-genotype comparison |
| `metaprofile_*_nestin_mut_vs_emx1_mut.png` | heatmaps_deeptools/ | Mutant genotype comparison |

## Running Individual Steps

```bash
# Step 1 only: Extract CREs
sbatch 1_extract_splicing_CREs_directional.sh

# Step 2 only: Compute matrices (requires Step 1 output)
sbatch 2_compute_signal_matrices.sh

# Step 3 only: Visualize (requires Step 2 output)
sbatch 3_visualize_directional_comparisons.sh

# Step 3 with options (via environment variables)
SKIP_INDIVIDUAL=1 sbatch 3_visualize_directional_comparisons.sh  # Fast mode (metaprofiles only)
MAX_INDIVIDUAL=100 sbatch 3_visualize_directional_comparisons.sh # More individual plots

# Step 3 with signal/fold-change thresholds for individual plots
# Default is MIN_SIGNAL=2.0 MIN_FC=2.0 (stricter filtering)
MIN_SIGNAL=3.0 MIN_FC=3.0 sbatch 3_visualize_directional_comparisons.sh  # Even stricter
MIN_SIGNAL=1.0 MIN_FC=1.5 sbatch 3_visualize_directional_comparisons.sh  # More permissive

# Combined options
MAX_INDIVIDUAL=100 sbatch 3_visualize_directional_comparisons.sh  # More plots with default thresholds

# Step 4 only: Differential CRE analysis (requires Step 2 output)
sbatch 4_differential_CRE_analysis.sh

# Step 4 with options (via environment variables)
DIFF_MIN_SIGNAL=1.5 DIFF_MIN_FC=2.0 sbatch 4_differential_CRE_analysis.sh  # More permissive
DIFF_MIN_SIGNAL=3.0 DIFF_MIN_FC=4.0 sbatch 4_differential_CRE_analysis.sh  # Stricter

# Step 4 with GABA-specific only (via USE_GABA=1)
USE_GABA=1 sbatch 4_differential_CRE_analysis.sh

# Skip differential analysis in main pipeline
SKIP_DIFFERENTIAL=1 sbatch 0_RUN_SPLICING_DIRECTIONAL_PIPELINE.sh
```

### Individual Plot Filtering Parameters

Individual CRE plots are filtered based on signal and fold-change thresholds:

| Parameter | Default | Description |
|-----------|---------|-------------|
| MIN_SIGNAL | 2.0 | Minimum max signal required (at least one condition must exceed) |
| MIN_FC | 2.0 | Minimum fold change required (Ctrl/Mut >= min_fc OR Mut/Ctrl >= min_fc) |
| MAX_INDIVIDUAL | 50 | Maximum number of individual plots to create per CRE type |
| INDIVIDUAL_DPI | 150 | DPI for individual plots (metaprofiles always 300) |

**Output directories include threshold suffix**: `individual_enhancer_minSig2.0_minFC2.0/`

### Differential CRE Analysis Parameters (Step 4)

**IMPORTANT**: Default thresholds are set for normalized BigWig signal values (typically 0.01-1.0 range).

| Parameter | Default | Description |
|-----------|---------|-------------|
| DIFF_MIN_SIGNAL | 0.05 | Minimum max signal required for differential detection |
| DIFF_MIN_FC | 1.5 | Minimum fold change required for significance |
| SKIP_DIFFERENTIAL | 0 | Set to 1 to skip differential analysis in main pipeline |
| DIFF_MAX_INDIVIDUAL | 20 | Maximum individual plots per comparison |
| USE_GABA | 0 | Set to 1 to analyze only GABA-specific CREs |
| DIAGNOSTICS_ONLY | 0 | Set to 1 to only show signal diagnostics |

**Differential CRE Detection Logic**:
```
Significant UP = (max_signal >= min_signal) AND (fc >= min_fc)
Significant DOWN = (max_signal >= min_signal) AND (fc <= 1/min_fc)
```

**Comparisons performed**:
- **nestin_ctrl_vs_mut**: Within-genotype mutation effect
- **nestin_ctrl_vs_emx1_mut**: Cross-genotype comparison
- **nestin_mut_vs_emx1_mut**: Mutant genotype comparison

**Output directories include threshold suffix**: `differential_CREs_minSig0.05_minFC1.5/`

## Comparison with Previous Pipelines

| Feature | Previous Pipelines | This Pipeline |
|---------|-------------------|---------------|
| **PCC filtering** | \|PCC\| > 0.2 (mixed) | PCC > 0.2 OR PCC < -0.2 (separated) |
| **CRE types** | Combined | Enhancers and Silencers separately |
| **Nestin-specific** | Not identified | Explicitly detected |
| **Differential CREs** | Not performed | Full differential analysis (Step 4) |
| **Cell type scope** | GABA only | Both GABA-specific AND all cell types |
| **Biological interpretation** | Ambiguous | Clear directional model |
| **Output organization** | Flat | Hierarchical with clear naming |