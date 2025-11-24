# Input and Output Files (Splicing Genes Workflow)

## Overview
This workflow analyzes CREs associated with splicing-related genes from Reactome pathways, using literature-based CRE-gene correlations (Table 16). The pipeline extracts CREs linked to 1,138 splicing genes and visualizes ATAC-seq signal patterns.

**Key Feature**: All CREs are from **published literature data** (Table 16), not from Signac peak-gene links.

---

## 1. Extract Splicing Gene CREs (`1_extract_splicing_gene_CREs.py`)

**Inputs:**
- `/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/splicing_genes/extracted_genes_final.csv`
  - Splicing genes from Reactome pathways (1,138 genes)
- `../data/table_16.txt`
  - Literature CRE-gene correlations (3,256,804 published links, 567 MB)

**Outputs:**
- `output/splicing_genes_analysis/splicing_genes_CREs_all_celltypes.tsv`
  - All CREs linked to splicing genes (any cell type)
  - ~832 unique CREs, multiple cell types per CRE
- `output/splicing_genes_analysis/splicing_genes_CREs_GABA.tsv`
  - CREs linked to splicing genes in GABA/hippocampal cell types
  - Filtered by SubType column in Table 16
- `output/splicing_genes_analysis/splicing_genes_CREs_GABA_specific.tsv`
  - GABA-exclusive CREs (appear ONLY in GABA/hippocampal SubTypes)
- `output/splicing_genes_analysis/SUMMARY_splicing_genes_CREs.txt`
  - Summary statistics and metadata

**Filters Applied:**
- FDR < 0.05
- |PCC| > 0.1 (Pearson correlation coefficient)
- Cell type classification based on SubType column from Table 16

---

## 2. Convert TSV to BED Format (`2_convert_splicing_CREs_to_bed.py`)

**Inputs:**
- `output/splicing_genes_analysis/splicing_genes_CREs_all_celltypes.tsv`
- `output/splicing_genes_analysis/splicing_genes_CREs_GABA.tsv`

**Outputs:**
- `output/splicing_genes_analysis/splicing_genes_CREs_all.bed` (BED6 format)
  - All CREs, all cell types (~832 CREs)
- `output/splicing_genes_analysis/splicing_genes_CREs_GABA.bed` (BED6 format)
  - GABA CREs only (~772 CREs)

**BED Format:** chr, start, end, cre_id, score (0), strand (.)

---

## 3. Create ATAC Signal Profiles (`3_create_splicing_profiles.sh`)

**Inputs:**
- `output/splicing_genes_analysis/splicing_genes_CREs_GABA.bed`
- `../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw`
  - GABA_Nestin-Ctrl.bw
  - GABA_Nestin-Mut.bw
  - GABA_Emx1-Ctrl.bw
  - GABA_Emx1-Mut.bw

**Outputs:**
- `output/splicing_genes_analysis/heatmaps_deeptools/` (directory)
  - **Matrices** (reusable for further analysis):
    - `matrix_GABA_nestin.gz` (Nestin Ctrl vs Mut)
    - `matrix_GABA_nestin.tab`
    - `matrix_GABA_emx1.gz` (Emx1 Ctrl vs Mut)
    - `matrix_GABA_emx1.tab`

  - **Metaprofiles** (always created, 300 DPI):
    - `profiles/metaprofile_nestin_ctrl_vs_mut.png`
    - `profiles/metaprofile_emx1_ctrl_vs_mut.png`

  - **Individual CRE plots** (optional, skipped by default):
    - `profiles/individual_Nestin_<Gene>_<CRE_ID>.png` (~772 plots)
    - `profiles/individual_Emx1_<Gene>_<CRE_ID>.png` (~772 plots)

  - **Documentation:**
    - `README.txt` (interpretation guide)
    - `computeMatrix_nestin.log`
    - `computeMatrix_emx1.log`

**Performance:**
- Default (metaprofiles only): 3-5 minutes
- Full mode (with individual plots): 30-50 minutes
- Full mode + parallelization: 8 minutes

**Usage:**
```bash
# Default (fast, metaprofiles only)
sbatch 3_create_splicing_profiles.sh

# Full mode with individual plots
SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 3_create_splicing_profiles.sh
```

---

## 4. Create Traditional Heatmaps (`4_create_heatmaps_splicing_genes.sh`)

**Inputs:**
- Same as Step 3

**Outputs:**
- `output/splicing_genes_analysis/heatmaps_deeptools/`
  - `heatmap_GABA_nestin.png` (traditional deepTools heatmap)
  - `heatmap_GABA_emx1.png`
  - Additional deepTools-generated plots

**Note:** Uses deepTools plotHeatmap/plotProfile commands (not custom Python)

---

## 5. Visualize BigWig Signal (`5_visualize_bigwig_signal.sh`)

**Inputs:**
- `output/splicing_genes_analysis/splicing_genes_CREs_GABA.bed`
- `output/splicing_genes_analysis/splicing_genes_CREs_all_celltypes.tsv` (for gene mapping)
- `../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw`

**Outputs:**
- `output/splicing_genes_analysis/bigwig_profiles/` (directory)
  - **Metaprofiles** (always created, 300 DPI):
    - `metaprofile_nestin_ctrl_vs_mut.png`
    - `metaprofile_emx1_ctrl_vs_mut.png`

  - **Individual CRE plots** (optional, skipped by default):
    - `individual_Nestin_<Gene>_<CRE_ID>.png` (~772 plots)
    - `individual_Emx1_<Gene>_<CRE_ID>.png` (~772 plots)

  - **Summary statistics:**
    - `summary_nestin.tsv`
    - `summary_emx1.tsv`

**Method:** Direct BigWig reading with pyBigWig (alternative to deepTools)

**Usage:**
```bash
# Default (fast, metaprofiles only)
sbatch 5_visualize_bigwig_signal.sh

# Full mode with individual plots
SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 5_visualize_bigwig_signal.sh
```

---

## 6. Create Custom Comparisons (`6_create_custom_comparisons.sh`)

**Inputs:**
- `output/splicing_genes_analysis/splicing_genes_CREs_all.bed` (ALL CREs, not just GABA)
- `../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw`
  - GABA_Nestin-Ctrl.bw
  - GABA_Nestin-Mut.bw
  - GABA_Emx1-Mut.bw (Emx1-Ctrl excluded - failed sample)

**Outputs:**
- `output/splicing_genes_analysis/custom_comparisons/` (directory)
  - **Matrices** (3 custom comparisons):
    - `matrix_nestin_ctrl_vs_mut.gz` / `.tab`
    - `matrix_nestin_ctrl_vs_emx1_mut.gz` / `.tab`
    - `matrix_nestin_mut_vs_emx1_mut.gz` / `.tab`

  - **Metaprofiles** (always created, 300 DPI):
    - `profiles/metaprofile_nestin_ctrl_vs_nestin_mut.png`
    - `profiles/metaprofile_nestin_ctrl_vs_emx1_mut.png`
    - `profiles/metaprofile_nestin_mut_vs_emx1_mut.png`

  - **Individual CRE plots** (optional, skipped by default):
    - `profiles/individual_comparison1_<Gene>_<CRE_ID>.png` (~832 plots)
    - `profiles/individual_comparison2_<Gene>_<CRE_ID>.png` (~832 plots)
    - `profiles/individual_comparison3_<Gene>_<CRE_ID>.png` (~832 plots)
    - Total: ~2,496 plots

  - **Statistics and documentation:**
    - `comparison_statistics.txt`
    - `README.txt`
    - `computeMatrix_*.log` (3 log files)

**Comparisons:**
1. **Nestin-Ctrl vs Nestin-Mut**: Within-genotype mutation effect
2. **Nestin-Ctrl vs Emx1-Mut**: Cross-genotype comparison
3. **Nestin-Mut vs Emx1-Mut**: Mutant-to-mutant genotype effect

**Usage:**
```bash
# Default (fast, metaprofiles only)
sbatch 6_create_custom_comparisons.sh

# Full mode with individual plots
SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 6_create_custom_comparisons.sh
```

---

## Master Script (`0_RUN_SPLICING_GENES_ANALYSIS.sh`)

Runs the complete workflow in sequence (Steps 1-4):

**Inputs:**
- All inputs from individual steps

**Outputs:**
- All outputs from steps 1-4
- Combined logs and error files

**Execution:**
```bash
sbatch 0_RUN_SPLICING_GENES_ANALYSIS.sh
```

**Note:** Steps 5 and 6 are optional and run separately for alternative visualizations.

---

## Key Differences from Other Workflows

### Compared to Specific_CREs Workflow:
1. **Data Source**: Uses literature CRE-gene correlations (Table 16) instead of Signac peak-gene links
2. **Gene Focus**: Specifically targets splicing machinery genes (1,138 genes from Reactome)
3. **CRE Count**: ~832 total CREs, ~772 GABA CREs (vs. thousands in Specific_CREs)
4. **Performance**: Default mode skips individual plots for faster analysis (3-5 min vs. 30-50 min)
5. **Comparisons**: Includes custom cross-genotype comparisons (Step 6)

### Compared to Exclusive_CREs Workflow:
1. **Exclusivity**: GABA CREs are NOT exclusive (overlap with other cell types allowed in Steps 3-6)
2. **GABA-Specific Analysis**: Step 1 generates GABA-specific CREs, but main analysis uses all GABA CREs
3. **Literature-Based**: All CREs from published correlations, not from our Signac analysis

---

## Workflow Execution Order

### Quick Analysis (Recommended):
```bash
# Step 1: Extract CREs from literature
python 1_extract_splicing_gene_CREs.py

# Step 2: Convert to BED format
python 2_convert_splicing_CREs_to_bed.py

# Step 3: Create metaprofiles (fast mode, 3-5 minutes)
sbatch 3_create_splicing_profiles.sh
```

### Full Analysis with Individual Plots:
```bash
# Steps 1-2 same as above

# Step 3: Full mode with parallelization
SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 3_create_splicing_profiles.sh

# Optional: Step 5 for alternative BigWig visualization
SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 5_visualize_bigwig_signal.sh

# Optional: Step 6 for custom cross-genotype comparisons
SKIP_INDIVIDUAL=0 PARALLEL_JOBS=8 sbatch 6_create_custom_comparisons.sh
```

### Complete Workflow (Steps 1-4):
```bash
sbatch 0_RUN_SPLICING_GENES_ANALYSIS.sh
```

---

## Performance Guide

See [PERFORMANCE_OPTIONS.md](PERFORMANCE_OPTIONS.md) for detailed performance optimization strategies.

**Summary:**
- **Default mode**: Metaprofiles only (3-5 minutes) - **NEW DEFAULT**
- **Full mode**: All plots, sequential (30-50 minutes)
- **Optimized mode**: All plots, parallel (8 minutes)

**Environment Variables:**
- `SKIP_INDIVIDUAL=0`: Enable individual CRE plots (default: skip)
- `PARALLEL_JOBS=8`: Use 8 parallel processes
- `INDIVIDUAL_DPI=150`: Set DPI for individual plots (default: 150)

---

## File Size Estimates

**Input Files:**
- `table_16.txt`: 567 MB (literature data)
- `extracted_genes_final.csv`: ~50 KB
- BigWig files: ~100-500 MB each

**Output Files:**
- TSV files (Step 1): ~5-10 MB total
- BED files (Step 2): ~50-100 KB
- deepTools matrices: ~10-20 MB per matrix
- Metaprofiles (300 DPI): ~500 KB - 1 MB each
- Individual plots (150 DPI): ~100-200 KB each
- Individual plots (300 DPI): ~200-400 KB each

**Total Storage (with all individual plots):**
- ~5-10 GB for full analysis with individual plots
- ~100-200 MB for metaprofiles only

---

## Notes

1. **Default Behavior Change**: As of version 1.3, individual plots are **skipped by default** to save time. Use `SKIP_INDIVIDUAL=0` to create them.

2. **Literature Data**: All CREs come from published CRE-gene correlations (Table 16), not from Signac peak-gene linkage analysis.

3. **Cell Type Classification**: Based solely on SubType column from Table 16 (independent of external BED files).

4. **CRE Classification Used in Analysis**:
   - **ALL CREs**: CREs linked to splicing genes from all cell types
   - **GABA CREs**: Subset of ALL CREs that appear in GABA/hippocampal cell types
   - **GABA-specific CREs**: CREs appearing ONLY in GABA/hippocampal SubTypes (generated but not used in main analysis)
   
   **Important Note**: Main analysis (Steps 3-6) uses **ALL CREs** with GABA BigWig files, not limited to GABA-specific CREs

5. **Failed Sample**: Emx1-Ctrl is excluded from analysis (failed sample quality).

6. **Parallel Processing**: Steps 3, 5, and 6 support multiprocessing for faster individual plot generation.

7. **Reusable Matrices**: deepTools matrices (.gz/.tab) can be reused for additional visualizations without recomputing signal.