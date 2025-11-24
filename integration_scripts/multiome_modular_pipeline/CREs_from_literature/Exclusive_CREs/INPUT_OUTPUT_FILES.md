# Input and Output Files (Exclusive CREs Workflow)

## Master Pipeline Script (`0_RUN_COMPLETE_EXCLUSIVE_PIPELINE.sh`)
**Purpose:** Master script that submits all pipeline steps with proper SLURM dependencies
**Inputs:** None (submits other scripts)
**Outputs:** SLURM job submissions with dependencies
- Job 1: `1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh` (extracts CREs and creates deepTools heatmaps)
- Job 2: `2_fold_enrichment_exclusive_analysis.sh` (depends on Job 1 - calculates fold enrichment)
- Job 3: `3_compare_cell_type_specific_CREs.sh` (depends on Job 1 - creates Python comparison plots)

## 1. Main Analysis (`1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh`)
**Purpose:** Complete cell-type specificity analysis with mutually exclusive CRE sets
**Inputs:**
- `../data/table_1.xlsx` (ENCODE data table 1 - sample and dissection summary)
- `../data/table_2.tsv` (ENCODE data table 2 - cell metadata with sample assignments)
- `../data/table_8.txt` (ENCODE data table 8 - cell type assignment of cCREs, 389MB)
- `../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw` (BigWig files with ATAC-seq signal data)

**Internal Scripts Called:**
- `1a_extract_cell_type_specific_CREs.py` (extracts mutually exclusive CRE sets)
- `1b_create_heatmaps_specific_CREs_deeptools.sh` (creates deepTools heatmaps and metaprofiles)

**Outputs:**
- `output/GABA_specific_CREs.bed` (GABA-specific CREs, mutually exclusive)
- `output/Excitatory_specific_CREs.bed` (Excitatory-specific CREs, mutually exclusive)
- `output/cell_type_specific_CREs_summary.txt` (summary report of CRE extraction)
- `logs/1_RUN_EXCLUSIVE_CRES_ANALYSIS.log` (execution log)
- `logs/1_RUN_EXCLUSIVE_CRES_ANALYSIS.err` (error log)

**deepTools Output Directory:** `output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/`
- `heatmap_GABA_specific.png` (POSITIVE CONTROL - expect STRONG red signal)
- `heatmap_Excitatory_specific.png` (NEGATIVE CONTROL - expect PALE/white signal)
- `metaprofile_GABA_specific.png` (POSITIVE CONTROL - expect HIGH curves)
- `metaprofile_Excitatory_specific.png` (NEGATIVE CONTROL - expect LOW/flat curves)
- `comparison_cell_type_specific.png` (side-by-side comparison with fold enrichment)
- `cell_type_specificity_statistics.txt` (quantitative metrics and interpretation)
- `matrix_GABA_specific.gz` and `matrix_Excitatory_specific.gz` (compressed signal matrices)
- `README.txt` (detailed methodology and interpretation guide)
- Genotype-specific analyses: `heatmap_GABA_specific_nestin.png`, `heatmap_GABA_specific_emx1.png`, etc.

## 2. Fold Enrichment Analysis (`2_fold_enrichment_exclusive_analysis.sh`)
**Purpose:** Post-processing analysis on deepTools matrices to calculate fold-enrichment statistics
**Inputs:**
- `output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/matrix_GABA_specific.gz`
- `output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/matrix_Excitatory_specific.gz`

**Internal Scripts Called:**
- `2_fold_enrichment_exclusive_analysis.py` (performs statistical analysis and creates figures)

**Outputs:**
- `output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/fold_enrichment_statistics.txt` (quantitative statistics with biological interpretation)
- `output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/comparison_fold_enrichment.png` (4-panel comparison figure)
- `output/GABA_DEG_analysis/heatmaps_specific_CREs_deeptools/fold_enrichment_by_condition.png` (detailed fold-enrichment bar chart)
- `logs/2_fold_enrichment_exclusive_analysis.log` (execution log)
- `logs/2_fold_enrichment_exclusive_analysis.err` (error log)

## 3. Python Comparison Analysis (`3_compare_cell_type_specific_CREs.sh`)
**Purpose:** Alternative Python-based analysis using pyBigWig for signal extraction
**Inputs:**
- `output/GABA_specific_CREs.bed` (from step 1)
- `output/Excitatory_specific_CREs.bed` (from step 1)
- `../../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_*.bw` (BigWig files)

**Internal Scripts Called:**
- `3_compare_cell_type_specific_CREs.py` (generates Python-based heatmaps and metaprofiles)

**Outputs:**
**Output Directory:** `output/GABA_DEG_analysis/heatmaps_specific_CREs/`
- `heatmap_GABA_specific.png` (POSITIVE CONTROL - expect HIGH signal)
- `heatmap_Excitatory_specific.png` (NEGATIVE CONTROL - expect LOW signal)
- `metaprofile_GABA_specific.png` (POSITIVE CONTROL - expect HIGH curves)
- `metaprofile_Excitatory_specific.png` (NEGATIVE CONTROL - expect LOW curves)
- `comparison_cell_type_specific.png` (combined comparison plot with metaprofiles, bar charts, and fold enrichment)
- `cell_type_specificity_statistics.txt` (detailed statistics file with signal comparisons and interpretations)
- `logs/3_compare_cell_type_specific_CREs.log` (execution log)
- `logs/3_compare_cell_type_specific_CREs.err` (error log)

## Key Features of This Pipeline:

1. **Mutually Exclusive CREs**: Creates GABA-specific and Excitatory-specific CREs with 0% overlap (vs previous 60% overlap)
2. **Dual Analysis Approaches**:
   - deepTools (fast, parallel C++ implementation) in step 1
   - Python/pyBigWig (alternative implementation) in step 3
3. **Positive/Negative Controls**:
   - GABA-specific CREs = POSITIVE CONTROL (expect HIGH signal in GABA samples)
   - Excitatory-specific CREs = NEGATIVE CONTROL (expect LOW signal in GABA samples)
4. **Comprehensive Statistics**: Fold-enrichment analysis with biological interpretation (≥3x = strong specificity)
5. **SLURM Integration**: Master script handles job dependencies automatically

## Data Requirements:

**ENCODE Data Tables:**
- `table_1.xlsx`: Sample and dissection summary
- `table_2.tsv`: Cell metadata with sample assignments
- `table_8.txt`: Cell type assignment of cCREs (389MB file)

**BigWig Files Required:**
- `GABA_Nestin-Ctrl.bw`
- `GABA_Nestin-Mut.bw`
- `GABA_Emx1-Ctrl.bw`
- `GABA_Emx1-Mut.bw`

## Workflow Order:
1. Run `0_RUN_COMPLETE_EXCLUSIVE_PIPELINE.sh` (master script) OR run steps individually:
2. Run `1_RUN_EXCLUSIVE_CRES_ANALYSIS.sh` to extract mutually exclusive CREs and generate deepTools outputs
3. Run `2_fold_enrichment_exclusive_analysis.sh` to calculate fold enrichment from deepTools matrices
4. Run `3_compare_cell_type_specific_CREs.sh` to generate Python-based comparison plots (optional)

## Expected Results:
- **GABA-specific CREs**: Strong red signal in heatmaps, high curves in metaprofiles
- **Excitatory-specific CREs**: Pale/white signal in heatmaps, low/flat curves in metaprofiles
- **Fold enrichment**: ≥3x enrichment indicating strong cell-type specificity
- **Clear visual difference** between the two CRE sets demonstrating cell-type specificity
