# Striatum_EP CRE Analysis Pipeline

This pipeline analyzes cis-regulatory elements (CREs) from `Striatum_EP.txt` that are associated with splicing-related genes.

**Key Features:**
- Uses CRE-gene linkages from `Striatum_EP.txt`
- Filters for 1,138 splicing-related genes
- Generates heatmaps and metaprofiles for GABA neurons

## Results Summary

Based on the latest run:
- **Input File**: `Striatum_EP.txt`
- **Splicing Genes**: 1,138
- **Linked CREs Found**: ~2,761

## Pipeline Steps

### Step 1: Extract CREs (`1_extract_Striatum_EP.py`)
- Parses `Striatum_EP.txt`
- Filters for splicing genes
- Outputs:
  - `output/Striatum_EP_splicing_genes.tsv`
  - `output/Striatum_EP_splicing_genes.bed`

### Step 4: Create Heatmaps (`4_create_heatmaps.sh`)
- Uses deepTools to generate heatmaps and metaprofiles
- Input: `Striatum_EP_splicing_genes.bed` + GABA BigWigs
- Output: `output/heatmaps_deeptools/`

### Step 5: Visualize Comparisons (`5_visualize_comparisons.py`)
- Generates publication-quality metaprofiles and statistics
- Output: `output/custom_comparisons/`

## Usage

```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/Striatum_EP_analysis

# Run complete pipeline
chmod +x 0_RUN_ANALYSIS.sh
./0_RUN_ANALYSIS.sh
```
