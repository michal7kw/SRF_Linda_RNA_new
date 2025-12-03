#!/usr/bin/env python3
"""
Extract CREs Associated with GABA Differentially Expressed Genes

This script identifies CREs (cis-regulatory elements) that are associated with
GABA-specific differentially expressed genes (DEGs), using literature CRE-gene
correlations from Table 16.

INPUT FILES:
1. GABA DEG lists:
   - Up-regulated DEGs: Cluster_GABA_vs_Rest_up_significant.csv
   - Down-regulated DEGs: Cluster_GABA_vs_Rest_down_significant.csv

2. CRE-gene correlations (LITERATURE DATA):
   ../data/table_16.txt (567 MB, literature CRE-gene correlations)

CLASSIFICATION:
- Cell type classification based on SubType column from Table 16
- GABA CREs: CREs that appear in hippocampal/GABAergic SubTypes

OUTPUT FILES:
1. GABA_DEGs_up_CREs.tsv - CREs linked to up-regulated DEGs
2. GABA_DEGs_down_CREs.tsv - CREs linked to down-regulated DEGs
3. GABA_DEGs_up_CREs.bed - BED format for up-regulated DEGs CREs
4. GABA_DEGs_down_CREs.bed - BED format for down-regulated DEGs CREs

Usage:
    python 1_extract_gaba_deg_CREs.py
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime

os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_GABA_DEGs_paper")

print("="*80)
print("EXTRACT CREs ASSOCIATED WITH GABA DIFFERENTIALLY EXPRESSED GENES")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# =============================================================================
# Configuration
# =============================================================================

# Input files - GABA DEGs
DEG_UP_FILE = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/biomarkers/sig_deg_lists/GABA/Cluster_GABA_vs_Rest_up_significant.csv"
DEG_DOWN_FILE = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/biomarkers/sig_deg_lists/GABA/Cluster_GABA_vs_Rest_down_significant.csv"
TABLE16_FILE = "../data/table_16.txt"  # Literature CRE-gene correlations

# Statistical filters (same as splicing genes analysis)
FDR_THRESHOLD = 0.05
PCC_THRESHOLD = 0.2

# Output directory
OUTPUT_DIR = "./output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# Step 1: Load GABA DEG Lists
# =============================================================================

print("="*80)
print("STEP 1: Loading GABA DEG lists")
print("-"*80)

# Load up-regulated DEGs
deg_up = pd.read_csv(DEG_UP_FILE)
deg_up_genes = set(deg_up['names'].dropna())
print(f"Up-regulated DEGs: {len(deg_up_genes)} genes")
print(f"  Sample genes: {', '.join(list(deg_up_genes)[:10])}")

# Load down-regulated DEGs
deg_down = pd.read_csv(DEG_DOWN_FILE)
deg_down_genes = set(deg_down['names'].dropna())
print(f"\nDown-regulated DEGs: {len(deg_down_genes)} genes")
print(f"  Sample genes: {', '.join(list(deg_down_genes)[:10])}")

# Convert to uppercase for case-insensitive matching
deg_up_genes_upper = {g.upper() for g in deg_up_genes}
deg_down_genes_upper = {g.upper() for g in deg_down_genes}
print()

# =============================================================================
# Step 2: Load CRE-Gene Correlations from Table 16 (LITERATURE DATA)
# =============================================================================

print("="*80)
print("STEP 2: Loading CRE-gene correlations from Table 16 (literature data)")
print("-"*80)
print("This may take 1-2 minutes (567 MB file)...")

table16 = pd.read_csv(TABLE16_FILE, sep='\t')
print(f"Loaded {len(table16):,} CRE-gene correlation entries from literature")

# Create uppercase version of Gene column for matching
table16['Gene_upper'] = table16['Gene'].str.upper()

print("\nColumns available:")
for col in table16.columns:
    print(f"  {col}")
print()

# =============================================================================
# Step 3: Define GABA/Hippocampal Cell Type Filter
# =============================================================================

print("="*80)
print("STEP 3: Setting up cell type classification")
print("-"*80)

# Import hippocampal GABAergic cell types from shared helper module
# This uses EXACT matching against 46 validated cell types (not keyword matching)
# Excludes glutamatergic neurons like CA1GL, CA3GL, DGGR
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'helpers'))
from gaba_cell_types import HIPPOCAMPAL_GABA_CELLTYPES, is_gaba_subtype

# Alias for compatibility with existing code
is_hippocampal_or_gaba = is_gaba_subtype

print(f"Using {len(HIPPOCAMPAL_GABA_CELLTYPES)} validated hippocampal GABAergic cell types")
print("(EXACT matching - excludes CA1GL, CA3GL, DGGR, etc.)")
print()

# =============================================================================
# Step 4: Extract CREs for Up-regulated DEGs
# =============================================================================

print("="*80)
print("STEP 4: Extracting CREs for up-regulated DEGs")
print("-"*80)

# Filter for up-regulated DEGs
cre_links_up = table16[table16['Gene_upper'].isin(deg_up_genes_upper)].copy()
print(f"Initial CRE-gene links for up-DEGs: {len(cre_links_up):,}")

# Apply statistical filters
cre_links_up = cre_links_up[
    (cre_links_up['FDR'] < FDR_THRESHOLD) &
    (cre_links_up['PCC'].abs() > PCC_THRESHOLD)
].copy()
print(f"After statistical filters (FDR < {FDR_THRESHOLD}, |PCC| > {PCC_THRESHOLD}): {len(cre_links_up):,}")

# Apply GABA cell type filter
cre_links_up['is_gaba'] = cre_links_up['SubType'].apply(is_hippocampal_or_gaba)
cre_links_up_gaba = cre_links_up[cre_links_up['is_gaba']].copy()
print(f"After GABA cell type filter: {len(cre_links_up_gaba):,}")

# Get statistics
unique_cres_up = cre_links_up_gaba['cCRE1'].nunique()
unique_genes_up = cre_links_up_gaba['Gene'].nunique()
print(f"\nUnique CREs (up-DEGs): {unique_cres_up:,}")
print(f"Unique genes with CREs (up-DEGs): {unique_genes_up:,}")
print()

# =============================================================================
# Step 5: Extract CREs for Down-regulated DEGs
# =============================================================================

print("="*80)
print("STEP 5: Extracting CREs for down-regulated DEGs")
print("-"*80)

# Filter for down-regulated DEGs
cre_links_down = table16[table16['Gene_upper'].isin(deg_down_genes_upper)].copy()
print(f"Initial CRE-gene links for down-DEGs: {len(cre_links_down):,}")

# Apply statistical filters
cre_links_down = cre_links_down[
    (cre_links_down['FDR'] < FDR_THRESHOLD) &
    (cre_links_down['PCC'].abs() > PCC_THRESHOLD)
].copy()
print(f"After statistical filters (FDR < {FDR_THRESHOLD}, |PCC| > {PCC_THRESHOLD}): {len(cre_links_down):,}")

# Apply GABA cell type filter
cre_links_down['is_gaba'] = cre_links_down['SubType'].apply(is_hippocampal_or_gaba)
cre_links_down_gaba = cre_links_down[cre_links_down['is_gaba']].copy()
print(f"After GABA cell type filter: {len(cre_links_down_gaba):,}")

# Get statistics
unique_cres_down = cre_links_down_gaba['cCRE1'].nunique()
unique_genes_down = cre_links_down_gaba['Gene'].nunique()
print(f"\nUnique CREs (down-DEGs): {unique_cres_down:,}")
print(f"Unique genes with CREs (down-DEGs): {unique_genes_down:,}")
print()

# =============================================================================
# Step 6: Convert to BED format and Save
# =============================================================================

print("="*80)
print("STEP 6: Converting to BED format and saving output files")
print("-"*80)

def parse_coordinate(coord_str):
    """Parse coordinate string to chr, start, end"""
    parts = coord_str.split('_')
    if len(parts) == 3:
        return parts[0], int(parts[1]), int(parts[2])
    else:
        raise ValueError(f"Invalid coordinate format: {coord_str}")

def create_bed_from_cre_links(cre_links_df, tsv_output, bed_output, description):
    """Create BED file from CRE links DataFrame"""
    print(f"\nProcessing {description}...")

    if len(cre_links_df) == 0:
        print(f"  WARNING: No CRE links found for {description}")
        return None

    # Save TSV file
    cre_links_df.to_csv(tsv_output, sep='\t', index=False)
    print(f"  Saved TSV: {tsv_output} ({len(cre_links_df):,} rows)")

    # Parse coordinates
    cre_links_df[['chr', 'start', 'end']] = cre_links_df['Coordinate1'].apply(
        lambda x: pd.Series(parse_coordinate(x))
    )

    # Extract unique CREs with coordinates
    unique_cres = cre_links_df[['chr', 'start', 'end', 'cCRE1']].drop_duplicates()

    # Sort by chromosome and start position
    chr_order = [f'chr{i}' for i in range(1, 20)] + ['chrX', 'chrY', 'chrM']
    unique_cres['chr'] = pd.Categorical(unique_cres['chr'], categories=chr_order, ordered=True)
    unique_cres = unique_cres.sort_values(['chr', 'start'])

    # Add score column (all 1000) and strand column (all .)
    unique_cres['score'] = 1000
    unique_cres['strand'] = '.'

    # Reorder columns for BED6 format
    bed_df = unique_cres[['chr', 'start', 'end', 'cCRE1', 'score', 'strand']]

    # Save BED file
    bed_df.to_csv(bed_output, sep='\t', header=False, index=False)
    print(f"  Saved BED: {bed_output} ({len(bed_df):,} unique CREs)")

    # Show first few entries
    print(f"\n  First 5 CREs:")
    print(bed_df.head().to_string(index=False, header=False))

    # Get genes
    genes = cre_links_df['Gene'].unique()
    print(f"\n  Genes represented: {len(genes)}")
    if len(genes) <= 20:
        print(f"  Gene list: {', '.join(sorted(genes))}")
    else:
        print(f"  Sample genes: {', '.join(sorted(genes)[:20])}...")

    return bed_df

# Process up-regulated DEGs
bed_up = create_bed_from_cre_links(
    cre_links_up_gaba,
    os.path.join(OUTPUT_DIR, "GABA_DEGs_up_CREs.tsv"),
    os.path.join(OUTPUT_DIR, "GABA_DEGs_up_CREs.bed"),
    "Up-regulated DEGs"
)

# Process down-regulated DEGs
bed_down = create_bed_from_cre_links(
    cre_links_down_gaba,
    os.path.join(OUTPUT_DIR, "GABA_DEGs_down_CREs.tsv"),
    os.path.join(OUTPUT_DIR, "GABA_DEGs_down_CREs.bed"),
    "Down-regulated DEGs"
)

# =============================================================================
# Step 7: Generate Summary Report
# =============================================================================

print("\n" + "="*80)
print("STEP 7: Creating summary report")
print("-"*80)

summary_file = os.path.join(OUTPUT_DIR, "SUMMARY_GABA_DEGs_CREs.txt")

with open(summary_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("SUMMARY: CREs ASSOCIATED WITH GABA DIFFERENTIALLY EXPRESSED GENES\n")
    f.write("="*80 + "\n\n")

    f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("INPUT FILES:\n")
    f.write("-"*80 + "\n")
    f.write(f"Up-regulated DEGs: {DEG_UP_FILE}\n")
    f.write(f"  Total genes: {len(deg_up_genes):,}\n\n")
    f.write(f"Down-regulated DEGs: {DEG_DOWN_FILE}\n")
    f.write(f"  Total genes: {len(deg_down_genes):,}\n\n")
    f.write(f"CRE-gene correlations: {TABLE16_FILE}\n")
    f.write(f"  Total links: {len(table16):,}\n\n")

    f.write("FILTERS APPLIED:\n")
    f.write("-"*80 + "\n")
    f.write(f"FDR threshold: < {FDR_THRESHOLD}\n")
    f.write(f"|PCC| threshold: > {PCC_THRESHOLD}\n")
    f.write(f"Cell type: GABA/Hippocampal (keywords: {', '.join(hippocampal_keywords)})\n\n")

    f.write("RESULTS:\n")
    f.write("="*80 + "\n\n")

    f.write("UP-REGULATED DEGs:\n")
    f.write("-"*80 + "\n")
    f.write(f"  Input DEGs: {len(deg_up_genes):,}\n")
    f.write(f"  DEGs with CRE links: {unique_genes_up:,} ({100*unique_genes_up/len(deg_up_genes):.1f}%)\n")
    f.write(f"  Unique CREs: {unique_cres_up:,}\n")
    f.write(f"  Total CRE-gene links: {len(cre_links_up_gaba):,}\n\n")

    f.write("DOWN-REGULATED DEGs:\n")
    f.write("-"*80 + "\n")
    f.write(f"  Input DEGs: {len(deg_down_genes):,}\n")
    f.write(f"  DEGs with CRE links: {unique_genes_down:,} ({100*unique_genes_down/len(deg_down_genes):.1f}%)\n")
    f.write(f"  Unique CREs: {unique_cres_down:,}\n")
    f.write(f"  Total CRE-gene links: {len(cre_links_down_gaba):,}\n\n")

    f.write("OUTPUT FILES:\n")
    f.write("-"*80 + "\n")
    f.write("  output/GABA_DEGs_up_CREs.tsv - CRE links for up-regulated DEGs\n")
    f.write("  output/GABA_DEGs_up_CREs.bed - BED format for up-regulated DEGs CREs\n")
    f.write("  output/GABA_DEGs_down_CREs.tsv - CRE links for down-regulated DEGs\n")
    f.write("  output/GABA_DEGs_down_CREs.bed - BED format for down-regulated DEGs CREs\n\n")

    f.write("="*80 + "\n")
    f.write(f"Analysis completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write("="*80 + "\n")

print(f"Saved summary: {summary_file}")
print()

# =============================================================================
# Completion
# =============================================================================

print("="*80)
print("ANALYSIS COMPLETE!")
print("="*80)
print()
print(f"Output directory: {OUTPUT_DIR}/")
print()
print("Generated files:")
print(f"  1. GABA_DEGs_up_CREs.tsv - {len(cre_links_up_gaba):,} CRE-gene links")
print(f"  2. GABA_DEGs_up_CREs.bed - {unique_cres_up:,} unique CREs")
print(f"  3. GABA_DEGs_down_CREs.tsv - {len(cre_links_down_gaba):,} CRE-gene links")
print(f"  4. GABA_DEGs_down_CREs.bed - {unique_cres_down:,} unique CREs")
print(f"  5. SUMMARY_GABA_DEGs_CREs.txt - Summary report")
print()
print("Key Statistics:")
print("-"*80)
print(f"Up-regulated DEGs: {len(deg_up_genes):,} genes -> {unique_genes_up:,} with CREs ({unique_cres_up:,} CREs)")
print(f"Down-regulated DEGs: {len(deg_down_genes):,} genes -> {unique_genes_down:,} with CREs ({unique_cres_down:,} CREs)")
print()
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
