#!/usr/bin/env python3
"""
Extract CREs Associated with Splicing Genes from Literature Data (CORRECTED)

This script identifies CREs (cis-regulatory elements) that are associated with
splicing-related genes from Reactome pathways, using literature CRE-gene
correlations from Table 16.

CORRECT APPROACH:
- Uses Table 16 (literature CRE-gene correlations), NOT Signac peak-gene links
- Filters for ALL splicing genes (not just genes_inter.txt)
- Follows same methodology as link_genes_to_cres.py

INPUT FILES:
1. Splicing genes list:
   /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/splicing_genes/extracted_genes_final.csv

2. CRE-gene correlations (LITERATURE DATA):
   ../data/table_16.txt (567 MB, literature CRE-gene correlations)

CLASSIFICATION:
- Cell type classification based ONLY on SubType column from Table 16
- GABA-specific CREs: CREs that appear exclusively in hippocampal/GABAergic SubTypes
- Independent from external BED files

OUTPUT FILES:
1. splicing_genes_CREs_all_celltypes.tsv
   - All CREs linked to splicing genes (any cell type)

2. splicing_genes_CREs_GABA.tsv
   - CREs linked to splicing genes in GABA/interneuron cell types

3. splicing_genes_CREs_GABA_specific.tsv
   - GABA-specific CREs linked to splicing genes

Usage:
    python 1_extract_splicing_gene_CREs_CORRECTED.py
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime

os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/splicing_genes")

print("="*80)
print("EXTRACT CREs ASSOCIATED WITH SPLICING GENES FROM LITERATURE (CORRECTED)")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# =============================================================================
# Configuration
# =============================================================================

# Input files
SPLICING_GENES_FILE = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/splicing_genes/extracted_genes_final.csv"
TABLE16_FILE = "../data/table_16.txt"  # Literature CRE-gene correlations

# Statistical filters (same as link_genes_to_cres.py)
FDR_THRESHOLD = 0.05
PCC_THRESHOLD = 0.2

# Output directory
OUTPUT_DIR = "./output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Output files
OUTPUT_ALL = os.path.join(OUTPUT_DIR, "splicing_genes_CREs_all_celltypes.tsv")
OUTPUT_GABA = os.path.join(OUTPUT_DIR, "splicing_genes_CREs_GABA.tsv")
OUTPUT_GABA_SPECIFIC = os.path.join(OUTPUT_DIR, "splicing_genes_CREs_GABA_specific.tsv")

# =============================================================================
# Step 1: Load Splicing Genes List
# =============================================================================

print("="*80)
print("STEP 1: Loading splicing genes list")
print("-"*80)

splicing_genes = pd.read_csv(SPLICING_GENES_FILE)
print(f"Loaded {len(splicing_genes)} splicing genes from Reactome/GO")

# Extract unique gene symbols (filter out non-gene entries like ChEBI)
splicing_gene_symbols = set(splicing_genes['gene_symbol'].dropna())
splicing_gene_symbols = {g for g in splicing_gene_symbols if not g.startswith('[')}

# Convert to uppercase for case-insensitive matching
splicing_gene_symbols_upper = {g.upper() for g in splicing_gene_symbols}

print(f"Unique gene symbols: {len(splicing_gene_symbols)}")
print(f"Sample genes: {', '.join(list(splicing_gene_symbols)[:10])}")
print()

# =============================================================================
# Step 2: Load CRE-Gene Correlations from Table 16 (LITERATURE DATA)
# =============================================================================

print("="*80)
print("STEP 2: Loading CRE-gene correlations from Table 16 (literature data)")
print("-"*80)
print("This may take 1-2 minutes (567 MB file)...")

table16 = pd.read_csv(TABLE16_FILE, sep='\t')
print(f"✓ Loaded {len(table16):,} CRE-gene correlation entries from literature")

print("\nColumns available:")
for col in table16.columns:
    print(f"  {col}")

print("\nSample data:")
print(table16.head(3))
print()

# =============================================================================
# Step 3: Filter for Splicing Genes
# =============================================================================

print("="*80)
print("STEP 3: Filtering for splicing genes")
print("-"*80)

# Create uppercase version of Gene column for matching
table16['Gene_upper'] = table16['Gene'].str.upper()
splicing_cre_links_all = table16[table16['Gene_upper'].isin(splicing_gene_symbols_upper)].copy()

print(f"CRE-gene links for splicing genes: {len(splicing_cre_links_all):,}")

# Get unique CREs and genes
unique_cres_before_filter = splicing_cre_links_all['cCRE1'].nunique()
unique_genes_before_filter = splicing_cre_links_all['Gene'].nunique()
print(f"Unique CREs (before filtering): {unique_cres_before_filter:,}")
print(f"Unique genes (before filtering): {unique_genes_before_filter:,}")
print()

# =============================================================================
# Step 4: Apply Statistical Filters (FDR < 0.05, |PCC| > 0.2)
# =============================================================================

print("="*80)
print("STEP 4: Applying statistical filters")
print("-"*80)

print(f"Before filtering: {len(splicing_cre_links_all):,} links")

# Check distribution of statistics
print("\nStatistical measures:")
print(f"  PCC (Pearson correlation): min={splicing_cre_links_all['PCC'].min():.3f}, "
      f"median={splicing_cre_links_all['PCC'].median():.3f}, max={splicing_cre_links_all['PCC'].max():.3f}")
print(f"  FDR: min={splicing_cre_links_all['FDR'].min():.6f}, "
      f"median={splicing_cre_links_all['FDR'].median():.6f}, max={splicing_cre_links_all['FDR'].max():.6f}")

# Apply filters (same as link_genes_to_cres.py)
splicing_cre_links_all = splicing_cre_links_all[
    (splicing_cre_links_all['FDR'] < FDR_THRESHOLD) &
    (splicing_cre_links_all['PCC'].abs() > PCC_THRESHOLD)
].copy()

print(f"\nAfter filtering (FDR < {FDR_THRESHOLD}, |PCC| > {PCC_THRESHOLD}):")
print(f"  Significant links: {len(splicing_cre_links_all):,}")

# Get unique CREs and genes after filtering
unique_cres_all = splicing_cre_links_all['cCRE1'].nunique()
unique_genes_all = splicing_cre_links_all['Gene'].nunique()
print(f"  Unique CREs: {unique_cres_all:,}")
print(f"  Unique genes: {unique_genes_all:,}")
print()

# =============================================================================
# Step 5: Filter by Cell Type (Hippocampal/GABAergic)
# =============================================================================

print("="*80)
print("STEP 5: Filtering by cell type (hippocampal/GABAergic)")
print("-"*80)

# Define hippocampal/GABAergic subtypes (same as link_genes_to_cres.py)
hippocampal_keywords = [
    'CA1', 'CA2', 'CA3', 'DG', 'DGNBL', 'GRC',  # Hippocampal
    'LAMP5', 'LAMP', 'VIP', 'SST', 'PV', 'PVGA', 'SSTGA', 'VIPGA', 'LAMGA',  # GABAergic interneurons
    'GABA', 'INH', 'interneuron'
]

# Filter for hippocampal/GABAergic cell types
def is_hippocampal_or_gaba(subtype):
    if pd.isna(subtype):
        return False
    subtype_str = str(subtype).upper()
    return any(keyword.upper() in subtype_str for keyword in hippocampal_keywords)

splicing_cre_links_all['is_hippo_gaba'] = splicing_cre_links_all['SubType'].apply(is_hippocampal_or_gaba)

print(f"Links in hippocampal/GABAergic cell types:")
print(f"  Total links: {splicing_cre_links_all['is_hippo_gaba'].sum():,} / {len(splicing_cre_links_all):,}")
print(f"  ({100 * splicing_cre_links_all['is_hippo_gaba'].sum() / len(splicing_cre_links_all):.1f}%)")

# Filter for hippocampal/GABAergic
splicing_cre_links_gaba = splicing_cre_links_all[splicing_cre_links_all['is_hippo_gaba']].copy()

unique_cres_gaba = splicing_cre_links_gaba['cCRE1'].nunique()
unique_genes_gaba = splicing_cre_links_gaba['Gene'].nunique()
print(f"\nGABA CRE-gene links:")
print(f"  Links: {len(splicing_cre_links_gaba):,}")
print(f"  Unique CREs: {unique_cres_gaba:,}")
print(f"  Unique genes: {unique_genes_gaba:,}")
print()

# =============================================================================
# Step 6: Identify GABA-Specific CREs (from Table 16 SubType data only)
# =============================================================================

print("="*80)
print("STEP 6: Identifying GABA-specific CREs from SubType data")
print("-"*80)

# GABA-specific CREs are those that appear ONLY in hippocampal/GABAergic cell types
# within the splicing genes CRE links from Table 16

# For each CRE, check if it appears in ANY non-GABA cell type
# Step 1: Get all CRE IDs that appear in splicing gene links
all_splicing_cres = set(splicing_cre_links_all['cCRE1'])

# Step 2: For each CRE, check if it appears in non-GABA cell types
# We'll look at the full filtered Table 16 data (all splicing gene CREs)
cre_in_non_gaba = set()
for cre_id in all_splicing_cres:
    cre_links = splicing_cre_links_all[splicing_cre_links_all['cCRE1'] == cre_id]
    # Check if any link is to a non-GABA cell type
    has_non_gaba = (~cre_links['is_hippo_gaba']).any()
    if has_non_gaba:
        cre_in_non_gaba.add(cre_id)

# GABA-specific CREs are those NOT in the non-GABA set
gaba_specific_cre_ids = all_splicing_cres - cre_in_non_gaba

print(f"Total CREs analyzed: {len(all_splicing_cres):,}")
print(f"CREs appearing in non-GABA cell types: {len(cre_in_non_gaba):,}")
print(f"CREs exclusive to GABA/hippocampal cell types: {len(gaba_specific_cre_ids):,}")
print()

# Filter GABA links for GABA-specific CREs only
splicing_cre_links_gaba_specific = splicing_cre_links_gaba[
    splicing_cre_links_gaba['cCRE1'].isin(gaba_specific_cre_ids)
].copy()

unique_cres_gaba_specific = splicing_cre_links_gaba_specific['cCRE1'].nunique()
unique_genes_gaba_specific = splicing_cre_links_gaba_specific['Gene'].nunique()

print(f"CRE-gene links for GABA-specific CREs:")
print(f"  Links: {len(splicing_cre_links_gaba_specific):,}")
print(f"  Unique CREs: {unique_cres_gaba_specific:,}")
print(f"  Unique genes: {unique_genes_gaba_specific:,}")
print()

# =============================================================================
# Step 7: Save Output Files
# =============================================================================

print("="*80)
print("STEP 7: Saving output files")
print("-"*80)

# Output 1: All cell types
splicing_cre_links_all.to_csv(OUTPUT_ALL, sep='\t', index=False)
print(f"✓ Saved: {OUTPUT_ALL}")
print(f"  Rows: {len(splicing_cre_links_all):,}")
print(f"  Unique CREs: {unique_cres_all:,}")
print(f"  Unique genes: {unique_genes_all:,}")
print()

# Output 2: GABA CREs
splicing_cre_links_gaba.to_csv(OUTPUT_GABA, sep='\t', index=False)
print(f"✓ Saved: {OUTPUT_GABA}")
print(f"  Rows: {len(splicing_cre_links_gaba):,}")
print(f"  Unique CREs: {unique_cres_gaba:,}")
print(f"  Unique genes: {unique_genes_gaba:,}")
print()

# Output 3: GABA-specific CREs
splicing_cre_links_gaba_specific.to_csv(OUTPUT_GABA_SPECIFIC, sep='\t', index=False)
print(f"✓ Saved: {OUTPUT_GABA_SPECIFIC}")
print(f"  Rows: {len(splicing_cre_links_gaba_specific):,}")
print(f"  Unique CREs: {unique_cres_gaba_specific:,}")
print(f"  Unique genes: {unique_genes_gaba_specific:,}")
print()

# =============================================================================
# Step 8: Generate Summary Statistics
# =============================================================================

print("="*80)
print("STEP 8: Summary statistics")
print("-"*80)

print("\nSplicing Genes Coverage:")
print("-"*40)
print(f"Total splicing genes in list: {len(splicing_gene_symbols):,}")
print(f"Genes with CREs (all cell types): {unique_genes_all:,} ({100*unique_genes_all/len(splicing_gene_symbols):.1f}%)")
print(f"Genes with GABA CREs: {unique_genes_gaba:,} ({100*unique_genes_gaba/len(splicing_gene_symbols):.1f}%)")
print(f"Genes with GABA-specific CREs: {unique_genes_gaba_specific:,}")

print("\nCRE Statistics:")
print("-"*40)
print(f"Total CREs linked to splicing genes (all): {unique_cres_all:,}")
if unique_cres_all > 0:
    print(f"GABA CREs linked to splicing genes: {unique_cres_gaba:,} ({100*unique_cres_gaba/unique_cres_all:.1f}%)")
    print(f"GABA-specific CREs linked to splicing genes: {unique_cres_gaba_specific:,}")

# Top genes by number of CRE links
print("\nTop 20 Splicing Genes by Number of CRE Links (All Cell Types):")
print("-"*40)
top_genes = splicing_cre_links_all.groupby('Gene')['cCRE1'].nunique().sort_values(ascending=False).head(20)
for gene, count in top_genes.items():
    print(f"  {gene:15s}: {count:,} CREs")

# Cell type distribution
print("\nCell Type Distribution (Top 15):")
print("-"*40)
if 'SubType' in splicing_cre_links_all.columns:
    subtype_counts = splicing_cre_links_all['SubType'].value_counts().head(15)
    for subtype, count in subtype_counts.items():
        print(f"  {subtype:20s}: {count:,} links")

# =============================================================================
# Step 9: Create Summary Report
# =============================================================================

print("\n" + "="*80)
print("STEP 9: Creating summary report")
print("-"*80)

summary_file = os.path.join(OUTPUT_DIR, "SUMMARY_splicing_genes_CREs.txt")

with open(summary_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("SUMMARY: CREs ASSOCIATED WITH SPLICING GENES (FROM LITERATURE DATA)\n")
    f.write("="*80 + "\n\n")

    f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("INPUT FILES:\n")
    f.write("-"*80 + "\n")
    f.write(f"Splicing genes: {SPLICING_GENES_FILE}\n")
    f.write(f"  Total genes: {len(splicing_gene_symbols):,}\n\n")

    f.write(f"CRE-gene correlations: {TABLE16_FILE}\n")
    f.write(f"  Total links (literature data): {len(table16):,}\n\n")

    f.write("CLASSIFICATION METHOD:\n")
    f.write("-"*80 + "\n")
    f.write("Cell type classification: Based ONLY on SubType column from Table 16\n")
    f.write("GABA-specific CREs: CREs appearing exclusively in hippocampal/GABAergic SubTypes\n")
    f.write("  (independent from external BED files)\n\n")

    f.write("FILTERS APPLIED:\n")
    f.write("-"*80 + "\n")
    f.write(f"FDR threshold: < {FDR_THRESHOLD}\n")
    f.write(f"|PCC| threshold: > {PCC_THRESHOLD}\n")
    f.write(f"Cell types: Hippocampal and GABAergic (for GABA analysis)\n\n")

    f.write("OUTPUT FILES:\n")
    f.write("="*80 + "\n\n")

    f.write(f"1. ALL CELL TYPES: {OUTPUT_ALL}\n")
    f.write("-"*80 + "\n")
    f.write(f"   CRE-gene links: {len(splicing_cre_links_all):,}\n")
    f.write(f"   Unique CREs: {unique_cres_all:,}\n")
    f.write(f"   Unique genes: {unique_genes_all:,}\n")
    f.write(f"   Coverage: {100*unique_genes_all/len(splicing_gene_symbols):.1f}% of splicing genes\n\n")

    f.write(f"2. GABA CREs: {OUTPUT_GABA}\n")
    f.write("-"*80 + "\n")
    f.write(f"   CRE-gene links: {len(splicing_cre_links_gaba):,}\n")
    f.write(f"   Unique CREs: {unique_cres_gaba:,}\n")
    f.write(f"   Unique genes: {unique_genes_gaba:,}\n")
    f.write(f"   Coverage: {100*unique_genes_gaba/len(splicing_gene_symbols):.1f}% of splicing genes\n\n")

    f.write(f"3. GABA-SPECIFIC CREs: {OUTPUT_GABA_SPECIFIC}\n")
    f.write("-"*80 + "\n")
    f.write(f"   CRE-gene links: {len(splicing_cre_links_gaba_specific):,}\n")
    f.write(f"   Unique CREs: {unique_cres_gaba_specific:,}\n")
    f.write(f"   Unique genes: {unique_genes_gaba_specific:,}\n\n")

    f.write("KEY FINDINGS:\n")
    f.write("="*80 + "\n\n")

    f.write("1. Gene Coverage:\n")
    f.write(f"   - {unique_genes_all:,} of {len(splicing_gene_symbols):,} splicing genes have associated CREs\n")
    f.write(f"   - {unique_genes_gaba:,} genes have GABA CREs\n")
    f.write(f"   - {unique_genes_gaba_specific:,} genes have GABA-specific CREs\n\n")

    f.write("2. CRE Distribution:\n")
    f.write(f"   - Total CREs for splicing genes: {unique_cres_all:,}\n")
    if unique_cres_all > 0:
        f.write(f"   - GABA CREs: {unique_cres_gaba:,} ({100*unique_cres_gaba/unique_cres_all:.1f}%)\n")
        f.write(f"   - GABA-specific: {unique_cres_gaba_specific:,}\n\n")

    f.write("3. Top Genes (by number of CREs):\n")
    for i, (gene, count) in enumerate(top_genes.items(), 1):
        f.write(f"   {i:2d}. {gene:15s}: {count:,} CREs\n")

    f.write("\n" + "="*80 + "\n")
    f.write(f"Analysis completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write("="*80 + "\n")

print(f"✓ Saved: {summary_file}")
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
print(f"  1. {os.path.basename(OUTPUT_ALL)}")
print(f"     → {len(splicing_cre_links_all):,} CRE-gene links (all cell types)")
print(f"     → {unique_cres_all:,} unique CREs")
print()
print(f"  2. {os.path.basename(OUTPUT_GABA)}")
print(f"     → {len(splicing_cre_links_gaba):,} CRE-gene links (GABA cell types)")
print(f"     → {unique_cres_gaba:,} unique CREs")
print()
print(f"  3. {os.path.basename(OUTPUT_GABA_SPECIFIC)}")
print(f"     → {len(splicing_cre_links_gaba_specific):,} CRE-gene links (GABA-specific)")
print()
print(f"  4. SUMMARY_splicing_genes_CREs.txt")
print(f"     → Summary report with statistics")
print()
print("Key Results:")
print("-"*80)
print(f"Splicing genes analyzed: {len(splicing_gene_symbols):,}")
print(f"Genes with CREs (all): {unique_genes_all:,} ({100*unique_genes_all/len(splicing_gene_symbols):.1f}%)")
print(f"Total CREs identified: {unique_cres_all:,}")
print(f"Total CRE-gene links: {len(splicing_cre_links_all):,}")
print()
print(f"DATA SOURCE: Literature CRE-gene correlations (Table 16)")
print(f"FILTERS: FDR < {FDR_THRESHOLD}, |PCC| > {PCC_THRESHOLD}")
print()
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
