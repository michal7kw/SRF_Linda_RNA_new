#!/usr/bin/env python3
"""
Convert ENCODE cCRE TSV Files to BED Format for Splicing Genes

This script converts the ENCODE cCRE TSV files to BED format for use with deepTools.
Creates separate BED files for:
1. All CREs combined
2. By CRE type (dELS, pELS, CTCF-only, etc.)

INPUT:
- output/splicing_encode_cCREs_all.tsv
- output/splicing_encode_cCREs_by_type.tsv

OUTPUT:
- output/splicing_encode_cCREs_all.bed (all CREs, BED6 format)
- output/splicing_encode_cCREs_GABA.bed (same, for GABA analysis)
- output/splicing_encode_cCREs_{type}.bed (type-specific BED files)

Usage:
    python 2_convert_to_bed.py
"""

import pandas as pd
import os
from datetime import datetime

os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/splicing_encode_cCREs")

print("="*80)
print("CONVERT ENCODE cCRE TSV FILES TO BED FORMAT")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# =============================================================================
# Configuration
# =============================================================================

OUTPUT_DIR = "./output"

# =============================================================================
# Helper Function
# =============================================================================

def convert_tsv_to_bed(df, bed_file, description, name_col='cCRE_id1'):
    """
    Convert TSV dataframe to BED format.

    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe with chr, start, end columns
    bed_file : str
        Output BED file path
    description : str
        Description for logging
    name_col : str
        Column to use for BED name field (default: 'cCRE_id1')
    """
    print(f"\n{'='*80}")
    print(f"Converting: {description}")
    print(f"{'='*80}")
    print(f"Output: {bed_file}")

    # Extract unique CREs
    unique_cres = df[['chr', 'start', 'end', name_col]].drop_duplicates()

    # Sort by chromosome and start position
    chr_order = [f'chr{i}' for i in range(1, 20)] + ['chrX', 'chrY', 'chrM']
    unique_cres['chr'] = pd.Categorical(unique_cres['chr'], categories=chr_order, ordered=True)
    unique_cres = unique_cres.sort_values(['chr', 'start'])

    # Add score and strand columns for BED6 format
    unique_cres['score'] = 1000
    unique_cres['strand'] = '.'

    # Reorder columns: chr, start, end, name, score, strand
    bed_df = unique_cres[['chr', 'start', 'end', name_col, 'score', 'strand']]

    # Save to BED file
    bed_df.to_csv(bed_file, sep='\t', header=False, index=False)

    print(f"Unique CREs: {len(bed_df):,}")
    print(f"\nFirst 5 entries:")
    print(bed_df.head().to_string(index=False, header=False))

    # Show genes if available in original df
    if 'Gene' in df.columns:
        genes = df['Gene'].unique()
        print(f"\nGenes represented: {len(genes)}")
        print(f"Sample genes: {', '.join(sorted(genes)[:10])}")

    return bed_df

# =============================================================================
# Step 1: Convert All CREs
# =============================================================================

print("="*80)
print("STEP 1: Converting all CREs")
print("-"*80)

tsv_all = os.path.join(OUTPUT_DIR, "splicing_encode_cCREs_all.tsv")
bed_all = os.path.join(OUTPUT_DIR, "splicing_encode_cCREs_all.bed")
bed_gaba = os.path.join(OUTPUT_DIR, "splicing_encode_cCREs_GABA.bed")

if os.path.exists(tsv_all):
    df_all = pd.read_csv(tsv_all, sep='\t')
    print(f"Loaded {len(df_all):,} CRE-gene links")
    bed_all_df = convert_tsv_to_bed(df_all, bed_all, "All CREs")

    # Copy to GABA file (same CREs, will be analyzed with GABA BigWigs)
    bed_all_df.to_csv(bed_gaba, sep='\t', header=False, index=False)
    print(f"\nAlso saved: {bed_gaba}")
else:
    print(f"WARNING: File not found: {tsv_all}")
    bed_all_df = None

# =============================================================================
# Step 2: Convert by CRE Type
# =============================================================================

print("\n" + "="*80)
print("STEP 2: Converting by CRE type")
print("-"*80)

tsv_by_type = os.path.join(OUTPUT_DIR, "splicing_encode_cCREs_by_type.tsv")

type_bed_files = {}

if os.path.exists(tsv_by_type):
    df_by_type = pd.read_csv(tsv_by_type, sep='\t')
    print(f"Loaded {len(df_by_type):,} CRE-gene links")

    # Get unique CRE types
    cre_types = df_by_type['cre_type'].unique()
    print(f"\nFound {len(cre_types)} CRE types:")
    for cre_type in sorted(cre_types):
        count = len(df_by_type[df_by_type['cre_type'] == cre_type])
        print(f"  {cre_type}: {count:,} links")

    # Create BED file for each type
    print(f"\nCreating type-specific BED files:")
    for cre_type in sorted(cre_types):
        # Create safe filename (replace special characters)
        safe_type = cre_type.replace(',', '_').replace(' ', '_').replace('-', '_')
        bed_file = os.path.join(OUTPUT_DIR, f"splicing_encode_cCREs_{safe_type}.bed")

        # Filter for this type
        df_type = df_by_type[df_by_type['cre_type'] == cre_type]

        # Convert to BED
        type_bed_df = convert_tsv_to_bed(
            df_type,
            bed_file,
            f"CRE type: {cre_type}"
        )

        type_bed_files[cre_type] = bed_file
else:
    print(f"WARNING: File not found: {tsv_by_type}")

# =============================================================================
# Step 3: Summary
# =============================================================================

print("\n" + "="*80)
print("CONVERSION COMPLETE")
print("="*80)

print("\nGenerated BED files:")
print("\nMain files:")
if bed_all_df is not None:
    print(f"  {bed_all} ({len(bed_all_df):,} unique CREs)")
    print(f"  {bed_gaba} ({len(bed_all_df):,} unique CREs)")

if type_bed_files:
    print("\nType-specific files:")
    for cre_type, bed_file in sorted(type_bed_files.items()):
        # Count lines in file
        with open(bed_file) as f:
            n_cres = sum(1 for _ in f)
        print(f"  {bed_file} ({n_cres:,} CREs)")

print("\nThese BED files are ready for use with deepTools:")
print("  - computeMatrix")
print("  - plotHeatmap")
print("  - plotProfile")

print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
