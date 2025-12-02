#!/usr/bin/env python3
"""
Convert GABA-Specific CRE TSV Files to BED Format

This script converts the GABA-specific CRE TSV files to BED format
for use with deepTools.

Creates TWO sets of BED files:

1. TABLE 16-ONLY CREs (no ENCODE intersection required):
   - GABA_specific_table16_cCREs.bed - All Table 16 GABA-specific regions

2. ENCODE-INTERSECTED CREs (Table 16 + mm10-cCREs.bed):
   - GABA_specific_encode_cCREs.bed - Intersected CREs
   - GABA_specific_encode_cCREs_{type}.bed - Type-specific BED files

INPUT:
- output/GABA_specific_table16_cCREs.tsv (Table 16-only)
- output/GABA_specific_encode_cCREs.tsv (intersected)
- output/GABA_specific_encode_cCREs_by_type.tsv (intersected by type)

OUTPUT:
- output/GABA_specific_table16_cCREs.bed (Table 16-only, BED6 format)
- output/GABA_specific_encode_cCREs.bed (intersected, BED6 format)
- output/GABA_specific_encode_cCREs_{type}.bed (type-specific BED files)
"""

import pandas as pd
import os
from datetime import datetime

# Change to script directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)

print("="*80)
print("CONVERT GABA-SPECIFIC ENCODE cCRE TSV FILES TO BED FORMAT")
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
# Step 1: Convert Table 16-Only CREs (no ENCODE intersection)
# =============================================================================

print("="*80)
print("STEP 1: Converting Table 16-only CREs to BED format")
print("-"*80)

tsv_table16 = os.path.join(OUTPUT_DIR, "GABA_specific_table16_cCREs.tsv")
bed_table16 = os.path.join(OUTPUT_DIR, "GABA_specific_table16_cCREs.bed")

bed_table16_df = None
if os.path.exists(tsv_table16):
    df_table16 = pd.read_csv(tsv_table16, sep='\t')
    print(f"Loaded {len(df_table16):,} CRE-gene links")

    # For Table 16-only, we need to create a unique ID from coordinates
    df_table16['cre_id'] = df_table16['chr'] + '_' + df_table16['start'].astype(str) + '_' + df_table16['end'].astype(str)

    # Convert to BED using coordinate-based ID
    bed_table16_df = convert_tsv_to_bed(df_table16, bed_table16, "Table 16-only GABA CREs", name_col='cre_id')
else:
    print(f"WARNING: File not found: {tsv_table16}")
    print("Table 16-only BED file will not be created.")

# =============================================================================
# Step 2: Convert ENCODE-Intersected CREs
# =============================================================================

print("\n" + "="*80)
print("STEP 2: Converting ENCODE-intersected CREs to BED format")
print("-"*80)

tsv_gaba = os.path.join(OUTPUT_DIR, "GABA_specific_encode_cCREs.tsv")
bed_gaba = os.path.join(OUTPUT_DIR, "GABA_specific_encode_cCREs.bed")

bed_gaba_df = None
if os.path.exists(tsv_gaba):
    df_gaba = pd.read_csv(tsv_gaba, sep='\t')
    if len(df_gaba) > 0:
        print(f"Loaded {len(df_gaba):,} CRE-gene links")
        bed_gaba_df = convert_tsv_to_bed(df_gaba, bed_gaba, "ENCODE-intersected GABA CREs")
    else:
        print("No ENCODE-intersected CREs found (empty file)")
else:
    print(f"WARNING: File not found: {tsv_gaba}")
    print("ENCODE-intersected BED file will not be created.")

# =============================================================================
# Step 3: Convert by CRE Type (ENCODE-intersected only)
# =============================================================================

print("\n" + "="*80)
print("STEP 3: Converting by CRE type (ENCODE-intersected only)")
print("-"*80)

tsv_by_type = os.path.join(OUTPUT_DIR, "GABA_specific_encode_cCREs_by_type.tsv")

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
        bed_file = os.path.join(OUTPUT_DIR, f"GABA_specific_encode_cCREs_{safe_type}.bed")

        # Filter for this type
        df_type = df_by_type[df_by_type['cre_type'] == cre_type]

        # Convert to BED
        type_bed_df = convert_tsv_to_bed(
            df_type,
            bed_file,
            f"GABA-specific CRE type: {cre_type}"
        )

        type_bed_files[cre_type] = bed_file
else:
    print(f"WARNING: File not found: {tsv_by_type}")

# =============================================================================
# Step 4: Summary
# =============================================================================

print("\n" + "="*80)
print("CONVERSION COMPLETE")
print("="*80)

print("\nGenerated BED files:")

print("\n1. TABLE 16-ONLY (no ENCODE intersection):")
if bed_table16_df is not None:
    print(f"   {bed_table16} ({len(bed_table16_df):,} unique CREs)")
else:
    print("   Not created (source file missing)")

print("\n2. ENCODE-INTERSECTED:")
if bed_gaba_df is not None:
    print(f"   {bed_gaba} ({len(bed_gaba_df):,} unique CREs)")
else:
    print("   Not created (no overlapping CREs or source file missing)")

if type_bed_files:
    print("\n3. TYPE-SPECIFIC (ENCODE-intersected):")
    for cre_type, bed_file in sorted(type_bed_files.items()):
        # Count lines in file
        with open(bed_file) as f:
            n_cres = sum(1 for _ in f)
        print(f"   {bed_file} ({n_cres:,} CREs)")

print("\nRecommendation:")
print("  - Use TABLE 16-ONLY for maximum coverage (all GABA CREs)")
print("  - Use ENCODE-INTERSECTED for high-confidence validated CREs")

print("\nThese BED files are ready for use with deepTools:")
print("  - computeMatrix")
print("  - plotHeatmap")
print("  - plotProfile")

print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
