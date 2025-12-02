#!/usr/bin/env python3
"""
Extract ENCODE cCREs Associated with GABA DEGs

This script identifies ENCODE cCREs (from mm10-cCREs.bed) that are associated with
GABA differentially expressed genes (up-regulated and down-regulated separately),
using literature CRE-gene correlations from Table 16.

WORKFLOW:
1. Load GABA DEG lists (up-regulated and down-regulated)
2. Query Table 16 for CRE-gene links for these DEGs
3. Apply statistical filters (FDR < 0.05, |PCC| > 0.2)
4. Filter for GABA/hippocampal cell types
5. Use bedtools intersect to find overlapping ENCODE cCREs
6. Create output files for up/down DEGs separately

INPUT:
- GABA DEG lists (up/down)
- Table 16 CRE-gene correlations
- ENCODE cCREs (mm10-cCREs.bed)

OUTPUT:
- ENCODE cCRE TSV and BED files for up-regulated DEGs
- ENCODE cCRE TSV and BED files for down-regulated DEGs
- Summary statistics

Usage:
    python 1_extract_encode_cCREs_for_DEGs.py
"""

import pandas as pd
import numpy as np
import os
import subprocess
from datetime import datetime
from collections import defaultdict

os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_GABA_DEGs_encode")

print("="*80)
print("EXTRACT ENCODE cCREs ASSOCIATED WITH GABA DEGs")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("", flush=True)

# =============================================================================
# Configuration
# =============================================================================

# Input files - GABA DEGs
DEG_UP_FILE = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/biomarkers/sig_deg_lists/GABA/Cluster_GABA_vs_Rest_up_significant.csv"
DEG_DOWN_FILE = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/biomarkers/sig_deg_lists/GABA/Cluster_GABA_vs_Rest_down_significant.csv"
TABLE16_FILE = "../data/table_16.txt"
ENCODE_CCRES_FILE = "../data/mm10-cCREs.bed"

# Statistical filters
FDR_THRESHOLD = 0.05
PCC_THRESHOLD = 0.2

# GABA/hippocampal cell types
GABA_SUBTYPES = {
    'Hippocampus', 'CA1', 'CA2', 'CA3', 'CA4', 'DG', 'DGNBL', 'GRC',
    'GABA', 'GABAergic', 'Interneuron', 'PV', 'SST', 'VIP',
    'LAMP5', 'LAMP', 'PVGA', 'SSTGA', 'VIPGA', 'LAMGA', 'INH'
}

# Output directory
OUTPUT_DIR = "./output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# Step 1: Load GABA DEG Lists
# =============================================================================

print("="*80)
print("STEP 1: Loading GABA DEG lists")
print("-"*80, flush=True)

# Load up-regulated DEGs
deg_up = pd.read_csv(DEG_UP_FILE)
deg_up_genes = set(deg_up['names'].dropna())
deg_up_genes_upper = {g.upper() for g in deg_up_genes}
print(f"Up-regulated DEGs: {len(deg_up_genes)} genes")

# Load down-regulated DEGs
deg_down = pd.read_csv(DEG_DOWN_FILE)
deg_down_genes = set(deg_down['names'].dropna())
deg_down_genes_upper = {g.upper() for g in deg_down_genes}
print(f"Down-regulated DEGs: {len(deg_down_genes)} genes")

# Combined for initial Table 16 filtering
all_deg_genes_upper = deg_up_genes_upper | deg_down_genes_upper
print(f"\nTotal unique DEGs: {len(all_deg_genes_upper)}")
print("", flush=True)

# =============================================================================
# Step 2: Load and Filter Table 16 (CHUNKED)
# =============================================================================

print("="*80)
print("STEP 2: Loading and filtering Table 16 for DEGs (chunked processing)")
print("-"*80)
print(f"Loading: {TABLE16_FILE}")
print("Processing in chunks to manage memory...", flush=True)

def is_gaba_subtype(subtype):
    """Check if SubType is GABA/hippocampal"""
    if pd.isna(subtype):
        return False
    subtype_str = str(subtype).upper()
    return any(gaba_type.upper() in subtype_str for gaba_type in GABA_SUBTYPES)

chunk_size = 500000
filtered_chunks = []
total_rows = 0
filtered_rows = 0

for i, chunk in enumerate(pd.read_csv(TABLE16_FILE, sep='\t', chunksize=chunk_size)):
    total_rows += len(chunk)

    # Filter for DEGs
    chunk['Gene_upper'] = chunk['Gene'].str.upper()
    deg_chunk = chunk[chunk['Gene_upper'].isin(all_deg_genes_upper)]

    # Apply statistical filters
    deg_chunk = deg_chunk[
        (deg_chunk['FDR'] < FDR_THRESHOLD) &
        (abs(deg_chunk['PCC']) > PCC_THRESHOLD)
    ]

    # Apply GABA cell type filter
    deg_chunk = deg_chunk[deg_chunk['SubType'].apply(is_gaba_subtype)]

    if len(deg_chunk) > 0:
        filtered_chunks.append(deg_chunk)
        filtered_rows += len(deg_chunk)

    print(f"  Chunk {i+1}: processed {total_rows:,} rows, kept {filtered_rows:,} DEG links", flush=True)

if len(filtered_chunks) == 0:
    print("\nERROR: No DEG links found after filtering!")
    exit(1)

deg_links = pd.concat(filtered_chunks, ignore_index=True)
print(f"\nTotal Table 16 rows processed: {total_rows:,}")
print(f"DEG links (filtered, GABA only): {len(deg_links):,}")

# Split into up and down
deg_links['is_up'] = deg_links['Gene_upper'].isin(deg_up_genes_upper)
deg_links['is_down'] = deg_links['Gene_upper'].isin(deg_down_genes_upper)

deg_links_up = deg_links[deg_links['is_up']].copy()
deg_links_down = deg_links[deg_links['is_down']].copy()

print(f"  Up-DEG links: {len(deg_links_up):,}")
print(f"  Down-DEG links: {len(deg_links_down):,}")
print("", flush=True)

# =============================================================================
# Step 3: Parse Coordinates
# =============================================================================

print("="*80)
print("STEP 3: Parsing coordinates from Table 16")
print("-"*80, flush=True)

def parse_coordinate(coord_str):
    """Parse coordinate string (chr11_40754370_40756166) to chr, start, end"""
    try:
        parts = coord_str.split('_')
        if len(parts) == 3:
            return parts[0], int(parts[1]), int(parts[2])
    except:
        pass
    return None, None, None

for df, name in [(deg_links_up, 'Up-DEGs'), (deg_links_down, 'Down-DEGs')]:
    coords = df['Coordinate1'].apply(lambda x: pd.Series(parse_coordinate(x)))
    df[['coord_chr', 'coord_start', 'coord_end']] = coords
    df.dropna(subset=['coord_chr', 'coord_start', 'coord_end'], inplace=True)
    df['coord_start'] = df['coord_start'].astype(int)
    df['coord_end'] = df['coord_end'].astype(int)
    print(f"  {name}: {len(df):,} genomic regions parsed")

print("", flush=True)

# =============================================================================
# Step 4: Find Overlapping ENCODE cCREs using bedtools
# =============================================================================

print("="*80)
print("STEP 4: Finding overlapping ENCODE cCREs (using bedtools)")
print("-"*80, flush=True)

def find_encode_overlaps(df, name):
    """Find ENCODE cCREs overlapping with Table 16 regions"""
    print(f"\nProcessing {name}...")

    if len(df) == 0:
        print(f"  WARNING: No regions for {name}")
        return pd.DataFrame()

    # Reset index for merging
    df = df.reset_index(drop=True)

    # Create temporary BED file
    temp_bed = os.path.join(OUTPUT_DIR, f"temp_{name.lower().replace('-', '_')}_regions.bed")
    with open(temp_bed, 'w') as f:
        for idx, row in df.iterrows():
            f.write(f"{row['coord_chr']}\t{int(row['coord_start'])}\t{int(row['coord_end'])}\t{idx}\n")

    print(f"  Created temp BED: {len(df):,} regions")

    # Run bedtools intersect
    cmd = [
        "bedtools", "intersect",
        "-a", ENCODE_CCRES_FILE,
        "-b", temp_bed,
        "-wa", "-wb"
    ]

    print(f"  Running bedtools intersect...", flush=True)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=1800)

        overlaps = []
        lines = result.stdout.strip().split('\n')

        for line in lines:
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 10:
                try:
                    overlaps.append({
                        'chr': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'cCRE_id1': parts[3],
                        'cCRE_id2': parts[4],
                        'cre_type': parts[5],
                        'table16_idx': int(parts[9])
                    })
                except (ValueError, IndexError):
                    continue

        print(f"  Found {len(overlaps):,} overlaps")

        # Clean up temp file
        os.unlink(temp_bed)

        if len(overlaps) == 0:
            return pd.DataFrame()

        # Create dataframe and merge with metadata
        overlaps_df = pd.DataFrame(overlaps)
        merged = overlaps_df.merge(
            df[['Gene', 'SubType', 'PCC', 'FDR', 'Coordinate1', 'cCRE1']],
            left_on='table16_idx',
            right_index=True,
            how='left'
        )
        merged = merged.drop(columns=['table16_idx'])

        print(f"  Unique ENCODE cCREs: {merged['cCRE_id1'].nunique():,}")
        print(f"  Unique genes: {merged['Gene'].nunique()}")

        return merged

    except subprocess.CalledProcessError as e:
        print(f"ERROR running bedtools: {e.stderr}")
        os.unlink(temp_bed)
        return pd.DataFrame()
    except FileNotFoundError:
        print("ERROR: bedtools not found. Please ensure bedtools is in PATH.")
        os.unlink(temp_bed)
        return pd.DataFrame()

# Find overlaps for up and down DEGs
encode_up = find_encode_overlaps(deg_links_up, 'Up-DEGs')
encode_down = find_encode_overlaps(deg_links_down, 'Down-DEGs')

print("", flush=True)

# =============================================================================
# Step 5: Create BED files
# =============================================================================

print("="*80)
print("STEP 5: Creating output files")
print("-"*80, flush=True)

def save_outputs(df, prefix, description):
    """Save TSV and BED files"""
    if len(df) == 0:
        print(f"  WARNING: No data for {description}")
        return 0

    # Save TSV
    tsv_file = os.path.join(OUTPUT_DIR, f"encode_cCREs_{prefix}.tsv")
    df.to_csv(tsv_file, sep='\t', index=False)
    print(f"  Saved TSV: {tsv_file} ({len(df):,} links)")

    # Create BED file (unique CREs)
    unique_cres = df[['chr', 'start', 'end', 'cCRE_id1', 'cre_type']].drop_duplicates()

    # Sort by chromosome and position
    chr_order = [f'chr{i}' for i in range(1, 20)] + ['chrX', 'chrY', 'chrM']
    unique_cres['chr'] = pd.Categorical(unique_cres['chr'], categories=chr_order, ordered=True)
    unique_cres = unique_cres.sort_values(['chr', 'start'])

    # BED6 format
    unique_cres['score'] = 1000
    unique_cres['strand'] = '.'
    bed_df = unique_cres[['chr', 'start', 'end', 'cCRE_id1', 'score', 'strand']]

    bed_file = os.path.join(OUTPUT_DIR, f"encode_cCREs_{prefix}.bed")
    bed_df.to_csv(bed_file, sep='\t', header=False, index=False)
    print(f"  Saved BED: {bed_file} ({len(bed_df):,} unique CREs)")

    # Show CRE type distribution
    print(f"  CRE type distribution for {description}:")
    type_counts = unique_cres['cre_type'].value_counts()
    for cre_type, count in type_counts.items():
        print(f"    {cre_type}: {count}")

    return len(bed_df)

n_up = save_outputs(encode_up, 'up_DEGs', 'Up-regulated DEGs')
print()
n_down = save_outputs(encode_down, 'down_DEGs', 'Down-regulated DEGs')

print("", flush=True)

# =============================================================================
# Step 6: Generate Summary Report
# =============================================================================

print("="*80)
print("STEP 6: Generating summary report")
print("-"*80, flush=True)

summary_file = os.path.join(OUTPUT_DIR, "SUMMARY_GABA_DEGs_encode_cCREs.txt")

with open(summary_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("ENCODE cCREs ASSOCIATED WITH GABA DEGs - SUMMARY\n")
    f.write("="*80 + "\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("INPUT DATA:\n")
    f.write("-"*80 + "\n")
    f.write(f"Up-regulated DEGs: {len(deg_up_genes)} genes\n")
    f.write(f"Down-regulated DEGs: {len(deg_down_genes)} genes\n")
    f.write(f"ENCODE cCREs database: mm10-cCREs.bed\n\n")

    f.write("FILTERS:\n")
    f.write(f"  FDR < {FDR_THRESHOLD}\n")
    f.write(f"  |PCC| > {PCC_THRESHOLD}\n")
    f.write(f"  Cell types: GABA/Hippocampal\n\n")

    f.write("RESULTS:\n")
    f.write("="*80 + "\n\n")

    f.write("UP-REGULATED DEGs:\n")
    f.write("-"*40 + "\n")
    if len(encode_up) > 0:
        f.write(f"  Unique ENCODE cCREs: {encode_up['cCRE_id1'].nunique():,}\n")
        f.write(f"  Unique genes with cCREs: {encode_up['Gene'].nunique()}\n")
        f.write(f"  Total cCRE-gene links: {len(encode_up):,}\n")
        f.write("\n  CRE Type Distribution:\n")
        type_counts = encode_up.drop_duplicates(subset=['cCRE_id1'])['cre_type'].value_counts()
        for cre_type, count in type_counts.items():
            f.write(f"    {cre_type}: {count}\n")
    else:
        f.write("  No overlapping cCREs found\n")

    f.write("\nDOWN-REGULATED DEGs:\n")
    f.write("-"*40 + "\n")
    if len(encode_down) > 0:
        f.write(f"  Unique ENCODE cCREs: {encode_down['cCRE_id1'].nunique():,}\n")
        f.write(f"  Unique genes with cCREs: {encode_down['Gene'].nunique()}\n")
        f.write(f"  Total cCRE-gene links: {len(encode_down):,}\n")
        f.write("\n  CRE Type Distribution:\n")
        type_counts = encode_down.drop_duplicates(subset=['cCRE_id1'])['cre_type'].value_counts()
        for cre_type, count in type_counts.items():
            f.write(f"    {cre_type}: {count}\n")
    else:
        f.write("  No overlapping cCREs found\n")

    f.write("\nOUTPUT FILES:\n")
    f.write("-"*80 + "\n")
    f.write("  encode_cCREs_up_DEGs.tsv - All links for up-regulated DEGs\n")
    f.write("  encode_cCREs_up_DEGs.bed - BED format (unique CREs)\n")
    f.write("  encode_cCREs_down_DEGs.tsv - All links for down-regulated DEGs\n")
    f.write("  encode_cCREs_down_DEGs.bed - BED format (unique CREs)\n")

print(f"Saved summary: {summary_file}")
print("", flush=True)

# =============================================================================
# Completion
# =============================================================================

print("="*80)
print("ANALYSIS COMPLETE")
print("="*80)
print(f"\nKey Results:")
print(f"  Up-regulated DEGs: {n_up} unique ENCODE cCREs")
print(f"  Down-regulated DEGs: {n_down} unique ENCODE cCREs")

if len(encode_up) > 0:
    print(f"\nTop up-DEGs by cCRE count:")
    top_up = encode_up.groupby('Gene').size().sort_values(ascending=False).head(5)
    for gene, count in top_up.items():
        print(f"  {gene}: {count} cCREs")

if len(encode_down) > 0:
    print(f"\nTop down-DEGs by cCRE count:")
    top_down = encode_down.groupby('Gene').size().sort_values(ascending=False).head(5)
    for gene, count in top_down.items():
        print(f"  {gene}: {count} cCREs")

print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
