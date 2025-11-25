#!/usr/bin/env python3
"""
Extract ENCODE cCREs Associated with Splicing Genes (OPTIMIZED VERSION)

This script identifies ENCODE cCREs (from mm10-cCREs.bed) that are associated with
splicing-related genes, using literature CRE-gene correlations from Table 16.

OPTIMIZATIONS:
- Uses bedtools intersect for fast overlap detection (vs naive O(n*m) loops)
- Chunked processing of large Table 16 file
- Memory-efficient data structures
- Progress tracking with flush

Usage:
    python 1_extract_encode_cCREs.py
"""

import pandas as pd
import numpy as np
import os
import subprocess
import tempfile
from datetime import datetime
from collections import defaultdict

os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/encode_cCREs")

print("="*80)
print("EXTRACT ENCODE cCREs ASSOCIATED WITH SPLICING GENES (OPTIMIZED)")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("", flush=True)

# =============================================================================
# Configuration
# =============================================================================

# Input files
SPLICING_GENES_FILE = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/splicing_genes/extracted_genes_final.csv"
TABLE16_FILE = "../data/table_16.txt"  # Literature CRE-gene correlations
ENCODE_CCRES_FILE = "../data/mm10-cCREs.bed"  # ENCODE cCREs

# Statistical filters
FDR_THRESHOLD = 0.05
PCC_THRESHOLD = 0.2

# Output directory
OUTPUT_DIR = "./output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Output files
OUTPUT_ALL = os.path.join(OUTPUT_DIR, "encode_cCREs_all_celltypes.tsv")
OUTPUT_GABA = os.path.join(OUTPUT_DIR, "encode_cCREs_GABA.tsv")
OUTPUT_BY_TYPE = os.path.join(OUTPUT_DIR, "encode_cCREs_by_type.tsv")
OUTPUT_SUMMARY = os.path.join(OUTPUT_DIR, "SUMMARY_encode_cCREs.txt")

# =============================================================================
# Step 1: Load Splicing Genes List
# =============================================================================

print("="*80)
print("STEP 1: Loading splicing genes list")
print("-"*80, flush=True)

splicing_genes = pd.read_csv(SPLICING_GENES_FILE)
print(f"Loaded {len(splicing_genes)} splicing genes from Reactome/GO")

# Extract unique gene symbols
splicing_gene_symbols = set(splicing_genes['gene_symbol'].dropna())
splicing_gene_symbols = {g for g in splicing_gene_symbols if not g.startswith('[')}
splicing_gene_symbols_upper = {g.upper() for g in splicing_gene_symbols}

print(f"Unique gene symbols: {len(splicing_gene_symbols)}")
print(f"Sample genes: {', '.join(list(splicing_gene_symbols)[:10])}")
print("", flush=True)

# =============================================================================
# Step 2: Load and Filter Table 16 (CHUNKED for memory efficiency)
# =============================================================================

print("="*80)
print("STEP 2: Loading and filtering Table 16 (chunked processing)")
print("-"*80)
print(f"Loading: {TABLE16_FILE}")
print("Processing in chunks to manage memory...", flush=True)

# Read in chunks and filter immediately to reduce memory
chunk_size = 500000
filtered_chunks = []
total_rows = 0
filtered_rows = 0

for i, chunk in enumerate(pd.read_csv(TABLE16_FILE, sep='\t', chunksize=chunk_size)):
    total_rows += len(chunk)

    # Filter for splicing genes
    chunk['Gene_upper'] = chunk['Gene'].str.upper()
    splicing_chunk = chunk[chunk['Gene_upper'].isin(splicing_gene_symbols_upper)]

    # Apply statistical filters
    splicing_chunk = splicing_chunk[
        (splicing_chunk['FDR'] < FDR_THRESHOLD) &
        (abs(splicing_chunk['PCC']) > PCC_THRESHOLD)
    ]

    if len(splicing_chunk) > 0:
        filtered_chunks.append(splicing_chunk)
        filtered_rows += len(splicing_chunk)

    print(f"  Chunk {i+1}: processed {total_rows:,} rows, kept {filtered_rows:,} splicing gene links", flush=True)

# Combine filtered chunks
if len(filtered_chunks) == 0:
    print("\nERROR: No splicing gene links found after filtering!")
    print("Check that Table 16 contains the expected genes.")
    exit(1)

splicing_links = pd.concat(filtered_chunks, ignore_index=True)
print(f"\nTotal rows processed: {total_rows:,}")
print(f"Splicing gene links (filtered): {len(splicing_links):,}")
print("", flush=True)

# =============================================================================
# Step 3: Parse Coordinates and Create BED file for Table 16 regions
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

# Parse coordinates
coords = splicing_links['Coordinate1'].apply(lambda x: pd.Series(parse_coordinate(x)))
splicing_links[['coord_chr', 'coord_start', 'coord_end']] = coords

# Remove invalid coordinates
splicing_links = splicing_links.dropna(subset=['coord_chr', 'coord_start', 'coord_end'])
splicing_links['coord_start'] = splicing_links['coord_start'].astype(int)
splicing_links['coord_end'] = splicing_links['coord_end'].astype(int)

print(f"Parsed {len(splicing_links):,} genomic regions")
print("", flush=True)

# =============================================================================
# Step 4: Use bedtools intersect for FAST overlap detection
# =============================================================================

print("="*80)
print("STEP 4: Finding overlapping ENCODE cCREs (using bedtools)")
print("-"*80, flush=True)

# Reset index and save it for later merging
splicing_links = splicing_links.reset_index(drop=True)

# Create temporary BED file for Table 16 regions
table16_bed = os.path.join(OUTPUT_DIR, "temp_table16_regions.bed")
with open(table16_bed, 'w') as f:
    for idx, row in splicing_links.iterrows():
        # BED format: chr, start, end, name (use index as identifier)
        f.write(f"{row['coord_chr']}\t{int(row['coord_start'])}\t{int(row['coord_end'])}\t{idx}\n")

print(f"Created temporary BED file: {table16_bed}")
print(f"  Contains {len(splicing_links):,} regions")

# Run bedtools intersect
print("\nRunning bedtools intersect...", flush=True)

cmd = [
    "bedtools", "intersect",
    "-a", ENCODE_CCRES_FILE,  # ENCODE cCREs
    "-b", table16_bed,         # Table 16 regions
    "-wa", "-wb"               # Report both A and B features
]

print(f"Command: {' '.join(cmd)}", flush=True)

try:
    result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=1800)  # 30 min timeout

    # Parse intersect output
    overlaps = []
    lines = result.stdout.strip().split('\n')
    print(f"\nParsing {len(lines):,} intersection results...", flush=True)

    for line in lines:
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) >= 10:
            # Parts 0-5: ENCODE cCRE (chr, start, end, id1, id2, type)
            # Parts 6-9: Table 16 region (chr, start, end, index)
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
            except (ValueError, IndexError) as e:
                continue  # Skip malformed lines

    print(f"Found {len(overlaps):,} overlaps")

except subprocess.TimeoutExpired:
    print("ERROR: bedtools intersect timed out after 30 minutes")
    raise
except subprocess.CalledProcessError as e:
    print(f"ERROR running bedtools: {e.stderr}")
    raise
except FileNotFoundError:
    print("ERROR: bedtools not found. Please ensure bedtools is in PATH.")
    print("Try: module load bedtools  or  conda activate <env with bedtools>")
    raise

# Clean up temp file
os.unlink(table16_bed)
print("Cleaned up temporary files", flush=True)
print("", flush=True)

# =============================================================================
# Step 5: Merge overlap data with Table 16 metadata
# =============================================================================

print("="*80)
print("STEP 5: Merging overlap data with Table 16 metadata")
print("-"*80, flush=True)

# Create dataframe from overlaps
overlaps_df = pd.DataFrame(overlaps)

if len(overlaps_df) == 0:
    print("WARNING: No overlapping cCREs found!")
    print("This may indicate:")
    print("  - No ENCODE cCREs overlap with splicing gene regions")
    print("  - Check that chromosome names match (chr1 vs 1)")
    print("  - Consider relaxing statistical thresholds")

    # Create empty output files
    empty_df = pd.DataFrame()
    empty_df.to_csv(OUTPUT_ALL, sep='\t', index=False)
    empty_df.to_csv(OUTPUT_GABA, sep='\t', index=False)
    empty_df.to_csv(OUTPUT_BY_TYPE, sep='\t', index=False)

    with open(OUTPUT_SUMMARY, 'w') as f:
        f.write("No overlapping cCREs found.\n")

    print("\nCreated empty output files.")
    exit(0)

# Merge with Table 16 metadata
merged = overlaps_df.merge(
    splicing_links[['Gene', 'SubType', 'PCC', 'FDR', 'Coordinate1', 'cCRE1']],
    left_on='table16_idx',
    right_index=True,
    how='left'
)

# Drop the temporary index column
merged = merged.drop(columns=['table16_idx'])

print(f"Merged {len(merged):,} cCRE-gene associations")
print(f"Unique ENCODE cCREs: {merged['cCRE_id1'].nunique():,}")
print(f"Unique genes: {merged['Gene'].nunique()}")
print("", flush=True)

# =============================================================================
# Step 6: Classify by Cell Type
# =============================================================================

print("="*80)
print("STEP 6: Classifying by cell type")
print("-"*80, flush=True)

# Define GABA/hippocampal cell types
GABA_SUBTYPES = {
    'Hippocampus', 'CA1', 'CA2', 'CA3', 'CA4', 'DG',
    'GABA', 'GABAergic', 'Interneuron', 'PV', 'SST', 'VIP'
}

def is_gaba_subtype(subtype):
    """Check if SubType is GABA/hippocampal"""
    if pd.isna(subtype):
        return False
    subtype_str = str(subtype).upper()
    return any(gaba_type.upper() in subtype_str for gaba_type in GABA_SUBTYPES)

merged['is_gaba'] = merged['SubType'].apply(is_gaba_subtype)
gaba_ccres = merged[merged['is_gaba']].copy()

print(f"Total cCRE-gene links: {len(merged):,}")
print(f"GABA/hippocampal links: {len(gaba_ccres):,}")
print("", flush=True)

# =============================================================================
# Step 7: CRE Type Analysis
# =============================================================================

print("="*80)
print("STEP 7: Analyzing CRE types")
print("-"*80, flush=True)

print("CRE type distribution (all cell types):")
type_counts_all = merged['cre_type'].value_counts()
for cre_type, count in type_counts_all.items():
    print(f"  {cre_type}: {count:,}")

print("\nCRE type distribution (GABA cell types):")
if len(gaba_ccres) > 0:
    type_counts_gaba = gaba_ccres['cre_type'].value_counts()
    for cre_type, count in type_counts_gaba.items():
        print(f"  {cre_type}: {count:,}")
else:
    print("  No GABA-specific cCREs found")
    type_counts_gaba = pd.Series(dtype=int)
print("", flush=True)

# =============================================================================
# Step 8: Save Output Files
# =============================================================================

print("="*80)
print("STEP 8: Saving output files")
print("-"*80, flush=True)

# Save all cell types
merged.to_csv(OUTPUT_ALL, sep='\t', index=False)
print(f"Saved: {OUTPUT_ALL}")
print(f"  {len(merged):,} CRE-gene links")

# Save GABA cell types
gaba_ccres.to_csv(OUTPUT_GABA, sep='\t', index=False)
print(f"Saved: {OUTPUT_GABA}")
print(f"  {len(gaba_ccres):,} CRE-gene links")

# Save by type (same as all, includes cre_type column)
merged.to_csv(OUTPUT_BY_TYPE, sep='\t', index=False)
print(f"Saved: {OUTPUT_BY_TYPE}")
print("", flush=True)

# =============================================================================
# Step 9: Generate Summary Report
# =============================================================================

print("="*80)
print("STEP 9: Generating summary report")
print("-"*80, flush=True)

with open(OUTPUT_SUMMARY, 'w') as f:
    f.write("="*80 + "\n")
    f.write("ENCODE cCREs ASSOCIATED WITH SPLICING GENES - SUMMARY\n")
    f.write("="*80 + "\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("INPUT DATA:\n")
    f.write(f"  Splicing genes: {len(splicing_gene_symbols)}\n")
    f.write(f"  Table 16 links (filtered): {len(splicing_links):,}\n\n")

    f.write("OVERLAPPING cCREs:\n")
    f.write(f"  Total cCRE-gene links: {len(merged):,}\n")
    f.write(f"  Unique cCREs: {merged['cCRE_id1'].nunique():,}\n")
    f.write(f"  Unique genes: {merged['Gene'].nunique()}\n\n")

    f.write("CELL TYPE CLASSIFICATION:\n")
    f.write(f"  GABA/hippocampal links: {len(gaba_ccres):,}\n")
    if len(gaba_ccres) > 0:
        f.write(f"  GABA unique cCREs: {gaba_ccres['cCRE_id1'].nunique():,}\n")
        f.write(f"  GABA unique genes: {gaba_ccres['Gene'].nunique()}\n\n")

    f.write("CRE TYPE DISTRIBUTION (ALL):\n")
    for cre_type, count in type_counts_all.items():
        f.write(f"  {cre_type}: {count:,}\n")

    if len(type_counts_gaba) > 0:
        f.write("\nCRE TYPE DISTRIBUTION (GABA):\n")
        for cre_type, count in type_counts_gaba.items():
            f.write(f"  {cre_type}: {count:,}\n")

    f.write("\nOUTPUT FILES:\n")
    f.write(f"  {OUTPUT_ALL}\n")
    f.write(f"  {OUTPUT_GABA}\n")
    f.write(f"  {OUTPUT_BY_TYPE}\n")

print(f"Saved summary: {OUTPUT_SUMMARY}")
print("", flush=True)

# =============================================================================
# Step 10: Display Results
# =============================================================================

print("="*80)
print("ANALYSIS COMPLETE")
print("="*80)
print(f"\nKey Results:")
print(f"  ENCODE cCREs linked to splicing genes: {merged['cCRE_id1'].nunique():,}")
print(f"  Splicing genes represented: {merged['Gene'].nunique()}")
print(f"  Total CRE-gene associations: {len(merged):,}")
print(f"\nGABA/Hippocampal subset:")
if len(gaba_ccres) > 0:
    print(f"  Unique cCREs: {gaba_ccres['cCRE_id1'].nunique():,}")
    print(f"  Genes: {gaba_ccres['Gene'].nunique()}")
    print(f"  Associations: {len(gaba_ccres):,}")
else:
    print("  No GABA-specific cCREs found")

print(f"\nTop genes by number of associated cCREs:")
top_genes = merged.groupby('Gene').size().sort_values(ascending=False).head(10)
for gene, count in top_genes.items():
    print(f"  {gene}: {count} cCREs")

print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
