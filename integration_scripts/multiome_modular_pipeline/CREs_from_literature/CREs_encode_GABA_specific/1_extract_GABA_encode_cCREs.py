#!/usr/bin/env python3
"""
Extract GABA Cell Type Specific ENCODE cCREs

This script identifies ENCODE cCREs (from mm10-cCREs.bed) that are specifically
associated with hippocampal GABAergic cell types, using literature CRE-gene
correlations from Table 16 with EXPLICIT cell type filtering.

This is an alternative to CREs_encode_paper_intersection that focuses EXCLUSIVELY
on GABA cell type specific CREs from the start, providing cleaner positive controls
for GABA ATAC-seq validation.

METHODOLOGY (based on data exploration in helpers/explore_data/):
1. Hippocampal samples identified from Table 1 (Major Region == 'HPF')
2. Cells from those samples filtered by Class == 'GABA' in Table 2
3. The 46 unique GABAergic CellTypes from hippocampus are used for filtering

This EXCLUDES glutamatergic cell types like CA1GL, CA3GL, DGGR that were
incorrectly matched by the previous keyword-based approach.

Included GABAergic cell types (46 total):
- PVGA1-7: Parvalbumin interneurons
- SSTGA1-10: Somatostatin interneurons
- VIPGA1-4: VIP interneurons
- LAMGA1-4: Lamp5 interneurons
- DGNBL1-2: Dentate gyrus neuroblasts
- Plus minor populations: MSGA, LSXGA, CNUGA, STRGA, OBNBL, D2MSN

Usage:
    python 1_extract_GABA_encode_cCREs.py
"""

import pandas as pd
import numpy as np
import os
import subprocess
import tempfile
from datetime import datetime
from collections import defaultdict

# Change to script directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)

print("="*80)
print("EXTRACT GABA CELL TYPE SPECIFIC ENCODE cCREs")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("", flush=True)

# =============================================================================
# Configuration
# =============================================================================

# Input files
TABLE16_FILE = "../data/table_16.txt"  # Literature CRE-gene correlations
ENCODE_CCRES_FILE = "../data/mm10-cCREs.bed"  # ENCODE cCREs

# Statistical filters
FDR_THRESHOLD = 0.05
PCC_THRESHOLD = 0.2

# Output directory
OUTPUT_DIR = "./output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Output files - Intersected (Table 16 + ENCODE)
OUTPUT_INTERSECTED = os.path.join(OUTPUT_DIR, "GABA_specific_encode_cCREs.tsv")
OUTPUT_INTERSECTED_BY_TYPE = os.path.join(OUTPUT_DIR, "GABA_specific_encode_cCREs_by_type.tsv")

# Output files - Table 16 only (no ENCODE intersection required)
OUTPUT_TABLE16_ONLY = os.path.join(OUTPUT_DIR, "GABA_specific_table16_cCREs.tsv")

# Summary
OUTPUT_SUMMARY = os.path.join(OUTPUT_DIR, "SUMMARY_GABA_specific_cCREs.txt")

# =============================================================================
# Hippocampal GABAergic Cell Types (from data exploration)
# =============================================================================

# Explicit list of 46 GABAergic cell types found in hippocampal samples
# Derived from: helpers/explore_data/identify_hippocampal_gaba_types.py
# Method: Table 1 (HPF samples) -> Table 2 (Class == 'GABA') -> CellType
#
# This EXCLUDES glutamatergic cell types that were incorrectly matched before:
# - CA1GL1-3, CA3GL1-6: Glutamatergic pyramidal neurons
# - DGGR: Dentate gyrus granule cells (glutamatergic)
# - PVM: Perivascular macrophages (non-neuronal)
# - RGDG: Radial glia

HIPPOCAMPAL_GABA_CELLTYPES = {
    # Parvalbumin interneurons (PV+)
    'PVGA1', 'PVGA2', 'PVGA3', 'PVGA4', 'PVGA5', 'PVGA6', 'PVGA7',
    # Somatostatin interneurons (SST+)
    'SSTGA1', 'SSTGA2', 'SSTGA3', 'SSTGA4', 'SSTGA5', 'SSTGA6',
    'SSTGA7', 'SSTGA8', 'SSTGA9', 'SSTGA10',
    # VIP interneurons
    'VIPGA1', 'VIPGA2', 'VIPGA3', 'VIPGA4',
    # Lamp5 interneurons
    'LAMGA1', 'LAMGA2', 'LAMGA3', 'LAMGA4',
    # Dentate gyrus neuroblasts (GABAergic immature neurons)
    'DGNBL1', 'DGNBL2',
    # Lateral septal complex GABAergic
    'LSXGA3', 'LSXGA4', 'LSXGA5', 'LSXGA7',
    # Medial septal GABAergic
    'MSGA1', 'MSGA2', 'MSGA4', 'MSGA6', 'MSGA7', 'MSGA8', 'MSGA9',
    'MSGA11', 'MSGA12',
    # Other minor populations
    'CNUGA',      # Cerebellar nuclei GABAergic (sparse in hippocampus)
    'OBNBL',      # Olfactory bulb neuroblasts
    'STRGA2', 'STRGA3',  # Striatal GABAergic
    'D2MSN2', 'D2MSN3',  # D2 medium spiny neurons (very sparse)
}

def is_gaba_subtype(subtype):
    """
    Check if SubType is a hippocampal GABAergic cell type.
    Uses EXACT matching against the curated list of 46 cell types.

    This excludes:
    - Glutamatergic neurons (CA1GL, CA3GL, DGGR, etc.)
    - Non-neuronal cells (PVM, etc.)
    """
    if pd.isna(subtype):
        return False
    return str(subtype) in HIPPOCAMPAL_GABA_CELLTYPES

# =============================================================================
# Step 1: Load and Filter Table 16 for GABA SubTypes
# =============================================================================

print("="*80)
print("STEP 1: Loading and filtering Table 16 for GABA cell types")
print("-"*80)
print(f"Loading: {TABLE16_FILE}")
print(f"Hippocampal GABAergic cell types: {len(HIPPOCAMPAL_GABA_CELLTYPES)}")
print(f"Sample cell types: {', '.join(sorted(list(HIPPOCAMPAL_GABA_CELLTYPES))[:10])}")
print("Using EXACT cell type matching (excludes CA1GL, CA3GL, DGGR, etc.)")
print("Processing in chunks to manage memory...", flush=True)

# Read in chunks and filter immediately
chunk_size = 500000
filtered_chunks = []
total_rows = 0
gaba_rows = 0

for i, chunk in enumerate(pd.read_csv(TABLE16_FILE, sep='\t', chunksize=chunk_size)):
    total_rows += len(chunk)

    # Filter for GABA SubTypes
    chunk['is_gaba'] = chunk['SubType'].apply(is_gaba_subtype)
    gaba_chunk = chunk[chunk['is_gaba']]

    # Apply statistical filters
    gaba_chunk = gaba_chunk[
        (gaba_chunk['FDR'] < FDR_THRESHOLD) &
        (abs(gaba_chunk['PCC']) > PCC_THRESHOLD)
    ]

    if len(gaba_chunk) > 0:
        filtered_chunks.append(gaba_chunk)
        gaba_rows += len(gaba_chunk)

    print(f"  Chunk {i+1}: processed {total_rows:,} rows, kept {gaba_rows:,} GABA links", flush=True)

# Combine filtered chunks
if len(filtered_chunks) == 0:
    print("\nERROR: No GABA cell type links found after filtering!")
    print("Check that Table 16 contains the expected SubTypes.")
    exit(1)

gaba_links = pd.concat(filtered_chunks, ignore_index=True)
print(f"\nTotal rows processed: {total_rows:,}")
print(f"GABA cell type links (filtered): {len(gaba_links):,}")
print("", flush=True)

# Show SubType distribution
print("SubType distribution in GABA links:")
subtype_counts = gaba_links['SubType'].value_counts().head(20)
for subtype, count in subtype_counts.items():
    print(f"  {subtype}: {count:,}")
print("", flush=True)

# =============================================================================
# Step 2: Parse Coordinates and Create BED file
# =============================================================================

print("="*80)
print("STEP 2: Parsing coordinates from Table 16")
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
coords = gaba_links['Coordinate1'].apply(lambda x: pd.Series(parse_coordinate(x)))
gaba_links[['coord_chr', 'coord_start', 'coord_end']] = coords

# Remove invalid coordinates
gaba_links = gaba_links.dropna(subset=['coord_chr', 'coord_start', 'coord_end'])
gaba_links['coord_start'] = gaba_links['coord_start'].astype(int)
gaba_links['coord_end'] = gaba_links['coord_end'].astype(int)

print(f"Parsed {len(gaba_links):,} genomic regions")
print("", flush=True)

# =============================================================================
# Step 3: Save Table 16-only CREs (before ENCODE intersection)
# =============================================================================

print("="*80)
print("STEP 3: Saving Table 16-only CREs (no ENCODE intersection)")
print("-"*80, flush=True)

# Create Table 16-only output with coordinates
table16_only = gaba_links[['coord_chr', 'coord_start', 'coord_end', 'Gene', 'SubType', 'PCC', 'FDR', 'Coordinate1', 'cCRE1']].copy()
table16_only.columns = ['chr', 'start', 'end', 'Gene', 'SubType', 'PCC', 'FDR', 'Coordinate1', 'cCRE1']

# Get unique regions (same coordinates can have multiple gene associations)
table16_unique_regions = table16_only[['chr', 'start', 'end']].drop_duplicates()

print(f"Table 16-only CREs:")
print(f"  Total CRE-gene links: {len(table16_only):,}")
print(f"  Unique regions: {len(table16_unique_regions):,}")
print(f"  Unique genes: {table16_only['Gene'].nunique()}")

# Save Table 16-only output
table16_only.to_csv(OUTPUT_TABLE16_ONLY, sep='\t', index=False)
print(f"\nSaved: {OUTPUT_TABLE16_ONLY}")
print("", flush=True)

# =============================================================================
# Step 4: Use bedtools intersect for FAST overlap detection
# =============================================================================

print("="*80)
print("STEP 4: Finding overlapping ENCODE cCREs (using bedtools)")
print("-"*80, flush=True)

# Reset index and save it for later merging
gaba_links = gaba_links.reset_index(drop=True)

# Create temporary BED file for Table 16 regions
table16_bed = os.path.join(OUTPUT_DIR, "temp_table16_GABA_regions.bed")
with open(table16_bed, 'w') as f:
    for idx, row in gaba_links.iterrows():
        f.write(f"{row['coord_chr']}\t{int(row['coord_start'])}\t{int(row['coord_end'])}\t{idx}\n")

print(f"Created temporary BED file: {table16_bed}")
print(f"  Contains {len(gaba_links):,} GABA-specific regions")

# Run bedtools intersect
print("\nRunning bedtools intersect...", flush=True)

cmd = [
    "bedtools", "intersect",
    "-a", ENCODE_CCRES_FILE,  # ENCODE cCREs
    "-b", table16_bed,         # Table 16 GABA regions
    "-wa", "-wb"               # Report both A and B features
]

print(f"Command: {' '.join(cmd)}", flush=True)

try:
    result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=1800)

    # Parse intersect output
    overlaps = []
    lines = result.stdout.strip().split('\n')
    print(f"\nParsing {len(lines):,} intersection results...", flush=True)

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

    print(f"Found {len(overlaps):,} overlaps")

except subprocess.TimeoutExpired:
    print("ERROR: bedtools intersect timed out after 30 minutes")
    raise
except subprocess.CalledProcessError as e:
    print(f"ERROR running bedtools: {e.stderr}")
    raise
except FileNotFoundError:
    print("ERROR: bedtools not found. Please ensure bedtools is in PATH.")
    raise

# Clean up temp file
os.unlink(table16_bed)
print("Cleaned up temporary files", flush=True)
print("", flush=True)

# =============================================================================
# Step 5: Merge overlap data with Table 16 metadata
# =============================================================================

print("="*80)
print("STEP 5: Merging overlap data with GABA metadata")
print("-"*80, flush=True)

# Create dataframe from overlaps
overlaps_df = pd.DataFrame(overlaps)

if len(overlaps_df) == 0:
    print("WARNING: No overlapping cCREs found!")

    # Create empty output files
    empty_df = pd.DataFrame()
    empty_df.to_csv(OUTPUT_GABA, sep='\t', index=False)
    empty_df.to_csv(OUTPUT_BY_TYPE, sep='\t', index=False)

    with open(OUTPUT_SUMMARY, 'w') as f:
        f.write("No overlapping cCREs found.\n")

    print("\nCreated empty output files.")
    exit(0)

# Merge with Table 16 metadata
merged = overlaps_df.merge(
    gaba_links[['Gene', 'SubType', 'PCC', 'FDR', 'Coordinate1', 'cCRE1']],
    left_on='table16_idx',
    right_index=True,
    how='left'
)

# Drop the temporary index column
merged = merged.drop(columns=['table16_idx'])

print(f"Merged {len(merged):,} GABA-specific cCRE-gene associations")
print(f"Unique ENCODE cCREs: {merged['cCRE_id1'].nunique():,}")
print(f"Unique genes: {merged['Gene'].nunique()}")
print("", flush=True)

# =============================================================================
# Step 6: CRE Type Analysis
# =============================================================================

print("="*80)
print("STEP 6: Analyzing CRE types (intersected CREs only)")
print("-"*80, flush=True)

print("CRE type distribution (GABA cell types):")
type_counts = merged['cre_type'].value_counts()
for cre_type, count in type_counts.items():
    print(f"  {cre_type}: {count:,}")
print("", flush=True)

# =============================================================================
# Step 7: Save Output Files
# =============================================================================

print("="*80)
print("STEP 7: Saving output files")
print("-"*80, flush=True)

# Save intersected GABA-specific CREs
merged.to_csv(OUTPUT_INTERSECTED, sep='\t', index=False)
print(f"Saved: {OUTPUT_INTERSECTED}")
print(f"  {len(merged):,} CRE-gene links (intersected with ENCODE)")

# Save by type (same data, includes cre_type column)
merged.to_csv(OUTPUT_INTERSECTED_BY_TYPE, sep='\t', index=False)
print(f"Saved: {OUTPUT_INTERSECTED_BY_TYPE}")

# Table 16-only already saved in Step 3
print(f"\nTable 16-only output (saved earlier): {OUTPUT_TABLE16_ONLY}")
print(f"  {len(table16_only):,} CRE-gene links (no ENCODE intersection)")
print("", flush=True)

# =============================================================================
# Step 8: Generate Summary Report
# =============================================================================

print("="*80)
print("STEP 8: Generating summary report")
print("-"*80, flush=True)

# Calculate overlap statistics
table16_unique_count = len(table16_unique_regions)
intersected_unique_count = merged['cCRE_id1'].nunique() if len(merged) > 0 else 0
overlap_rate = (intersected_unique_count / table16_unique_count * 100) if table16_unique_count > 0 else 0

with open(OUTPUT_SUMMARY, 'w') as f:
    f.write("="*80 + "\n")
    f.write("GABA CELL TYPE SPECIFIC CREs - SUMMARY\n")
    f.write("="*80 + "\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("INPUT DATA:\n")
    f.write(f"  Table 16 GABA links (filtered): {len(gaba_links):,}\n")
    f.write(f"  Hippocampal GABAergic cell types: {len(HIPPOCAMPAL_GABA_CELLTYPES)}\n")
    f.write(f"  Method: EXACT cell type matching (not keyword-based)\n\n")

    f.write("STATISTICAL FILTERS:\n")
    f.write(f"  FDR threshold: < {FDR_THRESHOLD}\n")
    f.write(f"  |PCC| threshold: > {PCC_THRESHOLD}\n\n")

    f.write("="*80 + "\n")
    f.write("TABLE 16-ONLY CREs (No ENCODE intersection required)\n")
    f.write("="*80 + "\n")
    f.write(f"  Total CRE-gene links: {len(table16_only):,}\n")
    f.write(f"  Unique regions: {table16_unique_count:,}\n")
    f.write(f"  Unique genes: {table16_only['Gene'].nunique()}\n")
    f.write(f"  Output file: {OUTPUT_TABLE16_ONLY}\n\n")

    f.write("="*80 + "\n")
    f.write("ENCODE-INTERSECTED CREs (Table 16 + mm10-cCREs.bed)\n")
    f.write("="*80 + "\n")
    if len(merged) > 0:
        f.write(f"  Total CRE-gene links: {len(merged):,}\n")
        f.write(f"  Unique ENCODE cCREs: {intersected_unique_count:,}\n")
        f.write(f"  Unique genes: {merged['Gene'].nunique()}\n")
    else:
        f.write("  No overlapping CREs found!\n")
    f.write(f"  Output file: {OUTPUT_INTERSECTED}\n\n")

    f.write("="*80 + "\n")
    f.write("OVERLAP STATISTICS\n")
    f.write("="*80 + "\n")
    f.write(f"  Table 16 unique regions: {table16_unique_count:,}\n")
    f.write(f"  Regions with ENCODE overlap: {intersected_unique_count:,}\n")
    f.write(f"  Overlap rate: {overlap_rate:.1f}%\n")
    f.write(f"  Regions WITHOUT ENCODE overlap: {table16_unique_count - intersected_unique_count:,}\n\n")

    if len(merged) > 0:
        f.write("CRE TYPE DISTRIBUTION (Intersected only):\n")
        for cre_type, count in type_counts.items():
            f.write(f"  {cre_type}: {count:,}\n")

        f.write("\nTOP GENES BY NUMBER OF ASSOCIATED cCREs (Intersected):\n")
        top_genes = merged.groupby('Gene').size().sort_values(ascending=False).head(20)
        for gene, count in top_genes.items():
            f.write(f"  {gene}: {count} cCREs\n")

    f.write("\nSUBTYPE DISTRIBUTION (Table 16):\n")
    for subtype, count in subtype_counts.items():
        f.write(f"  {subtype}: {count:,}\n")

    f.write("\nHIPPOCAMPAL GABAERGIC CELL TYPES USED:\n")
    f.write(f"  {', '.join(sorted(HIPPOCAMPAL_GABA_CELLTYPES))}\n")
    f.write("\nEXCLUDED (glutamatergic/non-neuronal):\n")
    f.write("  CA1GL1-3, CA3GL1-6, DGGR, RGDG, PVM\n")

    f.write("\nOUTPUT FILES:\n")
    f.write(f"  Table 16-only: {OUTPUT_TABLE16_ONLY}\n")
    f.write(f"  Intersected: {OUTPUT_INTERSECTED}\n")
    f.write(f"  Intersected by type: {OUTPUT_INTERSECTED_BY_TYPE}\n")

print(f"Saved summary: {OUTPUT_SUMMARY}")
print("", flush=True)

# =============================================================================
# Step 9: Display Results
# =============================================================================

print("="*80)
print("ANALYSIS COMPLETE")
print("="*80)

print(f"\nKey Results:")
print(f"\n  TABLE 16-ONLY CREs (no ENCODE intersection):")
print(f"    Unique regions: {table16_unique_count:,}")
print(f"    Genes represented: {table16_only['Gene'].nunique()}")
print(f"    CRE-gene associations: {len(table16_only):,}")

print(f"\n  ENCODE-INTERSECTED CREs:")
if len(merged) > 0:
    print(f"    Unique ENCODE cCREs: {merged['cCRE_id1'].nunique():,}")
    print(f"    Genes represented: {merged['Gene'].nunique()}")
    print(f"    CRE-gene associations: {len(merged):,}")
else:
    print(f"    No overlapping CREs found")

print(f"\n  OVERLAP STATISTICS:")
print(f"    Table 16 regions: {table16_unique_count:,}")
print(f"    With ENCODE overlap: {intersected_unique_count:,} ({overlap_rate:.1f}%)")
print(f"    Without ENCODE overlap: {table16_unique_count - intersected_unique_count:,}")

if len(merged) > 0:
    print(f"\nTop genes by number of associated cCREs (intersected):")
    top_genes = merged.groupby('Gene').size().sort_values(ascending=False).head(10)
    for gene, count in top_genes.items():
        print(f"  {gene}: {count} cCREs")

print(f"\nOutput files:")
print(f"  Table 16-only: {OUTPUT_TABLE16_ONLY}")
print(f"  Intersected: {OUTPUT_INTERSECTED}")

print(f"\nRecommendation:")
print(f"  - Use TABLE 16-ONLY for maximum coverage (all GABA CREs)")
print(f"  - Use INTERSECTED for high-confidence ENCODE-validated CREs")

print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
