#!/usr/bin/env python3
"""
Extract CREs (Cis-Regulatory Elements) for hippocampal Excitatory neurons

Purpose:
- Extract CREs active in Excitatory (Glutamatergic) neurons
- Use as NEGATIVE CONTROL for GABA ATAC signal specificity
- Expected: GABA samples should show LOW signal at these CREs

Approach: Same data-driven method as hippocampal GABA
- Uses hippocampal tissue samples (CA and DG) from Table 1
- Identifies Excitatory cell types (Class='Glutamate') from Table 2
- Filters for cell types with ≥50 cells
- Extracts corresponding CREs from Table 8
- Outputs merged, deduplicated BED file

INPUT FILES:
- data/table_1.xlsx: Sample and dissection summary
- data/table_2.tsv: Cell metadata with sample assignments
- data/table_8.txt: Cell type assignment of cCREs

OUTPUT FILES:
- output/excitatory_neuron_CREs.bed: BED file with excitatory neuron CREs
- output/excitatory_neuron_CREs_with_header.bed: Same with column headers
- output/excitatory_neuron_CREs_metadata.txt: Detailed metadata
- output/excitatory_neuron_CREs_summary.txt: Summary statistics
"""

import pandas as pd
import numpy as np
from datetime import datetime
import os

print("="*80)
print("EXCITATORY NEURON CRE EXTRACTION (NEGATIVE CONTROL)")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# ============================================================================
# STEP 1: Identify hippocampal Excitatory cell types
# ============================================================================
print("STEP 1: Identifying hippocampal Excitatory (Glutamatergic) cell types...")
print("-"*80)

# Read tables
table1 = pd.read_excel("../data/table_1.xlsx")
table2 = pd.read_csv("../data/table_2.tsv", sep="\t")

# Get hippocampal samples
hpf_samples = table1[table1['Major Region'] == 'HPF']['sample'].tolist()
print(f"Hippocampal samples (HPF region): {len(hpf_samples)}")

# Filter for hippocampal cells
hpf_cells = table2[table2['Sample'].isin(hpf_samples)]
print(f"Total cells from hippocampal samples: {len(hpf_cells):,}")

# Filter for Excitatory (Glutamate) cells
hpf_glut_cells = hpf_cells[hpf_cells['Class'] == 'Glutamate']
print(f"Excitatory (Glutamate) cells: {len(hpf_glut_cells):,}")

# Count cells per type
celltype_counts = hpf_glut_cells['CellType'].value_counts()

# Apply threshold: ≥50 cells (same as GABA)
CELL_THRESHOLD = 50
selected_celltypes = celltype_counts[celltype_counts >= CELL_THRESHOLD].index.tolist()

print(f"\nCell types with ≥{CELL_THRESHOLD} cells: {len(selected_celltypes)}")
print(f"Total cells represented: {celltype_counts[celltype_counts >= CELL_THRESHOLD].sum():,} / {len(hpf_glut_cells):,} ({100*celltype_counts[celltype_counts >= CELL_THRESHOLD].sum()/len(hpf_glut_cells):.1f}%)")

print("\nSelected cell types:")
for celltype in sorted(selected_celltypes):
    count = celltype_counts[celltype]
    print(f"  {celltype}: {count:,} cells")

# ============================================================================
# STEP 2: Extract CREs for selected cell types
# ============================================================================
print(f"\n{'='*80}")
print("STEP 2: Extracting CREs from Table 8...")
print("-"*80)

# Read Table 8
print("Loading Table 8 (this may take a moment - 389MB file)...")
table8 = pd.read_csv("../data/table_8.txt", sep="\t")
print(f"Table 8 loaded: {len(table8):,} total CRE-celltype assignments")

# Filter for selected cell types
print(f"\nFiltering for {len(selected_celltypes)} hippocampal excitatory types...")
excit_cres = table8[table8['CellType'].isin(selected_celltypes)]
print(f"Total CRE entries for excitatory neurons: {len(excit_cres):,}")

# Show breakdown by cell type
print("\nCRE counts by cell type:")
cre_counts = excit_cres['CellType'].value_counts()
for celltype in sorted(selected_celltypes):
    if celltype in cre_counts:
        print(f"  {celltype}: {cre_counts[celltype]:,} CREs")

# ============================================================================
# STEP 3: Parse genomic coordinates and create BED format
# ============================================================================
print(f"\n{'='*80}")
print("STEP 3: Parsing genomic coordinates...")
print("-"*80)

def parse_cre_coordinates(cre_string):
    """Parse CRE string into chromosome, start, end, ID"""
    try:
        coords, cre_id = cre_string.split('|')
        chrom, positions = coords.split(':')
        start, end = positions.split('-')
        return pd.Series({
            'chr': chrom,
            'start': int(start),
            'end': int(end),
            'cre_id': cre_id
        })
    except:
        return pd.Series({
            'chr': None,
            'start': None,
            'end': None,
            'cre_id': None
        })

# Parse coordinates
print("Parsing CRE coordinates...")
parsed_coords = excit_cres['cCREs'].apply(parse_cre_coordinates)
excit_cres_bed = pd.concat([excit_cres, parsed_coords], axis=1)

# Remove any failed parses
excit_cres_bed = excit_cres_bed.dropna(subset=['chr', 'start', 'end'])
print(f"Successfully parsed {len(excit_cres_bed):,} CRE coordinates")

# ============================================================================
# STEP 4: Merge and deduplicate CREs
# ============================================================================
print(f"\n{'='*80}")
print("STEP 4: Merging and deduplicating CREs...")
print("-"*80)

print("Removing duplicate CREs across cell types...")
print(f"Before deduplication: {len(excit_cres_bed):,} entries")

# Group by CRE ID
cre_grouped = excit_cres_bed.groupby('cre_id').agg({
    'chr': 'first',
    'start': 'first',
    'end': 'first',
    'CellType': lambda x: ','.join(sorted(set(x)))
}).reset_index()

cre_grouped.rename(columns={'CellType': 'cell_types'}, inplace=True)

print(f"After deduplication: {len(cre_grouped):,} unique CREs")

# Add number of cell types for each CRE
cre_grouped['n_cell_types'] = cre_grouped['cell_types'].apply(lambda x: len(x.split(',')))

print(f"\nCRE sharing statistics:")
print(f"  CREs active in 1 cell type: {(cre_grouped['n_cell_types'] == 1).sum():,}")
print(f"  CREs active in 2-5 cell types: {((cre_grouped['n_cell_types'] >= 2) & (cre_grouped['n_cell_types'] <= 5)).sum():,}")
print(f"  CREs active in 6-10 cell types: {((cre_grouped['n_cell_types'] >= 6) & (cre_grouped['n_cell_types'] <= 10)).sum():,}")
print(f"  CREs active in >10 cell types: {(cre_grouped['n_cell_types'] > 10).sum():,}")

# ============================================================================
# STEP 5: Create output files
# ============================================================================
print(f"\n{'='*80}")
print("STEP 5: Creating output files...")
print("-"*80)

# Sort by chromosome and position
def chr_sort_key(chr_name):
    """Custom sorting key for chromosome names"""
    chr_name = str(chr_name)
    if chr_name.startswith('chr'):
        chr_name = chr_name[3:]
    if chr_name.isdigit():
        return (0, int(chr_name))
    elif chr_name == 'X':
        return (1, 0)
    elif chr_name == 'Y':
        return (1, 1)
    elif chr_name == 'M' or chr_name == 'MT':
        return (1, 2)
    else:
        return (2, chr_name)

print("Sorting by genomic coordinates...")
cre_grouped['chr_sort'] = cre_grouped['chr'].apply(chr_sort_key)
cre_grouped = cre_grouped.sort_values(['chr_sort', 'start', 'end']).drop('chr_sort', axis=1)

# Create BED file
bed_output = cre_grouped[['chr', 'start', 'end', 'cre_id', 'n_cell_types', 'cell_types']].copy()

# Create output directory
os.makedirs('output', exist_ok=True)

# Save BED file
bed_file = 'output/excitatory_neuron_CREs.bed'
bed_output.to_csv(bed_file, sep='\t', index=False, header=False)
print(f"✓ BED file saved: {bed_file}")

# Save BED file with header
bed_file_header = 'output/excitatory_neuron_CREs_with_header.bed'
bed_output.to_csv(bed_file_header, sep='\t', index=False)
print(f"✓ BED file with header saved: {bed_file_header}")

# Save detailed metadata
metadata_file = 'output/excitatory_neuron_CREs_metadata.txt'
with open(metadata_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("EXCITATORY NEURON CRE EXTRACTION - METADATA\n")
    f.write("="*80 + "\n\n")
    f.write(f"Extraction date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("PURPOSE:\n")
    f.write(f"  Negative control for GABA ATAC signal specificity\n")
    f.write(f"  Expected: LOW signal in GABA samples at these CREs\n\n")

    f.write("APPROACH:\n")
    f.write(f"  Data-driven, conservative (≥{CELL_THRESHOLD} cells per type)\n")
    f.write(f"  Source tissue: Hippocampal formation (CA and DG)\n")
    f.write(f"  Cell class: Glutamatergic (Excitatory) neurons\n\n")

    f.write("STATISTICS:\n")
    f.write(f"  Hippocampal samples: {len(hpf_samples)}\n")
    f.write(f"  Total hippocampal Glutamate cells: {len(hpf_glut_cells):,}\n")
    f.write(f"  Cell types selected: {len(selected_celltypes)}\n")
    f.write(f"  Cells represented: {celltype_counts[celltype_counts >= CELL_THRESHOLD].sum():,} ({100*celltype_counts[celltype_counts >= CELL_THRESHOLD].sum()/len(hpf_glut_cells):.1f}%)\n")
    f.write(f"  Total unique CREs: {len(cre_grouped):,}\n\n")

    f.write("CELL TYPES INCLUDED:\n")
    for celltype in sorted(selected_celltypes):
        n_cells = celltype_counts[celltype]
        n_cres = cre_counts[celltype] if celltype in cre_counts else 0
        f.write(f"  {celltype}: {n_cells:,} cells, {n_cres:,} CRE entries\n")

print(f"✓ Metadata saved: {metadata_file}")

# Save summary statistics
summary_file = 'output/excitatory_neuron_CREs_summary.txt'
with open(summary_file, 'w') as f:
    f.write(f"Total unique CREs: {len(cre_grouped):,}\n")
    f.write(f"Genomic coverage: {(cre_grouped['end'] - cre_grouped['start']).sum():,} bp\n")
    f.write(f"Average CRE length: {(cre_grouped['end'] - cre_grouped['start']).mean():.1f} bp\n")
    f.write(f"Median CRE length: {(cre_grouped['end'] - cre_grouped['start']).median():.1f} bp\n\n")

    f.write("Chromosomal distribution:\n")
    chr_dist = cre_grouped['chr'].value_counts().sort_index()
    for chrom, count in chr_dist.items():
        f.write(f"  {chrom}: {count:,} CREs\n")

print(f"✓ Summary statistics saved: {summary_file}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print(f"\n{'='*80}")
print("EXTRACTION COMPLETE!")
print("="*80)
print(f"\nResults:")
print(f"  Unique CREs extracted: {len(cre_grouped):,}")
print(f"  Cell types included: {len(selected_celltypes)}")
print(f"  Cells represented: {celltype_counts[celltype_counts >= CELL_THRESHOLD].sum():,} / {len(hpf_glut_cells):,} ({100*celltype_counts[celltype_counts >= CELL_THRESHOLD].sum()/len(hpf_glut_cells):.1f}%)")
print(f"  Genomic coverage: {(cre_grouped['end'] - cre_grouped['start']).sum() / 1e6:.1f} Mb")
print(f"\nOutput files in: output/")
print(f"  - excitatory_neuron_CREs.bed")
print(f"  - excitatory_neuron_CREs_with_header.bed")
print(f"  - excitatory_neuron_CREs_metadata.txt")
print(f"  - excitatory_neuron_CREs_summary.txt")

print(f"\n{'='*80}")
print("EXPECTED USE:")
print("="*80)
print("These Excitatory neuron CREs should show LOW ATAC signal in GABA samples.")
print("Compare to GABA CREs (high signal) to demonstrate cell-type specificity.")
print("\nNext: Run heatmap analysis with both GABA and Excitatory CREs side-by-side")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
