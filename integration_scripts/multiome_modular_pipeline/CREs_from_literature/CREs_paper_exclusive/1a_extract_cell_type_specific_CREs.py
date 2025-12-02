#!/usr/bin/env python3
"""
Extract Cell-Type-SPECIFIC CREs (Mutually Exclusive)

Problem with previous approach:
- 60% of CREs overlap between GABA and Excitatory sets
- This causes identical signal patterns
- Not useful as negative control

New approach:
- Extract CREs that are EXCLUSIVE to each cell type
- GABA-specific: Active in GABA neurons but NOT in Excitatory neurons
- Excitatory-specific: Active in Excitatory neurons but NOT in GABA neurons
- This ensures they serve as true positive/negative controls

INPUT FILES:
- data/table_1.xlsx: Sample and dissection summary
- data/table_2.tsv: Cell metadata with sample assignments
- data/table_8.txt: Cell type assignment of cCREs

OUTPUT FILES:
- output/GABA_specific_CREs.bed: CREs ONLY active in GABA neurons
- output/Excitatory_specific_CREs.bed: CREs ONLY active in Excitatory neurons
- output/cell_type_specific_CREs_summary.txt: Overlap analysis
"""

import pandas as pd
import numpy as np
from datetime import datetime
import os

print("="*80)
print("CELL-TYPE-SPECIFIC CRE EXTRACTION")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# ============================================================================
# STEP 1: Load data and identify cell types
# ============================================================================
print("STEP 1: Loading data and identifying cell types...")
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

# Separate GABA and Excitatory cells
hpf_gaba_cells = hpf_cells[hpf_cells['Class'] == 'GABA']
hpf_glut_cells = hpf_cells[hpf_cells['Class'] == 'Glutamate']

print(f"\nGABA cells: {len(hpf_gaba_cells):,}")
print(f"Excitatory (Glutamate) cells: {len(hpf_glut_cells):,}")

# Get cell types with ≥50 cells
CELL_THRESHOLD = 50

gaba_celltype_counts = hpf_gaba_cells['CellType'].value_counts()
gaba_celltypes = gaba_celltype_counts[gaba_celltype_counts >= CELL_THRESHOLD].index.tolist()

glut_celltype_counts = hpf_glut_cells['CellType'].value_counts()
glut_celltypes = glut_celltype_counts[glut_celltype_counts >= CELL_THRESHOLD].index.tolist()

print(f"\nGABA cell types (≥{CELL_THRESHOLD} cells): {len(gaba_celltypes)}")
print(f"Excitatory cell types (≥{CELL_THRESHOLD} cells): {len(glut_celltypes)}")

# ============================================================================
# STEP 2: Extract CREs for each class
# ============================================================================
print(f"\n{'='*80}")
print("STEP 2: Extracting CREs from Table 8...")
print("-"*80)

# Read Table 8
print("Loading Table 8 (this may take a moment - 389MB file)...")
table8 = pd.read_csv("../data/table_8.txt", sep="\t")
print(f"Table 8 loaded: {len(table8):,} total CRE-celltype assignments")

# Filter for GABA and Excitatory cell types
gaba_cres = table8[table8['CellType'].isin(gaba_celltypes)]
glut_cres = table8[table8['CellType'].isin(glut_celltypes)]

print(f"\nCRE entries:")
print(f"  GABA: {len(gaba_cres):,}")
print(f"  Excitatory: {len(glut_cres):,}")

# ============================================================================
# STEP 3: Parse coordinates and identify unique CRE IDs
# ============================================================================
print(f"\n{'='*80}")
print("STEP 3: Identifying unique CRE IDs for each class...")
print("-"*80)

# Get unique CRE IDs for each class
gaba_cre_ids = set(gaba_cres['cCREs'].str.split('|').str[1])
glut_cre_ids = set(glut_cres['cCREs'].str.split('|').str[1])

print(f"\nUnique CRE IDs:")
print(f"  GABA: {len(gaba_cre_ids):,}")
print(f"  Excitatory: {len(glut_cre_ids):,}")

# Find overlapping and exclusive CREs
overlap_cre_ids = gaba_cre_ids & glut_cre_ids
gaba_specific_ids = gaba_cre_ids - glut_cre_ids
glut_specific_ids = glut_cre_ids - gaba_cre_ids

print(f"\nOverlap analysis:")
print(f"  Shared CREs: {len(overlap_cre_ids):,} ({100*len(overlap_cre_ids)/len(gaba_cre_ids):.1f}% of GABA, {100*len(overlap_cre_ids)/len(glut_cre_ids):.1f}% of Excitatory)")
print(f"  GABA-specific CREs: {len(gaba_specific_ids):,}")
print(f"  Excitatory-specific CREs: {len(glut_specific_ids):,}")

# ============================================================================
# STEP 4: Create BED files for specific CREs
# ============================================================================
print(f"\n{'='*80}")
print("STEP 4: Creating BED files for cell-type-specific CREs...")
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

# GABA-specific CREs
print("\nProcessing GABA-specific CREs...")
gaba_specific = gaba_cres[gaba_cres['cCREs'].str.split('|').str[1].isin(gaba_specific_ids)]
gaba_specific_parsed = gaba_specific['cCREs'].apply(parse_cre_coordinates)
gaba_specific_bed = pd.concat([gaba_specific, gaba_specific_parsed], axis=1)
gaba_specific_bed = gaba_specific_bed.dropna(subset=['chr', 'start', 'end'])

# Group by CRE ID and deduplicate
gaba_specific_grouped = gaba_specific_bed.groupby('cre_id').agg({
    'chr': 'first',
    'start': 'first',
    'end': 'first',
    'CellType': lambda x: ','.join(sorted(set(x)))
}).reset_index()
gaba_specific_grouped.rename(columns={'CellType': 'cell_types'}, inplace=True)
gaba_specific_grouped['n_cell_types'] = gaba_specific_grouped['cell_types'].apply(lambda x: len(x.split(',')))

print(f"  Unique GABA-specific CREs: {len(gaba_specific_grouped):,}")

# Excitatory-specific CREs
print("Processing Excitatory-specific CREs...")
glut_specific = glut_cres[glut_cres['cCREs'].str.split('|').str[1].isin(glut_specific_ids)]
glut_specific_parsed = glut_specific['cCREs'].apply(parse_cre_coordinates)
glut_specific_bed = pd.concat([glut_specific, glut_specific_parsed], axis=1)
glut_specific_bed = glut_specific_bed.dropna(subset=['chr', 'start', 'end'])

# Group by CRE ID and deduplicate
glut_specific_grouped = glut_specific_bed.groupby('cre_id').agg({
    'chr': 'first',
    'start': 'first',
    'end': 'first',
    'CellType': lambda x: ','.join(sorted(set(x)))
}).reset_index()
glut_specific_grouped.rename(columns={'CellType': 'cell_types'}, inplace=True)
glut_specific_grouped['n_cell_types'] = glut_specific_grouped['cell_types'].apply(lambda x: len(x.split(',')))

print(f"  Unique Excitatory-specific CREs: {len(glut_specific_grouped):,}")

# ============================================================================
# STEP 5: Sort and save BED files
# ============================================================================
print(f"\n{'='*80}")
print("STEP 5: Sorting and saving BED files...")
print("-"*80)

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

# Sort GABA-specific
gaba_specific_grouped['chr_sort'] = gaba_specific_grouped['chr'].apply(chr_sort_key)
gaba_specific_grouped = gaba_specific_grouped.sort_values(['chr_sort', 'start', 'end']).drop('chr_sort', axis=1)

# Sort Excitatory-specific
glut_specific_grouped['chr_sort'] = glut_specific_grouped['chr'].apply(chr_sort_key)
glut_specific_grouped = glut_specific_grouped.sort_values(['chr_sort', 'start', 'end']).drop('chr_sort', axis=1)

# Create output directory
os.makedirs('output', exist_ok=True)

# Save GABA-specific BED
gaba_bed_output = gaba_specific_grouped[['chr', 'start', 'end', 'cre_id', 'n_cell_types', 'cell_types']].copy()
gaba_bed_file = 'output/GABA_specific_CREs.bed'
gaba_bed_output.to_csv(gaba_bed_file, sep='\t', index=False, header=False)
print(f"✓ GABA-specific BED file saved: {gaba_bed_file}")

# Save Excitatory-specific BED
glut_bed_output = glut_specific_grouped[['chr', 'start', 'end', 'cre_id', 'n_cell_types', 'cell_types']].copy()
glut_bed_file = 'output/Excitatory_specific_CREs.bed'
glut_bed_output.to_csv(glut_bed_file, sep='\t', index=False, header=False)
print(f"✓ Excitatory-specific BED file saved: {glut_bed_file}")

# ============================================================================
# STEP 6: Save summary report
# ============================================================================
print(f"\n{'='*80}")
print("STEP 6: Creating summary report...")
print("-"*80)

summary_file = 'output/cell_type_specific_CREs_summary.txt'
with open(summary_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("CELL-TYPE-SPECIFIC CRE EXTRACTION - SUMMARY\n")
    f.write("="*80 + "\n\n")
    f.write(f"Extraction date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("APPROACH:\n")
    f.write("  Extract CREs that are EXCLUSIVE to each cell type\n")
    f.write("  GABA-specific: Active in GABA neurons but NOT in Excitatory neurons\n")
    f.write("  Excitatory-specific: Active in Excitatory neurons but NOT in GABA neurons\n\n")

    f.write("STATISTICS:\n")
    f.write(f"  Total CRE IDs in GABA neurons: {len(gaba_cre_ids):,}\n")
    f.write(f"  Total CRE IDs in Excitatory neurons: {len(glut_cre_ids):,}\n")
    f.write(f"  Shared (overlapping) CREs: {len(overlap_cre_ids):,}\n")
    f.write(f"    - {100*len(overlap_cre_ids)/len(gaba_cre_ids):.1f}% of GABA CREs\n")
    f.write(f"    - {100*len(overlap_cre_ids)/len(glut_cre_ids):.1f}% of Excitatory CREs\n\n")

    f.write("CELL-TYPE-SPECIFIC CREs:\n")
    f.write(f"  GABA-specific CREs: {len(gaba_specific_grouped):,}\n")
    f.write(f"  Excitatory-specific CREs: {len(glut_specific_grouped):,}\n\n")

    f.write("OUTPUT FILES:\n")
    f.write(f"  - {gaba_bed_file}\n")
    f.write(f"  - {glut_bed_file}\n\n")

    f.write("EXPECTED RESULTS:\n")
    f.write("  GABA samples should show:\n")
    f.write("    - HIGH signal at GABA-specific CREs (positive control)\n")
    f.write("    - LOW/NO signal at Excitatory-specific CREs (negative control)\n\n")

    f.write("NEXT STEPS:\n")
    f.write("  1. Run heatmap analysis using GABA_specific_CREs.bed and Excitatory_specific_CREs.bed\n")
    f.write("  2. Expect clear visual difference between the two sets\n")
    f.write("  3. Calculate fold-enrichment to quantify specificity\n")

print(f"✓ Summary report saved: {summary_file}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print(f"\n{'='*80}")
print("EXTRACTION COMPLETE!")
print("="*80)
print(f"\nResults:")
print(f"  GABA-specific CREs: {len(gaba_specific_grouped):,}")
print(f"  Excitatory-specific CREs: {len(glut_specific_grouped):,}")
print(f"  Overlap eliminated: {len(overlap_cre_ids):,} shared CREs excluded")
print(f"\nOutput files in: output/")
print(f"  - GABA_specific_CREs.bed")
print(f"  - Excitatory_specific_CREs.bed")
print(f"  - cell_type_specific_CREs_summary.txt")

print(f"\n{'='*80}")
print("EXPECTED USE:")
print("="*80)
print("These are MUTUALLY EXCLUSIVE CRE sets:")
print("  - GABA-specific: Should show HIGH signal in GABA ATAC samples")
print("  - Excitatory-specific: Should show LOW signal in GABA ATAC samples")
print("\nThis ensures a true positive/negative control comparison!")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
