#!/usr/bin/env python3
"""
Detailed exploration of hippocampal GABAergic cell types

INPUT FILES:
- data/table_1.xlsx: Sample and dissection summary
- data/table_3.xlsx: Cell cluster annotation
- data/table_8.txt: Cell type assignment of cCREs

OUTPUT FILES:
- None (console output only)
"""

import pandas as pd
import numpy as np

# Read Table 3 - Cell cluster annotation
print("="*80)
print("TABLE 3: Exploring hippocampal GABAergic cell types")
print("="*80)
table3 = pd.read_excel("data/table_3.xlsx")

# Filter for GABAergic cells
gaba_subclasses = ['LAMGA', 'VIPGA', 'PVGA', 'SSTGA', 'LSXGA', 'MSGA', 'CNUGA', 'STRGA']
gaba_cells = table3[table3['Subclass'].isin(gaba_subclasses)]

print(f"\nTotal GABAergic cell types: {len(gaba_cells)}")
print(f"\nColumns: {table3.columns.tolist()}")

# Look at probable locations
print("\n" + "="*80)
print("GABAergic cells by probable location:")
print("="*80)
location_counts = gaba_cells['Probable location (manually assigned)'].value_counts()
print(location_counts)

# Filter for hippocampal cells
print("\n" + "="*80)
print("Hippocampal GABAergic cell types:")
print("="*80)
hippocampal_keywords = ['hippocampus', 'hippo', 'CA1', 'CA2', 'CA3', 'DG', 'dentate']
hippocampal_gaba = gaba_cells[
    gaba_cells['Probable location (manually assigned)'].str.lower().str.contains(
        '|'.join(hippocampal_keywords), na=False
    )
]

print(f"\nFound {len(hippocampal_gaba)} hippocampal GABAergic cell types:")
print(hippocampal_gaba[['Subclass', 'Cell type', 'Probable location (manually assigned)', 'NucleiCounts']])

# Save the cell type IDs
hippocampal_cell_types = hippocampal_gaba['Cell type'].tolist()
print(f"\nHippocampal GABAergic cell type IDs:")
for ct in hippocampal_cell_types:
    print(f"  {ct}")

# Look at Table 1 to understand dissection regions
print("\n" + "="*80)
print("TABLE 1: Hippocampal dissection regions")
print("="*80)
table1 = pd.read_excel("data/table_1.xlsx")
print(f"\nColumns: {table1.columns.tolist()}")

# Check for hippocampus-related regions
hippocampal_regions = table1[
    table1['Region Name'].str.lower().str.contains(
        '|'.join(hippocampal_keywords), na=False
    )
]

if len(hippocampal_regions) > 0:
    print(f"\nFound {len(hippocampal_regions)} hippocampal dissection regions:")
    print(hippocampal_regions[['DissectionRegion', 'Region Name', 'Major Region', 'Sub Region', 'Detail Region']])
else:
    # Try other columns
    print("\nSearching in all region columns...")
    for col in ['Major Region', 'Sub Region', 'Detail Region', 'DissectionRegion']:
        if col in table1.columns:
            matches = table1[
                table1[col].astype(str).str.lower().str.contains(
                    '|'.join(hippocampal_keywords), na=False
                )
            ]
            if len(matches) > 0:
                print(f"\nMatches in {col}:")
                print(matches[['DissectionRegion', 'Region Name', col]].drop_duplicates())

# Look at all unique regions to manually identify hippocampus
print("\n" + "="*80)
print("All unique region names (to manually identify hippocampal regions):")
print("="*80)
for col in ['Region Name', 'Major Region', 'Sub Region', 'Detail Region']:
    if col in table1.columns:
        print(f"\n{col}:")
        print(table1[col].unique())

# Look at Table 8 structure
print("\n" + "="*80)
print("TABLE 8: Understanding cell type naming")
print("="*80)
table8 = pd.read_csv("data/table_8.txt", sep="\t")
print(f"\nUnique cell types in Table 8: {len(table8['CellType'].unique())}")
print("\nSample of cell types:")
print(sorted(table8['CellType'].unique())[:50])

# Check if hippocampal cell types from Table 3 match Table 8
print("\n" + "="*80)
print("Matching hippocampal cell types between Table 3 and Table 8:")
print("="*80)
table8_types = set(table8['CellType'].unique())
for ct in hippocampal_cell_types:
    # Check exact match
    if ct in table8_types:
        count = (table8['CellType'] == ct).sum()
        print(f"  ✓ {ct}: {count} cCREs")
    else:
        # Check if there are numbered variants (e.g., PVGA1, PVGA2 for PVGA)
        variants = [t for t in table8_types if t.startswith(ct)]
        if variants:
            total = sum((table8['CellType'] == v).sum() for v in variants)
            print(f"  ~ {ct} has {len(variants)} subtypes: {variants[:5]}... ({total} total cCREs)")
        else:
            print(f"  ✗ {ct}: NOT FOUND in Table 8")
