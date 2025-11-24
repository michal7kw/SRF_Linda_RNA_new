#!/usr/bin/env python3
"""
Explore different approaches to identify hippocampal interneuron CREs

INPUT FILES:
- data/table_1.xlsx: Sample and dissection summary
- data/table_3.xlsx: Cell cluster annotation
- data/table_8.txt: Cell type assignment of cCREs

OUTPUT FILES:
- None (console output only)
"""

import pandas as pd

# Read the tables
table1 = pd.read_excel("data/table_1.xlsx")
table3 = pd.read_excel("data/table_3.xlsx")
table8 = pd.read_csv("data/table_8.txt", sep="\t")

print("="*80)
print("APPROACH 1: Strict hippocampal annotation from Table 3")
print("="*80)
print("\nOnly LAMGA1 is explicitly marked as 'Hippocampus' in Table 3")
print(f"LAMGA1 has {(table8['CellType'] == 'LAMGA1').sum():,} cCREs")

print("\n" + "="*80)
print("APPROACH 2: All GABAergic interneuron types (pan-brain)")
print("="*80)
gaba_subclasses = ['LAMGA', 'VIPGA', 'PVGA', 'SSTGA', 'LSXGA', 'MSGA', 'CNUGA', 'STRGA']
gaba_cells = table3[table3['Subclass'].isin(gaba_subclasses)]
print(f"\nFound {len(gaba_cells)} GABAergic cell types across all brain regions")
print("\nBreakdown by subclass:")
for subclass in gaba_subclasses:
    cells = table3[table3['Subclass'] == subclass]
    # Find matching cell types in Table 8
    matching_types = [t for t in table8['CellType'].unique() if t.startswith(subclass)]
    n_cres = table8[table8['CellType'].isin(matching_types)].shape[0]
    print(f"  {subclass}: {len(cells)} cell types, {len(matching_types)} subtypes in Table 8, {n_cres:,} cCREs")

print("\n" + "="*80)
print("APPROACH 3: Hippocampal samples (using Table 1 dissections)")
print("="*80)
# Identify hippocampal dissections
hpf_samples = table1[table1['Major Region'] == 'HPF']
print(f"\nFound {len(hpf_samples)} samples from HPF (Hippocampal Formation):")
print(hpf_samples[['sample', 'DissectionRegion', 'Sub Region', 'Detail Region']].to_string(index=False))

print("\n" + "="*80)
print("APPROACH 4: Look for hippocampus-specific cell types in Table 8")
print("="*80)
# Cell types with CA (Cornu Ammonis) or DG (Dentate Gyrus) in name
hippo_keywords = ['CA1', 'CA2', 'CA3', 'DG']
hippo_cell_types = [ct for ct in table8['CellType'].unique()
                    if any(kw in ct for kw in hippo_keywords)]
print(f"\nFound {len(hippo_cell_types)} cell types with hippocampal markers:")
for ct in sorted(hippo_cell_types):
    n_cres = (table8['CellType'] == ct).sum()
    print(f"  {ct}: {n_cres:,} cCREs")

# Check which are GABAergic
hippo_gaba = [ct for ct in hippo_cell_types if 'GA' in ct or any(ct.startswith(g) for g in gaba_subclasses)]
if hippo_gaba:
    print(f"\nHippocampal GABAergic cell types: {hippo_gaba}")
else:
    print("\nNote: Hippocampal cell types (CA1GL, CA3GL, DGGR, DGNBL) are mostly glutamatergic")

print("\n" + "="*80)
print("RECOMMENDATION")
print("="*80)
print("""
Option 1 (STRICT): Use only LAMGA1 (~59k cCREs)
  - Most conservative
  - Only explicitly hippocampal interneurons

Option 2 (MODERATE): Use major hippocampal interneuron types
  - PVGA (Parvalbumin) - basket cells, chandelier cells
  - SSTGA (Somatostatin) - O-LM cells, bistratified cells
  - VIPGA (VIP) - interneuron-specific interneurons
  - LAMGA (various types)
  - These are well-characterized hippocampal interneuron markers

Option 3 (BROAD): Use all GABAergic subtypes
  - Includes all brain regions
  - Most comprehensive but includes non-hippocampal cells

Which approach would you prefer?
""")
