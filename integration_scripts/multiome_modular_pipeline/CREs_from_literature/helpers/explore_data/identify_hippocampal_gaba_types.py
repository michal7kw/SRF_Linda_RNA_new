#!/usr/bin/env python3
"""
Identify GABAergic cell types present in hippocampal samples

INPUT FILES:
- data/table_1.xlsx: Sample and dissection summary
- data/table_2.tsv: Cell metadata with sample assignments
- data/table_8.txt: Cell type assignment of cCREs

OUTPUT FILES:
- None (console output only)
"""
import pandas as pd

# Read tables
table1 = pd.read_excel("data/table_1.xlsx")
table2 = pd.read_csv("data/table_2.tsv", sep="\t")
table8 = pd.read_csv("data/table_8.txt", sep="\t")

# Get hippocampal samples from Table 1
hpf_samples = table1[table1['Major Region'] == 'HPF']['sample'].tolist()
print(f"{'='*80}")
print(f"Hippocampal samples from Table 1 (n={len(hpf_samples)}):")
print(f"{'='*80}")
for sample in hpf_samples:
    region_info = table1[table1['sample'] == sample][['DissectionRegion', 'Sub Region', 'Detail Region']].iloc[0]
    print(f"  {sample}: {region_info['Sub Region']} ({region_info['Detail Region']})")

# Filter Table 2 for hippocampal samples
print(f"\n{'='*80}")
print("Cells from hippocampal samples in Table 2:")
print(f"{'='*80}")
hpf_cells = table2[table2['Sample'].isin(hpf_samples)]
print(f"Total cells from hippocampal samples: {len(hpf_cells):,}")

# Breakdown by cell class
print(f"\nBreakdown by Class:")
class_counts = hpf_cells['Class'].value_counts()
for cls, count in class_counts.items():
    pct = 100 * count / len(hpf_cells)
    print(f"  {cls}: {count:,} cells ({pct:.1f}%)")

# Filter for GABAergic cells only
print(f"\n{'='*80}")
print("GABAergic cells from hippocampal samples:")
print(f"{'='*80}")
hpf_gaba_cells = hpf_cells[hpf_cells['Class'] == 'GABA']
print(f"Total GABAergic cells: {len(hpf_gaba_cells):,}")

# Get unique GABAergic cell types
print(f"\nUnique GABAergic SubClasses in hippocampus:")
subclass_counts = hpf_gaba_cells['SubClass'].value_counts()
for subclass, count in subclass_counts.items():
    pct = 100 * count / len(hpf_gaba_cells)
    print(f"  {subclass}: {count:,} cells ({pct:.1f}%)")

print(f"\nUnique GABAergic CellTypes in hippocampus:")
celltype_counts = hpf_gaba_cells['CellType'].value_counts()
for celltype, count in celltype_counts.items():
    pct = 100 * count / len(hpf_gaba_cells)
    print(f"  {celltype}: {count:,} cells ({pct:.1f}%)")

# Now check which of these cell types have CREs in Table 8
print(f"\n{'='*80}")
print("Matching hippocampal GABAergic cell types to CREs in Table 8:")
print(f"{'='*80}")
hippocampal_gaba_types = hpf_gaba_cells['CellType'].unique()
table8_types = set(table8['CellType'].unique())

total_cres = 0
matched_types = []
for celltype in sorted(hippocampal_gaba_types):
    if celltype in table8_types:
        n_cres = (table8['CellType'] == celltype).sum()
        total_cres += n_cres
        matched_types.append(celltype)
        cell_count = (hpf_gaba_cells['CellType'] == celltype).sum()
        print(f"  ✓ {celltype}: {n_cres:,} cCREs ({cell_count:,} cells in hippocampus)")
    else:
        print(f"  ✗ {celltype}: NOT FOUND in Table 8")

print(f"\n{'='*80}")
print(f"SUMMARY:")
print(f"{'='*80}")
print(f"Total hippocampal GABAergic cell types: {len(hippocampal_gaba_types)}")
print(f"Cell types with CREs in Table 8: {len(matched_types)}")
print(f"Total CREs for hippocampal GABAergic interneurons: {total_cres:,}")

print(f"\n{'='*80}")
print(f"Cell types to extract (for final script):")
print(f"{'='*80}")
print(matched_types)
