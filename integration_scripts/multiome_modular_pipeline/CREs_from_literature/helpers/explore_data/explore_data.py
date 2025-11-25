#!/usr/bin/env python3
"""
Explore supplementary tables to understand data structure

INPUT FILES:
- ../../data/table_1.xlsx: Sample and dissection summary
- ../../data/table_3.xlsx: Cell cluster annotation
- ../../data/table_7.txt: cCREs list
- ../../data/table_8.txt: Cell type assignment of cCREs

OUTPUT FILES:
- None (console output only)
"""

import pandas as pd
import numpy as np

# Read the Excel files
print("="*80)
print("TABLE 1: Sample and dissection summary")
print("="*80)
table1 = pd.read_excel("../../data/table_1.xlsx")
print(f"Shape: {table1.shape}")
print(f"Columns: {table1.columns.tolist()}")
print("\nFirst few rows:")
print(table1.head(10))
print("\nUnique dissection regions:")
if 'Dissection' in table1.columns:
    print(table1['Dissection'].unique())
elif 'Region' in table1.columns:
    print(table1['Region'].unique())
print("\n")

print("="*80)
print("TABLE 3: Cell cluster annotation")
print("="*80)
table3 = pd.read_excel("../../data/table_3.xlsx")
print(f"Shape: {table3.shape}")
print(f"Columns: {table3.columns.tolist()}")
print("\nFirst few rows:")
print(table3.head(10))
print("\nUnique cell types (sample):")
if 'CellType' in table3.columns:
    cell_types = table3['CellType'].unique()
    print(cell_types[:20])
elif 'Subclass' in table3.columns:
    cell_types = table3['Subclass'].unique()
    print(cell_types[:20])
print("\n")

# Read text files
print("="*80)
print("TABLE 7: cCREs list")
print("="*80)
table7 = pd.read_csv("../../data/table_7.txt", sep="\t")
print(f"Shape: {table7.shape}")
print(f"Columns: {table7.columns.tolist()}")
print("\nFirst few rows:")
print(table7.head(10))
print("\n")

print("="*80)
print("TABLE 8: Cell type assignment of cCREs")
print("="*80)
table8 = pd.read_csv("../../data/table_8.txt", sep="\t")
print(f"Shape: {table8.shape}")
print(f"Columns: {table8.columns.tolist()}")
print("\nFirst few rows:")
print(table8.head(10))
print("\nUnique cell types (sample):")
print(table8['CellType'].unique()[:30])
print("\n")

# Look for GABAergic cell types
print("="*80)
print("Searching for GABAergic/Interneuron cell types in Table 8")
print("="*80)
gaba_patterns = ['GA', 'GABA', 'interneuron', 'Inter', 'SST', 'PV', 'VIP', 'LAMP']
gaba_types = [ct for ct in table8['CellType'].unique() if any(pattern.lower() in ct.lower() for pattern in gaba_patterns)]
print(f"Found {len(gaba_types)} GABAergic-like cell types:")
for ct in sorted(gaba_types):
    count = (table8['CellType'] == ct).sum()
    print(f"  {ct}: {count} cCREs")
