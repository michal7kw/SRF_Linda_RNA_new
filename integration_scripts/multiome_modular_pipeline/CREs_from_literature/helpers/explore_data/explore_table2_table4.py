#!/usr/bin/env python3
"""
Explore Table 2 to link hippocampal samples to GABAergic cell types

INPUT FILES:
- ../../data/table_1.xlsx: Sample and dissection summary
- ../../data/table_2.tsv: Cell metadata with sample assignments
- ../../data/table_3.xlsx: Cell cluster annotation
- ../../data/table_4.xlsx: Cell type annotations

OUTPUT FILES:
- None (console output only)
"""
import pandas as pd

# Read tables
table1 = pd.read_excel("../../data/table_1.xlsx")
table2 = pd.read_csv("../../data/table_2.tsv", sep="\t")
table3 = pd.read_excel("../../data/table_3.xlsx")
table4 = pd.read_excel("../../data/table_4.xlsx")

print("="*80)
print("TABLE 2 Structure")
print("="*80)
print(f"Shape: {table2.shape}")
print(f"Columns: {table2.columns.tolist()}")
print("\nFirst 10 rows:")
print(table2.head(10))

print("\n" + "="*80)
print("TABLE 4 Structure")
print("="*80)
print(f"Shape: {table4.shape}")
print(f"Columns: {table4.columns.tolist()}")
print("\nFirst 10 rows:")
print(table4.head(10))

# Get hippocampal samples from Table 1
hpf_samples = table1[table1['Major Region'] == 'HPF']['sample'].tolist()
print(f"\n{'='*80}")
print(f"Hippocampal samples (n={len(hpf_samples)}):")
print(f"{'='*80}")
print(hpf_samples)

# Try to link samples to cell types via Table 2
print(f"\n{'='*80}")
print("Linking hippocampal samples to cell types via Table 2")
print(f"{'='*80}")
# Check if table2 has sample or cluster information
if 'sample' in table2.columns:
    hpf_table2 = table2[table2['sample'].isin(hpf_samples)]
    print(f"Found {len(hpf_table2)} entries from hippocampal samples")
    print(hpf_table2.head(20))
