#!/usr/bin/env python3
"""
Explore CRE-gene linkage tables (Tables 16-18)

INPUT FILES:
- data/table_16.*: Summary of gene-cCRE correlations (various formats)
- data/table_17.*: Association modules of enhancer-gene pairs (various formats)
- data/table_18.*: Association modules with individual putative enhancers (various formats)

OUTPUT FILES:
- None (console output only)
"""
import pandas as pd
import numpy as np

print("="*80)
print("TABLE 16: Summary of gene-cCRE correlations")
print("="*80)

# Try different file formats
import os
files_16 = [f for f in os.listdir('data') if 'table_16' in f.lower() or '16' in f]
print(f"Found files: {files_16}")

for file in files_16:
    filepath = f'data/{file}'
    print(f"\nTrying to read: {filepath}")
    try:
        if file.endswith('.xlsx'):
            table16 = pd.read_excel(filepath)
        elif file.endswith('.tsv') or file.endswith('.txt'):
            table16 = pd.read_csv(filepath, sep="\t")
        elif file.endswith('.csv'):
            table16 = pd.read_csv(filepath)
        else:
            # Try tab-separated by default
            table16 = pd.read_csv(filepath, sep="\t")

        print(f"✓ Successfully read {file}")
        print(f"Shape: {table16.shape}")
        print(f"Columns: {table16.columns.tolist()}")
        print("\nFirst 10 rows:")
        print(table16.head(10))

        # Check for gene names
        print("\nSample gene names:")
        gene_cols = [col for col in table16.columns if 'gene' in col.lower() or 'symbol' in col.lower()]
        if gene_cols:
            print(f"Gene columns: {gene_cols}")
            print(table16[gene_cols[0]].head(20).tolist())

        # Check for CRE identifiers
        print("\nSample CRE identifiers:")
        cre_cols = [col for col in table16.columns if 'cre' in col.lower() or 'ccre' in col.lower() or 'peak' in col.lower()]
        if cre_cols:
            print(f"CRE columns: {cre_cols}")
            print(table16[cre_cols[0]].head(20).tolist())

    except Exception as e:
        print(f"✗ Failed to read {file}: {e}")

print("\n" + "="*80)
print("TABLE 17: Association modules of enhancer-gene pairs")
print("="*80)

files_17 = [f for f in os.listdir('data') if 'table_17' in f.lower() or '17' in f]
print(f"Found files: {files_17}")

for file in files_17:
    filepath = f'data/{file}'
    print(f"\nTrying to read: {filepath}")
    try:
        if file.endswith('.xlsx'):
            table17 = pd.read_excel(filepath)
        elif file.endswith('.tsv') or file.endswith('.txt'):
            table17 = pd.read_csv(filepath, sep="\t")
        elif file.endswith('.csv'):
            table17 = pd.read_csv(filepath)
        else:
            table17 = pd.read_csv(filepath, sep="\t")

        print(f"✓ Successfully read {file}")
        print(f"Shape: {table17.shape}")
        print(f"Columns: {table17.columns.tolist()}")
        print("\nFirst 10 rows:")
        print(table17.head(10))

    except Exception as e:
        print(f"✗ Failed to read {file}: {e}")

print("\n" + "="*80)
print("TABLE 18: Association modules with individual putative enhancers")
print("="*80)

files_18 = [f for f in os.listdir('data') if 'table_18' in f.lower() or '18' in f]
print(f"Found files: {files_18}")

for file in files_18:
    filepath = f'data/{file}'
    print(f"\nTrying to read: {filepath}")
    try:
        if file.endswith('.xlsx'):
            table18 = pd.read_excel(filepath)
        elif file.endswith('.tsv') or file.endswith('.txt'):
            table18 = pd.read_csv(filepath, sep="\t")
        elif file.endswith('.csv'):
            table18 = pd.read_csv(filepath)
        else:
            table18 = pd.read_csv(filepath, sep="\t")

        print(f"✓ Successfully read {file}")
        print(f"Shape: {table18.shape}")
        print(f"Columns: {table18.columns.tolist()}")
        print("\nFirst 10 rows:")
        print(table18.head(10))

        # Check for gene and CRE info
        print("\nChecking for gene-CRE pairs...")
        gene_cols = [col for col in table18.columns if 'gene' in col.lower()]
        cre_cols = [col for col in table18.columns if 'cre' in col.lower() or 'ccre' in col.lower() or 'peak' in col.lower() or 'enhancer' in col.lower()]

        if gene_cols and cre_cols:
            print(f"\nGene columns: {gene_cols}")
            print(f"CRE columns: {cre_cols}")
            print("\nSample gene-CRE pairs:")
            print(table18[[gene_cols[0], cre_cols[0]]].head(20))

    except Exception as e:
        print(f"✗ Failed to read {file}: {e}")

# List all files in data directory
print("\n" + "="*80)
print("All files in data/ directory:")
print("="*80)
for file in sorted(os.listdir('data')):
    filepath = f'data/{file}'
    size = os.path.getsize(filepath) / (1024*1024)  # MB
    print(f"{file:40s} {size:8.1f} MB")
