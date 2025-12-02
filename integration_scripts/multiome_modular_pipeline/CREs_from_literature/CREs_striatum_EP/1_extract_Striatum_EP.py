#!/usr/bin/env python3
"""
Extract Striatum_EP CREs Associated with Splicing Genes

This script parses the Striatum_EP.txt file and filters for CREs associated with
splicing-related genes.

Input Format:
chr1:4775480-4775650_ENSMUSG00000025902$Sox17$chr1$4486494      0.600838

Usage:
    python 1_extract_Striatum_EP.py
"""

import pandas as pd
import os
import re
from datetime import datetime

# Configuration
BASE_DIR = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/Striatum_EP_analysis"
os.chdir(BASE_DIR)

INPUT_FILE = "../Striatum_EP.txt"
SPLICING_GENES_FILE = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/splicing_genes/extracted_genes_final.csv"

OUTPUT_DIR = "./output"
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_TSV = os.path.join(OUTPUT_DIR, "Striatum_EP_splicing_genes.tsv")
OUTPUT_BED = os.path.join(OUTPUT_DIR, "Striatum_EP_splicing_genes.bed")
OUTPUT_SUMMARY = os.path.join(OUTPUT_DIR, "SUMMARY_Striatum_EP.txt")

print("="*80)
print("EXTRACT STRIATUM_EP CREs ASSOCIATED WITH SPLICING GENES")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# 1. Load Splicing Genes
print("\nLoading splicing genes...")
splicing_genes_df = pd.read_csv(SPLICING_GENES_FILE)
splicing_gene_symbols = set(splicing_genes_df['gene_symbol'].dropna())
# Clean up symbols (remove brackets if any)
splicing_gene_symbols = {g.split('[')[0] for g in splicing_gene_symbols}
# Create uppercase set for case-insensitive matching
splicing_gene_symbols_upper = {g.upper() for g in splicing_gene_symbols}

print(f"Loaded {len(splicing_gene_symbols)} unique splicing genes.")

# 2. Process Striatum_EP.txt
print("\nProcessing Striatum_EP.txt...")

extracted_data = []
total_lines = 0
matched_lines = 0

with open(INPUT_FILE, 'r') as f:
    for line in f:
        total_lines += 1
        line = line.strip()
        if not line:
            continue
            
        # Split score from ID
        parts = line.split()
        if len(parts) < 2:
            continue
            
        full_id = parts[0]
        score = parts[1]
        
        # Parse ID: chr1:4775480-4775650_ENSMUSG00000025902$Sox17$chr1$4486494
        # 1. Region: chr1:4775480-4775650
        # 2. Gene Info: ENSMUSG00000025902$Sox17$chr1$4486494
        
        try:
            region_part, gene_info_part = full_id.split('_', 1)
            
            # Parse Region
            chrom, coords = region_part.split(':')
            start, end = coords.split('-')
            
            # Parse Gene Symbol (between first and second $)
            # ENSMUSG00000025902$Sox17$chr1$4486494
            gene_info_parts = gene_info_part.split('$')
            if len(gene_info_parts) >= 2:
                gene_symbol = gene_info_parts[1]
                
                # Check if gene is in splicing list
                if gene_symbol.upper() in splicing_gene_symbols_upper:
                    matched_lines += 1
                    extracted_data.append({
                        'chrom': chrom,
                        'start': int(start),
                        'end': int(end),
                        'gene_symbol': gene_symbol,
                        'score': score,
                        'full_id': full_id
                    })
                    
        except ValueError as e:
            # print(f"Skipping malformed line: {line}")
            continue

print(f"Processed {total_lines:,} lines.")
print(f"Found {matched_lines:,} CREs linked to splicing genes.")

# 3. Save Results
if not extracted_data:
    print("No matches found!")
    exit(1)

df = pd.DataFrame(extracted_data)

# Save TSV
df.to_csv(OUTPUT_TSV, sep='\t', index=False)
print(f"\nSaved TSV: {OUTPUT_TSV}")

# Save BED
# BED4: chr, start, end, name
# We'll use full_id as name to preserve all info, or construct a unique name
# Let's use: GeneSymbol_Region as name for readability, but full_id is safer for uniqueness if multiple genes link to same region?
# Actually, region+gene makes it unique.
df['bed_name'] = df['gene_symbol'] + "_" + df['chrom'] + ":" + df['start'].astype(str) + "-" + df['end'].astype(str)
df[['chrom', 'start', 'end', 'bed_name', 'score']].to_csv(OUTPUT_BED, sep='\t', header=False, index=False)
print(f"Saved BED: {OUTPUT_BED}")

# 4. Summary
unique_genes = df['gene_symbol'].nunique()
unique_regions = df.groupby(['chrom', 'start', 'end']).ngroups

with open(OUTPUT_SUMMARY, 'w') as f:
    f.write(f"Analysis Date: {datetime.now()}\n")
    f.write(f"Input File: {INPUT_FILE}\n")
    f.write(f"Total Lines Processed: {total_lines}\n")
    f.write(f"Splicing Gene Matches: {matched_lines}\n")
    f.write(f"Unique Splicing Genes: {unique_genes}\n")
    f.write(f"Unique Genomic Regions: {unique_regions}\n")

print(f"\nSummary:")
print(f"  Unique Genes: {unique_genes}")
print(f"  Unique Regions: {unique_regions}")
print("Done.")
