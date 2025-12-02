#!/usr/bin/env python3
"""
Extract ENCODE cCREs Associated with Splicing Genes (PROXIMITY VERSION)

This script identifies ENCODE cCREs (from mm10-cCREs.bed) that are associated with
splicing-related genes, using GENOMIC PROXIMITY (bedtools window).

This is an adaptation of the encode_cCREs pipeline specifically for splicing genes.

METHOD:
- Fetches gene coordinates from Ensembl BioMart (mm10/GRCm38)
- Links cCREs to genes within a specified genomic window
- Uses bedtools window for proximity-based linkage

INPUT FILES:
1. Splicing genes list:
   /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/splicing_genes/extracted_genes_final.csv

2. ENCODE cCREs:
   ../data/mm10-cCREs.bed

OUTPUT FILES:
1. CREs_splicing_genes_encode_all.tsv - All cCRE-gene associations
2. CREs_splicing_genes_encode_GABA.tsv - Associations (same, marked for GABA analysis)
3. CREs_splicing_genes_encode_by_type.tsv - Associations with CRE type information
4. SUMMARY_CREs_splicing_genes_encode.txt - Summary statistics

Usage:
    python 1_extract_encode_cCREs.py
"""

import pandas as pd
import numpy as np
import os
import subprocess
import urllib.request
import urllib.parse
import urllib.error
import time
from datetime import datetime

os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_genes_encode")

print("="*80)
print("EXTRACT ENCODE cCREs ASSOCIATED WITH SPLICING GENES (PROXIMITY)")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("", flush=True)

# =============================================================================
# Configuration
# =============================================================================

# Input files
SPLICING_GENES_FILE = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/splicing_genes/extracted_genes_final.csv"
ENCODE_CCRES_FILE = "../data/mm10-cCREs.bed"  # ENCODE cCREs

# Output directory
OUTPUT_DIR = "./output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Output files
OUTPUT_ALL = os.path.join(OUTPUT_DIR, "CREs_splicing_genes_encode_all.tsv")
OUTPUT_GABA = os.path.join(OUTPUT_DIR, "CREs_splicing_genes_encode_GABA.tsv")
OUTPUT_BY_TYPE = os.path.join(OUTPUT_DIR, "CREs_splicing_genes_encode_by_type.tsv")
OUTPUT_SUMMARY = os.path.join(OUTPUT_DIR, "SUMMARY_CREs_splicing_genes_encode.txt")

# BioMart Configuration (mm10 / GRCm38)
BIOMART_URL = "http://nov2020.archive.ensembl.org/biomart/martservice"
WINDOW_SIZE = 500000  # 500kb window for linking

# =============================================================================
# Step 1: Load Splicing Genes List
# =============================================================================

print("="*80)
print("STEP 1: Loading splicing genes list")
print("-"*80, flush=True)

CREs_splicing_genes_paper = pd.read_csv(SPLICING_GENES_FILE)
print(f"Loaded {len(CREs_splicing_genes_paper)} splicing genes from Reactome/GO")

# Extract unique gene symbols
splicing_gene_symbols = sorted(list(set(CREs_splicing_genes_paper['gene_symbol'].dropna())))
splicing_gene_symbols = [g for g in splicing_gene_symbols if not g.startswith('[')]

print(f"Unique gene symbols: {len(splicing_gene_symbols)}")
print(f"Sample genes: {', '.join(splicing_gene_symbols[:10])}")
print("", flush=True)

# =============================================================================
# Step 2: Fetch Gene Coordinates from BioMart (mm10)
# =============================================================================

print("="*80)
print("STEP 2: Fetching gene coordinates from BioMart (mm10)")
print("-"*80, flush=True)

def fetch_biomart_coords(gene_list):
    """Fetch coordinates for a list of genes from Ensembl BioMart (mm10) using urllib."""
    # Split into chunks to avoid URL length limits
    chunk_size = 200
    all_results = []

    print(f"Querying BioMart in chunks of {chunk_size} genes...")

    for i in range(0, len(gene_list), chunk_size):
        chunk = gene_list[i:i+chunk_size]
        print(f"  Fetching chunk {i//chunk_size + 1} ({len(chunk)} genes)...", flush=True)

        query_xml = f"""
        <Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">
            <Dataset name="mmusculus_gene_ensembl" interface="default">
                <Filter name="external_gene_name" value="{','.join(chunk)}"/>
                <Attribute name="external_gene_name"/>
                <Attribute name="chromosome_name"/>
                <Attribute name="start_position"/>
                <Attribute name="end_position"/>
                <Attribute name="strand"/>
            </Dataset>
        </Query>
        """

        params = urllib.parse.urlencode({'query': query_xml})
        url = f"{BIOMART_URL}?{params}"

        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                if response.status == 200:
                    text = response.read().decode('utf-8')
                    lines = text.strip().split('\n')
                    for line in lines:
                        if line.strip():
                            parts = line.split('\t')
                            if len(parts) >= 4:
                                all_results.append({
                                    'gene_symbol': parts[0],
                                    'chr': parts[1],
                                    'start': int(parts[2]),
                                    'end': int(parts[3]),
                                    'strand': parts[4] if len(parts) > 4 else '.'
                                })
                else:
                    print(f"    Error: Status {response.status}")
        except urllib.error.URLError as e:
            print(f"    URLError: {e}")
        except Exception as e:
            print(f"    Exception: {e}")

        time.sleep(1)  # Be nice to the server

    return pd.DataFrame(all_results)

# Fetch coordinates
gene_coords = fetch_biomart_coords(splicing_gene_symbols)

if len(gene_coords) == 0:
    print("ERROR: Could not fetch any gene coordinates!")
    exit(1)

# Filter for standard chromosomes
valid_chroms = [str(i) for i in range(1, 20)] + ['X', 'Y', 'MT']
gene_coords = gene_coords[gene_coords['chr'].isin(valid_chroms)].copy()

# Add 'chr' prefix if needed (Ensembl usually returns just numbers, UCSC/BED needs 'chr')
gene_coords['chr'] = 'chr' + gene_coords['chr']

print(f"Fetched coordinates for {gene_coords['gene_symbol'].nunique()} genes")
print(f"Total gene entries (including isoforms/duplicates): {len(gene_coords)}")
print("", flush=True)

# =============================================================================
# Step 3: Create BED file for Splicing Genes
# =============================================================================

print("="*80)
print("STEP 3: Creating temporary BED file for splicing genes")
print("-"*80, flush=True)

genes_bed = os.path.join(OUTPUT_DIR, "temp_splicing_genes.bed")
gene_coords = gene_coords.sort_values(['chr', 'start'])

with open(genes_bed, 'w') as f:
    for _, row in gene_coords.iterrows():
        # BED format: chr, start, end, name
        f.write(f"{row['chr']}\t{row['start']}\t{row['end']}\t{row['gene_symbol']}\n")

print(f"Created BED file: {genes_bed}")
print("", flush=True)

# =============================================================================
# Step 4: Find Closest cCREs (bedtools window)
# =============================================================================

print("="*80)
print("STEP 4: Finding ENCODE cCREs near splicing genes (bedtools window)")
print("-"*80, flush=True)

# Check if bedtools is available
try:
    subprocess.run(["bedtools", "--version"], check=True, capture_output=True)
except (FileNotFoundError, subprocess.CalledProcessError):
    print("ERROR: bedtools not found. Please load bedtools module or activate environment.")
    exit(1)

print(f"Using bedtools window with size +/- {WINDOW_SIZE} bp")

cmd = [
    "bedtools", "window",
    "-a", ENCODE_CCRES_FILE,  # cCREs (we want to keep these rows)
    "-b", genes_bed,          # Genes
    "-w", str(WINDOW_SIZE)
]

print(f"Command: {' '.join(cmd)}", flush=True)

overlaps = []
try:
    # Run bedtools
    # Output format:
    # A: chr, start, end, id1, id2, type
    # B: chr, start, end, gene_symbol
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    lines = result.stdout.strip().split('\n')
    print(f"Parsing {len(lines):,} associations...", flush=True)

    for line in lines:
        if not line: continue
        parts = line.split('\t')
        if len(parts) >= 10:
            # cCRE (A)
            ccre_chr = parts[0]
            ccre_start = int(parts[1])
            ccre_end = int(parts[2])
            ccre_id1 = parts[3]
            ccre_id2 = parts[4]
            cre_type = parts[5]

            # Gene (B) - parts[6] starts the B file
            # B file: chr, start, end, gene
            gene_symbol = parts[9]
            gene_start = int(parts[7])
            gene_end = int(parts[8])

            # Calculate distance (approximate)
            dist = 0
            if ccre_end < gene_start:
                dist = gene_start - ccre_end
            elif gene_end < ccre_start:
                dist = ccre_start - gene_end
            else:
                dist = 0 # Overlap

            overlaps.append({
                'chr': ccre_chr,
                'start': ccre_start,
                'end': ccre_end,
                'cCRE_id1': ccre_id1,
                'cCRE_id2': ccre_id2,
                'cre_type': cre_type,
                'Gene': gene_symbol,
                'Distance': dist
            })

except subprocess.CalledProcessError as e:
    print(f"ERROR running bedtools: {e.stderr}")
    exit(1)

print(f"Found {len(overlaps):,} cCRE-gene associations")

# Clean up
os.unlink(genes_bed)
print("Cleaned up temporary files", flush=True)
print("", flush=True)

# =============================================================================
# Step 5: Process and Save Results
# =============================================================================

print("="*80)
print("STEP 5: Processing and saving results")
print("-"*80, flush=True)

merged = pd.DataFrame(overlaps)

if len(merged) == 0:
    print("No associations found.")
    exit(0)

# Add metadata columns
merged['SubType'] = 'Proximity_Linked'
merged['is_gaba'] = True  # All will be analyzed with GABA BigWig files

# Save files
merged.to_csv(OUTPUT_ALL, sep='\t', index=False)
merged.to_csv(OUTPUT_GABA, sep='\t', index=False)
merged.to_csv(OUTPUT_BY_TYPE, sep='\t', index=False)

print(f"Saved: {OUTPUT_ALL}")
print(f"Saved: {OUTPUT_GABA}")
print(f"Saved: {OUTPUT_BY_TYPE}")
print("", flush=True)

# =============================================================================
# Step 6: Summary Statistics
# =============================================================================

print("="*80)
print("STEP 6: Generating summary statistics")
print("-"*80, flush=True)

# Statistics
unique_cres = merged['cCRE_id1'].nunique()
unique_genes = merged['Gene'].nunique()
total_links = len(merged)

print(f"\nSUMMARY STATISTICS:")
print(f"  Input splicing genes: {len(splicing_gene_symbols)}")
print(f"  Genes with coordinates: {gene_coords['gene_symbol'].nunique()}")
print(f"  Linked cCREs: {unique_cres:,}")
print(f"  Linked genes: {unique_genes}")
print(f"  Total associations: {total_links:,}")

# CRE type distribution
print(f"\nCRE TYPE DISTRIBUTION:")
for cre_type, count in merged['cre_type'].value_counts().items():
    print(f"  {cre_type}: {count:,}")

# Top genes by CRE count
print(f"\nTOP 20 GENES BY CRE COUNT:")
top_genes = merged.groupby('Gene')['cCRE_id1'].nunique().sort_values(ascending=False).head(20)
for gene, count in top_genes.items():
    print(f"  {gene}: {count:,} cCREs")

# Distance distribution
print(f"\nDISTANCE DISTRIBUTION:")
print(f"  Min: {merged['Distance'].min():,} bp")
print(f"  Median: {merged['Distance'].median():,.0f} bp")
print(f"  Mean: {merged['Distance'].mean():,.0f} bp")
print(f"  Max: {merged['Distance'].max():,} bp")

print("", flush=True)

# =============================================================================
# Step 7: Save Summary Report
# =============================================================================

print("="*80)
print("STEP 7: Saving summary report")
print("-"*80, flush=True)

with open(OUTPUT_SUMMARY, 'w') as f:
    f.write("="*80 + "\n")
    f.write("ENCODE cCREs ASSOCIATED WITH SPLICING GENES (PROXIMITY)\n")
    f.write("="*80 + "\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("METHOD:\n")
    f.write("-"*80 + "\n")
    f.write(f"  Genomic Proximity (Window +/- {WINDOW_SIZE:,} bp)\n")
    f.write(f"  Gene coordinates: Ensembl BioMart (mm10/GRCm38)\n")
    f.write(f"  cCRE source: ENCODE mm10-cCREs.bed\n\n")

    f.write("INPUT FILES:\n")
    f.write("-"*80 + "\n")
    f.write(f"  Splicing genes: {SPLICING_GENES_FILE}\n")
    f.write(f"  ENCODE cCREs: {ENCODE_CCRES_FILE}\n\n")

    f.write("RESULTS:\n")
    f.write("-"*80 + "\n")
    f.write(f"  Input splicing genes: {len(splicing_gene_symbols)}\n")
    f.write(f"  Genes with coordinates: {gene_coords['gene_symbol'].nunique()}\n")
    f.write(f"  Linked cCREs: {unique_cres:,}\n")
    f.write(f"  Linked genes: {unique_genes}\n")
    f.write(f"  Total associations: {total_links:,}\n\n")

    f.write("CRE TYPE DISTRIBUTION:\n")
    f.write("-"*80 + "\n")
    for cre_type, count in merged['cre_type'].value_counts().items():
        f.write(f"  {cre_type}: {count:,}\n")

    f.write("\nDISTANCE DISTRIBUTION:\n")
    f.write("-"*80 + "\n")
    f.write(f"  Min: {merged['Distance'].min():,} bp\n")
    f.write(f"  Median: {merged['Distance'].median():,.0f} bp\n")
    f.write(f"  Mean: {merged['Distance'].mean():,.0f} bp\n")
    f.write(f"  Max: {merged['Distance'].max():,} bp\n")

    f.write("\nTOP 20 GENES BY CRE COUNT:\n")
    f.write("-"*80 + "\n")
    for gene, count in top_genes.items():
        f.write(f"  {gene}: {count:,} cCREs\n")

    f.write("\n" + "="*80 + "\n")
    f.write("OUTPUT FILES:\n")
    f.write("="*80 + "\n")
    f.write(f"  1. {OUTPUT_ALL}\n")
    f.write(f"  2. {OUTPUT_GABA}\n")
    f.write(f"  3. {OUTPUT_BY_TYPE}\n")
    f.write(f"  4. {OUTPUT_SUMMARY}\n")

print(f"Saved summary: {OUTPUT_SUMMARY}")

# =============================================================================
# Completion
# =============================================================================

print("\n" + "="*80)
print("ANALYSIS COMPLETE!")
print("="*80)
print()
print(f"Output directory: {OUTPUT_DIR}/")
print()
print("Generated files:")
print(f"  1. CREs_splicing_genes_encode_all.tsv ({total_links:,} associations)")
print(f"  2. CREs_splicing_genes_encode_GABA.tsv ({total_links:,} associations)")
print(f"  3. CREs_splicing_genes_encode_by_type.tsv ({total_links:,} associations)")
print(f"  4. SUMMARY_CREs_splicing_genes_encode.txt (summary report)")
print()
print("Key Results:")
print("-"*80)
print(f"  Splicing genes analyzed: {len(splicing_gene_symbols)}")
print(f"  Genes with ENCODE cCREs: {unique_genes} ({100*unique_genes/len(splicing_gene_symbols):.1f}%)")
print(f"  Total cCREs linked: {unique_cres:,}")
print(f"  Total associations: {total_links:,}")
print()
print(f"DATA SOURCE: ENCODE cCREs (mm10-cCREs.bed)")
print(f"LINKAGE: Genomic proximity (+/- {WINDOW_SIZE:,} bp)")
print()
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
