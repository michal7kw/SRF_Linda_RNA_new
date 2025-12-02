#!/usr/bin/env python3
"""
Extract ENCODE cCREs Associated with Splicing Genes (PROXIMITY VERSION)

This script identifies ENCODE cCREs (from mm10-cCREs.bed) that are associated with
splicing-related genes, using GENOMIC PROXIMITY (bedtools closest).

CHANGES:
- Removed Table 16 dependency.
- Fetches gene coordinates from Ensembl BioMart (mm10/GRCm38).
- Links cCREs to the nearest splicing gene.
- Uses urllib (standard library) to avoid dependency issues.

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

os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_encode_all")

print("="*80)
print("EXTRACT ENCODE cCREs ASSOCIATED WITH SPLICING GENES (PROXIMITY)")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("", flush=True)

# =============================================================================
# Configuration
# =============================================================================

# Input files
SPLICING_GENES_FILE = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/CREs_splicing_genes_paper/extracted_genes_final.csv"
ENCODE_CCRES_FILE = "../data/mm10-cCREs.bed"  # ENCODE cCREs

# Output directory
OUTPUT_DIR = "./output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Output files
OUTPUT_ALL = os.path.join(OUTPUT_DIR, "encode_cCREs_linked.tsv")
OUTPUT_BY_TYPE = os.path.join(OUTPUT_DIR, "encode_cCREs_by_type.tsv")
OUTPUT_SUMMARY = os.path.join(OUTPUT_DIR, "SUMMARY_encode_cCREs.txt")

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
# Step 4: Find Closest cCREs (bedtools closest)
# =============================================================================

print("="*80)
print("STEP 4: Finding closest ENCODE cCREs (bedtools window)")
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

# Add "is_gaba" column (dummy for now)
merged['SubType'] = 'Proximity_Linked'
merged['is_gaba'] = True 

# Save files
# Save files
merged.to_csv(OUTPUT_ALL, sep='\t', index=False)
merged.to_csv(OUTPUT_BY_TYPE, sep='\t', index=False)

print(f"Saved: {OUTPUT_ALL}")
print(f"Saved: {OUTPUT_BY_TYPE}")
print("", flush=True)

# =============================================================================
# Step 6: Summary
# =============================================================================

print("="*80)
print("STEP 6: Generating summary")
print("-"*80, flush=True)

with open(OUTPUT_SUMMARY, 'w') as f:
    f.write("="*80 + "\n")
    f.write("ENCODE cCREs ASSOCIATED WITH SPLICING GENES (PROXIMITY)\n")
    f.write("="*80 + "\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write(f"Method: Genomic Proximity (Window +/- {WINDOW_SIZE} bp)\n")
    f.write(f"Genes Source: Ensembl BioMart (mm10)\n\n")
    
    f.write("RESULTS:\n")
    f.write(f"  Input Genes: {len(splicing_gene_symbols)}\n")
    f.write(f"  Genes with Coords: {gene_coords['gene_symbol'].nunique()}\n")
    f.write(f"  Linked cCREs: {merged['cCRE_id1'].nunique():,}\n")
    f.write(f"  Linked Genes: {merged['Gene'].nunique()}\n")
    f.write(f"  Total Associations: {len(merged):,}\n\n")
    
    f.write("CRE TYPE DISTRIBUTION:\n")
    for cre_type, count in merged['cre_type'].value_counts().items():
        f.write(f"  {cre_type}: {count:,}\n")

print(f"Saved summary: {OUTPUT_SUMMARY}")
print("="*80)
print("ANALYSIS COMPLETE")
print("="*80)
