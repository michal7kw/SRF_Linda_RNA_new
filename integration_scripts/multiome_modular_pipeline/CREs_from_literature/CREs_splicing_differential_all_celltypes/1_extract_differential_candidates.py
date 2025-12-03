#!/usr/bin/env python3
"""
Step 1: Extract Splicing Gene CREs (All Cell Types) with PCC Directionality

This script extracts CREs linked to splicing genes from Table 16,
SEPARATING by correlation direction:
- Enhancer-like CREs (PCC > 0): accessibility ↑ = expression ↑
- Silencer-like CREs (PCC < 0): accessibility ↑ = expression ↓

DIFFERENCE FROM PREVIOUS PIPELINES:
- NO filtering for GABA/Hippocampal cell types.
- Uses ALL available CRE-gene links for splicing genes.

Author: Claude Code
Date: December 2024
"""

import pandas as pd
import numpy as np
import os
import subprocess
from pathlib import Path
from collections import defaultdict

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths
BASE_DIR = Path("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature")
DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = BASE_DIR / "CREs_splicing_differential_all_celltypes" / "output"

# Input files
SPLICING_GENES_FILE = Path("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/splicing_genes/extracted_genes_final.csv")
TABLE16_FILE = DATA_DIR / "table_16.txt"
ENCODE_CCRES_FILE = DATA_DIR / "mm10-cCREs.bed"

# Statistical thresholds
FDR_THRESHOLD = 0.05
PCC_THRESHOLD = 0.2  # Will be applied directionally: PCC > 0.2 OR PCC < -0.2

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_splicing_genes():
    """Load and process splicing genes list."""
    print("=" * 80)
    print("LOADING SPLICING GENES")
    print("=" * 80)

    df = pd.read_csv(SPLICING_GENES_FILE)
    print(f"Loaded {len(df)} entries from splicing genes file")

    # Get unique gene symbols (convert to uppercase for matching)
    genes = df['gene_symbol'].str.upper().unique().tolist()
    print(f"Unique splicing genes: {len(genes)}")

    return set(genes)


def parse_coordinate(coord_str):
    """Parse coordinate string (chr_start_end) to BED format."""
    try:
        parts = coord_str.replace(':', '_').replace('-', '_').split('_')
        if len(parts) >= 3:
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            return chrom, start, end
    except:
        pass
    return None, None, None


def process_table16_chunks(CREs_splicing_genes_paper, chunk_size=500000):
    """Process Table 16 in chunks to handle large file."""
    print("\n" + "=" * 80)
    print("PROCESSING TABLE 16 (CRE-Gene Correlations)")
    print("=" * 80)

    # Initialize result containers
    enhancer_links = []  # PCC > threshold (positive correlation)
    silencer_links = []  # PCC < -threshold (negative correlation)

    total_rows = 0
    splicing_matches = 0

    # Process in chunks
    print(f"Processing Table 16 in chunks of {chunk_size:,} rows...")

    for chunk_num, chunk in enumerate(pd.read_csv(TABLE16_FILE, sep='\t', chunksize=chunk_size)):
        total_rows += len(chunk)

        # Filter for splicing genes (case-insensitive)
        chunk['Gene_upper'] = chunk['Gene'].str.upper()
        splicing_chunk = chunk[chunk['Gene_upper'].isin(CREs_splicing_genes_paper)].copy()

        if len(splicing_chunk) == 0:
            continue

        splicing_matches += len(splicing_chunk)

        # Apply FDR filter
        splicing_chunk = splicing_chunk[splicing_chunk['FDR'] < FDR_THRESHOLD]

        if len(splicing_chunk) == 0:
            continue

        # Separate by PCC direction
        # ENHANCERS: positive correlation (PCC > threshold)
        enhancers = splicing_chunk[splicing_chunk['PCC'] > PCC_THRESHOLD].copy()
        enhancers['CRE_type'] = 'enhancer'

        # SILENCERS: negative correlation (PCC < -threshold)
        silencers = splicing_chunk[splicing_chunk['PCC'] < -PCC_THRESHOLD].copy()
        silencers['CRE_type'] = 'silencer'

        if len(enhancers) > 0:
            enhancer_links.append(enhancers)
        if len(silencers) > 0:
            silencer_links.append(silencers)

        if (chunk_num + 1) % 5 == 0:
            print(f"  Processed {total_rows:,} rows, found {splicing_matches:,} splicing gene matches...")

    print(f"\nTotal rows processed: {total_rows:,}")
    print(f"Total splicing gene matches: {splicing_matches:,}")

    # Combine results
    enhancer_df = pd.concat(enhancer_links, ignore_index=True) if enhancer_links else pd.DataFrame()
    silencer_df = pd.concat(silencer_links, ignore_index=True) if silencer_links else pd.DataFrame()

    print(f"\nAfter FDR < {FDR_THRESHOLD} and PCC filtering:")
    print(f"  Enhancer-like CRE links (PCC > {PCC_THRESHOLD}): {len(enhancer_df):,}")
    print(f"  Silencer-like CRE links (PCC < -{PCC_THRESHOLD}): {len(silencer_df):,}")

    return enhancer_df, silencer_df


def extract_bed_coordinates(df, cre_type):
    """Extract unique CRE coordinates in BED format."""
    if len(df) == 0:
        return pd.DataFrame(), pd.DataFrame()

    # Parse coordinates
    coords = []
    for _, row in df.iterrows():
        chrom, start, end = parse_coordinate(row['Coordinate1'])
        if chrom is not None:
            coords.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'cre_id': row.get('cCRE1', f"{chrom}_{start}_{end}"),
                'score': abs(row['PCC']) * 1000,  # Score based on correlation strength
                'strand': '+' if row['PCC'] > 0 else '-',  # Encode direction in strand
                'gene': row['Gene'],
                'pcc': row['PCC'],
                'fdr': row['FDR'],
                'subtype': row.get('SubType', 'NA'),
                'cre_type': cre_type
            })

    bed_df = pd.DataFrame(coords)

    # Get unique CREs (same CRE can link to multiple genes)
    unique_cres = bed_df.drop_duplicates(subset=['chrom', 'start', 'end'])

    print(f"  Unique {cre_type} CREs: {len(unique_cres):,}")

    return bed_df, unique_cres


def save_outputs(enhancer_all, silencer_all):
    """Save all output files."""
    print("\n" + "=" * 80)
    print("SAVING OUTPUT FILES")
    print("=" * 80)

    # Save TSV files with all information
    outputs = [
        (enhancer_all, "enhancer_CREs_all_celltypes.tsv"),
        (silencer_all, "silencer_CREs_all_celltypes.tsv"),
    ]

    for df, filename in outputs:
        if len(df) > 0:
            outpath = OUTPUT_DIR / filename
            df.to_csv(outpath, sep='\t', index=False)
            print(f"Saved: {filename} ({len(df):,} rows)")

    # Create BED files for deepTools
    bed_outputs = []

    for df, name, cre_type in [
        (enhancer_all, "enhancer_CREs_all_celltypes", "enhancer"),
        (silencer_all, "silencer_CREs_all_celltypes", "silencer"),
    ]:
        if len(df) > 0:
            bed_df, unique_cres = extract_bed_coordinates(df, cre_type)

            # Save full BED with gene info
            full_bed_path = OUTPUT_DIR / f"{name}_full.bed"
            bed_df[['chrom', 'start', 'end', 'cre_id', 'score', 'strand']].to_csv(
                full_bed_path, sep='\t', index=False, header=False
            )

            # Save unique CREs BED (for deepTools)
            unique_bed_path = OUTPUT_DIR / f"{name}.bed"
            unique_cres[['chrom', 'start', 'end', 'cre_id', 'score', 'strand']].to_csv(
                unique_bed_path, sep='\t', index=False, header=False
            )

            print(f"Saved: {name}.bed ({len(unique_cres):,} unique CREs)")

            # Save gene linkage info
            gene_info_path = OUTPUT_DIR / f"{name}_gene_links.tsv"
            bed_df.to_csv(gene_info_path, sep='\t', index=False)

            bed_outputs.append((unique_bed_path, name, len(unique_cres)))

    return bed_outputs


def write_summary(CREs_splicing_genes_paper, enhancer_all, silencer_all, bed_outputs):
    """Write summary statistics file."""
    print("\n" + "=" * 80)
    print("WRITING SUMMARY")
    print("=" * 80)

    summary_path = OUTPUT_DIR / "SUMMARY_extraction.txt"

    with open(summary_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("SPLICING GENE CRE ANALYSIS - ALL CELL TYPES\n")
        f.write("=" * 80 + "\n\n")

        f.write("INPUT DATA:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Splicing genes: {len(CREs_splicing_genes_paper):,}\n")
        f.write(f"Table 16 (CRE-gene correlations): {TABLE16_FILE}\n")
        f.write(f"ENCODE cCREs: {ENCODE_CCRES_FILE}\n\n")

        f.write("STATISTICAL FILTERS:\n")
        f.write("-" * 80 + "\n")
        f.write(f"FDR threshold: < {FDR_THRESHOLD}\n")
        f.write(f"PCC threshold: > {PCC_THRESHOLD} (enhancers) OR < -{PCC_THRESHOLD} (silencers)\n")
        f.write(f"Cell type filter: NONE (All cell types included)\n\n")

        f.write("RESULTS:\n")
        f.write("-" * 80 + "\n")
        f.write("\nENHANCER-like CREs (PCC > 0.2):\n")
        f.write(f"  All cell types: {len(enhancer_all):,} CRE-gene links\n")

        f.write("\nSILENCER-like CREs (PCC < -0.2):\n")
        f.write(f"  All cell types: {len(silencer_all):,} CRE-gene links\n")

        f.write("\nBED FILES GENERATED:\n")
        for bed_path, name, count in bed_outputs:
            f.write(f"  {name}.bed: {count:,} unique CREs\n")

    print(f"Saved: {summary_path}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    print("\n" + "=" * 80)
    print("SPLICING GENE CRE ANALYSIS - ALL CELL TYPES")
    print("=" * 80)
    print("\nThis analysis extracts ALL CREs linked to splicing genes, without")
    print("cell type restrictions, separating by correlation direction.\n")

    # Load splicing genes
    CREs_splicing_genes_paper = load_splicing_genes()

    # Process Table 16 with PCC directionality
    enhancer_df, silencer_df = process_table16_chunks(CREs_splicing_genes_paper)

    # Save outputs (No GABA filtering step)
    bed_outputs = save_outputs(enhancer_df, silencer_df)

    # Write summary
    write_summary(CREs_splicing_genes_paper, enhancer_df, silencer_df, bed_outputs)

    print("\n" + "=" * 80)
    print("EXTRACTION COMPLETE")
    print("=" * 80)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print("\nNext: Run 2_compute_signal_matrices.sh")


if __name__ == "__main__":
    main()
