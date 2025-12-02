#!/usr/bin/env python3
"""
Step 1: Extract Splicing Gene CREs with PCC Directionality

This script extracts CREs linked to splicing genes from Table 16,
SEPARATING by correlation direction:
- Enhancer-like CREs (PCC > 0): accessibility ↑ = expression ↑
- Silencer-like CREs (PCC < 0): accessibility ↑ = expression ↓

This is critical for proper biological interpretation:
- Loss of accessibility at ENHANCERS → reduced gene expression
- Loss of accessibility at SILENCERS → increased gene expression

For splicing disturbances (likely reduced splicing gene expression),
we expect to see loss of accessibility at ENHANCER-like CREs.

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
OUTPUT_DIR = BASE_DIR / "CREs_splicing_directional" / "output"

# Input files
SPLICING_GENES_FILE = Path("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/CREs_splicing_genes_paper/extracted_genes_final.csv")
TABLE16_FILE = DATA_DIR / "table_16.txt"
ENCODE_CCRES_FILE = DATA_DIR / "mm10-cCREs.bed"

# Statistical thresholds
FDR_THRESHOLD = 0.05
PCC_THRESHOLD = 0.2  # Will be applied directionally: PCC > 0.2 OR PCC < -0.2

# GABA/Hippocampal cell type keywords (for SubType filtering)
GABA_KEYWORDS = [
    # Hippocampal excitatory (for comparison)
    'CA1', 'CA2', 'CA3', 'CA4', 'DG', 'DGNBL', 'GRC',
    # GABAergic interneurons
    'GABA', 'PV', 'SST', 'VIP', 'LAMP5', 'LAMP',
    'PVGA', 'SSTGA', 'VIPGA', 'LAMGA', 'INH',
    # Additional subtypes
    'CGE', 'MGE', 'LGE'
]

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


def is_gaba_subtype(subtype):
    """Check if a SubType matches GABA/hippocampal keywords."""
    if pd.isna(subtype):
        return False
    subtype_upper = str(subtype).upper()
    return any(kw.upper() in subtype_upper for kw in GABA_KEYWORDS)


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


def filter_gaba_subtypes(df, cre_type):
    """Filter for GABA/hippocampal cell types."""
    if len(df) == 0:
        return df, df

    print(f"\nFiltering {cre_type} CREs for GABA/hippocampal subtypes...")

    # All cell types
    all_celltypes = df.copy()

    # GABA-specific
    gaba_mask = df['SubType'].apply(is_gaba_subtype)
    gaba_specific = df[gaba_mask].copy()

    print(f"  All cell types: {len(all_celltypes):,} links")
    print(f"  GABA/hippocampal: {len(gaba_specific):,} links")

    return all_celltypes, gaba_specific


def extract_bed_coordinates(df, cre_type):
    """Extract unique CRE coordinates in BED format."""
    if len(df) == 0:
        return pd.DataFrame()

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


def intersect_with_encode(bed_file, output_file):
    """Intersect with ENCODE cCREs to get type annotations."""
    print(f"\nIntersecting with ENCODE cCREs...")

    # Run bedtools intersect
    cmd = f"bedtools intersect -a {bed_file} -b {ENCODE_CCRES_FILE} -wa -wb"

    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)

        # Parse output
        lines = result.stdout.strip().split('\n')
        if not lines or lines[0] == '':
            print("  No overlaps found with ENCODE cCREs")
            return pd.DataFrame()

        # Parse intersected results
        records = []
        for line in lines:
            parts = line.split('\t')
            if len(parts) >= 12:  # Our BED6 + ENCODE BED6
                records.append({
                    'chrom': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'cre_id': parts[3],
                    'score': float(parts[4]),
                    'strand': parts[5],
                    'encode_chrom': parts[6],
                    'encode_start': int(parts[7]),
                    'encode_end': int(parts[8]),
                    'encode_id': parts[9],
                    'encode_score': parts[10],
                    'encode_type': parts[11] if len(parts) > 11 else 'unknown'
                })

        intersect_df = pd.DataFrame(records)
        print(f"  Found {len(intersect_df):,} overlaps with ENCODE cCREs")

        # Count by ENCODE type
        if len(intersect_df) > 0:
            type_counts = intersect_df['encode_type'].value_counts()
            print("\n  ENCODE CRE type distribution:")
            for cre_type, count in type_counts.head(10).items():
                print(f"    {cre_type}: {count:,}")

        return intersect_df

    except subprocess.CalledProcessError as e:
        print(f"  bedtools error: {e.stderr}")
        return pd.DataFrame()


def save_outputs(enhancer_all, enhancer_gaba, silencer_all, silencer_gaba):
    """Save all output files."""
    print("\n" + "=" * 80)
    print("SAVING OUTPUT FILES")
    print("=" * 80)

    # Save TSV files with all information
    outputs = [
        (enhancer_all, "enhancer_CREs_all_celltypes.tsv"),
        (enhancer_gaba, "enhancer_CREs_GABA.tsv"),
        (silencer_all, "silencer_CREs_all_celltypes.tsv"),
        (silencer_gaba, "silencer_CREs_GABA.tsv"),
    ]

    for df, filename in outputs:
        if len(df) > 0:
            outpath = OUTPUT_DIR / filename
            df.to_csv(outpath, sep='\t', index=False)
            print(f"Saved: {filename} ({len(df):,} rows)")

    # Create BED files for deepTools
    bed_outputs = []

    for df, name, cre_type in [
        (enhancer_gaba, "enhancer_CREs_GABA", "enhancer"),
        (silencer_gaba, "silencer_CREs_GABA", "silencer"),
        (enhancer_all, "enhancer_CREs_all", "enhancer"),
        (silencer_all, "silencer_CREs_all", "silencer"),
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


def write_summary(CREs_splicing_genes_paper, enhancer_all, enhancer_gaba, silencer_all, silencer_gaba, bed_outputs):
    """Write summary statistics file."""
    print("\n" + "=" * 80)
    print("WRITING SUMMARY")
    print("=" * 80)

    summary_path = OUTPUT_DIR / "SUMMARY_splicing_CREs_directional.txt"

    with open(summary_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("SPLICING GENE CRE ANALYSIS - DIRECTIONAL\n")
        f.write("=" * 80 + "\n\n")

        f.write("RESEARCH QUESTION:\n")
        f.write("-" * 80 + "\n")
        f.write("Identify loss of chromatin accessibility at splicing gene regulatory elements\n")
        f.write("specifically in Nestin-MUT GABA neurons compared to:\n")
        f.write("  1. Nestin-Ctrl (within-genotype mutation effect)\n")
        f.write("  2. Emx1-Mut (genotype-specific response)\n\n")

        f.write("KEY INSIGHT - PCC DIRECTIONALITY:\n")
        f.write("-" * 80 + "\n")
        f.write("ENHANCER-like CREs (PCC > 0): accessibility UP = expression UP\n")
        f.write("  -> Loss of accessibility at enhancers -> REDUCED gene expression\n")
        f.write("  -> Relevant if splicing genes are DOWN-regulated in Nestin-MUT\n\n")
        f.write("SILENCER-like CREs (PCC < 0): accessibility UP = expression DOWN\n")
        f.write("  -> Loss of accessibility at silencers -> INCREASED gene expression\n")
        f.write("  -> Relevant if splicing genes are UP-regulated in Nestin-MUT\n\n")

        f.write("INPUT DATA:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Splicing genes: {len(CREs_splicing_genes_paper):,}\n")
        f.write(f"Table 16 (CRE-gene correlations): {TABLE16_FILE}\n")
        f.write(f"ENCODE cCREs: {ENCODE_CCRES_FILE}\n\n")

        f.write("STATISTICAL FILTERS:\n")
        f.write("-" * 80 + "\n")
        f.write(f"FDR threshold: < {FDR_THRESHOLD}\n")
        f.write(f"PCC threshold: > {PCC_THRESHOLD} (enhancers) OR < -{PCC_THRESHOLD} (silencers)\n")
        f.write(f"Cell type filter: GABA/Hippocampal keywords\n\n")

        f.write("RESULTS:\n")
        f.write("-" * 80 + "\n")
        f.write("\nENHANCER-like CREs (PCC > 0.2):\n")
        f.write(f"  All cell types: {len(enhancer_all):,} CRE-gene links\n")
        f.write(f"  GABA/hippocampal: {len(enhancer_gaba):,} CRE-gene links\n")

        f.write("\nSILENCER-like CREs (PCC < -0.2):\n")
        f.write(f"  All cell types: {len(silencer_all):,} CRE-gene links\n")
        f.write(f"  GABA/hippocampal: {len(silencer_gaba):,} CRE-gene links\n")

        f.write("\nBED FILES GENERATED:\n")
        for bed_path, name, count in bed_outputs:
            f.write(f"  {name}.bed: {count:,} unique CREs\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("BIOLOGICAL INTERPRETATION:\n")
        f.write("=" * 80 + "\n\n")

        f.write("For splicing disturbances in Nestin-MUT GABA neurons:\n\n")
        f.write("1. If splicing genes show REDUCED expression:\n")
        f.write("   -> Focus on ENHANCER-like CREs\n")
        f.write("   -> Look for LOSS of accessibility in Nestin-MUT\n")
        f.write("   -> Pattern: Nestin-Ctrl HIGH, Nestin-Mut LOW, Emx1-Mut HIGH/NORMAL\n\n")

        f.write("2. If splicing genes show INCREASED expression:\n")
        f.write("   -> Focus on SILENCER-like CREs\n")
        f.write("   -> Look for LOSS of accessibility in Nestin-MUT\n")
        f.write("   -> Pattern: Nestin-Ctrl HIGH, Nestin-Mut LOW, Emx1-Mut HIGH/NORMAL\n\n")

        f.write("NEXT STEPS:\n")
        f.write("-" * 80 + "\n")
        f.write("1. Run Step 2: Compute signal matrices (2_compute_signal_matrices.sh)\n")
        f.write("2. Run Step 3: Visualize comparisons and identify Nestin-specific loss\n")
        f.write("   (3_visualize_directional_comparisons.py)\n")

    print(f"Saved: {summary_path}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    print("\n" + "=" * 80)
    print("SPLICING GENE CRE ANALYSIS - WITH DIRECTIONALITY")
    print("=" * 80)
    print("\nThis analysis separates CREs by correlation direction to enable")
    print("proper biological interpretation of accessibility changes.\n")

    # Load splicing genes
    CREs_splicing_genes_paper = load_splicing_genes()

    # Process Table 16 with PCC directionality
    enhancer_df, silencer_df = process_table16_chunks(CREs_splicing_genes_paper)

    # Filter for GABA subtypes
    enhancer_all, enhancer_gaba = filter_gaba_subtypes(enhancer_df, "enhancer")
    silencer_all, silencer_gaba = filter_gaba_subtypes(silencer_df, "silencer")

    # Save outputs
    bed_outputs = save_outputs(enhancer_all, enhancer_gaba, silencer_all, silencer_gaba)

    # Write summary
    write_summary(CREs_splicing_genes_paper, enhancer_all, enhancer_gaba, silencer_all, silencer_gaba, bed_outputs)

    print("\n" + "=" * 80)
    print("EXTRACTION COMPLETE")
    print("=" * 80)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print("\nKey files for downstream analysis:")
    print("  - enhancer_CREs_GABA.bed (for loss of accessibility -> reduced expression)")
    print("  - silencer_CREs_GABA.bed (for loss of accessibility -> increased expression)")
    print("\nNext: Run 2_compute_signal_matrices.sh")


if __name__ == "__main__":
    main()
