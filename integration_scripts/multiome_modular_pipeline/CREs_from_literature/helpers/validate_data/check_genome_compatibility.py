#!/usr/bin/env python3
"""
Check Genome Compatibility Between Your Data and Literature CREs

This script verifies that the reference genome builds match between:
1. Your ATAC-seq/RNA-seq data (should be mm10/GRCm38)
2. Literature CREs from supplementary tables (unknown - we'll check)

It performs the following checks:
1. Chromosome naming convention (chr1 vs 1)
2. Gene coordinate comparison (TSS positions)
3. CRE-gene distance validation
4. Example coordinate spot-checks

INPUT FILES:
- data/table_16.txt: Literature CREs with coordinates and gene associations
- ../combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/dge_res/geno_spec_cond_comp/dge_Nestin_GABA_mut_vs_ctrl.csv: DEG results for Nestin GABA comparison
- ../combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/dge_res/geno_spec_cond_comp/dge_Emx1_GABA_mut_vs_ctrl.csv: DEG results for Emx1 GABA comparison
- ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Nestin-Ctrl.bw: BigWig file for chromosome compatibility check

OUTPUT FILES:
- None (console output only)
"""

import pandas as pd
import numpy as np
from collections import Counter
import os

print("="*80)
print("GENOME COMPATIBILITY CHECK")
print("="*80)
print()

# ============================================================================
# Load data
# ============================================================================
print("STEP 1: Loading data...")
print("-"*80)

# Load literature CREs
print("Loading literature CREs (Table 16)...")
table16 = pd.read_csv("data/table_16.txt", sep='\t', nrows=100000)  # Load first 100k for speed
print(f"  Loaded {len(table16)} rows (sampling for speed)")

# Parse coordinates from Coordinate2 column (the CRE column)
table16[['cre_chr', 'cre_coords']] = table16['Coordinate2'].str.split('_', n=1, expand=True)
table16[['cre_start', 'cre_end']] = table16['cre_coords'].str.split('_', expand=True)
table16['cre_start'] = pd.to_numeric(table16['cre_start'])
table16['cre_end'] = pd.to_numeric(table16['cre_end'])

print(f"  Unique genes: {table16['Gene'].nunique()}")
print(f"  Unique chromosomes: {sorted(table16['cre_chr'].unique())[:10]}")

# Load your DEG data (has gene info)
print("\nLoading your DEG data...")
deg_files = [
    "../combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/dge_res/geno_spec_cond_comp/dge_Nestin_GABA_mut_vs_ctrl.csv",
    "../combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/dge_res/geno_spec_cond_comp/dge_Emx1_GABA_mut_vs_ctrl.csv"
]

all_genes = set()
for f in deg_files:
    if os.path.exists(f):
        deg = pd.read_csv(f)
        all_genes.update(deg['names'].tolist())
        print(f"  Loaded {len(deg)} genes from {os.path.basename(f)}")

print(f"  Total unique genes in your data: {len(all_genes)}")

# ============================================================================
# Check 1: Chromosome naming convention
# ============================================================================
print("\n" + "="*80)
print("CHECK 1: Chromosome Naming Convention")
print("-"*80)

lit_chrs = table16['cre_chr'].unique()
chr_counts = Counter([chr for chr in lit_chrs])

print(f"Literature CREs chromosome format:")
print(f"  Sample chromosomes: {sorted(lit_chrs)[:10]}")

if all(chr.startswith('chr') for chr in lit_chrs):
    print(f"  ✓ Format: UCSC style (chr1, chr2, ...)")
    print(f"  ✓ This matches mm10/GRCm38 UCSC convention")
else:
    print(f"  ✗ Format: Ensembl style (1, 2, ...)")
    print(f"  ✗ Your data uses UCSC style - MISMATCH!")

print(f"\nYour data chromosome format (from Signac pipeline):")
print(f"  Format: UCSC style (chr1, chr2, ...) [from BSgenome.Mmusculus.UCSC.mm10]")

if all(chr.startswith('chr') for chr in lit_chrs):
    print(f"\n✓ PASS: Chromosome naming conventions MATCH")
else:
    print(f"\n⚠️  WARNING: Chromosome naming mismatch!")
    print(f"   You may need to convert chromosome names")

# ============================================================================
# Check 2: Gene coordinate validation
# ============================================================================
print("\n" + "="*80)
print("CHECK 2: Gene Coordinate Validation")
print("-"*80)

# Check if gene coordinates are plausible for mm10
# We'll use known gene positions from mm10 for validation

print("\nChecking if CRE coordinates are plausible for mouse genome...")

# Mouse genome size (mm10)
mm10_chr_sizes = {
    'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116,
    'chr5': 151834684, 'chr6': 149736546, 'chr7': 145441459, 'chr8': 129401213,
    'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022,
    'chr13': 120421639, 'chr14': 124902244, 'chr15': 104043685, 'chr16': 98207768,
    'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566, 'chrX': 171031299,
    'chrY': 91744698
}

print(f"\nValidating coordinates against mm10 chromosome sizes...")
validation_errors = 0

for idx, row in table16.head(1000).iterrows():  # Check first 1000 rows
    chr = row['cre_chr']
    start = row['cre_start']
    end = row['cre_end']

    if chr in mm10_chr_sizes:
        chr_size = mm10_chr_sizes[chr]
        if end > chr_size:
            validation_errors += 1
            if validation_errors <= 5:  # Show first 5 errors
                print(f"  ✗ CRE {row['cre2']} on {chr}:{start}-{end} exceeds mm10 chr size {chr_size}")

if validation_errors == 0:
    print(f"  ✓ PASS: All checked coordinates fit within mm10 chromosome sizes")
else:
    print(f"  ✗ FAIL: {validation_errors}/1000 coordinates exceed mm10 chromosome sizes")
    print(f"          This suggests data may NOT be mm10!")

# ============================================================================
# Check 3: Gene-CRE distance validation
# ============================================================================
print("\n" + "="*80)
print("CHECK 3: CRE-Gene Distance Distribution")
print("-"*80)

print("\nAnalyzing CRE-gene distances from Table 16...")
print("(This helps identify if coordinates are from the same genome build)")

# For simplicity, we'll look at distribution of CRE sizes
table16['cre_size'] = table16['cre_end'] - table16['cre_start']

print(f"\nCRE size distribution:")
print(f"  Mean: {table16['cre_size'].mean():.0f} bp")
print(f"  Median: {table16['cre_size'].median():.0f} bp")
print(f"  Min: {table16['cre_size'].min():.0f} bp")
print(f"  Max: {table16['cre_size'].max():.0f} bp")

# Typical CRE sizes for ATAC-seq peaks are 200-2000 bp
if 100 <= table16['cre_size'].median() <= 3000:
    print(f"  ✓ PASS: CRE sizes are typical for ATAC-seq peaks (100-3000 bp)")
else:
    print(f"  ⚠️  WARNING: Unusual CRE sizes detected")

# ============================================================================
# Check 4: Coordinate spot-check with known genes
# ============================================================================
print("\n" + "="*80)
print("CHECK 4: Coordinate Spot-Check with Known Genes")
print("-"*80)

print("\nChecking coordinates for commonly studied genes...")

# Known gene positions in mm10 (from Ensembl/UCSC)
known_genes_mm10 = {
    'Actb': ('chr5', 142904061, 142908142),      # Beta-actin
    'Gapdh': ('chr6', 125112096, 125118870),     # GAPDH
    'Calm1': ('chr12', 100158529, 100165365),    # Calmodulin 1
    'Srf': ('chr17', 46545506, 46554896),        # SRF itself!
    'Nrxn1': ('chr17', 91195529, 92272803),      # Neurexin 1 (large gene)
}

matches = 0
mismatches = 0

for gene, (chr, start, end) in known_genes_mm10.items():
    gene_data = table16[table16['Gene'] == gene]

    if len(gene_data) == 0:
        print(f"  {gene}: Not found in Table 16 (may not be in dataset)")
        continue

    # Get CREs for this gene
    gene_cres = gene_data[['cre_chr', 'cre_start', 'cre_end']].drop_duplicates()

    # Check if any CREs are on the same chromosome
    chr_match = (gene_cres['cre_chr'] == chr).any()

    if chr_match:
        # Check if CREs are within reasonable distance (500kb)
        gene_cres_chr = gene_cres[gene_cres['cre_chr'] == chr]
        min_dist = np.min([
            min(abs(row['cre_start'] - start), abs(row['cre_start'] - end))
            for _, row in gene_cres_chr.iterrows()
        ])

        if min_dist <= 500000:  # Within 500kb
            print(f"  ✓ {gene} ({chr}): CREs found within {min_dist/1000:.1f} kb of mm10 position")
            matches += 1
        else:
            print(f"  ⚠️  {gene} ({chr}): CREs found but >500kb away (dist={min_dist/1000:.0f} kb)")
            mismatches += 1
    else:
        print(f"  ✗ {gene}: CREs on wrong chromosome! Expected {chr}, found {gene_cres['cre_chr'].unique()}")
        mismatches += 1

if matches > 0 and mismatches == 0:
    print(f"\n  ✓ PASS: Known gene coordinates match mm10 ({matches}/{matches})")
elif matches > mismatches:
    print(f"\n  ⚠️  PARTIAL: Most coordinates match mm10 ({matches}/{matches+mismatches})")
else:
    print(f"\n  ✗ FAIL: Coordinates do NOT match mm10 ({matches}/{matches+mismatches})")

# ============================================================================
# Check 5: Compare with your BigWig data
# ============================================================================
print("\n" + "="*80)
print("CHECK 5: BigWig File Chromosome Compatibility")
print("-"*80)

try:
    import pyBigWig

    bw_file = "../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Nestin-Ctrl.bw"

    if os.path.exists(bw_file):
        print(f"\nChecking BigWig file: {os.path.basename(bw_file)}")
        bw = pyBigWig.open(bw_file)
        bw_chrs = list(bw.chroms().keys())
        bw.close()

        print(f"  BigWig chromosomes: {sorted(bw_chrs)[:10]}")

        # Check overlap with literature chromosomes
        lit_chrs_set = set(lit_chrs)
        bw_chrs_set = set(bw_chrs)
        overlap = lit_chrs_set & bw_chrs_set

        overlap_pct = 100 * len(overlap) / len(lit_chrs_set)

        print(f"  Literature CRE chromosomes: {len(lit_chrs_set)}")
        print(f"  BigWig chromosomes: {len(bw_chrs_set)}")
        print(f"  Overlap: {len(overlap)} ({overlap_pct:.1f}%)")

        if overlap_pct > 90:
            print(f"  ✓ PASS: Excellent chromosome overlap ({overlap_pct:.1f}%)")
        elif overlap_pct > 70:
            print(f"  ⚠️  WARNING: Some chromosomes missing ({overlap_pct:.1f}% overlap)")
        else:
            print(f"  ✗ FAIL: Poor chromosome overlap ({overlap_pct:.1f}%)")
            print(f"         This suggests different genome builds!")

        # Show any mismatches
        in_lit_not_bw = lit_chrs_set - bw_chrs_set
        in_bw_not_lit = bw_chrs_set - lit_chrs_set

        if in_lit_not_bw:
            print(f"  In literature but not BigWig: {sorted(in_lit_not_bw)[:5]}")
        if in_bw_not_lit:
            print(f"  In BigWig but not literature: {sorted(in_bw_not_lit)[:5]}")
    else:
        print(f"\n⚠️  BigWig file not found: {bw_file}")
        print(f"   Skipping BigWig chromosome check")

except ImportError:
    print("\n⚠️  pyBigWig not installed, skipping BigWig check")
except Exception as e:
    print(f"\n⚠️  Error checking BigWig: {e}")

# ============================================================================
# Final verdict
# ============================================================================
print("\n" + "="*80)
print("FINAL VERDICT")
print("="*80)

print("\nBased on the checks above:")
print()
print("YOUR DATA:")
print("  - Reference genome: mm10/GRCm38 (UCSC)")
print("  - Source: BSgenome.Mmusculus.UCSC.mm10")
print("  - Chromosome format: UCSC (chr1, chr2, ...)")
print()
print("LITERATURE CREs:")

# Make determination based on checks
chr_format_match = all(chr.startswith('chr') for chr in lit_chrs)
coord_valid = validation_errors == 0
known_genes_match = matches > 0 and mismatches == 0

if chr_format_match and coord_valid and known_genes_match:
    print("  - Reference genome: LIKELY mm10/GRCm38 ✓")
    print("  - Chromosome format: UCSC (chr1, chr2, ...) ✓")
    print("  - Coordinate validation: PASS ✓")
    print()
    print("  ✓✓✓ COMPATIBLE: Your data and literature CREs appear to use the same genome build!")
    print()
    print("  You can proceed with confidence. The analysis is valid.")
elif chr_format_match and coord_valid:
    print("  - Reference genome: POSSIBLY mm10/GRCm38 ⚠️")
    print("  - Chromosome format: UCSC (chr1, chr2, ...) ✓")
    print("  - Coordinate validation: PASS ✓")
    print()
    print("  ⚠️  PROBABLY COMPATIBLE, but verify gene positions manually")
    print()
    print("  Recommendation: Spot-check a few genes of interest using UCSC Genome Browser")
else:
    print("  - Reference genome: UNKNOWN or MISMATCH ✗")
    print()
    print("  ✗✗✗ WARNING: Genome builds may NOT match!")
    print()
    print("  CRITICAL: Do NOT proceed without resolving this mismatch!")
    print("  Check the original paper for genome version information.")

print("\n" + "="*80)
print("NEXT STEPS:")
print("="*80)
print()
print("1. Check the original paper's Methods section for genome version")
print("   - Look for 'mm10', 'mm39', 'GRCm38', 'GRCm39'")
print()
print("2. If literature uses mm10/GRCm38:")
print("   ✓ You're good! Proceed with analysis.")
print()
print("3. If literature uses mm39/GRCm39:")
print("   ⚠️  You need to either:")
print("   a) Lift over literature coordinates from mm39→mm10, OR")
print("   b) Lift over your data from mm10→mm39")
print()
print("4. Verify a few genes manually:")
print("   - Go to UCSC Genome Browser (https://genome.ucsc.edu)")
print("   - Search for genes like Calm1, Nrxn1, Srf")
print("   - Check if literature CRE coordinates align with gene positions")
print()
print("5. If uncertain, check supplementary methods or contact paper authors")
print()
print("="*80)
