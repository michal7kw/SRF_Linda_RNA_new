#!/usr/bin/env python3
"""
Investigate genes of interest across all cell types

Check:
1. Are these genes expressed in GABA cells?
2. Are they DE in other cell types?
3. Where are the literature CREs located?
4. Where are your DA peaks located?

INPUT FILES:
- genes_inter.txt: List of genes of interest
- ../signac_results_L1/celltype_results/DEG/*_DEGs.csv: Differential expression results for all cell types
- output/unique_CREs_for_genes_of_interest_with_header.bed: Literature CREs for genes of interest
- ../signac_results_L1/celltype_results/filtered/GABA_*_DA_peaks.csv: Differential accessibility peaks for GABA cells

OUTPUT FILES:
- None (console output only)
"""

import pandas as pd
import numpy as np
import os

print("="*80)
print("INVESTIGATE GENES ACROSS CELL TYPES")
print("="*80)

# Configuration
SIGNAC_RESULTS_DIR = "../signac_results_L1"
GENES_FILE = "genes_inter.txt"
LITERATURE_CRE_DIR = "output"

# Load genes
genes_of_interest = []
with open(GENES_FILE, 'r') as f:
    for line in f:
        gene = line.strip()
        if gene and not gene.startswith('#'):
            genes_of_interest.append(gene)

print(f"\nGenes of interest ({len(genes_of_interest)}):")
print(", ".join(genes_of_interest))

# ============================================================================
# 1. Check DEGs across all cell types
# ============================================================================
print(f"\n{'='*80}")
print("1. Checking DEGs across ALL cell types...")
print("-"*80)

deg_dir = os.path.join(SIGNAC_RESULTS_DIR, "celltype_results", "DEG")
deg_files = [f for f in os.listdir(deg_dir) if f.endswith('_DEGs.csv')]

print(f"\nFound {len(deg_files)} DEG files")

all_degs = []
for deg_file in deg_files:
    # Parse filename: CellType_Genotype_DEGs.csv
    parts = deg_file.replace('_DEGs.csv', '').split('_')
    if len(parts) == 2:
        celltype, genotype = parts
    else:
        continue

    # Load DEGs
    deg_path = os.path.join(deg_dir, deg_file)
    deg = pd.read_csv(deg_path)

    # Filter for genes of interest
    deg_genes = deg[deg['gene'].isin(genes_of_interest)]

    # Add metadata
    deg_genes['celltype'] = celltype
    deg_genes['genotype'] = genotype

    all_degs.append(deg_genes)

all_degs_df = pd.concat(all_degs, ignore_index=True)

# Filter for significant
sig_degs = all_degs_df[all_degs_df['p_val_adj'] < 0.05]

print(f"\nTotal gene-celltype-genotype combinations: {len(all_degs_df)}")
print(f"Significant DEGs: {len(sig_degs)}")

if len(sig_degs) > 0:
    print("\nSignificant DEGs by gene:")
    for gene in genes_of_interest:
        gene_degs = sig_degs[sig_degs['gene'] == gene]
        if len(gene_degs) > 0:
            print(f"\n  {gene}:")
            for _, row in gene_degs.iterrows():
                print(f"    {row['celltype']} {row['genotype']}: log2FC={row['avg_log2FC']:.3f}, "
                      f"padj={row['p_val_adj']:.2e}, {row['direction']}")
else:
    print("\n⚠️  NO significant DEGs found for any gene of interest in any cell type!")

# Show expression stats (all, not just significant)
print("\n" + "="*80)
print("Expression in GABA cells (including non-significant):")
print("-"*80)

gaba_degs = all_degs_df[all_degs_df['celltype'] == 'GABA']
print(f"\nGenes found in GABA DEG results: {gaba_degs['gene'].nunique()}/{len(genes_of_interest)}")

if len(gaba_degs) > 0:
    for gene in genes_of_interest:
        gene_gaba = gaba_degs[gaba_degs['gene'] == gene]
        if len(gene_gaba) > 0:
            print(f"\n  {gene}:")
            for _, row in gene_gaba.iterrows():
                sig = "***" if row['p_val_adj'] < 0.05 else ""
                print(f"    {row['genotype']}: log2FC={row['avg_log2FC']:.3f}, "
                      f"padj={row['p_val_adj']:.2e}, pct.1={row['pct.1']:.3f}, pct.2={row['pct.2']:.3f} {sig}")
        else:
            print(f"  {gene}: NOT FOUND (not expressed or filtered out)")

# ============================================================================
# 2. Check literature CRE locations
# ============================================================================
print(f"\n{'='*80}")
print("2. Literature CRE locations...")
print("-"*80)

lit_cres_file = os.path.join(LITERATURE_CRE_DIR, "unique_CREs_for_genes_of_interest_with_header.bed")
lit_cres = pd.read_csv(lit_cres_file, sep='\t')

print(f"\nLiterature CREs: {len(lit_cres)}")
print("\nCRE locations:")
for _, cre in lit_cres.iterrows():
    print(f"  {cre['cCRE1']}: {cre['chr']}:{cre['start']}-{cre['end']} "
          f"({(cre['end']-cre['start'])/1000:.1f}kb, {cre['n_genes']} genes: {cre['Gene']})")

# Chromosome distribution
print("\nChromosome distribution:")
chr_dist = lit_cres['chr'].value_counts().sort_index()
for chrom, count in chr_dist.items():
    print(f"  {chrom}: {count} CREs")

# ============================================================================
# 3. Check DA peak locations (GABA)
# ============================================================================
print(f"\n{'='*80}")
print("3. DA peak locations (GABA cells)...")
print("-"*80)

for genotype in ['Nestin', 'Emx1']:
    da_file = os.path.join(SIGNAC_RESULTS_DIR, "celltype_results", "filtered",
                          f"GABA_{genotype}_DA_peaks.csv")

    if os.path.exists(da_file):
        da = pd.read_csv(da_file)

        # Parse coordinates
        da[['chr', 'start', 'end']] = da['peak'].str.split('-', expand=True)
        da['start'] = pd.to_numeric(da['start'])
        da['end'] = pd.to_numeric(da['end'])

        print(f"\n{genotype}: {len(da)} DA peaks")

        # Chromosome distribution
        chr_dist = da['chr'].value_counts().sort_index()
        print(f"Chromosomes represented: {len(chr_dist)}")
        print("Top chromosomes:")
        for chrom, count in chr_dist.head(10).items():
            print(f"  {chrom}: {count} peaks")

        # Check if any peaks are on same chromosomes as literature CREs
        lit_cre_chrs = set(lit_cres['chr'].unique())
        da_chrs = set(da['chr'].unique())
        common_chrs = lit_cre_chrs.intersection(da_chrs)

        print(f"\nChromosomes in common with literature CREs: {len(common_chrs)}")
        if common_chrs:
            print(f"  {', '.join(sorted(common_chrs))}")

            # For common chromosomes, show ranges
            for chrom in sorted(common_chrs):
                lit_cre_chr = lit_cres[lit_cres['chr'] == chrom]
                da_chr = da[da['chr'] == chrom]

                print(f"\n  {chrom}:")
                print(f"    Literature CREs: {len(lit_cre_chr)} CREs, "
                      f"range {lit_cre_chr['start'].min():,}-{lit_cre_chr['end'].max():,}")
                print(f"    DA peaks: {len(da_chr)} peaks, "
                      f"range {da_chr['start'].min():,}-{da_chr['end'].max():,}")

# ============================================================================
# 4. Recommendations
# ============================================================================
print(f"\n{'='*80}")
print("RECOMMENDATIONS")
print("="*80)

print("\nBased on this analysis:")

if len(sig_degs) == 0:
    print("\n1. ⚠️  None of your genes of interest are significantly DE in ANY cell type")
    print("   → These splicing factors may not be affected by SRF deletion")
    print("   → Consider analyzing other gene sets (e.g., SRF targets, cell type markers)")

if len(gaba_degs) < len(genes_of_interest):
    missing = set(genes_of_interest) - set(gaba_degs['gene'].unique())
    print(f"\n2. ⚠️  {len(missing)} genes not found in GABA DEG results:")
    print(f"   {', '.join(sorted(missing))}")
    print("   → These may have very low expression in GABA cells")
    print("   → Check expression in other cell types (Excitatory, Astrocytes, etc.)")

print("\n3. ✓ Literature CREs are from hippocampal interneurons")
print("   → Your GABA cells may include different interneuron subtypes")
print("   → Consider analyzing L2 (fine-grained) cell types for more specificity")

print("\n4. Next steps:")
print("   a. Check if genes are expressed but not DE (examine raw expression)")
print("   b. Expand to all cell types to find where these genes ARE dysregulated")
print("   c. Try using L2 cell type annotations for more specific interneuron subtypes")
print("   d. Analyze a different gene set (e.g., known SRF targets)")

print("\n" + "="*80)
