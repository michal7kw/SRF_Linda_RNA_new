#!/usr/bin/env python3
"""
Link genes of interest to their regulatory CREs using literature data

Working dir: "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature"

Purpose:
- Extract CRE-gene links from Table 16 for genes in genes_inter.txt
- Filter for hippocampal/GABAergic cell types
- Output CREs that regulate these genes
- Prepare for intersection with your differential ATAC peaks

INPUT FILES:
- data/genes_inter.txt: List of genes of interest
- data/table_16.txt: Literature CRE-gene correlation data

OUTPUT FILES:
- output/CREs_linked_to_genes_of_interest.bed: BED file with CREs linked to genes of interest
- output/CREs_linked_to_genes_of_interest_with_header.bed: Same BED file with column headers
- output/CRE_gene_links_detailed.tsv: Detailed table with all correlation statistics
- output/gene_summary.tsv: Per-gene summary statistics
- output/unique_CREs_for_genes_of_interest.bed: Deduplicated CREs for all genes
- output/unique_CREs_for_genes_of_interest_with_header.bed: Same BED file with headers
- output/CRE_gene_links_metadata.txt: Complete metadata report
"""

import pandas as pd
import numpy as np
import re
from datetime import datetime
import os

os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature")

print("="*80)
print("CRE-GENE LINKAGE EXTRACTION FOR GENES OF INTEREST")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# ============================================================================
# STEP 1: Load genes of interest
# ============================================================================
print("STEP 1: Loading genes of interest...")
print("-"*80)

genes_of_interest = []
with open('data/genes_inter.txt', 'r') as f:
    for line in f:
        gene = line.strip()
        if gene and not gene.startswith('#'):
            genes_of_interest.append(gene)

print(f"Loaded {len(genes_of_interest)} genes of interest:")
for gene in genes_of_interest:
    print(f"  {gene}")

# ============================================================================
# STEP 2: Load CRE-gene correlations from Table 16
# ============================================================================
print(f"\n{'='*80}")
print("STEP 2: Loading CRE-gene correlations (Table 16)...")
print("-"*80)
print("This may take 1-2 minutes (567 MB file)...")

table16 = pd.read_csv('data/table_16.txt', sep='\t')
print(f"✓ Loaded {len(table16):,} CRE-gene correlation entries")

print("\nColumns available:")
for col in table16.columns:
    print(f"  {col}")

print("\nSample data:")
print(table16.head(3))

# ============================================================================
# STEP 3: Filter for genes of interest
# ============================================================================
print(f"\n{'='*80}")
print("STEP 3: Filtering for genes of interest...")
print("-"*80)

# Filter for genes of interest
genes_links = table16[table16['Gene'].isin(genes_of_interest)].copy()
print(f"Found {len(genes_links):,} CRE-gene links for genes of interest")

# Show distribution
print("\nLinks per gene:")
gene_counts = genes_links['Gene'].value_counts()
for gene in genes_of_interest:
    if gene in gene_counts:
        print(f"  {gene}: {gene_counts[gene]:,} CRE links")
    else:
        print(f"  {gene}: 0 CRE links (NOT FOUND)")

# ============================================================================
# STEP 4: Filter by cell type (hippocampal/GABAergic)
# ============================================================================
print(f"\n{'='*80}")
print("STEP 4: Filtering by cell type...")
print("-"*80)

# Check unique SubTypes
print(f"\nUnique SubTypes in Table 16 (sample):")
unique_subtypes = table16['SubType'].unique()
print(f"Total unique subtypes: {len(unique_subtypes)}")
print("Sample:")
print(unique_subtypes[:30])

# Define hippocampal/GABAergic subtypes
# Based on our previous analysis and Table 17
hippocampal_keywords = [
    'CA1', 'CA2', 'CA3', 'DG', 'DGNBL', 'GRC',  # Hippocampal
    'LAMP5', 'LAMP', 'VIP', 'SST', 'PV', 'PVGA', 'SSTGA', 'VIPGA', 'LAMGA',  # GABAergic interneurons
    'GABA', 'INH', 'interneuron'
]

# Filter for hippocampal/GABAergic cell types
def is_hippocampal_or_gaba(subtype):
    if pd.isna(subtype):
        return False
    subtype_str = str(subtype).upper()
    return any(keyword.upper() in subtype_str for keyword in hippocampal_keywords)

genes_links['is_hippo_gaba'] = genes_links['SubType'].apply(is_hippocampal_or_gaba)

print(f"\nLinks in hippocampal/GABAergic cell types:")
print(f"  Total links: {genes_links['is_hippo_gaba'].sum():,} / {len(genes_links):,}")
print(f"  ({100 * genes_links['is_hippo_gaba'].sum() / len(genes_links):.1f}%)")

# Check subtypes that matched
matched_subtypes = genes_links[genes_links['is_hippo_gaba']]['SubType'].unique()
print(f"\nMatched SubTypes ({len(matched_subtypes)}):")
for subtype in sorted(matched_subtypes)[:20]:
    count = (genes_links['SubType'] == subtype).sum()
    print(f"  {subtype}: {count:,} links")
if len(matched_subtypes) > 20:
    print(f"  ... and {len(matched_subtypes) - 20} more")

# Filter for hippocampal/GABAergic
genes_links_filtered = genes_links[genes_links['is_hippo_gaba']].copy()

# ============================================================================
# STEP 5: Apply statistical filters
# ============================================================================
print(f"\n{'='*80}")
print("STEP 5: Applying statistical filters...")
print("-"*80)

print(f"Before filtering: {len(genes_links_filtered):,} links")

# Check distribution of statistics
print("\nStatistical measures:")
print(f"  PCC (Pearson correlation): min={genes_links_filtered['PCC'].min():.3f}, "
      f"median={genes_links_filtered['PCC'].median():.3f}, max={genes_links_filtered['PCC'].max():.3f}")
print(f"  FDR: min={genes_links_filtered['FDR'].min():.6f}, "
      f"median={genes_links_filtered['FDR'].median():.6f}, max={genes_links_filtered['FDR'].max():.6f}")

# Apply filters
# Option 1: Lenient (FDR < 0.1, |PCC| > 0.1)
# Option 2: Moderate (FDR < 0.05, |PCC| > 0.2)
# Option 3: Stringent (FDR < 0.01, |PCC| > 0.3)

# Let's try moderate first
FDR_THRESHOLD = 0.05
PCC_THRESHOLD = 0.2

genes_links_sig = genes_links_filtered[
    (genes_links_filtered['FDR'] < FDR_THRESHOLD) &
    (genes_links_filtered['PCC'].abs() > PCC_THRESHOLD)
].copy()

print(f"\nAfter filtering (FDR < {FDR_THRESHOLD}, |PCC| > {PCC_THRESHOLD}):")
print(f"  Significant links: {len(genes_links_sig):,}")
print(f"  Genes with links: {genes_links_sig['Gene'].nunique()}")
print(f"  Unique CREs: {genes_links_sig['cCRE1'].nunique()}")

# Show per-gene breakdown
print("\nSignificant links per gene:")
gene_sig_counts = genes_links_sig['Gene'].value_counts()
for gene in genes_of_interest:
    if gene in gene_sig_counts:
        print(f"  {gene}: {gene_sig_counts[gene]:,} significant CRE links")
    else:
        print(f"  {gene}: 0 significant links")

# ============================================================================
# STEP 6: Parse CRE coordinates
# ============================================================================
print(f"\n{'='*80}")
print("STEP 6: Parsing CRE coordinates...")
print("-"*80)

def parse_coordinate(coord_str):
    """Parse coordinate format: chr11_40754370_40756166"""
    try:
        parts = coord_str.split('_')
        if len(parts) == 3:
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            return pd.Series({'chr': chrom, 'start': start, 'end': end})
        else:
            return pd.Series({'chr': None, 'start': None, 'end': None})
    except:
        return pd.Series({'chr': None, 'start': None, 'end': None})

print("Parsing Coordinate1 (CRE locations)...")
parsed_coords = genes_links_sig['Coordinate1'].apply(parse_coordinate)
genes_links_sig = pd.concat([genes_links_sig, parsed_coords], axis=1)

# Remove unparsed entries
genes_links_sig = genes_links_sig.dropna(subset=['chr', 'start', 'end'])
print(f"Successfully parsed {len(genes_links_sig):,} CRE coordinates")

# ============================================================================
# STEP 7: Create output files
# ============================================================================
print(f"\n{'='*80}")
print("STEP 7: Creating output files...")
print("-"*80)

import os
os.makedirs('output', exist_ok=True)

# Sort by chromosome and position
def chr_sort_key(chr_name):
    chr_name = str(chr_name)
    if chr_name.startswith('chr'):
        chr_name = chr_name[3:]
    if chr_name.isdigit():
        return (0, int(chr_name))
    elif chr_name == 'X':
        return (1, 0)
    elif chr_name == 'Y':
        return (1, 1)
    elif chr_name in ['M', 'MT']:
        return (1, 2)
    else:
        return (2, chr_name)

genes_links_sig['chr_sort'] = genes_links_sig['chr'].apply(chr_sort_key)
genes_links_sig = genes_links_sig.sort_values(['chr_sort', 'start', 'end']).drop('chr_sort', axis=1)

# Output 1: BED file of CREs linked to genes of interest
bed_output = genes_links_sig[['chr', 'start', 'end', 'cCRE1', 'Gene', 'PCC', 'FDR', 'SubType']].copy()
bed_file = 'output/CREs_linked_to_genes_of_interest.bed'
bed_output.to_csv(bed_file, sep='\t', index=False, header=False)
print(f"✓ BED file saved: {bed_file}")
print(f"  Format: chr, start, end, cre_id, gene, PCC, FDR, subtype")

# Output 2: BED file with header
bed_file_header = 'output/CREs_linked_to_genes_of_interest_with_header.bed'
bed_output.to_csv(bed_file_header, sep='\t', index=False)
print(f"✓ BED file with header saved: {bed_file_header}")

# Output 3: Detailed table with all information
detailed_output = genes_links_sig[[
    'Gene', 'chr', 'start', 'end', 'cCRE1', 'PCC', 'Pval', 'FDR',
    'SubType', 'Type', 'Coaccess', 'Conns'
]].copy()
detailed_file = 'output/CRE_gene_links_detailed.tsv'
detailed_output.to_csv(detailed_file, sep='\t', index=False)
print(f"✓ Detailed table saved: {detailed_file}")

# Output 4: Gene-centric summary
print("\nCreating gene-centric summary...")
gene_summary = []
for gene in genes_of_interest:
    gene_data = genes_links_sig[genes_links_sig['Gene'] == gene]
    n_cres = len(gene_data)
    n_unique_cres = gene_data['cCRE1'].nunique()
    if n_cres > 0:
        avg_pcc = gene_data['PCC'].mean()
        median_pcc = gene_data['PCC'].median()
        top_cre = gene_data.nlargest(1, 'PCC').iloc[0]
        top_cre_info = f"{top_cre['cCRE1']} ({top_cre['chr']}:{top_cre['start']}-{top_cre['end']}, PCC={top_cre['PCC']:.3f})"
    else:
        avg_pcc = np.nan
        median_pcc = np.nan
        top_cre_info = "N/A"

    gene_summary.append({
        'Gene': gene,
        'N_CRE_links': n_cres,
        'N_unique_CREs': n_unique_cres,
        'Avg_PCC': avg_pcc,
        'Median_PCC': median_pcc,
        'Top_CRE': top_cre_info
    })

gene_summary_df = pd.DataFrame(gene_summary)
gene_summary_file = 'output/gene_summary.tsv'
gene_summary_df.to_csv(gene_summary_file, sep='\t', index=False)
print(f"✓ Gene summary saved: {gene_summary_file}")

# Output 5: Unique CREs (deduplicated)
unique_cres = genes_links_sig.groupby('cCRE1').agg({
    'chr': 'first',
    'start': 'first',
    'end': 'first',
    'Gene': lambda x: ','.join(sorted(set(x))),
    'PCC': 'mean',
    'FDR': 'min',
    'SubType': lambda x: ','.join(sorted(set(x)))
}).reset_index()

unique_cres['n_genes'] = unique_cres['Gene'].apply(lambda x: len(x.split(',')))
unique_cres = unique_cres.sort_values(['chr_sort', 'start', 'end']) if 'chr_sort' in unique_cres.columns else unique_cres

# Sort properly
unique_cres['chr_sort'] = unique_cres['chr'].apply(chr_sort_key)
unique_cres = unique_cres.sort_values(['chr_sort', 'start', 'end']).drop('chr_sort', axis=1)

unique_cres_bed = unique_cres[['chr', 'start', 'end', 'cCRE1', 'n_genes', 'Gene', 'PCC', 'FDR']].copy()
unique_cres_file = 'output/unique_CREs_for_genes_of_interest.bed'
unique_cres_bed.to_csv(unique_cres_file, sep='\t', index=False, header=False)
print(f"✓ Unique CREs saved: {unique_cres_file}")

unique_cres_file_header = 'output/unique_CREs_for_genes_of_interest_with_header.bed'
unique_cres_bed.to_csv(unique_cres_file_header, sep='\t', index=False)
print(f"✓ Unique CREs with header saved: {unique_cres_file_header}")

# Output 6: Metadata
metadata_file = 'output/CRE_gene_links_metadata.txt'
with open(metadata_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("CRE-GENE LINKAGE EXTRACTION - METADATA\n")
    f.write("="*80 + "\n\n")
    f.write(f"Extraction date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("SOURCE DATA:\n")
    f.write(f"  Table 16: CRE-gene correlations from literature\n")
    f.write(f"  Genes of interest: genes_inter.txt ({len(genes_of_interest)} genes)\n\n")

    f.write("FILTERS APPLIED:\n")
    f.write(f"  Cell types: Hippocampal and GABAergic\n")
    f.write(f"  FDR threshold: < {FDR_THRESHOLD}\n")
    f.write(f"  |PCC| threshold: > {PCC_THRESHOLD}\n\n")

    f.write("RESULTS:\n")
    f.write(f"  Genes with significant links: {genes_links_sig['Gene'].nunique()} / {len(genes_of_interest)}\n")
    f.write(f"  Total CRE-gene links: {len(genes_links_sig):,}\n")
    f.write(f"  Unique CREs: {genes_links_sig['cCRE1'].nunique():,}\n")
    f.write(f"  Average links per gene: {len(genes_links_sig) / genes_links_sig['Gene'].nunique():.1f}\n\n")

    f.write("GENES WITH SIGNIFICANT CRE LINKS:\n")
    for gene in genes_of_interest:
        count = gene_sig_counts.get(gene, 0)
        f.write(f"  {gene}: {count:,} CRE links\n")

    f.write("\n" + "="*80 + "\n")
    f.write("OUTPUT FILES:\n")
    f.write("="*80 + "\n")
    f.write(f"\n{bed_file}\n")
    f.write(f"  BED file with all CRE-gene links\n")
    f.write(f"  Format: chr, start, end, cre_id, gene, PCC, FDR, subtype\n\n")

    f.write(f"{unique_cres_file}\n")
    f.write(f"  Unique CREs (deduplicated across genes)\n")
    f.write(f"  Format: chr, start, end, cre_id, n_genes, gene_list, avg_PCC, min_FDR\n\n")

    f.write(f"{detailed_file}\n")
    f.write(f"  Detailed table with all correlation statistics\n\n")

    f.write(f"{gene_summary_file}\n")
    f.write(f"  Per-gene summary statistics\n\n")

    f.write("\nNEXT STEPS:\n")
    f.write("="*80 + "\n")
    f.write("1. Intersect CREs with your differential ATAC peaks:\n")
    f.write("   bedtools intersect -a signac_results/celltype_results/DA/YourCellType_peaks.bed \\\n")
    f.write("                      -b output/unique_CREs_for_genes_of_interest.bed -wa -wb\n\n")
    f.write("2. Check if CRE accessibility changes correlate with gene expression changes\n\n")
    f.write("3. Use the detailed table to identify cell-type-specific regulatory relationships\n\n")

print(f"✓ Metadata saved: {metadata_file}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print(f"\n{'='*80}")
print("EXTRACTION COMPLETE!")
print("="*80)

print(f"\nGenes analyzed: {len(genes_of_interest)}")
print(f"Genes with CRE links: {genes_links_sig['Gene'].nunique()} ({100*genes_links_sig['Gene'].nunique()/len(genes_of_interest):.0f}%)")
print(f"\nTotal CRE-gene links: {len(genes_links_sig):,}")
print(f"Unique CREs: {genes_links_sig['cCRE1'].nunique():,}")
print(f"Average CRE correlation: {genes_links_sig['PCC'].mean():.3f} (median: {genes_links_sig['PCC'].median():.3f})")

print("\nTop genes by number of CRE links:")
for gene, count in gene_sig_counts.head(10).items():
    print(f"  {gene}: {count:,} CRE links")

print(f"\nOutput files saved in: output/")
print(f"  - CREs_linked_to_genes_of_interest.bed")
print(f"  - unique_CREs_for_genes_of_interest.bed")
print(f"  - CRE_gene_links_detailed.tsv")
print(f"  - gene_summary.tsv")
print(f"  - CRE_gene_links_metadata.txt")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)

# Print quick preview of gene summary
print("\nGENE SUMMARY:")
print("-"*80)
print(gene_summary_df.to_string(index=False))
