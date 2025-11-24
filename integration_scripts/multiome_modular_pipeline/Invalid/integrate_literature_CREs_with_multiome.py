#!/usr/bin/env python3
"""
Integrate Literature CREs with Multiome Results

Purpose:
- For genes in genes_inter.txt, check if their expression changes can be
  explained by accessibility changes in literature-validated CREs
- Compare Nestin Ctrl vs Mut and Emx1 Ctrl vs Mut separately
- Use GABA cell type (L1 broad annotation)

Analysis:
1. Load literature CRE-gene links
2. Load DEGs (GABA Nestin/Emx1)
3. Load DA peaks (GABA Nestin/Emx1)
4. Intersect literature CREs with DA peaks
5. For each gene: check if expression change is consistent with CRE accessibility

INPUT FILES:
- genes_inter.txt: List of genes of interest
- output/unique_CREs_for_genes_of_interest_with_header.bed: Literature CREs for genes of interest
- ../signac_results_L1/celltype_results/DEG/GABA_Nestin_DEGs.csv: Differential expression for Nestin genotype
- ../signac_results_L1/celltype_results/DEG/GABA_Emx1_DEGs.csv: Differential expression for Emx1 genotype
- ../signac_results_L1/celltype_results/filtered/GABA_Nestin_DA_peaks.csv: Differential accessibility for Nestin genotype
- ../signac_results_L1/celltype_results/filtered/GABA_Emx1_DA_peaks.csv: Differential accessibility for Emx1 genotype

OUTPUT FILES:
- output/multiome_integration/literature_CRE_support_GABA_*.tsv: Integration results for each genotype
- output/multiome_integration/CRE_peak_overlaps_GABA_*.tsv: CRE-DA peak overlap data
- output/multiome_integration/summary_literature_CRE_support.tsv: Summary table of support status
- output/multiome_integration/analysis_report.txt: Detailed analysis report
"""

import pandas as pd
import numpy as np
from datetime import datetime
import os

print("="*80)
print("INTEGRATE LITERATURE CREs WITH MULTIOME RESULTS")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# ============================================================================
# Configuration
# ============================================================================
SIGNAC_RESULTS_DIR = "../signac_results_L1"
LITERATURE_CRE_DIR = "output"
GENES_FILE = "genes_inter.txt"
OUTPUT_DIR = "output/multiome_integration"

CELL_TYPE = "GABA"
GENOTYPES = ["Nestin", "Emx1"]

# Statistical thresholds
PADJ_THRESHOLD_DEG = 0.05
PADJ_THRESHOLD_DA = 0.05
LOG2FC_THRESHOLD = 0.25  # Minimum |log2FC| for biological relevance

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# STEP 1: Load genes of interest
# ============================================================================
print("STEP 1: Loading genes of interest...")
print("-"*80)

genes_of_interest = []
with open(GENES_FILE, 'r') as f:
    for line in f:
        gene = line.strip()
        if gene and not gene.startswith('#'):
            genes_of_interest.append(gene)

print(f"Loaded {len(genes_of_interest)} genes of interest")

# ============================================================================
# STEP 2: Load literature CRE-gene links
# ============================================================================
print(f"\n{'='*80}")
print("STEP 2: Loading literature CRE-gene links...")
print("-"*80)

# Load the unique CREs file (deduplicated)
lit_cres_file = os.path.join(LITERATURE_CRE_DIR, "unique_CREs_for_genes_of_interest_with_header.bed")
if not os.path.exists(lit_cres_file):
    print(f"ERROR: Literature CRE file not found: {lit_cres_file}")
    print("Please run link_genes_to_cres.py first!")
    exit(1)

lit_cres = pd.read_csv(lit_cres_file, sep='\t')
print(f"Loaded {len(lit_cres)} unique literature CREs")
print(f"Covering {lit_cres['Gene'].apply(lambda x: len(x.split(','))).sum()} gene-CRE links")

# Parse gene list (can be comma-separated)
lit_cres['genes_list'] = lit_cres['Gene'].str.split(',')

# Create expanded version with one row per gene-CRE pair
lit_cres_expanded = lit_cres.explode('genes_list')
lit_cres_expanded = lit_cres_expanded.rename(columns={'genes_list': 'gene'})
print(f"Expanded to {len(lit_cres_expanded)} gene-CRE pairs")

# ============================================================================
# STEP 3: Load differential expression results
# ============================================================================
print(f"\n{'='*80}")
print("STEP 3: Loading differential expression (DEG) results...")
print("-"*80)

deg_results = {}
for genotype in GENOTYPES:
    deg_file = os.path.join(SIGNAC_RESULTS_DIR, "celltype_results", "DEG",
                            f"{CELL_TYPE}_{genotype}_DEGs.csv")

    if os.path.exists(deg_file):
        deg = pd.read_csv(deg_file)

        # Filter for significance
        deg_sig = deg[deg['p_val_adj'] < PADJ_THRESHOLD_DEG].copy()

        # Filter for genes of interest
        deg_sig_genes = deg_sig[deg_sig['gene'].isin(genes_of_interest)].copy()

        deg_results[genotype] = {
            'all': deg,
            'significant': deg_sig,
            'genes_of_interest': deg_sig_genes
        }

        print(f"\n{genotype}:")
        print(f"  Total DEGs: {len(deg):,}")
        print(f"  Significant DEGs (padj < {PADJ_THRESHOLD_DEG}): {len(deg_sig):,}")
        print(f"  Genes of interest (significant): {len(deg_sig_genes)}")

        if len(deg_sig_genes) > 0:
            print(f"  Significant genes of interest:")
            for _, row in deg_sig_genes.iterrows():
                print(f"    {row['gene']}: log2FC={row['avg_log2FC']:.3f}, padj={row['p_val_adj']:.2e}, {row['direction']}")
    else:
        print(f"\n{genotype}: DEG file not found: {deg_file}")
        deg_results[genotype] = None

# ============================================================================
# STEP 4: Load differential accessibility results
# ============================================================================
print(f"\n{'='*80}")
print("STEP 4: Loading differential accessibility (DA) results...")
print("-"*80)

da_results = {}
for genotype in GENOTYPES:
    da_file = os.path.join(SIGNAC_RESULTS_DIR, "celltype_results", "filtered",
                          f"{CELL_TYPE}_{genotype}_DA_peaks.csv")

    if os.path.exists(da_file):
        da = pd.read_csv(da_file)

        # Parse peak coordinates: chr-start-end
        da[['chr', 'start', 'end']] = da['peak'].str.split('-', expand=True)
        da['start'] = pd.to_numeric(da['start'])
        da['end'] = pd.to_numeric(da['end'])

        # Filter for significance
        da_sig = da[da['p_val_adj'] < PADJ_THRESHOLD_DA].copy()

        da_results[genotype] = {
            'all': da,
            'significant': da_sig
        }

        print(f"\n{genotype}:")
        print(f"  Total DA peaks: {len(da):,}")
        print(f"  Significant DA peaks (padj < {PADJ_THRESHOLD_DA}): {len(da_sig):,}")
        print(f"    More accessible: {(da_sig['direction'] == 'More_Accessible').sum():,}")
        print(f"    Less accessible: {(da_sig['direction'] == 'Less_Accessible').sum():,}")
    else:
        print(f"\n{genotype}: DA peaks file not found: {da_file}")
        da_results[genotype] = None

# ============================================================================
# STEP 5: Intersect literature CREs with DA peaks
# ============================================================================
print(f"\n{'='*80}")
print("STEP 5: Intersecting literature CREs with DA peaks...")
print("-"*80)

def find_overlaps(cres_df, peaks_df, min_overlap=1):
    """
    Find overlaps between CREs and peaks
    Returns DataFrame with overlapping CRE-peak pairs
    """
    overlaps = []

    for _, cre in cres_df.iterrows():
        cre_chr = cre['chr']
        cre_start = cre['start']
        cre_end = cre['end']

        # Find peaks on same chromosome that overlap
        peaks_chr = peaks_df[peaks_df['chr'] == cre_chr]

        for _, peak in peaks_chr.iterrows():
            peak_start = peak['start']
            peak_end = peak['end']

            # Check for overlap
            overlap_start = max(cre_start, peak_start)
            overlap_end = min(cre_end, peak_end)
            overlap_size = overlap_end - overlap_start

            if overlap_size >= min_overlap:
                overlaps.append({
                    'cre_id': cre['cCRE1'],
                    'cre_chr': cre_chr,
                    'cre_start': cre_start,
                    'cre_end': cre_end,
                    'cre_genes': cre['Gene'],
                    'peak': peak['peak'],
                    'peak_chr': peak['chr'],
                    'peak_start': peak_start,
                    'peak_end': peak_end,
                    'peak_log2fc': peak['avg_log2FC'],
                    'peak_padj': peak['p_val_adj'],
                    'peak_direction': peak['direction'],
                    'overlap_size': overlap_size
                })

    return pd.DataFrame(overlaps)

# Find overlaps for each genotype
overlap_results = {}
for genotype in GENOTYPES:
    if da_results[genotype] is not None:
        print(f"\n{genotype}:")
        print(f"  Intersecting {len(lit_cres)} literature CREs with {len(da_results[genotype]['significant']):,} DA peaks...")

        overlaps = find_overlaps(lit_cres, da_results[genotype]['significant'])
        overlap_results[genotype] = overlaps

        print(f"  Found {len(overlaps)} CRE-peak overlaps")

        if len(overlaps) > 0:
            print(f"  Unique CREs with DA peaks: {overlaps['cre_id'].nunique()}")
            print(f"  Genes affected: {overlaps['cre_genes'].str.split(',').explode().nunique()}")
    else:
        overlap_results[genotype] = pd.DataFrame()

# ============================================================================
# STEP 6: Integrate DEG + DA + Literature CREs
# ============================================================================
print(f"\n{'='*80}")
print("STEP 6: Integrating DEG + DA + Literature CREs...")
print("-"*80)

integration_results = {}

for genotype in GENOTYPES:
    print(f"\n{genotype}:")

    if deg_results[genotype] is None or da_results[genotype] is None:
        print("  Skipping (missing data)")
        continue

    # Get genes of interest with DEG results
    deg_genes = deg_results[genotype]['genes_of_interest']

    if len(deg_genes) == 0:
        print("  No significant DEGs for genes of interest")
        continue

    # For each DEG, find literature CREs and check DA status
    gene_reports = []

    for _, deg_row in deg_genes.iterrows():
        gene = deg_row['gene']
        gene_log2fc = deg_row['avg_log2FC']
        gene_padj = deg_row['p_val_adj']
        gene_direction = deg_row['direction']

        # Find literature CREs for this gene
        gene_cres = lit_cres_expanded[lit_cres_expanded['gene'] == gene]

        if len(gene_cres) == 0:
            gene_reports.append({
                'gene': gene,
                'gene_log2fc': gene_log2fc,
                'gene_padj': gene_padj,
                'gene_direction': gene_direction,
                'n_literature_cres': 0,
                'n_da_cres': 0,
                'da_cres_consistent': 0,
                'da_cres_inconsistent': 0,
                'support_status': 'NO_LITERATURE_CRES',
                'details': 'No literature CREs found for this gene'
            })
            continue

        # Check which CREs have DA peaks
        gene_cre_ids = gene_cres['cCRE1'].tolist()
        overlaps_for_gene = overlap_results[genotype][
            overlap_results[genotype]['cre_id'].isin(gene_cre_ids)
        ]

        n_cres_total = len(gene_cre_ids)
        n_cres_with_da = len(overlaps_for_gene)

        if n_cres_with_da == 0:
            gene_reports.append({
                'gene': gene,
                'gene_log2fc': gene_log2fc,
                'gene_padj': gene_padj,
                'gene_direction': gene_direction,
                'n_literature_cres': n_cres_total,
                'n_da_cres': 0,
                'da_cres_consistent': 0,
                'da_cres_inconsistent': 0,
                'support_status': 'NO_DA_PEAKS',
                'details': f'{n_cres_total} literature CREs, but none are DA'
            })
            continue

        # Check consistency: gene down + CRE less accessible = consistent
        n_consistent = 0
        n_inconsistent = 0
        cre_details = []

        for _, overlap in overlaps_for_gene.iterrows():
            cre_direction = overlap['peak_direction']

            # Check consistency
            if gene_direction == 'Down' and cre_direction == 'Less_Accessible':
                consistent = True
                n_consistent += 1
            elif gene_direction == 'Up' and cre_direction == 'More_Accessible':
                consistent = True
                n_consistent += 1
            else:
                consistent = False
                n_inconsistent += 1

            cre_details.append(
                f"{overlap['cre_id']}({cre_direction},log2FC={overlap['peak_log2fc']:.2f})"
            )

        # Determine support status
        if n_consistent > 0 and n_inconsistent == 0:
            support_status = 'STRONG_SUPPORT'
        elif n_consistent > n_inconsistent:
            support_status = 'PARTIAL_SUPPORT'
        elif n_consistent > 0:
            support_status = 'WEAK_SUPPORT'
        else:
            support_status = 'NO_SUPPORT'

        gene_reports.append({
            'gene': gene,
            'gene_log2fc': gene_log2fc,
            'gene_padj': gene_padj,
            'gene_direction': gene_direction,
            'n_literature_cres': n_cres_total,
            'n_da_cres': n_cres_with_da,
            'da_cres_consistent': n_consistent,
            'da_cres_inconsistent': n_inconsistent,
            'support_status': support_status,
            'details': f"{n_consistent}/{n_cres_with_da} DA CREs consistent; " + '; '.join(cre_details)
        })

    integration_df = pd.DataFrame(gene_reports)
    integration_results[genotype] = integration_df

    # Print summary
    print(f"  Genes analyzed: {len(integration_df)}")
    for status in ['STRONG_SUPPORT', 'PARTIAL_SUPPORT', 'WEAK_SUPPORT', 'NO_SUPPORT', 'NO_DA_PEAKS', 'NO_LITERATURE_CRES']:
        count = (integration_df['support_status'] == status).sum()
        if count > 0:
            print(f"    {status}: {count}")

# ============================================================================
# STEP 7: Create output files
# ============================================================================
print(f"\n{'='*80}")
print("STEP 7: Creating output files...")
print("-"*80)

# Save integration results
for genotype in GENOTYPES:
    if genotype in integration_results:
        output_file = os.path.join(OUTPUT_DIR, f"literature_CRE_support_{CELL_TYPE}_{genotype}.tsv")
        integration_results[genotype].to_csv(output_file, sep='\t', index=False)
        print(f"✓ Saved: {output_file}")

# Save CRE-peak overlaps
for genotype in GENOTYPES:
    if genotype in overlap_results and len(overlap_results[genotype]) > 0:
        output_file = os.path.join(OUTPUT_DIR, f"CRE_peak_overlaps_{CELL_TYPE}_{genotype}.tsv")
        overlap_results[genotype].to_csv(output_file, sep='\t', index=False)
        print(f"✓ Saved: {output_file}")

# Create combined summary
print("\nCreating combined summary...")
summary_rows = []

for genotype in GENOTYPES:
    if genotype in integration_results:
        df = integration_results[genotype]

        for status in ['STRONG_SUPPORT', 'PARTIAL_SUPPORT', 'WEAK_SUPPORT', 'NO_SUPPORT', 'NO_DA_PEAKS', 'NO_LITERATURE_CRES']:
            count = (df['support_status'] == status).sum()
            summary_rows.append({
                'genotype': genotype,
                'support_status': status,
                'n_genes': count
            })

summary_df = pd.DataFrame(summary_rows)
summary_file = os.path.join(OUTPUT_DIR, "summary_literature_CRE_support.tsv")
summary_df.to_csv(summary_file, sep='\t', index=False)
print(f"✓ Saved: {summary_file}")

# Create detailed report
report_file = os.path.join(OUTPUT_DIR, "analysis_report.txt")
with open(report_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("LITERATURE CRE - MULTIOME INTEGRATION REPORT\n")
    f.write("="*80 + "\n\n")
    f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("CONFIGURATION:\n")
    f.write(f"  Cell type: {CELL_TYPE}\n")
    f.write(f"  Genotypes: {', '.join(GENOTYPES)}\n")
    f.write(f"  Genes analyzed: {len(genes_of_interest)}\n")
    f.write(f"  Literature CREs: {len(lit_cres)}\n\n")

    f.write("RESULTS BY GENOTYPE:\n")
    f.write("-"*80 + "\n")

    for genotype in GENOTYPES:
        if genotype not in integration_results:
            continue

        f.write(f"\n{genotype}:\n")
        df = integration_results[genotype]

        # Overall stats
        f.write(f"  Total genes analyzed: {len(df)}\n")
        f.write(f"  Genes with literature CREs: {(df['n_literature_cres'] > 0).sum()}\n")
        f.write(f"  Genes with DA CREs: {(df['n_da_cres'] > 0).sum()}\n\n")

        # By support status
        f.write(f"  Support status:\n")
        for status in ['STRONG_SUPPORT', 'PARTIAL_SUPPORT', 'WEAK_SUPPORT', 'NO_SUPPORT', 'NO_DA_PEAKS', 'NO_LITERATURE_CRES']:
            genes = df[df['support_status'] == status]['gene'].tolist()
            if len(genes) > 0:
                f.write(f"    {status}: {len(genes)} genes\n")
                for gene in genes:
                    f.write(f"      - {gene}\n")

        f.write("\n  Detailed gene reports:\n")
        for _, row in df.iterrows():
            f.write(f"\n    {row['gene']}:\n")
            f.write(f"      Expression: {row['gene_direction']} (log2FC={row['gene_log2fc']:.3f}, padj={row['gene_padj']:.2e})\n")
            f.write(f"      Literature CREs: {row['n_literature_cres']}\n")
            f.write(f"      DA CREs: {row['n_da_cres']} ({row['da_cres_consistent']} consistent, {row['da_cres_inconsistent']} inconsistent)\n")
            f.write(f"      Support: {row['support_status']}\n")
            f.write(f"      Details: {row['details']}\n")

print(f"✓ Saved: {report_file}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print(f"\n{'='*80}")
print("ANALYSIS COMPLETE!")
print("="*80)

print("\nKey findings:")
for genotype in GENOTYPES:
    if genotype in integration_results:
        df = integration_results[genotype]
        n_strong = (df['support_status'] == 'STRONG_SUPPORT').sum()
        n_partial = (df['support_status'] == 'PARTIAL_SUPPORT').sum()
        n_total = len(df)

        print(f"\n{genotype}:")
        print(f"  {n_strong}/{n_total} genes with STRONG CRE support")
        print(f"  {n_partial}/{n_total} genes with PARTIAL CRE support")

        # Show genes with strong support
        strong_genes = df[df['support_status'] == 'STRONG_SUPPORT']['gene'].tolist()
        if strong_genes:
            print(f"  Strong support genes: {', '.join(strong_genes)}")

print(f"\nOutput files in: {OUTPUT_DIR}/")
print(f"  - literature_CRE_support_{{genotype}}.tsv")
print(f"  - CRE_peak_overlaps_{{genotype}}.tsv")
print(f"  - summary_literature_CRE_support.tsv")
print(f"  - analysis_report.txt")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
