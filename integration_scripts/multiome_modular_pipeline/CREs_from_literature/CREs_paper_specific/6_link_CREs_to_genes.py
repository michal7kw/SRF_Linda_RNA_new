#!/usr/bin/env python3
"""
Link Cell-Type CREs to Genes Using Table 16

This script links the hippocampal interneuron CREs and excitatory neuron CREs
to genes using published CRE-gene correlations from ENCODE Table 16.

Note: This pipeline uses OVERLAPPING CRE sets (60% overlap) - for mutually
exclusive CRE analysis, see CREs_paper_exclusive/

INPUT FILES:
- output/hippocampal_interneuron_CREs.bed: GABA neuron CREs
- output/excitatory_neuron_CREs.bed: Excitatory neuron CREs (if available)
- ../data/table_16.txt: Literature CRE-gene correlations

OUTPUT FILES:
- output/hippocampal_interneuron_CREs_genes.tsv: GABA CREs with gene linkage
- output/hippocampal_interneuron_CREs_genes_linked.bed: GABA CREs with genes (for viz)
- output/excitatory_neuron_CREs_genes.tsv: Excitatory CREs with gene linkage
- output/gene_linkage_summary.txt: Summary statistics

Usage:
    python 6_link_CREs_to_genes.py
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime

os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_paper_specific")

print("="*80)
print("LINK CELL-TYPE CREs TO GENES")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# =============================================================================
# Configuration
# =============================================================================

TABLE16_FILE = "../data/table_16.txt"
GABA_BED = "output/hippocampal_interneuron_CREs.bed"
EXCITATORY_BED = "output/excitatory_neuron_CREs.bed"
OUTPUT_DIR = "output"

# Statistical filters
FDR_THRESHOLD = 0.05
PCC_THRESHOLD = 0.2

# =============================================================================
# Step 1: Load CRE BED files
# =============================================================================

print("="*80)
print("STEP 1: Loading cell-type CRE BED files")
print("-"*80)

# Load GABA CREs
if not os.path.exists(GABA_BED):
    print(f"ERROR: GABA BED file not found: {GABA_BED}")
    print("Please run 1_extract_hippocampal_interneuron_CREs.py first")
    exit(1)

gaba_bed = pd.read_csv(GABA_BED, sep='\t', header=None,
                        names=['chr', 'start', 'end', 'cre_id', 'n_cell_types', 'cell_types'])
print(f"Loaded GABA CREs: {len(gaba_bed):,}")

# Load Excitatory CREs (optional)
excit_bed = None
if os.path.exists(EXCITATORY_BED):
    excit_bed = pd.read_csv(EXCITATORY_BED, sep='\t', header=None,
                             names=['chr', 'start', 'end', 'cre_id', 'n_cell_types', 'cell_types'])
    print(f"Loaded Excitatory CREs: {len(excit_bed):,}")
else:
    print(f"Note: Excitatory BED file not found: {EXCITATORY_BED}")
    print("  Run 2a_extract_excitatory_neuron_CREs.py if needed")
print()

# =============================================================================
# Step 2: Load Table 16 (CRE-Gene Correlations)
# =============================================================================

print("="*80)
print("STEP 2: Loading Table 16 (CRE-gene correlations)")
print("-"*80)
print("This may take 1-2 minutes (567 MB file)...")

table16 = pd.read_csv(TABLE16_FILE, sep='\t')
print(f"Loaded {len(table16):,} CRE-gene correlation entries")

print("\nColumns available:")
for col in table16.columns:
    print(f"  {col}")
print()

# =============================================================================
# Step 3: Apply statistical filters
# =============================================================================

print("="*80)
print("STEP 3: Applying statistical filters")
print("-"*80)

print(f"Initial entries: {len(table16):,}")

table16_filtered = table16[
    (table16['FDR'] < FDR_THRESHOLD) &
    (table16['PCC'].abs() > PCC_THRESHOLD)
].copy()

print(f"After FDR < {FDR_THRESHOLD}: {len(table16_filtered):,}")
print(f"After |PCC| > {PCC_THRESHOLD}: {len(table16_filtered):,}")
print()

# =============================================================================
# Step 4: Create lookup by cCRE ID
# =============================================================================

print("="*80)
print("STEP 4: Creating CRE-gene lookup")
print("-"*80)

cre_gene_lookup = table16_filtered.groupby('cCRE1').agg({
    'Gene': lambda x: ','.join(sorted(set(x))),
    'PCC': ['mean', 'max', 'min'],
    'FDR': 'min',
    'Coordinate1': 'first',
    'SubType': lambda x: ','.join(sorted(set(str(s) for s in x if pd.notna(s))))
}).reset_index()

cre_gene_lookup.columns = ['cCRE1', 'Genes', 'PCC_mean', 'PCC_max', 'PCC_min',
                            'FDR_min', 'Coordinate1', 'SubTypes']

cre_gene_lookup['n_genes'] = cre_gene_lookup['Genes'].apply(lambda x: len(x.split(',')))

print(f"Unique CREs with gene links: {len(cre_gene_lookup):,}")
print(f"Average genes per CRE: {cre_gene_lookup['n_genes'].mean():.2f}")
print()

# =============================================================================
# Step 5: Link GABA CREs to genes
# =============================================================================

print("="*80)
print("STEP 5: Linking GABA CREs to genes")
print("-"*80)

gaba_with_genes = gaba_bed.merge(
    cre_gene_lookup,
    left_on='cre_id',
    right_on='cCRE1',
    how='left'
)

gaba_with_links = gaba_with_genes[gaba_with_genes['Genes'].notna()]
gaba_no_links = gaba_with_genes[gaba_with_genes['Genes'].isna()]

print(f"Total GABA CREs: {len(gaba_bed):,}")
print(f"CREs with gene links: {len(gaba_with_links):,} ({100*len(gaba_with_links)/len(gaba_bed):.1f}%)")
print(f"CREs without gene links: {len(gaba_no_links):,}")

unique_genes_gaba = set()
if len(gaba_with_links) > 0:
    unique_genes_gaba = set(','.join(gaba_with_links['Genes'].dropna()).split(','))
    print(f"Unique genes linked: {len(unique_genes_gaba):,}")
    print(f"Sample genes: {', '.join(list(unique_genes_gaba)[:10])}")
print()

# =============================================================================
# Step 6: Link Excitatory CREs to genes (if available)
# =============================================================================

excit_with_links = None
unique_genes_excit = set()

if excit_bed is not None:
    print("="*80)
    print("STEP 6: Linking Excitatory CREs to genes")
    print("-"*80)

    excit_with_genes = excit_bed.merge(
        cre_gene_lookup,
        left_on='cre_id',
        right_on='cCRE1',
        how='left'
    )

    excit_with_links = excit_with_genes[excit_with_genes['Genes'].notna()]
    excit_no_links = excit_with_genes[excit_with_genes['Genes'].isna()]

    print(f"Total Excitatory CREs: {len(excit_bed):,}")
    print(f"CREs with gene links: {len(excit_with_links):,} ({100*len(excit_with_links)/len(excit_bed):.1f}%)")

    if len(excit_with_links) > 0:
        unique_genes_excit = set(','.join(excit_with_links['Genes'].dropna()).split(','))
        print(f"Unique genes linked: {len(unique_genes_excit):,}")
    print()

# =============================================================================
# Step 7: Save output files
# =============================================================================

print("="*80)
print("STEP 7: Saving output files")
print("-"*80)

# Save GABA CREs with genes
gaba_output = os.path.join(OUTPUT_DIR, "hippocampal_interneuron_CREs_genes.tsv")
gaba_with_genes.to_csv(gaba_output, sep='\t', index=False)
print(f"Saved: {gaba_output}")

# Save only linked GABA CREs
gaba_linked_output = os.path.join(OUTPUT_DIR, "hippocampal_interneuron_CREs_genes_linked.tsv")
gaba_with_links.to_csv(gaba_linked_output, sep='\t', index=False)
print(f"Saved: {gaba_linked_output}")

# Create BED file for linked CREs
gaba_bed_linked = gaba_with_links[['chr', 'start', 'end', 'cre_id', 'n_genes', 'Genes']].copy()
gaba_bed_linked_output = os.path.join(OUTPUT_DIR, "hippocampal_interneuron_CREs_genes_linked.bed")
gaba_bed_linked.to_csv(gaba_bed_linked_output, sep='\t', index=False, header=False)
print(f"Saved: {gaba_bed_linked_output}")

# Save Excitatory if available
if excit_with_links is not None:
    excit_output = os.path.join(OUTPUT_DIR, "excitatory_neuron_CREs_genes.tsv")
    excit_with_genes.to_csv(excit_output, sep='\t', index=False)
    print(f"Saved: {excit_output}")

    excit_linked_output = os.path.join(OUTPUT_DIR, "excitatory_neuron_CREs_genes_linked.tsv")
    excit_with_links.to_csv(excit_linked_output, sep='\t', index=False)
    print(f"Saved: {excit_linked_output}")

print()

# =============================================================================
# Step 8: Generate summary report
# =============================================================================

print("="*80)
print("STEP 8: Creating summary report")
print("-"*80)

summary_file = os.path.join(OUTPUT_DIR, "gene_linkage_summary.txt")

with open(summary_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("GENE LINKAGE SUMMARY: CELL-TYPE CREs\n")
    f.write("="*80 + "\n\n")

    f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("NOTE: This pipeline uses OVERLAPPING CRE sets (~60% overlap)\n")
    f.write("For mutually exclusive CREs, see CREs_paper_exclusive/\n\n")

    f.write("DATA SOURCES:\n")
    f.write("-"*80 + "\n")
    f.write(f"Table 16 (CRE-gene correlations): {TABLE16_FILE}\n")
    f.write(f"  Total entries: {len(table16):,}\n")
    f.write(f"  After filtering: {len(table16_filtered):,}\n\n")

    f.write("FILTERS APPLIED:\n")
    f.write("-"*80 + "\n")
    f.write(f"FDR threshold: < {FDR_THRESHOLD}\n")
    f.write(f"|PCC| threshold: > {PCC_THRESHOLD}\n\n")

    f.write("GABA CREs (Hippocampal Interneurons):\n")
    f.write("-"*80 + "\n")
    f.write(f"  Total CREs: {len(gaba_bed):,}\n")
    f.write(f"  CREs with gene links: {len(gaba_with_links):,} ({100*len(gaba_with_links)/len(gaba_bed):.1f}%)\n")
    f.write(f"  Unique genes linked: {len(unique_genes_gaba):,}\n\n")

    if excit_bed is not None:
        f.write("EXCITATORY CREs:\n")
        f.write("-"*80 + "\n")
        f.write(f"  Total CREs: {len(excit_bed):,}\n")
        if excit_with_links is not None:
            f.write(f"  CREs with gene links: {len(excit_with_links):,} ({100*len(excit_with_links)/len(excit_bed):.1f}%)\n")
            f.write(f"  Unique genes linked: {len(unique_genes_excit):,}\n")
        f.write("\n")

    f.write("NEXT STEPS:\n")
    f.write("-"*80 + "\n")
    f.write("1. Run 7_visualize_DA_CREs.py for individual CRE plots with minSig/minFC filtering\n")
    f.write("2. Use *_linked.bed files for visualization with gene annotations\n")
    f.write("\n")
    f.write("="*80 + "\n")

print(f"Saved: {summary_file}")
print()

# =============================================================================
# Completion
# =============================================================================

print("="*80)
print("GENE LINKAGE COMPLETE!")
print("="*80)
print()
print("Summary:")
print(f"  GABA CREs: {len(gaba_bed):,} total, {len(gaba_with_links):,} with gene links")
if excit_bed is not None:
    print(f"  Excitatory CREs: {len(excit_bed):,} total", end="")
    if excit_with_links is not None:
        print(f", {len(excit_with_links):,} with gene links")
    else:
        print()
print()
print(f"Output files saved to: {OUTPUT_DIR}/")
print()
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
