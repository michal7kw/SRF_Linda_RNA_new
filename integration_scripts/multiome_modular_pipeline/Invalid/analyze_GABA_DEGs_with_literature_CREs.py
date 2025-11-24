#!/usr/bin/env python3
"""
Analyze GABA DEGs with Literature CREs

Purpose:
- Load significant GABA DEGs (Nestin and Emx1 separately)
- Find literature-validated CREs for these genes
- Check if CRE accessibility changes explain gene expression changes
- Extract ATAC signal from BigWig files for visual comparison

Analysis approach:
1. Load GABA DEGs (padj < 0.05, |log2FC| > 0.25)
2. Search Table 16 for CRE-gene links (hippocampal/GABAergic filtered)
3. Check DA peak overlaps
4. Extract BigWig signal at CREs (Ctrl vs Mut)
5. Create hypothesis-generating visualizations


INPUTS:
- /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/dge_res/geno_spec_cond_comp/dge_Nestin_GABA_mut_vs_ctrl.csv
	- (Differential expression results for Nestin genotype GABA cells)
- /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/dge_res/geno_spec_cond_comp/dge_Emx1_GABA_mut_vs_ctrl.csv
	- (Differential expression results for Emx1 genotype GABA cells)
- data/table_16.txt
	- (Literature CRE-gene correlations)
- ../signac_results_L1/celltype_results/filtered/GABA_Nestin_DA_peaks.csv
	- (Differential accessibility peaks for Nestin genotype GABA cells)
- ../signac_results_L1/celltype_results/filtered/GABA_Emx1_DA_peaks.csv
	- (Differential accessibility peaks for Emx1 genotype GABA cells)
- ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Nestin-Ctrl.bw
	- (BigWig signal file for Nestin control)
- ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Nestin-Mut.bw
	- (BigWig signal file for Nestin mutant)
- ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Emx1-Ctrl.bw
	- (BigWig signal file for Emx1 control)
- ../signac_results_L1/bigwig_tracks_L1/by_celltype/GABA_Emx1-Mut.bw
	- (BigWig signal file for Emx1 mutant)

OUTPUTS:
- output/GABA_DEG_analysis/DEG_CRE_links_Nestin.tsv
	- (Detailed gene-CRE-DA links for Nestin genotype)
- output/GABA_DEG_analysis/DEG_CRE_links_Emx1.tsv
	- (Detailed gene-CRE-DA links for Emx1 genotype)
- output/GABA_DEG_analysis/gene_summary_Nestin.tsv
	- (Per-gene summary for Nestin genotype)
- output/GABA_DEG_analysis/gene_summary_Emx1.tsv
	- (Per-gene summary for Emx1 genotype)
- output/GABA_DEG_analysis/analysis_report.txt
	- (Comprehensive analysis report)

"""

import pandas as pd
import numpy as np
import os
from datetime import datetime

print("="*80)
print("ANALYZE GABA DEGs WITH LITERATURE CREs")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# ============================================================================
# Configuration
# ============================================================================
DEG_BASE = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/DEGs_cell_type_L1FC_0_25/dge_res/geno_spec_cond_comp"
DA_BASE = "../signac_results_L1/celltype_results/filtered"
BIGWIG_BASE = "../signac_results_L1/bigwig_tracks_L1/by_celltype"
LITERATURE_CRE_DIR = "data"
OUTPUT_DIR = "output/GABA_DEG_analysis"

# Statistical thresholds
PADJ_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 0.25

# Cell type and genotypes
CELL_TYPE = "GABA"
GENOTYPES = ["Nestin", "Emx1"]

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# STEP 1: Load and filter GABA DEGs
# ============================================================================
print("STEP 1: Loading GABA DEGs...")
print("-"*80)

deg_data = {}

for genotype in GENOTYPES:
    deg_file = os.path.join(DEG_BASE, f"dge_{genotype}_GABA_mut_vs_ctrl.csv")

    if not os.path.exists(deg_file):
        print(f"ERROR: {deg_file} not found!")
        continue

    # Load DEGs
    deg = pd.read_csv(deg_file)
    print(f"\n{genotype}:")
    print(f"  Total genes tested: {len(deg):,}")

    # Filter for significance
    deg_sig = deg[
        (deg['pvals_adj'] < PADJ_THRESHOLD) &
        (deg['logfoldchanges'].abs() > LOG2FC_THRESHOLD)
    ].copy()

    # Add direction
    deg_sig['direction'] = deg_sig['logfoldchanges'].apply(
        lambda x: 'Up' if x > 0 else 'Down'
    )

    # Sort by significance
    deg_sig = deg_sig.sort_values('pvals_adj')

    deg_data[genotype] = {
        'all': deg,
        'significant': deg_sig
    }

    print(f"  Significant DEGs: {len(deg_sig):,}")
    print(f"    Up-regulated: {(deg_sig['direction'] == 'Up').sum():,}")
    print(f"    Down-regulated: {(deg_sig['direction'] == 'Down').sum():,}")
    print(f"    Top 5 by significance:")
    for _, row in deg_sig.head(5).iterrows():
        print(f"      {row['names']}: log2FC={row['logfoldchanges']:.3f}, "
              f"padj={row['pvals_adj']:.2e}, {row['direction']}")

# ============================================================================
# STEP 2: Load literature CRE-gene correlations
# ============================================================================
print(f"\n{'='*80}")
print("STEP 2: Loading literature CRE-gene correlations...")
print("-"*80)

print("Loading Table 16 (may take 1-2 minutes)...")
table16 = pd.read_csv(os.path.join(LITERATURE_CRE_DIR, 'table_16.txt'), sep='\t')
print(f"✓ Loaded {len(table16):,} CRE-gene correlations")

# Filter for hippocampal/GABAergic cell types
hippocampal_keywords = [
    'CA1', 'CA2', 'CA3', 'DG', 'DGNBL', 'GRC',
    'LAMP5', 'LAMP', 'VIP', 'SST', 'PV', 'PVGA', 'SSTGA', 'VIPGA', 'LAMGA',
    'GABA', 'INH', 'interneuron'
]

def is_hippocampal_or_gaba(subtype):
    if pd.isna(subtype):
        return False
    subtype_str = str(subtype).upper()
    return any(keyword.upper() in subtype_str for keyword in hippocampal_keywords)

table16['is_hippo_gaba'] = table16['SubType'].apply(is_hippocampal_or_gaba)
table16_filtered = table16[table16['is_hippo_gaba']].copy()

print(f"Filtered to hippocampal/GABAergic: {len(table16_filtered):,} links")
print(f"Unique genes: {table16_filtered['Gene'].nunique():,}")

# Apply statistical filters (FDR < 0.05, |PCC| > 0.2)
table16_sig = table16_filtered[
    (table16_filtered['FDR'] < 0.05) &
    (table16_filtered['PCC'].abs() > 0.2)
].copy()

print(f"After statistical filtering: {len(table16_sig):,} significant links")
print(f"Unique genes with CRE links: {table16_sig['Gene'].nunique():,}")

# Parse CRE coordinates
def parse_coordinate(coord_str):
    try:
        parts = coord_str.split('_')
        if len(parts) == 3:
            return pd.Series({
                'cre_chr': parts[0],
                'cre_start': int(parts[1]),
                'cre_end': int(parts[2])
            })
    except:
        pass
    return pd.Series({'cre_chr': None, 'cre_start': None, 'cre_end': None})

coords = table16_sig['Coordinate1'].apply(parse_coordinate)
table16_sig = pd.concat([table16_sig, coords], axis=1)
table16_sig = table16_sig.dropna(subset=['cre_chr', 'cre_start', 'cre_end'])

print(f"Successfully parsed coordinates: {len(table16_sig):,} CRE-gene links")

# ============================================================================
# STEP 3: Link DEGs to literature CREs
# ============================================================================
print(f"\n{'='*80}")
print("STEP 3: Linking DEGs to literature CREs...")
print("-"*80)

deg_with_cres = {}

for genotype in GENOTYPES:
    if genotype not in deg_data:
        continue

    print(f"\n{genotype}:")
    deg_sig = deg_data[genotype]['significant']

    # Find DEGs with literature CREs
    deg_genes = set(deg_sig['names'].tolist())
    lit_genes = set(table16_sig['Gene'].tolist())
    genes_with_cres = deg_genes.intersection(lit_genes)

    print(f"  DEGs with literature CREs: {len(genes_with_cres)} / {len(deg_sig)}")

    # Create merged dataset
    deg_cre_links = []

    for gene in genes_with_cres:
        # Get DEG info
        deg_info = deg_sig[deg_sig['names'] == gene].iloc[0]

        # Get CRE links
        cre_links = table16_sig[table16_sig['Gene'] == gene]

        for _, cre in cre_links.iterrows():
            deg_cre_links.append({
                'gene': gene,
                'gene_log2fc': deg_info['logfoldchanges'],
                'gene_padj': deg_info['pvals_adj'],
                'gene_direction': deg_info['direction'],
                'cre_id': cre['cCRE1'],
                'cre_chr': cre['cre_chr'],
                'cre_start': cre['cre_start'],
                'cre_end': cre['cre_end'],
                'cre_pcc': cre['PCC'],
                'cre_fdr': cre['FDR'],
                'cre_subtype': cre['SubType']
            })

    deg_cre_df = pd.DataFrame(deg_cre_links)
    deg_with_cres[genotype] = deg_cre_df

    print(f"  Total gene-CRE links: {len(deg_cre_df)}")

    if len(genes_with_cres) > 0:
        print(f"  Top genes with most CREs:")
        gene_cre_counts = deg_cre_df['gene'].value_counts().head(10)
        for gene, count in gene_cre_counts.items():
            gene_info = deg_sig[deg_sig['names'] == gene].iloc[0]
            print(f"    {gene}: {count} CREs (log2FC={gene_info['logfoldchanges']:.3f}, "
                  f"{gene_info['direction']})")

# ============================================================================
# STEP 4: Check DA peak overlaps
# ============================================================================
print(f"\n{'='*80}")
print("STEP 4: Checking DA peak overlaps...")
print("-"*80)

def find_overlaps(cres_df, peaks_df, min_overlap=1):
    """Find overlaps between CREs and peaks"""
    overlaps = []

    for _, cre in cres_df.iterrows():
        cre_chr = cre['cre_chr']
        cre_start = cre['cre_start']
        cre_end = cre['cre_end']

        # Find peaks on same chromosome
        peaks_chr = peaks_df[peaks_df['chr'] == cre_chr]

        for _, peak in peaks_chr.iterrows():
            overlap_start = max(cre_start, peak['start'])
            overlap_end = min(cre_end, peak['end'])
            overlap_size = overlap_end - overlap_start

            if overlap_size >= min_overlap:
                overlaps.append({
                    'gene': cre['gene'],
                    'cre_id': cre['cre_id'],
                    'peak': peak['peak'],
                    'peak_log2fc': peak['avg_log2FC'],
                    'peak_padj': peak['p_val_adj'],
                    'peak_direction': peak['direction'],
                    'overlap_size': overlap_size
                })

    return pd.DataFrame(overlaps)

da_overlaps = {}

for genotype in GENOTYPES:
    if genotype not in deg_with_cres:
        continue

    print(f"\n{genotype}:")

    # Load DA peaks
    da_file = os.path.join(DA_BASE, f"GABA_{genotype}_DA_peaks.csv")

    if not os.path.exists(da_file):
        print(f"  DA peaks file not found: {da_file}")
        continue

    da = pd.read_csv(da_file)
    da[['chr', 'start', 'end']] = da['peak'].str.split('-', expand=True)
    da['start'] = pd.to_numeric(da['start'])
    da['end'] = pd.to_numeric(da['end'])

    print(f"  DA peaks: {len(da):,}")

    # Find overlaps
    overlaps = find_overlaps(deg_with_cres[genotype], da)
    da_overlaps[genotype] = overlaps

    print(f"  CRE-DA peak overlaps: {len(overlaps)}")

    if len(overlaps) > 0:
        print(f"  Genes with DA peak overlaps:")
        for gene in overlaps['gene'].unique():
            gene_overlaps = overlaps[overlaps['gene'] == gene]
            gene_info = deg_with_cres[genotype][deg_with_cres[genotype]['gene'] == gene].iloc[0]
            print(f"    {gene}: {len(gene_overlaps)} DA peaks, "
                  f"gene {gene_info['gene_direction']} (log2FC={gene_info['gene_log2fc']:.3f})")

# ============================================================================
# STEP 5: Check BigWig availability
# ============================================================================
print(f"\n{'='*80}")
print("STEP 5: Checking BigWig files...")
print("-"*80)

bigwig_files = {}
for genotype in GENOTYPES:
    for condition in ['Ctrl', 'Mut']:
        bw_file = os.path.join(BIGWIG_BASE, f"GABA_{genotype}-{condition}.bw")
        if os.path.exists(bw_file):
            bigwig_files[f"{genotype}_{condition}"] = bw_file
            print(f"✓ Found: GABA_{genotype}-{condition}.bw")
        else:
            print(f"✗ Missing: GABA_{genotype}-{condition}.bw")

if len(bigwig_files) == 0:
    print("\n⚠️  No BigWig files found. Will skip signal extraction.")
    print("    Analysis will continue with DA peaks only.")
    EXTRACT_BIGWIG = False
else:
    print(f"\nFound {len(bigwig_files)} BigWig files. Will extract signals.")
    EXTRACT_BIGWIG = True

# ============================================================================
# STEP 6: Create summary reports
# ============================================================================
print(f"\n{'='*80}")
print("STEP 6: Creating summary reports...")
print("-"*80)

for genotype in GENOTYPES:
    if genotype not in deg_with_cres:
        continue

    print(f"\n{genotype}:")

    deg_cre_df = deg_with_cres[genotype]
    overlaps_df = da_overlaps.get(genotype, pd.DataFrame())

    # Merge DEG-CRE with DA overlaps
    if len(overlaps_df) > 0:
        deg_cre_df = deg_cre_df.merge(
            overlaps_df[['gene', 'cre_id', 'peak_log2fc', 'peak_direction']],
            on=['gene', 'cre_id'],
            how='left'
        )
        deg_cre_df['has_da_peak'] = ~deg_cre_df['peak_direction'].isna()
    else:
        deg_cre_df['has_da_peak'] = False
        deg_cre_df['peak_log2fc'] = np.nan
        deg_cre_df['peak_direction'] = None

    # Check consistency
    def check_consistency(row):
        if not row['has_da_peak']:
            return 'NO_DA_PEAK'

        gene_dir = row['gene_direction']
        peak_dir = row['peak_direction']

        if gene_dir == 'Down' and peak_dir == 'Less_Accessible':
            return 'CONSISTENT'
        elif gene_dir == 'Up' and peak_dir == 'More_Accessible':
            return 'CONSISTENT'
        else:
            return 'INCONSISTENT'

    deg_cre_df['consistency'] = deg_cre_df.apply(check_consistency, axis=1)

    # Summary statistics
    print(f"  Total gene-CRE pairs: {len(deg_cre_df)}")
    print(f"  With DA peaks: {deg_cre_df['has_da_peak'].sum()}")
    print(f"    Consistent: {(deg_cre_df['consistency'] == 'CONSISTENT').sum()}")
    print(f"    Inconsistent: {(deg_cre_df['consistency'] == 'INCONSISTENT').sum()}")
    print(f"  Without DA peaks: {(deg_cre_df['consistency'] == 'NO_DA_PEAK').sum()}")

    # Save detailed report
    output_file = os.path.join(OUTPUT_DIR, f"DEG_CRE_links_{genotype}.tsv")
    deg_cre_df.to_csv(output_file, sep='\t', index=False)
    print(f"  ✓ Saved: {output_file}")

    # Create gene-level summary
    gene_summary = []
    for gene in deg_cre_df['gene'].unique():
        gene_data = deg_cre_df[deg_cre_df['gene'] == gene]
        gene_summary.append({
            'gene': gene,
            'gene_log2fc': gene_data['gene_log2fc'].iloc[0],
            'gene_padj': gene_data['gene_padj'].iloc[0],
            'gene_direction': gene_data['gene_direction'].iloc[0],
            'n_cres': len(gene_data),
            'n_cres_with_da': gene_data['has_da_peak'].sum(),
            'n_consistent': (gene_data['consistency'] == 'CONSISTENT').sum(),
            'n_inconsistent': (gene_data['consistency'] == 'INCONSISTENT').sum(),
            'cre_list': ','.join(gene_data['cre_id'].tolist())
        })

    gene_summary_df = pd.DataFrame(gene_summary)
    gene_summary_df = gene_summary_df.sort_values('gene_padj')

    summary_file = os.path.join(OUTPUT_DIR, f"gene_summary_{genotype}.tsv")
    gene_summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"  ✓ Saved: {summary_file}")

    # Show top genes with CRE support
    print(f"\n  Top genes with CRE support:")
    for _, row in gene_summary_df.head(10).iterrows():
        support = ""
        if row['n_consistent'] > 0:
            support = f"✓ {row['n_consistent']}/{row['n_cres_with_da']} consistent DA"
        elif row['n_cres_with_da'] > 0:
            support = f"✗ {row['n_inconsistent']}/{row['n_cres_with_da']} inconsistent DA"
        else:
            support = f"○ {row['n_cres']} CREs (no DA)"

        print(f"    {row['gene']}: log2FC={row['gene_log2fc']:.3f}, "
              f"{row['gene_direction']}, {support}")

# ============================================================================
# STEP 7: Create analysis report
# ============================================================================
print(f"\n{'='*80}")
print("STEP 7: Creating analysis report...")
print("-"*80)

report_file = os.path.join(OUTPUT_DIR, "analysis_report.txt")
with open(report_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("GABA DEG ANALYSIS WITH LITERATURE CREs\n")
    f.write("="*80 + "\n\n")
    f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("CONFIGURATION:\n")
    f.write(f"  Cell type: {CELL_TYPE}\n")
    f.write(f"  Genotypes: {', '.join(GENOTYPES)}\n")
    f.write(f"  DEG thresholds: padj < {PADJ_THRESHOLD}, |log2FC| > {LOG2FC_THRESHOLD}\n")
    f.write(f"  CRE source: Table 16 (hippocampal/GABAergic filtered)\n")
    f.write(f"  CRE filters: FDR < 0.05, |PCC| > 0.2\n\n")

    f.write("RESULTS:\n")
    f.write("-"*80 + "\n\n")

    for genotype in GENOTYPES:
        if genotype not in deg_data:
            continue

        f.write(f"{genotype}:\n")

        deg_sig = deg_data[genotype]['significant']
        f.write(f"  Total significant DEGs: {len(deg_sig):,}\n")
        f.write(f"    Up-regulated: {(deg_sig['direction'] == 'Up').sum():,}\n")
        f.write(f"    Down-regulated: {(deg_sig['direction'] == 'Down').sum():,}\n\n")

        if genotype in deg_with_cres:
            deg_cre_df = deg_with_cres[genotype]
            n_genes_with_cres = deg_cre_df['gene'].nunique()

            f.write(f"  DEGs with literature CREs: {n_genes_with_cres}\n")
            f.write(f"  Total gene-CRE pairs: {len(deg_cre_df)}\n\n")

            if 'consistency' in deg_cre_df.columns:
                f.write(f"  CRE-DA peak consistency:\n")
                f.write(f"    Consistent: {(deg_cre_df['consistency'] == 'CONSISTENT').sum()}\n")
                f.write(f"    Inconsistent: {(deg_cre_df['consistency'] == 'INCONSISTENT').sum()}\n")
                f.write(f"    No DA peak: {(deg_cre_df['consistency'] == 'NO_DA_PEAK').sum()}\n\n")

        f.write("\n")

    f.write("="*80 + "\n")
    f.write("OUTPUT FILES:\n")
    f.write("="*80 + "\n")
    f.write(f"  DEG_CRE_links_{{genotype}}.tsv - Detailed gene-CRE-DA links\n")
    f.write(f"  gene_summary_{{genotype}}.tsv - Per-gene summary\n")
    f.write(f"  analysis_report.txt - This file\n\n")

    f.write("NEXT STEPS:\n")
    f.write("-"*80 + "\n")
    f.write("1. Review genes with CONSISTENT CRE support\n")
    f.write("2. Extract BigWig signals for visual validation (see separate script)\n")
    f.write("3. Generate coverage plots for top candidates\n")
    f.write("4. Prioritize for experimental validation\n")

print(f"✓ Saved: {report_file}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print(f"\n{'='*80}")
print("ANALYSIS COMPLETE!")
print("="*80)

for genotype in GENOTYPES:
    if genotype in deg_with_cres:
        deg_cre_df = deg_with_cres[genotype]
        n_genes = deg_cre_df['gene'].nunique()
        n_consistent = (deg_cre_df['consistency'] == 'CONSISTENT').sum() if 'consistency' in deg_cre_df.columns else 0

        print(f"\n{genotype}:")
        print(f"  DEGs with literature CREs: {n_genes}")
        print(f"  Gene-CRE pairs with consistent DA: {n_consistent}")

print(f"\nOutput directory: {OUTPUT_DIR}/")
print(f"  - DEG_CRE_links_{{genotype}}.tsv")
print(f"  - gene_summary_{{genotype}}.tsv")
print(f"  - analysis_report.txt")

if EXTRACT_BIGWIG:
    print(f"\nBigWig files found: {len(bigwig_files)}")
    print("  Run the BigWig extraction script next for signal visualization")
else:
    print("\n⚠️  BigWig files not found - signal extraction skipped")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
