#!/usr/bin/env python3
"""
Step 3: Analyze Differential CREs (Pairwise Comparisons)

This script:
1. Loads deepTools matrices
2. Performs pairwise differential analysis:
   - Nestin-Ctrl vs Nestin-Mut (Within-genotype)
   - Nestin-Ctrl vs Emx1-Mut (Cross-genotype)
   - Nestin-Mut vs Emx1-Mut (Mutant comparison)
3. Identifies significant CREs based on Fold Change and Signal thresholds.
4. Generates visualizations (Volcano plots, Metaprofiles, Individual plots).

Author: Claude Code
Date: December 2024
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import gzip
import argparse
import warnings
import json

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = Path("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_differential_all_celltypes")
OUTPUT_DIR = BASE_DIR / "output"
MATRIX_DIR = OUTPUT_DIR / "matrices"
VIZ_DIR = OUTPUT_DIR / "visualizations"
DIFF_DIR = OUTPUT_DIR / "differential_lists"

# Default Parameters
DEFAULT_MIN_SIGNAL = 2.0
DEFAULT_MIN_FC = 3.0  # User requested 3.0

# Visualization
DPI = 300
COLORS = {
    'Nestin-Ctrl': '#3182bd',  # Blue
    'Nestin-Mut': '#de2d26',   # Red
    'Emx1-Mut': '#ff7f00',     # Orange
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_deeptools_matrix(matrix_file):
    """Load deepTools matrix file and extract signal data."""
    print(f"Loading matrix: {matrix_file}")

    with gzip.open(matrix_file, 'rt') as f:
        header_line = f.readline().strip()

    header_json = header_line.replace('@', '')
    header_data = json.loads(header_json)

    sample_labels = header_data.get('sample_labels', [])
    sample_boundaries = header_data.get('sample_boundaries', [])
    upstream = header_data.get('upstream', [2000])[0]
    downstream = header_data.get('downstream', [2000])[0]
    bin_size = header_data.get('bin size', [50])[0]

    n_bins = int((upstream + downstream) / bin_size)
    
    # Load data
    df = pd.read_csv(matrix_file, sep='\t', skiprows=1, header=None)
    n_regions = len(df)
    
    # Signal data starts at column 6
    signal_data = df.iloc[:, 6:].values.astype(float)

    sample_data = {}
    for i, label in enumerate(sample_labels):
        start_bin = sample_boundaries[i]
        end_bin = sample_boundaries[i + 1]
        sample_data[label] = signal_data[:, start_bin:end_bin]

    positions = np.linspace(-upstream, downstream, n_bins)
    
    return sample_data, positions, n_regions


def load_bed_with_genes(bed_file, gene_links_file):
    """Load BED file and associate with gene information."""
    bed_df = pd.read_csv(bed_file, sep='\t', header=None,
                         names=['chrom', 'start', 'end', 'cre_id', 'score', 'strand'])

    if gene_links_file.exists():
        gene_df = pd.read_csv(gene_links_file, sep='\t')
        cre_genes = gene_df.groupby(['chrom', 'start', 'end'])['gene'].apply(
            lambda x: ', '.join(sorted(set(x)))
        ).reset_index()
        cre_genes.columns = ['chrom', 'start', 'end', 'genes']
        bed_df = bed_df.merge(cre_genes, on=['chrom', 'start', 'end'], how='left')
        bed_df['genes'] = bed_df['genes'].fillna('unknown')
    else:
        bed_df['genes'] = 'unknown'

    return bed_df


def compute_region_stats(sample_data):
    """Compute mean signal for each region (central 500bp)."""
    region_stats = {}
    for sample, data in sample_data.items():
        n_bins = data.shape[1]
        center = n_bins // 2
        central_bins = slice(center - 5, center + 5)
        region_stats[sample] = np.nanmean(data[:, central_bins], axis=1)
    return pd.DataFrame(region_stats)


def perform_differential_analysis(df, comparison_name, cond1, cond2, min_signal, min_fc):
    """
    Identify differential CREs between two conditions.
    
    Returns:
        up_df: CREs significantly HIGHER in cond2 vs cond1
        down_df: CREs significantly LOWER in cond2 vs cond1
    """
    print(f"\nAnalyzing {comparison_name}: {cond1} vs {cond2}")
    
    # Calculate Fold Change (cond2 / cond1)
    # Add small constant to avoid division by zero
    epsilon = 0.01
    
    # Calculate FC
    # If cond2 > cond1 -> FC > 1
    # If cond1 > cond2 -> FC < 1
    
    # We'll use log2FC for symmetry in logic, but filter by linear FC
    df[f'log2fc_{comparison_name}'] = np.log2((df[cond2] + epsilon) / (df[cond1] + epsilon))
    
    # Linear FC for filtering
    # Upregulated in cond2: cond2 / cond1
    fc_up = (df[cond2] + epsilon) / (df[cond1] + epsilon)
    
    # Downregulated in cond2 (Upregulated in cond1): cond1 / cond2
    fc_down = (df[cond1] + epsilon) / (df[cond2] + epsilon)
    
    # Max signal check (at least one condition must be active)
    max_signal = df[[cond1, cond2]].max(axis=1)
    
    # Filter UPREGULATED (Higher in cond2)
    up_mask = (
        (max_signal >= min_signal) &
        (fc_up >= min_fc)
    )
    up_df = df[up_mask].copy()
    
    # Filter DOWNREGULATED (Lower in cond2)
    down_mask = (
        (max_signal >= min_signal) &
        (fc_down >= min_fc)
    )
    down_df = df[down_mask].copy()
    
    print(f"  Total CREs: {len(df)}")
    print(f"  Upregulated in {cond2} (FC >= {min_fc}): {len(up_df)}")
    print(f"  Downregulated in {cond2} (FC <= 1/{min_fc}): {len(down_df)}")
    
    return up_df, down_df


def create_volcano_plot(df, comparison_name, cond1, cond2, up_df, down_df, output_file):
    """Create a scatter plot highlighting differential CREs."""
    plt.figure(figsize=(8, 8))
    
    # Plot all points
    plt.scatter(df[cond1], df[cond2], c='gray', alpha=0.3, s=10, label='Unchanged')
    
    # Highlight Up
    if len(up_df) > 0:
        plt.scatter(up_df[cond1], up_df[cond2], c=COLORS.get(cond2, 'red'), 
                   alpha=0.7, s=20, label=f'Higher in {cond2}')
        
    # Highlight Down
    if len(down_df) > 0:
        plt.scatter(down_df[cond1], down_df[cond2], c=COLORS.get(cond1, 'blue'), 
                   alpha=0.7, s=20, label=f'Higher in {cond1}')
        
    # Diagonal line
    max_val = max(df[cond1].max(), df[cond2].max())
    plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    
    plt.xlabel(f'{cond1} Signal')
    plt.ylabel(f'{cond2} Signal')
    plt.title(f'{comparison_name}\nDifferential CREs')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=DPI)
    plt.close()


def create_metaprofile(sample_data, positions, indices, title, output_file):
    """Create metaprofile for a subset of CREs."""
    plt.figure(figsize=(8, 6))
    
    for sample, color in COLORS.items():
        if sample in sample_data:
            # Extract subset
            data = sample_data[sample][indices, :]
            if len(data) == 0:
                continue
                
            mean_signal = np.nanmean(data, axis=0)
            sem_signal = stats.sem(data, axis=0, nan_policy='omit')
            
            plt.plot(positions, mean_signal, color=color, label=sample, linewidth=2)
            plt.fill_between(positions, mean_signal - sem_signal, mean_signal + sem_signal,
                           color=color, alpha=0.2)
            
    plt.axvline(0, color='gray', linestyle='--', alpha=0.5)
    plt.xlabel('Distance from Center (bp)')
    plt.ylabel('ATAC Signal')
    plt.title(f"{title}\n(n={len(indices)})")
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=DPI)
    plt.close()


def create_individual_plot(sample_data, positions, cre_idx, gene, fc, output_file):
    """Create plot for a single CRE."""
    plt.figure(figsize=(6, 4))
    
    for sample, color in COLORS.items():
        if sample in sample_data:
            signal = sample_data[sample][cre_idx, :]
            plt.plot(positions, signal, color=color, label=sample, linewidth=2)
            
    plt.axvline(0, color='gray', linestyle='--', alpha=0.5)
    plt.title(f"Gene: {gene}\nFC: {fc:.2f}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Differential CRE Analysis')
    parser.add_argument('--min-signal', type=float, default=DEFAULT_MIN_SIGNAL)
    parser.add_argument('--min-fc', type=float, default=DEFAULT_MIN_FC)
    parser.add_argument('--max-individual', type=int, default=20)
    args = parser.parse_args()

    print("=" * 80)
    print("DIFFERENTIAL CRE ANALYSIS")
    print("=" * 80)
    print(f"Parameters: Min Signal={args.min_signal}, Min FC={args.min_fc}")

    # Create directories
    VIZ_DIR.mkdir(parents=True, exist_ok=True)
    DIFF_DIR.mkdir(parents=True, exist_ok=True)

    # Comparisons to perform
    comparisons = [
        ('Nestin_Ctrl_vs_Nestin_Mut', 'Nestin-Ctrl', 'Nestin-Mut'),
        ('Nestin_Ctrl_vs_Emx1_Mut', 'Nestin-Ctrl', 'Emx1-Mut'),
        ('Nestin_Mut_vs_Emx1_Mut', 'Nestin-Mut', 'Emx1-Mut')
    ]

    for cre_type in ['enhancer', 'silencer']:
        matrix_file = MATRIX_DIR / f"matrix_{cre_type}_all.gz"
        bed_file = OUTPUT_DIR / f"{cre_type}_CREs_all_celltypes.bed"
        gene_links_file = OUTPUT_DIR / f"{cre_type}_CREs_all_celltypes_gene_links.tsv"

        if not matrix_file.exists():
            print(f"Skipping {cre_type}: matrix not found")
            continue

        print(f"\nProcessing {cre_type.upper()} CREs...")
        
        # Load data
        sample_data, positions, n_regions = load_deeptools_matrix(matrix_file)
        bed_df = load_bed_with_genes(bed_file, gene_links_file)
        
        # Compute stats
        region_stats = compute_region_stats(sample_data)
        
        # Add metadata to stats for easier handling
        full_df = pd.concat([bed_df, region_stats], axis=1)
        full_df['cre_idx'] = range(len(full_df))

        # Perform comparisons
        for comp_name, cond1, cond2 in comparisons:
            up_df, down_df = perform_differential_analysis(
                full_df, comp_name, cond1, cond2, args.min_signal, args.min_fc
            )
            
            # Save lists
            if len(up_df) > 0:
                up_df.to_csv(DIFF_DIR / f"{cre_type}_{comp_name}_UP_in_{cond2}.tsv", sep='\t', index=False)
            if len(down_df) > 0:
                down_df.to_csv(DIFF_DIR / f"{cre_type}_{comp_name}_DOWN_in_{cond2}.tsv", sep='\t', index=False)
                
            # Visualizations
            create_volcano_plot(
                full_df, comp_name, cond1, cond2, up_df, down_df,
                VIZ_DIR / f"scatter_{cre_type}_{comp_name}.png"
            )
            
            # Metaprofiles of differential sets
            if len(up_df) > 0:
                create_metaprofile(
                    sample_data, positions, up_df['cre_idx'].values,
                    f"{cre_type}: Higher in {cond2} vs {cond1}",
                    VIZ_DIR / f"metaprofile_{cre_type}_{comp_name}_UP_in_{cond2}.png"
                )
            
            if len(down_df) > 0:
                create_metaprofile(
                    sample_data, positions, down_df['cre_idx'].values,
                    f"{cre_type}: Lower in {cond2} vs {cond1}",
                    VIZ_DIR / f"metaprofile_{cre_type}_{comp_name}_DOWN_in_{cond2}.png"
                )
                
            # Individual plots (Top hits)
            # Combine up and down for plotting
            combined_diff = pd.concat([up_df, down_df])
            if len(combined_diff) > 0:
                # Sort by max signal to show most robust peaks
                combined_diff['max_sig'] = combined_diff[[cond1, cond2]].max(axis=1)
                top_hits = combined_diff.sort_values('max_sig', ascending=False).head(args.max_individual)
                
                indiv_dir = VIZ_DIR / f"individual_{cre_type}_{comp_name}"
                indiv_dir.mkdir(exist_ok=True)
                
                for _, row in top_hits.iterrows():
                    gene = str(row['genes']).split(',')[0]
                    # Calculate specific FC for filename
                    fc = (row[cond2] + 0.01) / (row[cond1] + 0.01)
                    if fc < 1: fc = -1/fc
                    
                    create_individual_plot(
                        sample_data, positions, row['cre_idx'], 
                        gene, fc, 
                        indiv_dir / f"{gene}_fc{fc:.1f}.png"
                    )

    print("\nAnalysis Complete.")

if __name__ == "__main__":
    main()
