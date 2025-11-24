# %% [markdown]
# # Individual Sample Processing Script
# 
# This script loads raw count data (e.g., from Cell Ranger) for each sample individually,
# performs standard preprocessing, normalization, dimensionality reduction,
# and clustering on each sample separately, then saves each processed sample
# as a .h5ad file in its respective sample directory.

# %%
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import anndata as ad
from tqdm import tqdm
import random
from pathlib import Path
import seaborn as sns

# %%
# --- Configuration ---

# Set seeds for reproducibility
random_seed = 0
np.random.seed(random_seed)
random.seed(random_seed)

PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = os.path.join(PROJECT_DIR, "combine_data")
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

# Define Input Directory for Raw Data
# Assumes Cell Ranger structure: INPUT_RAW_DATA_BASE/SAMPLE_NAME/outs/filtered_feature_bc_matrix/
INPUT_RAW_DATA_BASE = os.path.join(PROJECT_DIR, "cellranger_final_count_data")

# Define Samples to Process
SAMPLES = ["Emx1_Ctrl", "Emx1_Mut", "Nestin_Ctrl", "Nestin_Mut"]

# Define QC Parameters
MIN_GENES_PER_CELL = 500 
MIN_CELLS_PER_GENE = 6
MAX_MITO_PERCENT = 20

# Define Processing Parameters
N_TOP_HVG = 3000
N_PCS = 50
LEIDEN_RESOLUTION = 0.4

# %% [markdown]
# ## Process Each Sample Individually

# %%
def process_single_sample(sample_name):
    """
    Process a single sample: load, QC, filter, normalize, scale, 
    dimensionality reduction, clustering, and save.
    """
    print(f"\n{'='*60}")
    print(f"Processing sample: {sample_name}")
    print(f"{'='*60}")
    
    # Path to the directory containing matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz
    sample_dir_path = os.path.join(INPUT_RAW_DATA_BASE, sample_name, "outs", "filtered_feature_bc_matrix")
    print(f"Loading {sample_name} from directory: {sample_dir_path}")
    
    if not os.path.isdir(sample_dir_path):
        print(f"Warning: Directory not found for sample {sample_name} at {sample_dir_path}. Skipping.")
        return None
    
    try:
        # Use read_10x_mtx for Matrix Market format
        adata = sc.read_10x_mtx(sample_dir_path, cache=True)
        adata.var_names_make_unique()  # Ensure unique gene names
        
        # Add sample information
        adata.obs['sample'] = sample_name
        # Add condition information
        if "Ctrl" in sample_name:
            adata.obs['condition'] = 'Control'
        elif "Mut" in sample_name:
            adata.obs['condition'] = 'Mutant'
        else:
            adata.obs['condition'] = 'Unknown'

        # Add genotype information
        if "Emx1" in sample_name:
            adata.obs['genotype'] = 'Emx1'
        elif "Nestin" in sample_name:
            adata.obs['genotype'] = 'Nestin'
        else:
            adata.obs['genotype'] = 'Unknown'

        print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes for {sample_name}")
        
        # Make a copy of raw counts
        adata.raw = adata.copy()
        
        # --- Quality Control and Filtering ---
        print(f"\nPerforming QC and Filtering for {sample_name}...")
        
        # Calculate mitochondrial gene percentage
        adata.var['mt'] = adata.var_names.str.startswith('mt-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        
        # Store data before filtering for comparison
        adata_before = adata.copy()
        
        # Apply filters
        print(f"Initial cell count: {adata.n_obs}")
        print(f"Initial gene count: {adata.n_vars}")
        
        # Calculate the 95th percentile for max_genes
        max_genes_percentile = np.percentile(adata.obs['n_genes_by_counts'], 95)
        print(f"Calculated 95th percentile for max genes: {max_genes_percentile:.0f}")
        
        # 1. Filter cells based on min number of genes
        n_obs_before = adata.n_obs
        adata = adata[np.array(adata.obs.n_genes_by_counts) > MIN_GENES_PER_CELL]
        print(f"Cells after min_genes filter (> {MIN_GENES_PER_CELL}): {adata.n_obs} (removed {n_obs_before - adata.n_obs})")
        
        # 2. Filter cells based on max number of genes (percentile)
        n_obs_before = adata.n_obs
        adata = adata[np.array(adata.obs.n_genes_by_counts) < int(max_genes_percentile)]
        print(f"Cells after max_genes filter (< {max_genes_percentile:.0f}): {adata.n_obs} (removed {n_obs_before - adata.n_obs})")
        
        # 3. Filter cells based on mitochondrial content
        n_obs_before = adata.n_obs
        adata = adata[np.array(adata.obs.pct_counts_mt) < MAX_MITO_PERCENT]
        print(f"Cells after mito filter (< {MAX_MITO_PERCENT}%): {adata.n_obs} (removed {n_obs_before - adata.n_obs})")
        
        # 4. Filter genes based on minimum number of cells expressing them
        n_vars_before = adata.n_vars
        sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
        print(f"Genes after min_cells filter ({MIN_CELLS_PER_GENE}): {adata.n_vars} (removed {n_vars_before - adata.n_vars})")
        
        print("QC and Filtering complete.")
        
        # --- Create QC Visualization ---
        print(f"\nCreating quality control visualizations for {sample_name}...")
        
        # Create sample-specific output directory
        sample_output_dir = os.path.join(INPUT_RAW_DATA_BASE, sample_name)
        os.makedirs(sample_output_dir, exist_ok=True)
        
        # Create plots directory
        plot_dir = os.path.join(sample_output_dir, "plots")
        os.makedirs(plot_dir, exist_ok=True)
        
        # Create QC comparison plot
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle(f'Quality Control: {sample_name} - Before vs After Filtering', 
                     fontsize=14, fontweight='bold')
        
        # Colors
        before_color = '#E74C3C'  # Red
        after_color = '#2ECC71'   # Green
        threshold_color = '#F39C12'  # Orange
        
        # Genes per cell
        axes[0,0].hist(adata_before.obs['n_genes_by_counts'], bins=50, alpha=0.7, 
                       color=before_color, label='Before', density=True)
        axes[0,0].hist(adata.obs['n_genes_by_counts'], bins=50, alpha=0.7, 
                       color=after_color, label='After', density=True)
        axes[0,0].axvline(MIN_GENES_PER_CELL, color=threshold_color, linestyle='--', linewidth=2)
        axes[0,0].axvline(max_genes_percentile, color=threshold_color, linestyle='--', linewidth=2)
        axes[0,0].set_title('Genes per Cell')
        axes[0,0].legend()
        axes[0,0].grid(True, alpha=0.3)
        
        # Total UMI counts
        axes[0,1].hist(adata_before.obs['total_counts'], bins=50, alpha=0.7, 
                       color=before_color, label='Before', density=True)
        axes[0,1].hist(adata.obs['total_counts'], bins=50, alpha=0.7, 
                       color=after_color, label='After', density=True)
        axes[0,1].set_title('Total UMI Counts')
        axes[0,1].legend()
        axes[0,1].grid(True, alpha=0.3)
        
        # Mitochondrial percentage
        axes[0,2].hist(adata_before.obs['pct_counts_mt'], bins=50, alpha=0.7, 
                       color=before_color, label='Before', density=True)
        axes[0,2].hist(adata.obs['pct_counts_mt'], bins=50, alpha=0.7, 
                       color=after_color, label='After', density=True)
        axes[0,2].axvline(MAX_MITO_PERCENT, color=threshold_color, linestyle='--', linewidth=2)
        axes[0,2].set_title('Mitochondrial %')
        axes[0,2].legend()
        axes[0,2].grid(True, alpha=0.3)
        
        # Scatter plots
        scatter1 = axes[1,0].scatter(adata_before.obs['total_counts'], 
                                    adata_before.obs['n_genes_by_counts'], 
                                    c=adata_before.obs['pct_counts_mt'], 
                                    cmap='Reds', alpha=0.6, s=8)
        axes[1,0].set_xlabel('Total UMI Counts')
        axes[1,0].set_ylabel('Genes per Cell')
        axes[1,0].set_title('Before: Genes vs UMI (colored by mito %)')
        plt.colorbar(scatter1, ax=axes[1,0], label='Mito %')
        
        scatter2 = axes[1,1].scatter(adata.obs['total_counts'], 
                                    adata.obs['n_genes_by_counts'], 
                                    c=adata.obs['pct_counts_mt'], 
                                    cmap='Greens', alpha=0.6, s=8)
        axes[1,1].set_xlabel('Total UMI Counts')
        axes[1,1].set_ylabel('Genes per Cell')
        axes[1,1].set_title('After: Genes vs UMI (colored by mito %)')
        plt.colorbar(scatter2, ax=axes[1,1], label='Mito %')
        
        # Summary stats
        axes[1,2].axis('off')
        summary_text = f"""
        FILTERING SUMMARY
        
        Cells: {adata_before.n_obs:,} → {adata.n_obs:,}
        Removed: {adata_before.n_obs - adata.n_obs:,} ({((adata_before.n_obs - adata.n_obs)/adata_before.n_obs*100):.1f}%)
        
        Genes: {adata_before.n_vars:,} → {adata.n_vars:,}
        Removed: {adata_before.n_vars - adata.n_vars:,} ({((adata_before.n_vars - adata.n_vars)/adata_before.n_vars*100):.1f}%)
        
        THRESHOLDS USED:
        • Min genes/cell: {MIN_GENES_PER_CELL}
        • Max genes/cell: {max_genes_percentile:.0f}
        • Max mito %: {MAX_MITO_PERCENT}
        • Min cells/gene: {MIN_CELLS_PER_GENE}
        """
        axes[1,2].text(0.1, 0.9, summary_text, transform=axes[1,2].transAxes, 
                       fontsize=9, verticalalignment='top',
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.8))
        
        plt.tight_layout()
        qc_plot_path = os.path.join(plot_dir, f'{sample_name}_qc_before_after_comparison.png')
        plt.savefig(qc_plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"QC visualization saved to {qc_plot_path}")
        
        # --- Identify Highly Variable Genes ---
        print(f"\nIdentifying top {N_TOP_HVG} highly variable genes for {sample_name}...")
        sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_HVG, flavor='seurat_v3')
        
        # Save HVG plot
        sc.settings.figdir = plot_dir
        sc.pl.highly_variable_genes(adata, save=f'_{sample_name}_hvg.png', show=False)
        
        # --- Normalization and Scaling ---
        print(f"\nNormalizing and Scaling data for {sample_name}...")
        
        # Normalize counts per cell (library size normalization)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Scale data to unit variance and zero mean
        sc.pp.scale(adata, max_value=10)
        
        print("Normalization and Scaling complete.")
        
        # --- Dimensionality Reduction ---
        print(f"\nPerforming Dimensionality Reduction (PCA and UMAP) for {sample_name}...")
        
        # PCA
        sc.tl.pca(adata, n_comps=N_PCS, use_highly_variable=True, svd_solver='arpack')
        
        # Save PCA variance plot
        sc.pl.pca_variance_ratio(adata, log=True, n_pcs=N_PCS, 
                                save=f'_{sample_name}_pca_variance.png', show=False)
        
        # Compute neighborhood graph
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=N_PCS)
        
        # UMAP
        sc.tl.umap(adata)
        
        print("Dimensionality Reduction complete.")
        
        # --- Clustering ---
        print(f"\nPerforming Clustering (Leiden) for {sample_name}...")
        
        # Leiden clustering
        sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, key_added=f'leiden_{LEIDEN_RESOLUTION}')
        
        print("Clustering complete.")
        
        # --- Visualization ---
        print(f"\nGenerating visualizations for {sample_name}...")
        
        # Plot UMAP colored by different features
        umap_plots = {
            'sample': 'sample',
            'condition': 'condition',
            'genotype': 'genotype',
            f'leiden_{LEIDEN_RESOLUTION}': f'leiden_{LEIDEN_RESOLUTION}'
        }
        
        for key, column in umap_plots.items():
            if column in adata.obs:
                sc.pl.umap(adata, color=column, legend_loc='right margin',
                          save=f'_{sample_name}_umap_{key}.png', show=False,
                          title=f'{sample_name} UMAP colored by {key.replace("_", " ").title()}')
        
        # --- Save Processed Data ---
        output_file = os.path.join(sample_output_dir, f'{sample_name}_processed.h5ad')
        print(f"\nSaving processed dataset to {output_file}")
        
        try:
            adata.write(Path(output_file))
            print(f"Successfully saved {sample_name} processed data.")
        except Exception as e:
            print(f"Error saving {sample_name} processed data: {e}")
            return None
        
        # --- Generate Sample Summary ---
        print(f"\nGenerating sample summary for {sample_name}...")
        
        summary_data = {
            'sample': sample_name,
            'n_cells_before_filtering': adata_before.n_obs,
            'n_genes_before_filtering': adata_before.n_vars,
            'n_cells_after_filtering': adata.n_obs,
            'n_genes_after_filtering': adata.n_vars,
            'cells_removed': adata_before.n_obs - adata.n_obs,
            'genes_removed': adata_before.n_vars - adata.n_vars,
            'percent_cells_removed': ((adata_before.n_obs - adata.n_obs)/adata_before.n_obs*100),
            'percent_genes_removed': ((adata_before.n_vars - adata.n_vars)/adata_before.n_vars*100),
            'median_genes_per_cell': adata.obs['n_genes_by_counts'].median(),
            'median_umi_per_cell': adata.obs['total_counts'].median(),
            'n_clusters': len(adata.obs[f'leiden_{LEIDEN_RESOLUTION}'].unique())
        }
        
        summary_df = pd.DataFrame([summary_data])
        summary_csv_path = os.path.join(sample_output_dir, f'{sample_name}_processing_summary.csv')
        summary_df.to_csv(summary_csv_path, index=False)
        print(f"Sample summary saved to {summary_csv_path}")
        
        print(f"\n✅ Successfully processed {sample_name}")
        print(f"📊 Final dataset: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
        print(f"🗂️  Files saved in: {sample_output_dir}")
        
        return adata
        
    except Exception as e:
        print(f"Error processing {sample_name}: {e}")
        return None

# %%
# Process all samples
print("Starting individual sample processing...")
print(f"Samples to process: {SAMPLES}")

processed_samples = {}
all_summaries = []

for sample in SAMPLES:
    adata_sample = process_single_sample(sample)
    if adata_sample is not None:
        processed_samples[sample] = adata_sample
        
        # Collect summary for combined report
        sample_output_dir = os.path.join(INPUT_RAW_DATA_BASE, sample)
        summary_csv_path = os.path.join(sample_output_dir, f'{sample}_processing_summary.csv')
        if os.path.exists(summary_csv_path):
            sample_summary = pd.read_csv(summary_csv_path)
            all_summaries.append(sample_summary)

# %%
# Create combined summary report
if all_summaries:
    print("\nCreating combined summary report...")
    combined_summary = pd.concat(all_summaries, ignore_index=True)
    
    # Create a general output directory for combined reports
    general_output_dir = os.path.join(WORKING_DIR, "individual_sample_results")
    os.makedirs(general_output_dir, exist_ok=True)
    
    combined_summary_path = os.path.join(general_output_dir, "all_samples_processing_summary.csv")
    combined_summary.to_csv(combined_summary_path, index=False)
    
    print(f"Combined summary saved to {combined_summary_path}")
    print("\n--- COMBINED PROCESSING SUMMARY ---")
    print(combined_summary.to_string(index=False))

print(f"\n🎉 Individual sample processing complete!")
print(f"📁 Processed {len(processed_samples)} out of {len(SAMPLES)} samples successfully")
print(f"💾 Each sample's .h5ad file saved in its respective directory under {INPUT_RAW_DATA_BASE}")

# %%
# Display final summary
print("\n" + "="*80)
print("FINAL PROCESSING SUMMARY")
print("="*80)

for sample_name, adata in processed_samples.items():
    print(f"{sample_name}:")
    print(f"  • Cells: {adata.n_obs:,}")
    print(f"  • Genes: {adata.n_vars:,}")
    print(f"  • Clusters: {len(adata.obs[f'leiden_{LEIDEN_RESOLUTION}'].unique())}")
    print(f"  • Output: {os.path.join(INPUT_RAW_DATA_BASE, sample_name, f'{sample_name}_processed.h5ad')}")
    print()

print("All individual sample processing completed successfully! 🚀")