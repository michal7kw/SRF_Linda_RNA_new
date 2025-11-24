# %% [markdown]
# # 0. Merge Raw scRNA-seq Data and Process
# 
# This script loads raw count data (e.g., from Cell Ranger) for multiple samples,
# merges them into a single AnnData object, adds sample information,
# and then performs standard preprocessing, normalization, dimensionality reduction,
# and clustering on the combined dataset.

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

# Define Output Directory
OUTPUT_DIR = os.path.join(WORKING_DIR, "results_from_raw")
os.makedirs(OUTPUT_DIR, exist_ok=True)
output_file = os.path.join(OUTPUT_DIR, 'merged_raw_processed.h5ad')

# Define Samples to Merge
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
# ## Load and Merge Raw Data

# %%
adata_list = []
print("Loading raw data for each sample...")

for sample in tqdm(SAMPLES, desc="Loading Samples"):
    # Path to the directory containing matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz
    sample_dir_path = os.path.join(INPUT_RAW_DATA_BASE, sample, "outs", "filtered_feature_bc_matrix")
    print(f"  Loading {sample} from directory: {sample_dir_path}")
    if not os.path.isdir(sample_dir_path):
        print(f"  Warning: Directory not found for sample {sample} at {sample_dir_path}. Skipping.")
        continue
    try:
        # Use read_10x_mtx for Matrix Market format
        adata_sample = sc.read_10x_mtx(sample_dir_path, cache=True) # cache=True can speed up subsequent loads
        adata_sample.var_names_make_unique()  # Ensure unique gene names

        # Add sample information
        adata_sample.obs['sample'] = sample
        # Add condition information
        if "Ctrl" in sample:
            adata_sample.obs['condition'] = 'Control'
        elif "Mut" in sample:
            adata_sample.obs['condition'] = 'Mutant'
        else:
            adata_sample.obs['condition'] = 'Unknown'

        # Add genotype information
        if "Emx1" in sample:
            adata_sample.obs['genotype'] = 'Emx1'
        elif "Nestin" in sample:
            adata_sample.obs['genotype'] = 'Nestin'
        else:
            adata_sample.obs['genotype'] = 'Unknown'

        adata_list.append(adata_sample)
        print(f"  Loaded {adata_sample.n_obs} cells and {adata_sample.n_vars} genes for {sample}")
    except Exception as e:
        print(f"  Error loading {sample}: {e}")

# %%
if not adata_list:
    print("Error: No datasets were loaded. Exiting.")
    sys.exit(1)

print("\nConcatenating datasets...")
# Ensure all AnnData objects have the same variables (genes) before concatenating
# Using 'outer' join keeps all genes, filling missing values with 0
# Using 'inner' join keeps only common genes
adata_merged = ad.concat(adata_list, join='outer', label='sample_batch', index_unique='-') # index_unique adds sample suffix to barcodes
adata_merged.obs_names_make_unique() # Ensure unique cell barcodes after merging


print(f"Total cells in merged dataset: {adata_merged.n_obs}")
print(f"Total genes in merged dataset: {adata_merged.n_vars}")

# Make a copy of raw counts
adata_merged.raw = adata_merged.copy()

# %% [markdown]
# ## Quality Control (QC) and Filtering

# %%

# Set aesthetic styling
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

print("\nPerforming QC and Filtering...")

# Calculate mitochondrial gene percentage
adata_merged.var['mt'] = adata_merged.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata_merged, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# === STORE DATA BEFORE FILTERING ===
adata_before = adata_merged.copy()
print("Data stored for before/after comparison")

# --- Filtering ---
print(f"Initial cell count: {adata_merged.n_obs}")
print(f"Initial gene count: {adata_merged.n_vars}")

# Calculate the 95th percentile for max_genes
max_genes_percentile = np.percentile(adata_merged.obs['n_genes_by_counts'], 95)
print(f"Calculated 95th percentile for max genes: {max_genes_percentile:.0f}")

# --- Apply Filters ---
# 1. Filter cells based on min number of genes
n_obs_before = adata_merged.n_obs
adata_merged = adata_merged[np.array(adata_merged.obs.n_genes_by_counts) > MIN_GENES_PER_CELL]
print(f"Cells after min_genes filter (> {MIN_GENES_PER_CELL}): {adata_merged.n_obs} (removed {n_obs_before - adata_merged.n_obs})")

# 2. Filter cells based on max number of genes (percentile)
n_obs_before = adata_merged.n_obs
adata_merged = adata_merged[np.array(adata_merged.obs.n_genes_by_counts) < int(max_genes_percentile)]
print(f"Cells after max_genes filter (< {max_genes_percentile:.0f}): {adata_merged.n_obs} (removed {n_obs_before - adata_merged.n_obs})")

# 3. Filter cells based on mitochondrial content
n_obs_before = adata_merged.n_obs
adata_merged = adata_merged[np.array(adata_merged.obs.pct_counts_mt) < MAX_MITO_PERCENT]
print(f"Cells after mito filter (< {MAX_MITO_PERCENT}%): {adata_merged.n_obs} (removed {n_obs_before - adata_merged.n_obs})")

# 4. Filter genes based on minimum number of cells expressing them
n_vars_before = adata_merged.n_vars
sc.pp.filter_genes(adata_merged, min_cells=MIN_CELLS_PER_GENE)
print(f"Genes after min_cells filter ({MIN_CELLS_PER_GENE}): {adata_merged.n_vars} (removed {n_vars_before - adata_merged.n_vars})")

print("QC and Filtering complete.")

# === CREATE AESTHETIC VISUALIZATIONS ===
print("\nCreating quality control visualizations...")

def plot_qc_metrics_comparison(adata_before, adata_after, figsize=(16, 12)):
    """Create aesthetic before/after QC comparison plots"""
    
    # Colors
    before_color = '#E74C3C'  # Red
    after_color = '#2ECC71'   # Green
    threshold_color = '#F39C12'  # Orange
    
    fig, axes = plt.subplots(3, 4, figsize=figsize)
    fig.suptitle('Single-cell RNA-seq Quality Control: Before vs After Filtering', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # Calculate thresholds
    max_genes = np.percentile(adata_before.obs['n_genes_by_counts'], 95)
    
    # Row 1: Histograms
    # Genes per cell
    axes[0,0].hist(adata_before.obs['n_genes_by_counts'], bins=50, alpha=0.7, 
                   color=before_color, label='Before', density=True)
    axes[0,0].hist(adata_after.obs['n_genes_by_counts'], bins=50, alpha=0.7, 
                   color=after_color, label='After', density=True)
    axes[0,0].axvline(MIN_GENES_PER_CELL, color=threshold_color, linestyle='--', linewidth=2)
    axes[0,0].axvline(max_genes, color=threshold_color, linestyle='--', linewidth=2)
    axes[0,0].set_title('Genes per Cell')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)
    
    # Total UMI counts
    axes[0,1].hist(adata_before.obs['total_counts'], bins=50, alpha=0.7, 
                   color=before_color, label='Before', density=True)
    axes[0,1].hist(adata_after.obs['total_counts'], bins=50, alpha=0.7, 
                   color=after_color, label='After', density=True)
    axes[0,1].set_title('Total UMI Counts')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    
    # Mitochondrial percentage
    axes[0,2].hist(adata_before.obs['pct_counts_mt'], bins=50, alpha=0.7, 
                   color=before_color, label='Before', density=True)
    axes[0,2].hist(adata_after.obs['pct_counts_mt'], bins=50, alpha=0.7, 
                   color=after_color, label='After', density=True)
    axes[0,2].axvline(MAX_MITO_PERCENT, color=threshold_color, linestyle='--', linewidth=2)
    axes[0,2].set_title('Mitochondrial %')
    axes[0,2].legend()
    axes[0,2].grid(True, alpha=0.3)
    
    # Summary stats
    axes[0,3].axis('off')
    summary_text = f"""
    FILTERING SUMMARY
    
    Cells: {adata_before.n_obs:,} → {adata_after.n_obs:,}
    Removed: {adata_before.n_obs - adata_after.n_obs:,} ({((adata_before.n_obs - adata_after.n_obs)/adata_before.n_obs*100):.1f}%)
    
    Genes: {adata_before.n_vars:,} → {adata_after.n_vars:,}
    Removed: {adata_before.n_vars - adata_after.n_vars:,} ({((adata_before.n_vars - adata_after.n_vars)/adata_before.n_vars*100):.1f}%)
    
    THRESHOLDS USED:
    • Min genes/cell: {MIN_GENES_PER_CELL}
    • Max genes/cell: {max_genes:.0f}
    • Max mito %: {MAX_MITO_PERCENT}
    • Min cells/gene: {MIN_CELLS_PER_GENE}
    """
    axes[0,3].text(0.1, 0.9, summary_text, transform=axes[0,3].transAxes, 
                   fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.8))
    
    # Row 2: Box plots
    # Prepare data for box plots
    df_before = pd.DataFrame({
        'n_genes': adata_before.obs['n_genes_by_counts'],
        'total_counts': adata_before.obs['total_counts'],
        'pct_mt': adata_before.obs['pct_counts_mt'],
        'stage': 'Before'
    })
    df_after = pd.DataFrame({
        'n_genes': adata_after.obs['n_genes_by_counts'],
        'total_counts': adata_after.obs['total_counts'],
        'pct_mt': adata_after.obs['pct_counts_mt'],
        'stage': 'After'
    })
    df_combined = pd.concat([df_before, df_after])
    
    # Box plots
    sns.boxplot(data=df_combined, x='stage', y='n_genes', ax=axes[1,0],
                hue='stage', palette=[before_color, after_color], legend=False)
    axes[1,0].set_title('Genes per Cell Distribution')
    axes[1,0].grid(True, alpha=0.3)
    
    sns.boxplot(data=df_combined, x='stage', y='total_counts', ax=axes[1,1],
                hue='stage', palette=[before_color, after_color], legend=False)
    axes[1,1].set_title('Total UMI Distribution')
    axes[1,1].grid(True, alpha=0.3)
    
    sns.boxplot(data=df_combined, x='stage', y='pct_mt', ax=axes[1,2],
                hue='stage', palette=[before_color, after_color], legend=False)
    axes[1,2].set_title('Mitochondrial % Distribution')
    axes[1,2].grid(True, alpha=0.3)
    
    # Pie chart showing filtering impact
    removed_cells = adata_before.n_obs - adata_after.n_obs
    retained_cells = adata_after.n_obs
    sizes = [retained_cells, removed_cells]
    labels = [f'Retained\n{retained_cells:,}', f'Removed\n{removed_cells:,}']
    colors = [after_color, '#FFB3B3']
    
    axes[1,3].pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', 
                  startangle=90)
    axes[1,3].set_title('Cell Filtering Impact')
    
    # Row 3: Scatter plots
    # Before filtering scatter
    scatter1 = axes[2,0].scatter(adata_before.obs['total_counts'], 
                                adata_before.obs['n_genes_by_counts'], 
                                c=adata_before.obs['pct_counts_mt'], 
                                cmap='Reds', alpha=0.6, s=8)
    axes[2,0].set_xlabel('Total UMI Counts')
    axes[2,0].set_ylabel('Genes per Cell')
    axes[2,0].set_title('Before: Genes vs UMI (colored by mito %)')
    plt.colorbar(scatter1, ax=axes[2,0], label='Mito %')
    
    # After filtering scatter
    scatter2 = axes[2,1].scatter(adata_after.obs['total_counts'], 
                                adata_after.obs['n_genes_by_counts'], 
                                c=adata_after.obs['pct_counts_mt'], 
                                cmap='Greens', alpha=0.6, s=8)
    axes[2,1].set_xlabel('Total UMI Counts')
    axes[2,1].set_ylabel('Genes per Cell')
    axes[2,1].set_title('After: Genes vs UMI (colored by mito %)')
    plt.colorbar(scatter2, ax=axes[2,1], label='Mito %')
    
    # Mitochondrial vs UMI before
    axes[2,2].scatter(adata_before.obs['total_counts'], 
                      adata_before.obs['pct_counts_mt'], 
                      color=before_color, alpha=0.6, s=8)
    axes[2,2].axhline(MAX_MITO_PERCENT, color=threshold_color, linestyle='--', linewidth=2)
    axes[2,2].set_xlabel('Total UMI Counts')
    axes[2,2].set_ylabel('Mitochondrial %')
    axes[2,2].set_title('Before: UMI vs Mito %')
    
    # Mitochondrial vs UMI after
    axes[2,3].scatter(adata_after.obs['total_counts'], 
                      adata_after.obs['pct_counts_mt'], 
                      color=after_color, alpha=0.6, s=8)
    axes[2,3].set_xlabel('Total UMI Counts')
    axes[2,3].set_ylabel('Mitochondrial %')
    axes[2,3].set_title('After: UMI vs Mito %')
    
    plt.tight_layout()
    return fig

# Create the comprehensive visualization
fig = plot_qc_metrics_comparison(adata_before, adata_merged)

# Save the figure
plt.savefig('qc_before_after_comparison.png', dpi=300, bbox_inches='tight')
plt.savefig('qc_before_after_comparison.pdf', bbox_inches='tight')

print("✅ Quality control visualizations saved as 'qc_before_after_comparison.png/pdf'")
plt.show()


print("\n🎨 Visualization complete! Generated beautiful before/after QC plots.")
print(f"📊 Original dataset: {adata_before.n_obs:,} cells, {adata_before.n_vars:,} genes")
print(f"📊 Filtered dataset: {adata_merged.n_obs:,} cells, {adata_merged.n_vars:,} genes")
print(f"🗑️  Removed: {adata_before.n_obs - adata_merged.n_obs:,} cells ({((adata_before.n_obs - adata_merged.n_obs)/adata_before.n_obs*100):.1f}%) and {adata_before.n_vars - adata_merged.n_vars:,} genes ({((adata_before.n_vars - adata_merged.n_vars)/adata_before.n_vars*100):.1f}%)")


# %%
# Check the data type
print("Data type:", adata_merged.X.dtype) # type: ignore

# Look at some actual values
print("Sample values:")
print(adata_merged.X[:5, :5].toarray() if hasattr(adata_merged.X, 'toarray') else adata_merged.X[:5, :5]) # type: ignore

# Check if all values are integers
if hasattr(adata_merged.X, 'toarray'):
    sample_data = adata_merged.X[:1000, :100].toarray() # type: ignore
else:
    sample_data = adata_merged.X[:1000, :100] # type: ignore

print("Are all sampled values integers?", np.all(sample_data == sample_data.astype(int)))
print("Min value:", sample_data.min())
print("Max value:", sample_data.max())

# %%
# Identify Highly Variable Genes (HVGs)
print(f"Identifying top {N_TOP_HVG} highly variable genes...")
sc.pp.highly_variable_genes(adata_merged, n_top_genes=N_TOP_HVG, flavor='seurat_v3')
sc.pl.highly_variable_genes(adata_merged, show=True)


# %% [markdown]
# ## Normalization and Scaling

# %%
print("\nNormalizing and Scaling data...")

# Normalize counts per cell (library size normalization)
sc.pp.normalize_total(adata_merged, target_sum=1e4)
sc.pp.log1p(adata_merged)

# %%
# Subset data to HVGs (optional, but often done before scaling and PCA)
# adata_merged = adata_merged[:, adata_merged.var.highly_variable]
# print(f"Subsetting to {adata_merged.n_vars} HVGs.")

# %%
# Scale data to unit variance and zero mean (regress out total counts and mito % if desired)
print("Scaling data...")
# sc.pp.scale(adata_merged, max_value=10, vars_to_regress=['total_counts', 'pct_counts_mt']) # Example regression
sc.pp.scale(adata_merged, max_value=10) # Simpler scaling

print("Normalization and Scaling complete.")

# %%
print(adata_merged.layers)

# %% [markdown]
# ## Dimensionality Reduction

# %%
print("\nPerforming Dimensionality Reduction (PCA and UMAP)...")

# Use HVGs for PCA if subsetting was done, otherwise uses all genes if not subsetted before scaling
print(f"Running PCA using {adata_merged.n_vars} genes...")
sc.tl.pca(adata_merged, n_comps=N_PCS, use_highly_variable=True, svd_solver='arpack') # Ensure HVGs are used if identified

sc.pl.pca_variance_ratio(adata_merged, log=True, n_pcs=N_PCS, show=True)

print(f"Computing neighborhood graph using {N_PCS} PCs...")
sc.pp.neighbors(adata_merged, n_neighbors=15, n_pcs=N_PCS)

print("Running UMAP...")
sc.tl.umap(adata_merged)

print("Dimensionality Reduction complete.")

# %% [markdown]
# ## Clustering

# %%
adata_merged.write(Path('adata_merged_intermediate.h5ad'))

# %%
adata_merged = ad.read_h5ad('adata_merged_intermediate.h5ad')

# %%
print("\nPerforming Clustering (Leiden)...")

# Leiden clustering
print(f"Running Leiden clustering with resolution {LEIDEN_RESOLUTION}...")
sc.tl.leiden(adata_merged, resolution=LEIDEN_RESOLUTION, key_added=f'leiden_{LEIDEN_RESOLUTION}')

print("Clustering complete.")

# %% [markdown]
# ## Visualization and Saving

# %%
print("\nGenerating visualizations...")

# Create directory for plots if it doesn't exist
plot_dir = os.path.join(OUTPUT_DIR, "plots")
os.makedirs(plot_dir, exist_ok=True)
sc.settings.figdir = plot_dir

# Plot UMAP colored by sample, condition, genotype, and Leiden clusters
umap_plots = {
    'sample': 'sample',
    'condition': 'condition',
    'genotype': 'genotype',
    f'leiden_{LEIDEN_RESOLUTION}': f'leiden_{LEIDEN_RESOLUTION}'
}

for key, column in umap_plots.items():
    if column in adata_merged.obs:
        print(f"Plotting UMAP colored by {key}...")
        sc.pl.umap(adata_merged, color=column, legend_loc='right margin',
                   save=f"_umap_{key}.png", show=True,
                   title=f'UMAP colored by {key.replace("_", " ").title()}')
    else:
        print(f"Warning: Column '{column}' not found in adata_merged.obs. Skipping UMAP plot for {key}.")

# %%
# Save the processed merged datasetoutput_file = os.path.join(OUTPUT_DIR, 'merged_raw_processed.h5ad')
print(f"\nSaving processed merged dataset to {output_file}")
try:
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    adata_merged.write(Path(output_file))
    print("Successfully saved the AnnData object.")
except Exception as e:
    print(f"Error saving AnnData object: {e}")

print("\nScript finished!")


# # %%
# adata_merged

# %% [markdown]
# ## Post-Filtering Sample Summary &amp; Transcriptomic Complexity

# %%
print("\nGenerating post-filtering sample summary and complexity plots...")

# --- Create a summary DataFrame ---
# This contains the number of cells and genes for each sample after filtering
summary_df = adata_merged.obs.groupby('sample').agg(
    n_cells=('sample', 'size'),
    median_genes_per_cell=('n_genes_by_counts', 'median'),
    median_umi_per_cell=('total_counts', 'median')
).reset_index()

print("\n--- Summary of Cells and Genes per Sample (Post-Filtering) ---")
print(summary_df)
summary_csv_path = os.path.join(OUTPUT_DIR, "post_filtering_sample_summary.csv")
summary_df.to_csv(summary_csv_path, index=False)
print(f"Saved post-filtering summary to {summary_csv_path}")


# --- Create Histograms for Transcriptomic Complexity ---
print("\nCreating transcriptomic complexity histograms for each sample...")

# Get the list of samples
samples = adata_merged.obs['sample'].cat.categories.tolist()
n_samples = len(samples)

# Create a figure with subplots
# We'll have two rows of plots: one for genes per cell, one for UMI counts per cell
fig, axes = plt.subplots(2, n_samples, figsize=(5 * n_samples, 8), sharex='row', sharey='row')
fig.suptitle('Transcriptomic Complexity per Sample (Post-Filtering)', fontsize=16, fontweight='bold', y=1.02)

for i, sample in enumerate(samples):
    sample_obs = adata_merged.obs[adata_merged.obs['sample'] == sample]
    
    # Plot n_genes_by_counts
    sns.histplot(sample_obs['n_genes_by_counts'], bins=50, ax=axes[0, i], kde=True)
    axes[0, i].set_title(f'{sample}\n({sample_obs.shape[0]} cells)')
    axes[0, i].set_xlabel('Genes per Cell')
    axes[0, i].grid(True, alpha=0.3)
    if i == 0:
        axes[0, i].set_ylabel('Density')
    else:
        axes[0, i].set_ylabel('')


    # Plot total_counts
    sns.histplot(sample_obs['total_counts'], bins=50, ax=axes[1, i], kde=True)
    axes[1, i].set_xlabel('UMI Counts per Cell')
    axes[1, i].grid(True, alpha=0.3)
    if i == 0:
        axes[1, i].set_ylabel('Density')
    else:
        axes[1, i].set_ylabel('')

plt.tight_layout(rect=(0, 0, 1, 0.98))

# Save the figure
complexity_plot_path = os.path.join(plot_dir, 'sample_complexity_histograms.png')
plt.savefig(complexity_plot_path, dpi=300, bbox_inches='tight')
print(f"✅ Sample complexity histograms saved to {complexity_plot_path}")
plt.show()

print("\nSummary generation complete!")