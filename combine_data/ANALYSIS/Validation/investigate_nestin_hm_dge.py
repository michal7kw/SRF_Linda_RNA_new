#%% [markdown]
# # Investigation of DGE Skew in Nestin-HM Comparison
#
# This notebook investigates potential reasons for the observed downregulation
# skew in the Differential Gene Expression (DGE) results when comparing
# Nestin_Mut vs Nestin_Ctrl within a specific cell cluster (e.g., 'HM').
# We will examine QC metrics and expression distributions.

#%%
# ## Environment Setup
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import anndata as ad
import warnings

# Set plotting style
sns.set_theme(style="ticks")
plt.rcParams['figure.figsize'] = (6, 5) # Default figure size

# Ignore common warnings for cleaner output
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=ad.ImplicitModificationWarning)


#%%
# ## Configuration
# --- User Configuration ---
TARGET_CLUSTER = 'HM' # <<< IMPORTANT: Verify this is a valid cluster name in 'mapmycells_second_layer'
TARGET_GENOTYPE = 'Nestin'
GROUPING_KEY = 'mapmycells_second_layer' # Column containing the cluster names
CONDITION_KEY = 'condition' # Column containing 'Control'/'Mutant'
GENOTYPE_KEY = 'genotype' # Column containing 'Emx1'/'Nestin'
MITO_PERCENT_KEY = 'percent_mt' # <<< IMPORTANT: Verify this is the correct column name for mitochondrial percentage in adata.obs, or set to None if not available
HOUSEKEEPING_GENES = ['Actb', 'Gapdh', 'Rpl10'] # Example housekeeping genes (verify they exist in data)
# --- End User Configuration ---

# Project Directories (adjust if needed)
# PROJECT_DIR = "/home/michal/Github/SRF_Linda_RNA" # Linux path
PROJECT_DIR = "D:/Github/SRF_Linda_RNA" # Windows path
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

# Input AnnData file (output from script 16)
INPUT_DIR = f"{WORKING_DIR}/results_from_raw/final_annotation"
ADATA_FILENAME = 'merged_raw_final_annotated_simple_mapmycells_merged_GC.h5ad'
adata_path = os.path.join(INPUT_DIR, ADATA_FILENAME)

# Output directory for plots from this notebook
INVESTIGATION_PLOT_DIR = os.path.join(INPUT_DIR, f"investigation_plots_{TARGET_GENOTYPE}_{TARGET_CLUSTER}")
os.makedirs(INVESTIGATION_PLOT_DIR, exist_ok=True)
print(f"Project Directory: {PROJECT_DIR}")
print(f"Working Directory: {WORKING_DIR}")
print(f"Input AnnData Path: {adata_path}")
print(f"Investigation Plot Directory: {INVESTIGATION_PLOT_DIR}")
print(f"Target Cluster ({GROUPING_KEY}): {TARGET_CLUSTER}")
print(f"Target Genotype ({GENOTYPE_KEY}): {TARGET_GENOTYPE}")


#%%
# ## Load Data
print(f"Loading AnnData object from: {adata_path}")
if os.path.exists(adata_path):
    adata = sc.read_h5ad(adata_path)
    print("AnnData loaded successfully.")
    print(adata)
else:
    print(f"Error: AnnData file not found at {adata_path}")
    adata = None # Set to None if loading fails

# Display available metadata columns and layers
if adata:
    print("\nAvailable adata.obs columns:")
    print(list(adata.obs.columns))
    print("\nAvailable adata.layers:")
    print(list(adata.layers.keys()))
    print(f"\nChecking for '{GROUPING_KEY}' column...")
    if GROUPING_KEY not in adata.obs.columns:
        print(f"Error: Grouping key '{GROUPING_KEY}' not found in adata.obs!")
    else:
        print(f"'{GROUPING_KEY}' found. Unique values:")
        print(list(adata.obs[GROUPING_KEY].unique()))
        if TARGET_CLUSTER not in adata.obs[GROUPING_KEY].unique():
             print(f"Warning: Target cluster '{TARGET_CLUSTER}' not found in '{GROUPING_KEY}' column!")

    print(f"\nChecking for '{GENOTYPE_KEY}' column...")
    if GENOTYPE_KEY not in adata.obs.columns:
        print(f"Error: Genotype key '{GENOTYPE_KEY}' not found in adata.obs!")
    else:
        print(f"'{GENOTYPE_KEY}' found. Unique values:")
        print(list(adata.obs[GENOTYPE_KEY].unique()))

    print(f"\nChecking for '{CONDITION_KEY}' column...")
    if CONDITION_KEY not in adata.obs.columns:
        print(f"Error: Condition key '{CONDITION_KEY}' not found in adata.obs!")
    else:
        print(f"'{CONDITION_KEY}' found. Unique values:")
        print(list(adata.obs[CONDITION_KEY].unique()))

    if MITO_PERCENT_KEY and MITO_PERCENT_KEY not in adata.obs.columns:
        print(f"Warning: Mito percent key '{MITO_PERCENT_KEY}' not found in adata.obs. Mitochondrial plots will be skipped.")
        MITO_PERCENT_KEY = None # Disable mito plots if key not found

    # Check if raw data is available for pre-normalization counts
    if adata.raw is None:
        print("\nWarning: adata.raw not found. Cannot plot pre-normalization total counts directly from .raw.")
        # Consider adding logic here to check for an alternative raw counts layer if needed


#%% [markdown]
# ## Subset Data
# Focus on the specific cluster and genotype of interest.

#%%
if adata and GROUPING_KEY in adata.obs.columns and GENOTYPE_KEY in adata.obs.columns:
    print(f"Subsetting data for {GROUPING_KEY} == '{TARGET_CLUSTER}' and {GENOTYPE_KEY} == '{TARGET_GENOTYPE}'...")

    # Create boolean masks
    cluster_mask = adata.obs[GROUPING_KEY] == TARGET_CLUSTER
    genotype_mask = adata.obs[GENOTYPE_KEY] == TARGET_GENOTYPE

    # Combine masks
    subset_mask = cluster_mask & genotype_mask

    # Perform subsetting
    adata_sub = adata[subset_mask].copy()

    # Check counts in the subset
    print(f"\nSubset created with {adata_sub.n_obs} cells.")
    if adata_sub.n_obs > 0:
        print("Cell counts by condition in subset:")
        print(adata_sub.obs[CONDITION_KEY].value_counts())

        # Ensure condition is categorical for plotting order
        if not pd.api.types.is_categorical_dtype(adata_sub.obs[CONDITION_KEY]):
             adata_sub.obs[CONDITION_KEY] = adata_sub.obs[CONDITION_KEY].astype('category')
        # Set order if possible
        if 'Control' in adata_sub.obs[CONDITION_KEY].cat.categories and 'Mutant' in adata_sub.obs[CONDITION_KEY].cat.categories:
             adata_sub.obs[CONDITION_KEY] = adata_sub.obs[CONDITION_KEY].cat.reorder_categories(['Control', 'Mutant'], ordered=True)

    else:
        print("Warning: Subset is empty. Check TARGET_CLUSTER and TARGET_GENOTYPE.")
        adata_sub = None # Set to None if subset is empty
else:
    print("Skipping subsetting due to missing data or keys.")
    adata_sub = None


#%% [markdown]
# ## QC Metric Comparison (Violin Plots)
# Compare key QC metrics between Control and Mutant cells within the subset.

#%%
# --- Plotting Function ---
def plot_qc_violin(adata_subset, qc_key, title, ylabel, output_dir, filename_suffix):
    """Helper function to generate and save QC violin plots."""
    if adata_subset is None or qc_key not in adata_subset.obs.columns:
        print(f"Skipping plot: {title} (data or key '{qc_key}' missing)")
        return

    plt.figure(figsize=(5, 4))
    sns.violinplot(data=adata_subset.obs, x=CONDITION_KEY, y=qc_key, inner='quartile', cut=0)
    sns.stripplot(data=adata_subset.obs, x=CONDITION_KEY, y=qc_key, color='black', size=2, alpha=0.3)
    plt.title(f"{title}\n({TARGET_GENOTYPE} - {TARGET_CLUSTER})")
    plt.xlabel("Condition")
    plt.ylabel(ylabel)
    sns.despine()
    plt.tight_layout()
    output_path = os.path.join(output_dir, f"qc_{filename_suffix}_{TARGET_GENOTYPE}_{TARGET_CLUSTER}.png")
    plt.savefig(output_path, dpi=150)
    print(f"Saved plot: {output_path}")
    plt.close() # Close figure

# --- Calculate Metrics if Needed ---
if adata_sub:
    # 1. Total Counts (Pre-Normalization) - Requires raw counts
    if adata_sub.raw is not None:
        print("Calculating total counts from adata_sub.raw...")
        adata_sub.obs['total_counts_raw'] = np.ravel(adata_sub.raw.X.sum(axis=1))
        plot_qc_violin(adata_sub, 'total_counts_raw', 'Total UMI Counts (Raw)', 'Total Counts', INVESTIGATION_PLOT_DIR, 'total_counts_raw')
    else:
        print("Skipping raw total counts plot (adata.raw not found).")
        # Optional: Check for a raw counts layer if .raw is not used
        # if 'counts' in adata_sub.layers:
        #     adata_sub.obs['total_counts_layer'] = np.ravel(adata_sub.layers['counts'].sum(axis=1))
        #     plot_qc_violin(adata_sub, 'total_counts_layer', ...)

    # 2. Total Counts (Post-Normalization) - Assuming DGE used log1p(normalized) data
    #    We sum the *non-logarithmized* normalized counts if possible, or plot sums from the log layer.
    #    Let's check if 'for_DEGs' layer exists, otherwise use .X
    target_layer_for_dge_sum = None
    if 'for_DEGs' in adata_sub.layers:
        target_layer_for_dge_sum = 'for_DEGs'
    elif adata_sub.X is not None:
         target_layer_for_dge_sum = 'X' # Fallback to .X

    if target_layer_for_dge_sum:
        print(f"Calculating sum of expression from layer: '{target_layer_for_dge_sum}' (Log-Transformed)")
        # Note: Summing log-transformed values isn't the same as total counts, but shows relative magnitude.
        adata_sub.obs[f'sum_expr_{target_layer_for_dge_sum}'] = np.ravel(adata_sub.layers[target_layer_for_dge_sum].sum(axis=1)) if target_layer_for_dge_sum in adata_sub.layers else np.ravel(adata_sub.X.sum(axis=1))
        plot_qc_violin(adata_sub, f'sum_expr_{target_layer_for_dge_sum}', f'Sum of Expression ({target_layer_for_dge_sum})', 'Sum(LogExpr)', INVESTIGATION_PLOT_DIR, f'sum_expr_{target_layer_for_dge_sum}')
    else:
        print("Skipping post-normalization sum plot (no suitable layer found).")


    # 3. Number of Genes Detected
    #    Recalculate based on the layer used for DGE (or raw if preferred)
    if target_layer_for_dge_sum:
        print(f"Calculating genes detected from layer: '{target_layer_for_dge_sum}'")
        layer_data = adata_sub.layers[target_layer_for_dge_sum] if target_layer_for_dge_sum in adata_sub.layers else adata_sub.X
        adata_sub.obs['n_genes_detected'] = np.ravel((layer_data > 0).sum(axis=1))
        plot_qc_violin(adata_sub, 'n_genes_detected', 'Number of Genes Detected', 'Genes Detected', INVESTIGATION_PLOT_DIR, 'n_genes_detected')
    elif adata_sub.raw is not None: # Fallback to raw if layer not found
         print("Calculating genes detected from adata_sub.raw...")
         adata_sub.obs['n_genes_detected_raw'] = np.ravel((adata_sub.raw.X > 0).sum(axis=1))
         plot_qc_violin(adata_sub, 'n_genes_detected_raw', 'Number of Genes Detected (Raw)', 'Genes Detected', INVESTIGATION_PLOT_DIR, 'n_genes_detected_raw')
    else:
        print("Skipping genes detected plot (no suitable data found).")


    # 4. Mitochondrial Percentage
    if MITO_PERCENT_KEY:
        plot_qc_violin(adata_sub, MITO_PERCENT_KEY, 'Mitochondrial Gene Percentage', '% Mito', INVESTIGATION_PLOT_DIR, 'percent_mt')
    else:
        print(f"Skipping mitochondrial percentage plot (key '{MITO_PERCENT_KEY}' not found or not specified).")

else:
    print("Skipping QC plots because subset is empty or adata was not loaded.")


#%% [markdown]
# ## Calculate and Save Summary Statistics
# Calculate descriptive statistics for the QC metrics grouped by condition and save to a text file.

#%%
if adata_sub:
    print("Calculating summary statistics...")
    summary_stats = {}
    qc_metrics_to_summarize = []

    # Check which metrics were successfully calculated and added to obs
    if 'total_counts_raw' in adata_sub.obs.columns:
        qc_metrics_to_summarize.append('total_counts_raw')
    # Find the sum_expr column name dynamically
    sum_expr_col = next((col for col in adata_sub.obs.columns if col.startswith('sum_expr_')), None)
    if sum_expr_col:
        qc_metrics_to_summarize.append(sum_expr_col)
    # Find the n_genes_detected column name dynamically
    n_genes_col = next((col for col in adata_sub.obs.columns if col.startswith('n_genes_detected')), None)
    if n_genes_col:
        qc_metrics_to_summarize.append(n_genes_col)
    if MITO_PERCENT_KEY and MITO_PERCENT_KEY in adata_sub.obs.columns:
        qc_metrics_to_summarize.append(MITO_PERCENT_KEY)

    if qc_metrics_to_summarize:
        # Group by condition and calculate describe() stats
        grouped_stats = adata_sub.obs.groupby(CONDITION_KEY)[qc_metrics_to_summarize].describe()
        summary_stats['qc_metrics'] = grouped_stats
        print("\nSummary Statistics (Grouped by Condition):")
        print(grouped_stats)

        # --- Save to File ---
        summary_filename = os.path.join(INVESTIGATION_PLOT_DIR, 'summary.txt')
        print(f"\nSaving summary statistics to: {summary_filename}")
        try:
            with open(summary_filename, 'w') as f:
                f.write(f"Summary Statistics for Cluster: {TARGET_CLUSTER}, Genotype: {TARGET_GENOTYPE}\n")
                f.write("="*40 + "\n\n")
                f.write("Cell Counts:\n")
                f.write(adata_sub.obs[CONDITION_KEY].value_counts().to_string())
                f.write("\n\n" + "="*40 + "\n\n")
                f.write("QC Metric Summary (Grouped by Condition):\n")
                f.write(grouped_stats.to_string())
                f.write("\n\n" + "="*40 + "\n")
            print("Successfully saved summary statistics.")
        except Exception as e:
            print(f"Error saving summary statistics: {e}")

    else:
        print("No QC metrics found in adata_sub.obs to calculate statistics for.")

else:
    print("Skipping summary statistics calculation because subset is empty.")


#%% [markdown]
# ## Housekeeping Gene Expression
# Visualize expression of common housekeeping genes to see if they are also affected.

#%%
if adata_sub:
    print("Plotting housekeeping gene expression...")
    # Use the layer that was likely used for DGE analysis
    layer_to_plot = 'for_DEGs' if 'for_DEGs' in adata_sub.layers else None # Use None to default to .X if layer missing

    # Filter genes to those present in the subset
    genes_present = [g for g in HOUSEKEEPING_GENES if g in adata_sub.var_names]
    missing_genes = [g for g in HOUSEKEEPING_GENES if g not in adata_sub.var_names]

    if missing_genes:
        print(f"  Warning: Housekeeping genes not found in subset: {missing_genes}")

    if genes_present:
        # Explicitly set use_raw=False when layer is specified
        sc.pl.violin(adata_sub, keys=genes_present, groupby=CONDITION_KEY, layer=layer_to_plot, use_raw=False,
                     stripplot=True, multi_panel=True, show=False, cut=0)
        plt.suptitle(f"Housekeeping Gene Expression\n({TARGET_GENOTYPE} - {TARGET_CLUSTER})", y=1.02)
        output_path = os.path.join(INVESTIGATION_PLOT_DIR, f"expr_housekeeping_{TARGET_GENOTYPE}_{TARGET_CLUSTER}.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved plot: {output_path}")
        plt.close()
    else:
        print("  Skipping housekeeping plot: None of the specified genes found in the subset.")

else:
    print("Skipping housekeeping gene plots because subset is empty.")


#%% [markdown]
# ---
# # End of Investigation Script
# Review the generated plots and `summary.txt` in the `investigation_plots_...` directory.
# - Differences in **Total UMI Counts (Raw)** might indicate capture efficiency issues.
# - Differences in **Sum of Expression (Log-Transformed)** or **Number of Genes Detected** post-normalization could reflect normalization artifacts or real biological effects.
# - Differences in **Mitochondrial Gene Percentage** could indicate quality differences.
# - Consistent downregulation of **Housekeeping Genes** would strongly suggest a technical or normalization issue rather than specific biological regulation.
# ---