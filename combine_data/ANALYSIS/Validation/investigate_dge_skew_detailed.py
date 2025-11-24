#%% [markdown]
# # Detailed Investigation of DGE Skew
#
# This notebook provides a more in-depth investigation into the observed DGE skew,
# building upon the initial QC checks. We will examine:
# 1. Top highly expressed genes in Control vs. Mutant.
# 2. Percentage of ribosomal protein gene expression.
# 3. Expression of hemoglobin genes (potential contamination).
# 4. Potential batch effects.
# 5. Overall expression distributions.

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
from scipy.sparse import issparse

# Set plotting style
sns.set_theme(style="ticks")
plt.rcParams['figure.figsize'] = (6, 5) # Default figure size

# Ignore common warnings for cleaner output
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=ad.ImplicitModificationWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)


#%%
# ## Configuration
# --- User Configuration ---
TARGET_CLUSTER = 'HM' # <<< IMPORTANT: Verify this is a valid cluster name in 'mapmycells_second_layer'
TARGET_GENOTYPE = 'Nestin'
GROUPING_KEY = 'mapmycells_second_layer' # Column containing the cluster names
CONDITION_KEY = 'condition' # Column containing 'Control'/'Mutant'
GENOTYPE_KEY = 'genotype' # Column containing 'Emx1'/'Nestin'
MITO_PERCENT_KEY = 'percent_mt' # <<< IMPORTANT: Verify this is the correct column name for mitochondrial percentage
# List of potential batch keys to check in adata.obs
POTENTIAL_BATCH_KEYS = ['sample_id', 'batch', 'sequencing_run', 'library_prep_batch'] # Add more as needed
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
INVESTIGATION_PLOT_DIR_BASE = os.path.join(INPUT_DIR, f"investigation_plots_{TARGET_GENOTYPE}_{TARGET_CLUSTER}")
INVESTIGATION_PLOT_DIR_DETAILED = os.path.join(INPUT_DIR, f"investigation_plots_detailed_{TARGET_GENOTYPE}_{TARGET_CLUSTER}")
os.makedirs(INVESTIGATION_PLOT_DIR_DETAILED, exist_ok=True)

print(f"Project Directory: {PROJECT_DIR}")
print(f"Working Directory: {WORKING_DIR}")
print(f"Input AnnData Path: {adata_path}")
print(f"Base Investigation Plot Directory: {INVESTIGATION_PLOT_DIR_BASE}") # For summary.txt reference
print(f"Detailed Investigation Plot Directory: {INVESTIGATION_PLOT_DIR_DETAILED}")
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
    sys.exit(f"Exiting: AnnData file not found at {adata_path}")

# Display available metadata columns and layers
if adata:
    print("\nAvailable adata.obs columns:")
    print(list(adata.obs.columns))
    print("\nAvailable adata.layers:")
    print(list(adata.layers.keys()))
    print(f"\nChecking for '{GROUPING_KEY}' column...")
    if GROUPING_KEY not in adata.obs.columns:
        sys.exit(f"Error: Grouping key '{GROUPING_KEY}' not found in adata.obs!")
    else:
        print(f"'{GROUPING_KEY}' found. Unique values:")
        print(list(adata.obs[GROUPING_KEY].unique()))
        if TARGET_CLUSTER not in adata.obs[GROUPING_KEY].unique():
             print(f"Warning: Target cluster '{TARGET_CLUSTER}' not found in '{GROUPING_KEY}' column!")

    if GENOTYPE_KEY not in adata.obs.columns:
        sys.exit(f"Error: Genotype key '{GENOTYPE_KEY}' not found in adata.obs!")
    if CONDITION_KEY not in adata.obs.columns:
        sys.exit(f"Error: Condition key '{CONDITION_KEY}' not found in adata.obs!")
    if MITO_PERCENT_KEY and MITO_PERCENT_KEY not in adata.obs.columns:
        print(f"Warning: Mito percent key '{MITO_PERCENT_KEY}' not found. Mitochondrial plots might be affected.")
        # MITO_PERCENT_KEY = None # Optionally disable if strictly needed elsewhere

    if adata.raw is None and 'counts' not in adata.layers:
        print("\nWarning: adata.raw is None and 'counts' layer not found. Raw count based metrics might be unavailable.")


#%% [markdown]
# ## Subset Data
# Focus on the specific cluster and genotype of interest.

#%%
if adata and GROUPING_KEY in adata.obs.columns and GENOTYPE_KEY in adata.obs.columns:
    print(f"Subsetting data for {GROUPING_KEY} == '{TARGET_CLUSTER}' and {GENOTYPE_KEY} == '{TARGET_GENOTYPE}'...")
    cluster_mask = adata.obs[GROUPING_KEY] == TARGET_CLUSTER
    genotype_mask = adata.obs[GENOTYPE_KEY] == TARGET_GENOTYPE
    subset_mask = cluster_mask & genotype_mask
    adata_sub = adata[subset_mask].copy()

    print(f"\nSubset created with {adata_sub.n_obs} cells.")
    if adata_sub.n_obs > 0:
        print("Cell counts by condition in subset:")
        print(adata_sub.obs[CONDITION_KEY].value_counts())
        if not pd.api.types.is_categorical_dtype(adata_sub.obs[CONDITION_KEY]):
             adata_sub.obs[CONDITION_KEY] = adata_sub.obs[CONDITION_KEY].astype('category')
        if 'Control' in adata_sub.obs[CONDITION_KEY].cat.categories and 'Mutant' in adata_sub.obs[CONDITION_KEY].cat.categories:
             adata_sub.obs[CONDITION_KEY] = adata_sub.obs[CONDITION_KEY].cat.reorder_categories(['Control', 'Mutant'], ordered=True)
    else:
        sys.exit("Error: Subset is empty. Check TARGET_CLUSTER and TARGET_GENOTYPE.")
else:
    sys.exit("Skipping analysis due to missing data or keys for subsetting.")

# Determine layer for DGE-like analysis (normalized, log-transformed)
DGE_LAYER = None
if 'for_DEGs' in adata_sub.layers:
    DGE_LAYER = 'for_DEGs'
    print(f"Using layer '{DGE_LAYER}' for DGE-related analyses.")
elif adata_sub.X is not None: # Fallback to .X if 'for_DEGs' is not present
    DGE_LAYER = None # sc.pl functions use .X by default if layer is None
    print("Layer 'for_DEGs' not found. Using adata_sub.X for DGE-related analyses.")
else:
    print("Warning: Neither 'for_DEGs' layer nor adata_sub.X found. Some analyses might fail.")

# Determine layer for raw counts (for percentages like ribosomal, mito)
RAW_COUNT_LAYER_DATA = None
if adata_sub.raw is not None:
    RAW_COUNT_LAYER_DATA = adata_sub.raw.X
    adata_sub.obs['total_counts_raw'] = np.ravel(RAW_COUNT_LAYER_DATA.sum(axis=1))
    print("Using adata_sub.raw.X for raw count based calculations ('total_counts_raw' added to obs).")
elif 'counts' in adata_sub.layers: # Check for a 'counts' layer
    RAW_COUNT_LAYER_DATA = adata_sub.layers['counts']
    adata_sub.obs['total_counts_raw'] = np.ravel(RAW_COUNT_LAYER_DATA.sum(axis=1))
    print("Using adata_sub.layers['counts'] for raw count based calculations ('total_counts_raw' added to obs).")
else:
    print("Warning: No raw count data (adata_sub.raw or 'counts' layer) found. Raw count based metrics will be skipped or may be inaccurate if derived from processed layers.")
    # If 'total_counts_raw' was expected from the previous script's summary.txt, ensure it's present or recalculated.
    if 'total_counts_raw' not in adata_sub.obs.columns:
        print("   'total_counts_raw' not found in obs. Downstream raw count metrics will be affected.")


#%% [markdown]
# ## 1. Top Highly Expressed Genes Analysis
# Compare the most abundant transcripts between Control and Mutant groups.

#%%
if adata_sub and DGE_LAYER is not None or adata_sub.X is not None:
    print("\n--- Analyzing Top Highly Expressed Genes ---")
    n_top_genes = 30

    # Calculate mean expression for each group using the DGE layer or .X
    mean_expr_dfs = []
    for condition_val in adata_sub.obs[CONDITION_KEY].cat.categories:
        adata_cond = adata_sub[adata_sub.obs[CONDITION_KEY] == condition_val]
        if DGE_LAYER:
            data_matrix = adata_cond.layers[DGE_LAYER]
        else:
            data_matrix = adata_cond.X

        if issparse(data_matrix):
            mean_expr = pd.Series(np.ravel(data_matrix.mean(axis=0)), index=adata_cond.var_names, name=f"mean_expr_{condition_val}")
        else:
            mean_expr = pd.Series(data_matrix.mean(axis=0), index=adata_cond.var_names, name=f"mean_expr_{condition_val}")
        mean_expr_dfs.append(mean_expr.sort_values(ascending=False).head(n_top_genes))

    if len(mean_expr_dfs) == 2:
        print(f"\nTop {n_top_genes} genes in Control ({TARGET_CLUSTER} - {TARGET_GENOTYPE}):")
        print(mean_expr_dfs[0])
        print(f"\nTop {n_top_genes} genes in Mutant ({TARGET_CLUSTER} - {TARGET_GENOTYPE}):")
        print(mean_expr_dfs[1])

        # Get union of top genes for heatmap
        union_top_genes = pd.concat([mean_expr_dfs[0], mean_expr_dfs[1]]).index.unique().tolist()
        print(f"\nPlotting heatmap for {len(union_top_genes)} unique top genes.")

        # Ensure genes are in adata_sub.var_names
        union_top_genes = [g for g in union_top_genes if g in adata_sub.var_names]

        if union_top_genes:
            plt.figure() # Ensure a new figure context for sc.pl
            sc.pl.heatmap(adata_sub, var_names=union_top_genes, groupby=CONDITION_KEY, layer=DGE_LAYER,
                          cmap='viridis', dendrogram=True, standard_scale='var', show=True)
            plt.suptitle(f"Top Expressed Genes Heatmap\n({TARGET_GENOTYPE} - {TARGET_CLUSTER})", y=1.05)
            output_path = os.path.join(INVESTIGATION_PLOT_DIR_DETAILED, f"heatmap_top_genes_{TARGET_GENOTYPE}_{TARGET_CLUSTER}.png")
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            print(f"Saved heatmap: {output_path}")
            plt.close()
        else:
            print("No common top genes found to plot heatmap.")
    else:
        print("Could not generate top gene lists for both conditions.")
else:
    print("Skipping top highly expressed genes analysis: DGE layer or .X not available or subset empty.")


#%% [markdown]
# ## 2. Ribosomal Protein Gene Percentage
# High or variable ribosomal content can sometimes indicate quality issues or specific biological states.

#%%
if adata_sub and RAW_COUNT_LAYER_DATA is not None:
    print("\n--- Analyzing Ribosomal Protein Gene Percentage ---")
    # Identify ribosomal genes (common prefixes for mouse/human)
    # Mitochondrial ribosomal proteins (Mrp) are often excluded unless specifically investigating mitochondria

    # Determine var_names corresponding to RAW_COUNT_LAYER_DATA
    if RAW_COUNT_LAYER_DATA is (adata_sub.raw.X if adata_sub.raw else None): # Check identity
        current_var_names = adata_sub.raw.var_names
        print(f"Deriving ribosomal mask from adata_sub.raw.var_names ({len(current_var_names)} genes)")
    elif 'counts' in adata_sub.layers and RAW_COUNT_LAYER_DATA is adata_sub.layers['counts']: # Check identity
        current_var_names = adata_sub.var_names # Assumes 'counts' layer uses adata_sub.var_names
        print(f"Deriving ribosomal mask from adata_sub.var_names ({len(current_var_names)} genes for 'counts' layer)")
    else:
        # Fallback or error, but RAW_COUNT_LAYER_DATA should be one of these if set
        print("Warning: Could not definitively determine var_names for RAW_COUNT_LAYER_DATA. Using adata_sub.var_names as a fallback.")
        current_var_names = adata_sub.var_names

    rp_gene_mask_for_raw_layer = current_var_names.str.contains(r'^(Rps|Rpl)[0-9]', case=False, regex=True)
    num_rp_genes = rp_gene_mask_for_raw_layer.sum()
    print(f"Found {num_rp_genes} ribosomal genes in the source for raw counts.")

    if num_rp_genes > 0:
        # Manual calculation
        # Ensure ribo_sum_per_cell is a dense array for division
        ribo_sum_per_cell = np.asarray(RAW_COUNT_LAYER_DATA[:, rp_gene_mask_for_raw_layer].sum(axis=1)).ravel()
        
        # Ensure 'total_counts_raw' is available and correctly calculated from RAW_COUNT_LAYER_DATA
        # Also ensure it's a dense array for division
        if 'total_counts_raw' not in adata_sub.obs or not np.allclose(np.asarray(adata_sub.obs['total_counts_raw'].values), np.asarray(RAW_COUNT_LAYER_DATA.sum(axis=1)).ravel()):
            print("Recalculating 'total_counts_raw' from RAW_COUNT_LAYER_DATA for consistency.")
            adata_sub.obs['total_counts_raw'] = np.asarray(RAW_COUNT_LAYER_DATA.sum(axis=1)).ravel()
        
        total_sum_per_cell = adata_sub.obs['total_counts_raw'].values # This is already a dense numpy array

        adata_sub.obs['pct_counts_ribosomal'] = np.divide(ribo_sum_per_cell, total_sum_per_cell,
                                                          out=np.zeros_like(ribo_sum_per_cell, dtype=float),
                                                          where=total_sum_per_cell != 0) * 100
        print("Calculated 'pct_counts_ribosomal' manually.")
        
        # Add 'ribosomal' to adata_sub.var for completeness if needed by other parts, though not used by this manual calc directly for indexing RAW_COUNT_LAYER_DATA
        # This definition is based on adata_sub.var_names, so it's for the subsetted view.
        adata_sub.var['ribosomal'] = adata_sub.var_names.str.contains(r'^(Rps|Rpl)[0-9]', case=False, regex=True)


        if 'pct_counts_ribosomal' in adata_sub.obs.columns:
            plt.figure(figsize=(5, 4))
            sns.violinplot(data=adata_sub.obs, x=CONDITION_KEY, y='pct_counts_ribosomal', inner='quartile', cut=0)
            sns.stripplot(data=adata_sub.obs, x=CONDITION_KEY, y='pct_counts_ribosomal', color='black', size=2, alpha=0.3)
            plt.title(f"Ribosomal Gene Percentage\n({TARGET_GENOTYPE} - {TARGET_CLUSTER})")
            plt.xlabel("Condition")
            plt.ylabel("% Ribosomal Counts")
            sns.despine()
            plt.tight_layout()
            output_path = os.path.join(INVESTIGATION_PLOT_DIR_DETAILED, f"qc_ribo_percent_{TARGET_GENOTYPE}_{TARGET_CLUSTER}.png")
            plt.savefig(output_path, dpi=150)
            print(f"Saved plot: {output_path}")
            plt.close()

            # Add to summary statistics if possible (re-run original summary or append)
            print("Ribosomal percentage stats:")
            print(adata_sub.obs.groupby(CONDITION_KEY)['pct_counts_ribosomal'].describe())
        else:
            print("Could not calculate 'pct_counts_ribosomal'.")
    else:
        print("No ribosomal genes found with the specified pattern.")
else:
    print("Skipping ribosomal gene percentage: Raw count data not available or subset empty.")


#%% [markdown]
# ## 3. Hemoglobin Gene Expression
# Check for hemoglobin genes, which can indicate blood contamination if not expected.

#%%
if adata_sub:
    print("\n--- Analyzing Hemoglobin Gene Expression ---")
    # Common hemoglobin gene prefixes for mouse
    hb_genes_patterns = [r'^Hba-', r'^Hbb-', r'^Hbg-', r'^Hbe-', r'^Hbz-']
    # Use np.full for boolean array initialization to avoid DeprecationWarning
    hb_gene_mask = np.full(adata_sub.n_vars, False, dtype=bool)
    for pattern in hb_genes_patterns:
        hb_gene_mask |= adata_sub.var_names.str.contains(pattern, case=False, regex=True)
    
    actual_hb_genes = adata_sub.var_names[hb_gene_mask].tolist()
    print(f"Found {len(actual_hb_genes)} potential hemoglobin genes: {actual_hb_genes}")

    if actual_hb_genes:
        plt.figure()
        sc.pl.violin(adata_sub, keys=actual_hb_genes, groupby=CONDITION_KEY, layer=DGE_LAYER,
                      stripplot=True, multi_panel=True, show=True, cut=0, use_raw=False) # use_raw=False as DGE_LAYER is specified or .X is used
        plt.suptitle(f"Hemoglobin Gene Expression\n({TARGET_GENOTYPE} - {TARGET_CLUSTER})", y=1.02)
        output_path = os.path.join(INVESTIGATION_PLOT_DIR_DETAILED, f"expr_hemoglobin_{TARGET_GENOTYPE}_{TARGET_CLUSTER}.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved plot: {output_path}")
        plt.close()
    else:
        print("No hemoglobin genes found with the specified patterns.")
else:
    print("Skipping hemoglobin gene analysis: Subset empty.")


#%% [markdown]
# ## 4. Potential Batch Effects
# If batch information is available, visualize its impact.

#%%
if adata_sub:
    print("\n--- Checking for Potential Batch Effects ---")
    found_batch_keys = [key for key in POTENTIAL_BATCH_KEYS if key in adata_sub.obs.columns]

    if not found_batch_keys:
        print("No pre-defined batch keys found in adata_sub.obs.")
    else:
        for batch_key in found_batch_keys:
            print(f"\nAnalyzing batch key: '{batch_key}'")
            if adata_sub.obs[batch_key].nunique() < 2:
                print(f"  Batch key '{batch_key}' has only one unique value. Skipping.")
                continue

            # Ensure batch key is categorical
            if not pd.api.types.is_categorical_dtype(adata_sub.obs[batch_key]):
                adata_sub.obs[batch_key] = adata_sub.obs[batch_key].astype('category')

            # 1. Cell counts per batch and condition
            plt.figure(figsize=(8, 6))
            count_df = adata_sub.obs.groupby([batch_key, CONDITION_KEY]).size().unstack(fill_value=0)
            count_df.plot(kind='bar', stacked=True, ax=plt.gca())
            plt.title(f"Cell Counts by '{batch_key}' and Condition\n({TARGET_GENOTYPE} - {TARGET_CLUSTER})")
            plt.ylabel("Number of Cells")
            plt.xticks(rotation=45, ha='right')
            sns.despine()
            plt.tight_layout()
            output_path = os.path.join(INVESTIGATION_PLOT_DIR_DETAILED, f"batch_{batch_key}_counts_{TARGET_GENOTYPE}_{TARGET_CLUSTER}.png")
            plt.savefig(output_path, dpi=150)
            print(f"Saved plot: {output_path}")
            plt.close()

            # 2. Key QC metrics per batch
            qc_to_plot_by_batch = []
            if 'total_counts_raw' in adata_sub.obs: qc_to_plot_by_batch.append('total_counts_raw')
            if 'n_genes_by_counts' in adata_sub.obs: qc_to_plot_by_batch.append('n_genes_by_counts') # from sc.pp.calculate_qc_metrics
            elif 'n_genes_detected' in adata_sub.obs: qc_to_plot_by_batch.append('n_genes_detected') # from original script
            
            if MITO_PERCENT_KEY and MITO_PERCENT_KEY in adata_sub.obs: qc_to_plot_by_batch.append(MITO_PERCENT_KEY)
            if 'pct_counts_ribosomal' in adata_sub.obs: qc_to_plot_by_batch.append('pct_counts_ribosomal')


            for qc_metric in qc_to_plot_by_batch:
                plt.figure(figsize=(max(6, adata_sub.obs[batch_key].nunique() * 1.5), 5)) # Adjust width for many batches
                sns.violinplot(data=adata_sub.obs, x=batch_key, y=qc_metric, hue=CONDITION_KEY, inner='quartile', cut=0, dodge=True)
                plt.title(f"{qc_metric} by '{batch_key}' and Condition\n({TARGET_GENOTYPE} - {TARGET_CLUSTER})")
                plt.xticks(rotation=45, ha='right')
                sns.despine()
                plt.tight_layout()
                output_path = os.path.join(INVESTIGATION_PLOT_DIR_DETAILED, f"batch_{batch_key}_{qc_metric}_{TARGET_GENOTYPE}_{TARGET_CLUSTER}.png")
                plt.savefig(output_path, dpi=150)
                print(f"Saved plot: {output_path}")
                plt.close()
else:
    print("Skipping batch effect check: Subset empty.")


#%% [markdown]
# ## 5. Overall Expression Distributions
# Plot density of average gene expressions and total counts per cell.

#%%
if adata_sub:
    print("\n--- Analyzing Overall Expression Distributions ---")

    # 1. Density of average gene expression (using DGE layer or .X)
    plt.figure(figsize=(7, 5))
    for condition_val in adata_sub.obs[CONDITION_KEY].cat.categories:
        adata_cond = adata_sub[adata_sub.obs[CONDITION_KEY] == condition_val]
        if DGE_LAYER:
            data_matrix = adata_cond.layers[DGE_LAYER]
        else:
            data_matrix = adata_cond.X
        
        # Ensure mean_gene_expr is dense for kdeplot
        if issparse(data_matrix):
            mean_gene_expr = np.asarray(data_matrix.mean(axis=0)).ravel()
        else:
            mean_gene_expr = data_matrix.mean(axis=0) # Already dense
        
        sns.kdeplot(mean_gene_expr, label=condition_val, fill=True, alpha=0.3)
    plt.title(f"Density of Mean Gene Expression (per gene)\n({TARGET_GENOTYPE} - {TARGET_CLUSTER}, Layer: {DGE_LAYER if DGE_LAYER else '.X'})")
    plt.xlabel("Mean Expression (log-scale if applicable)")
    plt.ylabel("Density")
    plt.legend()
    sns.despine()
    plt.tight_layout()
    output_path = os.path.join(INVESTIGATION_PLOT_DIR_DETAILED, f"dist_mean_gene_expr_{TARGET_GENOTYPE}_{TARGET_CLUSTER}.png")
    plt.savefig(output_path, dpi=150)
    print(f"Saved plot: {output_path}")
    plt.close()

    # 2. Density of total counts per cell (Raw)
    if 'total_counts_raw' in adata_sub.obs.columns:
        plt.figure(figsize=(7, 5))
        for condition_val in adata_sub.obs[CONDITION_KEY].cat.categories:
            sns.kdeplot(adata_sub.obs.loc[adata_sub.obs[CONDITION_KEY] == condition_val, 'total_counts_raw'],
                        label=condition_val, fill=True, alpha=0.3)
        plt.title(f"Density of Total Raw Counts per Cell\n({TARGET_GENOTYPE} - {TARGET_CLUSTER})")
        plt.xlabel("Total Raw Counts")
        plt.ylabel("Density")
        plt.legend()
        sns.despine()
        plt.tight_layout()
        output_path = os.path.join(INVESTIGATION_PLOT_DIR_DETAILED, f"dist_total_counts_raw_{TARGET_GENOTYPE}_{TARGET_CLUSTER}.png")
        plt.savefig(output_path, dpi=150)
        print(f"Saved plot: {output_path}")
        plt.close()
    else:
        print("Skipping density plot of total raw counts: 'total_counts_raw' not in obs.")

    # 3. Density of sum of expression from DGE layer (if different from raw counts)
    # This is often sum of log-transformed normalized counts.
    sum_expr_col_dge = None
    if DGE_LAYER:
        sum_expr_col_dge = f'sum_expr_{DGE_LAYER}'
        # Ensure dense array for obs assignment
        adata_sub.obs[sum_expr_col_dge] = np.asarray(adata_sub.layers[DGE_LAYER].sum(axis=1)).ravel()
    elif adata_sub.X is not None: # Using .X
        sum_expr_col_dge = 'sum_expr_X'
        adata_sub.obs[sum_expr_col_dge] = np.asarray(adata_sub.X.sum(axis=1)).ravel()

    if sum_expr_col_dge and sum_expr_col_dge in adata_sub.obs.columns:
        plt.figure(figsize=(7, 5))
        for condition_val in adata_sub.obs[CONDITION_KEY].cat.categories:
            sns.kdeplot(adata_sub.obs.loc[adata_sub.obs[CONDITION_KEY] == condition_val, sum_expr_col_dge],
                        label=condition_val, fill=True, alpha=0.3)
        plt.title(f"Density of Summed Expression per Cell (from {DGE_LAYER if DGE_LAYER else '.X'})\n({TARGET_GENOTYPE} - {TARGET_CLUSTER})")
        plt.xlabel(f"Summed Expression ({DGE_LAYER if DGE_LAYER else '.X'})")
        plt.ylabel("Density")
        plt.legend()
        sns.despine()
        plt.tight_layout()
        output_path = os.path.join(INVESTIGATION_PLOT_DIR_DETAILED, f"dist_sum_expr_dge_layer_{TARGET_GENOTYPE}_{TARGET_CLUSTER}.png")
        plt.savefig(output_path, dpi=150)
        print(f"Saved plot: {output_path}")
        plt.close()
    else:
        print(f"Skipping density plot of summed DGE layer expression: '{sum_expr_col_dge}' not calculated.")

else:
    print("Skipping overall expression distribution analysis: Subset empty.")


#%% [markdown]
# ---
# # End of Detailed Investigation Script
# Review the generated plots in the `investigation_plots_detailed_...` directory.
# These plots, along with the initial QC, should provide a clearer picture of potential technical issues.
#
# Key things to look for:
# - **Top Expressed Genes:** Are they vastly different, or are the same genes present but at lower levels in Mutant?
# - **Ribosomal Percentage:** Is it unusually high or different between conditions?
# - **Hemoglobin:** Any unexpected expression?
# - **Batch Effects:** Do any of your `POTENTIAL_BATCH_KEYS` correlate with condition and QC metrics?
# - **Expression Distributions:** Do the overall shapes of distributions for counts/expression differ significantly, beyond a simple shift?
#
# The large discrepancy in `total_counts_raw` and `n_genes_detected` from `summary.txt` is the primary concern.
# These detailed plots will help to characterize this difference further.
# ---
print("\nDetailed investigation script finished.")