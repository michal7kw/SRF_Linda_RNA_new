# %%
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import os
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches # Needed for custom legend

transgene_name = "Rosa26_SBP1"
annotation_cols = ['cell_type', 'mapmycells_first_layer', 'mapmycells_second_layer']
sel_annotation = annotation_cols[1] # Example: 'mapmycells_first_layer'

# %%
# Define file paths
transgene_adata_path = "counts_trans/Emx1_Mut.h5ad"
annotated_merged_adata_path = "combine_data/results_from_raw/final_annotation/merged_raw_final_annotated_simple_mapmycells.h5ad"
output_adata_path = "counts_trans/Emx1_Mut_annotated_final_original_umap.h5ad" # Updated output name

print(f"Transgene AnnData path: {transgene_adata_path}")
print(f"Annotated merged AnnData path: {annotated_merged_adata_path}")
print(f"Output AnnData path: {output_adata_path}")

# %%
# Load AnnData objects
print("\nLoading AnnData objects...")
try:
    adata_transgene = sc.read_h5ad(transgene_adata_path)
    print(f"Loaded adata_transgene with shape: {adata_transgene.shape}")
    if transgene_name not in adata_transgene.var_names:
        print(f"ERROR: Transgene '{transgene_name}' not found in adata_transgene.var_names.")
        exit()
except FileNotFoundError:
    print(f"ERROR: File not found at {transgene_adata_path}")
    exit()
except Exception as e:
    print(f"ERROR loading adata_transgene: {e}")
    exit()

try:
    adata_annotated_full = sc.read_h5ad(annotated_merged_adata_path)
    print(f"Loaded adata_annotated_full with shape: {adata_annotated_full.shape}")
except FileNotFoundError:
    print(f"ERROR: File not found at {annotated_merged_adata_path}")
    exit()
except Exception as e:
    print(f"ERROR loading adata_annotated_full: {e}")
    exit()

# %%
# Filter for the target sample and prepare base AnnData
print("\nFiltering for target sample ('Emx1_Mut') and preparing base AnnData...")
if 'sample' not in adata_annotated_full.obs.columns:
    print("ERROR: 'sample' column not found in adata_annotated_full.obs")
    exit()

adata_main = adata_annotated_full[adata_annotated_full.obs['sample'] == 'Emx1_Mut'].copy()
print(f"Shape of adata_main after sample filtering: {adata_main.shape}")

if adata_main.n_obs == 0:
    print("ERROR: No cells found for sample 'Emx1_Mut'. Check the sample name.")
    exit()

# Verify necessary embeddings exist in adata_main
required_obsm = ['X_umap', 'X_pca']
for req_obsm in required_obsm:
    if req_obsm not in adata_main.obsm:
        print(f"ERROR: Required embedding '{req_obsm}' not found in adata_main.obsm.")
        print("Make sure the input annotated AnnData has been processed with PCA and UMAP.")
        exit()
print(f"Found required embeddings in adata_main.obsm: {list(adata_main.obsm.keys())}")


# Ensure 'original_barcode' column exists in adata_main.obs for mapping
if 'original_barcode' not in adata_main.obs.columns:
    print("Creating 'original_barcode' column in adata_main.obs")
    adata_main.obs['original_barcode'] = adata_main.obs_names.str.replace(r'-\d+$', '', regex=True)
else:
    print("'original_barcode' column already exists in adata_main.obs")

# Determine the source of raw counts for adata_main to be augmented
if adata_main.raw is not None:
    print("Using adata_main.raw.X as the base for original host raw counts.")
    original_host_raw_X = adata_main.raw.X.copy()
    original_host_raw_var = adata_main.raw.var.copy()
else:
    print("Using adata_main.X as the base for original host raw counts (adata_main.raw was None).")
    print("WARNING: adata_main.X is assumed to be raw counts for this step. If it's processed, raw transgene integration might be inconsistent.")
    original_host_raw_X = adata_main.X.copy()
    original_host_raw_var = adata_main.var.copy()

if transgene_name in original_host_raw_var.index:
    print(f"ERROR: Transgene name '{transgene_name}' already exists in the host dataset's var_names. This script expects it to be a new gene.")
    exit()

# %%
# Integrate Transgene Raw Counts
print("\nIntegrating transgene raw counts...")

tg_gene_loc = adata_transgene.var_names.get_loc(transgene_name)
transgene_expression_values_from_tg_adata = adata_transgene.X[:, tg_gene_loc]
if sp.issparse(transgene_expression_values_from_tg_adata):
    transgene_expression_values_from_tg_adata = transgene_expression_values_from_tg_adata.toarray()
transgene_expression_values_from_tg_adata = transgene_expression_values_from_tg_adata.flatten()

transgene_barcode_to_expression_map = pd.Series(
    transgene_expression_values_from_tg_adata,
    index=adata_transgene.obs_names
)

mapped_transgene_expression_for_main = adata_main.obs['original_barcode'].map(transgene_barcode_to_expression_map).fillna(0).values
mapped_transgene_expression_sparse = sp.csr_matrix(mapped_transgene_expression_for_main.reshape(-1, 1))

combined_raw_X = sp.hstack([original_host_raw_X, mapped_transgene_expression_sparse], format='csr')

new_var_row_data = {col: pd.NA for col in original_host_raw_var.columns}
if 'gene_ids' in new_var_row_data and 'gene_ids' not in original_host_raw_var.columns:
    new_var_row_data['gene_ids'] = transgene_name
if 'feature_types' in new_var_row_data:
    new_var_row_data['feature_types'] = 'Gene Expression'

new_var_entry_df = pd.DataFrame(
    [new_var_row_data],
    index=pd.Index([transgene_name], name=original_host_raw_var.index.name if original_host_raw_var.index.name else 'gene_name')
)
combined_raw_var = pd.concat([original_host_raw_var, new_var_entry_df])

# Create the AnnData object. Its .X will initially be the combined raw counts.
adata_to_process = sc.AnnData(X=combined_raw_X.copy(), obs=adata_main.obs.copy(), var=combined_raw_var.copy())

# Store this combined raw data in its own .raw attribute
adata_to_process.raw = adata_to_process.copy() # .raw.X contains combined raw counts
print(f"Shape of combined raw data (adata_to_process.raw.X): {adata_to_process.raw.X.shape}")
print(f"Transgene '{transgene_name}' added to var_names: {transgene_name in adata_to_process.var_names}")

# %%
# Normalize the combined dataset (host + transgene) for adata_to_process.X
print("\nNormalizing combined data (host + transgene) for .X layer...")
# adata_to_process.X currently holds raw combined counts from the previous step
sc.pp.normalize_total(adata_to_process, target_sum=1e4)
sc.pp.log1p(adata_to_process)
print(f"adata_to_process.X is now normalized and log-transformed combined data. Shape: {adata_to_process.X.shape}")

# %%
# Copy embeddings and neighbor graphs from adata_main
print("\nCopying existing UMAP, PCA, and neighbor information from adata_main...")
adata_to_process.obsm['X_umap'] = adata_main.obsm['X_umap'].copy()
adata_to_process.obsm['X_pca'] = adata_main.obsm['X_pca'].copy()
print(f"Copied X_umap (shape: {adata_to_process.obsm['X_umap'].shape}) and X_pca (shape: {adata_to_process.obsm['X_pca'].shape})")

# Copy other .obsm if needed, e.g., other dimensionality reductions
for key in adata_main.obsm:
    if key not in adata_to_process.obsm: # Avoid re-copying X_umap, X_pca
        adata_to_process.obsm[key] = adata_main.obsm[key].copy()
        print(f"Copied additional obsm key: {key}")

if adata_main.obsp:
    adata_to_process.obsp = adata_main.obsp.copy()
    print(f"Copied .obsp (neighbor graphs): {list(adata_to_process.obsp.keys())}")
else:
    print("adata_main.obsp is empty, nothing to copy for neighbor graphs.")
    
# Copy .uns if it contains relevant info like 'neighbors' settings or UMAP params, or color schemes
if adata_main.uns:
    adata_to_process.uns = adata_main.uns.copy() # Be cautious with this, copy selectively if needed
    print(f"Copied .uns. Keys include: {list(adata_to_process.uns.keys())[:5]}")


print(f"Shape of final adata_to_process: {adata_to_process.shape}")
print(f"UMAP coordinates available and copied: {'X_umap' in adata_to_process.obsm}")

# %%
# Annotations are already in adata_to_process.obs from adata_main
print("\nVerifying annotations and transgene presence...")
# print("Available annotation columns (subset):", [col for col in annotation_cols if col in adata_to_process.obs][:5])
if sel_annotation not in adata_to_process.obs:
    print(f"ERROR: Selected annotation column '{sel_annotation}' not found in adata_to_process.obs.")
    if annotation_cols[0] in adata_to_process.obs:
        sel_annotation = annotation_cols[0]
        print(f"Falling back to '{sel_annotation}' for plotting.")
    else:
        print("No suitable annotation column found for plotting. Exiting plot section.")
        sel_annotation = None
else:
    print(f"Using '{sel_annotation}' for plotting.")

# Add a 'transgene_positive' column based on raw counts for clarity in plots/analysis
if transgene_name in adata_to_process.raw.var_names:
    tg_expr_raw = adata_to_process.raw[:, transgene_name].X
    if sp.issparse(tg_expr_raw):
        tg_expr_raw = tg_expr_raw.toarray()
    adata_to_process.obs[f'{transgene_name}_raw_counts'] = tg_expr_raw.flatten()
    adata_to_process.obs[f'{transgene_name}_positive_in_raw'] = adata_to_process.obs[f'{transgene_name}_raw_counts'] > 0
    print(f"Added '{transgene_name}_raw_counts' and '{transgene_name}_positive_in_raw' to .obs")
else:
    print(f"Warning: Transgene '{transgene_name}' not found in adata_to_process.raw.var_names for positive flag calculation.")

# %%
# Plotting
if sel_annotation and transgene_name in adata_to_process.raw.var_names:
    print("\nPreparing data for plotting transgene expression (violin plot using raw counts)...")
    valid_cell_types_for_plot = sorted([
        ct for ct in adata_to_process.obs[sel_annotation].unique()
        if pd.notna(ct)
    ])
    print(f"Found valid cell types for plot using '{sel_annotation}': {valid_cell_types_for_plot}")

    if not valid_cell_types_for_plot:
        print("No valid cell types found for plotting.")
    else:
        gene_idx_plot = adata_to_process.raw.var_names.get_loc(transgene_name)
        plot_data_list = []

        for cell_type_plot in valid_cell_types_for_plot:
            cell_mask_series_plot = adata_to_process.obs[sel_annotation] == cell_type_plot
            if not cell_mask_series_plot.any():
                continue
            cell_mask_numpy_plot = cell_mask_series_plot.to_numpy()

            expression_values_for_gene_plot = adata_to_process.raw.X[cell_mask_numpy_plot, gene_idx_plot]
            if sp.issparse(expression_values_for_gene_plot):
                dense_expression_plot = expression_values_for_gene_plot.toarray().flatten()
            else:
                dense_expression_plot = np.asarray(expression_values_for_gene_plot).flatten()

            non_zero_expression_plot = dense_expression_plot[dense_expression_plot > 0]

            if non_zero_expression_plot.size > 0:
                temp_df_plot = pd.DataFrame({
                    'Cell Type': [cell_type_plot] * len(non_zero_expression_plot),
                    'Expression': non_zero_expression_plot
                 })
                plot_data_list.append(temp_df_plot)

        if plot_data_list:
            df_plot_final = pd.concat(plot_data_list, ignore_index=True)
            print("\nHead of plotting DataFrame (df_plot_final):")
            print(df_plot_final.head())

            plot_order_dynamic = sorted(df_plot_final['Cell Type'].unique())
            current_plot_order = plot_order_dynamic

            plt.figure(figsize=(max(8, len(current_plot_order) * 1.5), 6))
            ax = sns.violinplot(
                x='Cell Type',
                y='Expression',
                data=df_plot_final,
                order=current_plot_order,
                palette='viridis',
                hue='Cell Type',
                hue_order=current_plot_order,
                inner='quartile',
                cut=0,
                dodge=False
            )
            if not ax.get_legend() is None:
                 ax.get_legend().remove()

            plt.title(f'Distribution of Non-Zero Raw {transgene_name} Expression by {sel_annotation}', fontsize=16)
            plt.xlabel(sel_annotation, fontsize=14)
            plt.ylabel(f'Non-Zero Raw {transgene_name} Expression', fontsize=14)
            plt.xticks(rotation=45, ha='right', fontsize=11)
            plt.yticks(fontsize=11)
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            plt.tight_layout()
            print("Displaying violin plot...")
            plt.show()
        else:
            print(f"No non-zero expression data found for {transgene_name} across specified cell types for violin plot.")
else:
    print("\nSkipping violin plotting section: sel_annotation not available or transgene not found in raw data.")

# Example of plotting UMAP with transgene expression (using processed data from .X)
if 'X_umap' in adata_to_process.obsm and transgene_name in adata_to_process.var_names:
    print(f"\nPlotting UMAP colored by processed '{transgene_name}' expression...")
    plt.figure()
    sc.pl.umap(adata_to_process, color=transgene_name, show=False,
               title=f'{transgene_name} (Processed Expression) on Original UMAP')
    plt.show()
else:
    print(f"\nSkipping UMAP plot for '{transgene_name}': UMAP or transgene not available in expected locations.")


# %%
# Save the result
print("\nSaving the processed AnnData object...")
output_dir = os.path.dirname(output_adata_path)
if output_dir and not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Created output directory: {output_dir}")

try:
    # Ensure .raw is correctly set up if it was derived from an X that changed
    # In our case, adata_to_process.raw was set from an AnnData whose .X was combined_raw_X
    # and .var was combined_raw_var. This should be fine.
    adata_to_process.write_h5ad(output_adata_path)
    print(f"Successfully saved processed AnnData to: {output_adata_path}")
except Exception as e:
    print(f"ERROR saving AnnData object: {e}")

print("\nScript finished.")
