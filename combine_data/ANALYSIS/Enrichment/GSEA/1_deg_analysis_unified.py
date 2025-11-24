import scanpy as sc
import numpy as np
import pandas as pd
import os
import sys
import pathlib
import argparse
import json
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

# --- Configuration ---
PROJECT_DIR = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

# Set up directories
BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")
INPUT_DIR = BASE_RESULTS_DIR
ADATA_PATH = os.path.join(INPUT_DIR, "annotation_final.h5ad")

# --- DEG Specific Configuration ---
CLUSTER_OF_INTEREST = 'Mature GC' # Initial cluster, will be updated if merging
CLUSTER_COLUMN = 'cell_type_L2_new'

MIN_CELLS  = 50

# --- Main Execution Block ---

def plot_volcano(de_results, p_val_col='pvals_adj', log2fc_col='logfoldchanges', gene_col='names', p_val_threshold=0.05, log2fc_threshold=0.5, title="Volcano Plot", output_path=None):
    """
    Generate and save a volcano plot.
    """
    # Calculate -log10 p-value
    de_results['-log10_p_val'] = -np.log10(de_results[p_val_col])

    # Determine significance and regulation status
    de_results['significant'] = (de_results[p_val_col] < p_val_threshold) & (abs(de_results[log2fc_col]) > log2fc_threshold)
    de_results['regulation'] = 'Not significant'
    de_results.loc[(de_results['significant']) & (de_results[log2fc_col] > 0), 'regulation'] = 'Upregulated'
    de_results.loc[(de_results['significant']) & (de_results[log2fc_col] < 0), 'regulation'] = 'Downregulated'

    # Create plot
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=de_results,
        x=log2fc_col,
        y='-log10_p_val',
        hue='regulation',
        palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not significant': 'grey'},
        edgecolor=None,
        s=20
    )

    # Add threshold lines
    plt.axhline(y=-np.log10(p_val_threshold), color='r', linestyle='--', linewidth=0.8)
    plt.axvline(x=log2fc_threshold, color='r', linestyle='--', linewidth=0.8)
    plt.axvline(x=-log2fc_threshold, color='r', linestyle='--', linewidth=0.8)

    # Add labels and title
    plt.title(title, fontsize=16)
    plt.xlabel("Log2 Fold Change", fontsize=12)
    plt.ylabel("-log10(Adjusted p-value)", fontsize=12)
    plt.legend(title='Regulation')

    # Save or show plot
    if output_path:
        plt.savefig(output_path, bbox_inches='tight', dpi=300)
        print(f"Volcano plot saved to {output_path}")
        plt.close()
    else:
        plt.show()


def main():
    """
    Main function to run the Differential Expression Gene (DEG) analysis pipeline.
    """
    parser = argparse.ArgumentParser(description="Run DEG analysis for specific genotypes and cluster merging options.")
    parser.add_argument('--genotype', type=str, required=True, choices=['Emx1', 'Nestin', 'both'],
                        help="Genotype to filter by (e.g., 'Emx1', 'Nestin', 'both').")
    parser.add_argument('--merge_clusters', type=str, required=True, choices=['True', 'False'],
                        help="Whether to merge Immature GC and Mature GC clusters ('True' or 'False').")
    
    args = parser.parse_args()
    
    genotype_filter = args.genotype
    merge_gc_clusters = args.merge_clusters == 'True'

    # --- Load Data ---
    print(f"Loading AnnData from {ADATA_PATH}...")
    try:
        adata = sc.read_h5ad(ADATA_PATH)
        print("AnnData object loaded.")
    except Exception as e:
        print(f"Error loading AnnData file: {e}")
        sys.exit(1)
        
    # --- Apply Genotype Filter ---
    if genotype_filter != "both":
        print(f"Filtering data for '{genotype_filter}' genotype...")
        adata = adata[adata.obs["genotype"] == genotype_filter].copy()
    else:
        print("No genotype filter applied (using 'both' genotypes).")

    if adata.shape[0] == 0:
        print(f"No cells found for genotype '{genotype_filter}'. Exiting.")
        sys.exit(0)

    # Using a local variable for the cluster of interest to handle merging
    target_cluster = CLUSTER_OF_INTEREST

    # --- Optional: Merge Clusters ---
    if merge_gc_clusters:
        print("\n--- Merging 'Immature GC' and 'Mature GC' into 'Combined GC' ---")
        clusters_to_merge = ['Immature GC', 'Mature GC']
        new_cluster_name = 'Combined GC'
        if not pd.api.types.is_categorical_dtype(adata.obs[CLUSTER_COLUMN]):
            adata.obs[CLUSTER_COLUMN] = adata.obs[CLUSTER_COLUMN].astype('category')
        present_clusters = [c for c in clusters_to_merge if c in adata.obs[CLUSTER_COLUMN].cat.categories]
        if present_clusters:
            if new_cluster_name not in adata.obs[CLUSTER_COLUMN].cat.categories:
                adata.obs[CLUSTER_COLUMN] = adata.obs[CLUSTER_COLUMN].cat.add_categories([new_cluster_name])
            adata.obs.loc[adata.obs[CLUSTER_COLUMN].isin(present_clusters), CLUSTER_COLUMN] = new_cluster_name
            adata.obs[CLUSTER_COLUMN] = adata.obs[CLUSTER_COLUMN].cat.remove_unused_categories()
            target_cluster = new_cluster_name # Update target cluster
            print(f"Successfully merged clusters. New target cluster: '{target_cluster}'")
        else:
            print("Warning: None of the clusters to merge were found. Skipping merge.")
    else:
        print("Cluster merging is disabled.")

    # --- Prepare for DEG Analysis ---
    # Construct output directory based on arguments
    genotype_dir_name = genotype_filter.lower() if genotype_filter != 'both' else 'both_genotypes'
    
    DEG_OUTPUT_DIR_BASE_DYNAMIC = os.path.join(BASE_RESULTS_DIR, 'deg_analysis_between_conditions', genotype_dir_name)

    sanitized_cluster_name = target_cluster.replace(' ', '_').replace('/', '_')
    comparison_name = "Mutant_vs_Control"
    cluster_output_dir = pathlib.Path(DEG_OUTPUT_DIR_BASE_DYNAMIC) / f"{sanitized_cluster_name}_{comparison_name}"
    cluster_output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nDEG analysis results will be saved to: {cluster_output_dir}")
    
    # Filter for the cluster of interest
    print(f"\nFiltering data for cluster: '{target_cluster}'")
    adata_cluster = adata[adata.obs[CLUSTER_COLUMN] == target_cluster].copy()

    # --- Filter lowly expressed genes ---
    print(f"\nFiltering genes for cluster: '{target_cluster}'")
    n_genes_before = adata_cluster.n_vars
    sc.pp.filter_genes(adata_cluster, min_cells=MIN_CELLS)
    n_genes_after = adata_cluster.n_vars
    print(f"  - Genes before filtering: {n_genes_before}")
    print(f"  - Genes after filtering (min_cells={MIN_CELLS}): {n_genes_after}")
    print(f"  - Genes removed: {n_genes_before - n_genes_after}")

    # --- Get and Save Cluster Sizes ---
    print(f"\nCalculating group sizes for cluster: '{target_cluster}'")
    condition_counts = adata_cluster.obs['condition'].value_counts()
    mutant_size = condition_counts.get('Mutant', 0)
    control_size = condition_counts.get('Control', 0)
    print(f"  - Mutant group size: {mutant_size}")
    print(f"  - Control group size: {control_size}")

    cluster_sizes = {
        'cluster_name': target_cluster,
        'mutant_count': int(mutant_size),
        'control_count': int(control_size)
    }
    sizes_file_path = cluster_output_dir / "cluster_sizes.json"
    with open(sizes_file_path, 'w') as f:
        json.dump(cluster_sizes, f, indent=4)
    print(f"Saved cluster sizes to {sizes_file_path}")

    # --- Prepare data for DEG analysis from raw counts ---
    print("\nPreparing data for DEG analysis from raw counts...")
    use_layer = None
    if hasattr(adata_cluster, 'raw') and adata_cluster.raw is not None and adata_cluster.raw.X is not None:
        print("Found .raw attribute. Creating a layer for DEG analysis.")
        
        # Create a temporary AnnData with raw counts to perform normalization and log1p
        adata_for_dge = ad.AnnData(adata_cluster.raw.X.copy())
        adata_for_dge.obs_names = adata_cluster.obs_names
        adata_for_dge.var_names = adata_cluster.raw.var_names

        # Print example raw counts for verification
        print("\n--- Example Raw Counts (before normalization) ---")
        try:
            # If it's a sparse matrix, convert to dense for printing
            if hasattr(adata_for_dge.X, "toarray"):
                print(adata_for_dge.X[:5, :5].toarray())
            else:
                print(adata_for_dge.X[:5, :5])
        except Exception as e:
            print(f"Could not print example counts: {e}")
        print("--------------------------------------------------")


        print("Normalizing total counts (target_sum=1e4)...")
        sc.pp.normalize_total(adata_for_dge, target_sum=1e4)

        print("Applying log1p transformation...")
        sc.pp.log1p(adata_for_dge)

        # Ensure the gene order matches the main adata_cluster object
        adata_for_dge = adata_for_dge[:, adata_cluster.var_names].copy()

        print("Storing result in adata_cluster.layers['for_DEGs']...")
        adata_cluster.layers['for_DEGs'] = adata_for_dge.X.copy()
        use_layer = 'for_DEGs'
        print("'for_DEGs' layer created.")
    else:
        print("Warning: adata_cluster.raw.X not found. Using adata_cluster.X for DEG analysis.")

    # --- Differential Expression Analysis ---
    print(f"\nRunning differential expression for '{target_cluster}' between Mutant and Control conditions...")
    sc.tl.rank_genes_groups(
        adata_cluster,
        groupby='condition',
        groups=['Mutant'],
        reference='Control',
        method='t-test',
        use_raw=False,
        layer=use_layer
    )
    
    # --- Prepare and Save Ranked Gene List ---
    print("Preparing ranked gene list and saving DEG results...")
    de_results = sc.get.rank_genes_groups_df(adata_cluster, group='Mutant')
    de_results = de_results.copy()
    
    # Save the full DEG results
    deg_output_file_excel = cluster_output_dir / "DEG_Results_Full.xlsx"
    deg_output_file_csv = cluster_output_dir / "DEG_Results_Full.csv"
    de_results.to_excel(deg_output_file_excel, index=False)
    de_results.to_csv(deg_output_file_csv, index=False)
    print(f"Full DEG results saved to: {deg_output_file_excel} and .csv")

    # --- Generate and Save Volcano Plot ---
    print("\nGenerating volcano plot...")
    volcano_plot_path = cluster_output_dir / "Volcano_Plot.png"
    
    # Check if necessary columns exist
    if 'pvals_adj' in de_results.columns and 'logfoldchanges' in de_results.columns:
        plot_volcano(
            de_results=de_results.copy(), # Pass a copy to avoid SettingWithCopyWarning
            p_val_col='pvals_adj',
            log2fc_col='logfoldchanges',
            gene_col='names',
            p_val_threshold=0.05,
            log2fc_threshold=0.25,
            title=f'Volcano Plot for {target_cluster}\n(Mutant: {mutant_size} cells vs. Control: {control_size} cells)',
            output_path=volcano_plot_path
        )
    else:
        print("Warning: 'pvals_adj' or 'logfoldchanges' not found in results. Skipping volcano plot.")

    print("\n--- DEG analysis complete. ---")

if __name__ == '__main__':
    main()