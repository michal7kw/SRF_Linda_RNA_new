# %%
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import pathlib
import gseapy as gp
import multiprocessing
from functools import partial
import traceback
import json
import argparse
import anndata as ad
import re

# --- Configuration ---
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

ORGANISM = 'Mouse'

# Set up directories
BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw", "gsea_analysis_between_conditions")
INPUT_DIR = os.path.join(WORKING_DIR, "results_from_raw")
ADATA_PATH = os.path.join(INPUT_DIR, "annotation_final.h5ad")

# --- GSEA Specific Configuration (defaults, will be overridden by args) ---
CLUSTER_OF_INTEREST = 'Mature GC' # Initial cluster, will be updated if merging
CLUSTER_COLUMN = 'cell_type_L2_new'
GSEA_GENE_SETS = ['KEGG_2019_Mouse', 'GO_Biological_Process_2023', 'Reactome_2022', 'MSigDB_Hallmark_2020']
N_TOP_TERMS_PLOT = 15
GSEA_PADJ_THRESHOLD = 0.25
MIN_CELLS = 10

# --- Worker Function for Parallel GSEA ---

def run_gsea_for_gene_set(gene_set, rank_df, organism, cluster_output_dir, adata_var_names, gsea_padj_threshold, n_top_terms_plot):
    """
    Performs GSEA analysis for a single gene set.
    This function is designed to be called by a multiprocessing Pool.
    """
    print(f"--- Starting GSEA for gene set: {gene_set} ---")
    
    gsea_lib_output_dir = cluster_output_dir / gene_set
    gsea_lib_output_dir.mkdir(exist_ok=True)

    try:
        pre_res = gp.prerank(rnk=rank_df,
                             gene_sets=gene_set,
                             organism=organism,
                             outdir=str(gsea_lib_output_dir),
                             no_plot=True,
                             verbose=True,
                             background=adata_var_names,
                             min_size=15,
                             max_size=500)

        # --- Compatibility fix for older gseapy versions ---
        if 'fdr' not in pre_res.res2d.columns and 'FDR q-val' in pre_res.res2d.columns:
            pre_res.res2d.rename(columns={'FDR q-val': 'fdr'}, inplace=True)
            for term in pre_res.results:
                if 'FDR q-val' in pre_res.results[term]:
                    pre_res.results[term]['fdr'] = pre_res.results[term]['FDR q-val']
        if 'nes' not in pre_res.res2d.columns and 'NES' in pre_res.res2d.columns:
            pre_res.res2d.rename(columns={'NES': 'nes'}, inplace=True)
            for term in pre_res.results:
                if 'NES' in pre_res.results[term]:
                    pre_res.results[term]['nes'] = pre_res.results[term]['NES']

        results_df = pre_res.res2d
        
        # --- Generate Plots for Significant Terms ---
        if 'fdr' in results_df.columns and 'nes' in results_df.columns:
            significant_results = results_df[results_df['fdr'] < gsea_padj_threshold]
            if not significant_results.empty:
                print(f"Found {len(significant_results)} significant terms in {gene_set}. Visualizing top {n_top_terms_plot}.")
                
                top_positive = significant_results[significant_results['nes'] > 0].sort_values('nes', ascending=False).head(n_top_terms_plot)
                top_negative = significant_results[significant_results['nes'] < 0].sort_values('nes', ascending=True).head(n_top_terms_plot)
                
                for term_df, direction in [(top_positive, 'positive'), (top_negative, 'negative')]:
                    if term_df.empty: continue
                    
                    term_name_col = 'Term' if 'Term' in term_df.columns else 'term'
                    for _, row in term_df.iterrows():
                        term_name = row[term_name_col]
                        try:
                            # Sanitize term name for filename compatibility
                            sanitized_term_name = re.sub(r'[^a-zA-Z0-9_.-]', '_', str(term_name))[:100]
                            plot_filename = gsea_lib_output_dir / f"{direction}_{sanitized_term_name}.png"
                            
                            rank_metric = pre_res.rank_metric if hasattr(pre_res, 'rank_metric') else pre_res.ranking
                            gp.gseaplot(rank_metric=rank_metric, term=term_name, **pre_res.results[term_name], ofname=str(plot_filename))
                            plt.close()
                        except Exception as plot_err:
                            print(f"  Could not generate plot for term '{term_name}' in {gene_set}. Error: {plot_err}")
        
        print(f"--- Finished GSEA for gene set: {gene_set} ---")
        return results_df

    except Exception as e:
        print(f"An error occurred during GSEA for gene set '{gene_set}': {e}")
        traceback.print_exc()
        return None


# --- Main Execution Block ---

def main():
    """
    Main function to run the GSEA analysis pipeline.
    """
    parser = argparse.ArgumentParser(description="Run GSEA analysis for specific genotypes and cluster merging options.")
    parser.add_argument('--genotype', type=str, required=True, choices=['Emx1', 'Nestin', 'both'],
                        help="Genotype to filter by (e.g., 'Emx1', 'Nestin', 'both').")
    parser.add_argument('--merge_clusters', type=str, required=True, choices=['True', 'False'],
                        help="Whether to merge Immature GC and Mature GC clusters ('True' or 'False'). This is ignored if --cluster_name is provided.")
    parser.add_argument('--cluster_name', type=str, help="Name of the cluster of interest (e.g., 'GABA'). Overrides merge_clusters.")
    parser.add_argument('--cluster_column', type=str, help="Column in adata.obs for the cluster (e.g., 'cell_type_L1').")
    
    args = parser.parse_args()
    
    genotype_filter = args.genotype
    merge_gc_clusters = args.merge_clusters == 'True'
    cluster_name_arg = args.cluster_name
    cluster_column_arg = args.cluster_column

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
    cluster_column = CLUSTER_COLUMN

    # --- Optional: Merge Clusters or Select Custom Cluster ---
    if cluster_name_arg and cluster_column_arg:
        print(f"\n--- Selecting custom cluster '{cluster_name_arg}' from column '{cluster_column_arg}' ---")
        target_cluster = cluster_name_arg
        cluster_column = cluster_column_arg
        print(f"Target cluster set to '{target_cluster}' from '{cluster_column}'. Merging logic is skipped.")
    elif merge_gc_clusters:
        print("\n--- Merging 'Immature GC' and 'Mature GC' into 'Combined GC' ---")
        clusters_to_merge = ['Immature GC', 'Mature GC']
        new_cluster_name = 'Combined GC'
        if not pd.api.types.is_categorical_dtype(adata.obs[cluster_column]):
            adata.obs[cluster_column] = adata.obs[cluster_column].astype('category')
        present_clusters = [c for c in clusters_to_merge if c in adata.obs[cluster_column].cat.categories]
        if present_clusters:
            if new_cluster_name not in adata.obs[cluster_column].cat.categories:
                adata.obs[cluster_column] = adata.obs[cluster_column].cat.add_categories([new_cluster_name])
            adata.obs.loc[adata.obs[cluster_column].isin(present_clusters), cluster_column] = new_cluster_name
            adata.obs[cluster_column] = adata.obs[cluster_column].cat.remove_unused_categories()
            target_cluster = new_cluster_name # Update target cluster
            print(f"Successfully merged clusters. New target cluster: '{target_cluster}'")
        else:
            print("Warning: None of the clusters to merge were found. Skipping merge.")
    else:
        print("Cluster merging is disabled.")

    # --- Prepare for GSEA ---
    # Construct output directory based on arguments
    genotype_dir_name = genotype_filter.lower() if genotype_filter != 'both' else 'both_genotypes'
    
    GSEA_OUTPUT_DIR_BASE_DYNAMIC = os.path.join(BASE_RESULTS_DIR, genotype_dir_name)

    sanitized_cluster_name = target_cluster.replace(' ', '_').replace('/', '_')
    comparison_name = "Mutant_vs_Control"
    cluster_output_dir = pathlib.Path(GSEA_OUTPUT_DIR_BASE_DYNAMIC) / f"{sanitized_cluster_name}_{comparison_name}"
    cluster_output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nGSEA output will be saved to: {cluster_output_dir}")
    
    # Filter for the cluster of interest
    print(f"\nFiltering data for cluster: '{target_cluster}'")
    adata_cluster = adata[adata.obs[cluster_column] == target_cluster].copy()

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
    
    # --- Prepare Ranked Gene List ---
    print("Preparing ranked gene list for GSEA...")
    de_results = sc.get.rank_genes_groups_df(adata_cluster, group='Mutant')
    de_results = de_results.copy()
    de_results['pvals_adj'] = de_results['pvals_adj'].replace(0, 1e-300)
    de_results['logfoldchanges'] = de_results['logfoldchanges'].fillna(0)
    de_results['rank_metric'] = np.sign(de_results['logfoldchanges']) * (-np.log10(de_results['pvals_adj']))
    rank_df = de_results.set_index('names')['rank_metric'].sort_values(ascending=False)
    if rank_df.index.duplicated().any():
        rank_df = rank_df[~rank_df.index.duplicated()]
    print(f"Created ranked list with {len(rank_df)} genes.")

    # --- Run GSEA in Parallel ---
    print("\n--- Starting Parallel GSEA Analysis ---")
    
    # Use functools.partial to create a function with fixed arguments
    worker_func = partial(run_gsea_for_gene_set,
                          rank_df=rank_df,
                          organism=ORGANISM,
                          cluster_output_dir=cluster_output_dir,
                          adata_var_names=list(adata.var_names),
                          gsea_padj_threshold=GSEA_PADJ_THRESHOLD,
                          n_top_terms_plot=N_TOP_TERMS_PLOT)

    # Determine number of processes
    num_processes = min(len(GSEA_GENE_SETS), os.cpu_count() - 1 if os.cpu_count() > 1 else 1)
    print(f"Using {num_processes} processes to run GSEA for {len(GSEA_GENE_SETS)} gene sets.")

    with multiprocessing.Pool(processes=num_processes) as pool:
        # pool.map only takes one iterable argument, so we use partial
        all_results_dfs = pool.map(worker_func, GSEA_GENE_SETS)
        
    print("\n--- Parallel GSEA analysis complete. Generating summary files. ---")

    # --- Process and Save Results ---
    summary_results = []
    full_summary_results = []

    for i, df in enumerate(all_results_dfs):
        if df is not None:
            gene_set = GSEA_GENE_SETS[i]
            df['Gene_Set'] = gene_set
            df['Direction'] = df['nes'].apply(lambda x: 'Up' if x > 0 else 'Down')
            full_summary_results.append(df.copy())
            
            if 'fdr' in df.columns:
                sig_results = df[df['fdr'] < GSEA_PADJ_THRESHOLD].copy()
                if not sig_results.empty:
                    summary_results.append(sig_results)

    # --- Save significant results summary ---
    summary_file = cluster_output_dir / "GSEA_Summary_All_Genesets.xlsx"
    if summary_results:
        all_sig_results = pd.concat(summary_results, ignore_index=True)
        all_sig_results.to_excel(summary_file, index=False)
        all_sig_results.to_csv(summary_file.with_suffix('.csv'), index=False)
        print(f"\nSummary table with {len(all_sig_results)} significant terms saved to: {summary_file} and .csv")
    else:
        print("\nNo significant results found. Creating empty summary files.")
        pd.DataFrame().to_excel(summary_file, index=False)
        pd.DataFrame().to_csv(summary_file.with_suffix('.csv'), index=False)
    
    # --- Save full report ---
    full_summary_file = cluster_output_dir / "GSEA_Full_Report_All_Genesets.xlsx"
    if full_summary_results:
        all_full_results = pd.concat(full_summary_results, ignore_index=True)
        all_full_results.to_excel(full_summary_file, index=False)
        all_full_results.to_csv(full_summary_file.with_suffix('.csv'), index=False)
        print(f"\nFull report with all {len(all_full_results)} terms saved to: {full_summary_file} and .csv")
    else:
        print("\nNo full report generated as no result files were found.")

if __name__ == '__main__':
    # This is crucial for multiprocessing to work correctly
    multiprocessing.freeze_support() # For Windows compatibility
    main()