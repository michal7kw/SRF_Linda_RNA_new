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
import argparse

# --- Configuration ---
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

ORGANISM = 'Mouse'

# Set up directories
BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")
INPUT_DIR = BASE_RESULTS_DIR
ADATA_PATH = os.path.join(INPUT_DIR, "annotation_final.h5ad")

# --- GSEA Specific Configuration ---
CLUSTER_OF_INTEREST = 'Mature GC' # Initial cluster, will be updated if merging
CLUSTER_COLUMN = 'cell_type_L2_new'
GSEA_GENE_SETS = ['KEGG_2019_Mouse', 'GO_Biological_Process_2023', 'Reactome_2022', 'MSigDB_Hallmark_2020']
N_TOP_TERMS_PLOT = 15
GSEA_PADJ_THRESHOLD = 0.25

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
                            sanitized_term_name = str(term_name).replace(' ', '_').replace('/', '_')[:100]
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
    parser = argparse.ArgumentParser(description="Run GSEA analysis for specific conditions and genotypes.")
    parser.add_argument('--condition', type=str, required=True, choices=['Control', 'Mutant', 'both'],
                        help="Condition to filter by (e.g., 'Control', 'Mutant', 'both').")
    parser.add_argument('--genotype', type=str, required=True, choices=['Emx1', 'Nestin', 'both'],
                        help="Genotype to filter by (e.g., 'Emx1', 'Nestin', 'both').")
    parser.add_argument('--merge_clusters', type=str, required=True, choices=['True', 'False'],
                        help="Whether to merge Immature GC and Mature GC clusters ('True' or 'False').")
    
    args = parser.parse_args()
    
    condition_filter = args.condition
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
        
    # --- Apply Filters ---
    if condition_filter != "both":
        print(f"Filtering for condition: {condition_filter}")
        adata = adata[adata.obs["condition"] == condition_filter]
    else:
        print("No condition filter applied (using 'both').")

    if genotype_filter != "both":
        print(f"Filtering for genotype: {genotype_filter}")
        adata = adata[adata.obs["genotype"] == genotype_filter]
    else:
        print("No genotype filter applied (using 'both').")

    if adata.shape[0] == 0:
        print(f"No cells found for condition '{condition_filter}' and genotype '{genotype_filter}'. Exiting.")
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

    # --- Prepare for GSEA ---
    # Construct output directory based on arguments
    condition_dir_name = 'ctrl' if condition_filter == 'Control' else 'mut' if condition_filter == 'Mutant' else 'both_conditions'
    genotype_dir_name = genotype_filter.lower() if genotype_filter != 'both' else 'both_genotypes'
    
    GSEA_OUTPUT_DIR_BASE_DYNAMIC = os.path.join(BASE_RESULTS_DIR, 'gsea_analysis_between_clusters', condition_dir_name, genotype_dir_name)

    sanitized_cluster_name = target_cluster.replace(' ', '_').replace('/', '_')
    cluster_output_dir = pathlib.Path(GSEA_OUTPUT_DIR_BASE_DYNAMIC) / sanitized_cluster_name
    cluster_output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nGSEA output will be saved to: {cluster_output_dir}")
    
    # --- Differential Expression Analysis ---
    print(f"\nRunning differential expression for '{target_cluster}' vs all other cells...")
    sc.tl.rank_genes_groups(adata, groupby=CLUSTER_COLUMN, groups=[target_cluster], reference='rest', method='t-test', use_raw=False)
    
    # --- Prepare Ranked Gene List ---
    print("Preparing ranked gene list for GSEA...")
    de_results = sc.get.rank_genes_groups_df(adata, group=target_cluster)
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