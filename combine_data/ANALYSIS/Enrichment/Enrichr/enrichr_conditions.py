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

# --- Configuration ---
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

ORGANISM = 'Mouse'

# Set up directories
BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw", "enrichr_between_conditions") # New output dir
INPUT_DIR = os.path.join(WORKING_DIR, "results_from_raw")
ADATA_PATH = os.path.join(INPUT_DIR, "annotation_final.h5ad")

# --- Enrichr Specific Configuration ---
CLUSTER_OF_INTEREST = 'Mature GC'
CLUSTER_COLUMN = 'cell_type_L2_new'
GENE_SETS = ['KEGG_2019_Mouse', 'GO_Biological_Process_2023', 'Reactome_2022', 'MSigDB_Hallmark_2020']
PADJ_THRESHOLD = 0.05
LOGFC_THRESHOLD = 0.25
MIN_CELLS = 50

# --- Worker Function for Parallel Enrichr ---

def run_enrichr_for_gene_list(gene_list, gene_set_name, organism, output_dir, direction):
    """
    Performs Enrichr analysis for a single gene list and gene set.
    """
    print(f"--- Starting Enrichr for {direction} genes in gene set: {gene_set_name} ---")
    
    try:
        enr_res = gp.enrichr(gene_list=gene_list,
                             gene_sets=[gene_set_name],
                             organism=organism,
                             outdir=None, # Don't write to disk inside the function
                             verbose=True)
        
        if enr_res.results.empty:
            print(f"No enrichment results for {gene_set_name}.")
            return None

        results_df = enr_res.results
        
        # Filter for significant results before returning
        significant_results = results_df[results_df['Adjusted P-value'] < PADJ_THRESHOLD].copy()
        
        if significant_results.empty:
            print(f"No significant enrichment for {gene_set_name}.")
            return None
            
        significant_results['Gene_Set_Library'] = gene_set_name
        significant_results['Direction'] = direction
        
        print(f"--- Finished Enrichr for gene set: {gene_set_name} ---")
        return significant_results

    except Exception as e:
        print(f"An error occurred during Enrichr for gene set '{gene_set_name}': {e}")
        traceback.print_exc()
        return None

# --- Main Execution Block ---

def main():
    """
    Main function to run the Enrichr analysis pipeline.
    """
    parser = argparse.ArgumentParser(description="Run Enrichr analysis for specific genotypes and cluster merging options.")
    parser.add_argument('--genotype', type=str, required=True, choices=['Emx1', 'Nestin', 'both'],
                        help="Genotype to filter by (e.g., 'Emx1', 'Nestin', 'both').")
    parser.add_argument('--merge_clusters', type=str, required=True, choices=['True', 'False'],
                        help="Whether to merge Immature GC and Mature GC clusters ('True' or 'False').")
    parser.add_argument('--cluster_name', type=str, help="Name of the cluster of interest (e.g., 'GABA'). Overrides merge_clusters.")
    parser.add_argument('--cluster_column', type=str, help="Column in adata.obs for the cluster (e.g., 'cell_type_L1').")
    
    args = parser.parse_args()
    
    genotype_filter = args.genotype
    merge_gc_clusters = args.merge_clusters == 'True'
    cluster_name_arg = args.cluster_name
    cluster_column_arg = args.cluster_column

    # --- Load Data ---
    print(f"Loading AnnData from {ADATA_PATH}...")
    adata = sc.read_h5ad(ADATA_PATH)
    
    # --- Apply Genotype Filter ---
    if genotype_filter != "both":
        adata = adata[adata.obs["genotype"] == genotype_filter].copy()

    if adata.shape[0] == 0:
        print(f"No cells for genotype '{genotype_filter}'. Exiting.")
        sys.exit(0)

    target_cluster = CLUSTER_OF_INTEREST
    cluster_column = CLUSTER_COLUMN

    # --- Cluster Merging/Selection ---
    if cluster_name_arg and cluster_column_arg:
        target_cluster = cluster_name_arg
        cluster_column = cluster_column_arg
    elif merge_gc_clusters:
        # (Same merging logic as before)
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
            target_cluster = new_cluster_name

    # --- Prepare Output Directory ---
    genotype_dir_name = genotype_filter.lower() if genotype_filter != 'both' else 'both_genotypes'
    sanitized_cluster_name = target_cluster.replace(' ', '_').replace('/', '_')
    output_dir = pathlib.Path(BASE_RESULTS_DIR) / genotype_dir_name / f"{sanitized_cluster_name}_Mutant_vs_Control"
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nEnrichr output will be saved to: {output_dir}")

    # --- Filter for Cluster and Genes ---
    adata_cluster = adata[adata.obs[cluster_column] == target_cluster].copy()
    sc.pp.filter_genes(adata_cluster, min_cells=MIN_CELLS)

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

    # --- DEG Analysis ---
    print(f"\nRunning DEG analysis for '{target_cluster}'...")
    sc.tl.rank_genes_groups(adata_cluster, groupby='condition', groups=['Mutant'], reference='Control', method='wilcoxon', use_raw=False, layer=use_layer)
    
    de_results = sc.get.rank_genes_groups_df(adata_cluster, group='Mutant')
    
    # --- Get Up and Down Regulated Gene Lists ---
    degs_sig = de_results[de_results['pvals_adj'] < PADJ_THRESHOLD]
    up_genes = degs_sig[degs_sig['logfoldchanges'] > LOGFC_THRESHOLD]['names'].tolist()
    down_genes = degs_sig[degs_sig['logfoldchanges'] < -LOGFC_THRESHOLD]['names'].tolist()

    print(f"Found {len(up_genes)} upregulated genes and {len(down_genes)} downregulated genes in Mutant vs Control.")

    # --- Run Enrichr Analysis ---
    all_results = []
    
    if up_genes:
        print("\n--- Running Enrichr for Upregulated Genes ---")
        for gene_set in GENE_SETS:
            res = run_enrichr_for_gene_list(up_genes, gene_set, ORGANISM, output_dir, 'Upregulated_in_Mutant')
            if res is not None: all_results.append(res)

    if down_genes:
        print("\n--- Running Enrichr for Downregulated Genes ---")
        for gene_set in GENE_SETS:
            res = run_enrichr_for_gene_list(down_genes, gene_set, ORGANISM, output_dir, 'Upregulated_in_Control')
            if res is not None: all_results.append(res)

    # --- Save Combined Results ---
    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        # Sort by p-value and direction
        final_df = final_df.sort_values(by=['Direction', 'Adjusted P-value'], ascending=[False, True])
        
        summary_file = output_dir / "Enrichr_Summary_All_Genesets.xlsx"
        final_df.to_excel(summary_file, index=False)
        final_df.to_csv(summary_file.with_suffix('.csv'), index=False)
        print(f"\nEnrichr summary table saved to: {summary_file} and .csv")
    else:
        print("\nNo significant enrichment results found.")

if __name__ == '__main__':
    main()