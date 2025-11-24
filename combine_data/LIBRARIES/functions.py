import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import anndata as ad
from tqdm import tqdm
import argparse
import seaborn as sns

def load_annotated_data(adata_paths):
    """Load annotated AnnData objects from paths"""
    adata_dict = {}
    for sample_name, file_path in adata_paths.items():
        print(f"Loading {sample_name} from {file_path}")
        try:
            adata = sc.read_h5ad(file_path)
            # Add sample identifier
            adata.obs['sample'] = sample_name
            # Add condition and genotype metadata
            adata.obs['condition'] = 'Control' if 'Ctrl' in sample_name else 'Mutant'
            adata.obs['genotype'] = sample_name.split('_')[0]  # Emx1 or Nestin
            adata_dict[sample_name] = adata
            print(f"  Loaded {adata.n_obs} cells and {adata.n_vars} genes")
        except Exception as e:
            print(f"Error loading {sample_name}: {e}")
    return adata_dict

def get_common_genes(adata_dict):
    """Find genes common to all datasets"""
    if not adata_dict:
        return []

    # Get sets of gene names from each dataset
    gene_sets = [set(adata.var_names) for adata in adata_dict.values()]

    # Find intersection of all gene sets
    common_genes = set.intersection(*gene_sets)

    print(f"Found {len(common_genes)} genes common to all datasets")
    return list(common_genes)

def merge_datasets(adata_dict, common_genes=None):
    """Merge datasets using only common genes"""
    if not adata_dict:
        print("No datasets to merge")
        return None

    # List of datasets to merge
    adatas_to_merge = []

    # Process each dataset
    for sample_name, adata in adata_dict.items():
        # Make a copy to avoid modifying the original
        adata_copy = adata.copy()

        # Subset to common genes if provided
        if common_genes:
            adata_copy = adata_copy[:, common_genes].copy()

        # Convert any layer matrices to dense if needed
        for layer_key in adata_copy.layers:
            if not isinstance(adata_copy.layers[layer_key], np.ndarray):
                adata_copy.layers[layer_key] = adata_copy.layers[layer_key].toarray()

        # Add to list for merging
        adatas_to_merge.append(adata_copy)

    # Merge datasets
    print("Merging datasets...")
    adata_merged = adatas_to_merge[0].concatenate(
        adatas_to_merge[1:],
        join='inner',  # Use only features found in all objects
        batch_key='batch',  # Add batch annotation
        batch_categories=[sample for sample in adata_dict.keys()],
        index_unique='-'  # Add suffix to cell names to make them unique
    )

    print(f"Merged dataset: {adata_merged.n_obs} cells, {adata_merged.n_vars} genes")
    return adata_merged

def process_merged_data(adata_merged):
    """Perform basic processing on merged data"""
    # Create a copy of the raw counts for later use
    # adata_merged.layers['raw_counts'] = adata_merged.raw.X.copy()

    # Calculate highly variable genes
    print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(adata_merged, min_mean=0.0125, max_mean=3, min_disp=0.5)
    print(f"Found {sum(adata_merged.var.highly_variable)} highly variable genes")

    # Run PCA
    print("Running PCA...")
    sc.pp.scale(adata_merged, max_value=10)
    sc.tl.pca(adata_merged, svd_solver='arpack')

    # Run UMAP and clustering
    print("Computing neighborhood graph...")
    sc.pp.neighbors(adata_merged, n_neighbors=15, n_pcs=30)

    print("Running UMAP...")
    sc.tl.umap(adata_merged)

    print("Running Leiden clustering...")
    # Try with default flavor and more iterations to potentially improve stability
    for resolution in [0.3, 0.5, 0.8]:
        print(f"  Running Leiden with resolution {resolution}...")
        try:
            sc.tl.leiden(adata_merged, resolution=resolution, key_added=f'leiden_{resolution}', n_iterations=5, random_state=0) # Added random_state for consistency
        except ValueError as e:
            print(f"    Error during Leiden clustering with resolution {resolution}: {e}")
            # Optionally, decide if you want to continue with other resolutions or stop
            # raise e # Uncomment to stop execution on error
            continue # Continue to the next resolution if one fails

    return adata_merged

def analyze_merged_data(adata_merged, output_dir, display_figures=False):
    """Generate analysis plots for merged data

    Parameters
    ----------
    adata_merged : AnnData
        Merged AnnData object
    output_dir : str
        Directory to save results
    display_figures : bool
        Whether to display figures in notebook/display environment
    """
    os.makedirs(output_dir, exist_ok=True)

    # Create plots directory
    plots_dir = os.path.join(output_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    # Plot sample distribution
    sc.pl.umap(adata_merged, color=['sample', 'condition', 'genotype'],
              wspace=0.4, ncols=2, save=None if display_figures else '_samples_merged.png')

    # Plot cell type distributions from both models
    if 'DG_majority_voting' in adata_merged.obs.columns:
        sc.pl.umap(adata_merged, color=['DG_majority_voting'],
                  save=None if display_figures else '_DG_celltypes_merged.png')
    if 'ISO_majority_voting' in adata_merged.obs.columns:
        sc.pl.umap(adata_merged, color=['ISO_majority_voting'],
                  save=None if display_figures else '_ISO_celltypes_merged.png')

    # Plot clustering results for specific resolutions
    leiden_cols = [f'leiden_{res}' for res in [0.3, 0.5, 0.8] if f'leiden_{res}' in adata_merged.obs.columns]
    sc.pl.umap(adata_merged, color=leiden_cols, wspace=0.4, ncols=len(leiden_cols),
               save=None if display_figures else '_clusters_merged.png')

    # Sample count statistics
    sample_counts = adata_merged.obs['sample'].value_counts()
    condition_counts = adata_merged.obs['condition'].value_counts()
    genotype_counts = adata_merged.obs['genotype'].value_counts()

    # Create summary dataframe
    summary_df = pd.DataFrame({
        'Sample': sample_counts.index,
        'Cells': sample_counts.values,
        'Percentage': (sample_counts.values / sample_counts.sum() * 100).round(2)
    })

    # Save summary statistics
    summary_df.to_csv(os.path.join(output_dir, 'sample_statistics.csv'), index=False)

    # Analyze cell types by condition for both models
    for model in ['DG', 'ISO']:
        col_name = f'{model}_majority_voting'
        if col_name in adata_merged.obs.columns:
            # Contingency table of cell types by condition
            celltype_by_condition = pd.crosstab(
                adata_merged.obs[col_name],
                adata_merged.obs['condition'],
                normalize='index'
            ) * 100  # Convert to percentage

            # Save cell type statistics
            celltype_by_condition.to_csv(os.path.join(output_dir, f'{model}_celltype_by_condition.csv'))

            # Plot cell type distribution
            plt.figure(figsize=(14, 10))
            celltype_by_condition.plot(kind='bar', stacked=False)
            plt.title(f'{model} Cell Type Distribution by Condition')
            plt.ylabel('Percentage')
            plt.xticks(rotation=90)
            plt.tight_layout()

            if not display_figures:
                plt.savefig(os.path.join(plots_dir, f'{model}_celltype_condition_distribution.png'), dpi=150)
                plt.close()

    print(f"Analysis complete. Results saved to {output_dir}")
    return

def compare_annotations_with_leiden(adata, leiden_key, models, output_dir, annotation_suffix='majority_voting'):
    """
    Compares CellTypist annotations with Leiden clusters and generates stacked bar plots and heatmaps.

    Args:
        adata (anndata.AnnData): AnnData object containing the data.
        leiden_key (str): Key in adata.obs containing Leiden cluster assignments.
        models (dict): Dictionary of models to compare against.
        output_dir (str): Directory to save the plots.
        annotation_suffix (str): Suffix for the annotation column (default: 'majority_voting').
    """
    print(f"\nComparing CellTypist annotations with Leiden clusters ('{leiden_key}')...")
    if leiden_key not in adata.obs:
        print(f"Leiden key '{leiden_key}' not found. Skipping comparison with Leiden clusters.")
    else:
        for model_name in models.keys():
            prefix = 'ISO_' if 'Isocortex' in model_name else 'DG_'
            pred_col = f'{prefix}{annotation_suffix}'

            if pred_col in adata.obs.columns:
                print(f"\n--- Comparison for {model_name} vs {leiden_key} --- ")
                # Create a cross-tabulation (contingency table)
                ct = pd.crosstab(adata.obs[leiden_key], adata.obs[pred_col])

                # Plot as a stacked bar chart (showing proportion)
                try:
                    ct_prop = ct.apply(lambda r: r/r.sum()*100, axis=1) # Normalize each Leiden cluster to 100%

                    fig, ax = plt.subplots(figsize=(12, 8))
                    ct_prop.plot(kind='bar', stacked=True, ax=ax, colormap='tab20') # Use a colormap like tab20
                    plt.title(f'{leiden_key} Composition by {model_name} Predictions (Raw Workflow)')
                    plt.xlabel(f'{leiden_key} Cluster')
                    plt.ylabel('Percentage of Cells')
                    plt.legend(title=f'{model_name} Labels (MV)', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
                    plt.tight_layout(rect=[0, 0, 0.85, 1]) # Adjust layout to make space for legend
                    # stackedbar_file = os.path.join(output_dir, "plots", f"stackedbar_{leiden_key}_vs_{prefix}raw.png")
                    # plt.savefig(stackedbar_file, dpi=150)
                    # plt.close(fig)
                    plt.show()
                    print(f"Generated stacked bar plot")
                except Exception as e:
                    print(f"Could not generate stacked bar plot for {model_name}: {e}")

                # Optional: Heatmap of the proportions
                try:
                    plt.figure(figsize=(max(8, ct_prop.shape[1]*0.5), max(6, ct_prop.shape[0]*0.4))) # Adjust size based on dimensions
                    sns.heatmap(ct_prop, annot=True, fmt=".1f", cmap="viridis", linewidths=.5)
                    plt.title(f'Heatmap: {leiden_key} vs {model_name} Predictions (%) (Raw Workflow)')
                    plt.xlabel(f'{model_name} Predicted Label')
                    plt.ylabel(f'{leiden_key} Cluster')
                    plt.tight_layout()
                    # heatmap_leiden_file = os.path.join(output_dir, "plots", f"heatmap_{leiden_key}_vs_{prefix}raw.png")
                    # plt.savefig(heatmap_leiden_file, dpi=150)
                    # plt.close()
                    plt.show()
                    print(f"Generated Leiden comparison heatmap")
                except Exception as e:
                    print(f"Could not generate Leiden comparison heatmap for {model_name}: {e}")
            else:
                print(f"Warning: Predicted label column '{pred_col}' not found for model {model_name}. Skipping Leiden comparison.")
