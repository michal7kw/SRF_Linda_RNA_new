import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import pathlib
import numpy as np
import argparse

# --- Configuration ---
PROJECT_DIR = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)

BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")

# Define both processing directories to process
PROCESSING_CONFIGS = [
    {
        'name': 'default',
        'gsea_dir': os.path.join(BASE_RESULTS_DIR, 'gsea_analysis_between_conditions_default'),
        'output_dir': os.path.join(BASE_RESULTS_DIR, 'gsea_summary_visuals_default')
    },
    {
        'name': 'non-default',
        'gsea_dir': os.path.join(BASE_RESULTS_DIR, 'gsea_analysis_between_conditions'),
        'output_dir': os.path.join(BASE_RESULTS_DIR, 'gsea_summary_visuals')
    }
]

# Create output directories
for config in PROCESSING_CONFIGS:
    pathlib.Path(config['output_dir']).mkdir(exist_ok=True)


# Define the analyses to include in the heatmap
# The key is the column name in the heatmap, the value is the directory name for the genotype
ANALYSES = {
    'Emx1': 'emx1',
    'Nestin': 'nestin',
    'Both': 'both_genotypes'
}

def load_gsea_results(analysis_name, genotype_dir, cluster_type, input_file, gsea_base_dir):
    """Loads filtered GSEA results for a given analysis."""
    sanitized_cluster_name = cluster_type.replace(' ', '_')
    comparison_name = "Mutant_vs_Control"

    file_path = pathlib.Path(gsea_base_dir) / genotype_dir / f"{sanitized_cluster_name}_{comparison_name}" / input_file

    if not file_path.exists():
        print(f"Warning: File not found for {analysis_name}: {file_path}")
        return None

    print(f"Loading results for {analysis_name} from {file_path}")
    df = pd.read_excel(file_path)
    
    if 'Term' not in df.columns or 'nes' not in df.columns or 'fdr' not in df.columns:
        print(f"Warning: 'Term', 'nes', or 'fdr' column not found in {file_path}")
        return None
        
    # Handle potential duplicate terms by keeping the one with the max absolute NES
    if df['Term'].duplicated().any():
        print(f"  - Found duplicate terms in {analysis_name}. Resolving by max absolute NES.")
        df['abs_nes'] = df['nes'].abs()
        df = df.loc[df.groupby('Term')['abs_nes'].idxmax()]
        df = df.drop(columns=['abs_nes'])

    # Keep only Term, nes, and fdr, and rename columns
    df = df[['Term', 'nes', 'fdr']].copy()
    df.rename(columns={
        'nes': f"{analysis_name}_nes",
        'fdr': f"{analysis_name}_fdr"
    }, inplace=True)
    df.set_index('Term', inplace=True)
    
    return df

def generate_heatmaps_for_config(cluster_type, input_file, config):
    """
    Generate heatmaps for a specific configuration (default or non-default).
    """
    gsea_base_dir = config['gsea_dir']
    output_dir = config['output_dir']
    config_name = config['name']

    print(f"\n{'='*70}")
    print(f"Processing {config_name} configuration: {cluster_type}")
    print(f"GSEA directory: {gsea_base_dir}")
    print(f"Output directory: {output_dir}")
    print(f"{'='*70}")

    all_results = []
    for analysis_name, genotype_dir in ANALYSES.items():
        results_df = load_gsea_results(analysis_name, genotype_dir, cluster_type, input_file, gsea_base_dir)
        if results_df is not None:
            all_results.append(results_df)

    if not all_results:
        print(f"No GSEA result files found for cluster type '{cluster_type}' in {config_name}. Skipping.")
        return

    # Merge all dataframes on the 'Term' index
    merged_df = pd.concat(all_results, axis=1, join='outer')

    if merged_df.empty:
        print(f"No GSEA result files found or they were empty for {config_name}. Skipping.")
        return

    print(f"Found {len(merged_df)} total terms across all analyses.")

    # --- Term Selection for Plotting ---
    nes_cols = [f"{name}_nes" for name in ANALYSES.keys()]
    merged_df['mean_nes'] = merged_df[nes_cols].mean(axis=1)

    up_regulated = merged_df[merged_df['mean_nes'] > 0].sort_values('mean_nes', ascending=False)
    down_regulated = merged_df[merged_df['mean_nes'] <= 0].sort_values('mean_nes', ascending=True)

    print(f"Found {len(up_regulated)} up-regulated terms and {len(down_regulated)} down-regulated terms.")

    # top_up = up_regulated.head(15)
    # top_down = down_regulated.head(15)
    top_up = up_regulated
    top_down = down_regulated

    combined_terms = pd.concat([top_up, top_down])
    combined_terms.sort_values(by=['mean_nes'], ascending=False, inplace=True)

    plot_df = combined_terms[nes_cols].copy()
    plot_df.columns = ANALYSES.keys()

    if plot_df.empty:
        print(f"No data to plot after filtering for {config_name}. Skipping.")
        return

    # --- Plotting with FDR ---
    # Create Annotation Dataframe with FDR values
    annot_df = pd.DataFrame(index=plot_df.index, columns=plot_df.columns)
    for term in plot_df.index:
        for analysis in ANALYSES.keys():
            fdr_val = merged_df.loc[term, f"{analysis}_fdr"]
            nes_val = merged_df.loc[term, f"{analysis}_nes"]
            if pd.notna(fdr_val):
                annot_df.loc[term, analysis] = f"{nes_val:.2f}; {fdr_val:.2f}"
            else:
                annot_df.loc[term, analysis] = ""

    # --- Plotting ---
    vmax = np.nanmax(np.abs(plot_df.values))
    vmax = max(vmax, 0.1)
    vmin = -vmax

    plt.figure(figsize=(14, max(8, len(plot_df) * 0.4)))

    ax = sns.heatmap(
        plot_df,
        cmap='RdBu_r',
        vmin=vmin,
        vmax=vmax,
        annot=annot_df,
        fmt="",
        linewidths=.5,
        cbar_kws={'label': 'Normalized Enrichment Score (NES)'}
    )

    ax.set_title(f'GSEA Focused Results Comparison ({cluster_type.replace("_", " ")})', fontsize=16, pad=20)
    ax.set_xlabel('Genotype Comparison (Mutant vs Control)', fontsize=12)
    ax.set_ylabel('Gene Set Term', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(fontsize=10)

    # Adjust layout to make space for the colorbar and prevent labels from being cut off
    plt.tight_layout(rect=[0, 0, 0.9, 1]) # Adjust rect to make space for colorbar

    # Save the plot
    output_filename = f"GSEA_heatmap_{cluster_type}.png"
    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

    print(f"\nHeatmap saved to: {output_path}")
    plt.close()

    # Generate no-FDR version
    print("\n--- Generating no-FDR version ---")

    # Create annotation without FDR
    annot_df_no_fdr = pd.DataFrame(index=plot_df.index, columns=plot_df.columns)
    annot_df_no_fdr = annot_df_no_fdr.fillna("")

    # Create new plot
    plt.figure(figsize=(14, max(8, len(plot_df) * 0.4)))

    ax = sns.heatmap(
        plot_df,
        cmap='RdBu_r',
        vmin=vmin,
        vmax=vmax,
        annot=annot_df_no_fdr,
        fmt="",
        linewidths=.5,
        cbar_kws={'label': 'Normalized Enrichment Score (NES)'}
    )

    ax.set_title(f'GSEA Focused Results Comparison ({cluster_type.replace("_", " ")})', fontsize=16, pad=20)
    ax.set_xlabel('Genotype Comparison (Mutant vs Control)', fontsize=12)
    ax.set_ylabel('Gene Set Term', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(fontsize=10)

    plt.tight_layout(rect=[0, 0, 0.9, 1])

    output_filename_no_fdr = f"GSEA_heatmap_{cluster_type}_no_fdr.png"
    output_path_no_fdr = os.path.join(output_dir, output_filename_no_fdr)
    plt.savefig(output_path_no_fdr, dpi=300, bbox_inches='tight')

    print(f"No-FDR heatmap saved to: {output_path_no_fdr}")
    plt.close()


def main():
    """
    Generates heatmaps comparing GSEA results across different genotypes.
    Processes BOTH default and non-default configurations automatically.
    """
    parser = argparse.ArgumentParser(description="Generate heatmaps to compare GSEA results across different analyses.")
    parser.add_argument('--cluster_type', type=str, required=True, choices=['Mature_GC', 'Combined_GC', 'GABA'],
                        help="The cluster type to generate the heatmap for ('Mature_GC' or 'Combined_GC' or 'GABA').")
    parser.add_argument('--input_file', type=str, required=True, default="GSEA_Focus_Neuro_CellDeath.xlsx")

    args = parser.parse_args()
    cluster_type = args.cluster_type
    input_file = args.input_file

    print("="*70)
    print(f"GSEA Heatmap Generator for cluster type: {cluster_type}")
    print("Processing BOTH default and non-default configurations")
    print("="*70)

    # Process both configurations
    for config in PROCESSING_CONFIGS:
        generate_heatmaps_for_config(cluster_type, input_file, config)

    print("\n" + "="*70)
    print("All heatmap generation completed!")
    print("="*70)


if __name__ == '__main__':
    main()