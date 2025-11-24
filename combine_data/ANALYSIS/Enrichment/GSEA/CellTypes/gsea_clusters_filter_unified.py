# %%
import pandas as pd
import os
import pathlib
import sys
import argparse

# --- Configuration ---
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

# --- GSEA Results Location ---
BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")

# --- Configuration to match the main analysis script ---
# IMPORTANT: This must match the setting in the script that generated the results.
CLUSTER_OF_INTEREST = 'Mature GC' # The *original* cluster name.

# --- Focused Analysis Configuration ---
# Keywords to search for in GSEA results
SEARCH_KEYWORDS = [
    'neurodegeneration', 'neuronal death', 'axonopathy', 'neuron death',
    'apoptosis', 'apoptotic', 'necroptosis', 'pyroptosis', 'ferroptosis', 
    'autophagy', 'autophagic', 'parthanatos', 'cell death'
]
# --- End Configuration ---

def main():
    parser = argparse.ArgumentParser(description="Filter GSEA analysis results for specific conditions and genotypes.")
    parser.add_argument('--condition', type=str, required=True, choices=['Control', 'Mutant', 'both'],
                        help="Condition used for the GSEA analysis (e.g., 'Control', 'Mutant', 'both').")
    parser.add_argument('--genotype', type=str, required=True, choices=['Emx1', 'Nestin', 'both'],
                        help="Genotype used for the GSEA analysis (e.g., 'Emx1', 'Nestin', 'both').")
    parser.add_argument('--merge_clusters', type=str, required=True, choices=['True', 'False'],
                        help="Whether clusters were merged in the GSEA analysis ('True' or 'False').")
    
    args = parser.parse_args()
    
    condition_filter = args.condition
    genotype_filter = args.genotype
    merge_gc_clusters = args.merge_clusters == 'True'

    # Determine the correct cluster name based on whether merging was done
    if merge_gc_clusters:
        effective_cluster_name = 'Combined GC'
        print("Assuming clusters were merged. Looking for results in 'Combined_GC' directory.")
    else:
        effective_cluster_name = CLUSTER_OF_INTEREST

    # Construct input GSEA_OUTPUT_DIR based on arguments
    condition_dir_name = 'ctrl' if condition_filter == 'Control' else 'mut' if condition_filter == 'Mutant' else 'both_conditions'
    genotype_dir_name = genotype_filter.lower() if genotype_filter != 'both' else 'both_genotypes'
    
    GSEA_OUTPUT_DIR = os.path.join(BASE_RESULTS_DIR, 'gsea_analysis_between_clusters', condition_dir_name, genotype_dir_name)

    # Construct paths
    sanitized_cluster_name = effective_cluster_name.replace(' ', '_').replace('/', '_')
    cluster_output_dir = pathlib.Path(GSEA_OUTPUT_DIR) / sanitized_cluster_name
    # Point to the new, full report file for the focus analysis
    summary_file = cluster_output_dir / "GSEA_Full_Report_All_Genesets.xlsx"
    focused_summary_file = cluster_output_dir / "GSEA_Focus_Neuro_CellDeath.xlsx"
    
    print("--- Focused Analysis for Neurodegeneration and Cell Death ---")
    print(f"This script expects GSEA results from 'gsea_clusters_unified.py' to be present.")
    print(f"Looking for full report file: {summary_file}")
    
    # Check if the summary file from the main analysis exists
    if not summary_file.exists():
        print(f"\nError: The GSEA full report file was not found at {summary_file}.")
        print(f"Please run 'gsea_clusters_unified.py' for condition '{condition_filter}' and genotype '{genotype_filter}' first to generate the necessary results.")
        sys.exit(1)
    
    # --- Load and Filter Results ---
    try:
        print(f"\nLoading full report from {summary_file}...")
        all_results = pd.read_excel(summary_file)
    
        # Check if the dataframe is empty or lacks the required 'Term' column.
        if all_results.empty or 'Term' not in all_results.columns:
            print("\nThe GSEA report file is empty or invalid. No terms to analyze.")
            print("This is unexpected if the main GSEA script ran correctly.")
            sys.exit(0)
    
        print(f"Successfully loaded {len(all_results)} total terms from the report.")
    
        # Create a regex pattern to find any of the keywords (case-insensitive)
        pattern = '|'.join(SEARCH_KEYWORDS)
        print(f"Searching for terms matching: {pattern}")
        
        # Filter the 'Term' column
        focused_results = all_results[all_results['Term'].str.contains(pattern, case=False, na=False)].copy()
        
        if not focused_results.empty:
            print(f"\nFound {len(focused_results)} terms related to neurodegeneration and cell death.")
            
            # Sort by significance and effect size for better ranking
            focused_results['abs_nes'] = focused_results['nes'].abs()
            focused_results.sort_values(by=['fdr', 'abs_nes'], ascending=[True, False], inplace=True)
            
            # Save to new files (Excel and CSV)
            focused_results.to_excel(focused_summary_file, index=False)
            focused_results.to_csv(focused_summary_file.with_suffix('.csv'), index=False)
            print(f"Focused summary saved to: {focused_summary_file} and .csv")
            
            # Print the focused results to the console
            print("\n--- Top Hits for Neurodegeneration & Cell Death ---")
            for _, row in focused_results.head(20).iterrows(): # Print top 20 hits
                print(f"  - Term: {row['Term']}")
                print(f"    NES: {row['nes']:.2f}, FDR: {row['fdr']:.4f}, Gene Set: {row['Gene_Set']}, Direction: {row['Direction']}")
    
        else:
            print("\nNo significant terms found matching the specified keywords.")
    
    except Exception as e:
        print(f"An error occurred during the analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    main()

# --- Visualization of Focused Results ---
if 'focused_results' in locals() and not focused_results.empty:
    print("\n--- Visualizing Focused GSEA Results ---")

    # Take the top 20 terms for visualization and reverse for plotting
    # (so the most significant is at the top)
    plot_data = focused_results.head(20).copy().iloc[::-1]
    
    # Create a new column for -log10(FDR) for plotting
    # Add a small constant to FDR to avoid log(0) to prevent -inf values
    plot_data['-log10(FDR)'] = -plot_data['fdr'].apply(lambda x: np.log10(x) if x > 0 else np.log10(1e-300))

    # Create a continuous color mapping based on NES values
    # Use a diverging colormap from red (negative NES) to blue (positive NES)
    nes_values = plot_data['nes'].values
    nes_min, nes_max = nes_values.min(), nes_values.max()
    
    # Normalize NES values to [0, 1] for color mapping
    # Center the colormap around 0
    nes_abs_max = max(abs(nes_min), abs(nes_max))
    normalized_nes = (nes_values + nes_abs_max) / (2 * nes_abs_max)
    
    # Create color mapping using a diverging colormap
    import matplotlib.cm as cm
    colormap = cm.RdBu_r  # Red-Blue diverging colormap (reversed)
    bar_colors = [colormap(norm_val) for norm_val in normalized_nes]

    plt.figure(figsize=(12, 10))
    
    # Create the horizontal bar plot
    ax = sns.barplot(
        x='-log10(FDR)',
        y='Term',
        data=plot_data,
        palette=bar_colors,
        orient='h'
    )
    
    plt.xlabel('-log10(FDR q-value)', fontsize=14)
    plt.ylabel('Gene Set Term', fontsize=14)
    
    # --- Create Dynamic Plot Title ---
    plot_title = f'Top 20 GSEA Terms in {effective_cluster_name}\n(Neurodegeneration & Cell Death Focus)'
    # The cluster_size_info is not available in this script, so we omit it from the title.
    # if cluster_size_info:
    #     mutant_n = cluster_size_info.get('mutant_count', "N/A")
    #     control_n = cluster_size_info.get('control_count', "N/A")
    #     plot_title += f"\n(Mutant n={mutant_n} vs Control n={control_n})"

    plt.title(plot_title, fontsize=16, pad=20)
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    
    # Add a vertical line for the significance threshold (FDR = 0.25)
    significance_threshold = -np.log10(0.25)
    ax.axvline(x=significance_threshold, color='green', linestyle='--', linewidth=2)
    
    # Add a custom legend to explain colors and the significance line
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    
    # Create a colorbar to show the NES scale
    sm = ScalarMappable(cmap=colormap, norm=Normalize(vmin=-nes_abs_max, vmax=nes_abs_max))
    sm.set_array([])
    
    legend_elements = [
        Line2D([0], [0], color='green', linestyle='--', linewidth=2, label='FDR = 0.25 Threshold')
    ]
    
    # Add colorbar
    cbar = plt.colorbar(sm, ax=ax, shrink=0.8, aspect=20)
    cbar.set_label('Normalized Enrichment Score (NES)', fontsize=12)
    
    ax.legend(handles=legend_elements, title='Legend', loc='lower right')

    # Adjust layout to make sure labels are not cut off
    plt.tight_layout()

    # Save the plot
    plot_file = focused_summary_file.with_suffix('.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to: {plot_file}")
    
    # Display the plot
    # plt.show()

else:
    print("\nNo focused results to visualize.")