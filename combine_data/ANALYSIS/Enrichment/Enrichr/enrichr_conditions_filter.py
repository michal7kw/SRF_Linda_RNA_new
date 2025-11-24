# %%
import pandas as pd
import os
import pathlib
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import json
import argparse

# --- Configuration ---
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

# --- Enrichr Results Location ---
BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")
EXPERIMENT = "enrichr_between_conditions"

# --- Configuration to match the main analysis script ---
CLUSTER_OF_INTEREST = 'Mature GC'

# --- Focused Analysis Configuration ---
FOCUS_CONFIG = {
    'neuro': {
        'keywords': [
            'neurodegeneration', 'neuronal death', 'axonopathy', 'neuron death',
            'apoptosis', 'apoptotic', 'necroptosis', 'pyroptosis', 'ferroptosis',
            'autophagy', 'autophagic', 'parthanatos', 'cell death'
        ],
        'file_suffix': 'Neuro_CellDeath',
        'title_name': 'Neurodegeneration & Cell Death'
    },
    'splicing': {
        'keywords': [
            'splicing', 'spliceosome', 'splice', 'rna processing',
            'alternative splicing', 'splice variant', 'isoform',
            'pre-mrna', 'intron', 'exon',
            'splice site', 'exon skipping', 'intron retention',
            'sr protein', 'hnrnp', 'splicing factor',
            'aberrant splicing', 'spliceopathy',
            'nonsense-mediated decay', 'nmd'
        ],
        'file_suffix': 'RNA_Splicing',
        'title_name': 'RNA Splicing'
    }
}
# --- End Configuration ---

def main():
    parser = argparse.ArgumentParser(description="Filter Enrichr analysis results for specific conditions and genotypes.")
    parser.add_argument('--genotype', type=str, required=True, choices=['Emx1', 'Nestin', 'both'],
                        help="Genotype used for the Enrichr analysis.")
    parser.add_argument('--merge_clusters', type=str, required=True, choices=['True', 'False'],
                        help="Whether clusters were merged. Ignored if --cluster_name is provided.")
    parser.add_argument('--cluster_name', type=str, help="Name of the cluster of interest. Overrides merge_clusters.")
    parser.add_argument('--focus', type=str, required=True, choices=list(FOCUS_CONFIG.keys()),
                        help="The focus area for the analysis.")
    
    args = parser.parse_args()
    
    genotype_filter = args.genotype
    merge_gc_clusters = args.merge_clusters == 'True'
    cluster_name_arg = args.cluster_name
    focus_key = args.focus

    focus_config = FOCUS_CONFIG[focus_key]
    search_keywords = focus_config['keywords']
    file_suffix = focus_config['file_suffix']
    title_name = focus_config['title_name']

    # Determine the correct cluster name
    if cluster_name_arg:
        effective_cluster_name = cluster_name_arg
    elif merge_gc_clusters:
        effective_cluster_name = 'Combined GC'
    else:
        effective_cluster_name = CLUSTER_OF_INTEREST

    # Construct input directory
    genotype_dir_name = genotype_filter.lower() if genotype_filter != 'both' else 'both_genotypes'
    sanitized_cluster_name = effective_cluster_name.replace(' ', '_').replace('/', '_')
    comparison_name = "Mutant_vs_Control"
    enrichr_output_dir = pathlib.Path(BASE_RESULTS_DIR) / EXPERIMENT / genotype_dir_name / f"{sanitized_cluster_name}_{comparison_name}"
    
    summary_file = enrichr_output_dir / "Enrichr_Summary_All_Genesets.xlsx"
    focused_summary_file = enrichr_output_dir / f"Enrichr_Focus_{file_suffix}.xlsx"
    
    print(f"--- Focused Analysis for {title_name} on Enrichr Results ---")
    print(f"Looking for summary file: {summary_file}")

    if not summary_file.exists():
        print(f"\nError: The Enrichr summary file was not found at {summary_file}.")
        print(f"Please run 'enrichr_conditions.py' for the corresponding settings first.")
        sys.exit(1)
    
    # --- Load and Filter Results ---
    try:
        all_results = pd.read_excel(summary_file)
        if all_results.empty or 'Term' not in all_results.columns:
            print("\nThe Enrichr report file is empty or invalid.")
            sys.exit(0)
        
        print(f"Successfully loaded {len(all_results)} total terms.")
        
        pattern = '|'.join(search_keywords)
        print(f"Searching for terms matching: {pattern}")
        
        focused_results = all_results[all_results['Term'].str.contains(pattern, case=False, na=False)].copy()
        
        if not focused_results.empty:
            print(f"\nFound {len(focused_results)} terms related to {title_name}.")
            
            focused_results.sort_values(by=['Direction', 'Adjusted P-value'], ascending=[False, True], inplace=True)
            
            focused_results.to_excel(focused_summary_file, index=False)
            focused_results.to_csv(focused_summary_file.with_suffix('.csv'), index=False)
            print(f"Focused summary saved to: {focused_summary_file} and .csv")
            
            print(f"\n--- Top Hits for {title_name} ---")
            for _, row in focused_results.head(20).iterrows():
                print(f"  - Term: {row['Term']}")
                print(f"    Adj P-val: {row['Adjusted P-value']:.4f}, Direction: {row['Direction']}, Gene Set: {row['Gene_Set_Library']}")
        else:
            print("\nNo significant terms found matching the specified keywords.")
    
    except Exception as e:
        print(f"An error occurred during the analysis: {e}")
        traceback.print_exc()

    # --- Visualization of Focused Results ---
    if 'focused_results' in locals() and not focused_results.empty:
        print("\n--- Visualizing Focused Enrichr Results ---")

        # The input is already significant, so we just take the top terms for the plot
        plot_data = focused_results.head(20).copy().iloc[::-1]
        
        plot_data['-log10(Adj. P-value)'] = -plot_data['Adjusted P-value'].apply(lambda x: np.log10(x) if x > 0 else np.log10(1e-300))
        
        color_map = {'Upregulated_in_Mutant': 'red', 'Upregulated_in_Control': 'blue'}

        plt.figure(figsize=(12, 10))
        
        ax = sns.barplot(
            x='-log10(Adj. P-value)',
            y='Term',
            data=plot_data,
            hue='Direction',
            palette=color_map,
            orient='h',
            dodge=False
        )
        
        plt.xlabel('-log10(Adjusted P-value)', fontsize=14)
        plt.ylabel('Gene Set Term', fontsize=14)
        
        plot_title = f'Top 20 Enrichr Terms in {effective_cluster_name}\n({title_name} Focus)'
        plt.title(plot_title, fontsize=16, pad=20)
        plt.grid(axis='x', linestyle='--', alpha=0.7)
        
        significance_threshold = -np.log10(0.05)
        ax.axvline(x=significance_threshold, color='green', linestyle='--', linewidth=2)
        
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='blue', label='Upregulated in Mutant'),
            Patch(facecolor='red', label='Upregulated in Control'),
            plt.Line2D([0], [0], color='green', linestyle='--', lw=2, label='Adj. P-value = 0.05')
        ]
        ax.legend(handles=legend_elements, title='Legend', loc='lower right')

        plt.tight_layout()

        plot_file = focused_summary_file.with_suffix('.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to: {plot_file}")
        
    else:
        print("\nNo focused results to visualize.")

if __name__ == '__main__':
    main()