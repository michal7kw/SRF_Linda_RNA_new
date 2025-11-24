# %%
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import pathlib
import gseapy as gp

# --- Configuration ---
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

ORGANISM = 'Mouse'

# Set up directories
REMOVE_DOUBLETS = False
FIX_TRESHOLD = True

if FIX_TRESHOLD:
    BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")
else:
    if REMOVE_DOUBLETS:
        BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw_percentile_threshold", "doublets_removed")
    else:
        BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw_percentile_threshold")

INPUT_DIR = BASE_RESULTS_DIR
ADATA_PATH = os.path.join(INPUT_DIR, "annotation_final.h5ad")

DGE_RESULTS_DIR_BASE = os.path.join(INPUT_DIR,"DEGs_cell_type_L2")

CUSTOM_ANALYSIS =  "FC_0_25"

DEG_BY = "cell_type_L2"
DEGs_folder = "DEGs_cell_type_L2FC_0_25"

if CUSTOM_ANALYSIS is not None:
    DGE_RESULTS_DIR = DGE_RESULTS_DIR_BASE + CUSTOM_ANALYSIS
else:
    DGE_RESULTS_DIR = DGE_RESULTS_DIR_BASE

DGE_RESULTS_DIR = os.path.join(DGE_RESULTS_DIR,"biomarkers")
GO_OUTPUT_DIR = os.path.join(DGE_RESULTS_DIR, 'go_enrichment_combined') # New output directory
GENE_SETS = ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023']

N_TOP_TERMS_PLOT = 10
GO_PADJ_THRESHOLD = 0.05

# --- End Configuration ---

# Create base output directory if it doesn't exist
GO_OUTPUT_DIR = pathlib.Path(GO_OUTPUT_DIR)
GO_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"DGE results input directory: {DGE_RESULTS_DIR}")
print(f"GO results output directory: {GO_OUTPUT_DIR}")
print(f"Organism: {ORGANISM}")
print(f"Gene sets: {GENE_SETS}")

# %%
# Check available gene sets for the organism
try:
    print(f"Checking available Enrichr libraries for organism: {ORGANISM}...")
    available_sets = gp.get_library_name(organism=ORGANISM)
    print(f"Found {len(available_sets)} libraries.")
    for gs in GENE_SETS:
        if gs not in available_sets:
            print(f"WARNING: Specified gene set '{gs}' not found in available Enrichr libraries for {ORGANISM}!")
except Exception as e:
    print(f"Could not retrieve available gene sets from Enrichr: {e}")

# %%
# Load the AnnData object to get the full list of genes (gene universe)
print(f"Loading AnnData from {ADATA_PATH} to get gene universe...")
try:
    adata = sc.read_h5ad(ADATA_PATH)
    print(adata)
    gene_universe = adata.var_names.tolist()
    print(f"Using {len(gene_universe)} genes as the background universe.")
    del adata
except FileNotFoundError:
    print(f"Error: AnnData file not found at {ADATA_PATH}. Cannot determine gene universe.")
    sys.exit(1)
except Exception as e:
    print(f"Error loading AnnData file: {e}")
    sys.exit(1)

# %%
# Find all significant DEG list files recursively
deg_list_files = list(pathlib.Path(DGE_RESULTS_DIR).rglob('*_significant.csv'))
print(f"Found {len(deg_list_files)} significant DEG list files to process.")

if not deg_list_files:
    print(f"Error: No significant DEG list files found matching '*_significant.csv' within {DGE_RESULTS_DIR}")
    sys.exit(1)

# Group DEG files by cell type and comparison name
grouped_deg_files = {}
for deg_file_path in deg_list_files:
    parts = deg_file_path.parts
    try:
        base_index = parts.index(DEGs_folder)
        if base_index + 1 < len(parts) and parts[base_index + 1] == 'biomarkers':
            biomarkers_index = base_index + 1
        else:
            continue

        if biomarkers_index + 1 >= len(parts) or parts[biomarkers_index + 1] != 'sig_deg_lists':
            continue

        context_start_index = biomarkers_index + 2
        sub_path_parts = parts[context_start_index:]

        if len(sub_path_parts) < 2:
            continue

        cell_type_name = sub_path_parts[0]
        filename = deg_file_path.name
        
        # Determine comparison_name and direction based on filename
        if '_up_significant.csv' in filename:
            direction = 'up'
            comparison_name = filename.replace('_up_significant.csv', '')
        elif '_down_significant.csv' in filename:
            direction = 'down'
            comparison_name = filename.replace('_down_significant.csv', '')
        else:
            direction = 'unknown'
            comparison_name = filename.replace('_significant.csv', '')

        key = (cell_type_name, comparison_name)
        if key not in grouped_deg_files:
            grouped_deg_files[key] = {}
        grouped_deg_files[key][direction] = deg_file_path
    except ValueError:
        continue
    except Exception as e:
        print(f"Error grouping file {deg_file_path}: {e}")
        continue

print(f"Found {len(grouped_deg_files)} unique comparisons to process.")

# %%
print(grouped_deg_files)

# %%
# Process each grouped DEG list - iterate through all cell type and comparison combinations
for (cell_type_name, comparison_name), files in grouped_deg_files.items():
    print(f"\n--- Processing comparison: {comparison_name} for cell type: {cell_type_name} ---")

    # Initialize variables to store gene lists and processing metadata
    gene_list_combined = []
    processing_type = ""
    output_file_suffix = ""
    
    # Extract UP and DOWN regulated gene files for this comparison
    up_file = files.get('up')
    down_file = files.get('down')

    # Initialize empty DataFrames for UP and DOWN gene data
    df_up = pd.DataFrame()
    df_down = pd.DataFrame()

    # Read and process UP-regulated genes if file exists
    if up_file:
        try:
            df_up = pd.read_csv(up_file)
            # Extract gene names from 'names' column and add to combined list
            if 'names' in df_up.columns:
                gene_list_combined.extend(df_up['names'].dropna().astype(str).tolist())
            print(f"  Read {len(df_up)} genes from {up_file.name} (UP)")
        except Exception as e:
            print(f"  Error reading UP file {up_file.name}: {e}")

    # Read and process DOWN-regulated genes if file exists
    if down_file:
        try:
            df_down = pd.read_csv(down_file)
            # Extract gene names from 'names' column and add to combined list
            if 'names' in df_down.columns:
                gene_list_combined.extend(df_down['names'].dropna().astype(str).tolist())
            print(f"  Read {len(df_down)} genes from {down_file.name} (DOWN)")
        except Exception as e:
            print(f"  Error reading DOWN file {down_file.name}: {e}")

    # Determine the type of analysis based on available files and set appropriate output suffix
    if up_file and down_file:
        processing_type = "Combined UP & DOWN"
        output_file_suffix = "_combined"
    elif up_file:
        processing_type = "UP only"
        output_file_suffix = "_up"
    elif down_file:
        processing_type = "DOWN only"
        output_file_suffix = "_down"
    else:
        print("  No UP or DOWN files found for this comparison. Skipping.")
        continue

    # Skip processing if no genes were found in the files
    if not gene_list_combined:
        print(f"  Combined gene list for {comparison_name} is empty. Skipping enrichment analysis.")
        continue

    # Remove duplicate genes and report final count
    gene_list_combined = list(set(gene_list_combined)) # Remove duplicates
    print(f"  Total {len(gene_list_combined)} unique genes for {processing_type} analysis.")

    # --- Prepare Output Directory and File Paths ---
    # Sanitize cell type name for file system compatibility (replace '/' with '_')
    sanitized_cell_type = cell_type_name.replace('/', '_')
    # Create output directory specific to this cell type
    specific_go_output_dir = pathlib.Path(os.path.join(GO_OUTPUT_DIR, sanitized_cell_type))
    specific_go_output_dir.mkdir(parents=True, exist_ok=True)

    # Define output file names for Excel results and dotplot visualization
    output_file_base = f"{comparison_name}{output_file_suffix}"
    excel_output = os.path.join(specific_go_output_dir, f"{output_file_base}_go_enrichment.xlsx")
    dotplot_output_file = os.path.join(specific_go_output_dir, f"{output_file_base}_go_enrichment_dotplot.png")

    # Create descriptive title for plots
    plot_title = f"GO Enrichment: {sanitized_cell_type} Biomarkers\n{comparison_name} ({processing_type})"

    # --- Run GO Enrichment Analysis ---
    print(f"  Running GO enrichment for {len(gene_list_combined)} genes ({processing_type})...")
    try:
        # Perform enrichment analysis using gseapy.enrichr
        enr = gp.enrichr(gene_list=gene_list_combined,
                         gene_sets=GENE_SETS,
                         organism=ORGANISM,
                         background=gene_universe,
                         outdir=None,
                         cutoff=0.1,
                         verbose=False)
    except Exception as enrich_err:
         print(f"  Error during gseapy.enrichr execution for {comparison_name} ({processing_type}): {enrich_err}")
         continue

    # --- Process and Save Enrichment Results ---
    # Check if enrichment analysis returned valid results
    if enr and hasattr(enr, 'results') and isinstance(enr.results, pd.DataFrame) and not enr.results.empty:
        print("  Filtering and saving enrichment results...")
        # Filter results to only include significantly enriched terms (adjusted p-value < threshold)
        filtered_results = enr.results[enr.results['Adjusted P-value'] < GO_PADJ_THRESHOLD].copy()

        if not filtered_results.empty:
             print(f"  Found {len(filtered_results)} significant terms (Adj P < {GO_PADJ_THRESHOLD}).")
             
             # Calculate gene count for each enriched term if 'Genes' column exists
             if 'Genes' in filtered_results.columns:
                 # Count number of genes in each term by splitting the semicolon-separated gene list
                 filtered_results['Gene_Count'] = filtered_results['Genes'].apply(
                     lambda x: len(str(x).split(';')) if pd.notna(x) and str(x) else 0
                 )
                 filtered_results['Gene_Count'] = pd.to_numeric(filtered_results['Gene_Count'], errors='coerce').fillna(1.0)
             else:
                 print("  Warning: 'Genes' column not found. Using default size.")
                 filtered_results['Gene_Count'] = 1.0

             # Save filtered results to Excel file
             try:
                 filtered_results.to_excel(excel_output, index=False)
                 print(f"  Saved significant results to {excel_output}")
             except Exception as excel_err:
                 print(f"  Error saving Excel file: {excel_err}")

             # --- Generate Visualization ---
             try:
                 # Select top terms for visualization (limited by N_TOP_TERMS_PLOT * number of gene sets)
                 plot_data = filtered_results.sort_values("Adjusted P-value").head(N_TOP_TERMS_PLOT * len(GENE_SETS)).copy()

                 if not plot_data.empty:
                     # Convert numeric columns to proper data types for plotting
                     for col in ['P-value', 'Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Gene_Count']:
                         if col in plot_data.columns:
                             plot_data[col] = pd.to_numeric(plot_data[col], errors='coerce')
                     
                     # Remove rows with missing adjusted p-values
                     plot_data.dropna(subset=['Adjusted P-value'], inplace=True)

                     if not plot_data.empty:
                         print("  Data types before plotting:")
                         print(plot_data.dtypes)
                         print(f"  Generating visualizations for {len(plot_data)} enriched terms...")
                         
                         try:
                             print(f"  Creating dotplot visualization: {dotplot_output_file}")
                             
                             # Select top 15 most significant terms for the plot
                             top_overall = plot_data.sort_values('Adjusted P-value').head(15)
                             
                             if not top_overall.empty:
                                 # Clean up term names by removing GO IDs and replacing underscores
                                 top_overall['Clean_Term'] = top_overall['Term'].str.replace(r' \(GO:[0-9]+\)', '', regex=True)
                                 top_overall['Clean_Term'] = top_overall['Clean_Term'].str.replace('_', ' ')
                                 
                                 # Calculate gene ratio (normalized by maximum gene count)
                                 top_overall['GeneRatio'] = top_overall['Gene_Count'] / top_overall['Gene_Count'].max()
                                 
                                 # Sort by adjusted p-value for better visualization
                                 top_overall = top_overall.sort_values('Adjusted P-value', ascending=False)
                                 
                                 # Create figure with dynamic height based on number of terms
                                 plt.figure(figsize=(10, max(6, len(top_overall) * 0.4)))
                                 
                                 # Create scatter plot with constant dot size
                                 scatter = plt.scatter(
                                     top_overall['GeneRatio'],
                                     top_overall['Clean_Term'],
                                     s=100,  # Constant dot size
                                     c=-np.log10(top_overall['Adjusted P-value']),  # Color by significance
                                     cmap='Reds',
                                     alpha=0.8
                                 )
                                 
                                 # Add gene count text next to each dot
                                 for idx, row in top_overall.iterrows():
                                     plt.text(row['GeneRatio'] + 0.02, row['Clean_Term'], 
                                             f"({int(row['Gene_Count'])})", 
                                             va='center', ha='left', fontsize=9)
                                 
                                 # Set plot labels and title
                                 plt.xlabel('GeneRatio')
                                 plt.ylabel('')
                                 plt.title(f"{plot_title}\nTop GO Terms")
                                 plt.grid(True, alpha=0.3, axis='x')
                                 
                                 cbar = plt.colorbar(scatter)
                                 cbar.set_label('-log10(Adj. P-value)')
                                 
                                 valid_counts = top_overall['Gene_Count'].dropna().unique()
                                 if len(valid_counts) > 0:
                                     size_steps = sorted(list(valid_counts))
                                     if len(size_steps) > 3:
                                         size_steps = [min(size_steps),
                                                       np.percentile(size_steps, 50),
                                                       max(size_steps)]
                                 else:
                                     size_steps = []
                                 # Remove the legend for gene count (no plt.scatter for legend, no plt.legend)
                                 # for count in size_steps:
                                 #     plt.scatter([], [], s=np.minimum(count * 5, 400), c='black',
                                 #               label=f'{int(count)}')
                                 # plt.legend(title="Gene Count", loc='lower right', frameon=True)
                                 
                                 plt.tight_layout()
                                 plt.savefig(dotplot_output_file, dpi=150, bbox_inches='tight')
                                 plt.close()
                                 
                                 print(f"  Dotplot saved to {dotplot_output_file}")
                         except Exception as dotplot_err:
                             print(f"  Error generating dotplot: {dotplot_err}")
                     else:
                         print("  No significant terms left after filtering for visualization.")
                 else:
                     print("  No significant terms left after filtering for visualization.")
             except Exception as viz_err:
                 print(f"  Warning: Could not generate visualizations. Error: {viz_err}")
        else:
            print(f"  No significantly enriched terms found after filtering (Adj P < {GO_PADJ_THRESHOLD}). No Excel file or visualizations saved.")
    else:
         print("  Enrichment analysis did not return results or results table was empty.")

print("\nGO enrichment analysis for DEG lists complete.")
print(f"Results saved in subdirectories under: {GO_OUTPUT_DIR}")


