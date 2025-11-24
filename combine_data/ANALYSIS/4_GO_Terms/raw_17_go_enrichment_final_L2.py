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
# REMOVE_DOUBLETS = True
REMOVE_DOUBLETS = False

FIX_TRESHOLD = True
# FIX_TRESHOLD = False

if FIX_TRESHOLD:
    BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")
else:
    if REMOVE_DOUBLETS:
        BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw_percentile_threshold", "doublets_removed")
    else:
        BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw_percentile_threshold")

INPUT_DIR = BASE_RESULTS_DIR
ADATA_PATH = os.path.join(INPUT_DIR, "annotation_final.h5ad")

DGE_RESULTS_DIR = os.path.join(INPUT_DIR,"DEGs_cell_type_L2")

# CUSTOM_ANALYSIS =  None
CUSTOM_ANALYSIS =  "FC_0_25"

DEG_BY = "cell_type_L2"
DEGs_folder = "DEGs_cell_type_L2FC_0_25"

if CUSTOM_ANALYSIS is not None:
    DGE_RESULTS_DIR = DGE_RESULTS_DIR + CUSTOM_ANALYSIS

DGE_RESULTS_DIR = os.path.join(DGE_RESULTS_DIR,"biomarkers")
GO_OUTPUT_DIR = os.path.join(DGE_RESULTS_DIR, 'go_enrichment')
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
    # You can optionally print all available_sets if needed for debugging:
    # print(available_sets) 
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

# %%
# Process each DEG list file
for deg_file_path in deg_list_files:
    print(f"\n--- Processing {deg_file_path.name} ---")
    print(f"Full path: {deg_file_path}")

    try:
        # --- Parse file path to get context ---
        parts = deg_file_path.parts
        # Find the index of the 'biomarkers' directory within the path
        try:
            # Find the index of DEG_BY first
            base_index = parts.index(DEGs_folder)
            # Check if 'biomarkers' is the next part
            if base_index + 1 < len(parts) and parts[base_index + 1] == 'biomarkers':
                biomarkers_index = base_index + 1
            else:
                print(f"Warning: Could not find 'biomarkers' directory after f'{DEGs_folder}' in path {deg_file_path}. Skipping.")
                continue
        except ValueError:
            print(f"Warning: Could not find f'{DEGs_folder}' in path {deg_file_path}. Skipping.")
            continue

        # Check for 'sig_deg_lists' after 'biomarkers'
        if biomarkers_index + 1 >= len(parts) or parts[biomarkers_index + 1] != 'sig_deg_lists':
            print(f"Warning: Expected 'sig_deg_lists' directory after 'biomarkers' not found in path {deg_file_path}. Skipping.")
            continue

        # Extract context based on the biomarker structure
        # Example path: .../*/biomarkers/sig_deg_lists/Astro/Cluster_Astro_vs_Rest_up_significant.csv
        context_start_index = biomarkers_index + 2 # Start after 'sig_deg_lists'
        sub_path_parts = parts[context_start_index:] # Parts after sig_deg_lists

        # Expected structure: <cell_type>/<filename> -> 2 parts minimum
        if len(sub_path_parts) < 2:
            print(f"Warning: Unexpected path structure after 'sig_deg_lists' for {deg_file_path}. Parts: {sub_path_parts}. Skipping.")
            continue

        analysis_type = 'biomarkers' # This is specifically for biomarkers
        cell_type_name = sub_path_parts[0] # The directory name is the cell type
        comparison_details = 'Marker Genes' # Implicitly marker genes

        filename = deg_file_path.name
        direction = 'up' if '_up_' in filename else 'down' if '_down_' in filename else 'unknown'
        # Extract comparison name from filename (e.g., Cluster_Astro_vs_Rest)
        comparison_name = filename.replace(f'_{direction}_significant.csv', '')

        # Basic validation for extracted cell_type_name
        if not cell_type_name:
            print(f"Warning: Could not extract cell type name for {deg_file_path}. Skipping.")
            continue

        print(f"  Context: Type={analysis_type}, Details={comparison_details or 'N/A'}, CellType={cell_type_name}, Comparison={comparison_name}, Direction={direction}")

        # --- Read DEG list ---
        df = pd.read_csv(deg_file_path)
        print(f"  Read {len(df)} genes from {deg_file_path.name}")

        if 'names' not in df.columns:
            print(f"  Error: 'names' column not found in {deg_file_path.name}. Skipping.")
            continue

        gene_list = df['names'].dropna().astype(str).tolist()

        if not gene_list:
            print("  Gene list is empty. Skipping enrichment analysis.")
            continue

        # --- Prepare Output ---
        # Create specific output directory for this result
        # Structure: GO_OUTPUT_DIR / <sanitized_cell_type> / ...
        sanitized_cell_type = cell_type_name.replace('/', '_') # Sanitize cell type name for path
        specific_go_output_dir = pathlib.Path(os.path.join(GO_OUTPUT_DIR, sanitized_cell_type))

        specific_go_output_dir.mkdir(parents=True, exist_ok=True)

        # Define output file names (using sanitized cell_type_name)
        output_file_base = f"{comparison_name}_{direction}" # Keep original comparison name here
        excel_output = os.path.join(specific_go_output_dir, f"{output_file_base}_go_enrichment.xlsx")
        
        # Define visualization output files - only keep bubble plots and dotplots
        dotplot_output_file = os.path.join(specific_go_output_dir, f"{output_file_base}_go_enrichment_dotplot.png")
        
        # Construct a more informative title for biomarkers
        plot_title = f"GO Enrichment: {sanitized_cell_type} Biomarkers\n{comparison_name} ({direction.upper()})"

        # --- Run Enrichment ---
        print(f"  Running GO enrichment for {len(gene_list)} genes...")
        try:
            enr = gp.enrichr(gene_list=gene_list,
                             gene_sets=GENE_SETS,
                             organism=ORGANISM,
                             background=gene_universe, # Provide the background gene list
                             outdir=None, # Prevent gseapy from creating its own directory structure
                             cutoff=0.1, # Use a slightly less stringent cutoff for enrichr, filter later
                             verbose=False)
        except Exception as enrich_err:
             print(f"  Error during gseapy.enrichr execution: {enrich_err}")
             continue # Skip to next file if enrichment fails

        # --- Process and Save Results ---
        if enr and hasattr(enr, 'results') and isinstance(enr.results, pd.DataFrame) and not enr.results.empty:
            print("Filtering and saving enrichment results...")
            # Filter results based on adjusted p-value
            filtered_results = enr.results[enr.results['Adjusted P-value'] < GO_PADJ_THRESHOLD].copy() # Use .copy()

            if not filtered_results.empty:
                 print(f"  Found {len(filtered_results)} significant terms (Adj P < {GO_PADJ_THRESHOLD}).")
                 # Calculate gene count from 'Genes' column for size mapping
                 if 'Genes' in filtered_results.columns:
                     # Ensure robust splitting even if genes are missing or format varies
                     filtered_results['Gene_Count'] = filtered_results['Genes'].apply(
                         lambda x: len(str(x).split(';')) if pd.notna(x) and str(x) else 0
                     )
                     # Immediately convert to float to ensure numeric type from the start
                     filtered_results['Gene_Count'] = pd.to_numeric(filtered_results['Gene_Count'], errors='coerce').fillna(1.0)
                 else:
                     print("  Warning: 'Genes' column not found. Using default size.")
                     filtered_results['Gene_Count'] = 1.0  # Assign a default size as float

                 # Save filtered results to Excel
                 try:
                     filtered_results.to_excel(excel_output, index=False)
                     print(f"  Saved significant results to {excel_output}")
                 except Exception as excel_err:
                     print(f"  Error saving Excel file: {excel_err}")


                 # Generate visualizations for top terms
                 try:
                     # Sort by Adj P-value for plotting top terms
                     plot_data = filtered_results.sort_values("Adjusted P-value").head(N_TOP_TERMS_PLOT * len(GENE_SETS)).copy() # Use .copy()

                     if not plot_data.empty:
                         # Ensure relevant columns are numeric before plotting
                         for col in ['P-value', 'Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Gene_Count']:
                             if col in plot_data.columns:
                                 plot_data[col] = pd.to_numeric(plot_data[col], errors='coerce')
                         
                         # Drop rows where conversion failed (NaN values)
                         plot_data.dropna(subset=['Adjusted P-value'], inplace=True)

                         if not plot_data.empty:
                             print("  Data types before plotting:")
                             print(plot_data.dtypes)
                             print(f"  Generating visualizations for {len(plot_data)} enriched terms...")
                             
                             # 1. Create a dotplot (similar to second example)
                             try:
                                 print(f"  Creating dotplot visualization: {dotplot_output_file}")
                                 
                                 # Get top terms across all categories (most significant)
                                 top_overall = plot_data.sort_values('Adjusted P-value').head(15)
                                 
                                 if not top_overall.empty:
                                     # Process the terms
                                     top_overall['Clean_Term'] = top_overall['Term'].str.replace(r' \(GO:[0-9]+\)', '', regex=True)
                                     top_overall['Clean_Term'] = top_overall['Clean_Term'].str.replace('_', ' ')
                                     
                                     # Calculate gene ratio for x-axis
                                     top_overall['GeneRatio'] = top_overall['Gene_Count'] / top_overall['Gene_Count'].max()
                                     
                                     # Sort by p-value for intuitive display
                                     top_overall = top_overall.sort_values('Adjusted P-value', ascending=False)
                                     
                                     # Plot
                                     plt.figure(figsize=(10, max(6, len(top_overall) * 0.4))) # Dynamic height
                                     
                                     # Create a scatter plot with constant dot size
                                     scatter = plt.scatter(
                                         top_overall['GeneRatio'],
                                         top_overall['Clean_Term'],
                                         s=100,  # Constant dot size
                                         c=-np.log10(top_overall['Adjusted P-value']),  # Color by p-value
                                         cmap='Reds',  # Red colormap
                                         alpha=0.8
                                     )
                                     
                                     # Add gene count text next to each dot
                                     for idx, row in top_overall.iterrows():
                                         plt.text(row['GeneRatio'] + 0.02, row['Clean_Term'], 
                                                 f"({int(row['Gene_Count'])})", 
                                                 va='center', ha='left', fontsize=9)
                                     
                                     # Customize
                                     plt.xlabel('GeneRatio')
                                     plt.ylabel('')
                                     plt.title(f"{plot_title}\nTop GO Terms")
                                     plt.grid(True, alpha=0.3, axis='x')
                                     
                                     # Add colorbar
                                     cbar = plt.colorbar(scatter)
                                     cbar.set_label('-log10(Adj. P-value)') # Updated label
                                     
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
    except FileNotFoundError:
        print(f"Error: DEG list file not found (should not happen with glob): {deg_file_path}")
    except pd.errors.EmptyDataError:
        print(f"Error: DEG list file is empty: {deg_file_path.name}")
    except KeyError as e:
        print(f"Error: Missing expected column in {deg_file_path.name}: {e}. Ensure 'names' column exists.")
    except Exception as e:
        print(f"An unexpected error occurred processing {deg_file_path.name}: {e}")
        import traceback
        traceback.print_exc() # Print detailed traceback for debugging

print("\nGO enrichment analysis for DEG lists complete.")
print(f"Results saved in subdirectories under: {GO_OUTPUT_DIR}")

# %%



