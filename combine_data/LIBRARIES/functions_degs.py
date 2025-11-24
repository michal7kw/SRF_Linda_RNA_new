import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re # Added for sanitizing filenames
import seaborn as sns
import pathlib # For creating paths robustly
try:
    import openpyxl # For writing .xlsx files
except ImportError:
    openpyxl = None # Flag that it's not available

# --- DGE Significance Thresholds ---
LOGFC_THRESHOLD = 0.25

# --- Data Quality & Filtering Thresholds ---
MIN_CELLS_FOR_DGE_AND_PLOTTING = 50  # Minimum cells per group for a DGE comparison to run.
MIN_PERCENTAGE_OF_CELLS = 0.10      # A gene must be in at least this percentage of cells in a subgroup.
MIN_GENE_EXPRESSION_CELLS_ABSOLUTE = 3 # Fallback absolute minimum number of cells for gene expression.
MIN_CELLS_GLOBAL_FILTER = 3         # A gene must be in at least this many cells across the entire dataset.

def sanitize_filename(name):
    """Sanitizes a string to be used as a valid filename."""
    # Convert to string and replace spaces and slashes with underscores
    s = str(name)
    s = re.sub(r'[\\/\s]+', '_', s)
    # Remove any characters that are not alphanumeric, underscore, or hyphen
    s = re.sub(r'[^a-zA-Z0-9_-]', '', s)
    # Replace multiple underscores with a single one
    s = re.sub(r'_+', '_', s)
    # Remove leading/trailing underscores
    s = s.strip('_')
    # If the name becomes empty after sanitization, return a default
    if not s:
        return "unnamed_group"
    return s


def save_significant_degs(result_df, base_output_dir, analysis_type, comparison_details, group_name,
                          comparison_name, padj_threshold=0.05, logfc_threshold=LOGFC_THRESHOLD):
    """
    Filters DGE results for significant up/down regulated genes and saves them to text files.

    Args:
        result_df (pd.DataFrame): DataFrame with DGE results (needs 'names', 'logfoldchanges', 'pvals_adj').
        base_output_dir (str): The main output directory for the specific DGE run (e.g., overall_dge_dir).
        analysis_type (str): Type of analysis (e.g., 'overall', 'genotype_specific', 'genotype_comparison').
        comparison_details (str): Subdirectory structure specific to the analysis (e.g., '', 'Emx1', 'Control').
        group_name (str): Name of the cell type or group being analyzed.
        comparison_name (str): Name for the file prefix (e.g., 'Mutant_vs_Control', 'Nestin_vs_Emx1').
        padj_threshold (float): Adjusted p-value threshold for significance.
        logfc_threshold (float): Absolute log2 fold change threshold for significance.
    """
    if result_df is None or result_df.empty:
        print(f"    Skipping DEG list saving for {analysis_type}/{comparison_details}/{group_name}/{comparison_name} (No results).")
        return

    sanitized_group_name = sanitize_filename(group_name)

    # Define output directory for the lists
    list_output_dir = pathlib.Path(base_output_dir) / 'sig_deg_lists' # Shortened
    if comparison_details:
        list_output_dir = list_output_dir / str(comparison_details) # Ensure comparison_details is also string
    list_output_dir = list_output_dir / sanitized_group_name

    print(f"    DEBUG: Attempting to create list_output_dir: '{list_output_dir}'")
    try:
        list_output_dir.mkdir(parents=True, exist_ok=True)
        if list_output_dir.is_dir():
            print(f"    DEBUG: Successfully created/confirmed directory: '{list_output_dir}'")
            # Try to create a temporary file to check writability
            temp_file_path = list_output_dir / "temp_writetest.txt"
            try:
                with open(temp_file_path, "w") as f_test:
                    f_test.write("test")
                os.remove(temp_file_path) # Clean up the temp file
                print(f"    DEBUG: Directory '{list_output_dir}' is writable.")
            except Exception as e_write:
                print(f"    DEBUG: ERROR - Directory '{list_output_dir}' IS NOT WRITABLE or temp file creation failed: {e_write}")
                # If directory is not writable, we should not proceed with saving files.
                return # Stop further execution in this function call
        else:
            print(f"    DEBUG: ERROR - Failed to create directory: '{list_output_dir}' (is_dir() is false after mkdir).")
            # This is a critical failure, the directory doesn't exist after mkdir.
            return # Stop further execution in this function call
    except Exception as e_mkdir:
        print(f"    DEBUG: ERROR - Exception during mkdir for '{list_output_dir}': {e_mkdir}")
        # This means mkdir itself failed.
        return # Stop further execution in this function call

    # Filter for significant genes
    sig_mask = result_df['pvals_adj'] < padj_threshold
    up_mask = sig_mask & (result_df['logfoldchanges'] > logfc_threshold)
    down_mask = sig_mask & (result_df['logfoldchanges'] < -logfc_threshold)

    up_df = result_df.loc[up_mask].copy()
    down_df = result_df.loc[down_mask].copy()

    # Save UP list as CSV
    up_filename = list_output_dir / f"{comparison_name}_up_significant.csv"
    if not up_df.empty:
        up_df.to_csv(up_filename, index=False)
        print(f"    Saved {len(up_df)} significant UP genes with stats to {up_filename}")
    else:
        # Create empty file with header to indicate no significant genes found
        pd.DataFrame(columns=result_df.columns).to_csv(up_filename, index=False)
        print(f"    No significant UP genes found for {comparison_name}. Empty CSV file with header created: {up_filename}")


    # Save DOWN list as CSV
    down_filename = list_output_dir / f"{comparison_name}_down_significant.csv"
    if not down_df.empty:
        down_df.to_csv(down_filename, index=False)
        print(f"    Saved {len(down_df)} significant DOWN genes with stats to {down_filename}")
    else:
        # Create empty file with header
        pd.DataFrame(columns=result_df.columns).to_csv(down_filename, index=False)
        print(f"    No significant DOWN genes found for {comparison_name}. Empty CSV file with header created: {down_filename}")

# Define a function for creating volcano plots
def plot_volcano(result_df, group_name, n1, n2, output_dir, filename_prefix, group1_name="Mutant", group2_name="Control", title_prefix="",
                 padj_threshold=0.05, logfc_threshold=LOGFC_THRESHOLD, num_genes_to_label=5,
                 logfc_filter_tolerance=0.05, up_color='red', down_color='blue'):
    """
    Generates and saves a volcano plot for differential gene expression results.

    Args:
        result_df (pd.DataFrame): DataFrame with DGE results (needs 'names', 'logfoldchanges', 'pvals_adj').
        group_name (str): Name of the cell type or group being plotted.
        n1 (int): Number of cells in group 1 (e.g., Mutant).
        n2 (int): Number of cells in group 2 (e.g., Control).
        output_dir (str): Directory to save the plot.
        filename_prefix (str): Prefix for the output plot filename.
        group1_name (str): Name of the group corresponding to positive logFC.
        group2_name (str): Name of the group corresponding to negative logFC (reference).
        title_prefix (str): Prefix for the plot title (e.g., "Genotype: Emx1 - ").
        padj_threshold (float): Adjusted p-value threshold for significance.
        logfc_threshold (float): Absolute log2 fold change threshold for significance.
        num_genes_to_label (int): Number of top up/down regulated genes to label.
        logfc_filter_tolerance (float): Tolerance for filtering near-zero log fold changes.
        up_color (str): Color for significantly upregulated genes.
        down_color (str): Color for significantly downregulated genes.
    """
    # Add a check for minimum cell count in either group
    min_cells_for_plot = MIN_CELLS_FOR_DGE_AND_PLOTTING
    if n1 < min_cells_for_plot or n2 < min_cells_for_plot:
        print(f"  Skipping volcano plot for {title_prefix}{group_name}: Not enough cells in one or both groups (< {min_cells_for_plot}). {group2_name}={n2}, {group1_name}={n1}.")
        return

    if result_df is None:
        print(f"  Skipping volcano plot for {title_prefix}{group_name} (no DGE results, Counts: {group2_name}={n2}, {group1_name}={n1}).")
        return

    if result_df.empty:
        print(f"  Skipping volcano plot for {title_prefix}{group_name} (empty DGE results DataFrame, Counts: {group2_name}={n2}, {group1_name}={n1}).")
        return

    sanitized_group_name = sanitize_filename(group_name)
    output_filename = f"{filename_prefix}_{sanitized_group_name}_volcano.png" # Use group_name in filename
    # Ensure output_dir exists before joining path
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_filename)

    print(f"  Generating volcano plot for: {title_prefix}{group_name} (Counts: {group2_name}={n2}, {group1_name}={n1}) -> {output_filename}")

    # Prepare data for plotting
    df = result_df.copy()
    # Filter out genes with logfoldchange very close to 0
    df = df[df['logfoldchanges'].abs() > logfc_filter_tolerance].copy()

    # If filtering results in an empty dataframe, skip
    if df.empty:
        print(f"    Skipping volcano plot for {title_prefix}{group_name} after filtering logFC near 0 (Counts: {group2_name}={n2}, {group1_name}={n1}).")
        return

    df['-log10p_adj'] = -np.log10(df['pvals_adj'])

    # Handle infinite values resulting from p_adj = 0
    max_finite_logp = df.loc[np.isfinite(df['-log10p_adj']), '-log10p_adj'].max()
    if pd.isna(max_finite_logp) or max_finite_logp == 0: # Handle case where all p-values might be 1 or NaN
        max_finite_logp = 300 # Default max if no finite values or max is 0
    inf_replacement = max_finite_logp * 1.1 if max_finite_logp > 0 else 10 # Ensure non-zero replacement
    df['-log10p_adj'] = df['-log10p_adj'].replace([np.inf, -np.inf], inf_replacement)
    df['-log10p_adj'] = df['-log10p_adj'].replace(np.nan, 0) # Replace NaN with 0

    # Replace potential NaN logfoldchanges with neutral values for plotting
    df['logfoldchanges'] = df['logfoldchanges'].fillna(0)


    # Determine significance and color
    df['color'] = 'grey' # Default: not significant
    sig_up_mask = (df['pvals_adj'] < padj_threshold) & (df['logfoldchanges'] > logfc_threshold)
    sig_down_mask = (df['pvals_adj'] < padj_threshold) & (df['logfoldchanges'] < -logfc_threshold)
    df.loc[sig_up_mask, 'color'] = up_color # Significant Up
    df.loc[sig_down_mask, 'color'] = down_color # Significant Down

    # Count significant DEGs
    n_degs_up = sig_up_mask.sum()
    n_degs_down = sig_down_mask.sum()
    n_degs_total = n_degs_up + n_degs_down

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 7))
    # Use rasterized=True for large number of points to keep vector file size manageable
    ax.scatter(df['logfoldchanges'], df['-log10p_adj'], c=df['color'], alpha=0.6, s=10, rasterized=True)

    # Add threshold lines
    ax.axhline(-np.log10(padj_threshold), color='grey', linestyle='--', lw=1)
    ax.axvline(logfc_threshold, color='grey', linestyle='--', lw=1)
    ax.axvline(-logfc_threshold, color='grey', linestyle='--', lw=1)

    # Add labels for top genes
    # Sort by significance first, then by fold change magnitude
    # Need boolean masks on the original df index before filtering for labeling
    df_label_candidates = result_df.copy() # Use original for labeling candidates
    df_label_candidates['-log10p_adj'] = -np.log10(df_label_candidates['pvals_adj'])
    max_finite_logp_label = df_label_candidates.loc[np.isfinite(df_label_candidates['-log10p_adj']), '-log10p_adj'].max()
    if pd.isna(max_finite_logp_label) or max_finite_logp_label == 0: 
        max_finite_logp_label = 300
    inf_replacement_label = max_finite_logp_label * 1.1 if max_finite_logp_label > 0 else 10
    df_label_candidates['-log10p_adj'] = df_label_candidates['-log10p_adj'].replace([np.inf, -np.inf], inf_replacement_label)
    df_label_candidates['-log10p_adj'] = df_label_candidates['-log10p_adj'].replace(np.nan, 0) # Replace NaN with 0
    df_label_candidates['logfoldchanges'] = df_label_candidates['logfoldchanges'].fillna(0) # Fill NaN logFC

    orig_sig_up_mask = (df_label_candidates['pvals_adj'] < padj_threshold) & (df_label_candidates['logfoldchanges'] > logfc_threshold)
    orig_sig_down_mask = (df_label_candidates['pvals_adj'] < padj_threshold) & (df_label_candidates['logfoldchanges'] < -logfc_threshold)
    df_label_candidates['abs_logfc'] = df_label_candidates['logfoldchanges'].abs()


    # Sort candidates for labeling based on p-value and fold change
    df_label_sorted = df_label_candidates[(orig_sig_up_mask | orig_sig_down_mask)].sort_values(
        by=['-log10p_adj', 'abs_logfc'], ascending=[False, False]
    )

    # Select top N up and down regulated genes from the sorted significant list
    top_up_labels = df_label_sorted[df_label_sorted['logfoldchanges'] > 0].head(num_genes_to_label)
    top_down_labels = df_label_sorted[df_label_sorted['logfoldchanges'] < 0].head(num_genes_to_label)
    genes_to_label = pd.concat([top_up_labels, top_down_labels])

    # Ensure logfoldchange is not near zero for labeled genes if they survived the main filter
    # Re-check against original result_df as df_label_candidates might have imputed NaNs
    genes_to_label = genes_to_label[result_df.loc[genes_to_label.index, 'logfoldchanges'].abs() > logfc_filter_tolerance]


    for idx, row in genes_to_label.iterrows():
        # Use the potentially updated -log10p_adj for plotting position
        plot_y = row['-log10p_adj']
        # Ensure we don't try to plot text for NaNs that might have slipped through
        if pd.notna(row['logfoldchanges']) and pd.notna(plot_y) and pd.notna(row['names']):
             ax.text(row['logfoldchanges'], plot_y, row['names'], fontsize=8)

    # Set labels and title
    ax.set_xlabel(f"Log2 Fold Change ({group1_name} vs {group2_name})")
    ax.set_ylabel("-log10 (Adjusted P-value)")
    ax.set_title(f"Volcano Plot: {title_prefix}{group_name}\n" # Use group_name in title
                 f"({group2_name}: {n2} cells, {group1_name}: {n1} cells)\n"
                 f"DEGs: {n_degs_total} (Up in {group1_name}: {n_degs_up}, Down in {group1_name}: {n_degs_down})")
    ax.grid(False)
    sns.despine(ax=ax)

    plt.tight_layout()
    # Save the figure
    try:
        fig.savefig(output_path, bbox_inches='tight', dpi=300)
    except Exception as e:
        print(f"    Error saving volcano plot for {group_name}: {e}")
    # plt.show() # Removed show
    plt.close(fig) # Close figure to free memory
    print("---")

def run_overall_dge(adata, grouping_key, dge_output_dir, plot_output_dir, layer=None):
    """
    Performs DGE analysis (Mutant vs Control) across groups defined by grouping_key for both genomes.

    Args:
        adata (ad.AnnData): Annotated data object.
        grouping_key (str): The column name in adata.obs to use for grouping (e.g., 'cell_type', 'leiden_0.8_alt').
        dge_output_dir (str): Base directory to save DGE result tables.
        plot_output_dir (str): Base directory to save volcano plots.
        layer (str, optional): Layer in adata to use for DGE. If None, uses adata.X. Defaults to None.
                             WARNING: If layer is None, adata.X is used. If adata.X was generated using
                             normalize_total, results may be unreliable for comparing groups with
                             large differences in raw counts/library size.

    Returns:
        dict: Dictionary containing DGE results per group.
              Key: group_value (e.g., specific cell type or cluster ID)
              Value: tuple (results_df, n_cells_control, n_cells_mutant)
    """
    # Define and create specific output subdirectories
    overall_dge_dir = os.path.join(dge_output_dir, 'both_geno_cond_comp') # Shortened
    overall_plot_dir = os.path.join(plot_output_dir, 'both_geno_cond_comp') # Shortened
    os.makedirs(overall_dge_dir, exist_ok=True)
    os.makedirs(overall_plot_dir, exist_ok=True)

    print("\nStarting both genomes Differential Gene Expression analysis (Mutant vs Control)...")

    # Ensure 'condition' is categorical for rank_genes_groups
    if not pd.api.types.is_categorical_dtype(adata.obs['condition']): # type: ignore
        adata.obs['condition'] = adata.obs['condition'].astype('category')

    # Set 'Control' as the reference category if not already
    if 'Control' in adata.obs['condition'].cat.categories and 'Mutant' in adata.obs['condition'].cat.categories:
        adata.obs['condition'] = adata.obs['condition'].cat.reorder_categories(['Control', 'Mutant'], ordered=True)
    else:
        print("Warning: Could not set 'Control' as reference. Categories might be different.")

    dge_results = {} # Store tuple: (results_df, n_cells_control, n_cells_mutant)
    # Check if grouping_key exists
    if grouping_key not in adata.obs.columns:
        raise ValueError(f"Grouping key '{grouping_key}' not found in adata.obs. Available columns: {list(adata.obs.columns)}")
    groups = adata.obs[grouping_key].unique()
    skipped_dge = []

    for group_value in groups:
        # min_cells_for_gene = 3 # Minimum number of cells a gene must be expressed in to be kept for DGE - Replaced by percentage
        print(f"  Running DGE for group: {group_value} (using key '{grouping_key}')")

        # Subset AnnData to the current group
        adata_subset = adata[adata.obs[grouping_key] == group_value].copy()

        # Get cell counts per condition for this subset
        condition_counts = adata_subset.obs['condition'].value_counts()
        n_ctrl = condition_counts.get('Control', 0)
        n_mut = condition_counts.get('Mutant', 0)

        # Skip if too few cells in either condition for comparison (minimum 50 for DGE calculation)
        if n_ctrl < MIN_CELLS_FOR_DGE_AND_PLOTTING or n_mut < MIN_CELLS_FOR_DGE_AND_PLOTTING:
            print(f"    Skipping DGE for group {group_value}: Not enough cells (< {MIN_CELLS_FOR_DGE_AND_PLOTTING}). Control: {n_ctrl}, Mutant: {n_mut}")
            dge_results[group_value] = (None, n_ctrl, n_mut)
            skipped_dge.append({grouping_key: group_value, 'reason': f'Low cell count (Ctrl={n_ctrl}, Mut={n_mut})'})
            continue

        # Filter genes expressed in too few cells within this subset (require expression in >= 10% of cells)
        print(f"    Initial genes: {adata_subset.n_vars}")
        
        calculated_min_cells_by_percentage = int(np.ceil(adata_subset.n_obs * MIN_PERCENTAGE_OF_CELLS))
        effective_min_cells = max(MIN_GENE_EXPRESSION_CELLS_ABSOLUTE, calculated_min_cells_by_percentage)
        
        print(f"    Filtering genes: Requiring expression in >= {effective_min_cells} cells (max of {MIN_GENE_EXPRESSION_CELLS_ABSOLUTE} absolute and {MIN_PERCENTAGE_OF_CELLS*100}% of {adata_subset.n_obs})...")
        sc.pp.filter_genes(adata_subset, min_cells=effective_min_cells)
        print(f"    Genes after filtering: {adata_subset.n_vars}")
        if adata_subset.n_vars == 0:
            print(f"    Skipping group {group_value}: No genes remaining after filtering.")
            dge_results[group_value] = (None, n_ctrl, n_mut)
            skipped_dge.append({grouping_key: group_value, 'reason': f'No genes left after filtering (min_cells={effective_min_cells}, 10%)'})
            continue

        # Run DGE comparing Mutant vs Control (Control is reference)
        # Uses specified layer or adata.X.
        if layer is None:
             print(f"    WARNING: Running DGE on adata.X for group {group_value}. Ensure this data representation is appropriate for DGE (see function docstring warning).")

        try:
            sc.tl.rank_genes_groups(
                adata_subset,
                groupby='condition',
                method='wilcoxon', # Changed from 't-test'
                layer=layer, # Use specified layer
                use_raw=False, # layer takes precedence over use_raw
                n_genes=adata_subset.n_vars, # Rank all genes to get full stats if needed
                pts=True # Calculate fraction of cells expressing the gene
            )

            # Store results in the dictionary
            # Extracting results for the 'Mutant' group compared to 'Control'
            result_df = sc.get.rank_genes_groups_df(adata_subset, group='Mutant')

            # --- Add mean expression columns ---
            if not result_df.empty:
                genes_in_result = result_df['names'].tolist()
                # Ensure genes exist in the subset's var_names
                valid_genes = [g for g in genes_in_result if g in adata_subset.var_names]
                if not valid_genes:
                    print("      Warning: No genes from DGE results found in the subset anndata. Skipping mean expression calculation.")
                else:
                    print(f"      Calculating mean expression for {len(valid_genes)} genes...")
                    # Get indices for control and mutant cells in the current subset
                    control_indices = adata_subset.obs.index[adata_subset.obs['condition'] == 'Control']
                    mutant_indices = adata_subset.obs.index[adata_subset.obs['condition'] == 'Mutant']

                    # Select the data layer (use the same layer as DGE)
                    data_source = adata_subset[:, valid_genes]
                    if layer and layer in data_source.layers:
                        data_matrix = data_source.layers[layer]
                        layer_name = layer
                    else:
                        data_matrix = data_source.X
                        layer_name = 'X'
                    print(f"      Using data from layer: {layer_name}")

                    # Ensure data_matrix is dense for indexing if it's sparse
                    if hasattr(data_matrix, "toarray"):
                         data_matrix = data_matrix.toarray()

                    # Create a DataFrame for easier indexing
                    expr_df = pd.DataFrame(data_matrix, index=adata_subset.obs_names, columns=valid_genes)

                    # Calculate mean expression for each group
                    mean_expr_control = expr_df.loc[control_indices].mean(axis=0)
                    mean_expr_mutant = expr_df.loc[mutant_indices].mean(axis=0)

                    # Add to result_df, aligning by gene name ('names' column)
                    result_df[f'mean_expr_Control ({layer_name})'] = result_df['names'].map(mean_expr_control)
                    result_df[f'mean_expr_Mutant ({layer_name})'] = result_df['names'].map(mean_expr_mutant)
                    print("Added mean expression columns.")
            # --- End mean expression ---

            dge_results[group_value] = (result_df, n_ctrl, n_mut) # Store potentially modified df

            # Save DGE results table (CSV) - Use specific subdir
            sanitized_group_value = sanitize_filename(group_value)
            file_base = f'dge_both_genomes_{sanitized_group_value}_mut_vs_ctrl' # Use group_value in filename
            dge_csv_path = os.path.join(overall_dge_dir, f'{file_base}.csv') # Use overall_dge_dir
            result_df.to_csv(dge_csv_path, index=False)
            print(f"    DGE results saved to {dge_csv_path}")

            # Save DGE results table (XLSX) - Use specific subdir
            if openpyxl:
                try:
                    dge_xlsx_path = os.path.join(overall_dge_dir, f'{file_base}.xlsx') # Use overall_dge_dir
                    result_df.to_excel(dge_xlsx_path, index=False, engine='openpyxl')
                    print(f"    DGE results saved to {dge_xlsx_path}")
                except Exception as e:
                    print(f"    Error saving DGE results to Excel for group {group_value}: {e}")
            else:
                 print("    Warning: 'openpyxl' not installed. Cannot save to .xlsx. Please install it (`pip install openpyxl`).")

            # --- Save significant DEG lists ---
            save_significant_degs(
                result_df=result_df,
                base_output_dir=overall_dge_dir, # Use the specific dir for this analysis type
                analysis_type='both_genomes',
                comparison_details='', # No genotype/condition subdir needed here
                group_name=group_value, # Pass the specific group value
                comparison_name='Mutant_vs_Control'
            )
            # --- End save significant DEG lists ---

            # print top N genes upregulated in Mutant
            top_n = 5 # Reduced for brevity
            print(f"    Top {top_n} genes upregulated in Mutant vs Control for group {group_value}: {result_df.head(top_n)['names'].tolist()}")
            print("---")
        except ValueError as e:
            print(f"    Skipping group {group_value} due to error during rank_genes_groups: {e}")
            dge_results[group_value] = (None, n_ctrl, n_mut)
            skipped_dge.append({grouping_key: group_value, 'reason': f'rank_genes_groups error: {e}'})
            print("---")

    # Save summary of skipped DGE runs - Use specific subdir
    if skipped_dge:
        skipped_df = pd.DataFrame(skipped_dge)
        skipped_output_path = os.path.join(overall_dge_dir, f'dge_overall_skipped_{grouping_key}.csv') # Use grouping_key in filename
        skipped_df.to_csv(skipped_output_path, index=False)
        print(f"\nSummary of skipped groups (using key '{grouping_key}') saved to {skipped_output_path}")

    print("\nBoth genomes DGE analysis finished.")

    # Visualize results
    print("\nGenerating Volcano Plots for Overall DGE results...")
    for group_value, (result_df, n_ctrl, n_mut) in dge_results.items():
        plot_volcano(
            result_df=result_df,
            group_name=group_value, # Pass the specific group value
            n1=n_mut,
            n2=n_ctrl,
            output_dir=overall_plot_dir, # Use overall_plot_dir
            filename_prefix="dge_conditions_comparison",
            group1_name="Mutant",
            group2_name="Control",
            title_prefix="Both genomes - "
        )
    print("\nBoth genomes DGE visualization finished.")

    return dge_results

def run_genotype_specific_dge(adata, grouping_key, dge_output_dir, plot_output_dir, layer=None):
    """
    Performs genotype-specific DGE analysis (Control vs Mutant within each genotype)
    across groups defined by grouping_key.

    Args:
        adata (ad.AnnData): Annotated data object.
        grouping_key (str): The column name in adata.obs to use for grouping (e.g., 'cell_type', 'leiden_0.8_alt').
        dge_output_dir (str): Base directory to save DGE result tables.
        plot_output_dir (str): Base directory to save volcano plots.
        layer (str, optional): Layer in adata to use for DGE. If None, uses adata.X. Defaults to None.
                             WARNING: If layer is None, adata.X is used. If adata.X was generated using
                             normalize_total, results may be unreliable for comparing groups with
                             large differences in raw counts/library size (e.g., Nestin_Mut vs Nestin_Ctrl).

    Returns:
        dict: Dictionary containing DGE results per genotype, then per group.
              Key: genotype ('Emx1', 'Nestin')
              Value: dict {group_value: tuple (results_df, n_ctrl, n_mut)}
    """
    # Define and create specific output subdirectories
    gs_dge_dir = os.path.join(dge_output_dir, 'geno_spec_cond_comp') # Shortened
    gs_plot_dir = os.path.join(plot_output_dir, 'geno_spec_cond_comp') # Shortened
    os.makedirs(gs_dge_dir, exist_ok=True)
    os.makedirs(gs_plot_dir, exist_ok=True)

    print("\nStarting Genotype-Specific DGE analysis...")

    # 1. Create combined genotype_condition column if it doesn't exist
    if 'genotype_condition' not in adata.obs.columns:
        adata.obs['genotype_condition'] = adata.obs['genotype'].astype(str) + '_' + adata.obs['condition'].astype(str)
        # Ensure it's categorical
        if not pd.api.types.is_categorical_dtype(adata.obs['genotype_condition']): # type: ignore
            adata.obs['genotype_condition'] = adata.obs['genotype_condition'].astype('category')
        print(f"Created 'genotype_condition' column with categories: {list(adata.obs['genotype_condition'].cat.categories)}")
    else:
        print("'genotype_condition' column already exists.")
        # Ensure it's categorical
        if not pd.api.types.is_categorical_dtype(adata.obs['genotype_condition']):  # type: ignore
            adata.obs['genotype_condition'] = adata.obs['genotype_condition'].astype('category')
        print(f"Existing categories: {list(adata.obs['genotype_condition'].cat.categories)}")


    # 2. Perform DGE for each genotype separately
    dge_by_genotype = {'Emx1': {}, 'Nestin': {}}
    # Check if grouping_key exists
    if grouping_key not in adata.obs.columns:
        raise ValueError(f"Grouping key '{grouping_key}' not found in adata.obs. Available columns: {list(adata.obs.columns)}")
    groups = adata.obs[grouping_key].unique()
    skipped_dge_genotype = {'Emx1': [], 'Nestin': []}

    for group_value in groups:
        # min_cells_for_gene = 3 # Minimum number of cells a gene must be expressed in to be kept for DGE - Replaced by percentage
        print(f"\nProcessing group: {group_value} (using key '{grouping_key}')")
        adata_group_orig = adata[adata.obs[grouping_key] == group_value].copy() # Keep original for Nestin
        counts = adata_group_orig.obs['genotype_condition'].value_counts()

        # --- Emx1 Genotype ---
        genotype = 'Emx1'
        print(f"  Running DGE for {genotype} genotype...")
        group_ctrl = f'{genotype}_Control'
        group_mut = f'{genotype}_Mutant'

        # Check counts
        n_ctrl = counts.get(group_ctrl, 0)
        n_mut = counts.get(group_mut, 0)

        if n_ctrl < MIN_CELLS_FOR_DGE_AND_PLOTTING or n_mut < MIN_CELLS_FOR_DGE_AND_PLOTTING:
            print(f"    Skipping {genotype} DGE for group {group_value}: Not enough cells (< {MIN_CELLS_FOR_DGE_AND_PLOTTING}). {group_ctrl}: {n_ctrl}, {group_mut}: {n_mut}")
            dge_by_genotype[genotype][group_value] = (None, n_ctrl, n_mut)
            skipped_dge_genotype[genotype].append({grouping_key: group_value, 'reason': f'Low cell count ({group_ctrl}={n_ctrl}, {group_mut}={n_mut})'})
        else:
            # Filter genes for this specific comparison (require expression in >= 10% of cells)
            adata_group_emx1 = adata_group_orig[adata_group_orig.obs['genotype'] == genotype].copy() # Subset further for filtering
            print(f"    Filtering genes for {genotype} (Initial: {adata_group_emx1.n_vars})")
            
            calculated_min_cells_by_percentage_emx1 = int(np.ceil(adata_group_emx1.n_obs * MIN_PERCENTAGE_OF_CELLS))
            effective_min_cells_emx1 = max(MIN_GENE_EXPRESSION_CELLS_ABSOLUTE, calculated_min_cells_by_percentage_emx1)

            print(f"    Filtering genes: Requiring expression in >= {effective_min_cells_emx1} cells (max of {MIN_GENE_EXPRESSION_CELLS_ABSOLUTE} absolute and {MIN_PERCENTAGE_OF_CELLS*100}% of {adata_group_emx1.n_obs})...")
            sc.pp.filter_genes(adata_group_emx1, min_cells=effective_min_cells_emx1)
            print(f"    Genes after filtering: {adata_group_emx1.n_vars}")
            if adata_group_emx1.n_vars == 0:
                print(f"    Skipping {genotype} for group {group_value}: No genes remaining after filtering.")
                dge_by_genotype[genotype][group_value] = (None, n_ctrl, n_mut)
                skipped_dge_genotype[genotype].append({grouping_key: group_value, 'reason': f'No genes left after filtering (min_cells={effective_min_cells_emx1}, 10%)'})
                # No continue here, need to let Nestin run below
            else:
                try:
                    if layer is None:
                        print(f"    WARNING: Running DGE on adata.X for {genotype}, group {group_value}. Ensure this data representation is appropriate for DGE (see function docstring warning).")
                    key_added = f'rank_genes_{genotype}'
                    sc.tl.rank_genes_groups(
                        adata_group_emx1, # Use filtered data
                        groupby='condition', # Compare Control vs Mutant within this genotype subset
                        groups=['Mutant'], # Group to compare
                        reference='Control', # Reference group
                        method='wilcoxon', # Changed from 't-test'
                        layer=layer, # Use specified layer
                        use_raw=False, # layer takes precedence over use_raw
                        n_genes=adata_group_emx1.n_vars, # Use filtered number of genes
                        pts=True, # Calculate fraction of cells expressing the gene
                        key_added=key_added
                )
                # Extract results for the specific comparison (Mutant vs Control)
                    result_df = sc.get.rank_genes_groups_df(adata_group_emx1, group='Mutant', key=key_added) # Get results for Mutant group

                    # --- Add mean expression columns ---
                    if not result_df.empty:
                        genes_in_result = result_df['names'].tolist()
                        valid_genes = [g for g in genes_in_result if g in adata_group_emx1.var_names]
                        if not valid_genes:
                             print(f"      Warning ({genotype}, group {group_value}): No genes from DGE results found in the subset anndata. Skipping mean expression calculation.")
                        else:
                            print(f"      Calculating mean expression for {len(valid_genes)} genes ({genotype}, group {group_value})...")
                            control_indices = adata_group_emx1.obs.index[adata_group_emx1.obs['condition'] == 'Control']
                            mutant_indices = adata_group_emx1.obs.index[adata_group_emx1.obs['condition'] == 'Mutant']

                            data_source = adata_group_emx1[:, valid_genes]
                            if layer and layer in data_source.layers:
                                data_matrix = data_source.layers[layer]
                                layer_name = layer
                            else:
                                data_matrix = data_source.X
                                layer_name = 'X'
                            print(f"      Using data from layer: {layer_name}")

                            if hasattr(data_matrix, "toarray"): 
                                data_matrix = data_matrix.toarray()
                            expr_df = pd.DataFrame(data_matrix, index=adata_group_emx1.obs_names, columns=valid_genes)

                            mean_expr_control = expr_df.loc[control_indices].mean(axis=0)
                            mean_expr_mutant = expr_df.loc[mutant_indices].mean(axis=0)

                            result_df[f'mean_expr_{genotype}_Control ({layer_name})'] = result_df['names'].map(mean_expr_control)
                            result_df[f'mean_expr_{genotype}_Mutant ({layer_name})'] = result_df['names'].map(mean_expr_mutant)
                            print(f"      Added mean expression columns ({genotype}, group {group_value}).")
                    # --- End mean expression ---

                    dge_by_genotype[genotype][group_value] = (result_df, n_ctrl, n_mut) # Store potentially modified df
                    print(f"    Finished DGE for {genotype}, group {group_value}. Found {len(result_df)} ranked genes.")

                    # Save DGE results table (CSV) - Use specific subdir
                    sanitized_group_value = sanitize_filename(group_value)
                    file_base = f'dge_{genotype}_{sanitized_group_value}_mut_vs_ctrl' # Use group_value in filename
                    dge_csv_path = os.path.join(gs_dge_dir, f'{file_base}.csv') # Use gs_dge_dir
                    result_df.to_csv(dge_csv_path, index=False)
                    print(f"    DGE results saved to {dge_csv_path}")

                    # Save DGE results table (XLSX) - Use specific subdir
                    if openpyxl:
                        try:
                            dge_xlsx_path = os.path.join(gs_dge_dir, f'{file_base}.xlsx') # Use gs_dge_dir
                            result_df.to_excel(dge_xlsx_path, index=False, engine='openpyxl')
                            print(f"    DGE results saved to {dge_xlsx_path}")
                        except Exception as e:
                            print(f"    Error saving DGE results to Excel for {genotype}, group {group_value}: {e}")
                    else:
                        print(f"    Warning ({genotype}, group {group_value}): 'openpyxl' not installed. Cannot save to .xlsx.")

                    # --- Save significant DEG lists ---
                    save_significant_degs(
                        result_df=result_df,
                        base_output_dir=gs_dge_dir, # Use the specific dir for this analysis type
                        analysis_type='genotype_specific',
                        comparison_details=genotype, # Subdir for genotype
                        group_name=group_value, # Pass the specific group value
                        comparison_name='Mutant_vs_Control'
                    )
                    # --- End save significant DEG lists ---

                except ValueError as e:
                     print(f"    Skipping {genotype} for group {group_value} due to error during rank_genes_groups: {e}")
                     dge_by_genotype[genotype][group_value] = (None, n_ctrl, n_mut)
                     skipped_dge_genotype[genotype].append({grouping_key: group_value, 'reason': f'rank_genes_groups error: {e}'})

        # --- Nestin Genotype ---
        genotype = 'Nestin'
        print(f"  Running DGE for {genotype} genotype...")
        group_ctrl = f'{genotype}_Control'
        group_mut = f'{genotype}_Mutant'

        # Check counts (reuse counts variable from above)
        n_ctrl = counts.get(group_ctrl, 0)
        n_mut = counts.get(group_mut, 0)

        if n_ctrl < MIN_CELLS_FOR_DGE_AND_PLOTTING or n_mut < MIN_CELLS_FOR_DGE_AND_PLOTTING:
            print(f"    Skipping {genotype} DGE for group {group_value}: Not enough cells (< {MIN_CELLS_FOR_DGE_AND_PLOTTING}). {group_ctrl}: {n_ctrl}, {group_mut}: {n_mut}")
            dge_by_genotype[genotype][group_value] = (None, n_ctrl, n_mut)
            skipped_dge_genotype[genotype].append({grouping_key: group_value, 'reason': f'Low cell count ({group_ctrl}={n_ctrl}, {group_mut}={n_mut})'})

        else:
            # Filter genes for this specific comparison (require expression in >= 10% of cells)
            adata_group_nes = adata_group_orig[adata_group_orig.obs['genotype'] == genotype].copy() # Subset further for filtering
            print(f"    Filtering genes for {genotype} (Initial: {adata_group_nes.n_vars})")

            calculated_min_cells_by_percentage_nes = int(np.ceil(adata_group_nes.n_obs * MIN_PERCENTAGE_OF_CELLS))
            effective_min_cells_nes = max(MIN_GENE_EXPRESSION_CELLS_ABSOLUTE, calculated_min_cells_by_percentage_nes)
            
            print(f"    Filtering genes: Requiring expression in >= {effective_min_cells_nes} cells (max of {MIN_GENE_EXPRESSION_CELLS_ABSOLUTE} absolute and {MIN_PERCENTAGE_OF_CELLS*100}% of {adata_group_nes.n_obs})...")
            sc.pp.filter_genes(adata_group_nes, min_cells=effective_min_cells_nes)
            print(f"    Genes after filtering: {adata_group_nes.n_vars}")
            if adata_group_nes.n_vars == 0:
                print(f"    Skipping {genotype} for group {group_value}: No genes remaining after filtering.")
                dge_by_genotype[genotype][group_value] = (None, n_ctrl, n_mut)
                skipped_dge_genotype[genotype].append({grouping_key: group_value, 'reason': f'No genes left after filtering (min_cells={effective_min_cells_nes}, 10%)'})
                # Continue to next group
            else:
                try:
                    if layer is None:
                         print(f"    WARNING: Running DGE on adata.X for {genotype}, group {group_value}. Ensure this data representation is appropriate for DGE (see function docstring warning).")
                    # Use different key to avoid overwrite
                    key_added = f'rank_genes_{genotype}'
                    # Need to re-run rank_genes_groups for the Nestin comparison on the same subset
                    sc.tl.rank_genes_groups(
                        adata_group_nes, # Use filtered data
                        groupby='condition', # Compare Control vs Mutant within this genotype subset
                        groups=['Mutant'], # Group to compare
                        reference='Control', # Reference group
                        method='wilcoxon', # Changed from 't-test'
                        layer=layer, # Use specified layer
                        use_raw=False, # layer takes precedence over use_raw
                        n_genes=adata_group_nes.n_vars, # Use filtered number of genes
                        pts=True, # Calculate fraction of cells expressing the gene
                        key_added=key_added
                )
                # Extract results for the specific comparison (Mutant vs Control)
                    result_df = sc.get.rank_genes_groups_df(adata_group_nes, group='Mutant', key=key_added) # Get results for Mutant group

                    # --- Add mean expression columns ---
                    if not result_df.empty:
                        genes_in_result = result_df['names'].tolist()
                        valid_genes = [g for g in genes_in_result if g in adata_group_nes.var_names]
                        if not valid_genes:
                             print(f"      Warning ({genotype}, group {group_value}): No genes from DGE results found in the subset anndata. Skipping mean expression calculation.")
                        else:
                            print(f"      Calculating mean expression for {len(valid_genes)} genes ({genotype}, group {group_value})...")
                            control_indices = adata_group_nes.obs.index[adata_group_nes.obs['condition'] == 'Control']
                            mutant_indices = adata_group_nes.obs.index[adata_group_nes.obs['condition'] == 'Mutant']

                            data_source = adata_group_nes[:, valid_genes]
                            if layer and layer in data_source.layers:
                                data_matrix = data_source.layers[layer]
                                layer_name = layer
                            else:
                                data_matrix = data_source.X
                                layer_name = 'X'
                            print(f"      Using data from layer: {layer_name}")

                            if hasattr(data_matrix, "toarray"): 
                                data_matrix = data_matrix.toarray()
                            expr_df = pd.DataFrame(data_matrix, index=adata_group_nes.obs_names, columns=valid_genes)

                            mean_expr_control = expr_df.loc[control_indices].mean(axis=0)
                            mean_expr_mutant = expr_df.loc[mutant_indices].mean(axis=0)

                            result_df[f'mean_expr_{genotype}_Control ({layer_name})'] = result_df['names'].map(mean_expr_control)
                            result_df[f'mean_expr_{genotype}_Mutant ({layer_name})'] = result_df['names'].map(mean_expr_mutant)
                            print(f"      Added mean expression columns ({genotype}, group {group_value}).")
                    # --- End mean expression ---

                    dge_by_genotype[genotype][group_value] = (result_df, n_ctrl, n_mut) # Store potentially modified df
                    print(f"    Finished DGE for {genotype}, group {group_value}. Found {len(result_df)} ranked genes.")

                    # Save DGE results table (CSV) - Use specific subdir
                    sanitized_group_value = sanitize_filename(group_value)
                    file_base = f'dge_{genotype}_{sanitized_group_value}_mut_vs_ctrl' # Use group_value in filename
                    dge_csv_path = os.path.join(gs_dge_dir, f'{file_base}.csv') # Use gs_dge_dir
                    result_df.to_csv(dge_csv_path, index=False)
                    print(f"    DGE results saved to {dge_csv_path}")

                    # Save DGE results table (XLSX) - Use specific subdir
                    if openpyxl:
                        try:
                            dge_xlsx_path = os.path.join(gs_dge_dir, f'{file_base}.xlsx') # Use gs_dge_dir
                            result_df.to_excel(dge_xlsx_path, index=False, engine='openpyxl')
                            print(f"    DGE results saved to {dge_xlsx_path}")
                        except Exception as e:
                            print(f"    Error saving DGE results to Excel for {genotype}, group {group_value}: {e}")
                    else:
                        print(f"    Warning ({genotype}, group {group_value}): 'openpyxl' not installed. Cannot save to .xlsx.")

                    # --- Save significant DEG lists ---
                    save_significant_degs(
                        result_df=result_df,
                        base_output_dir=gs_dge_dir, # Use the specific dir for this analysis type
                        analysis_type='genotype_specific',
                        comparison_details=genotype, # Subdir for genotype
                        group_name=group_value, # Pass the specific group value
                        comparison_name='Mutant_vs_Control'
                    )
                    # --- End save significant DEG lists ---

                except ValueError as e:
                     print(f"    Skipping {genotype} for group {group_value} due to error during rank_genes_groups: {e}")
                     dge_by_genotype[genotype][group_value] = (None, n_ctrl, n_mut)
                     skipped_dge_genotype[genotype].append({grouping_key: group_value, 'reason': f'rank_genes_groups error: {e}'})

    # Save summary of skipped DGE runs for each genotype
    for genotype, skipped_list in skipped_dge_genotype.items():
        if skipped_list:
            skipped_df = pd.DataFrame(skipped_list)
            skipped_output_path = os.path.join(gs_dge_dir, f'dge_{genotype}_skipped_{grouping_key}.csv') # Use grouping_key in filename
            skipped_df.to_csv(skipped_output_path, index=False)
            print(f"\nSummary of skipped groups (using key '{grouping_key}') for {genotype} saved to {skipped_output_path}")

    print("\nGenotype-Specific DGE analysis finished.")

    # Visualize results
    print("\nGenerating Volcano Plots for Genotype-Specific DGE results...")
    for genotype in dge_by_genotype:
        print(f"\n--- Visualizing for Genotype: {genotype} ---")
        genotype_results = dge_by_genotype[genotype]

        for group_value, (result_df, n_ctrl, n_mut) in genotype_results.items():
             plot_volcano(
                result_df=result_df,
                group_name=group_value, # Pass the specific group value
                n1=n_mut,
                n2=n_ctrl,
                output_dir=gs_plot_dir, # Use gs_plot_dir
                filename_prefix=f"dge_{genotype}",
                group1_name=f"{genotype}_Mut",
                group2_name=f"{genotype}_Ctrl",
                title_prefix=f"{genotype} - "
             )
    print("\nGenotype-Specific DGE visualization finished.")

    return dge_by_genotype

def run_genotype_comparison_dge(adata, grouping_key, dge_output_dir, plot_output_dir, layer=None):
    """
    Performs DGE analysis comparing genotypes (Nestin vs Emx1) within each condition,
    across groups defined by grouping_key.

    Args:
        adata (ad.AnnData): Annotated data object. Must have 'genotype_condition' column.
        grouping_key (str): The column name in adata.obs to use for grouping (e.g., 'cell_type', 'leiden_0.8_alt').
        dge_output_dir (str): Base directory to save DGE result tables.
        plot_output_dir (str): Base directory to save volcano plots.
        layer (str, optional): Layer in adata to use for DGE. If None, uses adata.X. Defaults to None.
                             WARNING: If layer is None, adata.X is used. If adata.X was generated using
                             normalize_total, results may be unreliable for comparing groups with
                             large differences in raw counts/library size. Assess count differences
                             between genotypes within conditions before interpreting results.

    Returns:
        dict: Dictionary containing DGE results per condition, then per group.
              Key: condition ('Control', 'Mutant')
              Value: dict {group_value: tuple (results_df, n_ref_genotype, n_comp_genotype)}
    """
    # Define and create specific output subdirectories
    gc_dge_dir = os.path.join(dge_output_dir, 'cond_spec_geno_comp')
    gc_plot_dir = os.path.join(plot_output_dir, 'cond_spec_geno_comp')
    os.makedirs(gc_dge_dir, exist_ok=True)
    os.makedirs(gc_plot_dir, exist_ok=True)

    print("\nStarting DGE analysis: Genotype comparison within conditions...")

    # Ensure 'genotype_condition' column exists (should be created by run_genotype_specific_dge if needed)
    if 'genotype_condition' not in adata.obs.columns:
         adata.obs['genotype_condition'] = adata.obs['genotype'].astype(str) + '_' + adata.obs['condition'].astype(str)
         if not pd.api.types.is_categorical_dtype(adata.obs['genotype_condition']): # type: ignore
            adata.obs['genotype_condition'] = adata.obs['genotype_condition'].astype('category')
         print(f"Created 'genotype_condition' column (was missing). Categories: {list(adata.obs['genotype_condition'].cat.categories)}")


    # Result structure: {condition: {group_value: (df, n_genotype1, n_genotype2)}}
    dge_genotype_within_condition = {
        'Control': {},
        'Mutant': {}
    }
    skipped_dge_within_condition = {'Control': [], 'Mutant': []}

    # Define reference for consistency (e.g., compare Nestin against Emx1)
    ref_genotype = 'Emx1'
    comp_genotype = 'Nestin'

    # Check if grouping_key exists
    if grouping_key not in adata.obs.columns:
        raise ValueError(f"Grouping key '{grouping_key}' not found in adata.obs. Available columns: {list(adata.obs.columns)}")
    groups = adata.obs[grouping_key].unique()

    for group_value in groups:
        # min_cells_for_gene = 3 # Minimum number of cells a gene must be expressed in to be kept for DGE - Replaced by percentage
        print(f"\nProcessing group: {group_value} (using key '{grouping_key}')")
        adata_group_orig = adata[adata.obs[grouping_key] == group_value].copy() # Keep original for Mutant condition
        counts = adata_group_orig.obs['genotype_condition'].value_counts()

        # --- Control Condition ---
        condition = 'Control'
        print(f"  Running DGE for {condition} condition ({comp_genotype} vs {ref_genotype})...")
        group_ref = f'{ref_genotype}_{condition}' # e.g., Emx1_Control
        group_comp = f'{comp_genotype}_{condition}' # e.g., Nestin_Control

        n_ref = counts.get(group_ref, 0)
        n_comp = counts.get(group_comp, 0)

        if n_ref < MIN_CELLS_FOR_DGE_AND_PLOTTING or n_comp < MIN_CELLS_FOR_DGE_AND_PLOTTING:
            print(f"    Skipping DGE for {condition} condition, group {group_value}: Not enough cells (< {MIN_CELLS_FOR_DGE_AND_PLOTTING}). {group_ref}: {n_ref}, {group_comp}: {n_comp}")
            dge_genotype_within_condition[condition][group_value] = (None, n_ref, n_comp)
            skipped_dge_within_condition[condition].append({grouping_key: group_value, 'reason': f'Low cell count ({group_ref}={n_ref}, {group_comp}={n_comp})'})
        else:
            # Subset for the specific condition *before* filtering
            adata_group_ctrl = adata_group_orig[adata_group_orig.obs['condition'] == condition].copy()

            # Filter genes expressed in too few cells within this subset (require expression in >= 10% of cells)
            print(f"    Initial genes for {condition} condition: {adata_group_ctrl.n_vars}")
            
            calculated_min_cells_by_percentage_ctrl = int(np.ceil(adata_group_ctrl.n_obs * MIN_PERCENTAGE_OF_CELLS))
            effective_min_cells_ctrl = max(MIN_GENE_EXPRESSION_CELLS_ABSOLUTE, calculated_min_cells_by_percentage_ctrl)

            print(f"    Filtering genes for {condition}: Requiring expression in >= {effective_min_cells_ctrl} cells (max of {MIN_GENE_EXPRESSION_CELLS_ABSOLUTE} absolute and {MIN_PERCENTAGE_OF_CELLS*100}% of {adata_group_ctrl.n_obs})...")
            sc.pp.filter_genes(adata_group_ctrl, min_cells=effective_min_cells_ctrl)
            print(f"    Genes after filtering for {condition}: {adata_group_ctrl.n_vars}")

            if adata_group_ctrl.n_vars == 0:
                print(f"    Skipping {condition} for group {group_value}: No genes remaining after filtering.")
                dge_genotype_within_condition[condition][group_value] = (None, n_ref, n_comp)
                skipped_dge_within_condition[condition].append({grouping_key: group_value, 'reason': f'No genes left after filtering (min_cells={effective_min_cells_ctrl}, 10%)'})
                # Don't continue here, still need to process Mutant condition below
            else:
                try:
                    if layer is None:
                        print(f"    WARNING: Running DGE on adata.X for {condition} condition (Control), group {group_value}. Ensure this data representation is appropriate for DGE (see function docstring warning).")
                    # Use a unique key for this comparison type
                    key_added_ctrl = f'rank_genes_geno_{condition.lower()}'
                    sc.tl.rank_genes_groups(
                        adata_group_ctrl, # Use filtered data
                        groupby='genotype', # Compare Nestin vs Emx1 within this condition subset
                        groups=[comp_genotype], # Group to compare (Nestin)
                        reference=ref_genotype,  # Reference group (Emx1)
                        method='wilcoxon', # Changed from 't-test'
                        layer=layer, # Use specified layer
                        use_raw=False, # layer takes precedence over use_raw
                        n_genes=adata_group_ctrl.n_vars, # Use filtered number of genes
                        pts=True, # Calculate fraction of cells expressing the gene
                        key_added=key_added_ctrl
                )
                    result_df = sc.get.rank_genes_groups_df(adata_group_ctrl, group=comp_genotype, key=key_added_ctrl) # Get results for comp_genotype

                    # --- Add mean expression columns ---
                    if not result_df.empty:
                        genes_in_result = result_df['names'].tolist()
                        valid_genes = [g for g in genes_in_result if g in adata_group_ctrl.var_names]
                        if not valid_genes:
                             print(f"      Warning ({condition} condition, group {group_value}): No genes from DGE results found. Skipping mean expression calculation.")
                        else:
                            print(f"      Calculating mean expression for {len(valid_genes)} genes ({condition} condition, group {group_value})...")
                            ref_indices = adata_group_ctrl.obs.index[adata_group_ctrl.obs['genotype'] == ref_genotype]
                            comp_indices = adata_group_ctrl.obs.index[adata_group_ctrl.obs['genotype'] == comp_genotype]

                            data_source = adata_group_ctrl[:, valid_genes]
                            if layer and layer in data_source.layers:
                                data_matrix = data_source.layers[layer]
                                layer_name = layer
                            else:
                                data_matrix = data_source.X
                                layer_name = 'X'
                            print(f"      Using data from layer: {layer_name}")

                            if hasattr(data_matrix, "toarray"): 
                                data_matrix = data_matrix.toarray()
                            expr_df = pd.DataFrame(data_matrix, index=adata_group_ctrl.obs_names, columns=valid_genes)

                            mean_expr_ref = expr_df.loc[ref_indices].mean(axis=0)
                            mean_expr_comp = expr_df.loc[comp_indices].mean(axis=0)

                            result_df[f'mean_expr_{ref_genotype}_{condition} ({layer_name})'] = result_df['names'].map(mean_expr_ref)
                            result_df[f'mean_expr_{comp_genotype}_{condition} ({layer_name})'] = result_df['names'].map(mean_expr_comp)
                            print(f"      Added mean expression columns ({condition} condition, group {group_value}).")
                    # --- End mean expression ---

                    dge_genotype_within_condition[condition][group_value] = (result_df, n_ref, n_comp) # Store potentially modified df
                    print(f"    Finished DGE for {condition}, group {group_value}. Found {len(result_df)} ranked genes.")

                    # Save DGE results table (CSV) - Use specific subdir
                    sanitized_group_value = sanitize_filename(group_value)
                    file_base = f'dge_{condition}_cond_{sanitized_group_value}_{comp_genotype}_vs_{ref_genotype}' # Use group_value in filename
                    dge_csv_path = os.path.join(gc_dge_dir, f'{file_base}.csv') # Use gc_dge_dir
                    result_df.to_csv(dge_csv_path, index=False)
                    print(f"    DGE results saved to {dge_csv_path}")

                    # Save DGE results table (XLSX) - Use specific subdir
                    if openpyxl:
                        try:
                            dge_xlsx_path = os.path.join(gc_dge_dir, f'{file_base}.xlsx') # Use gc_dge_dir
                            result_df.to_excel(dge_xlsx_path, index=False, engine='openpyxl')
                            print(f"    DGE results saved to {dge_xlsx_path}")
                        except Exception as e:
                            print(f"    Error saving DGE results to Excel for {condition} condition, group {group_value}: {e}")
                    else:
                        print(f"    Warning ({condition} condition, group {group_value}): 'openpyxl' not installed. Cannot save to .xlsx.")

                    # --- Save significant DEG lists ---
                    save_significant_degs(
                        result_df=result_df,
                        base_output_dir=gc_dge_dir, # Use the specific dir for this analysis type
                        analysis_type='genotype_comparison',
                        comparison_details=condition, # Subdir for condition
                        group_name=group_value, # Pass the specific group value
                        comparison_name=f'{comp_genotype}_vs_{ref_genotype}' # e.g., Nestin_vs_Emx1
                    )
                    # --- End save significant DEG lists ---

                except ValueError as e:
                     print(f"    Skipping {condition} for group {group_value} due to error during rank_genes_groups: {e}")
                     dge_genotype_within_condition[condition][group_value] = (None, n_ref, n_comp)
                     skipped_dge_within_condition[condition].append({grouping_key: group_value, 'reason': f'rank_genes_groups error: {e}'})

        # --- Mutant Condition ---
        condition = 'Mutant'
        print(f"  Running DGE for {condition} condition ({comp_genotype} vs {ref_genotype})...")
        group_ref = f'{ref_genotype}_{condition}' # e.g., Emx1_Mutant
        group_comp = f'{comp_genotype}_{condition}' # e.g., Nestin_Mutant

        # Reuse counts calculated above
        n_ref = counts.get(group_ref, 0)
        n_comp = counts.get(group_comp, 0)

        if n_ref < MIN_CELLS_FOR_DGE_AND_PLOTTING or n_comp < MIN_CELLS_FOR_DGE_AND_PLOTTING:
            print(f"    Skipping DGE for {condition} condition, group {group_value}: Not enough cells (< {MIN_CELLS_FOR_DGE_AND_PLOTTING}). {group_ref}: {n_ref}, {group_comp}: {n_comp}")
            dge_genotype_within_condition[condition][group_value] = (None, n_ref, n_comp)
            skipped_dge_within_condition[condition].append({grouping_key: group_value, 'reason': f'Low cell count ({group_ref}={n_ref}, {group_comp}={n_comp})'})
        else:
            # Subset for the specific condition *before* filtering
            adata_group_mut = adata_group_orig[adata_group_orig.obs['condition'] == condition].copy()

            # Filter genes expressed in too few cells within this subset (require expression in >= 10% of cells)
            print(f"    Initial genes for {condition} condition: {adata_group_mut.n_vars}")
            
            calculated_min_cells_by_percentage_mut = int(np.ceil(adata_group_mut.n_obs * MIN_PERCENTAGE_OF_CELLS))
            effective_min_cells_mut = max(MIN_GENE_EXPRESSION_CELLS_ABSOLUTE, calculated_min_cells_by_percentage_mut)

            print(f"    Filtering genes for {condition}: Requiring expression in >= {effective_min_cells_mut} cells (max of {MIN_GENE_EXPRESSION_CELLS_ABSOLUTE} absolute and {MIN_PERCENTAGE_OF_CELLS*100}% of {adata_group_mut.n_obs})...")
            sc.pp.filter_genes(adata_group_mut, min_cells=effective_min_cells_mut)
            print(f"    Genes after filtering for {condition}: {adata_group_mut.n_vars}")

            if adata_group_mut.n_vars == 0:
                print(f"    Skipping {condition} for group {group_value}: No genes remaining after filtering.")
                dge_genotype_within_condition[condition][group_value] = (None, n_ref, n_comp)
                skipped_dge_within_condition[condition].append({grouping_key: group_value, 'reason': f'No genes left after filtering (min_cells={effective_min_cells_mut}, 10%)'})
                continue # Skip to next group
            else:
                try:
                    if layer is None:
                        print(f"    WARNING: Running DGE on adata.X for {condition} condition (Mutant), group {group_value}. Ensure this data representation is appropriate for DGE (see function docstring warning).")
                    # Use a unique key for this comparison type
                    key_added_mut = f'rank_genes_geno_{condition.lower()}'
                    sc.tl.rank_genes_groups(
                        adata_group_mut, # Use filtered data
                        groupby='genotype', # Compare Nestin vs Emx1 within this condition subset
                        groups=[comp_genotype], # Group to compare (Nestin)
                        reference=ref_genotype,  # Reference group (Emx1)
                        method='wilcoxon', # Changed from 't-test'
                        layer=layer, # Use specified layer
                        use_raw=False, # layer takes precedence over use_raw
                        n_genes=adata_group_mut.n_vars, # Use filtered number of genes
                        pts=True, # Calculate fraction of cells expressing the gene
                        key_added=key_added_mut
                )
                    result_df = sc.get.rank_genes_groups_df(adata_group_mut, group=comp_genotype, key=key_added_mut) # Get results for comp_genotype

                    # --- Add mean expression columns ---
                    if not result_df.empty:
                        genes_in_result = result_df['names'].tolist()
                        valid_genes = [g for g in genes_in_result if g in adata_group_mut.var_names]
                        if not valid_genes:
                             print(f"      Warning ({condition} condition, group {group_value}): No genes from DGE results found. Skipping mean expression calculation.")
                        else:
                            print(f"      Calculating mean expression for {len(valid_genes)} genes ({condition} condition, group {group_value})...")
                            ref_indices = adata_group_mut.obs.index[adata_group_mut.obs['genotype'] == ref_genotype]
                            comp_indices = adata_group_mut.obs.index[adata_group_mut.obs['genotype'] == comp_genotype]

                            data_source = adata_group_mut[:, valid_genes]
                            if layer and layer in data_source.layers:
                                data_matrix = data_source.layers[layer]
                                layer_name = layer
                            else:
                                data_matrix = data_source.X
                                layer_name = 'X'
                            print(f"      Using data from layer: {layer_name}")

                            if hasattr(data_matrix, "toarray"): 
                                data_matrix = data_matrix.toarray()
                            expr_df = pd.DataFrame(data_matrix, index=adata_group_mut.obs_names, columns=valid_genes)

                            mean_expr_ref = expr_df.loc[ref_indices].mean(axis=0)
                            mean_expr_comp = expr_df.loc[comp_indices].mean(axis=0)

                            result_df[f'mean_expr_{ref_genotype}_{condition} ({layer_name})'] = result_df['names'].map(mean_expr_ref)
                            result_df[f'mean_expr_{comp_genotype}_{condition} ({layer_name})'] = result_df['names'].map(mean_expr_comp)
                            print(f"      Added mean expression columns ({condition} condition, group {group_value}).")
                    # --- End mean expression ---

                    dge_genotype_within_condition[condition][group_value] = (result_df, n_ref, n_comp) # Store potentially modified df
                    print(f"    Finished DGE for {condition}, group {group_value}. Found {len(result_df)} ranked genes.")

                    # Save DGE results table (CSV) - Use specific subdir
                    sanitized_group_value = sanitize_filename(group_value)
                    file_base = f'dge_{condition}_cond_{sanitized_group_value}_{comp_genotype}_vs_{ref_genotype}' # Use group_value in filename
                    dge_csv_path = os.path.join(gc_dge_dir, f'{file_base}.csv') # Use gc_dge_dir
                    result_df.to_csv(dge_csv_path, index=False)
                    print(f"    DGE results saved to {dge_csv_path}")

                    # Save DGE results table (XLSX) - Use specific subdir
                    if openpyxl:
                        try:
                            dge_xlsx_path = os.path.join(gc_dge_dir, f'{file_base}.xlsx') # Use gc_dge_dir
                            result_df.to_excel(dge_xlsx_path, index=False, engine='openpyxl')
                            print(f"    DGE results saved to {dge_xlsx_path}")
                        except Exception as e:
                            print(f"    Error saving DGE results to Excel for {condition} condition, group {group_value}: {e}")
                    else:
                        print(f"    Warning ({condition} condition, group {group_value}): 'openpyxl' not installed. Cannot save to .xlsx.")

                    # --- Save significant DEG lists ---
                    save_significant_degs(
                        result_df=result_df,
                        base_output_dir=gc_dge_dir, # Use the specific dir for this analysis type
                        analysis_type='genotype_comparison',
                        comparison_details=condition, # Subdir for condition
                        group_name=group_value, # Pass the specific group value
                        comparison_name=f'{comp_genotype}_vs_{ref_genotype}' # e.g., Nestin_vs_Emx1
                    )
                    # --- End save significant DEG lists ---

                except ValueError as e:
                     print(f"    Skipping {condition} for group {group_value} due to error during rank_genes_groups: {e}")
                     dge_genotype_within_condition[condition][group_value] = (None, n_ref, n_comp)
                     skipped_dge_within_condition[condition].append({grouping_key: group_value, 'reason': f'rank_genes_groups error: {e}'})

    # Save summary of skipped DGE runs for each condition comparison
    for condition, skipped_list in skipped_dge_within_condition.items():
        if skipped_list:
            skipped_df = pd.DataFrame(skipped_list)
            skipped_output_path = os.path.join(gc_dge_dir, f'dge_{condition}_cond_comparison_skipped_{grouping_key}.csv') # Use grouping_key in filename
            skipped_df.to_csv(skipped_output_path, index=False)
            print(f"\nSummary of skipped groups (using key '{grouping_key}') for {condition} condition comparison saved to {skipped_output_path}")

    print("\nDGE analysis (Genotype comparison within conditions) finished.")

    # Visualize results
    print("\nGenerating Volcano Plots for Genotype comparison within conditions...")
    for condition in dge_genotype_within_condition:
        print(f"\n--- Visualizing for Condition: {condition} ({comp_genotype} vs {ref_genotype}) ---")
        condition_results = dge_genotype_within_condition[condition]

        for group_value, (result_df, n_ref, n_comp) in condition_results.items():
            plot_volcano(
                result_df=result_df,
                group_name=group_value, # Pass the specific group value
                n1=n_comp, # Nestin count
                n2=n_ref,  # Emx1 count
                output_dir=gc_plot_dir, # Use gc_plot_dir
                filename_prefix=f"dge_{condition}_cond_{comp_genotype}_vs_{ref_genotype}",
                group1_name=f"{comp_genotype}_{condition}", # e.g. Nestin_Control
                group2_name=f"{ref_genotype}_{condition}",  # e.g. Emx1_Control
                title_prefix=f"{condition} Cond - ",
                up_color='purple', # Custom colors for this comparison
                down_color='orange'
            )
    print("\nDGE visualization (Genotype comparison within conditions) finished.")

    return dge_genotype_within_condition

# --- Cluster Comparison DGE ---

def run_cluster_comparison_dge(adata, grouping_key, dge_output_dir, plot_output_dir, layer=None, method='t-test'):
    """
    Performs DGE analysis comparing each cluster defined by grouping_key against all other cells.
    This identifies marker genes for each cluster. Saves top 50 markers and generates a summary plot.

    Args:
        adata (ad.AnnData): Annotated data object.
        grouping_key (str): The column name in adata.obs to use for grouping (e.g., 'leiden_0.8_alt').
        dge_output_dir (str): Base directory to save DGE result tables.
        plot_output_dir (str): Base directory to save plots. Passed to sc.pl.rank_genes_groups.
        layer (str, optional): Layer in adata to use for DGE. If None, uses adata.X. Defaults to None.
                             WARNING: If layer is None, adata.X is used. If adata.X was generated using
                             normalize_total, results may be less suitable for quantitative fold-change
                             interpretation compared to methods using raw counts, but still useful for
                             identifying potential marker genes.
        method (str): Method for DGE ('t-test' or 'wilcoxon'). Defaults to 't-test'.

    Returns:
        dict: Dictionary containing DGE results per cluster.
              Key: cluster_id
              Value: tuple (results_df, n_cells_cluster, n_cells_rest)
              The results_df returned contains ALL ranked genes, even though only top 50 are saved to files.
    """
    # Define and create specific output subdirectories
    # The dge_output_dir passed from raw_16 is already the target 'biomarkers' directory
    cluster_dge_dir = dge_output_dir
    cluster_plot_dir = plot_output_dir # Use the plot_output_dir passed from raw_16
    os.makedirs(cluster_dge_dir, exist_ok=True)
    os.makedirs(cluster_plot_dir, exist_ok=True)

    print(f"\nStarting Cluster Comparison DGE analysis (Markers for '{grouping_key}') into {cluster_dge_dir}...")

    # Check if grouping_key exists
    if grouping_key not in adata.obs.columns:
        raise ValueError(f"Grouping key '{grouping_key}' not found in adata.obs. Available columns: {list(adata.obs.columns)}")

    # --- Filter Genes Globally ---
    # Filter genes expressed in too few cells across the *entire dataset* before running rank_genes_groups
    # This reduces computational load and focuses on more reliable markers.
    adata_filt = adata.copy() # Work on a copy to avoid modifying the original adata object here
    print(f"  Initial genes: {adata_filt.n_vars}")
    # Filter genes expressed in < 3 cells globally
    # Alternatively, use a percentage threshold like 1% or 5% of total cells
    # min_cells_global = int(np.ceil(adata_filt.n_obs * 0.01)) # Example: 1% threshold
    print(f"  Filtering genes globally: Requiring expression in >= {MIN_CELLS_GLOBAL_FILTER} cells...")
    sc.pp.filter_genes(adata_filt, min_cells=MIN_CELLS_GLOBAL_FILTER)
    print(f"  Genes after global filtering: {adata_filt.n_vars}")
    if adata_filt.n_vars == 0:
        print(f"  Skipping cluster comparison DGE: No genes remaining after global filtering (min_cells={MIN_CELLS_GLOBAL_FILTER}).")
        return {}
    # --- End Global Gene Filtering ---

    # Run DGE comparing each cluster vs the rest
    # This modifies the adata_filt object by adding the results to .uns
    # Filter out groups with too few cells before running rank_genes_groups
    min_cells_for_cluster_dge = 3 # Minimum cells required for a cluster to be included in DGE
    initial_groups = adata_filt.obs[grouping_key].unique()
    valid_groups = []
    removed_groups = []

    for group in initial_groups:
        n_cells_in_group = (adata_filt.obs[grouping_key] == group).sum()
        if n_cells_in_group >= min_cells_for_cluster_dge:
            valid_groups.append(group)
        else:
            removed_groups.append(group)
            print(f"  Warning: Skipping group '{group}' from DGE as it only contains {n_cells_in_group} cells (less than {min_cells_for_cluster_dge}).")

    if not valid_groups:
        print(f"  Skipping cluster comparison DGE: No groups remaining after filtering for minimum cell count ({min_cells_for_cluster_dge}).")
        return {}

    # Filter adata_filt to only include cells from valid groups
    adata_filt = adata_filt[adata_filt.obs[grouping_key].isin(valid_groups)].copy()
    # Re-categorize the grouping_key to remove empty categories
    adata_filt.obs[grouping_key] = adata_filt.obs[grouping_key].cat.remove_unused_categories()

    print(f"  Proceeding with DGE for {len(valid_groups)} groups. Removed groups: {removed_groups}")

    print(f"  Running sc.tl.rank_genes_groups (groupby='{grouping_key}', method='{method}')...")
    try:
        if layer is None:
             print("WARNING: Running cluster marker DGE on adata.X. Ensure this data representation is appropriate (see function docstring warning).")
        sc.tl.rank_genes_groups(
            adata_filt, # Use the globally filtered AnnData
            groupby=grouping_key,
            method=method,# type: ignore
            layer=layer,
            use_raw=False, # layer takes precedence over use_raw
            n_genes=adata_filt.n_vars, # Rank all remaining genes
            pts=True # Calculate fraction of cells expressing the gene in cluster vs rest
        )
    except Exception as e:
         print(f"  Error during rank_genes_groups for cluster comparison: {e}")
         return {} # Return empty dict on error

    print("  rank_genes_groups finished.")
    # Extract results
    dge_results_clusters = {}
    rank_genes_key = 'rank_genes_groups' # Default key where results are stored
    if rank_genes_key not in adata_filt.uns:
        print(f"  Error: '{rank_genes_key}' not found in adata_filt.uns after running sc.tl.rank_genes_groups.")
        return {}
    groups = adata_filt.uns[rank_genes_key]['names'].dtype.names # Get cluster names from results
    skipped_clusters = []

    for group_value in groups:
        # Extract results for this specific cluster
        try:
            result_df = sc.get.rank_genes_groups_df(adata_filt, group=str(group_value), key=rank_genes_key) # Ensure group_value is string
            if result_df is None or result_df.empty:
                print(f"    Skipping cluster {group_value}: No DGE results found.")
                skipped_clusters.append({grouping_key: group_value, 'reason': 'No DGE results from rank_genes_groups_df'})
                dge_results_clusters[group_value] = (None, 0, 0) # Indicate no results
                continue

            # --- Calculate mean expression for cluster vs rest ---
            print(f"    Calculating mean expression for Cluster {group_value} vs Rest...")
            # Get cell indices for the current cluster and the rest
            cluster_indices = adata_filt.obs.index[adata_filt.obs[grouping_key] == group_value]
            rest_indices = adata_filt.obs.index[adata_filt.obs[grouping_key] != group_value]

            n_cluster = len(cluster_indices)
            n_rest = len(rest_indices)

            if n_cluster == 0 or n_rest == 0:
                print(f"    Skipping cluster {group_value}: Zero cells in cluster ({n_cluster}) or rest ({n_rest}).")
                skipped_clusters.append({grouping_key: group_value, 'reason': f'Zero cells (Cluster={n_cluster}, Rest={n_rest})'})
                dge_results_clusters[group_value] = (result_df, n_cluster, n_rest) # Store counts even if skipping means
                continue


            genes_in_result = result_df['names'].tolist()
            valid_genes = [g for g in genes_in_result if g in adata_filt.var_names]

            if not valid_genes:
                print(f"      Warning (Cluster {group_value}): No genes from DGE results found in the AnnData. Skipping mean expression calculation.")
                mean_expr_cluster = pd.Series(index=result_df['names'], dtype=float) # Empty series
                mean_expr_rest = pd.Series(index=result_df['names'], dtype=float)    # Empty series
            else:
                # Select the data layer (use the same layer as DGE)
                data_source = adata_filt[:, valid_genes] # Use the filtered anndata and valid genes
                if layer and layer in data_source.layers:
                    data_matrix = data_source.layers[layer]
                    layer_name = layer
                else:
                    data_matrix = data_source.X
                    layer_name = 'X'
                print(f"      Using data from layer: {layer_name}")

                # Ensure data_matrix is dense for indexing if it's sparse
                if hasattr(data_matrix, "toarray"):
                     data_matrix = data_matrix.toarray()

                # Create a DataFrame for easier indexing using the filtered anndata's obs_names
                expr_df = pd.DataFrame(data_matrix, index=adata_filt.obs_names, columns=valid_genes)

                # Calculate mean expression for cluster and rest
                mean_expr_cluster = expr_df.loc[cluster_indices].mean(axis=0)
                mean_expr_rest = expr_df.loc[rest_indices].mean(axis=0)

            # Add to result_df, aligning by gene name ('names' column)
            # Use group_value directly in column name
            result_df[f'mean_expr_Cluster_{group_value} ({layer_name})'] = result_df['names'].map(mean_expr_cluster)
            result_df[f'mean_expr_Rest ({layer_name})'] = result_df['names'].map(mean_expr_rest)
            print(f"      Added mean expression columns for Cluster {group_value}.")
            # --- End mean expression calculation ---

            dge_results_clusters[group_value] = (result_df, n_cluster, n_rest)

            # --- Save Top 50 Markers ---
            top_50_df = result_df.head(50).copy()
            sanitized_group_value = sanitize_filename(group_value)
            file_base_top50 = f'dge_cluster_comparison_{sanitized_group_value}_vs_Rest_top50' # Filename for top 50

            # Save Top 50 DGE results table (CSV)
            dge_csv_path_top50 = os.path.join(cluster_dge_dir, f'{file_base_top50}.csv')
            top_50_df.to_csv(dge_csv_path_top50, index=False)
            print(f"    Top 50 Cluster marker results saved to {dge_csv_path_top50}")

            # Save Top 50 DGE results table (XLSX)
            if openpyxl:
                try:
                    dge_xlsx_path_top50 = os.path.join(cluster_dge_dir, f'{file_base_top50}.xlsx')
                    top_50_df.to_excel(dge_xlsx_path_top50, index=False, engine='openpyxl')
                    print(f"    Top 50 Cluster marker results saved to {dge_xlsx_path_top50}")
                except Exception as e:
                    print(f"    Error saving Top 50 cluster marker results to Excel for cluster {group_value}: {e}")
            else:
                print(f"    Warning (Cluster {group_value}): 'openpyxl' not installed. Cannot save Top 50 to .xlsx.")
            # --- End Save Top 50 Markers ---


            # --- Save significant marker gene lists (based on FULL results) ---
            save_significant_degs(
                result_df=result_df, # Use full results for significance filtering
                base_output_dir=cluster_dge_dir, # Pass the corrected base biomarker directory
                analysis_type='cluster_comparison',
                comparison_details='', # Group name is already used below
                group_name=group_value, # Pass the specific cluster ID
                comparison_name=f'Cluster_{group_value}_vs_Rest' # Name indicates cluster vs rest
            )
            # --- End save significant DEG lists ---

        except KeyError as e:
            print(f"    Skipping cluster {group_value}: KeyError during result extraction - {e}. Might indicate issues with rank_genes_groups output structure.")
            skipped_clusters.append({grouping_key: group_value, 'reason': f'KeyError during extraction: {e}'})
            dge_results_clusters[group_value] = (None, 0, 0) # Indicate no results
        except Exception as e: # Catch other potential errors
            print(f"    Skipping cluster {group_value}: Error during processing - {e}")
            skipped_clusters.append({grouping_key: group_value, 'reason': f'General processing error: {e}'})
            dge_results_clusters[group_value] = (None, 0, 0) # Indicate no results


    # Save summary of skipped clusters
    if skipped_clusters:
        skipped_df = pd.DataFrame(skipped_clusters)
        skipped_output_path = os.path.join(cluster_dge_dir, f'dge_cluster_comparison_skipped_{grouping_key}.csv')
        skipped_df.to_csv(skipped_output_path, index=False)
        print(f"\nSummary of skipped clusters (markers for '{grouping_key}') saved to {skipped_output_path}")

    print("\nCluster Comparison DGE analysis (Marker identification) finished.")

    # --- Generate and Save Marker Plot ---
    # Use the plot_output_dir passed to the function
    # Check if rank_genes_groups results exist before plotting
    if rank_genes_key in adata_filt.uns and groups: # Check if results exist and groups were found
        print(f"\nGenerating marker plot for '{grouping_key}' clusters...")
        original_figdir = sc.settings.figdir
        try:
            # Set figdir temporarily to ensure plot saves in the correct location (cluster_plot_dir)
            sc.settings.figdir = cluster_plot_dir
            plot_filename_suffix = f'_cluster_markers_{grouping_key}.png'
            # The filename saved by scanpy automatically includes 'rank_genes_groups_' prefix
            expected_plot_filename = f"rank_genes_groups_{plot_filename_suffix[1:]}" # Remove leading '_'
            sc.pl.rank_genes_groups(adata_filt, n_genes=25, sharey=False, key=rank_genes_key, save=plot_filename_suffix, show=False) # type: ignore
            print(f"  Marker plot saved to {os.path.join(cluster_plot_dir, expected_plot_filename)}")
        except Exception as e:
            print(f"  Error generating marker plot: {e}")
        finally:
            sc.settings.figdir = original_figdir # Restore original setting
    else:
        print("\nSkipping marker plot generation: No rank_genes_groups results found in AnnData.")
    # --- End Marker Plot Generation ---


    return dge_results_clusters

# --- End Cluster Comparison DGE ---
