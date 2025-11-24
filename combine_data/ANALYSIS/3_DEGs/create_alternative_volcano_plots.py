# %% [markdown]
# # Alternative Volcano Plot Generator
# 
# This script reads existing DEG results and creates alternative volcano plots with custom styling:
# - Upregulated genes: #FDAC5D (orange/yellow)
# - Downregulated genes: #7C3294 (purple)
# - Background grid and data point counts in regions
# - No gene labels (cleaner appearance)

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import seaborn as sns
from pathlib import Path
import re

# Set up project directories
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = os.path.join(PROJECT_DIR, "combine_data")
os.chdir(WORKING_DIR)

# Configuration
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
CUSTOM_ANALYSIS = "FC_0_25"

DEG_BY = "cell_type_L2_new"
PARENT_OUTPUT_DIR = os.path.join(INPUT_DIR, "DEGs_cell_type_L2_new")

# DEG_BY = "cell_type_L1"
# PARENT_OUTPUT_DIR = os.path.join(INPUT_DIR, "DEGs_cell_type_L1")

if CUSTOM_ANALYSIS is not None:
    PARENT_OUTPUT_DIR = PARENT_OUTPUT_DIR + CUSTOM_ANALYSIS

DGE_INPUT_DIR = os.path.join(PARENT_OUTPUT_DIR, 'dge_res')
PLOT_ALT_OUTPUT_DIR = os.path.join(PARENT_OUTPUT_DIR, 'plots_alt')

# Create alternative plots output directory
os.makedirs(PLOT_ALT_OUTPUT_DIR, exist_ok=True)

print(f"Input DGE directory: {DGE_INPUT_DIR}")
print(f"Alternative plots output directory: {PLOT_ALT_OUTPUT_DIR}")

# %%
# Custom colors for alternative volcano plots
UP_COLOR = '#FDAC5D'      # Upregulated genes (orange/yellow)
DOWN_COLOR = '#7C3294'    # Downregulated genes (purple)
NS_COLOR = '#CCCCCC'      # Non-significant genes (light gray)

# Significance thresholds
PADJ_THRESHOLD = 0.05
LOGFC_THRESHOLD = 0.25

def sanitize_filename(name):
    """Sanitizes a string to be used as a valid filename."""
    s = str(name)
    s = re.sub(r'[\\/\s]+', '_', s)
    s = re.sub(r'[^a-zA-Z0-9_-]', '', s)
    s = re.sub(r'_+', '_', s)
    s = s.strip('_')
    if not s:
        return "unnamed_group"
    return s

def plot_alternative_volcano(result_df, group_name, output_dir, filename_prefix, 
                           group1_name="Mutant", group2_name="Control", title_prefix="",
                           padj_threshold=PADJ_THRESHOLD, logfc_threshold=LOGFC_THRESHOLD,
                           up_color=UP_COLOR, down_color=DOWN_COLOR, ns_color=NS_COLOR):
    """
    Creates an alternative volcano plot with custom styling and data point counts.
    
    Args:
        result_df (pd.DataFrame): DataFrame with DGE results
        group_name (str): Name of the cell type or group
        output_dir (str): Directory to save the plot
        filename_prefix (str): Prefix for the output filename
        group1_name (str): Name of group 1 (positive logFC)
        group2_name (str): Name of group 2 (reference)
        title_prefix (str): Prefix for plot title
        padj_threshold (float): Adjusted p-value threshold
        logfc_threshold (float): Log fold change threshold
        up_color (str): Color for upregulated genes
        down_color (str): Color for downregulated genes
        ns_color (str): Color for non-significant genes
    """
    if result_df is None or result_df.empty:
        print(f"  Skipping alternative volcano plot for {title_prefix}{group_name} (no data)")
        return
    
    # Prepare data
    df = result_df.copy()
    
    # Filter out genes with very small log fold changes
    logfc_filter_tolerance = 0.05
    df = df[df['logfoldchanges'].abs() > logfc_filter_tolerance].copy()
    
    if df.empty:
        print(f"  Skipping alternative volcano plot for {title_prefix}{group_name} (no data after filtering)")
        return
    
    # Calculate -log10(p_adj)
    df['-log10p_adj'] = -np.log10(df['pvals_adj'])
    
    # Handle infinite values
    max_finite_logp = df.loc[np.isfinite(df['-log10p_adj']), '-log10p_adj'].max()
    if pd.isna(max_finite_logp) or max_finite_logp == 0:
        max_finite_logp = 300
    inf_replacement = max_finite_logp * 1.1 if max_finite_logp > 0 else 10
    df['-log10p_adj'] = df['-log10p_adj'].replace([np.inf, -np.inf], inf_replacement)
    df['-log10p_adj'] = df['-log10p_adj'].replace(np.nan, 0)
    
    # Replace NaN logfoldchanges
    df['logfoldchanges'] = df['logfoldchanges'].fillna(0)
    
    # Determine significance and assign colors
    df['color'] = ns_color  # Default: not significant
    sig_up_mask = (df['pvals_adj'] < padj_threshold) & (df['logfoldchanges'] > logfc_threshold)
    sig_down_mask = (df['pvals_adj'] < padj_threshold) & (df['logfoldchanges'] < -logfc_threshold)
    df.loc[sig_up_mask, 'color'] = up_color
    df.loc[sig_down_mask, 'color'] = down_color
    
    # Count genes in different regions
    total_genes = len(df)
    n_up = sig_up_mask.sum()
    n_down = sig_down_mask.sum()
    n_ns = total_genes - n_up - n_down
    
    # Count genes in threshold regions (for display on plot)
    # Upper left: high significance, downregulated
    upper_left = ((df['pvals_adj'] < padj_threshold) &
                  (df['logfoldchanges'] < -logfc_threshold)).sum()
    
    # Upper middle: high significance, not differentially expressed
    upper_middle = ((df['pvals_adj'] < padj_threshold) &
                    (df['logfoldchanges'].abs() <= logfc_threshold)).sum()
    
    # Upper right: high significance, upregulated
    upper_right = ((df['pvals_adj'] < padj_threshold) &
                   (df['logfoldchanges'] > logfc_threshold)).sum()
    
    # Lower left: not significant, downregulated direction
    lower_left = ((df['pvals_adj'] >= padj_threshold) &
                  (df['logfoldchanges'] < -logfc_threshold)).sum()
    
    # Lower middle: not significant, neutral
    lower_middle = ((df['pvals_adj'] >= padj_threshold) &
                    (df['logfoldchanges'].abs() <= logfc_threshold)).sum()
    
    # Lower right: not significant, upregulated direction
    lower_right = ((df['pvals_adj'] >= padj_threshold) &
                   (df['logfoldchanges'] > logfc_threshold)).sum()
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(6, 5))
    
    # Set style for clean appearance
    plt.style.use('default')
    
    # Create scatter plot
    scatter = ax.scatter(df['logfoldchanges'], df['-log10p_adj'], 
                        c=df['color'], alpha=0.6, s=12, rasterized=True, edgecolors='none')
    
    # Add threshold lines
    h_threshold_line = -np.log10(padj_threshold)
    ax.axhline(h_threshold_line, color='black', linestyle='--', lw=1, alpha=0.7)
    ax.axvline(logfc_threshold, color='black', linestyle='--', lw=1, alpha=0.7)
    ax.axvline(-logfc_threshold, color='black', linestyle='--', lw=1, alpha=0.7)
    
    # Add background grid
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # Set labels and title
    ax.set_xlabel(f'log₂(FC)', fontsize=12, fontweight='bold')
    ax.set_ylabel('-log₁₀(P value)', fontsize=12, fontweight='bold')
    
    # Create title with gene counts
    title = f'{group_name} volcano plot\nTotal genes: {total_genes:,}'
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    # Extend axis limits to provide boundaries and better positioning
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    # Extend limits by 10% on each side for better boundaries
    x_range = xlim[1] - xlim[0]
    y_range = ylim[1] - ylim[0]
    new_xlim = [xlim[0] - 0.1 * x_range, xlim[1] + 0.1 * x_range]
    new_ylim = [ylim[0] - 0.1 * y_range, ylim[1] + 0.1 * y_range]
    
    ax.set_xlim(new_xlim[0], new_xlim[1])
    ax.set_ylim(new_ylim[0], new_ylim[1])
    
    # Update limits for positioning
    xlim = new_xlim
    ylim = new_ylim
    
    # Position count annotations in the six corner/edge positions: top left, top center, top right, bottom left, bottom center, bottom right
    # Calculate positions for the six regions
    edge_offset_y = 0.02 * (ylim[1] - ylim[0])  # Small offset from top/bottom edges
    edge_offset_x = 0.02 * (xlim[1] - xlim[0])  # Small offset from left/right edges
    
    top_y = ylim[1] - edge_offset_y
    bottom_y = ylim[0] + edge_offset_y
    left_x = xlim[0] + edge_offset_x
    right_x = xlim[1] - edge_offset_x
    center_x = 0
    
    # Top left corner (upper left region)
    if upper_left > 0:
        ax.text(left_x, top_y,
                str(upper_left),
                bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.9, edgecolor='black'),
                fontsize=10, ha='left', va='top')
    
    # Top center (upper middle region)
    if upper_middle > 0:
        ax.text(center_x, top_y,
                str(upper_middle),
                bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.9, edgecolor='black'),
                fontsize=10, ha='center', va='top')
    
    # Top right corner (upper right region)
    if upper_right > 0:
        ax.text(right_x, top_y,
                str(upper_right),
                bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.9, edgecolor='black'),
                fontsize=10, ha='right', va='top')
    
    # Bottom left corner (lower left region)
    if lower_left > 0:
        ax.text(left_x, bottom_y,
                str(lower_left),
                bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.9, edgecolor='black'),
                fontsize=10, ha='left', va='bottom')
    
    # Bottom center (lower middle region)
    if lower_middle > 0:
        ax.text(center_x, bottom_y,
                str(lower_middle),
                bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.9, edgecolor='black'),
                fontsize=10, ha='center', va='bottom')
    
    # Bottom right corner (lower right region)
    if lower_right > 0:
        ax.text(right_x, bottom_y,
                str(lower_right),
                bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.9, edgecolor='black'),
                fontsize=10, ha='right', va='bottom')
    
    # Style the plot - show all boundaries/frames
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    
    # Set all spine linewidths for visible boundaries
    ax.spines['top'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    
    # Set spine colors to black for clear boundaries
    ax.spines['top'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    
    # Set tick parameters
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    plt.tight_layout()
    
    # Save the plot
    sanitized_group_name = sanitize_filename(group_name)
    output_filename = f"{filename_prefix}_{sanitized_group_name}_volcano_alt.png"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_filename)
    
    try:
        fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  Alternative volcano plot saved: {output_filename}")
    except Exception as e:
        print(f"  Error saving alternative volcano plot for {group_name}: {e}")
    
    plt.close(fig)

# %%
def process_dge_results_directory(dge_subdir, plot_subdir, analysis_type):
    """
    Process all CSV files in a DGE results subdirectory and create alternative volcano plots.
    
    Args:
        dge_subdir (str): Path to DGE results subdirectory
        plot_subdir (str): Path to plots output subdirectory  
        analysis_type (str): Type of analysis for labeling
    """
    print(f"\nProcessing {analysis_type} results...")
    
    # Find all CSV files in the directory
    csv_pattern = os.path.join(dge_subdir, "*.csv")
    csv_files = glob.glob(csv_pattern)
    
    # Filter out skipped files
    csv_files = [f for f in csv_files if 'skipped' not in os.path.basename(f)]
    
    if not csv_files:
        print(f"  No CSV files found in {dge_subdir}")
        return
    
    print(f"  Found {len(csv_files)} CSV files to process")
    
    for csv_file in csv_files:
        try:
            # Read the CSV file
            df = pd.read_csv(csv_file)
            
            if df.empty:
                print(f"  Skipping empty file: {os.path.basename(csv_file)}")
                continue
            
            # Extract group name from filename
            filename = os.path.basename(csv_file)
            
            # Parse different filename patterns
            if analysis_type == "both_geno_cond_comp":
                # Pattern: dge_both_genomes_GROUPNAME_mut_vs_ctrl.csv
                match = re.search(r'dge_both_genomes_(.+)_mut_vs_ctrl\.csv', filename)
                if match:
                    group_name = match.group(1)
                    title_prefix = "Both genomes - "
                    filename_prefix = "dge_conditions_comparison"
                else:
                    continue
                    
            elif analysis_type == "geno_spec_cond_comp":
                # Pattern: dge_GENOTYPE_GROUPNAME_mut_vs_ctrl.csv
                match = re.search(r'dge_(Emx1|Nestin)_(.+)_mut_vs_ctrl\.csv', filename)
                if match:
                    genotype = match.group(1)
                    group_name = match.group(2)
                    title_prefix = f"{genotype} - "
                    filename_prefix = f"dge_{genotype}"
                else:
                    continue
                    
            elif analysis_type == "cond_spec_geno_comp":
                # Pattern: dge_CONDITION_cond_GROUPNAME_Nestin_vs_Emx1.csv
                match = re.search(r'dge_(Control|Mutant)_cond_(.+)_Nestin_vs_Emx1\.csv', filename)
                if match:
                    condition = match.group(1)
                    group_name = match.group(2)
                    title_prefix = f"{condition} Cond - "
                    filename_prefix = f"dge_{condition}_cond_Nestin_vs_Emx1"
                else:
                    continue
            else:
                print(f"  Unknown analysis type: {analysis_type}")
                continue
            
            # Create alternative volcano plot
            plot_alternative_volcano(
                result_df=df,
                group_name=group_name,
                output_dir=plot_subdir,
                filename_prefix=filename_prefix,
                title_prefix=title_prefix
            )
            
        except Exception as e:
            print(f"  Error processing {os.path.basename(csv_file)}: {e}")

# %%
# Process all three types of DGE analyses
analysis_types = [
    ("both_geno_cond_comp", "Both Genomes Condition Comparison"),
    ("geno_spec_cond_comp", "Genotype-Specific Condition Comparison"), 
    ("cond_spec_geno_comp", "Condition-Specific Genotype Comparison")
]

for analysis_dir, analysis_name in analysis_types:
    dge_subdir = os.path.join(DGE_INPUT_DIR, analysis_dir)
    plot_subdir = os.path.join(PLOT_ALT_OUTPUT_DIR, analysis_dir)
    
    if os.path.exists(dge_subdir):
        process_dge_results_directory(dge_subdir, plot_subdir, analysis_dir)
    else:
        print(f"\nSkipping {analysis_name}: Directory not found - {dge_subdir}")

print("\n" + "="*60)
print("Alternative volcano plot generation completed!")
print(f"Plots saved in: {PLOT_ALT_OUTPUT_DIR}")
print("="*60)