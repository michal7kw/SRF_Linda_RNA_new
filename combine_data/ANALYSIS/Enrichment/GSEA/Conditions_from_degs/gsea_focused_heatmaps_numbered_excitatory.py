#!/usr/bin/env python3
"""
Complete GSEA Pipeline for Excitatory Neurons (cell_type_L1 level)

This script performs a comprehensive analysis pipeline:
1. Loads scRNA-seq data (annotation_final.h5ad)
2. Subsets by cell_type_L1 == 'Excitatory' (combines ALL excitatory subtypes)
3. Performs differential gene expression analysis (Mutant vs Control)
4. Runs focused GSEA on DEGs using custom pathway keywords
5. Creates numbered heatmaps for visualization

This is a true L1-level analysis that treats all excitatory neurons as a single group,
rather than analyzing individual L2 subtypes separately.

Output:
- GSEA results for each genotype (Emx1, Nestin)
- Focused pathway results based on keyword filtering
- Publication-ready numbered heatmaps with separate legend tables
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
import os
import sys
import pathlib
import gseapy as gp
import multiprocessing
from functools import partial
import traceback
import json
import anndata as ad
import re
import warnings
from matplotlib import rcParams

warnings.filterwarnings('ignore')

# --- Configuration ---
PROJECT_DIR = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)

ORGANISM = 'Mouse'

# Set up directories
BASE_RESULTS_DIR = os.path.join(WORKING_DIR, "results_from_raw")
INPUT_DIR = BASE_RESULTS_DIR
ADATA_PATH = os.path.join(INPUT_DIR, "annotation_final.h5ad")

# Output directories
GSEA_BASE_DIR = os.path.join(BASE_RESULTS_DIR, 'gsea_analysis_excitatory_L1')
OUTPUT_DIR = os.path.join(BASE_RESULTS_DIR, 'gsea_focused_heatmaps_excitatory_L1')

pathlib.Path(GSEA_BASE_DIR).mkdir(exist_ok=True)
pathlib.Path(OUTPUT_DIR).mkdir(exist_ok=True)

# GSEA Configuration
GSEA_GENE_SETS = ['KEGG_2019_Mouse', 'GO_Biological_Process_2023', 'Reactome_2022', 'MSigDB_Hallmark_2020']
GSEA_PADJ_THRESHOLD = 0.25
MIN_CELLS = 10
N_TOP_TERMS_PLOT = 15

# Cell type filtering
CELL_TYPE_L1_FILTER = 'Excitatory'

# Genotypes to analyze separately
GENOTYPES = {
    'Emx1': 'Emx1',
    'Nestin': 'Nestin'
}

# Publication settings for heatmaps
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
rcParams['font.size'] = 6.5
rcParams['axes.labelsize'] = 7
rcParams['axes.titlesize'] = 8
rcParams['xtick.labelsize'] = 6.5
rcParams['ytick.labelsize'] = 6.5
rcParams['legend.fontsize'] = 5.5
rcParams['figure.titlesize'] = 8
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

# Define focus configurations for excitatory neurons
FOCUS_CONFIG = {
    'cell_death': {
        'keywords': [
            'neurodegeneration', 'neuronal death', 'axonopathy', 'neuron death',
            'apoptosis', 'apoptotic', 'necroptosis', 'pyroptosis', 'ferroptosis',
            'autophagy', 'autophagic', 'parthanatos', 'cell death'
        ],
        'file_suffix': 'Neuro_CellDeath',
        'title_name': 'Cell Death (Excitatory L1)'
    },
    'synaptic': {
        'keywords': [
            'synapse', 'synaptic', 'neurotransmitter', 'glutamate',
            'synaptic vesicle', 'synaptic transmission', 'synaptic plasticity',
            'excitatory postsynaptic', 'epsp', 'glutamatergic',
            'ampa', 'nmda', 'kainate', 'receptor',
            'long term potentiation', 'ltp', 'long term depression', 'ltd',
            'dendritic spine', 'postsynaptic density'
        ],
        'file_suffix': 'Synaptic_Function',
        'title_name': 'Synaptic Function (Excitatory L1)'
    },
    'development': {
        'keywords': [
            'neurogenesis', 'neuronal differentiation', 'neuron development',
            'axon', 'dendrite', 'axonogenesis', 'dendritogenesis',
            'neurite', 'growth cone', 'axon guidance',
            'neuronal migration', 'neuronal maturation',
            'neuron projection', 'neuron morphogenesis'
        ],
        'file_suffix': 'Neuro_Development',
        'title_name': 'Neuronal Development (Excitatory L1)'
    },
    'metabolism': {
        'keywords': [
            'mitochondr', 'oxidative phosphorylation', 'respiratory chain',
            'electron transport', 'atp synthase', 'citric acid cycle',
            'glycolysis', 'glucose metabolism', 'energy metabolism',
            'metabolic', 'oxidative stress', 'reactive oxygen species'
        ],
        'file_suffix': 'Metabolism_Energy',
        'title_name': 'Metabolism & Energy (Excitatory L1)'
    }
}


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

        # Compatibility fix for older gseapy versions
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

        # Generate plots for significant terms
        if 'fdr' in results_df.columns and 'nes' in results_df.columns:
            significant_results = results_df[results_df['fdr'] < gsea_padj_threshold]
            if not significant_results.empty:
                print(f"Found {len(significant_results)} significant terms in {gene_set}.")

        print(f"--- Finished GSEA for gene set: {gene_set} ---")
        return results_df

    except Exception as e:
        print(f"An error occurred during GSEA for gene set '{gene_set}': {e}")
        traceback.print_exc()
        return None


def filter_gsea_by_keywords(gsea_results_df, keywords, focus_name):
    """
    Filter GSEA results based on keyword matching.

    Parameters:
    -----------
    gsea_results_df : pd.DataFrame
        GSEA results with 'Term', 'nes', 'fdr' columns
    keywords : list
        List of keywords to search for
    focus_name : str
        Name of the focus area

    Returns:
    --------
    pd.DataFrame
        Filtered results matching keywords
    """
    if gsea_results_df is None or gsea_results_df.empty:
        return pd.DataFrame()

    term_col = 'Term' if 'Term' in gsea_results_df.columns else 'term'

    # Filter by keywords
    mask = pd.Series(False, index=gsea_results_df.index)
    for keyword in keywords:
        mask |= gsea_results_df[term_col].str.lower().str.contains(keyword.lower(), na=False)

    filtered_df = gsea_results_df[mask].copy()
    print(f"  Filtered {len(filtered_df)} terms for '{focus_name}' from {len(gsea_results_df)} total terms")

    return filtered_df


def perform_deg_and_gsea_for_genotype(adata_excitatory, genotype_name, genotype_label):
    """
    Perform DEG analysis and GSEA for a specific genotype.

    Parameters:
    -----------
    adata_excitatory : AnnData
        AnnData object with excitatory neurons only
    genotype_name : str
        Name for directory (e.g., 'Emx1', 'Nestin')
    genotype_label : str
        Label in the data (e.g., 'Emx1', 'Nestin')

    Returns:
    --------
    dict
        Dictionary with GSEA results for each gene set
    """
    print(f"\n{'='*70}")
    print(f"Processing Genotype: {genotype_name}")
    print(f"{'='*70}")

    # Filter by genotype
    adata_geno = adata_excitatory[adata_excitatory.obs['genotype'] == genotype_label].copy()
    print(f"Cells for {genotype_name}: {adata_geno.n_obs}")

    if adata_geno.n_obs == 0:
        print(f"No cells found for genotype '{genotype_name}'. Skipping.")
        return None

    # Check condition distribution
    condition_counts = adata_geno.obs['condition'].value_counts()
    print(f"  - Control: {condition_counts.get('Control', 0)}")
    print(f"  - Mutant: {condition_counts.get('Mutant', 0)}")

    # Create output directory
    genotype_dir = pathlib.Path(GSEA_BASE_DIR) / genotype_name.lower()
    comparison_dir = genotype_dir / "Excitatory_L1_Mutant_vs_Control"
    comparison_dir.mkdir(parents=True, exist_ok=True)

    # Filter genes
    print(f"\nFiltering genes...")
    n_genes_before = adata_geno.n_vars
    sc.pp.filter_genes(adata_geno, min_cells=MIN_CELLS)
    n_genes_after = adata_geno.n_vars
    print(f"  - Genes before: {n_genes_before}")
    print(f"  - Genes after (min_cells={MIN_CELLS}): {n_genes_after}")

    # Save cluster sizes
    cluster_sizes = {
        'cluster_name': 'Excitatory_L1',
        'genotype': genotype_name,
        'mutant_count': int(condition_counts.get('Mutant', 0)),
        'control_count': int(condition_counts.get('Control', 0))
    }
    with open(comparison_dir / "cluster_sizes.json", 'w') as f:
        json.dump(cluster_sizes, f, indent=4)

    # Prepare data for DEG from raw counts
    print("\nPreparing data for DEG analysis from raw counts...")
    use_layer = None
    if hasattr(adata_geno, 'raw') and adata_geno.raw is not None and adata_geno.raw.X is not None:
        print("Found .raw attribute. Creating normalized layer for DEG analysis.")

        adata_for_dge = ad.AnnData(adata_geno.raw.X.copy())
        adata_for_dge.obs_names = adata_geno.obs_names
        adata_for_dge.var_names = adata_geno.raw.var_names

        sc.pp.normalize_total(adata_for_dge, target_sum=1e4)
        sc.pp.log1p(adata_for_dge)

        adata_for_dge = adata_for_dge[:, adata_geno.var_names].copy()
        adata_geno.layers['for_DEGs'] = adata_for_dge.X.copy()
        use_layer = 'for_DEGs'
        print("'for_DEGs' layer created.")
    else:
        print("Warning: adata_geno.raw.X not found. Using adata_geno.X for DEG analysis.")

    # Differential Expression Analysis
    print(f"\nRunning differential expression for Excitatory L1 neurons ({genotype_name})...")
    sc.tl.rank_genes_groups(
        adata_geno,
        groupby='condition',
        groups=['Mutant'],
        reference='Control',
        method='t-test',
        use_raw=False,
        layer=use_layer
    )

    # Prepare ranked gene list
    print("Preparing ranked gene list for GSEA...")
    de_results = sc.get.rank_genes_groups_df(adata_geno, group='Mutant')
    de_results = de_results.copy()
    de_results['pvals_adj'] = de_results['pvals_adj'].replace(0, 1e-300)
    de_results['logfoldchanges'] = de_results['logfoldchanges'].fillna(0)
    de_results['rank_metric'] = np.sign(de_results['logfoldchanges']) * (-np.log10(de_results['pvals_adj']))
    rank_df = de_results.set_index('names')['rank_metric'].sort_values(ascending=False)

    if rank_df.index.duplicated().any():
        rank_df = rank_df[~rank_df.index.duplicated()]

    print(f"Created ranked list with {len(rank_df)} genes.")

    # Save DEG results
    de_results.to_csv(comparison_dir / "DEG_results.csv", index=False)
    rank_df.to_csv(comparison_dir / "ranked_genes.csv")
    print(f"Saved DEG results to {comparison_dir}")

    # Run GSEA in parallel
    print("\n--- Starting Parallel GSEA Analysis ---")

    worker_func = partial(run_gsea_for_gene_set,
                          rank_df=rank_df,
                          organism=ORGANISM,
                          cluster_output_dir=comparison_dir,
                          adata_var_names=list(adata_excitatory.var_names),
                          gsea_padj_threshold=GSEA_PADJ_THRESHOLD,
                          n_top_terms_plot=N_TOP_TERMS_PLOT)

    num_processes = min(len(GSEA_GENE_SETS), os.cpu_count() - 1 if os.cpu_count() > 1 else 1)
    print(f"Using {num_processes} processes for {len(GSEA_GENE_SETS)} gene sets.")

    with multiprocessing.Pool(processes=num_processes) as pool:
        all_results_dfs = pool.map(worker_func, GSEA_GENE_SETS)

    print("\n--- Parallel GSEA complete. Processing results. ---")

    # Process and save results
    summary_results = []
    full_summary_results = []
    gsea_results_dict = {}

    for i, df in enumerate(all_results_dfs):
        if df is not None:
            gene_set = GSEA_GENE_SETS[i]
            df['Gene_Set'] = gene_set
            df['Direction'] = df['nes'].apply(lambda x: 'Up' if x > 0 else 'Down')
            full_summary_results.append(df.copy())
            gsea_results_dict[gene_set] = df.copy()

            if 'fdr' in df.columns:
                sig_results = df[df['fdr'] < GSEA_PADJ_THRESHOLD].copy()
                if not sig_results.empty:
                    summary_results.append(sig_results)

    # Save summary files
    if summary_results:
        all_sig_results = pd.concat(summary_results, ignore_index=True)
        all_sig_results.to_csv(comparison_dir / "GSEA_Summary_All_Genesets.csv", index=False)
        print(f"\nSaved summary with {len(all_sig_results)} significant terms")

    if full_summary_results:
        all_full_results = pd.concat(full_summary_results, ignore_index=True)
        all_full_results.to_csv(comparison_dir / "GSEA_Full_Report_All_Genesets.csv", index=False)
        print(f"Saved full report with {len(all_full_results)} total terms")

    # Apply focused filtering for each focus area
    print(f"\n--- Applying Focused Filtering for {genotype_name} ---")
    for focus_key, focus_config in FOCUS_CONFIG.items():
        print(f"\nProcessing focus area: {focus_key}")

        focus_results = []
        for gene_set, df in gsea_results_dict.items():
            filtered_df = filter_gsea_by_keywords(df, focus_config['keywords'], focus_key)
            if not filtered_df.empty:
                focus_results.append(filtered_df)

        if focus_results:
            focus_combined = pd.concat(focus_results, ignore_index=True)
            focus_file = comparison_dir / f"GSEA_Focus_{focus_config['file_suffix']}.csv"
            focus_combined.to_csv(focus_file, index=False)
            print(f"  Saved {len(focus_combined)} focused terms to {focus_file}")

    return gsea_results_dict


def shorten_term_name(term, max_length=60):
    """Intelligently shorten pathway/GO term names for legend table."""
    prefixes_to_remove = [
        'GOBP_', 'GOCC_', 'GOMF_', 'KEGG_', 'REACTOME_', 'HALLMARK_',
        'GO_', 'WP_', 'BIOCARTA_', 'PID_'
    ]

    for prefix in prefixes_to_remove:
        if term.startswith(prefix):
            term = term[len(prefix):]
            break

    term = term.replace('_', ' ').title()

    abbreviations = {
        'Regulation Of': 'Reg.',
        'Positive Regulation': '+Reg.',
        'Negative Regulation': '-Reg.',
        'Programmed Cell Death': 'PCD',
        'Processing': 'Proc.',
        'Independent': 'Indep.',
        'Dependent': 'Dep.',
        'Mediated': 'Med.',
        'Activation': 'Activ.',
        'Recognition': 'Recog.',
        'Mitochondrial': 'Mito.',
    }

    for full, abbrev in sorted(abbreviations.items(), key=lambda x: len(x[0]), reverse=True):
        term = term.replace(full, abbrev)

    term = ' '.join(term.split())

    if len(term) > max_length:
        term = term[:max_length-3] + '...'

    return term


def create_term_legend_table(ax, term_dict, fontsize=6):
    """Create a clean table showing term numbers and their names."""
    ax.axis('off')

    table_data = []
    for num, term in term_dict.items():
        table_data.append([f"{num}. {term}"])

    table = ax.table(
        cellText=table_data,
        cellLoc='left',
        loc='center',
        bbox=[0, 0, 1, 1]
    )

    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)

    row_height = 1.0 / len(table_data)

    for i in range(len(table_data)):
        cell = table[(i, 0)]
        if i % 2 == 0:
            cell.set_facecolor('#F8F8F8')
        else:
            cell.set_facecolor('white')
        cell.set_text_props(fontsize=fontsize)
        cell.set_edgecolor('lightgray')
        cell.set_linewidth(0.3)
        cell.set_height(row_height)

    return table


def create_numbered_heatmap_from_files(focus_key, max_terms=12, fdr_threshold=0.25, no_fdr=False):
    """
    Create numbered heatmap by loading focused GSEA results from files.

    Parameters:
    -----------
    focus_key : str
        Focus area key
    max_terms : int
        Maximum number of terms
    fdr_threshold : float
        FDR threshold
    no_fdr : bool
        If True, create version without FDR stars
    """
    focus_config = FOCUS_CONFIG[focus_key]
    title_name = focus_config['title_name']
    file_suffix = focus_config['file_suffix']

    print(f"\n{'='*70}")
    print(f"Creating numbered heatmap for {title_name}{' (no FDR)' if no_fdr else ''}")
    print(f"{'='*70}")

    # Load data from both genotypes
    all_results = []

    for genotype_name in GENOTYPES.keys():
        genotype_dir = genotype_name.lower()
        file_path = pathlib.Path(GSEA_BASE_DIR) / genotype_dir / "Excitatory_L1_Mutant_vs_Control" / f"GSEA_Focus_{file_suffix}.csv"

        if not file_path.exists():
            print(f"Warning: File not found for {genotype_name}: {file_path}")
            continue

        print(f"Loading results from {file_path}")
        df = pd.read_csv(file_path)

        if 'Term' not in df.columns or 'nes' not in df.columns or 'fdr' not in df.columns:
            print(f"Warning: Required columns not found in {file_path}")
            continue

        df_renamed = df[['Term', 'nes', 'fdr']].copy()
        df_renamed.rename(columns={
            'nes': f"{genotype_name}_nes",
            'fdr': f"{genotype_name}_fdr"
        }, inplace=True)

        # Handle duplicate terms by keeping the one with the lowest FDR
        df_renamed = df_renamed.sort_values(f"{genotype_name}_fdr")
        df_renamed = df_renamed[~df_renamed['Term'].duplicated(keep='first')]

        df_renamed.set_index('Term', inplace=True)
        all_results.append(df_renamed)

    if len(all_results) != len(GENOTYPES):
        print(f"\nError: Could not load results for all genotypes")
        return

    merged_df = pd.concat(all_results, axis=1, join='outer')

    if merged_df.empty:
        print(f"No data found for {title_name}. Exiting.")
        return

    print(f"\nFound {len(merged_df)} total terms across both genotypes.")

    # Filter and select terms
    nes_cols = [f"{name}_nes" for name in GENOTYPES.keys()]
    fdr_cols = [f"{name}_fdr" for name in GENOTYPES.keys()]

    merged_df['min_fdr'] = merged_df[fdr_cols].min(axis=1)
    merged_df['mean_abs_nes'] = merged_df[nes_cols].abs().mean(axis=1)

    significant_df = merged_df[merged_df['min_fdr'] < fdr_threshold].copy()

    if len(significant_df) == 0:
        print(f"  Warning: No terms pass FDR threshold {fdr_threshold}. Using all terms.")
        significant_df = merged_df.copy()

    top_terms = significant_df.nlargest(max_terms, 'mean_abs_nes')

    print(f"Selected top {len(top_terms)} terms by mean absolute NES")

    # Create term number mapping
    term_dict = {}
    numbered_terms = []
    for i, term in enumerate(top_terms.index, 1):
        term_dict[i] = shorten_term_name(term)
        numbered_terms.append(str(i))

    # Prepare plot data
    plot_df = top_terms[nes_cols].copy()
    plot_df.columns = GENOTYPES.keys()
    plot_df.index = numbered_terms

    # Calculate dimensions
    n_terms = len(plot_df)
    n_cols = len(plot_df.columns)

    heatmap_width = min(1.8 + (n_cols * 0.5), 3.5)
    heatmap_height = min(1.5 + (n_terms * 0.2), 6)
    legend_width = 4.5
    total_width = heatmap_width + legend_width + 0.5

    # Create figure
    fig = plt.figure(figsize=(total_width, heatmap_height))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[heatmap_width, legend_width], wspace=0.2)

    ax_heatmap = fig.add_subplot(gs[0])
    ax_legend = fig.add_subplot(gs[1])

    # Color scale
    vmax = np.nanmax(np.abs(plot_df.values))
    vmax = max(vmax, 0.1)
    vmin = -vmax

    # Create significance annotations
    annot_df = pd.DataFrame(index=plot_df.index, columns=plot_df.columns, dtype=str)

    if not no_fdr:
        for idx, (num, orig_term) in enumerate(zip(numbered_terms, top_terms.index)):
            for genotype_name in GENOTYPES.keys():
                fdr_col = f"{genotype_name}_fdr"

                if orig_term in top_terms.index:
                    fdr_val = top_terms.loc[orig_term, fdr_col]

                    if pd.notna(fdr_val):
                        if fdr_val < 0.001:
                            sig = '***'
                        elif fdr_val < 0.01:
                            sig = '**'
                        elif fdr_val < fdr_threshold:
                            sig = '*'
                        else:
                            sig = ''
                        annot_df.loc[num, genotype_name] = sig
                    else:
                        annot_df.loc[num, genotype_name] = ''
    else:
        annot_df = annot_df.fillna('')

    # Create heatmap
    sns.heatmap(
        plot_df,
        cmap='RdBu_r',
        vmin=vmin,
        vmax=vmax,
        center=0,
        annot=annot_df,
        fmt='',
        linewidths=0.3,
        linecolor='lightgray',
        cbar_kws={
            'label': 'NES',
            'shrink': 0.7,
            'aspect': 12,
            'pad': 0.02
        },
        ax=ax_heatmap
    )

    # Format heatmap
    ax_heatmap.set_xlabel('Genotype', fontsize=7, fontweight='bold')
    ax_heatmap.set_ylabel('Pathway #', fontsize=7, fontweight='bold')
    ax_heatmap.set_title(f'{title_name}\n(Mutant vs Control)',
                         fontsize=8, fontweight='bold', pad=8)

    ax_heatmap.tick_params(axis='x', rotation=45, labelsize=6.5)
    ax_heatmap.tick_params(axis='y', rotation=0, labelsize=6.5)

    # Add significance legend
    if not no_fdr:
        legend_elements = [
            mpatches.Patch(facecolor='white', edgecolor='black', label='* FDR<0.25'),
            mpatches.Patch(facecolor='white', edgecolor='black', label='** FDR<0.01'),
            mpatches.Patch(facecolor='white', edgecolor='black', label='*** FDR<0.001')
        ]

        ax_heatmap.legend(
            handles=legend_elements,
            bbox_to_anchor=(0.5, -0.20),
            loc='upper center',
            frameon=True,
            fontsize=5.5,
            ncol=3,
            handlelength=0.8,
            handletextpad=0.4,
            columnspacing=0.8
        )

    # Create legend table
    create_term_legend_table(ax_legend, term_dict, fontsize=6)

    # Adjust layout
    if not no_fdr:
        plt.tight_layout(rect=[0, 0.08, 1, 1])
    else:
        plt.tight_layout()

    # Save outputs
    fdr_suffix = '_no_fdr' if no_fdr else ''
    base_filename = f"GSEA_Focused_Excitatory_L1_{focus_key}_numbered{fdr_suffix}"

    png_path = os.path.join(OUTPUT_DIR, f"{base_filename}.png")
    plt.savefig(png_path, dpi=600, bbox_inches='tight', facecolor='white', pad_inches=0.15)
    print(f"\n✓ Saved PNG (600 DPI): {png_path}")

    pdf_path = os.path.join(OUTPUT_DIR, f"{base_filename}.pdf")
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white', pad_inches=0.15)
    print(f"✓ Saved PDF (vector): {pdf_path}")

    # Save legend CSV
    term_df = pd.DataFrame(list(term_dict.items()), columns=['Number', 'Pathway'])
    csv_path = os.path.join(OUTPUT_DIR, f"{base_filename}_legend.csv")
    term_df.to_csv(csv_path, index=False)
    print(f"✓ Saved legend CSV: {csv_path}")

    plt.close()


def main():
    """Main function to run the complete pipeline."""
    print("=" * 70)
    print("Complete GSEA Pipeline for Excitatory Neurons (L1 Level)")
    print("=" * 70)
    print(f"GSEA output: {GSEA_BASE_DIR}")
    print(f"Heatmap output: {OUTPUT_DIR}")

    # Load data
    if not os.path.exists(ADATA_PATH):
        print(f"Error: Data file not found at {ADATA_PATH}")
        sys.exit(1)

    print(f"\nLoading scRNA data from: {ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)
    print(f"Loaded data: {adata.n_obs} cells x {adata.n_vars} genes")

    # Check for required columns
    if 'cell_type_L1' not in adata.obs.columns:
        print(f"Error: 'cell_type_L1' column not found")
        sys.exit(1)

    print(f"\nAvailable cell_type_L1 categories: {list(adata.obs['cell_type_L1'].unique())}")

    # Filter for excitatory neurons
    print(f"\nFiltering for excitatory neurons (cell_type_L1 == '{CELL_TYPE_L1_FILTER}')...")
    adata_excitatory = adata[adata.obs['cell_type_L1'] == CELL_TYPE_L1_FILTER].copy()
    print(f"Filtered data: {adata_excitatory.n_obs} excitatory neurons ({adata_excitatory.n_obs/adata.n_obs*100:.1f}% of total)")

    if adata_excitatory.n_obs == 0:
        print(f"Error: No excitatory neurons found")
        sys.exit(1)

    # Show genotype and condition distribution
    print(f"\nGenotype distribution:")
    for geno, count in adata_excitatory.obs['genotype'].value_counts().items():
        print(f"  - {geno}: {count} cells")

    print(f"\nCondition distribution:")
    for cond, count in adata_excitatory.obs['condition'].value_counts().items():
        print(f"  - {cond}: {count} cells")

    # Perform DEG and GSEA for each genotype
    for genotype_name, genotype_label in GENOTYPES.items():
        perform_deg_and_gsea_for_genotype(adata_excitatory, genotype_name, genotype_label)

    # Create numbered heatmaps for each focus area
    print("\n" + "=" * 70)
    print("Creating Numbered Heatmaps")
    print("=" * 70)

    for focus_key in FOCUS_CONFIG.keys():
        # Version with FDR
        create_numbered_heatmap_from_files(
            focus_key,
            max_terms=12,
            fdr_threshold=0.25,
            no_fdr=False
        )
        # Version without FDR
        create_numbered_heatmap_from_files(
            focus_key,
            max_terms=12,
            fdr_threshold=0.25,
            no_fdr=True
        )

    print("\n" + "=" * 70)
    print("✓ Complete pipeline finished successfully!")
    print(f"✓ GSEA results: {GSEA_BASE_DIR}")
    print(f"✓ Heatmaps: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == '__main__':
    # Required for multiprocessing
    multiprocessing.freeze_support()
    main()
