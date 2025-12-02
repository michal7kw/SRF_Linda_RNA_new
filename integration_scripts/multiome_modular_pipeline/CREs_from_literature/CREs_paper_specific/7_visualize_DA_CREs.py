#!/usr/bin/env python3
"""
Visualize Differentially Accessible CREs with minSig/minFC Filtering

This script creates individual CRE plots for cell-type CREs that show
significant differential accessibility between conditions.

Note: This pipeline uses OVERLAPPING CRE sets (~60% overlap).
For mutually exclusive analysis, see CREs_paper_exclusive/

FILTERING CRITERIA:
- min_signal: Minimum max signal required (default: 1.0)
- min_fc: Minimum fold change required (default: 1.5)

INPUT FILES:
- output/GABA_DEG_analysis/heatmaps_deeptools/matrix_GABA_all_conditions.gz
- output/hippocampal_interneuron_CREs_genes.tsv (gene linkage)
- output/hippocampal_interneuron_CREs.bed

OUTPUT FILES:
- Individual CRE plots with gene annotations
- Metaprofiles for each comparison
- Statistics summary

Usage:
    python 7_visualize_DA_CREs.py --skip-individual        # Fast mode
    python 7_visualize_DA_CREs.py --parallel 8             # Full with parallelization
    python 7_visualize_DA_CREs.py --min-signal 2.0 --min-fc 2.0  # Stricter filtering
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import argparse
from multiprocessing import Pool
from datetime import datetime

# Set plotting defaults
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 11
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
sns.set_style("ticks")


def read_deeptools_matrix(matrix_file):
    """Read deepTools matrix file and extract signal data."""
    print(f"Reading matrix: {matrix_file}")

    matrix_tab = matrix_file.replace('.gz', '.tab')

    if not os.path.exists(matrix_tab):
        print(f"  ERROR: Matrix tab file not found: {matrix_tab}")
        return None, None, None

    with open(matrix_tab, 'r') as f:
        lines = f.readlines()

    n_regions = int(lines[0].strip().split(':')[1])
    print(f"  Regions: {n_regions}")

    params = lines[1].strip().split('\t')
    bin_size = None
    downstream = None
    upstream = None
    for param in params:
        if 'bin size:' in param:
            bin_size = int(param.split(':')[1])
        elif 'downstream:' in param:
            downstream = int(param.split(':')[1])
        elif 'upstream:' in param:
            upstream = int(param.split(':')[1])

    print(f"  Window: +/- {upstream} bp, Bin size: {bin_size} bp")

    header = lines[2].strip().split('\t')
    sample_names_all = header[1:]

    sample_names_unique = []
    prev_name = None
    bin_count = 0
    temp_bins = []

    for name in sample_names_all:
        if name != prev_name:
            if prev_name is not None:
                sample_names_unique.append(prev_name)
                temp_bins.append(bin_count)
            prev_name = name
            bin_count = 1
        else:
            bin_count += 1

    if prev_name is not None:
        sample_names_unique.append(prev_name)
        temp_bins.append(bin_count)

    n_bins = temp_bins[0]
    print(f"  Samples: {sample_names_unique}")
    print(f"  Bins per sample: {n_bins}")

    data_lines = lines[3:]
    data_matrix = []

    for line in data_lines:
        if line.strip():
            parts = line.strip().split('\t')
            data_matrix.append([float(x) for x in parts])

    data_array = np.array(data_matrix)
    print(f"  Data shape: {data_array.shape}")

    regions = [f"Region_{i+1}" for i in range(len(data_array))]

    data = {}
    for i, sample in enumerate(sample_names_unique):
        start_col = i * n_bins
        end_col = (i + 1) * n_bins
        data[sample] = data_array[:, start_col:end_col]

    bin_labels = np.linspace(-upstream, downstream, n_bins)

    return data, regions, bin_labels


def filter_significant_cres(data, regions, gene_info, sample1_name, sample2_name,
                            min_signal=2.0, min_fc=2.0):
    """Filter CREs based on signal and fold change thresholds."""
    significant_indices = []
    stats_list = []

    for i in range(len(regions)):
        sample1_signal = data[sample1_name][i, :]
        sample2_signal = data[sample2_name][i, :]

        mean1 = np.mean(sample1_signal)
        mean2 = np.mean(sample2_signal)
        max_signal = max(np.max(sample1_signal), np.max(sample2_signal))

        if mean1 > 0.01:
            fc = mean2 / mean1
        elif mean2 > 0.01:
            fc = 100.0
        else:
            fc = 1.0

        if max_signal < min_signal:
            continue

        if not (fc >= min_fc or fc <= (1.0/min_fc)):
            continue

        significant_indices.append(i)

        if i < len(gene_info):
            gene = gene_info.iloc[i].get('Genes', 'Unknown')
            cre_id = gene_info.iloc[i].get('cre_id', f'Region_{i}')
        else:
            gene = 'Unknown'
            cre_id = f'Region_{i}'

        stats_list.append({
            'index': i,
            'cre_id': cre_id,
            'gene': gene,
            'mean1': mean1,
            'mean2': mean2,
            'max_signal': max_signal,
            'fold_change': fc
        })

    return significant_indices, stats_list


def plot_comparison_metaprofile(data, bin_labels, output_file, title,
                                sample1_name, sample2_name, color1, color2,
                                n_cres=None):
    """Create metaprofile comparing two samples."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10),
                                    gridspec_kw={'height_ratios': [3, 1]})

    sample1_data = data[sample1_name]
    sample2_data = data[sample2_name]

    if n_cres is None:
        n_cres = sample1_data.shape[0]

    sample1_mean = np.mean(sample1_data, axis=0)
    sample2_mean = np.mean(sample2_data, axis=0)
    sample1_sem = stats.sem(sample1_data, axis=0)
    sample2_sem = stats.sem(sample2_data, axis=0)

    ax1.plot(bin_labels, sample1_mean, label=sample1_name, color=color1,
             linewidth=2.5, alpha=0.9)
    ax1.plot(bin_labels, sample2_mean, label=sample2_name, color=color2,
             linewidth=2.5, alpha=0.9)

    ax1.fill_between(bin_labels,
                     sample1_mean - sample1_sem,
                     sample1_mean + sample1_sem,
                     color=color1, alpha=0.2)
    ax1.fill_between(bin_labels,
                     sample2_mean - sample2_sem,
                     sample2_mean + sample2_sem,
                     color=color2, alpha=0.2)

    ax1.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax1.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Mean ATAC Signal', fontsize=12, fontweight='bold')
    ax1.set_title(f"{title}\n(n={n_cres} CREs)", fontsize=14, fontweight='bold', pad=20)
    ax1.legend(loc='best', frameon=True, fontsize=11)
    ax1.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax1)

    difference = sample2_mean - sample1_mean
    ax2.plot(bin_labels, difference, color='#6C464F', linewidth=2, alpha=0.9)
    ax2.fill_between(bin_labels, 0, difference, where=(difference > 0),
                     color=color2, alpha=0.3, label=f'{sample2_name} > {sample1_name}')
    ax2.fill_between(bin_labels, 0, difference, where=(difference < 0),
                     color=color1, alpha=0.3, label=f'{sample1_name} > {sample2_name}')
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1.5)
    ax2.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

    ax2.set_xlabel('Distance from CRE Center (bp)', fontsize=12, fontweight='bold')
    ax2.set_ylabel(f'Delta ATAC Signal', fontsize=11, fontweight='bold')
    ax2.legend(loc='best', frameon=True, fontsize=10)
    ax2.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax2)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    plt.close()


def _plot_single_cre(args):
    """Helper for parallel individual CRE plotting."""
    (i, cre_id, gene, chrom, start, end, data, bin_labels,
     sample1_name, sample2_name, color1, color2, comparison_name,
     output_dir, dpi, cre_type) = args

    fig, ax = plt.subplots(figsize=(8, 5))

    sample1_signal = data[sample1_name][i, :]
    sample2_signal = data[sample2_name][i, :]

    ax.plot(bin_labels, sample1_signal, label=sample1_name,
            color=color1, linewidth=2.5, alpha=0.9)
    ax.plot(bin_labels, sample2_signal, label=sample2_name,
            color=color2, linewidth=2.5, alpha=0.9)

    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Distance from CRE Center (bp)', fontsize=11, fontweight='bold')
    ax.set_ylabel('ATAC Signal', fontsize=11, fontweight='bold')

    gene_display = gene if len(str(gene)) < 50 else str(gene)[:47] + "..."
    ax.set_title(f'{cre_type}: {cre_id}\nGene(s): {gene_display}\n{comparison_name}',
                 fontsize=11, fontweight='bold')
    ax.legend(loc='best', frameon=True)
    ax.grid(True, alpha=0.3, linestyle=':')
    sns.despine(ax=ax)

    plt.tight_layout()
    safe_cre_id = cre_id.replace('/', '_').replace(':', '_')
    output_file = f"{output_dir}/individual_{cre_type}_{comparison_name}_{safe_cre_id}.png"
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()

    return output_file


def create_individual_plots(data, regions, bin_labels, gene_info, output_dir,
                            sample1_name, sample2_name, color1, color2,
                            comparison_name, cre_type, n_processes=1, dpi=150,
                            min_signal=2.0, min_fc=2.0):
    """Create individual plots for significant CREs."""

    sig_indices, sig_stats = filter_significant_cres(
        data, regions, gene_info, sample1_name, sample2_name,
        min_signal, min_fc
    )

    print(f"  Significant CREs: {len(sig_indices)} / {len(regions)} "
          f"(min_signal={min_signal}, min_fc={min_fc})")

    if len(sig_indices) == 0:
        print(f"  No significant CREs found for {comparison_name}")
        return [], []

    plot_args = []
    for stat in sig_stats:
        i = stat['index']
        cre_id = stat['cre_id']
        gene = stat['gene']

        if i < len(gene_info):
            chrom = gene_info.iloc[i].get('chr', 'unknown')
            start = gene_info.iloc[i].get('start', 0)
            end = gene_info.iloc[i].get('end', 0)
        else:
            chrom, start, end = 'unknown', 0, 0

        args = (i, cre_id, gene, chrom, start, end, data, bin_labels,
                sample1_name, sample2_name, color1, color2, comparison_name,
                output_dir, dpi, cre_type)
        plot_args.append(args)

    if n_processes > 1:
        print(f"  Using {n_processes} parallel processes...")
        with Pool(n_processes) as pool:
            output_files = pool.map(_plot_single_cre, plot_args)
    else:
        print(f"  Creating {len(plot_args)} plots sequentially...")
        output_files = []
        for args in plot_args:
            output_file = _plot_single_cre(args)
            output_files.append(output_file)

    return output_files, sig_stats


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Visualize differentially accessible CREs with minSig/minFC filtering'
    )

    parser.add_argument('--skip-individual', action='store_true',
                        help='Skip individual CRE plots')
    parser.add_argument('--parallel', type=int, default=1, metavar='N',
                        help='Number of parallel processes')
    parser.add_argument('--individual-dpi', type=int, default=150, metavar='N',
                        help='DPI for individual plots')
    parser.add_argument('--min-signal', type=float, default=1.0, metavar='F',
                        help='Minimum max signal required')
    parser.add_argument('--min-fc', type=float, default=1.5, metavar='F',
                        help='Minimum fold change required')

    args = parser.parse_args()

    print("="*80)
    print("VISUALIZE DIFFERENTIALLY ACCESSIBLE CREs")
    print("="*80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    if args.skip_individual:
        print("FAST MODE: Skipping individual CRE plots")
    else:
        print(f"FULL MODE: Creating individual plots")
        print(f"  Parallel: {args.parallel}, DPI: {args.individual_dpi}")
        print(f"  Min Signal: {args.min_signal}, Min FC: {args.min_fc}")
    print()

    os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_paper_specific")

    matrix_dir = "output/GABA_DEG_analysis/heatmaps_deeptools"
    suffix = f"_minSig{args.min_signal}_minFC{args.min_fc}"
    output_dir = f"output/DA_profiles{suffix}"
    os.makedirs(output_dir, exist_ok=True)

    comparisons = [
        {
            'name': 'nestin_ctrl_vs_mut',
            'sample1': 'Nestin-Ctrl',
            'sample2': 'Nestin-Mut',
            'color1': '#2E86AB',
            'color2': '#C73E1D',
            'title': 'Nestin: Ctrl vs Mut'
        },
        {
            'name': 'nestin_ctrl_vs_emx1_mut',
            'sample1': 'Nestin-Ctrl',
            'sample2': 'Emx1-Mut',
            'color1': '#2E86AB',
            'color2': '#FF8C42',
            'title': 'Nestin-Ctrl vs Emx1-Mut'
        },
        {
            'name': 'nestin_mut_vs_emx1_mut',
            'sample1': 'Nestin-Mut',
            'sample2': 'Emx1-Mut',
            'color1': '#C73E1D',
            'color2': '#FF8C42',
            'title': 'Nestin-Mut vs Emx1-Mut'
        }
    ]

    cre_types = [
        {
            'name': 'GABA',
            'matrix': f"{matrix_dir}/matrix_GABA_all_conditions.gz",
            'bed': "output/hippocampal_interneuron_CREs.bed",
            'genes_tsv': "output/hippocampal_interneuron_CREs_genes.tsv",
            'label': 'GABA (Hippocampal Interneurons)'
        },
        {
            'name': 'Excitatory',
            'matrix': f"{matrix_dir}/matrix_Excitatory_all_conditions.gz",
            'bed': "output/excitatory_neuron_CREs.bed",
            'genes_tsv': "output/excitatory_neuron_CREs_genes.tsv",
            'label': 'Excitatory'
        }
    ]

    all_stats = {}

    for cre_type in cre_types:
        print("\n" + "#"*80)
        print(f"# PROCESSING {cre_type['label'].upper()} CREs")
        print("#"*80)

        if not os.path.exists(cre_type['matrix'].replace('.gz', '.tab')):
            print(f"  WARNING: Matrix not found: {cre_type['matrix']}")
            print(f"  Please run 2b_create_heatmaps_deeptools.sh first")
            continue

        data, regions, bin_labels = read_deeptools_matrix(cre_type['matrix'])
        if data is None:
            continue

        if os.path.exists(cre_type['genes_tsv']):
            gene_info = pd.read_csv(cre_type['genes_tsv'], sep='\t')
            print(f"  Loaded gene info: {len(gene_info)} entries")
        else:
            print(f"  WARNING: Gene info not found: {cre_type['genes_tsv']}")
            print(f"  Please run 6_link_CREs_to_genes.py first")
            gene_info = pd.DataFrame()

        all_stats[cre_type['name']] = {}

        for comp in comparisons:
            print(f"\n{'='*80}")
            print(f"COMPARISON: {comp['title']}")
            print("="*80)

            if comp['sample1'] not in data or comp['sample2'] not in data:
                print(f"  WARNING: Samples not found in matrix")
                continue

            print("\nCreating metaprofile...")
            plot_comparison_metaprofile(
                data, bin_labels,
                f"{output_dir}/metaprofile_{cre_type['name']}_{comp['name']}.png",
                f"{comp['title']}\nATAC Signal at {cre_type['label']} CREs",
                comp['sample1'], comp['sample2'],
                comp['color1'], comp['color2']
            )

            if not args.skip_individual:
                print(f"\nCreating individual CRE plots...")
                output_files, sig_stats = create_individual_plots(
                    data, regions, bin_labels, gene_info, output_dir,
                    comp['sample1'], comp['sample2'],
                    comp['color1'], comp['color2'],
                    comp['name'], cre_type['name'],
                    n_processes=args.parallel,
                    dpi=args.individual_dpi,
                    min_signal=args.min_signal,
                    min_fc=args.min_fc
                )
                all_stats[cre_type['name']][comp['name']] = sig_stats

    # Save statistics
    stats_file = f"{output_dir}/DA_statistics.txt"
    with open(stats_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("DIFFERENTIALLY ACCESSIBLE CRE STATISTICS\n")
        f.write("="*80 + "\n\n")
        f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Min Signal: {args.min_signal}\n")
        f.write(f"Min FC: {args.min_fc}\n\n")

        for cre_type_name, comp_stats in all_stats.items():
            f.write(f"\n# {cre_type_name.upper()}\n")
            for comp_name, stats_list in comp_stats.items():
                f.write(f"  {comp_name}: {len(stats_list)} significant CREs\n")

    print(f"\nSaved: {stats_file}")

    print("\n" + "="*80)
    print("VISUALIZATION COMPLETE!")
    print("="*80)
    print(f"\nOutput: {output_dir}/")
    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
