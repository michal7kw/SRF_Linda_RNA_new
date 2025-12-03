#!/usr/bin/env python3
"""
8_differential_motif_analysis.py

PURPOSE:
Compare TFBS motif enrichment between experimental conditions to identify
transcription factors differentially enriched in mutant vs control samples
at promoter-like ATAC peaks.

DESCRIPTION:
This script parses HOMER output files and performs differential enrichment
analysis to identify:
1. TF motifs enriched in Nestin-Mut vs Nestin-Ctrl (mutation effect in Nestin)
2. TF motifs enriched in Emx1-Mut vs Nestin-Ctrl (mutation effect in Emx1)
3. TF motifs differentially enriched between genotypes (Nestin-Mut vs Emx1-Mut)

The analysis helps identify which transcription factors show altered binding
at promoters of splicing genes due to SRF mutation.

INPUT:
- HOMER knownResults.txt files from Step 7

OUTPUT:
- Differential enrichment tables (TSV)
- TF motif comparison matrices
- Summary statistics

AUTHOR: Claude Code
DATE: December 2024
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import logging

# =============================================================================
# CONFIGURATION
# =============================================================================

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Default paths
SCRIPT_DIR = Path(__file__).parent.resolve()
HOMER_DIR = SCRIPT_DIR / "output" / "tfbs_analysis" / "homer_results"
OUTPUT_DIR = SCRIPT_DIR / "output" / "tfbs_analysis" / "differential_analysis"

# Sample names
SAMPLES = ["Nestin_Ctrl", "Nestin_Mut", "Emx1_Mut"]

# Statistical thresholds
P_VALUE_THRESHOLD = 0.05
LOG_PVALUE_THRESHOLD = 2.0  # -log10(0.01)
FOLD_CHANGE_THRESHOLD = 1.5


# =============================================================================
# FUNCTIONS
# =============================================================================

def parse_homer_results(filepath: Path) -> pd.DataFrame:
    """
    Parse HOMER knownResults.txt file into a pandas DataFrame.

    The HOMER output contains:
    - Motif Name
    - Consensus
    - P-value
    - Log P-value
    - q-value (FDR)
    - # of Target Sequences with Motif
    - % of Target Sequences with Motif
    - # of Background Sequences with Motif
    - % of Background Sequences with Motif
    """
    logger.info(f"Parsing HOMER results: {filepath}")

    if not filepath.exists():
        logger.error(f"File not found: {filepath}")
        return None

    try:
        # HOMER output is tab-separated with header
        df = pd.read_csv(filepath, sep='\t', header=0)

        # Standardize column names
        df.columns = [
            'Motif_Name', 'Consensus', 'P_value', 'Log_P_value', 'Q_value',
            'Num_Target', 'Pct_Target', 'Num_Background', 'Pct_Background'
        ]

        # Clean percentage columns (remove % sign if present)
        for col in ['Pct_Target', 'Pct_Background']:
            if df[col].dtype == object:
                df[col] = df[col].str.rstrip('%').astype(float)

        # Extract TF name from motif name (usually in format "TF_Name(Family)/Source")
        df['TF_Name'] = df['Motif_Name'].apply(lambda x: x.split('(')[0].strip())
        df['TF_Family'] = df['Motif_Name'].apply(
            lambda x: x.split('(')[1].split(')')[0] if '(' in x else 'Unknown'
        )

        logger.info(f"  Loaded {len(df)} motifs")
        return df

    except Exception as e:
        logger.error(f"Error parsing {filepath}: {e}")
        return None


def calculate_enrichment_score(pct_target: float, pct_background: float) -> float:
    """
    Calculate fold enrichment of motif in target vs background.
    """
    if pct_background == 0:
        return float('inf') if pct_target > 0 else 1.0
    return pct_target / pct_background


def compare_conditions(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    name1: str,
    name2: str,
    output_dir: Path
) -> pd.DataFrame:
    """
    Compare motif enrichment between two conditions.

    Returns DataFrame with differential enrichment statistics.
    """
    logger.info(f"Comparing {name1} vs {name2}...")

    # Merge on motif name
    merged = pd.merge(
        df1[['Motif_Name', 'TF_Name', 'TF_Family', 'Consensus',
             'Log_P_value', 'Pct_Target', 'Pct_Background']],
        df2[['Motif_Name', 'Log_P_value', 'Pct_Target', 'Pct_Background']],
        on='Motif_Name',
        suffixes=(f'_{name1}', f'_{name2}'),
        how='outer'
    )

    # Fill NaN for motifs only in one condition
    merged = merged.fillna(0)

    # Calculate enrichment scores for each condition
    merged[f'Enrichment_{name1}'] = merged.apply(
        lambda x: calculate_enrichment_score(
            x[f'Pct_Target_{name1}'], x[f'Pct_Background_{name1}']
        ), axis=1
    )
    merged[f'Enrichment_{name2}'] = merged.apply(
        lambda x: calculate_enrichment_score(
            x[f'Pct_Target_{name2}'], x[f'Pct_Background_{name2}']
        ), axis=1
    )

    # Calculate differential enrichment (log2 fold change)
    # Positive = more enriched in condition 2 (typically Mut)
    # Negative = more enriched in condition 1 (typically Ctrl)
    merged['Log2_Enrichment_Ratio'] = np.log2(
        (merged[f'Enrichment_{name2}'] + 0.01) /
        (merged[f'Enrichment_{name1}'] + 0.01)
    )

    # Calculate differential significance (-log10 p-value difference)
    merged['Delta_LogP'] = (
        merged[f'Log_P_value_{name2}'] - merged[f'Log_P_value_{name1}']
    )

    # Determine differential status
    def classify_differential(row):
        if row['Log2_Enrichment_Ratio'] > np.log2(FOLD_CHANGE_THRESHOLD):
            if row[f'Log_P_value_{name2}'] > LOG_PVALUE_THRESHOLD:
                return f'Enriched_in_{name2}'
        elif row['Log2_Enrichment_Ratio'] < -np.log2(FOLD_CHANGE_THRESHOLD):
            if row[f'Log_P_value_{name1}'] > LOG_PVALUE_THRESHOLD:
                return f'Enriched_in_{name1}'
        return 'No_Change'

    merged['Differential_Status'] = merged.apply(classify_differential, axis=1)

    # Sort by absolute enrichment ratio
    merged = merged.sort_values('Log2_Enrichment_Ratio', key=abs, ascending=False)

    # Save full results
    output_file = output_dir / f"differential_motifs_{name1}_vs_{name2}.tsv"
    merged.to_csv(output_file, sep='\t', index=False)
    logger.info(f"  Saved full results to: {output_file}")

    # Create summary of significant differential motifs
    sig_mask = merged['Differential_Status'] != 'No_Change'
    sig_df = merged[sig_mask].copy()

    summary_file = output_dir / f"significant_differential_motifs_{name1}_vs_{name2}.tsv"
    sig_df.to_csv(summary_file, sep='\t', index=False)
    logger.info(f"  Significant differential motifs: {len(sig_df)}")
    logger.info(f"    Enriched in {name1}: {(sig_df['Differential_Status'] == f'Enriched_in_{name1}').sum()}")
    logger.info(f"    Enriched in {name2}: {(sig_df['Differential_Status'] == f'Enriched_in_{name2}').sum()}")

    return merged


def create_comparison_matrix(
    results_dict: dict,
    samples: list,
    output_dir: Path
):
    """
    Create a matrix of TF motif enrichment across all samples.
    """
    logger.info("Creating TF enrichment comparison matrix...")

    # Collect all unique motifs
    all_motifs = set()
    for sample, df in results_dict.items():
        if df is not None:
            all_motifs.update(df['Motif_Name'].tolist())

    # Create matrix
    matrix_data = []
    for motif in sorted(all_motifs):
        row = {'Motif_Name': motif}
        for sample in samples:
            df = results_dict.get(sample)
            if df is not None and motif in df['Motif_Name'].values:
                motif_row = df[df['Motif_Name'] == motif].iloc[0]
                row[f'{sample}_LogP'] = motif_row['Log_P_value']
                row[f'{sample}_PctTarget'] = motif_row['Pct_Target']
                row[f'{sample}_Enrichment'] = calculate_enrichment_score(
                    motif_row['Pct_Target'], motif_row['Pct_Background']
                )
            else:
                row[f'{sample}_LogP'] = 0
                row[f'{sample}_PctTarget'] = 0
                row[f'{sample}_Enrichment'] = 1.0
        matrix_data.append(row)

    matrix_df = pd.DataFrame(matrix_data)

    # Extract TF name and family
    matrix_df['TF_Name'] = matrix_df['Motif_Name'].apply(lambda x: x.split('(')[0].strip())
    matrix_df['TF_Family'] = matrix_df['Motif_Name'].apply(
        lambda x: x.split('(')[1].split(')')[0] if '(' in x else 'Unknown'
    )

    # Reorder columns
    cols = ['Motif_Name', 'TF_Name', 'TF_Family'] + \
           [c for c in matrix_df.columns if c not in ['Motif_Name', 'TF_Name', 'TF_Family']]
    matrix_df = matrix_df[cols]

    # Save matrix
    output_file = output_dir / "tf_enrichment_matrix_all_samples.tsv"
    matrix_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved enrichment matrix: {output_file}")

    return matrix_df


def generate_summary_report(
    results_dict: dict,
    comparisons: list,
    output_dir: Path
):
    """
    Generate a summary report of the differential motif analysis.
    """
    logger.info("Generating summary report...")

    report_lines = [
        "# Differential TFBS Motif Analysis Summary",
        f"# Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Overview",
        "",
        "This analysis identifies transcription factor binding site (TFBS) motifs",
        "that are differentially enriched in ATAC peaks at promoter-like elements (PLS)",
        "between control and SRF mutant conditions.",
        "",
        "## Samples Analyzed",
        ""
    ]

    for sample, df in results_dict.items():
        if df is not None:
            report_lines.append(f"- **{sample}**: {len(df)} motifs analyzed")
            sig_count = (df['Log_P_value'] > LOG_PVALUE_THRESHOLD).sum()
            report_lines.append(f"  - Significantly enriched (p < 0.01): {sig_count}")
        else:
            report_lines.append(f"- **{sample}**: Data not available")
    report_lines.append("")

    report_lines.extend([
        "## Comparison Results",
        "",
        f"Significance thresholds:",
        f"- Log10(p-value) > {LOG_PVALUE_THRESHOLD}",
        f"- Fold change > {FOLD_CHANGE_THRESHOLD}",
        ""
    ])

    for comp in comparisons:
        comp_file = output_dir / f"significant_differential_motifs_{comp}.tsv"
        if comp_file.exists():
            df = pd.read_csv(comp_file, sep='\t')
            name1, name2 = comp.split('_vs_')

            report_lines.extend([
                f"### {comp.replace('_', ' ')}",
                "",
                f"- Total differential motifs: {len(df)}",
                f"- Enriched in {name1}: {(df['Differential_Status'] == f'Enriched_in_{name1}').sum()}",
                f"- Enriched in {name2}: {(df['Differential_Status'] == f'Enriched_in_{name2}').sum()}",
                ""
            ])

            # Top differential motifs
            if len(df) > 0:
                report_lines.append("**Top 10 differentially enriched TFs:**")
                report_lines.append("")
                top_df = df.head(10)[['TF_Name', 'TF_Family', 'Log2_Enrichment_Ratio', 'Differential_Status']]
                for _, row in top_df.iterrows():
                    direction = "+" if row['Log2_Enrichment_Ratio'] > 0 else ""
                    report_lines.append(
                        f"  - {row['TF_Name']} ({row['TF_Family']}): "
                        f"{direction}{row['Log2_Enrichment_Ratio']:.2f} log2FC - {row['Differential_Status']}"
                    )
                report_lines.append("")

    report_lines.extend([
        "## Output Files",
        "",
        "- `tf_enrichment_matrix_all_samples.tsv`: Complete enrichment matrix",
        "- `differential_motifs_*.tsv`: Full comparison results",
        "- `significant_differential_motifs_*.tsv`: Filtered significant results",
        "",
        "## Interpretation",
        "",
        "- **Positive Log2_Enrichment_Ratio**: Motif more enriched in mutant condition",
        "- **Negative Log2_Enrichment_Ratio**: Motif more enriched in control condition",
        "",
        "TFs with increased accessibility in mutants may indicate:",
        "1. Compensatory transcriptional programs",
        "2. Altered chromatin state due to SRF loss",
        "3. Secondary effects on splicing gene regulation",
        "",
        "TFs with decreased accessibility in mutants may indicate:",
        "1. Direct or indirect SRF targets",
        "2. Co-regulatory relationships with SRF",
        "3. Affected transcriptional networks",
        ""
    ])

    # Save report
    report_file = output_dir / "DIFFERENTIAL_MOTIF_ANALYSIS_SUMMARY.md"
    with open(report_file, 'w') as f:
        f.write('\n'.join(report_lines))
    logger.info(f"Saved summary report: {report_file}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Compare TFBS motif enrichment between conditions'
    )
    parser.add_argument(
        '--input-dir', '-i',
        type=Path,
        default=HOMER_DIR,
        help=f'Directory containing HOMER results (default: {HOMER_DIR})'
    )
    parser.add_argument(
        '--output-dir', '-o',
        type=Path,
        default=OUTPUT_DIR,
        help=f'Output directory (default: {OUTPUT_DIR})'
    )
    args = parser.parse_args()

    logger.info("=" * 60)
    logger.info("Differential TFBS Motif Analysis")
    logger.info("=" * 60)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load HOMER results for each sample
    results_dict = {}
    for sample in SAMPLES:
        homer_file = args.input_dir / sample / "knownResults.txt"
        results_dict[sample] = parse_homer_results(homer_file)

    # Check if we have any valid data
    valid_samples = [s for s in SAMPLES if results_dict[s] is not None]
    if len(valid_samples) < 2:
        logger.error("Need at least 2 samples with valid HOMER results for comparison")
        sys.exit(1)

    logger.info(f"Valid samples for comparison: {valid_samples}")

    # Perform pairwise comparisons
    comparisons = []

    # 1. Nestin-Ctrl vs Nestin-Mut (mutation effect in Nestin)
    if results_dict.get("Nestin_Ctrl") is not None and results_dict.get("Nestin_Mut") is not None:
        compare_conditions(
            results_dict["Nestin_Ctrl"],
            results_dict["Nestin_Mut"],
            "Nestin_Ctrl",
            "Nestin_Mut",
            args.output_dir
        )
        comparisons.append("Nestin_Ctrl_vs_Nestin_Mut")

    # 2. Nestin-Ctrl vs Emx1-Mut (mutation effect in Emx1)
    if results_dict.get("Nestin_Ctrl") is not None and results_dict.get("Emx1_Mut") is not None:
        compare_conditions(
            results_dict["Nestin_Ctrl"],
            results_dict["Emx1_Mut"],
            "Nestin_Ctrl",
            "Emx1_Mut",
            args.output_dir
        )
        comparisons.append("Nestin_Ctrl_vs_Emx1_Mut")

    # 3. Nestin-Mut vs Emx1-Mut (genotype comparison)
    if results_dict.get("Nestin_Mut") is not None and results_dict.get("Emx1_Mut") is not None:
        compare_conditions(
            results_dict["Nestin_Mut"],
            results_dict["Emx1_Mut"],
            "Nestin_Mut",
            "Emx1_Mut",
            args.output_dir
        )
        comparisons.append("Nestin_Mut_vs_Emx1_Mut")

    # Create comparison matrix
    create_comparison_matrix(results_dict, SAMPLES, args.output_dir)

    # Generate summary report
    generate_summary_report(results_dict, comparisons, args.output_dir)

    logger.info("=" * 60)
    logger.info("Differential Motif Analysis Complete")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
