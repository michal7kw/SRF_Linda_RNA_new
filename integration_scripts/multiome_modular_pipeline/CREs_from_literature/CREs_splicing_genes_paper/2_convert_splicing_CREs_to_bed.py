#!/usr/bin/env python3
"""
Convert Splicing Gene CRE TSV Files to BED Format

This script converts the CRE-gene linkage TSV files (from Table 16 literature data)
to BED format for use with deepTools visualization tools.

The input TSV files contain genomic coordinates in Table 16 format:
- Coordinate1: chr_start_end format (e.g., "chr11_40754370_40756166")
- cCRE1: CRE identifier

This script parses the Coordinate1 column and creates standard BED6 format files.

INPUT:
- output/splicing_genes_analysis/splicing_genes_CREs_all_celltypes.tsv
- output/splicing_genes_analysis/splicing_genes_CREs_GABA.tsv

OUTPUT:
- output/splicing_genes_analysis/splicing_genes_CREs_all.bed (BED6 format)
- output/splicing_genes_analysis/splicing_genes_CREs_GABA.bed (BED6 format)
"""

import pandas as pd
import os

os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/CREs_splicing_genes_paper")

def convert_tsv_to_bed(tsv_file, bed_file, description):
    """
    Convert TSV file to BED format by extracting unique CREs.

    Parameters:
    -----------
    tsv_file : str
        Input TSV file path
    bed_file : str
        Output BED file path
    description : str
        Description for logging
    """
    print(f"\n{'='*80}")
    print(f"Converting {description}")
    print(f"{'='*80}")
    print(f"Input:  {tsv_file}")
    print(f"Output: {bed_file}")

    # Read TSV file
    df = pd.read_csv(tsv_file, sep='\t')
    print(f"\nTotal CRE-gene links: {len(df):,}")

    # Parse Coordinate1 column (format: chr11_40754370_40756166)
    # Extract chr, start, end from Coordinate1
    print("\nParsing genomic coordinates from Coordinate1 column...")

    def parse_coordinate(coord_str):
        """Parse coordinate string to chr, start, end"""
        parts = coord_str.split('_')
        if len(parts) == 3:
            return parts[0], int(parts[1]), int(parts[2])
        else:
            raise ValueError(f"Invalid coordinate format: {coord_str}")

    # Apply parsing to Coordinate1 column
    df[['chr', 'start', 'end']] = df['Coordinate1'].apply(
        lambda x: pd.Series(parse_coordinate(x))
    )

    # Extract unique CREs with coordinates
    # Each CRE is defined by: chr, start, end, cCRE1 (ID)
    unique_cres = df[['chr', 'start', 'end', 'cCRE1']].drop_duplicates()

    # Sort by chromosome and start position
    # Convert chr to categorical for proper sorting
    chr_order = [f'chr{i}' for i in range(1, 20)] + ['chrX', 'chrY', 'chrM']
    unique_cres['chr'] = pd.Categorical(unique_cres['chr'], categories=chr_order, ordered=True)
    unique_cres = unique_cres.sort_values(['chr', 'start'])

    # Add score column (all 1000 for now)
    unique_cres['score'] = 1000

    # Add strand column (all . for unstranded)
    unique_cres['strand'] = '.'

    # Reorder columns for BED6 format: chr, start, end, name, score, strand
    bed_df = unique_cres[['chr', 'start', 'end', 'cCRE1', 'score', 'strand']]

    # Save to BED file
    bed_df.to_csv(bed_file, sep='\t', header=False, index=False)

    print(f"\nUnique CREs extracted: {len(bed_df):,}")
    print(f"\nBED file saved to: {bed_file}")

    # Show first few entries
    print(f"\nFirst 5 CREs:")
    print(bed_df.head().to_string(index=False, header=False))

    # Get genes associated with these CREs
    genes = df['Gene'].unique()
    print(f"\nGenes represented: {len(genes)}")
    print(f"Gene list: {', '.join(sorted(genes))}")

    return bed_df

def main():
    """Main execution function"""
    print("="*80)
    print("CONVERT SPLICING GENE CRE TSV FILES TO BED FORMAT")
    print("="*80)

    # Define input/output paths
    output_dir = "./output"

    tsv_files = {
        'all_celltypes': {
            'tsv': f"{output_dir}/splicing_genes_CREs_all_celltypes.tsv",
            'bed': f"{output_dir}/splicing_genes_CREs_all.bed",
            'description': "All cell types"
        },
        'gaba': {
            'tsv': f"{output_dir}/splicing_genes_CREs_GABA.tsv",
            'bed': f"{output_dir}/splicing_genes_CREs_GABA.bed",
            'description': "GABA cell types"
        }
    }

    # Convert each file
    results = {}
    for key, files in tsv_files.items():
        if os.path.exists(files['tsv']):
            results[key] = convert_tsv_to_bed(
                files['tsv'],
                files['bed'],
                files['description']
            )
        else:
            print(f"\nWARNING: File not found: {files['tsv']}")

    # Summary
    print(f"\n{'='*80}")
    print("CONVERSION COMPLETE")
    print(f"{'='*80}")
    print(f"\nGenerated BED files:")
    for key, files in tsv_files.items():
        if key in results:
            print(f"  {files['bed']}")

    print(f"\nThese BED files are ready for use with deepTools:")
    print(f"  - computeMatrix")
    print(f"  - plotHeatmap")
    print(f"  - plotProfile")
    print()

if __name__ == "__main__":
    main()
