#!/usr/bin/env python3
"""
Analyze CRE lengths from hippocampal interneuron BED file.
Calculates mean length and plots histogram of CRE length distribution.

INPUT FILES:
- hippocampal_interneuron_CREs_with_header.bed: BED file containing CRE coordinates with header

OUTPUT FILES:
- CRE_length_analysis.png: Comprehensive 4-panel visualization (histogram, log-scale, box plot, cumulative distribution)
- CRE_length_histogram_simple.png: Simple single histogram visualization
- CRE_length_statistics.txt: Text file with detailed length statistics and categories
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

# Input and output paths
input_file = Path("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/output/hippocampal_interneuron_CREs_with_header.bed")
output_dir = Path("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature/output")

print(f"Reading BED file: {input_file}")

# Read BED file
df = pd.read_csv(input_file, sep='\t')

print(f"Total number of CREs: {len(df):,}")
print(f"Columns: {df.columns.tolist()}")

# Calculate CRE lengths
df['length'] = df['end'] - df['start']

# Calculate statistics
mean_length = df['length'].mean()
median_length = df['length'].median()
std_length = df['length'].std()
min_length = df['length'].min()
max_length = df['length'].max()
q25 = df['length'].quantile(0.25)
q75 = df['length'].quantile(0.75)

print("\n" + "="*70)
print("CRE LENGTH STATISTICS")
print("="*70)
print(f"Mean length:       {mean_length:>12,.2f} bp")
print(f"Median length:     {median_length:>12,.2f} bp")
print(f"Std deviation:     {std_length:>12,.2f} bp")
print(f"Min length:        {min_length:>12,} bp")
print(f"Max length:        {max_length:>12,} bp")
print(f"25th percentile:   {q25:>12,.2f} bp")
print(f"75th percentile:   {q75:>12,.2f} bp")
print("="*70)

# Create comprehensive visualization
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 1. Main histogram
ax1 = axes[0, 0]
n, bins, patches = ax1.hist(df['length'], bins=50, color='steelblue',
                             edgecolor='black', alpha=0.7)
ax1.axvline(mean_length, color='red', linestyle='--', linewidth=2,
            label=f'Mean: {mean_length:,.0f} bp')
ax1.axvline(median_length, color='orange', linestyle='--', linewidth=2,
            label=f'Median: {median_length:,.0f} bp')
ax1.set_xlabel('CRE Length (bp)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Frequency', fontsize=12, fontweight='bold')
ax1.set_title('Distribution of Hippocampal Interneuron CRE Lengths',
              fontsize=14, fontweight='bold', pad=20)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# 2. Log-scale histogram for better visualization of tail
ax2 = axes[0, 1]
ax2.hist(df['length'], bins=50, color='teal', edgecolor='black', alpha=0.7)
ax2.set_xlabel('CRE Length (bp)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Frequency (log scale)', fontsize=12, fontweight='bold')
ax2.set_title('CRE Length Distribution (Log Scale)',
              fontsize=14, fontweight='bold', pad=20)
ax2.set_yscale('log')
ax2.grid(True, alpha=0.3, which='both')

# 3. Box plot
ax3 = axes[1, 0]
bp = ax3.boxplot(df['length'], vert=True, patch_artist=True,
                  widths=0.5, showmeans=True,
                  meanprops=dict(marker='D', markerfacecolor='red',
                                markeredgecolor='red', markersize=8),
                  medianprops=dict(color='orange', linewidth=2),
                  boxprops=dict(facecolor='lightblue', alpha=0.7),
                  flierprops=dict(marker='o', markerfacecolor='gray',
                                 markersize=4, alpha=0.5))
ax3.set_ylabel('CRE Length (bp)', fontsize=12, fontweight='bold')
ax3.set_title('CRE Length Distribution (Box Plot)',
              fontsize=14, fontweight='bold', pad=20)
ax3.set_xticklabels(['Hippocampal\nInterneuron CREs'])
ax3.grid(True, alpha=0.3, axis='y')

# 4. Cumulative distribution
ax4 = axes[1, 1]
sorted_lengths = np.sort(df['length'])
cumulative = np.arange(1, len(sorted_lengths) + 1) / len(sorted_lengths) * 100
ax4.plot(sorted_lengths, cumulative, linewidth=2, color='darkgreen')
ax4.axhline(50, color='orange', linestyle='--', linewidth=1.5,
            label=f'50th percentile: {median_length:,.0f} bp')
ax4.axhline(25, color='gray', linestyle='--', linewidth=1, alpha=0.7,
            label=f'25th percentile: {q25:,.0f} bp')
ax4.axhline(75, color='gray', linestyle='--', linewidth=1, alpha=0.7,
            label=f'75th percentile: {q75:,.0f} bp')
ax4.set_xlabel('CRE Length (bp)', fontsize=12, fontweight='bold')
ax4.set_ylabel('Cumulative Percentage (%)', fontsize=12, fontweight='bold')
ax4.set_title('Cumulative Distribution of CRE Lengths',
              fontsize=14, fontweight='bold', pad=20)
ax4.legend(fontsize=9)
ax4.grid(True, alpha=0.3)

plt.tight_layout()

# Save figure
output_file = output_dir / "CRE_length_analysis.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nFigure saved to: {output_file}")

# Create a simple single histogram for quick reference
fig2, ax = plt.subplots(figsize=(10, 6))
ax.hist(df['length'], bins=50, color='steelblue', edgecolor='black', alpha=0.7)
ax.axvline(mean_length, color='red', linestyle='--', linewidth=2,
           label=f'Mean: {mean_length:,.0f} bp')
ax.axvline(median_length, color='orange', linestyle='--', linewidth=2,
           label=f'Median: {median_length:,.0f} bp')
ax.set_xlabel('CRE Length (bp)', fontsize=12, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
ax.set_title(f'Hippocampal Interneuron CRE Lengths (n={len(df):,})',
             fontsize=14, fontweight='bold', pad=20)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

output_file_simple = output_dir / "CRE_length_histogram_simple.png"
plt.tight_layout()
plt.savefig(output_file_simple, dpi=300, bbox_inches='tight')
print(f"Simple histogram saved to: {output_file_simple}")

# Save statistics to file
stats_file = output_dir / "CRE_length_statistics.txt"
with open(stats_file, 'w') as f:
    f.write("="*70 + "\n")
    f.write("HIPPOCAMPAL INTERNEURON CRE LENGTH STATISTICS\n")
    f.write("="*70 + "\n\n")
    f.write(f"Input file: {input_file}\n")
    f.write(f"Total number of CREs: {len(df):,}\n\n")
    f.write(f"Mean length:       {mean_length:12,.2f} bp\n")
    f.write(f"Median length:     {median_length:12,.2f} bp\n")
    f.write(f"Std deviation:     {std_length:12,.2f} bp\n")
    f.write(f"Min length:        {min_length:12,} bp\n")
    f.write(f"Max length:        {max_length:12,} bp\n")
    f.write(f"25th percentile:   {q25:12,.2f} bp\n")
    f.write(f"75th percentile:   {q75:12,.2f} bp\n")
    f.write("="*70 + "\n\n")

    # Length categories
    f.write("LENGTH CATEGORIES:\n")
    f.write("-"*70 + "\n")
    categories = [
        ("Very short (<500 bp)", (df['length'] < 500).sum()),
        ("Short (500-1000 bp)", ((df['length'] >= 500) & (df['length'] < 1000)).sum()),
        ("Medium (1000-2000 bp)", ((df['length'] >= 1000) & (df['length'] < 2000)).sum()),
        ("Long (2000-5000 bp)", ((df['length'] >= 2000) & (df['length'] < 5000)).sum()),
        ("Very long (>5000 bp)", (df['length'] >= 5000).sum())
    ]

    for category, count in categories:
        percentage = count / len(df) * 100
        f.write(f"{category:25s}: {count:6,} ({percentage:5.1f}%)\n")
    f.write("="*70 + "\n")

print(f"Statistics saved to: {stats_file}")

# Print length categories to console
print("\nLENGTH CATEGORIES:")
print("-"*70)
for category, count in categories:
    percentage = count / len(df) * 100
    print(f"{category:25s}: {count:6,} ({percentage:5.1f}%)")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
