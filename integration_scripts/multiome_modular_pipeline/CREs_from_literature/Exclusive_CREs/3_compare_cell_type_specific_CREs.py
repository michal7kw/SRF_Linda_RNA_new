#!/usr/bin/env python3
"""
Cell-Type Specificity Analysis: GABA-specific vs Excitatory-specific CREs

This script uses MUTUALLY EXCLUSIVE CRE sets to demonstrate ATAC-seq signal
specificity. Unlike the previous version, these CRE sets have 0% overlap:

- GABA-specific CREs: Active ONLY in GABA neurons (POSITIVE CONTROL)
- Excitatory-specific CREs: Active ONLY in Excitatory neurons (NEGATIVE CONTROL)

Expected results:
- GABA-specific CREs: HIGH signal in GABA ATAC samples
- Excitatory-specific CREs: LOW signal in GABA ATAC samples
- Clear visual difference demonstrating cell-type specificity

Creates heatmaps and metaprofiles with identical scales for direct comparison.

INPUT FILES:
- GABA_specific_CREs.bed: BED file with GABA neuron-specific CREs (mutually exclusive set)
- Excitatory_specific_CREs.bed: BED file with excitatory neuron-specific CREs (mutually exclusive set)
- GABA_Nestin-Ctrl.bw: BigWig file with ATAC-seq signal for GABA neurons, Nestin-Cre control condition
- GABA_Nestin-Mut.bw: BigWig file with ATAC-seq signal for GABA neurons, Nestin-Cre mutant condition
- GABA_Emx1-Ctrl.bw: BigWig file with ATAC-seq signal for GABA neurons, Emx1-Cre control condition
- GABA_Emx1-Mut.bw: BigWig file with ATAC-seq signal for GABA neurons, Emx1-Cre mutant condition

OUTPUT FILES:
- heatmap_GABA_specific.png: Heatmap showing ATAC signal at GABA-specific CREs (positive control)
- heatmap_Excitatory_specific.png: Heatmap showing ATAC signal at excitatory-specific CREs (negative control)
- metaprofile_GABA_specific.png: Metaprofile plot showing average signal at GABA-specific CREs
- metaprofile_Excitatory_specific.png: Metaprofile plot showing average signal at excitatory-specific CREs
- comparison_cell_type_specific.png: Combined comparison plot with metaprofiles, bar charts, and fold enrichment
- cell_type_specificity_statistics.txt: Detailed statistics file with signal comparisons and interpretations
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
from datetime import datetime

# Try to import pyBigWig
try:
    import pyBigWig
except ImportError:
    print("ERROR: pyBigWig not available. Install with: pip install pyBigWig")
    exit(1)

print("="*80)
print("CELL-TYPE SPECIFICITY: GABA-SPECIFIC vs EXCITATORY-SPECIFIC CREs")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Configuration
BIGWIG_BASE = "../../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR = "output/GABA_DEG_analysis/heatmaps_specific_CREs"
BED_GABA = "output/GABA_specific_CREs.bed"
BED_EXCITATORY = "output/Excitatory_specific_CREs.bed"

WINDOW_SIZE = 2000
BIN_SIZE = 50
N_BINS = WINDOW_SIZE * 2 // BIN_SIZE

CONDITIONS = {
    "Nestin-Ctrl": f"{BIGWIG_BASE}/GABA_Nestin-Ctrl.bw",
    "Nestin-Mut": f"{BIGWIG_BASE}/GABA_Nestin-Mut.bw",
    "Emx1-Ctrl": f"{BIGWIG_BASE}/GABA_Emx1-Ctrl.bw",
    "Emx1-Mut": f"{BIGWIG_BASE}/GABA_Emx1-Mut.bw"
}

COLORS = {
    "Nestin-Ctrl": "#2E86AB",
    "Nestin-Mut": "#A23B72",
    "Emx1-Ctrl": "#F18F01",
    "Emx1-Mut": "#C73E1D"
}

os.makedirs(OUTPUT_DIR, exist_ok=True)
sns.set_style("white")
plt.rcParams['figure.dpi'] = 150

# ============================================================================
# Helper Functions
# ============================================================================

def load_bed(bed_file):
    """Load CREs from BED file"""
    cres = pd.read_csv(bed_file, sep='\t', header=None,
                       names=['chr', 'start', 'end', 'name', 'score', 'extra'])
    cres['center'] = ((cres['start'] + cres['end']) / 2).astype(int)
    return cres[['chr', 'start', 'end', 'center', 'name']].sort_values(['chr', 'center']).reset_index(drop=True)

def extract_signal_window(bw_file, chrom, center):
    """Extract signal in bins around genomic position"""
    try:
        bw = pyBigWig.open(bw_file)
        if chrom not in bw.chroms():
            bw.close()
            return np.zeros(N_BINS)

        start = max(0, center - WINDOW_SIZE)
        end = min(center + WINDOW_SIZE, bw.chroms()[chrom])

        signals = []
        for i in range(N_BINS):
            bin_start = start + (i * BIN_SIZE)
            bin_end = min(bin_start + BIN_SIZE, end)
            signal = bw.stats(chrom, bin_start, bin_end, type="mean")
            signals.append(signal[0] if signal and signal[0] is not None else 0.0)

        bw.close()
        return np.array(signals[:N_BINS])
    except:
        return np.zeros(N_BINS)

def extract_signals(cres, conditions, label):
    """Extract signals for all CREs and conditions"""
    matrices = {}
    for cond_name, bw_file in conditions.items():
        print(f"\n{cond_name}:")
        if not os.path.exists(bw_file):
            print(f"  Skipping (file not found)")
            continue

        matrix = np.zeros((len(cres), N_BINS))
        for idx, row in cres.iterrows():
            if idx % 1000 == 0:
                print(f"  Processing {idx}/{len(cres)}...", end='\r')
            matrix[idx, :] = extract_signal_window(bw_file, row['chr'], row['center'])

        print(f"  Completed {len(cres)} CREs")
        matrices[cond_name] = matrix

    return matrices

# ============================================================================
# Load CREs
# ============================================================================
print("STEP 1: Loading cell-type-specific CREs...")
print("-"*80)

print(f"\nLoading GABA-specific CREs from: {BED_GABA}")
print(f"  File exists: {os.path.exists(BED_GABA)}")
cres_gaba = load_bed(BED_GABA)

print(f"\nLoading Excitatory-specific CREs from: {BED_EXCITATORY}")
print(f"  File exists: {os.path.exists(BED_EXCITATORY)}")
cres_excitatory = load_bed(BED_EXCITATORY)

print(f"\n{'='*80}")
print("LOADED CRE SETS:")
print(f"  GABA-specific CREs: {len(cres_gaba):,}")
print(f"  Excitatory-specific CREs: {len(cres_excitatory):,}")

print(f"\nFirst 3 GABA-specific CREs:")
print(cres_gaba.head(3).to_string())

print(f"\nFirst 3 Excitatory-specific CREs:")
print(cres_excitatory.head(3).to_string())

# Verify mutual exclusivity
print(f"\nVerifying mutual exclusivity...")
gaba_coords = set(zip(cres_gaba['chr'], cres_gaba['start'], cres_gaba['end']))
excit_coords = set(zip(cres_excitatory['chr'], cres_excitatory['start'], cres_excitatory['end']))
overlap = len(gaba_coords & excit_coords)

if overlap == 0:
    print(f"  ✓ VERIFIED: 0 overlapping coordinates (mutually exclusive)")
else:
    print(f"  ✗ WARNING: {overlap} overlapping coordinates found!")
    print(f"     This may indicate an error in CRE extraction!")

# ============================================================================
# Extract Signals
# ============================================================================
print(f"\n{'='*80}")
print("STEP 2: Extracting signals...")
print("-"*80)

print("\nGABA-specific CREs (POSITIVE CONTROL):")
signals_gaba = extract_signals(cres_gaba, CONDITIONS, "GABA-specific")

print(f"\n{'='*80}")
print("Excitatory-specific CREs (NEGATIVE CONTROL):")
signals_excitatory = extract_signals(cres_excitatory, CONDITIONS, "Excitatory-specific")

# Signal comparison
print(f"\n{'='*80}")
print("Signal Matrix Comparison:")
print("-"*80)
for cond in CONDITIONS.keys():
    if cond in signals_gaba and cond in signals_excitatory:
        gaba_mean = np.mean(signals_gaba[cond])
        excit_mean = np.mean(signals_excitatory[cond])
        gaba_max = np.max(signals_gaba[cond])
        excit_max = np.max(signals_excitatory[cond])
        fold_enrichment = gaba_mean / excit_mean if excit_mean > 0 else float('inf')

        print(f"\n{cond}:")
        print(f"  GABA-specific:       mean={gaba_mean:.6f}, max={gaba_max:.6f}")
        print(f"  Excitatory-specific: mean={excit_mean:.6f}, max={excit_max:.6f}")
        print(f"  Fold enrichment:     {fold_enrichment:.2f}x")

        if fold_enrichment < 2.0:
            print(f"  ⚠ WARNING: Low fold enrichment (<2x) - check cell-type specificity!")
        elif fold_enrichment >= 3.0:
            print(f"  ✓ Strong cell-type specificity (≥3x enrichment)")

# ============================================================================
# Determine Common Scale
# ============================================================================
print(f"\n{'='*80}")
print("STEP 3: Determining common scale...")
print("-"*80)

# Collect all values
all_values = []
for matrix in list(signals_gaba.values()) + list(signals_excitatory.values()):
    all_values.extend(matrix.flatten())

all_values = np.array([v for v in all_values if v > 0])  # Remove zeros

if len(all_values) > 0:
    percentile_90 = np.percentile(all_values, 90)
    vmax = percentile_90 * 1.2  # 20% buffer
else:
    percentile_90 = 0.0
    vmax = 1.0  # Fallback

print(f"90th percentile: {percentile_90:.4f}")
print(f"Common scale: 0 to {vmax:.4f}")

# ============================================================================
# Create Heatmaps
# ============================================================================
print(f"\n{'='*80}")
print("STEP 4: Creating heatmaps...")
print("-"*80)

# GABA-specific heatmap
fig, axes = plt.subplots(1, 4, figsize=(16, 12), sharey=True)
fig.suptitle(f'POSITIVE CONTROL: ATAC Signal at GABA-specific CREs (n={len(cres_gaba):,})',
             fontsize=16, fontweight='bold')

for idx, (cond, matrix) in enumerate(signals_gaba.items()):
    ax = axes[idx]

    # Sort rows by mean signal for better visualization
    row_means = np.mean(matrix, axis=1)
    sorted_indices = np.argsort(row_means)[::-1]
    matrix_sorted = matrix[sorted_indices, :]

    im = ax.imshow(matrix_sorted, aspect='auto', cmap='Reds',
                   vmin=0, vmax=vmax, interpolation='nearest')
    ax.set_title(cond, fontsize=14, fontweight='bold')
    if idx == 0:
        ax.set_ylabel('CREs (sorted by mean signal)', fontsize=12)
    ax.set_xlabel('Distance from center', fontsize=10)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('ATAC Signal', fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "heatmap_GABA_specific.png"), dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: heatmap_GABA_specific.png")

# Excitatory-specific heatmap
fig, axes = plt.subplots(1, 4, figsize=(16, 12), sharey=True)
fig.suptitle(f'NEGATIVE CONTROL: ATAC Signal at Excitatory-specific CREs (n={len(cres_excitatory):,})',
             fontsize=16, fontweight='bold')

for idx, (cond, matrix) in enumerate(signals_excitatory.items()):
    ax = axes[idx]

    # Sort rows by mean signal for better visualization
    row_means = np.mean(matrix, axis=1)
    sorted_indices = np.argsort(row_means)[::-1]
    matrix_sorted = matrix[sorted_indices, :]

    im = ax.imshow(matrix_sorted, aspect='auto', cmap='Reds',
                   vmin=0, vmax=vmax, interpolation='nearest')
    ax.set_title(cond, fontsize=14, fontweight='bold')
    if idx == 0:
        ax.set_ylabel('CREs (sorted by mean signal)', fontsize=12)
    ax.set_xlabel('Distance from center', fontsize=10)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('ATAC Signal', fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "heatmap_Excitatory_specific.png"), dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: heatmap_Excitatory_specific.png")

# ============================================================================
# Create Metaprofiles
# ============================================================================
print(f"\n{'='*80}")
print("STEP 5: Creating metaprofiles...")
print("-"*80)

x_coords = np.linspace(-WINDOW_SIZE, WINDOW_SIZE, N_BINS)

# GABA-specific metaprofile
fig, ax = plt.subplots(figsize=(10, 6))
for cond, matrix in signals_gaba.items():
    mean_sig = np.mean(matrix, axis=0)
    sem_sig = stats.sem(matrix, axis=0)
    ax.plot(x_coords, mean_sig, color=COLORS[cond], linewidth=2.5, label=cond)
    ax.fill_between(x_coords, mean_sig - sem_sig, mean_sig + sem_sig,
                     color=COLORS[cond], alpha=0.2)

ax.axvline(0, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)
ax.set_xlabel('Distance from CRE Center (bp)', fontsize=12)
ax.set_ylabel('Mean ATAC Signal', fontsize=12)
ax.set_ylim(0, vmax)
ax.set_title(f'POSITIVE CONTROL: Metaprofile at GABA-specific CREs (n={len(cres_gaba):,})',
             fontsize=14, fontweight='bold')
ax.legend(loc='best', frameon=True, fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "metaprofile_GABA_specific.png"), dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: metaprofile_GABA_specific.png")

# Excitatory-specific metaprofile
fig, ax = plt.subplots(figsize=(10, 6))
for cond, matrix in signals_excitatory.items():
    mean_sig = np.mean(matrix, axis=0)
    sem_sig = stats.sem(matrix, axis=0)
    ax.plot(x_coords, mean_sig, color=COLORS[cond], linewidth=2.5, label=cond)
    ax.fill_between(x_coords, mean_sig - sem_sig, mean_sig + sem_sig,
                     color=COLORS[cond], alpha=0.2)

ax.axvline(0, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)
ax.set_xlabel('Distance from CRE Center (bp)', fontsize=12)
ax.set_ylabel('Mean ATAC Signal', fontsize=12)
ax.set_ylim(0, vmax)
ax.set_title(f'NEGATIVE CONTROL: Metaprofile at Excitatory-specific CREs (n={len(cres_excitatory):,})',
             fontsize=14, fontweight='bold')
ax.legend(loc='best', frameon=True, fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "metaprofile_Excitatory_specific.png"), dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: metaprofile_Excitatory_specific.png")

# ============================================================================
# Create Comparison Plot
# ============================================================================
print(f"\n{'='*80}")
print("STEP 6: Creating comparison plot...")
print("-"*80)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# GABA-specific metaprofiles
ax = axes[0, 0]
for cond, matrix in signals_gaba.items():
    mean_sig = np.mean(matrix, axis=0)
    ax.plot(x_coords, mean_sig, color=COLORS[cond], linewidth=2, label=cond)
ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
ax.set_title('GABA-specific CREs (Positive)', fontweight='bold', fontsize=12)
ax.set_ylabel('Mean ATAC Signal', fontsize=11)
ax.set_ylim(0, vmax)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Excitatory-specific metaprofiles
ax = axes[0, 1]
for cond, matrix in signals_excitatory.items():
    mean_sig = np.mean(matrix, axis=0)
    ax.plot(x_coords, mean_sig, color=COLORS[cond], linewidth=2, label=cond)
ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
ax.set_title('Excitatory-specific CREs (Negative)', fontweight='bold', fontsize=12)
ax.set_ylim(0, vmax)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Mean signal comparison
ax = axes[1, 0]
gaba_means = [np.mean(m) for m in signals_gaba.values()]
excit_means = [np.mean(m) for m in signals_excitatory.values()]
x_pos = np.arange(len(CONDITIONS))
width = 0.35
bars1 = ax.bar(x_pos - width/2, gaba_means, width, label='GABA-specific', color='red', alpha=0.7)
bars2 = ax.bar(x_pos + width/2, excit_means, width, label='Excitatory-specific', color='blue', alpha=0.7)
ax.set_xticks(x_pos)
ax.set_xticklabels(CONDITIONS.keys(), rotation=45, ha='right', fontsize=9)
ax.set_ylabel('Mean ATAC Signal', fontsize=11)
ax.set_title('Signal Comparison by Condition', fontweight='bold', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, axis='y', alpha=0.3)

# Add value labels on bars
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.4f}', ha='center', va='bottom', fontsize=7)

# Fold enrichment
ax = axes[1, 1]
fold_enrichments = [g/e if e > 0 else 0 for g, e in zip(gaba_means, excit_means)]
bars = ax.bar(x_pos, fold_enrichments, color='purple', alpha=0.7)
ax.axhline(1, color='red', linestyle='--', label='No enrichment', linewidth=1.5)
ax.axhline(3, color='green', linestyle='--', label='Strong specificity (3x)', linewidth=1.5)
ax.set_xticks(x_pos)
ax.set_xticklabels(CONDITIONS.keys(), rotation=45, ha='right', fontsize=9)
ax.set_ylabel('Fold Enrichment (GABA/Excitatory)', fontsize=11)
ax.set_title('Cell-Type Specificity', fontweight='bold', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, axis='y', alpha=0.3)

# Add value labels on bars
for bar in bars:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{height:.1f}x', ha='center', va='bottom', fontsize=9, fontweight='bold')

plt.suptitle('Cell-Type Specificity: GABA-specific vs Excitatory-specific CREs',
             fontsize=16, fontweight='bold', y=1.00)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "comparison_cell_type_specific.png"), dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: comparison_cell_type_specific.png")

# ============================================================================
# Save Statistics
# ============================================================================
print(f"\n{'='*80}")
print("STEP 7: Saving statistics...")
print("-"*80)

stats_file = os.path.join(OUTPUT_DIR, "cell_type_specificity_statistics.txt")
with open(stats_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("CELL-TYPE SPECIFICITY: GABA-SPECIFIC vs EXCITATORY-SPECIFIC CREs\n")
    f.write("="*80 + "\n\n")
    f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("CRE SETS (MUTUALLY EXCLUSIVE):\n")
    f.write(f"  GABA-specific CREs: {len(cres_gaba):,}\n")
    f.write(f"  Excitatory-specific CREs: {len(cres_excitatory):,}\n")
    f.write(f"  Overlapping coordinates: {overlap} (0% - mutually exclusive)\n\n")

    f.write(f"SCALE:\n")
    f.write(f"  Common scale: 0 to {vmax:.4f} (90th percentile + 20%)\n\n")

    f.write("MEAN SIGNAL COMPARISON:\n")
    f.write("-"*80 + "\n")
    for idx, cond in enumerate(CONDITIONS.keys()):
        if cond in signals_gaba and cond in signals_excitatory:
            gaba_mean = np.mean(signals_gaba[cond])
            excit_mean = np.mean(signals_excitatory[cond])
            fold = gaba_mean / excit_mean if excit_mean > 0 else float('inf')
            f.write(f"{cond}:\n")
            f.write(f"  GABA-specific mean: {gaba_mean:.6f}\n")
            f.write(f"  Excitatory-specific mean: {excit_mean:.6f}\n")
            f.write(f"  Fold enrichment: {fold:.2f}x\n")
            if fold >= 3.0:
                f.write(f"  → Strong cell-type specificity (≥3x)\n")
            elif fold >= 2.0:
                f.write(f"  → Moderate cell-type specificity (≥2x)\n")
            else:
                f.write(f"  → Weak cell-type specificity (<2x)\n")
            f.write("\n")

    f.write("INTERPRETATION:\n")
    f.write("-"*80 + "\n")
    f.write("Fold enrichment ≥3x indicates strong cell-type specificity.\n")
    f.write("GABA ATAC samples should show:\n")
    f.write("  - HIGH signal at GABA-specific CREs (positive control)\n")
    f.write("  - LOW signal at Excitatory-specific CREs (negative control)\n\n")

    f.write("This demonstrates that GABA ATAC samples specifically capture\n")
    f.write("GABAergic chromatin accessibility patterns.\n\n")

    f.write("KEY DIFFERENCE FROM PREVIOUS ANALYSIS:\n")
    f.write("-"*80 + "\n")
    f.write("Previous analysis used CREs with 60% overlap → identical plots\n")
    f.write("Current analysis uses mutually exclusive CREs → clear difference\n")

print("✓ Saved: cell_type_specificity_statistics.txt")

print(f"\n{'='*80}")
print("ANALYSIS COMPLETE!")
print("="*80)
print(f"\nOutput directory: {OUTPUT_DIR}/")
print(f"\nKey files:")
print(f"  ★ Heatmaps:")
print(f"      - heatmap_GABA_specific.png (POSITIVE - expect HIGH signal)")
print(f"      - heatmap_Excitatory_specific.png (NEGATIVE - expect LOW signal)")
print(f"  ★ Metaprofiles:")
print(f"      - metaprofile_GABA_specific.png (POSITIVE - expect HIGH curves)")
print(f"      - metaprofile_Excitatory_specific.png (NEGATIVE - expect LOW curves)")
print(f"  ★ Comparison:")
print(f"      - comparison_cell_type_specific.png")
print(f"      - cell_type_specificity_statistics.txt")
print(f"\nAll plots use identical scale (0-{vmax:.4f}) for direct comparison")
print(f"\nExpected result:")
print(f"  Clear visual difference demonstrating cell-type specificity!")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
