#!/usr/bin/env python3
"""
Cell-Type Specificity Analysis: GABA vs Excitatory CREs

Demonstrates ATAC-seq signal specificity by comparing:
- GABA CREs (POSITIVE CONTROL): Expect HIGH signal in GABA ATAC samples
- Excitatory CREs (NEGATIVE CONTROL): Expect LOW signal in GABA ATAC samples

Creates heatmaps and metaprofiles with identical scales for direct comparison.

INPUT FILES:
- hippocampal_interneuron_CREs.bed: BED file containing GABA neuron CRE coordinates
- excitatory_neuron_CREs.bed: BED file containing excitatory neuron CRE coordinates
- GABA_Nestin-Ctrl.bw: BigWig file with ATAC-seq signal for GABA neurons, Nestin-Cre control condition
- GABA_Nestin-Mut.bw: BigWig file with ATAC-seq signal for GABA neurons, Nestin-Cre mutant condition
- GABA_Emx1-Ctrl.bw: BigWig file with ATAC-seq signal for GABA neurons, Emx1-Cre control condition
- GABA_Emx1-Mut.bw: BigWig file with ATAC-seq signal for GABA neurons, Emx1-Cre mutant condition

OUTPUT FILES:
- heatmap_GABA_all_conditions.png: Heatmap showing ATAC signal at GABA CREs across all conditions
- heatmap_Excitatory_all_conditions.png: Heatmap showing ATAC signal at excitatory CREs across all conditions
- metaprofile_GABA_all_CREs.png: Metaprofile plot showing average signal at GABA CREs
- metaprofile_Excitatory_all_CREs.png: Metaprofile plot showing average signal at excitatory CREs
- comparison_GABA_vs_Excitatory.png: Combined comparison plot with metaprofiles, bar charts, and fold enrichment
- enrichment_statistics.txt: Statistics file with quantitative signal comparisons and fold enrichment values
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
print("CELL-TYPE SPECIFICITY: GABA vs Excitatory CREs")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Configuration
BIGWIG_BASE = "../signac_results_L1/bigwig_tracks_L1/by_celltype"
OUTPUT_DIR = "output/GABA_DEG_analysis/heatmaps_metaprofiles"
BED_GABA = "output/hippocampal_interneuron_CREs.bed"
BED_EXCITATORY = "output/excitatory_neuron_CREs.bed"

WINDOW_SIZE = 2000
BIN_SIZE = 50
N_BINS = WINDOW_SIZE * 2 // BIN_SIZE

# NOTE: Emx1-Ctrl is excluded (failed sample)
CONDITIONS = {
    "Nestin-Ctrl": f"{BIGWIG_BASE}/GABA_Nestin-Ctrl.bw",
    "Nestin-Mut": f"{BIGWIG_BASE}/GABA_Nestin-Mut.bw",
    "Emx1-Mut": f"{BIGWIG_BASE}/GABA_Emx1-Mut.bw"
}

COLORS = {
    "Nestin-Ctrl": "#2E86AB",
    "Nestin-Mut": "#A23B72",
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

def extract_signals(cres, conditions):
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
print("STEP 1: Loading CREs...")
print("-"*80)

print(f"\nDEBUG: Loading GABA CREs from: {BED_GABA}")
print(f"  File exists: {os.path.exists(BED_GABA)}")
cres_gaba = load_bed(BED_GABA)

print(f"\nDEBUG: Loading Excitatory CREs from: {BED_EXCITATORY}")
print(f"  File exists: {os.path.exists(BED_EXCITATORY)}")
cres_excitatory = load_bed(BED_EXCITATORY)

print(f"\n{'='*80}")
print("LOADED CRE SETS:")
print(f"  GABA CREs: {len(cres_gaba)}")
print(f"  Excitatory CREs: {len(cres_excitatory)}")

print(f"\nDEBUG: First 3 GABA CREs:")
print(cres_gaba.head(3).to_string())

print(f"\nDEBUG: First 3 Excitatory CREs:")
print(cres_excitatory.head(3).to_string())

print(f"\nDEBUG: Comparing CRE sets...")
if len(cres_gaba) == len(cres_excitatory):
    # Check if they're identical
    gaba_coords = set(zip(cres_gaba['chr'], cres_gaba['start'], cres_gaba['end']))
    excit_coords = set(zip(cres_excitatory['chr'], cres_excitatory['start'], cres_excitatory['end']))
    overlap = len(gaba_coords & excit_coords)
    print(f"  WARNING: Same number of CREs ({len(cres_gaba)})")
    print(f"  Overlapping coordinates: {overlap}/{len(cres_gaba)} ({100*overlap/len(cres_gaba):.1f}%)")
    if overlap > len(cres_gaba) * 0.9:
        print(f"  ERROR: >90% overlap - CRE sets are nearly identical!")
        print(f"  This indicates the wrong BED file may have been loaded!")
else:
    print(f"  ✓ Different number of CREs - sets are distinct")

# ============================================================================
# Extract Signals
# ============================================================================
print(f"\n{'='*80}")
print("STEP 2: Extracting signals...")
print("-"*80)

print("\nGABA CREs (POSITIVE CONTROL):")
signals_gaba = extract_signals(cres_gaba, CONDITIONS)

print(f"\n{'='*80}")
print("Excitatory CREs (NEGATIVE CONTROL):")
signals_excitatory = extract_signals(cres_excitatory, CONDITIONS)

# DEBUG: Compare signal matrices
print(f"\n{'='*80}")
print("DEBUG: Signal Matrix Comparison")
print("-"*80)
for cond in CONDITIONS.keys():
    if cond in signals_gaba and cond in signals_excitatory:
        gaba_mean = np.mean(signals_gaba[cond])
        excit_mean = np.mean(signals_excitatory[cond])
        gaba_max = np.max(signals_gaba[cond])
        excit_max = np.max(signals_excitatory[cond])

        print(f"\n{cond}:")
        print(f"  GABA:  mean={gaba_mean:.6f}, max={gaba_max:.6f}, shape={signals_gaba[cond].shape}")
        print(f"  Excit: mean={excit_mean:.6f}, max={excit_max:.6f}, shape={signals_excitatory[cond].shape}")

        if abs(gaba_mean - excit_mean) < 0.001:
            print(f"  WARNING: Means are nearly identical!")
        else:
            print(f"  ✓ Means differ by {abs(gaba_mean - excit_mean):.6f}")

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

percentile_90 = np.percentile(all_values, 90)
vmax = percentile_90 * 1.2  # 20% buffer

print(f"90th percentile: {percentile_90:.4f}")
print(f"Common scale: 0 to {vmax:.4f}")

# ============================================================================
# Create Heatmaps
# ============================================================================
print(f"\n{'='*80}")
print("STEP 4: Creating heatmaps...")
print("-"*80)

# GABA heatmap
# NOTE: 3 conditions (Emx1-Ctrl excluded)
fig, axes = plt.subplots(1, 3, figsize=(12, 12), sharey=True)
fig.suptitle(f'POSITIVE CONTROL: ATAC Signal at GABA CREs (n={len(cres_gaba)})',
             fontsize=16, fontweight='bold')

for idx, (cond, matrix) in enumerate(signals_gaba.items()):
    ax = axes[idx]
    im = ax.imshow(matrix, aspect='auto', cmap='Reds',
                   vmin=0, vmax=vmax, interpolation='nearest')
    ax.set_title(cond, fontsize=14, fontweight='bold')
    if idx == 0:
        ax.set_ylabel('CREs', fontsize=12)
    plt.colorbar(im, ax=ax, label='ATAC Signal')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "heatmap_GABA_all_conditions.png"), dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: heatmap_GABA_all_conditions.png")

# Excitatory heatmap
# NOTE: 3 conditions (Emx1-Ctrl excluded)
fig, axes = plt.subplots(1, 3, figsize=(12, 12), sharey=True)
fig.suptitle(f'NEGATIVE CONTROL: ATAC Signal at Excitatory CREs (n={len(cres_excitatory)})',
             fontsize=16, fontweight='bold')

for idx, (cond, matrix) in enumerate(signals_excitatory.items()):
    ax = axes[idx]
    im = ax.imshow(matrix, aspect='auto', cmap='Reds',
                   vmin=0, vmax=vmax, interpolation='nearest')
    ax.set_title(cond, fontsize=14, fontweight='bold')
    if idx == 0:
        ax.set_ylabel('CREs', fontsize=12)
    plt.colorbar(im, ax=ax, label='ATAC Signal')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "heatmap_Excitatory_all_conditions.png"), dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: heatmap_Excitatory_all_conditions.png")

# ============================================================================
# Create Metaprofiles
# ============================================================================
print(f"\n{'='*80}")
print("STEP 5: Creating metaprofiles...")
print("-"*80)

x_coords = np.linspace(-WINDOW_SIZE, WINDOW_SIZE, N_BINS)

# GABA metaprofile
fig, ax = plt.subplots(figsize=(10, 6))
for cond, matrix in signals_gaba.items():
    mean_sig = np.mean(matrix, axis=0)
    sem_sig = stats.sem(matrix, axis=0)
    ax.plot(x_coords, mean_sig, color=COLORS[cond], linewidth=2.5, label=cond)
    ax.fill_between(x_coords, mean_sig - sem_sig, mean_sig + sem_sig,
                     color=COLORS[cond], alpha=0.2)

ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Distance from CRE Center (bp)', fontsize=12)
ax.set_ylabel('Mean ATAC Signal', fontsize=12)
ax.set_ylim(0, vmax)
ax.set_title(f'POSITIVE CONTROL: Metaprofile at GABA CREs (n={len(cres_gaba)})',
             fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "metaprofile_GABA_all_CREs.png"), dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: metaprofile_GABA_all_CREs.png")

# Excitatory metaprofile
fig, ax = plt.subplots(figsize=(10, 6))
for cond, matrix in signals_excitatory.items():
    mean_sig = np.mean(matrix, axis=0)
    sem_sig = stats.sem(matrix, axis=0)
    ax.plot(x_coords, mean_sig, color=COLORS[cond], linewidth=2.5, label=cond)
    ax.fill_between(x_coords, mean_sig - sem_sig, mean_sig + sem_sig,
                     color=COLORS[cond], alpha=0.2)

ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Distance from CRE Center (bp)', fontsize=12)
ax.set_ylabel('Mean ATAC Signal', fontsize=12)
ax.set_ylim(0, vmax)
ax.set_title(f'NEGATIVE CONTROL: Metaprofile at Excitatory CREs (n={len(cres_excitatory)})',
             fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "metaprofile_Excitatory_all_CREs.png"), dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: metaprofile_Excitatory_all_CREs.png")

# ============================================================================
# Create Comparison Plot
# ============================================================================
print(f"\n{'='*80}")
print("STEP 6: Creating comparison plot...")
print("-"*80)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# GABA metaprofiles
ax = axes[0, 0]
for cond, matrix in signals_gaba.items():
    mean_sig = np.mean(matrix, axis=0)
    ax.plot(x_coords, mean_sig, color=COLORS[cond], linewidth=2, label=cond)
ax.set_title('GABA CREs (Positive Control)', fontweight='bold')
ax.set_ylabel('Mean ATAC Signal')
ax.set_ylim(0, vmax)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Excitatory metaprofiles
ax = axes[0, 1]
for cond, matrix in signals_excitatory.items():
    mean_sig = np.mean(matrix, axis=0)
    ax.plot(x_coords, mean_sig, color=COLORS[cond], linewidth=2, label=cond)
ax.set_title('Excitatory CREs (Negative Control)', fontweight='bold')
ax.set_ylim(0, vmax)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Mean signal comparison
ax = axes[1, 0]
gaba_means = [np.mean(m) for m in signals_gaba.values()]
excit_means = [np.mean(m) for m in signals_excitatory.values()]
x_pos = np.arange(len(CONDITIONS))
width = 0.35
ax.bar(x_pos - width/2, gaba_means, width, label='GABA CREs', color='red', alpha=0.7)
ax.bar(x_pos + width/2, excit_means, width, label='Excitatory CREs', color='blue', alpha=0.7)
ax.set_xticks(x_pos)
ax.set_xticklabels(CONDITIONS.keys(), rotation=45, ha='right')
ax.set_ylabel('Mean ATAC Signal')
ax.set_title('Signal Comparison by Condition', fontweight='bold')
ax.legend()
ax.grid(True, axis='y', alpha=0.3)

# Fold enrichment
ax = axes[1, 1]
fold_enrichments = [g/e if e > 0 else 0 for g, e in zip(gaba_means, excit_means)]
ax.bar(x_pos, fold_enrichments, color='purple', alpha=0.7)
ax.axhline(1, color='red', linestyle='--', label='No enrichment')
ax.set_xticks(x_pos)
ax.set_xticklabels(CONDITIONS.keys(), rotation=45, ha='right')
ax.set_ylabel('Fold Enrichment (GABA / Excitatory)')
ax.set_title('Cell-Type Specificity', fontweight='bold')
ax.legend()
ax.grid(True, axis='y', alpha=0.3)

plt.suptitle('Cell-Type Specificity: GABA vs Excitatory CREs', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "comparison_GABA_vs_Excitatory.png"), dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: comparison_GABA_vs_Excitatory.png")

# ============================================================================
# Save Statistics
# ============================================================================
print(f"\n{'='*80}")
print("STEP 7: Saving statistics...")
print("-"*80)

stats_file = os.path.join(OUTPUT_DIR, "enrichment_statistics.txt")
with open(stats_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("CELL-TYPE SPECIFICITY: GABA vs Excitatory CREs\n")
    f.write("="*80 + "\n\n")
    f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write(f"GABA CREs: {len(cres_gaba)}\n")
    f.write(f"Excitatory CREs: {len(cres_excitatory)}\n")
    f.write(f"Common scale: 0 to {vmax:.4f} (90th percentile + 20%)\n\n")

    f.write("MEAN SIGNAL COMPARISON:\n")
    f.write("-"*80 + "\n")
    for idx, cond in enumerate(CONDITIONS.keys()):
        if cond in signals_gaba and cond in signals_excitatory:
            gaba_mean = np.mean(signals_gaba[cond])
            excit_mean = np.mean(signals_excitatory[cond])
            fold = gaba_mean / excit_mean if excit_mean > 0 else 0
            f.write(f"{cond}:\n")
            f.write(f"  GABA mean: {gaba_mean:.6f}\n")
            f.write(f"  Excitatory mean: {excit_mean:.6f}\n")
            f.write(f"  Fold enrichment: {fold:.2f}x\n\n")

    f.write("INTERPRETATION:\n")
    f.write("-"*80 + "\n")
    f.write("Fold enrichment >3x indicates strong cell-type specificity.\n")
    f.write("GABA ATAC samples should show HIGH signal at GABA CREs and\n")
    f.write("LOW signal at Excitatory CREs, demonstrating specificity.\n")

print("✓ Saved: enrichment_statistics.txt")

print(f"\n{'='*80}")
print("ANALYSIS COMPLETE!")
print("="*80)
print(f"\nOutput directory: {OUTPUT_DIR}/")
print(f"\nKey files:")
print(f"  - heatmap_GABA_all_conditions.png (POSITIVE)")
print(f"  - heatmap_Excitatory_all_conditions.png (NEGATIVE)")
print(f"  - metaprofile_GABA_all_CREs.png (POSITIVE)")
print(f"  - metaprofile_Excitatory_all_CREs.png (NEGATIVE)")
print(f"  - comparison_GABA_vs_Excitatory.png (COMPARISON)")
print(f"  - enrichment_statistics.txt (QUANTITATIVE)")
print(f"\nAll plots use identical scale (0-{vmax:.4f}) for direct comparison")
print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
