# %% [markdown]
# # Environment Setup and Data Loading

# %%
import scanpy as sc
import numpy as np
import pandas as pd
import os
import sys
import anndata as ad
import random
import importlib

PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = os.path.join(PROJECT_DIR, "combine_data")
LIB_DIRS = os.path.join(WORKING_DIR, "LIBRARIES")
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)
sys.path.insert(0, LIB_DIRS)

import functions_degs

importlib.reload(functions_degs)

# Set seeds for all random number generators
random_seed = 0
np.random.seed(random_seed)
random.seed(random_seed)

# %%
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
adata_path = os.path.join(INPUT_DIR, 'annotation_final.h5ad')

PARENT_OUTPUT_DIR = os.path.join(INPUT_DIR, "DEGs_cell_type_L2_new")

# CUSTOM_ANALYSIS =  None
CUSTOM_ANALYSIS =  "FC_0_25"

DEG_BY = "cell_type_L2_new"

if CUSTOM_ANALYSIS is not None:
    PARENT_OUTPUT_DIR = PARENT_OUTPUT_DIR + CUSTOM_ANALYSIS

PLOT_OUTPUT_DIR = os.path.join(PARENT_OUTPUT_DIR, 'plots')
DGE_OUTPUT_DIR = os.path.join(PARENT_OUTPUT_DIR, 'dge_res')
BIOMARKER_OUTPUT_DIR = os.path.join(PARENT_OUTPUT_DIR, 'biomarkers')

# Create output directories if they don't exist
os.makedirs(PLOT_OUTPUT_DIR, exist_ok=True)
os.makedirs(DGE_OUTPUT_DIR, exist_ok=True)
os.makedirs(BIOMARKER_OUTPUT_DIR, exist_ok=True)

print(f"Input directory: {INPUT_DIR}")
print(f"Plot output directory (Overall/Genotype): {PLOT_OUTPUT_DIR}")
print(f"DGE output directory (Overall/Genotype): {DGE_OUTPUT_DIR}")
print(f"Biomarker output directory: {BIOMARKER_OUTPUT_DIR}")

# Configure scanpy settings
sc.settings.figdir = PLOT_OUTPUT_DIR
sc.settings.set_figure_params(dpi=150, facecolor='white')

# %%
# Load data
print(f"\nLoading dataset from {adata_path}")
adata = sc.read_h5ad(adata_path)


# %% [markdown]

# %% [markdown]
# # Visualize Apoptosis Pathway Gene Expression

# %%
# Gene set for WP_APOPTOSIS (Mus musculus) from WikiPathways WP1254
# https://www.wikipathways.org/pathways/WP1254

with open('WP_APOPTOSIS.v2025.1.Mm.grp', 'r') as f:
    # Skip the first line (pathway name) and read the rest of the genes
    apoptosis_genes = [line.strip() for line in f.readlines()[1:] if line.strip()]

# Check which genes are present in the anndata object
apoptosis_genes_in_adata = [gene for gene in apoptosis_genes if gene in adata.var_names]
print(f"Found {len(apoptosis_genes_in_adata)} out of {len(apoptosis_genes)} apoptosis genes in the dataset.")
print(apoptosis_genes)

# %%
# Calculate the score for the apoptosis gene set
ctrl_size = min(50, len(apoptosis_genes_in_adata) - 1)
sc.tl.score_genes(adata, gene_list=apoptosis_genes_in_adata, score_name='apoptosis_score', ctrl_size=ctrl_size)
# Normalize and scale the apoptosis score (z-score normalization)
score = adata.obs["apoptosis_score"]
# adata.obs["apoptosis_score"] = (score - score.mean()) / score.std()
# adata.obs["apoptosis_score"] = np.sign(score) * np.power(score, 4)
max_score = adata.obs.apoptosis_score.max()
min_score = adata.obs.apoptosis_score.min()
print(f"max apoptosis score: {max_score}")
print(f"min apoptosis score: {min_score}")

# %%
import matplotlib.pyplot as plt
plt.hist(adata.obs.apoptosis_score, bins=30, edgecolor='black')
plt.xlabel("Apoptosis Score")
plt.ylabel("Frequency")
plt.title("Histogram of Apoptosis Score")
plt.show()

# %%
# Plot the score on the UMAP for Control cells
sc.pl.umap(
    adata[adata.obs.condition == "Control"],
    color='apoptosis_score',
    cmap='viridis',
    vmin=min_score,
    vmax=max_score,
    title='Apoptosis Pathway Score (Control)'
)

# Plot the score on the UMAP for Mutant cells
sc.pl.umap(
    adata[adata.obs.condition == "Mutant"],
    color='apoptosis_score',
    cmap='viridis',
    vmin=min_score,
    vmax=max_score,
    title='Apoptosis Pathway Score (Mutant)'
)

# %%
# Compare the apoptosis score between Control and Mutant in Mature GC cell types by plotting them side‐by‐side.
gc_mature_data = adata[adata.obs.cell_type_L2_new == "Mature GC"]

import matplotlib.pyplot as plt

# Create subplots for Control and Mutant conditions
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# Filter the GC_Mature cells by condition
control_data = gc_mature_data[gc_mature_data.obs.condition == "Control"]
mutant_data = gc_mature_data[gc_mature_data.obs.condition == "Mutant"]

# Plot UMAP for Control cells
sc.pl.umap(
    control_data,
    color='apoptosis_score',
    cmap='viridis',
    vmin=min_score,
    vmax=max_score,
    title='Apoptosis Score in Mature GC (Control)',
    ax=axes[0],
    show=False
)

# Plot UMAP for Mutant cells
sc.pl.umap(
    mutant_data,
    color='apoptosis_score',
    cmap='viridis',
    vmin=min_score,
    vmax=max_score,
    title='Apoptosis Score in Mature GC (Mutant)',
    ax=axes[1],
    show=False
)

plt.tight_layout()
plt.savefig('apoptosis_score_gc_mature.png')
plt.show()

# %%
sc.pl.violin(adata, 'apoptosis_score', groupby='condition', save='apoptosis_score_violin.png')
