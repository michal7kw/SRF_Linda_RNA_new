# %%
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import anndata as ad
import warnings
import seaborn as sns
import matplotlib.patches as mpatches

# PROJECT_DIR = "/home/michal/Github/SRF_Linda_RNA"
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)
    
# %%
INPUT_DIR = f"{WORKING_DIR}/results_from_raw/final_annotation"
adata_path = os.path.join(INPUT_DIR, 'merged_raw_final_annotated_simple.h5ad')
print(f"\nLoading merged dataset from {adata_path}")
adata = sc.read_h5ad(adata_path)

# %%
adata

# %%
# Check for specific genes
target_genes = ["Xist", "Tsix", "Eif2s3x", "Ddx3x", "Kdm5c", "Kdm6a", "Rps6ka3", "Uty", "Ddx3y", "Eif2s3y", "Kdm5d", "Rps4y1", "Jarid1d"]
adata_genes = adata.var_names.tolist()

print("\nChecking for target genes:")
for gene in target_genes:
    if gene in adata_genes:
        print(f"- Found: {gene}")
    else:
        print(f"- Not found: {gene}")

# %%
# %%
INPUT_DIR = os.path.join(WORKING_DIR, "results_from_raw", "celltypist_annotations") # Use output from script 1
adata_path = os.path.join(INPUT_DIR, 'merged_raw_annotated.h5ad') # Input is the annotated file from script 1
print(f"\nLoading annotated dataset from {adata_path}")
adata = sc.read_h5ad(adata_path)

# %%
adata

# %%
# Check for specific genes
target_genes = ["Xist", "Tsix", "Eif2s3x", "Ddx3x", "Kdm5c", "Kdm6a", "Rps6ka3", "Uty", "Ddx3y", "Eif2s3y", "Kdm5d", "Rps4y1", "Jarid1d"]
adata_genes = adata.var_names.tolist()

print("\nChecking for target genes:")
for gene in target_genes:
    if gene in adata_genes:
        print(f"- Found: {gene}")
    else:
        print(f"- Not found: {gene}")

# %%
