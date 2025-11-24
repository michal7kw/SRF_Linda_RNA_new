# %%
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import anndata as ad
import warnings
from celltypist import models, annotate
import seaborn as sns
from scipy.sparse import csr_matrix
import matplotlib.patches as mpatches

# PROJECT_DIR = "/home/michal/Github/SRF_Linda_RNA"
PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = f"{PROJECT_DIR}/combine_data"
os.chdir(WORKING_DIR)
sys.path.insert(0, WORKING_DIR)

import random

# Set seeds for all random number generators
random_seed = 0
np.random.seed(random_seed)
random.seed(random_seed)
    
# %%
# Set up directories
INPUT_DIR = f"{WORKING_DIR}/results/final_annotation"
adata_path = os.path.join(INPUT_DIR, 'all_samples_annotated.h5ad')
print(f"\nLoading merged dataset from {adata_path}")
adata_1 = sc.read_h5ad(adata_path)

# %%
# Set up directories
INPUT_DIR = f"{WORKING_DIR}/results_from_raw/final_annotation"
adata_path = os.path.join(INPUT_DIR, 'merged_raw_final_annotated.h5ad')
print(f"\nLoading merged dataset from {adata_path}")
adata_2 = sc.read_h5ad(adata_path)

# %%
adata_1
# %%
adata_2

# %%
print(adata_1.layers)
print(adata_2.layers)

# %%
print(adata_1.X[:5][:5])
print(adata_2.X[:5][:5])
# %%
