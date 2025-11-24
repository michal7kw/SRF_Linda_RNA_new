# %%
%%capture
import statsmodels.api as sm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#import ssl; ssl._create_default_https_context = ssl._create_unverified_context
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline

from sciRED import ensembleFCA as efca
from sciRED import glm
from sciRED import rotations as rot
from sciRED import metrics as met

from sciRED.utils import preprocess as proc
from sciRED.utils import visualize as vis
from sciRED.utils import corr
from sciRED.examples import ex_preprocess as exproc
from sciRED.examples import ex_visualize as exvis

# %%
import os

PROJECT_DIR = "D:/Github/SRF_Linda_RNA"
WORKING_DIR = os.path.join(PROJECT_DIR, "combine_data/ANALYSIS/Enrichment/sciRED")
DATA_DIR = os.path.join(PROJECT_DIR, "combine_data", "results_from_raw")

os.chdir(WORKING_DIR)

# %%
np.random.seed(10)
NUM_COMPONENTS = 30
NUM_GENES = 2000
NUM_COMP_TO_VIS = 3

# %%
%%capture
data_file_path = os.path.join(DATA_DIR, "annotation_final.h5ad")
data = exproc.import_AnnData(data_file_path)

# %%
print(f"Before filtering: {data.shape}")
data = data[:, ~data.var_names.str.startswith("mt-")].copy()
print(f"After filtering: {data.shape}")

# %%
data

# %%
print(data.obs["genotype"].unique())
print(data.obs["condition"].unique())
print(data.obs["cell_type_L2_new"].unique())

# %%
data = data[data.obs["genotype"]=="Emx1"]
print(data.shape)
data = data[data.obs["cell_type_L2_new"]=="Mature GC"]
print(data.shape)

# Let's check the distribution of conditions after filtering to be sure both are present
print("\nValue counts for 'condition' after filtering:")
print(data.obs['condition'].value_counts())

# %%
data, gene_idx = proc.get_sub_data(data, num_genes=NUM_GENES) # subset the data to num_genes HVGs
y, genes, num_cells, num_genes = proc.get_data_array(data)

# %%
data

# %%
print(f"gene_idx: {gene_idx[:10]} \n")
print(f"y[:5][:5] {y[:5][:5]}, y.shape: {y.shape} \n")
print(f"genes[:5] {genes[:5]} \n")
print(f"num_cells: {num_cells}\n") 
print(f"num_genes: {num_genes} \n")

# %% [markdown]
# **Step 1: Factor discovery:**
# 
# One option for the design matrix is to only regress out library size: 
# ```x = proc.get_library_design_mat(data, lib_size='total_counts')```

# %%
data.obs["total_counts"].head()

# %%
#### Design matrix - including library size
x = data.obs["total_counts"]
x = sm.add_constant(x) ## adding the intercept

# %%
print(x.shape)
print(x[:5])

# %%
glm_fit_dict = glm.poissonGLM(y, x)
resid_pearson = glm_fit_dict['resid_pearson'] 
print('pearson residuals: ', resid_pearson.shape) # numpy array of shape (num_genes, num_cells)
print('y shape: ', y.shape) # (num_cells, num_genes)
y = resid_pearson.T # (num_cells, num_genes)
print('y shape: ', y.shape) # (num_cells, num_genes)

# %% [markdown]
# We then apply PCA to the extracted residuals. PCA factors are then rotated (varimax or promax) to improve interpretibility. 

# %%
################################################
#### Running PCA on the pearson residual ######
################################################

### using pipeline to scale the gene expression data first
pipeline = Pipeline([('scaling', StandardScaler()), ('pca', PCA(n_components=NUM_COMPONENTS))])
pca_scores = pipeline.fit_transform(y)
pca = pipeline.named_steps['pca']
pca_loading = pca.components_
pca_loading.shape
plt.plot(pca.explained_variance_ratio_)

# %%
import matplotlib.pyplot as plt

# Create a colors vector from the 'genotype' column:
# 1. Get the unique genotype values.
unique_condition= data.obs["condition"].unique()

# 2. Use a discrete colormap (here 'tab10') to assign a distinct color to each genotype.
cmap = plt.get_cmap("tab10")
condition_color_dict = {condition: cmap(i % cmap.N) for i, condition in enumerate(unique_condition)}

# 3. Map the genotype column to its corresponding colors using a list comprehension to avoid MultiIndex issues.
condition_colors = [condition_color_dict[con] for con in data.obs["condition"]]

# %%
plt_legend_condition= exvis.get_legend_patch(data.obs["condition"], condition_colors )

# %%
title = 'PCA of pearson residuals - lib size removed'
### make a dictionary of colors for each sample in y_sample
vis.plot_pca(pca_scores, NUM_COMP_TO_VIS, 
               cell_color_vec= condition_colors, 
               legend_handles=True,
               title=title,
               plt_legend_list=plt_legend_condition)    

# %%
#### plot the loadings of the factors
vis.plot_factor_loading(pca_loading.T, genes, 0, 1, fontsize=10, 
                    num_gene_labels=2,
                    title='Scatter plot of the loading vectors', 
                    label_x=True, label_y=True)

# %%
vis.plot_umap(pca_scores, 
              title='UMAP',
              cell_color_vec= condition_colors, 
               legend_handles=True, plt_legend_list=plt_legend_condition)

# %%
######## Applying varimax rotation to the factor scores
rotation_results_varimax = rot.varimax(pca_loading.T)
varimax_loading = rotation_results_varimax['rotloading']
pca_scores_varimax = rot.get_rotated_scores(pca_scores, rotation_results_varimax['rotmat'])

# %%
print(f"rotation_results_varimax.keys(): {rotation_results_varimax.keys()}")

print(f"varimax_loading.shape: {varimax_loading.shape}")
print(f"pca_scores_varimax.shape: {pca_scores_varimax.shape}")

# %%
title = 'Varimax PCA of pearson residuals '

vis.plot_pca(pca_scores_varimax, NUM_COMP_TO_VIS+1, 
               cell_color_vec= condition_colors, 
               legend_handles=True,
               title=title,
               plt_legend_list=plt_legend_condition)

# %%
##### Setting the factor scores an loadings to be used in step-2 based on Varimax factors
factor_loading = rotation_results_varimax['rotloading']
factor_scores = pca_scores_varimax

# %%
print(factor_scores.shape)
print(factor_scores[:2][:2])

# %% [markdown]
# **Step 2: Factor-Covariate Association**: To identify factors that explain a specific covariate, sciRED employs an ensemble classifier as a second step. We apply four machine learning classifiers to predict covariate labels based on the cell-specific factor weights. Feature importance scores are obtained from each classifier are then scaled based on one out of three scaling methods, and averaged to generate a consensus association score. The consensus scores for every combination of covariate and factor are aggregated into the factor-covariate association table (FCAT) and visualized in a heatmap. 
# The FCAT function takes-in the cell-level labels for each covariate. The resulting tables for each covariate are then concatenated and visualized as a heatmap.  

# %%
print(data.obs["condition"].shape)
print(type(data.obs["condition"]))

# %%
%%capture
####################################
#### FCAT score calculation ######
####################################

### FCAT needs to be calculated for each covariate separately
# Create a lowercase version of the labels. The sciRED library likely
# gives special treatment to a 'control' level (lowercase) as a baseline.
fcat_labels = data.obs['condition'].str.lower()
fcat_condition = efca.FCAT(fcat_labels, factor_scores, scale='standard', mean='arithmatic')

### concatenate FCAT table for protocol and cell line
fcat = fcat_condition

# %%
### visualize the first 15 factors
vis.plot_FCAT(fcat.iloc[:,0:15],title='', color='coolwarm',x_axis_fontsize=35, 
              y_axis_fontsize=35, title_fontsize=35,
              x_axis_tick_fontsize=32, y_axis_tick_fontsize=34)

# %% [markdown]
# **Step 2a: Visualize Factor Directionality**
# 
# Since the FCAT shows identical rows for the binary 'condition' covariate, we need to check the direction of the effect for our factor of interest (Factor 7). This boxplot shows the distribution of factor scores for each condition.

# %%
INTERESTING_FACTOR_ID = 0

# Create a DataFrame for easy plotting with seaborn
plot_df = pd.DataFrame({
    f'Factor_{INTERESTING_FACTOR_ID + 1}_Scores': factor_scores[:, INTERESTING_FACTOR_ID],
    'Condition': data.obs['condition'].values
})

# Create the boxplot
plt.figure(figsize=(6, 6))
sns.boxplot(x='Condition', y=f'Factor_{INTERESTING_FACTOR_ID + 1}_Scores', data=plot_df)
# Add individual data points for better visualization
sns.stripplot(x='Condition', y=f'Factor_{INTERESTING_FACTOR_ID + 1}_Scores', data=plot_df, color=".25", size=3)
plt.title(f'Distribution of Scores for Factor {INTERESTING_FACTOR_ID + 1} by Condition')
plt.ylabel('Factor Score')
plt.xlabel('Condition')
plt.show()

# %% [markdown]
# Significant vs non-significant associations between factors and covariates are determined using a threshold automatically obtained using Otsu’s method. This threshold can assist users in defining the number of inferred factors (K). For example, if a considerable proportion of factors fail to align with any covariates, it may prompt the user to reduce K. 

# %%

## rownames of the FCAT table
all_covariate_levels = fcat.index.values

### using Otsu's method to calculate the threshold
threshold = efca.get_otsu_threshold(fcat.values.flatten())

vis.plot_histogram(fcat.values.flatten(),
                   xlabel='Factor-Covariate Association scores',
                   title='FCAT score distribution',
                   threshold=threshold)



# %%

matched_factor_dist, percent_matched_fact = efca.get_percent_matched_factors(fcat, threshold)
matched_covariate_dist, percent_matched_cov = efca.get_percent_matched_covariates(fcat, threshold=threshold)

print('percent_matched_fact: ', percent_matched_fact)
print('percent_matched_cov: ', percent_matched_cov)
vis.plot_matched_factor_dist(matched_factor_dist)
vis.plot_matched_covariate_dist(matched_covariate_dist, 
                                covariate_levels=all_covariate_levels)



# %% [markdown]
# We can check the correlation between the factors and the library size which was regressed out

# %%
factor_libsize_correlation = corr.get_factor_libsize_correlation(factor_scores, library_size = data.obs["total_counts"])
vis.plot_factor_cor_barplot(factor_libsize_correlation, 
             title='Correlation of factors with library size', 
             y_label='Correlation', x_label='Factors')

# %% [markdown]
# **Step 3: Interpretability scores:** The third step of sciRED involves quantifying the interpretability of identified factors. We defined four categories of metrics: separability, effect size, specificity, and homogeneity which are presented as the FIST table. 

# %%
%%capture
####################################
#### Bimodality scores
silhouette_score = met.kmeans_bimodal_score(factor_scores, time_eff=True)
bimodality_index = met.bimodality_index(factor_scores)
bimodality_score = np.mean([silhouette_score, bimodality_index], axis=0)

#### Effect size
factor_variance = met.factor_variance(factor_scores)

## Specificity
simpson_fcat = met.simpson_diversity_index(fcat)

### Homogeneity (how well-mixed the factors are given the covariate)
asv_condition = met.average_scaled_var(factor_scores, data.obs["condition"], mean_type='arithmetic')

#### plot the ralative variance table
svt_condition = met.scaled_var_table(factor_scores, data.obs["condition"])
svt = pd.concat([svt_condition], axis=0)
vis.plot_relativeVar(svt.iloc[:,0:15], title='Relative variance score table')

# %%
########### create factor-interpretibility score table (FIST) ######
metrics_dict = {'Bimodality':bimodality_score, 
                    'Specificity':simpson_fcat,
                    'Effect size': factor_variance,
                    'Homogeneity (condition)':asv_condition}

fist = met.FIST(metrics_dict)
### subset the first 15 factors of fist dataframe
vis.plot_FIST(fist.iloc[0:15,:])

# %% [markdown]
# **Step 4: Biological Interpretation of a Selected Factor**
# 
# After running the analysis, inspect the FCAT heatmap generated in Step 2. Identify a factor that shows a high association score with your 'condition' covariate (i.e., differentiates Control vs. Mutant). Let's assume Factor `k` is the one you're interested in.
# 
# The code below will extract and display the top genes for that factor. These genes represent the biological program captured by the factor. You can then perform pathway analysis on this gene list.

# %%
# Replace with the factor number you identified from the FCAT heatmap
INTERESTING_FACTOR_ID = 0

# Get the loadings for the factor of interest
factor_loadings_for_factor_k = factor_loading[:, INTERESTING_FACTOR_ID]

# Create a pandas Series for easy sorting and viewing
gene_loadings = pd.Series(factor_loadings_for_factor_k, index=genes)

# Sort genes by their loading values
sorted_gene_loadings = gene_loadings.sort_values(ascending=False)

print(f"--- Top genes for Factor {INTERESTING_FACTOR_ID + 1} ---")

print("\n--- Top 20 Positive-Loading Genes (up-regulated in one condition) ---")
print(sorted_gene_loadings.head(20))

print("\n--- Top 20 Negative-Loading Genes (up-regulated in the other condition) ---")
print(sorted_gene_loadings.tail(20))

# You can then take these gene lists for GO/pathway enrichment analysis.

# %%
gene_of_interest = "Calm2"

# Extract expression for the gene (ensure it's a numpy array, handling sparse matrices if necessary)
expr = data[:, gene_of_interest].X
if hasattr(expr, "toarray"):
    expr = expr.toarray().flatten()
else:
    expr = np.array(expr).flatten()

umap_coords = data.obsm["X_umap"]

# Plot UMAP colored by gene expression
plt.figure(figsize=(8, 6))
sc = plt.scatter(umap_coords[:, 0], umap_coords[:, 1], c=expr, cmap='Reds', s=10)
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title(f"UMAP Projection Colored by {gene_of_interest} Expression")
plt.colorbar(sc, label=f"{gene_of_interest} Expression")
plt.show()

# Plot UMAP colored by condition
plt.figure(figsize=(8, 6))
conditions = data.obs["condition"]
for cond in conditions.unique():
    idx = conditions == cond
    plt.scatter(umap_coords[idx, 0], umap_coords[idx, 1], s=10, label=cond)
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title("UMAP Projection Colored by Condition")
plt.legend(title="Condition")
plt.show()

# %%



