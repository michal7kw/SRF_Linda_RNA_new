# %% [markdown]
# # Env

# %%
%%capture
import statsmodels.api as sm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import warnings

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

from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.inspection import permutation_importance

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

# %% [markdown]
# # Load Data

# %%
%%capture
data_file_path = os.path.join(DATA_DIR, "annotation_final.h5ad")
data = exproc.import_AnnData(data_file_path)

# %%
data

# %%
# print(data.obs["genotype"].unique())
# print(data.obs["condition"].unique())

# %%
data = data[data.obs["cell_type_L2_new"]=="Mature GC"]
print(data.shape)

# %%
data, gene_idx = proc.get_sub_data(data, num_genes=NUM_GENES) # subset the data to num_genes HVGs
y, genes, num_cells, num_genes = proc.get_data_array(data)

# %%
# data

# %%
# print(f"gene_idx: {gene_idx[:10]} \n")
# print(f"y[:5][:5] {y[:5][:5]}, y.shape: {y.shape} \n")
# print(f"genes[:5] {genes[:5]} \n")
# print(f"num_cells: {num_cells}\n") 
# print(f"num_genes: {num_genes} \n")

# %% [markdown]
# # **Step 1: Factor discovery:**

# %%
data.obs["total_counts"].head()

# %%
#### Design matrix - including library size
x = data.obs["total_counts"]
x = sm.add_constant(x) ## adding the intercept
print(x[:5])

# %%
glm_fit_dict = glm.poissonGLM(y, x)
resid_pearson = glm_fit_dict['resid_pearson'] 
print('pearson residuals: ', resid_pearson.shape)
y = resid_pearson.T 
print('y shape: ', y.shape)

# %% [markdown]
# Apply PCA to the extracted residuals. PCA factors are then rotated (varimax or promax) to improve interpretibility. 

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
_ = plt.plot(pca.explained_variance_ratio_)

# %%
unique_condition= data.obs["condition"].unique()
cmap = plt.get_cmap("tab10")
condition_color_dict = {condition: cmap(i % cmap.N) for i, condition in enumerate(unique_condition)}
condition_colors = [condition_color_dict[con] for con in data.obs["condition"]]

# %%
unique_genotype= data.obs["genotype"].unique()
cmap = plt.get_cmap("tab10")
genotype_color_dict = {genotype: cmap(i % cmap.N) for i, genotype in enumerate(unique_genotype)}
genotype_colors = [genotype_color_dict[genotype] for genotype in data.obs["genotype"]]

# %%
plt_legend_genotype = exvis.get_legend_patch(data.obs["genotype"], genotype_colors )
plt_legend_condition= exvis.get_legend_patch(data.obs["condition"], condition_colors )

# %%
title = 'PCA of pearson residuals - lib size removed'
vis.plot_pca(pca_scores, NUM_COMP_TO_VIS, 
               cell_color_vec= genotype_colors, 
               legend_handles=True,
               title=title,
               plt_legend_list=plt_legend_genotype)

# %%
title = 'PCA of pearson residuals - lib size removed'
vis.plot_pca(pca_scores, NUM_COMP_TO_VIS, 
               cell_color_vec= condition_colors, 
               legend_handles=True,
               title=title,
               plt_legend_list=plt_legend_condition)    

# %%
#### plot the loadings of the factors
vis.plot_factor_loading(pca_loading.T, genes, 0, 2, fontsize=10, 
                    num_gene_labels=2,
                    title='Scatter plot of the loading vectors', 
                    label_x=True, label_y=True)

# %%
vis.plot_umap(pca_scores, 
              title='UMAP',
              cell_color_vec= condition_colors, 
               legend_handles=True, plt_legend_list=plt_legend_condition)

vis.plot_umap(pca_scores, 
              title='UMAP',
              cell_color_vec= genotype_colors, 
               legend_handles=True, plt_legend_list=plt_legend_genotype)

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
               cell_color_vec= genotype_colors, 
               legend_handles=True,
               title=title,
               plt_legend_list=plt_legend_genotype)

# %%
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
# print(factor_scores.shape)
# print(factor_scores[:2][:2])

# %% [markdown]
# # **Step 2: Factor-Covariate Association**:
# 
# - To identify factors that explain a specific covariate, sciRED employs an ensemble classifier as a second step. 
# - Apply four machine learning classifiers to predict covariate labels based on the cell-specific factor weights. 
# - Feature importance scores are obtained from each classifier are then scaled based on one out of three scaling methods, and averaged to generate a consensus association score. 
# - The consensus scores for every combination of covariate and factor are aggregated into the factor-covariate association table (FCAT) and visualized in a heatmap. 
# - The FCAT function takes-in the cell-level labels for each covariate. 
# - The resulting tables for each covariate are then concatenated and visualized as a heatmap.  

# %%
# print(data.obs["condition"].shape)
# print(type(data.obs["condition"]))

# print(data.obs["genotype"].shape)
# print(type(data.obs["genotype"]))

# %%
%%capture
####################################
#### FCAT score calculation ######
####################################

### FCAT needs to be calculated for each covariate separately
fcat_condition = efca.FCAT(data.obs["condition"], factor_scores, scale='standard', mean='arithmatic')
fcat_genotype = efca.FCAT(data.obs["genotype"], factor_scores, scale='standard', mean='arithmatic')

### concatenate FCAT table for protocol and cell line
fcat = pd.concat([fcat_condition, fcat_genotype], axis=0)

# %%
### visualize the first 15 factors
vis.plot_FCAT(fcat.iloc[:,0:15],title='', color='coolwarm',x_axis_fontsize=35, 
              y_axis_fontsize=35, title_fontsize=35,
              x_axis_tick_fontsize=32, y_axis_tick_fontsize=34)

# %% [markdown]
# # **Step 2a: Visualize Factor Directionality**
# 
# With the original implementation FCAT, for the binary 'condition' covariate, we need to check the direction of the effect for the factors of interest.   
# This boxplot shows the distribution of factor scores for each condition.

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

# %%
# Create a DataFrame for easy plotting with seaborn
plot_df = pd.DataFrame({
    f'Factor_{INTERESTING_FACTOR_ID + 1}_Scores': factor_scores[:, INTERESTING_FACTOR_ID],
    'Genotype': data.obs['genotype'].values
})

# Create the boxplot
plt.figure(figsize=(6, 6))
sns.boxplot(x='Genotype', y=f'Factor_{INTERESTING_FACTOR_ID + 1}_Scores', data=plot_df)
# Add individual data points for better visualization
sns.stripplot(x='Genotype', y=f'Factor_{INTERESTING_FACTOR_ID + 1}_Scores', data=plot_df, color=".25", size=3)
plt.title(f'Distribution of Scores for Factor {INTERESTING_FACTOR_ID + 1} by Genotype')
plt.ylabel('Factor Score')
plt.xlabel('Condition')
plt.show()

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
# Define the interesting factor IDs
INTERESTING_FACTOR_IDS = [0, 2, 10, 12]

num_factors = len(INTERESTING_FACTOR_IDS)
# Create a 2-row subplot: first row for Condition, second row for Genotype
fig, axes = plt.subplots(2, num_factors, figsize=(6 * num_factors, 12))

for idx, factor_id in enumerate(INTERESTING_FACTOR_IDS):
    # ----- Plot for Condition -----
    # Create a DataFrame for the current factor using condition data
    plot_df_cond = pd.DataFrame({
        f'Factor_{factor_id + 1}_Scores': factor_scores[:, factor_id],
        'Condition': data.obs['condition'].values
    })
    
    # Calculate mean and standard error for each condition
    summary_stats_cond = plot_df_cond.groupby('Condition', observed=False)[f'Factor_{factor_id + 1}_Scores'].agg(['mean', 'sem']).reset_index()
    
    ax_cond = axes[0, idx]
    bars_cond = ax_cond.bar(summary_stats_cond['Condition'], summary_stats_cond['mean'],
                            yerr=summary_stats_cond['sem'], capsize=5, alpha=0.7, color='steelblue')
    
    # Customize condition plot
    ax_cond.set_title(f'Mean Scores for Factor {factor_id + 1} by Condition', fontsize=14, fontweight='bold')
    ax_cond.set_ylabel('Mean Factor Score', fontsize=12)
    ax_cond.set_xlabel('Condition', fontsize=12)
    ax_cond.grid(axis='y', alpha=0.3)
    
    # Add value labels on top of condition bars
    for bar, mean_val in zip(bars_cond, summary_stats_cond['mean']):
        ax_cond.text(bar.get_x() + bar.get_width()/2, bar.get_height() + summary_stats_cond['sem'].max()*0.1,
                     f'{mean_val:.3f}', ha='center', va='bottom', fontsize=10)
    
    # ----- Plot for Genotype -----
    # Create a DataFrame for the current factor using genotype data
    plot_df_geno = pd.DataFrame({
        f'Factor_{factor_id + 1}_Scores': factor_scores[:, factor_id],
        'Genotype': data.obs['genotype'].values
    })
    
    # Calculate mean and standard error for each genotype
    summary_stats_geno = plot_df_geno.groupby('Genotype', observed=False)[f'Factor_{factor_id + 1}_Scores'].agg(['mean', 'sem']).reset_index()
    
    ax_geno = axes[1, idx]
    bars_geno = ax_geno.bar(summary_stats_geno['Genotype'], summary_stats_geno['mean'],
                            yerr=summary_stats_geno['sem'], capsize=5, alpha=0.7, color='steelblue')
    
    # Customize genotype plot
    ax_geno.set_title(f'Mean Scores for Factor {factor_id + 1} by Genotype', fontsize=14, fontweight='bold')
    ax_geno.set_ylabel('Mean Factor Score', fontsize=12)
    ax_geno.set_xlabel('Genotype', fontsize=12)
    ax_geno.grid(axis='y', alpha=0.3)
    
    # Add value labels on top of genotype bars
    for bar, mean_val in zip(bars_geno, summary_stats_geno['mean']):
        ax_geno.text(bar.get_x() + bar.get_width()/2, bar.get_height() + summary_stats_geno['sem'].max()*0.1,
                     f'{mean_val:.3f}', ha='center', va='bottom', fontsize=10)

plt.tight_layout()
plt.show()

# %% [markdown]
# Significant vs non-significant associations between factors and covariates are determined using a threshold automatically obtained using Otsu’s method.   
# This threshold can assist in defining the number of inferred factors (K).   
# For example, if a considerable proportion of factors fail to align with any covariates, it may prompt the to reduce K.   

# %%
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning, module="numpy.core.fromnumeric")
    matched_factor_dist, percent_matched_fact = efca.get_percent_matched_factors(fcat, threshold)
    matched_covariate_dist, percent_matched_cov = efca.get_percent_matched_covariates(fcat, threshold=threshold)

print('percent_matched_fact: ', percent_matched_fact)
print('percent_matched_cov: ', percent_matched_cov)

# %%
def plot_matched_factor_dist(matched_factor_dist, title='', save=False, save_path='./file.pdf', fontsize=18):
    """
    Plot the distribution of the number of matched covariate levels for each factor.

    Parameters
    ----------
    matched_factor_dist : array-like or Pandas Series
        Distribution of the number of matched covariate levels for each factor.
    title : str, optional
        Title of the plot. If not provided, a default title is used.
    save : bool, optional
        Whether to save the plot to a file.
    save_path : str, optional
        Path to save the plot if 'save' is True.
    """
    # Determine dynamic figure width based on the number of factors
    fig_width = np.round(len(matched_factor_dist) / 3)
    fig, ax = plt.subplots(figsize=(fig_width, 4))
    
    # Plot the bar chart
    x_positions = np.arange(len(matched_factor_dist))
    ax.bar(x_positions, matched_factor_dist, color='steelblue', edgecolor='black')
    
    # Set x-axis ticks and labels as F1, F2, ..., etc.
    labels = [f'F{i}' for i in range(1, len(matched_factor_dist) + 1)]
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, rotation=90, fontsize=fontsize)
    
    # Set labels and title similar to file_context_0 formatting
    ax.set_ylabel('Number of matched covariate levels', fontsize=fontsize, labelpad=12)
    ax.set_xlabel('Factors', fontsize=fontsize)
    default_title = 'Number of matched covariate levels per factor'
    ax.set_title(title if title else default_title, fontsize=fontsize, fontweight='bold')
    
    # Set y-axis ticks to show integer values
    ax.set_yticks(np.arange(0, max(matched_factor_dist) + 1, 1))
    ax.tick_params(axis='y', labelsize=11)
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save_path, bbox_inches='tight')
    plt.show()

plot_matched_factor_dist(matched_factor_dist, fontsize=12)

# %%
def plot_matched_covariate_dist(matched_covariate_dist, covariate_levels , title='',
                                save=False, save_path='./file.pdf', fontsize=12):
      """
        plot the distribution of the number of matched factors for each covariate level
        matched_covariate_dist: the distribution of the number of matched factors for each covariate level
        covariate_levels: the covariate levels
        title: the title of the plot
        save: whether to save the plot
        save_path: the path to save the plot

      """
      plt.figure(figsize=(np.round(len(matched_covariate_dist)/3),4))
      plt.bar(np.arange(len(matched_covariate_dist)), matched_covariate_dist)
      ### add covariate levels to the xticks
      plt.xticks(np.arange(len(matched_covariate_dist)), covariate_levels)

      ### make the xticks vertical and set the fontsize to 14
      plt.xticks(rotation=90, fontsize=fontsize)
      #plt.xlabel('Number of matched factors')
      ## set y ticks as digits and remove the decimal points and half points
      plt.yticks(np.arange(0, max(matched_covariate_dist)+1, 1), fontsize=fontsize) 
      plt.ylabel('Number of matched factors', fontsize=fontsize)
      plt.title(title)
      
      if save:
        plt.savefig(save_path, bbox_inches='tight')
      plt.show()

plot_matched_covariate_dist(
    matched_covariate_dist, 
    covariate_levels=all_covariate_levels,
    fontsize=12
)

# %% [markdown]
# We can check the correlation between the factors and the library size which was regressed out

# %%
factor_libsize_correlation = corr.get_factor_libsize_correlation(factor_scores, library_size = data.obs["total_counts"])
vis.plot_factor_cor_barplot(factor_libsize_correlation, 
             title='Correlation of factors with library size', 
             y_label='Correlation', x_label='Factors')

# %% [markdown]
# # **Step 3: Interpretability scores:** 
# The third step of sciRED involves quantifying the interpretability of identified factors.    
# Defined four categories of metrics: separability, effect size, specificity, and homogeneity which are presented as the FIST table. 

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
asv_genotype = met.average_scaled_var(factor_scores, data.obs["genotype"], mean_type='arithmetic')

#### plot the ralative variance table
svt_condition = met.scaled_var_table(factor_scores, data.obs["condition"])
svt_genotype = met.scaled_var_table(factor_scores, data.obs["genotype"])
svt = pd.concat([svt_condition, svt_genotype], axis=0)
vis.plot_relativeVar(svt.iloc[:,0:15], title='Relative variance score table')

# %%
########### create factor-interpretibility score table (FIST) ######
metrics_dict = {'Bimodality':bimodality_score, 
                    'Specificity':simpson_fcat,
                    'Effect size': factor_variance,
                    'Homogeneity (condition)':asv_condition,
                    'Homogeneity (genotype)':asv_genotype}

fist = met.FIST(metrics_dict)
### subset the first 15 factors of fist dataframe
vis.plot_FIST(fist.iloc[0:15,:])

# %% [markdown]
# # **Step 4: Biological Interpretation of a Selected Factor**

# %%
INTERESTING_FACTOR_ID = 0

# Get the loadings for the factor of interest
factor_loadings_for_factor_k = factor_loading[:, INTERESTING_FACTOR_ID]
gene_loadings = pd.Series(factor_loadings_for_factor_k, index=genes)
sorted_gene_loadings = gene_loadings.sort_values(ascending=False)

print(f"--- Top genes for Factor {INTERESTING_FACTOR_ID + 1} ---")

print("\n--- Top 20 Positive-Loading Genes (up-regulated in one condition) ---")
print(sorted_gene_loadings.head(20))

print("\n--- Top 20 Negative-Loading Genes (up-regulated in the other condition) ---")
print(sorted_gene_loadings.tail(20))

# %%
def plot_gene_umap(gene_of_interest):
    """
    Plots UMAP visualizations for a given gene:
      - Colors the UMAP by gene expression.
      - Colors the UMAP by condition.
      - Colors the UMAP by genotype.
    
    Parameters:
        gene_of_interest (str): Name of the gene to plot.
    """
    # Extract expression for the gene (ensure it's a numpy array, handling sparse matrices if necessary)
    expr = data[:, gene_of_interest].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray().flatten()
    else:
        expr = np.array(expr).flatten()
    
    umap_coords = data.obsm["X_umap"]
    
    # Create a figure with 3 subplots side by side (1 row, 3 columns)
    fig, axes = plt.subplots(1, 3, figsize=(24, 6))
    
    # Plot UMAP colored by gene expression on the first subplot
    sc = axes[0].scatter(umap_coords[:, 0], umap_coords[:, 1], c=expr, cmap='Reds', s=10)
    axes[0].set_xlabel("UMAP1")
    axes[0].set_ylabel("UMAP2")
    axes[0].set_title(f"UMAP Colored by {gene_of_interest} Expression")
    fig.colorbar(sc, ax=axes[0], label=f"{gene_of_interest} Expression")
    
    # Plot UMAP colored by condition on the second subplot
    conditions = data.obs["condition"]
    for cond in conditions.unique():
        idx = conditions == cond
        axes[1].scatter(umap_coords[idx, 0], umap_coords[idx, 1], s=10, label=cond)
    axes[1].set_xlabel("UMAP1")
    axes[1].set_ylabel("UMAP2")
    axes[1].set_title("UMAP Colored by Condition")
    axes[1].legend(title="Condition")
    
    # Plot UMAP colored by genotype on the third subplot
    genotypes = data.obs["genotype"]
    for genotype in genotypes.unique():
        idx = genotypes == genotype
        axes[2].scatter(umap_coords[idx, 0], umap_coords[idx, 1], s=10, label=genotype)
    axes[2].set_xlabel("UMAP1")
    axes[2].set_ylabel("UMAP2")
    axes[2].set_title("UMAP Colored by Genotype")
    axes[2].legend(title="Genotype")
    
    plt.tight_layout()
    plt.show()

# %%
INTERESTING_FACTOR_IDS = [0, 2, 10, 12]
combined_loading_tables = []  # List to store individual loading tables

for factor_id in INTERESTING_FACTOR_IDS:
    # Extract loadings for current factor from factor_loading array
    factor_loadings_for_factor = factor_loading[:, factor_id]
    
    # Create a pandas Series with gene names as the index
    gene_loadings = pd.Series(factor_loadings_for_factor, index=genes)
    
    # Sort the gene loadings in descending order
    sorted_gene_loadings = gene_loadings.sort_values(ascending=False)
    
    # Create a table (DataFrame) with two columns: one for gene names and one for their loading
    loading_table = pd.DataFrame({
        "gene": sorted_gene_loadings.index,
        "loading": sorted_gene_loadings.values
    })
    
    # Add a column to indicate the factor number (1-indexed for display)
    loading_table["factor"] = factor_id + 1
    
    print(f"--- Table for Factor {factor_id + 1} ---")
    print(loading_table.head(5))
    print("--------------------------------")
    print(loading_table.tail(5))
    print("\n")
    
    combined_loading_tables.append(loading_table)

# %%
# Example usage: Plotting UMAPs for different genes
for gene in ["Ubb", "Mir99ahg", "Dpp10"]:
    plot_gene_umap(gene)

# %% [markdown]
# # Modify sciRED functions to maintain directionality information

# %%
def get_AUC_all_factors_alevel_directional(factor_scores, a_binary_cov) -> np.ndarray:
    '''
    calculate the AUC of all the factors for a covariate level, preserving directionality.
    return a list of AUCs for all the factors
    factor_scores: a matrix of factor scores
    a_binary_cov: a binary vector of the covariate
    '''
    AUC_alevel_factors = []
    wilcoxon_pvalue_alevel_factors = []
    for i in range(factor_scores.shape[1]):
        a_factor = factor_scores[:,i]
        AUC, wilcoxon_pvalue = efca.get_AUC_alevel(a_factor, a_binary_cov)

        # Do NOT take absolute value or scale to 0-1 range
        AUC_alevel_factors.append(AUC)
        wilcoxon_pvalue_alevel_factors.append(wilcoxon_pvalue)
    ### convert to numpy array
    AUC_alevel_factors = np.asarray(AUC_alevel_factors)
    return AUC_alevel_factors


def get_importance_df_directional(factor_scores, a_binary_cov, time_eff=True) -> pd.DataFrame:
    '''
    calculate the importance of each factor for each covariate level, preserving directionality.
    factor_scores: numpy array of the factor scores for all the cells (n_cells, n_factors)
    a_binary_cov: numpy array of the binary covariate for a covariate level (n_cells, )
    time_eff: if True, skip RandomForest which is time consuming
    '''

    models = {'LogisticRegression': LogisticRegression(solver='lbfgs', max_iter=500),
              'DecisionTree': DecisionTreeClassifier(),
              'RandomForest': RandomForestClassifier(),
              'XGB': XGBClassifier(),
              'KNeighbors_permute': KNeighborsClassifier()}
    
    if time_eff:
        ### remove RandomForest, KN from the models dictionary
        models.pop('RandomForest')
        models.pop('KNeighbors_permute')
        
    importance_dict = {}

    for model_name, model in models.items():
        X, y = factor_scores, a_binary_cov
        model.fit(X, y)

        if model_name == 'LogisticRegression':
            # Use raw coefficients for directionality
            importance_dict[model_name] = model.coef_[0]

        elif model_name in ['DecisionTree', 'RandomForest', 'XGB']:
            # get importance values (these are inherently non-directional, but included for completeness)
            importance_dict[model_name] = model.feature_importances_

        elif model_name == 'KNeighbors_permute':
            # perform permutation importance (these are inherently non-directional, but included for completeness)
            perm_results = permutation_importance(model, X, y, scoring='accuracy')
            importance_dict[model_name] = perm_results.importances_mean # Keep raw mean, not absolute
    
    #### adding AUC as a measure of importance, preserving directionality
    AUC_alevel = get_AUC_all_factors_alevel_directional(factor_scores, a_binary_cov)
    importance_dict['AUC'] = AUC_alevel

    importance_df = pd.DataFrame.from_dict(importance_dict, orient='index',
                                           columns=['F'+str(i) for i in range(1, factor_scores.shape[1]+1)])
    return importance_df


def get_mean_importance_level_directional(importance_df_a_level, scale, mean) -> np.ndarray:
    '''
    calculate the mean importance of one level of a given covariate, preserving directionality.
    importance_df_a_level: a dataframe of the importance of each factor for a given covariate level
    scale: 'standard', 'minmax' or 'rank'
    mean: 'arithmatic' or 'geometric'
    '''
    importance_df_np = np.asarray(importance_df_a_level)
    ### normalize the importance score of each classifier in importance_df_np matrix
    if scale == 'standard':
        ### scale each row of the importance_df_np ( a model's importance results) to have zero mean and unit variance
        importance_df_np = (importance_df_np - importance_df_np.mean(axis=1, keepdims=True))/importance_df_np.std(axis=1, keepdims=True)

    if scale == 'minmax':
        ### scale each row of the importance_df_np to be between 0 and 1
        # For minmax scaling with directionality, we need to handle negative values correctly.
        # This scaling will map the range [min, max] to [0, 1].
        # If you need to preserve negative values, a different scaling approach might be needed,
        # or simply skip minmax scaling for directional results.
        importance_df_np = (importance_df_np - importance_df_np.min(axis=1, keepdims=True))/(importance_df_np.max(axis=1, keepdims=True) - importance_df_np.min(axis=1, keepdims=True))
    
    if scale == 'rank':
        ### replace each row of the importance_df_np with its rank
        # Rank data is inherently non-directional, but included for consistency.
        importance_df_np = np.apply_along_axis(ss.rankdata, 1, importance_df_np)
        ### for each row, devide ranks to n (number of factors) to get a value between 0 and 1
        importance_df_np = importance_df_np/importance_df_np.shape[1]

    ### calculate the mean of the importance_df_np matrix
    if mean == 'arithmatic':
        # Do NOT take absolute value here
        importance_df = np.mean(importance_df_np, axis=0)

    if mean == 'geometric':
        # Geometric mean requires positive values. If directional scores include negatives,
        # geometric mean is not appropriate. For now, we'll keep the original behavior
        # of replacing zeros/negatives with a small positive value, but this will
        # lose true negative directionality for geometric mean.
        # A warning or error might be appropriate here in a production system.
        #### if any value in a column is equal to zero, add a small value
        if np.any(importance_df_np == 0):
                importance_df_np[importance_df_np == 0] = 1e-10
        ### if any value is less than zero, replace with absolute value (this will lose directionality for geometric mean)
        if np.any(importance_df_np < 0):
            importance_df_np[importance_df_np < 0] = np.abs(importance_df_np[importance_df_np < 0])
        ### calculate the geometric mean of each column
        importance_df = ss.gmean(importance_df_np, axis=0)

    return importance_df


def FCAT_directional(covariate_vec, factor_scores,
                     scale='standard', mean='arithmatic', time_eff=True) -> pd.DataFrame:
    '''
    calculate the mean importance of all levels of a given covariate, preserving directionality.
    Returns a dataframe of size (num_levels, num_components).
    covariate_vec: numpy array of the covariate vector (n_cells, )
    factor_scores: numpy array of the factor scores for all the cells (n_cells, n_factors)
    '''

    mean_importance_df = pd.DataFrame(columns=['F'+str(i) for i in range(1, factor_scores.shape[1]+1)])

    for covariate_level in np.unique(covariate_vec):
        a_binary_cov = efca.get_binary_covariate(covariate_vec, covariate_level)
        importance_df_a_level = get_importance_df_directional(factor_scores, a_binary_cov, time_eff=time_eff)
        mean_importance_a_level = get_mean_importance_level_directional(importance_df_a_level, scale, mean)

        mean_importance_df.loc[covariate_level] = mean_importance_a_level

    return mean_importance_df

# %%
%%capture
####################################
#### FCAT score calculation ######
####################################

### FCAT needs to be calculated for each covariate separately
fcat_condition_directional = FCAT_directional(data.obs["condition"], factor_scores, scale='standard', mean='arithmatic')
fcat_genotype_directional = FCAT_directional(data.obs["genotype"], factor_scores, scale='standard', mean='arithmatic')

### concatenate FCAT table for protocol and cell line
fcat_directional = pd.concat([fcat_condition_directional, fcat_genotype_directional], axis=0)

# %%
### visualize the first 15 factors
vis.plot_FCAT(fcat_directional.iloc[:,0:15],title='Directional FCAT Scores', color='coolwarm',x_axis_fontsize=35, 
              y_axis_fontsize=35, title_fontsize=35,
              x_axis_tick_fontsize=32, y_axis_tick_fontsize=34)

# %% [markdown]
# # GSEA on the top factors

# %%
import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
from gseapy import prerank, enrichr
import os
import concurrent.futures
import warnings
warnings.filterwarnings('ignore')

# %%
INTERESTING_FACTOR_IDS = [0, 2, 10, 12]

# Define gene sets of interest
GENE_SETS_OF_INTEREST = [
    'GO_Biological_Process_2023',  # Note: 2025 might not exist yet, using 2023
    'GO_Cellular_Component_2023',
    'GO_Molecular_Function_2023',
    'KEGG_2019_Mouse',
    'Reactome_2022',  # Note: 2024 might not exist yet, using 2022
    'WikiPathways_2024_Mouse'  # Note: 2024 might not exist yet, using 2019
]

# %%
# Store all results
all_results = {
    'gsea_prerank': {},
    'enrichr': {},
    'summary': {}
}

# %%
def check_and_update_gene_sets():
    """Check available gene sets and update list with actually available ones."""
    print("Checking available gene sets for Mouse...")
    available_libs = gp.get_library_name(organism='Mouse')
    
    # Create mapping of requested to available
    gene_set_mapping = {}
    updated_sets = []
    
    for requested_set in GENE_SETS_OF_INTEREST:
        # Find exact match first
        if requested_set in available_libs:
            gene_set_mapping[requested_set] = requested_set
            updated_sets.append(requested_set)
        else:
            # Try to find closest match
            base_name = requested_set.rsplit('_', 1)[0]  # Remove year
            matches = [lib for lib in available_libs if base_name in lib]
            if matches:
                # Get the most recent version
                best_match = sorted(matches)[-1]
                gene_set_mapping[requested_set] = best_match
                updated_sets.append(best_match)
                print(f"  {requested_set} -> {best_match}")
            else:
                print(f"  {requested_set} -> NOT FOUND")
    
    return updated_sets, gene_set_mapping

available_sets, mapping = check_and_update_gene_sets()
print(available_sets)

# %%
def prepare_ranking_for_gsea(loading_table):
    """Prepare gene rankings from factor loadings for GSEA."""
    ranking = pd.Series(
        data=loading_table['loading'].values,
        index=loading_table['gene'].values
    )
    ranking = ranking.sort_values(ascending=False)
    ranking = ranking[~ranking.index.duplicated(keep='first')]
    return ranking

def run_gsea_prerank(loading_table, factor_name, gene_sets, output_dir='gsea_results'):
    """Run GSEA preranked analysis."""
    ranking = prepare_ranking_for_gsea(loading_table)
    
    results = {}
    for gene_set in gene_sets:
        print(f"  Running prerank with {gene_set}...")
        
        try:
            # Create output directory
            outdir = f'{output_dir}/{factor_name}/{gene_set}'
            os.makedirs(outdir, exist_ok=True)
            
            # Try to download gene set if needed
            gmt_file = f'{gene_set}.gmt'
            if not os.path.exists(gmt_file):
                try:
                    gp.parser.download_library(name=gene_set, organism='Mouse')
                except:
                    pass
            
            # Run prerank
            pre_res = gp.prerank(
                rnk=ranking,
                gene_sets=gene_set,
                threads=4,
                min_size=10,
                max_size=500,
                permutation_num=100,  # Reduced for speed, increase to 1000 for publication
                outdir=outdir,
                format='pdf',
                seed=42,
                verbose=False
            )
            
            if pre_res.res2d is not None and not pre_res.res2d.empty:
                results[gene_set] = pre_res.res2d
                print(f"    Success! Found {len(pre_res.res2d)} pathways")
                
        except Exception as e:
            print(f"    Error: {str(e)[:50]}...")
            results[gene_set] = None
    
    return results

# %%
available_sets[:3]

# %%
def process_factor(i, factor_id):
    if i >= len(combined_loading_tables):
        return None
    loading_table = combined_loading_tables[i]
    factor_name = f"Factor_{factor_id + 1}"
    
    prerank_results = run_gsea_prerank(loading_table, factor_name, available_sets[:3])  # Limit to GO terms for speed
    return factor_name, prerank_results

with concurrent.futures.ThreadPoolExecutor() as executor:
    futures = [executor.submit(process_factor, i, factor_id) for i, factor_id in enumerate(INTERESTING_FACTOR_IDS)]
    for future in concurrent.futures.as_completed(futures):
        result = future.result()
        if result is not None:
            factor_name, prerank_results = result
            all_results['gsea_prerank'][factor_name] = prerank_results

# %%
def run_enrichr_analysis(loading_table, factor_name, gene_sets, top_n=150):
    """Run Enrichr analysis."""
    # Get top and bottom genes
    top_genes = loading_table.nlargest(top_n, 'loading')['gene'].tolist()
    bottom_genes = loading_table.nsmallest(top_n, 'loading')['gene'].tolist()
    
    results = {'top': {}, 'bottom': {}}
    
    for direction, gene_list in [('top', top_genes), ('bottom', bottom_genes)]:
        print(f"  Analyzing {direction} {len(gene_list)} genes...")
        
        for gene_set in gene_sets:
            try:
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=gene_set,
                    organism='Mouse',
                    outdir=f'enrichr/{factor_name}/{direction}_{gene_set}',
                    cutoff=0.05
                )
                
                if enr.results is not None and not enr.results.empty:
                    results[direction][gene_set] = enr.results
                    print(f"    {gene_set}: {len(enr.results)} pathways")
                    
            except Exception as e:
                print(f"    {gene_set} error: {str(e)[:50]}...")
                results[direction][gene_set] = None
    
    return results

# %%
# %%capture
# def process_enrichr(idx, factor_id):
#     if idx >= len(combined_loading_tables):
#         return None
#     loading_table = combined_loading_tables[idx]
#     factor_name = f"Factor_{factor_id + 1}"
    
#     enrichr_results = run_enrichr_analysis(loading_table, factor_name, available_sets, top_n=25)
#     return factor_name, enrichr_results

# with concurrent.futures.ThreadPoolExecutor() as executor:
#     futures = [executor.submit(process_enrichr, i, factor_id) for i, factor_id in enumerate(INTERESTING_FACTOR_IDS)]
#     for future in concurrent.futures.as_completed(futures):
#         result = future.result()
#         if result is not None:
#             factor_name, enrichr_results = result
#             all_results['enrichr'][factor_name] = enrichr_results

# %%
# Sequential execution (simplified)
for i, factor_id in enumerate(INTERESTING_FACTOR_IDS):
    if i >= len(combined_loading_tables):
        continue
        
    loading_table = combined_loading_tables[i]
    factor_name = f"Factor_{factor_id + 1}"
    
    enrichr_results = run_enrichr_analysis(loading_table, factor_name, available_sets, top_n=25)
    
    if enrichr_results is not None:
        all_results['enrichr'][factor_name] = enrichr_results

# %%
def create_factor_summary(prerank_results, enrichr_results, factor_name):
    """Create a summary of enrichment results for a factor."""
    summary = {
        'factor': factor_name,
        'prerank_pathways': 0,
        'enrichr_top_pathways': 0,
        'enrichr_bottom_pathways': 0,
        'top_pathways': []
    }
    
    # Count prerank results
    for gene_set, results in prerank_results.items():
        if results is not None:
            sig_pathways = results[results['FDR q-val'] < 0.05]
            summary['prerank_pathways'] += len(sig_pathways)
            
            # Add top pathway
            if len(sig_pathways) > 0:
                top_pathway = sig_pathways.iloc[0]
                summary['top_pathways'].append({
                    'method': 'prerank',
                    'gene_set': gene_set,
                    'pathway': top_pathway.name,
                    'NES': top_pathway['NES'],
                    'FDR': top_pathway['FDR q-val']
                })
    
    # Count enrichr results
    for direction in ['top', 'bottom']:
        for gene_set, results in enrichr_results[direction].items():
            if results is not None:
                sig_pathways = results[results['Adjusted P-value'] < 0.05]
                if direction == 'top':
                    summary['enrichr_top_pathways'] += len(sig_pathways)
                else:
                    summary['enrichr_bottom_pathways'] += len(sig_pathways)
                
                # Add top pathway
                if len(sig_pathways) > 0:
                    top_pathway = sig_pathways.iloc[0]
                    summary['top_pathways'].append({
                        'method': f'enrichr_{direction}',
                        'gene_set': gene_set,
                        'pathway': top_pathway['Term'],
                        'pvalue': top_pathway['P-value'],
                        'adj_pvalue': top_pathway['Adjusted P-value']
                    })
    # Define a sorting key function
    def get_sort_key(pathway_dict):
        # Use FDR for prerank, Adjusted P-value for enrichr
        # Lower is better
        if 'FDR' in pathway_dict:
            return pathway_dict.get('FDR', 1.0)
        elif 'adj_pvalue' in pathway_dict:
            return pathway_dict.get('adj_pvalue', 1.0)
        return 1.0

    # Sort all collected "top" pathways to find the true top one
    summary['top_pathways'].sort(key=get_sort_key)

    return summary

# %%
for i, factor_id in enumerate(INTERESTING_FACTOR_IDS):
    factor_name = f"Factor_{factor_id + 1}"
    summary = create_factor_summary(all_results['gsea_prerank'][factor_name], all_results['enrichr'][factor_name], factor_name)
    all_results['summary'][factor_name] = summary

# %%
def plot_comprehensive_results(all_results, factor_name):
    """Create comprehensive visualization of results."""
    fig, ax = plt.subplots(1, 2, figsize=(15, 12))
    fig.suptitle(f'Comprehensive Enrichment Analysis - {factor_name}', fontsize=16)
    
    # Plot 1: Top Enrichr results (top genes)
    enrichr_top = all_results['enrichr'][factor_name]['top']
    plot_data = []
    
    for gene_set, results in enrichr_top.items():
        if results is not None and not results.empty:
            top_terms = results.head(5)
            for _, row in top_terms.iterrows():
                plot_data.append({
                    'Term': row['Term'][:40] + '...' if len(row['Term']) > 40 else row['Term'],
                    '-log10(p)': -np.log10(row['P-value']),
                    'Gene Set': gene_set.split('_')[0]
                })
    
    if plot_data:
        plot_df = pd.DataFrame(plot_data)
        plot_df = plot_df.sort_values('-log10(p)', ascending=True).tail(15)
        
        colors = {'GO': 'skyblue', 'KEGG': 'lightcoral', 'Reactome': 'lightgreen', 'Wiki': 'plum'}
        bar_colors = [colors.get(row['Gene Set'].split('_')[0], 'gray') for _, row in plot_df.iterrows()]
        
        ax[0].barh(range(len(plot_df)), plot_df['-log10(p)'], color=bar_colors)
        ax[0].set_yticks(range(len(plot_df)))
        ax[0].set_yticklabels(plot_df['Term'])
        ax[0].set_xlabel('-log10(p-value)')
        ax[0].set_title('Top Enriched Pathways (High Loading Genes)')
        ax[0].axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
    
    # Plot 2: Top Enrichr results (bottom genes)
    enrichr_bottom = all_results['enrichr'][factor_name]['bottom']
    plot_data = []
    
    for gene_set, results in enrichr_bottom.items():
        if results is not None and not results.empty:
            top_terms = results.head(5)
            for _, row in top_terms.iterrows():
                plot_data.append({
                    'Term': row['Term'][:40] + '...' if len(row['Term']) > 40 else row['Term'],
                    '-log10(p)': -np.log10(row['P-value']),
                    'Gene Set': gene_set.split('_')[0]
                })
    
    if plot_data:
        plot_df = pd.DataFrame(plot_data)
        plot_df = plot_df.sort_values('-log10(p)', ascending=True).tail(15)
        
        bar_colors = [colors.get(row['Gene Set'].split('_')[0], 'gray') for _, row in plot_df.iterrows()]
        
        ax[1].barh(range(len(plot_df)), plot_df['-log10(p)'], color=bar_colors)
        ax[1].set_yticks(range(len(plot_df)))
        ax[1].set_yticklabels(plot_df['Term'])
        ax[1].set_xlabel('-log10(p-value)')
        ax[1].set_title('Top Enriched Pathways (Low Loading Genes)')
        ax[1].axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(f'{factor_name}_Enrichr_summary.pdf', dpi=300, bbox_inches='tight')
    plt.show()

# %%
# Create visualizations
for factor_name in all_results['summary'].keys():
    plot_comprehensive_results(all_results, factor_name)

# %%
def save_comprehensive_results(all_results, available_sets):
    """Save all results to files."""
    print("\n" + "="*60)
    print("SAVING RESULTS")
    print("="*60)
    
    # Create results directory
    os.makedirs('enrichment_results_summary', exist_ok=True)
    
    # Save summary
    summary_df = pd.DataFrame([s for s in all_results['summary'].values()])
    summary_df.to_csv('enrichment_results_summary/summary_all_factors.csv', index=False)
    print("Saved summary_all_factors.csv")
    
    # Save detailed results for each factor
    for factor_name in all_results['enrichr'].keys():
        # Save Enrichr results
        for direction in ['top', 'bottom']:
            for gene_set, results in all_results['enrichr'][factor_name][direction].items():
                if results is not None:
                    filename = f'enrichment_results_summary/{factor_name}_enrichr_{direction}_{gene_set}.csv'
                    results.to_csv(filename, index=False)
        
        # Save GSEA results
        for gene_set, results in all_results['gsea_prerank'][factor_name].items():
            if results is not None:
                filename = f'enrichment_results_summary/{factor_name}_gsea_{gene_set}.csv'
                results.to_csv(filename)
    
    print("All results saved to 'enrichment_results_summary' directory")

# %%
save_comprehensive_results(all_results, available_sets)

# %%
# Print final summary
print("\n" + "="*60)
print("ANALYSIS COMPLETE - SUMMARY")
print("="*60)

for factor_name, summary in all_results['summary'].items():
    print(f"\n{factor_name}:")
    print(f"  GSEA significant pathways: {summary['prerank_pathways']}")
    print(f"  Enrichr top genes pathways: {summary['enrichr_top_pathways']}")
    print(f"  Enrichr bottom genes pathways: {summary['enrichr_bottom_pathways']}")
    
    i = 0
    if summary['top_pathways']:
        top_pathway_str = str(summary['top_pathways'][i].get('pathway', ''))
        while top_pathway_str=="0":
            i += 1
            top_pathway_str = str(summary['top_pathways'][i].get('pathway', ''))
        print(f"  Top pathway ({i}): {top_pathway_str[:60]}...")

# %% [markdown]
# # Format for STRING database

# %%
def extract_genes_for_string(combined_loading_tables, 
                           method='all', 
                           top_n=None, 
                           threshold=None,
                           factors=None,
                           unique_only=True):
    """
    Extract gene symbols from loading tables formatted for STRING database.
    
    Args:
        combined_loading_tables: List of DataFrames with gene loadings
        method: 'all', 'top_n', 'threshold', or 'bottom_n'
        top_n: Number of top genes to extract per factor (for 'top_n' method)
        threshold: Loading threshold (for 'threshold' method)
        factors: List of factor numbers to include (1-indexed), None for all
        unique_only: If True, return only unique genes across all factors
    
    Returns:
        String of gene symbols, one per line
    """
    genes = []
    
    # Filter factors if specified
    if factors:
        # Convert to 0-indexed
        factor_indices = [f - 1 for f in factors]
        tables_to_process = [combined_loading_tables[i] for i in factor_indices 
                           if i < len(combined_loading_tables)]
    else:
        tables_to_process = combined_loading_tables
    
    for table in tables_to_process:
        if method == 'all':
            # Get all genes
            genes.extend(table['gene'].tolist())
            
        elif method == 'top_n' and top_n:
            # Get top N genes by loading
            genes.extend(table.nlargest(top_n, 'loading')['gene'].tolist())
            
        elif method == 'bottom_n' and top_n:
            # Get bottom N genes by loading
            genes.extend(table.nsmallest(top_n, 'loading')['gene'].tolist())
            
        elif method == 'threshold' and threshold is not None:
            # Get genes above threshold
            genes.extend(table[table['loading'] > threshold]['gene'].tolist())
            
        elif method == 'abs_threshold' and threshold is not None:
            # Get genes with absolute loading above threshold
            genes.extend(table[abs(table['loading']) > threshold]['gene'].tolist())
    
    # Remove duplicates if requested
    if unique_only:
        genes = list(dict.fromkeys(genes))  # Preserves order
    
    # Return as string, one gene per line
    return '\n'.join(genes)

def get_top_genes_per_factor(combined_loading_tables, n=10):
    """Get top N genes from each factor."""
    print(f"Extracting top {n} genes from each factor...")
    for i, table in enumerate(combined_loading_tables):
        factor_num = table['factor'].iloc[0]
        print(f"\n--- Factor {factor_num} (top {n} genes) ---")
        top_genes = table.nlargest(n, 'loading')['gene'].tolist()
        print('\n'.join(top_genes))
    
    # Also return all unique genes
    all_genes = extract_genes_for_string(combined_loading_tables, 
                                       method='top_n', 
                                       top_n=n, 
                                       unique_only=True)
    print(f"\n--- All unique genes from top {n} of each factor ---")
    print(f"Total: {len(all_genes.split())} unique genes")
    print(all_genes)
    return all_genes

def get_high_loading_genes(combined_loading_tables, threshold=0.05):
    """Get genes with loading above threshold from all factors."""
    genes = extract_genes_for_string(combined_loading_tables, 
                                   method='threshold', 
                                   threshold=threshold,
                                   unique_only=True)
    print(f"Genes with loading > {threshold}:")
    print(f"Total: {len(genes.split())} unique genes")
    print(genes)
    return genes

def get_extreme_genes(combined_loading_tables, n=10):
    """Get top and bottom N genes from each factor."""
    all_genes = []
    
    for table in combined_loading_tables:
        factor_num = table['factor'].iloc[0]
        # Get top N
        top_genes = table.nlargest(n, 'loading')['gene'].tolist()
        # Get bottom N
        bottom_genes = table.nsmallest(n, 'loading')['gene'].tolist()
        all_genes.extend(top_genes + bottom_genes)
    
    # Remove duplicates
    unique_genes = list(dict.fromkeys(all_genes))
    genes_string = '\n'.join(unique_genes)
    
    print(f"Top {n} and bottom {n} genes from each factor:")
    print(f"Total: {len(unique_genes)} unique genes")
    print(genes_string)
    return genes_string
 

# %%
# Option 1: Get top 20 genes from each factor
print("OPTION 1: Top 20 genes per factor")
top_genes = get_top_genes_per_factor(combined_loading_tables, n=20)

# %%
# Option 2: Get all genes with loading > 0.05
print("OPTION 2: High loading genes")
high_loading = get_high_loading_genes(combined_loading_tables, threshold=0.05)

# %%
print("OPTION 3: Extreme genes")
extreme = get_extreme_genes(combined_loading_tables, n=15)

# %%
print("OPTION 4: Custom - Factors 1 & 2 only, top 30")
custom = extract_genes_for_string(combined_loading_tables, 
                                method='top_n',
                                top_n=30,
                                factors=[1, 2],
                                unique_only=True)
print(f"Total: {len(custom.split())} unique genes")
print(custom)

# %%
# Save to file for easy copying
with open('genes_for_string.txt', 'w') as f:
    f.write(top_genes)
print("\nGenes saved to 'genes_for_string.txt'")

# %%
# Quick function to get all genes from a specific factor
def get_genes_from_factor(combined_loading_tables, factor_num, n=None):
    """Get genes from a specific factor."""
    # Find the table for this factor
    for table in combined_loading_tables:
        if table['factor'].iloc[0] == factor_num:
            if n:
                genes = table.nlargest(n, 'loading')['gene'].tolist()
            else:
                genes = table['gene'].tolist()
            return '\n'.join(genes)
    return ""

# %%
# Quick copy-paste ready output
def quick_string_list(combined_loading_tables, n_per_factor=25):
    """Quick function to get a copy-paste ready list for STRING."""
    genes = extract_genes_for_string(combined_loading_tables, 
                                   method='top_n', 
                                   top_n=n_per_factor,
                                   unique_only=True)
    print("Copy and paste this into STRING:\n")
    print(genes)
    print(f"\n({len(genes.split())} genes total)")
    return genes

# %% [markdown]
# # Get Proteins

# %%
import pandas as pd
from typing import Dict, List
import mygene

# Using mygene package
def convert_genes_to_proteins_mygene(gene_symbols: List[str]) -> Dict[str, str]:
    """
    Convert mouse gene symbols to protein names using mygene API.
    
    Args:
        gene_symbols: List of mouse gene symbols
    
    Returns:
        Dictionary mapping gene symbols to protein names
    """
    mg = mygene.MyGeneInfo()
    
    # Query for mouse genes (taxid 10090 is for mouse)
    results = mg.querymany(gene_symbols, 
                            scopes='symbol', 
                            fields='name,symbol,uniprot', 
                            species='mouse',
                            returnall=True)
    
    gene_to_protein = {}
    for result in results['out']:
        if 'symbol' in result and 'name' in result:
            gene_to_protein[result['symbol']] = result['name']
    
    return gene_to_protein


# Main function to update your loading tables
def update_loading_tables_with_proteins(combined_loading_tables: List[pd.DataFrame]) -> List[pd.DataFrame]:
    """
    Update loading tables by adding protein names for each gene.
    Args:
        combined_loading_tables: List of DataFrames with gene loadings
    Returns:
        Updated list of DataFrames with protein names added
    """
    # Collect all unique gene symbols
    all_genes = set()
    for table in combined_loading_tables:
        all_genes.update(table['gene'].tolist())
    
    print(f"Converting {len(all_genes)} unique gene symbols to protein names...")
    
    # Get gene to protein mapping
    gene_to_protein = convert_genes_to_proteins_mygene(list(all_genes))

    print(f"Successfully mapped {len(gene_to_protein)} genes to proteins")
    
    # Update each loading table
    updated_tables = []
    for table in combined_loading_tables:
        # Create a copy to avoid modifying the original
        updated_table = table.copy()
        
        # Add protein name column
        updated_table['protein_name'] = updated_table['gene'].map(
            lambda x: gene_to_protein.get(x, f"Unknown protein for {x}")
        )
        
        # Reorder columns for better readability
        updated_table = updated_table[['factor', 'gene', 'protein_name', 'loading']]
        
        updated_tables.append(updated_table)
    
    return updated_tables

# %%
# # Alternative: If you just want to add protein names to existing tables in-place
# def add_protein_names_to_tables(combined_loading_tables: List[pd.DataFrame]) -> None:
#     """
#     Add protein names to existing loading tables in-place.
#     """
#     # Get all unique genes
#     all_genes = set()
#     for table in combined_loading_tables:
#         all_genes.update(table['gene'].tolist())
    
#     # Get mappings
#     gene_to_protein = convert_genes_to_proteins_mygene(list(all_genes))
    
#     # Update each table in-place
#     for table in combined_loading_tables:
#         table['protein_name'] = table['gene'].map(
#             lambda x: gene_to_protein.get(x, f"Unknown protein for {x}")
#         )

# %%
# Update tables with protein names
updated_loading_tables = update_loading_tables_with_proteins(combined_loading_tables)

# %%
# Display the updated tables
for i, (factor_id, table) in enumerate(zip(INTERESTING_FACTOR_IDS, updated_loading_tables)):
    print(f"\n--- Updated Table for Factor {factor_id + 1} ---")
    print("Top 5 genes/proteins:")
    print(table.head(5).to_string(index=False))
    print("\nBottom 5 genes/proteins:")
    print(table.tail(5).to_string(index=False))
    print("-" * 80)

# %%
# Save to files
for factor_id, table in zip(INTERESTING_FACTOR_IDS, updated_loading_tables):
    filename = f"loadings_with_proteins/factor_{factor_id+1}_loadings_with_proteins.csv"
    table.to_csv(filename, index=False)
    print(f"Saved to {filename}")


