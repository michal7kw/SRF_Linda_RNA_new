# %%
import scanpy as sc
import pandas as pd
import numpy as np
import os
import scipy.sparse as sp
import matplotlib.pyplot as plt

# %%
# --- Configuration ---
SAMPLES= ["Emx1_Ctrl", "Emx1_Mut"]
annotated_merged_adata_path = f"combine_data/results_from_raw_percentile_threshold/final_annotation/merged_raw_final_annotated_simple_mapmycells.h5ad"
output_adata_path = f"counts_trans/merged_with_trans_gene_normalized.h5ad"
transgene_name = "Rosa26_SBP1"

print(f"Annotated merged AnnData path: {annotated_merged_adata_path}")
print(f"Output AnnData path: {output_adata_path}")
print(f"Transgene name: {transgene_name}")

# %%
# --- Load AnnData objects (annotated_merged_adata_path) ---
try:
    adata_annotated_full = sc.read_h5ad(annotated_merged_adata_path)
    print(f"Loaded adata_annotated_full with shape: {adata_annotated_full.shape}")
except FileNotFoundError:
    print(f"ERROR: File not found at {annotated_merged_adata_path}")
    exit()
except Exception as e:
    print(f"ERROR loading {annotated_merged_adata_path}: {e}")
    exit()

# %%
# --- Prepare Barcodes in Annotated Data ---
print("\nPreparing barcodes in the annotated data...")
# Ensure they are strings
adata_annotated_full.obs_names = adata_annotated_full.obs_names.astype(str)
adata_annotated_full.obs['sample'] = adata_annotated_full.obs['sample'].astype(str)
adata_annotated_full.obs['default_barcodes'] = adata_annotated_full.obs_names

# Create 'original_barcode' column by removing potential suffixes (e.g., "-1")
try:
    adata_annotated_full.obs['original_barcode'] = adata_annotated_full.obs_names.str.replace(r'-\d+$', '', regex=True)
    print("Created 'original_barcode' column in adata_annotated_full.obs.")
    print(f"Head of original barcodes:\n{adata_annotated_full.obs['original_barcode'].head()}")
except Exception as e:
    print(f"ERROR creating 'original_barcode' column: {e}")
    exit()

# Create sample specific 'sample_specific_barcode' column by addoing sample name (e.g., "Emx1_Ctrl")
try:
    adata_annotated_full.obs['sample_specific_barcode'] = adata_annotated_full.obs['original_barcode'] + "-" + adata_annotated_full.obs['sample']
    print("Created 'sample_specific_barcode' column in adata_annotated_full.obs.")
    print(f"Head of sample_specific_barcode barcodes:\n{adata_annotated_full.obs['sample_specific_barcode'].head()}")
except Exception as e:
    print(f"ERROR creating 'sample_specific_barcode' column: {e}")
    exit()

# Set the 'sample_specific_barcode' as the new index (obs_names)
try:
    adata_annotated_full.obs_names = adata_annotated_full.obs['sample_specific_barcode']
    # Ensure index is unique, which is a requirement for AnnData
    if not adata_annotated_full.obs_names.is_unique:
        print(f"WARNING: Index is not unique after setting. Making unique...")
        adata_annotated_full.obs_names_make_unique()
    print("Set 'sample_specific_barcode' as index (obs_names) for adata_annotated_full.")
except Exception as e:
    print(f"ERROR setting 'sample_specific_barcode' as index for adata_annotated_full: {e}")
    exit()


# %%
# --- Load AnnData objects (transgene_adata_path) ---
adatas_transgene = dict() 
for sample in SAMPLES:
    transgene_adata_path = f"counts_trans/{sample}_normalized.h5ad"
    print(f"Transgene AnnData path: {transgene_adata_path}")
    print("\nLoading AnnData objects...")
    try:
        adatas_transgene[sample] = sc.read_h5ad(transgene_adata_path)
        print(f"Loaded adata_transgene with shape: {adatas_transgene[sample].shape}")
    except FileNotFoundError:
        print(f"ERROR: File not found at {transgene_adata_path}")
        exit()
    except Exception as e:
        print(f"ERROR loading {transgene_adata_path}: {e}")
        exit()

    # --- Prepare Barcodes in Annotated Data ---
    print("\nPreparing barcodes in the annotated data...")
    # Ensure they are strings
    adatas_transgene[sample].obs_names = adatas_transgene[sample].obs_names.astype(str)
    adatas_transgene[sample].obs["sample"] = adatas_transgene[sample].obs["sample"].astype(str)
    adatas_transgene[sample].obs["default_barcodes"] = adatas_transgene[sample].obs_names

    # Create sample specific 'sample_specific_barcode' column by addoing sample name (e.g., "Emx1_Ctrl")
    try:
        adatas_transgene[sample].obs['sample_specific_barcode'] = adatas_transgene[sample].obs_names + "-" + adatas_transgene[sample].obs["sample"]
        print(f"Created 'sample_specific_barcode' column in adatas_transgene[{sample}].obs.")
        print(f"Head of sample_specific_barcode barcodes:\n{adatas_transgene[sample].obs['sample_specific_barcode'].head()}")
    except Exception as e:
        print(f"ERROR creating 'sample_specific_barcode' column: {e}")
        exit()

    # Set the 'sample_specific_barcode' as the new index (obs_names)
    try:
        adatas_transgene[sample].obs_names = adatas_transgene[sample].obs['sample_specific_barcode'].astype(str)
        # Ensure index is unique, which is a requirement for AnnData
        if not adatas_transgene[sample].obs_names.is_unique:
            print(f"WARNING: Index for sample {sample} is not unique after setting. Making unique...")
            adatas_transgene[sample].obs_names_make_unique()
        print(f"Set 'sample_specific_barcode' as index (obs_names) for adatas_transgene[{sample}].")
    except Exception as e:
        print(f"ERROR setting 'sample_specific_barcode' as index for adatas_transgene[{sample}]: {e}")
        exit()

# %%
# --- Extract Transgene Expression ---
individual_expression_series = []  # To store Series for each sample
sample_order_for_concat = []     # To store the order of samples for pd.concat keys

print(f"\nExtracting expression for transgene: {transgene_name}...")
for sample, adata_transgene in adatas_transgene.items():
    # Ensure transgene exists in var_names
    if transgene_name not in adata_transgene.var_names:
        print(f"ERROR: Transgene '{transgene_name}' not found in {sample} adata_transgene.var_names.")
        print(f"Available var names head: {adata_transgene.var_names[:10].tolist()}")
        exit()

    # Get the expression vector for the transgene
    try:
        gene_idx = adata_transgene.var_names.get_loc(transgene_name)
        # Extract the column - result might be sparse or dense
        expression_vector = adata_transgene.X[:, gene_idx]

        # Convert to dense numpy array if it's sparse
        if sp.issparse(expression_vector):
            dense_expression = expression_vector.toarray().flatten()
        else:
            dense_expression = np.asarray(expression_vector).flatten()
        
        print(f"Successfully extracted expression for {transgene_name}. Shape: {dense_expression.shape}")
        print(f"--- Debug: Raw values for {transgene_name} from adata_transgene.X ---")
        print(f"Min: {np.min(dense_expression)}, Max: {np.max(dense_expression)}, Mean: {np.mean(dense_expression)}")
        # Check if all are integers
        is_integer = np.all(dense_expression == np.floor(dense_expression))
        print(f"Are all extracted values integers? {is_integer}")
        non_zero_values = dense_expression[dense_expression > 0]
        if len(non_zero_values) > 0:
            print(f"Number of non-zero values: {len(non_zero_values)}")
            print(f"Example non-zero values: {non_zero_values[:10]}")
        else:
            print("All extracted values are zero.")
        print(f"--------------------------------------------------------------------")

    except Exception as e:
        print(f"ERROR extracting expression for '{transgene_name}': {e}")
        exit()

    # Create a mapping from barcode to expression
    if len(adata_transgene.obs_names) != len(dense_expression):
        print(f"ERROR: Mismatch between number of barcodes ({len(adata_transgene.obs_names)}) and expression values ({len(dense_expression)}).")
        exit()

    current_series = pd.Series(dense_expression, index=adata_transgene.obs_names)
    individual_expression_series.append(current_series)
    sample_order_for_concat.append(sample)
    print(f"Processed sample {sample}. Series head:\n{current_series.head()}")

# %%
# Concatenate all individual series into one big expression map
concatenated_expression_map = pd.concat(individual_expression_series)
print(f"\nSuccessfully concatenated all expression series.")
print(f"Shape of concatenated_expression_map: {concatenated_expression_map.shape}")
print(f"Head of concatenated_expression_map:\n{concatenated_expression_map.head()}")
print(f"Tail of concatenated_expression_map:\n{concatenated_expression_map.tail()}")


# %%
# --- Align Transgene Expression Vector ---
print("\nAligning transgene expression vector to annotated data...")
# Create a new vector aligned with adata_annotated_full.obs_names
# Use the 'sample_specific_barcode' to map values from the expression_map Series
# Fill missing values with 0

aligned_expression = adata_annotated_full.obs_names.map(concatenated_expression_map).fillna(0).values
print(f"Created aligned expression vector with shape: {aligned_expression.shape}")

# Reshape to be a column vector (n_obs x 1)
aligned_expression_col = aligned_expression.reshape(-1, 1)

# Convert to sparse matrix format (CSC is efficient for column operations)
sparse_expression_col = sp.csc_matrix(aligned_expression_col)
print(f"Converted aligned expression to sparse matrix with shape: {sparse_expression_col.shape}")

# %%
# --- Add Transgene as a New Feature ---
print(f"\nAdding '{transgene_name}' as a new feature...")

# Create a minimal AnnData object for the transgene feature
print("Creating temporary AnnData object for the transgene feature...")
transgene_var = pd.DataFrame(index=[transgene_name])
# Ensure obs names match the main object for concatenation compatibility
# Use adata_annotated_full.obs_names as index for the temporary obs
temp_obs = pd.DataFrame(index=adata_annotated_full.obs_names)
adata_transgene_feature = sc.AnnData(X=sparse_expression_col, var=transgene_var, obs=temp_obs)
print(f"Temporary transgene AnnData shape: {adata_transgene_feature.shape}")

# %%
# --- Modify AnnData object to include transgene ---
if adata_annotated_full.raw is not None:
    print("Found existing .raw attribute. Modifying .raw to include transgene.")
    current_raw_adata = adata_annotated_full.raw

    if transgene_name in current_raw_adata.var_names:
        print(f"Feature '{transgene_name}' already exists in .raw.var_names. Overwriting its data in .raw.X.")
        existing_idx = current_raw_adata.var_names.get_loc(transgene_name)
        
        if not sp.isspmatrix_csc(current_raw_adata.X):
            raw_X_modifiable = current_raw_adata.X.tocsc()
        else:
            raw_X_modifiable = current_raw_adata.X
            
        raw_X_modifiable[:, existing_idx] = sparse_expression_col
        current_raw_adata.X = raw_X_modifiable.tocsr()
        # Ensure var_names are unique, although overwriting shouldn't change them
        current_raw_adata.var_names_make_unique()
        adata_annotated_full.raw = current_raw_adata # Assign back the modified Raw object
        print(f"Updated existing feature '{transgene_name}' in .raw.X. .raw shape: {adata_annotated_full.raw.shape}")
    else:
        print(f"Adding new feature '{transgene_name}' to .raw.")
        # Create a temporary AnnData from existing .raw components for concatenation
        temp_raw_for_concat = sc.AnnData(
            X=current_raw_adata.X.copy(),
            var=current_raw_adata.var.copy(),
            obs=adata_annotated_full.obs.copy() # Use main .obs for consistent indexing
        )
        temp_raw_for_concat.obs_names = adata_annotated_full.obs_names
        adata_transgene_feature.obs_names = adata_annotated_full.obs_names

        try:
            concatenated_raw_adata = sc.concat(
                [temp_raw_for_concat, adata_transgene_feature],
                axis=1, join='outer', merge='unique', uns_merge='unique', label="source"
            )
            concatenated_raw_adata.var_names_make_unique() # Ensure unique var names after concat
            # Replace the .raw attribute with this new, fully formed AnnData object
            adata_annotated_full.raw = concatenated_raw_adata
            print(f"Replaced .raw with concatenated data. New .raw.shape: {adata_annotated_full.raw.shape}")
        except Exception as e:
            print(f"ERROR during anndata.concat for .raw: {e}")
            exit()
else: # .raw does not exist
    print("No .raw attribute found. Adding transgene to main .X and creating .raw from it.")
    
    if transgene_name in adata_annotated_full.var_names:
        print(f"Feature '{transgene_name}' already exists in main AnnData. Overwriting its data in .X.")
        existing_idx = adata_annotated_full.var_names.get_loc(transgene_name)
        if not sp.isspmatrix_csc(adata_annotated_full.X):
            main_X_modifiable = adata_annotated_full.X.tocsc()
        else:
            main_X_modifiable = adata_annotated_full.X
        main_X_modifiable[:, existing_idx] = sparse_expression_col
        adata_annotated_full.X = main_X_modifiable.tocsr()
        adata_annotated_full.var_names_make_unique() # Ensure unique var names
    else:
        print(f"Adding new feature '{transgene_name}' to main AnnData .X and .var.")
        try:
            adata_transgene_feature.obs_names = adata_annotated_full.obs_names
            adata_annotated_full = sc.concat(
                [adata_annotated_full, adata_transgene_feature],
                axis=1, join='outer', merge='unique', uns_merge='unique', label="source"
            )
            adata_annotated_full.var_names_make_unique() # Ensure unique var names after concat
            print(f"Concatenated new feature to main AnnData. New shape: {adata_annotated_full.shape}")
        except Exception as e:
            print(f"ERROR during anndata.concat for main AnnData: {e}")
            exit()
            
    # Create .raw from the modified main AnnData
    print("Creating .raw attribute from the modified main AnnData...")
    adata_annotated_full.raw = adata_annotated_full.copy()
    adata_annotated_full.raw.var_names_make_unique() # Also ensure .raw has unique var names
    print(f".raw created. .raw.shape: {adata_annotated_full.raw.shape}")


# %%
# --- Verification ---
print("\nVerification:")
target_adata_updated = adata_annotated_full.raw if adata_annotated_full.raw is not None else adata_annotated_full

if transgene_name in target_adata_updated.var_names:
    print(f"  Feature '{transgene_name}' successfully added/updated in target var_names.")
    print(f"  Target matrix shape: {target_adata_updated.shape}")
    # Verify expression values for a few cells
    gene_idx_verify = target_adata_updated.var_names.get_loc(transgene_name)
    print(f"  Example expression values for '{transgene_name}': {target_adata_updated.X[:5, gene_idx_verify].toarray().flatten()}")
else:
    print(f"ERROR: Feature '{transgene_name}' NOT found in target var_names after update.")

# Remove the temporary original_barcode column
if 'original_barcode' in adata_annotated_full.obs.columns:
    del adata_annotated_full.obs['original_barcode']
    print("Removed temporary 'original_barcode' column from .obs.")

# %%
# Restore original barcodes
print("\nRestoring original barcodes...")
# Set the 'default_barcodes' as the new index (obs_names)
try:
    adata_annotated_full.obs_names = adata_annotated_full.obs['default_barcodes']
    # Ensure index is unique, which is a requirement for AnnData
    if not adata_annotated_full.obs_names.is_unique:
        print(f"WARNING: Index is not unique after setting. Making unique...")
        adata_annotated_full.obs_names_make_unique()
    print("Set 'sample_specific_barcode' as index (obs_names) for adata_annotated_full.")
except Exception as e:
    print(f"ERROR setting 'sample_specific_barcode' as index for adata_annotated_full: {e}")
    exit()

# %%
# --- Save Result ---
print("\nSaving the result...")
# Ensure the output directory exists
output_dir = os.path.dirname(output_adata_path)
if output_dir and not os.path.exists(output_dir):
    try:
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    except Exception as e:
        print(f"ERROR creating output directory '{output_dir}': {e}")
        exit()

# Save the modified AnnData object
try:
    # Remove the temporary original_barcode column if not needed
    adata_annotated_full.write_h5ad(output_adata_path, compression="gzip")
    print(f"Successfully saved merged AnnData with transgene feature to: {output_adata_path}")
except Exception as e:
    print(f"ERROR saving AnnData object to '{output_adata_path}': {e}")

print("\nScript finished.")
# %%

# --- Visualize Transgene Expression on UMAP (Optional) ---
# Visualization now needs to fetch data from the matrix, not .obs
print("\nGenerating UMAP plot colored by transgene expression (fetching from matrix)...")
target_adata_final = adata_annotated_full.raw if adata_annotated_full.raw is not None else adata_annotated_full

if 'X_umap_alt' in adata_annotated_full.obsm and transgene_name in target_adata_final.var_names:
    try:
        # Add the expression data temporarily to obs for plotting
        adata_annotated_full.obs[f'{transgene_name}_viz'] = target_adata_final[:, transgene_name].X.toarray().flatten()
        
        plot_filename = f"_{transgene_name}_umap_expression.png"
        sc.pl.embedding(
            adata_annotated_full,
            color=f'{transgene_name}_viz',
            basis='X_umap_alt',
            show=True,
            cmap='Reds',
            title=f'{transgene_name} Expression',
            size=20.0,  # Set point size
            save=plot_filename
        )
        print(f"Saved UMAP plot for {transgene_name} to {plot_filename}")
        
        # Clean up temporary column
        del adata_annotated_full.obs[f'{transgene_name}_viz']

    except Exception as e:
        print(f"ERROR generating UMAP plot: {e}")
        # Clean up temporary column in case of error
        if f'{transgene_name}_viz' in adata_annotated_full.obs:
            del adata_annotated_full.obs[f'{transgene_name}_viz']
elif 'X_umap_alt' not in adata_annotated_full.obsm:
    print("WARNING: 'X_umap_alt' not found in adata_annotated_full.obsm. Skipping UMAP plot.")
elif transgene_name not in target_adata_final.var_names:
    print(f"WARNING: Transgene '{transgene_name}' not found in target matrix var_names. Skipping UMAP plot.")

print("\nScript finished completely.")
# %%
