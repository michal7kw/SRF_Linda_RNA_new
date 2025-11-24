
# # %% [markdown]
# # ## Add Alternative UMAP Embedding

# # %%
# alt_leiden_key = "leiden_0.8"
# print("\nLoading alternative UMAP embedding...")
# alt_adata_path = os.path.join(PROJECT_DIR, "combine_data", "Archive", "from_preprocessed_counts", "results", "all_samples_merged.h5ad")
# alt_umap_key = 'X_umap_alt' # Key to store the alternative UMAP

# # Load obs_names, UMAP coordinates, and leiden_0.8 clustering from the alternative file
# print(f"Loading obs_names, X_umap, obs['leiden_0.8'], DG_majority_voting, ISO_majority_voting, DG_conf_score, and ISO_conf_score from {alt_adata_path}")
# adata_alt_subset = ad.read_h5ad(alt_adata_path, backed='r') # Read in backed mode initially

# # %%
# adata_alt_subset.obs.columns
# adata_alt_subset.obs.ISO_majority_voting

# # %%

# # --- Explicitly load required .obs columns into memory ---
# print("Loading required .obs columns from alternative data into memory...")
# obs_cols_to_load = []
# if 'leiden_0.8' in adata_alt_subset.obs.columns:
#     obs_cols_to_load.append('leiden_0.8')
# if 'DG_majority_voting' in adata_alt_subset.obs.columns: 
#     obs_cols_to_load.append('DG_majority_voting')
# if 'ISO_majority_voting' in adata_alt_subset.obs.columns: 
#      obs_cols_to_load.append('ISO_majority_voting')
# if 'DG_conf_score' in adata_alt_subset.obs.columns: 
#      obs_cols_to_load.append('DG_conf_score')
# if 'ISO_conf_score' in adata_alt_subset.obs.columns: 
#     obs_cols_to_load.append('ISO_conf_score')

# if obs_cols_to_load:
#     adata_alt_obs_in_memory = adata_alt_subset.obs[obs_cols_to_load].copy()
#     # We also need the original obs_names
#     alt_obs_names_orig = adata_alt_subset.obs_names.to_list()
#     print(f"Loaded {len(obs_cols_to_load)} .obs columns and obs_names into memory.")
# else:
#     print("Warning: No required .obs columns found to load into memory.")
#     adata_alt_obs_in_memory = pd.DataFrame(index=adata_alt_subset.obs_names) # Empty df with index
#     alt_obs_names_orig = adata_alt_subset.obs_names.to_list()
# # Also load UMAP coordinates if present
# alt_umap_coords = None
# if 'X_umap' in adata_alt_subset.obsm_keys():
#     umap_data = adata_alt_subset.obsm['X_umap']
#     if hasattr(umap_data, 'toarray'):  # Check if it's a sparse matrix
#         alt_umap_coords = umap_data.toarray().copy()
#     else:  # Assume it's already a dense numpy array
#         alt_umap_coords = umap_data.copy()
#     print("Loaded X_umap coordinates into memory.")
# # Close the backed file connection if possible (optional, may help release resources)
# # adata_alt_subset.file.close()
# # ----------------------------------------------------------

# # Check if required data exists in the alternative file (using in-memory data now)
# alt_umap_present = alt_umap_coords is not None
# alt_leiden_key = 'leiden_0.8' # Keep original key name
# alt_leiden_present = alt_leiden_key in adata_alt_obs_in_memory.columns
# alt_dg_present = 'DG_majority_voting' in adata_alt_obs_in_memory.columns
# alt_iso_present = 'ISO_majority_voting' in adata_alt_obs_in_memory.columns
# alt_dg_conf_present = 'DG_conf_score' in adata_alt_obs_in_memory.columns
# alt_iso_conf_present = 'ISO_conf_score' in adata_alt_obs_in_memory.columns

# if not alt_umap_present:
#         print(f"Warning: 'X_umap' not found in {alt_adata_path}. Skipping addition of alternative UMAP.")
# if not alt_leiden_present:
#         print(f"Warning: '{alt_leiden_key}' not found in {alt_adata_path}. Skipping addition of alternative Leiden clustering.")

# if not alt_dg_present:
#     print(f"Warning: 'DG_majority_voting' not found in {alt_adata_path}. Skipping addition of alternative DG majority voting.")

# if not alt_iso_present:
#     print(f"Warning: 'ISO_majority_voting' not found in {alt_adata_path}. Skipping addition of alternative ISO majority voting.")
    
# if not alt_dg_conf_present:
#     print(f"Warning: 'DG_conf_score' not found in {alt_adata_path}. Skipping addition of alternative DG confidence score.")
    
# if not alt_iso_conf_present:
#     print(f"Warning: 'ISO_conf_score' not found in {alt_adata_path}. Skipping addition of alternative ISO confidence score.")

# if alt_umap_present or alt_leiden_present or alt_dg_present or alt_iso_present or alt_dg_conf_present or alt_iso_conf_present:
#     # --- Barcode Matching Adjustment ---
#     # 1. Get barcodes and required data from alternative data (using in-memory data)
#     # alt_obs_names_orig is already loaded
#     print(f"Processing {len(alt_obs_names_orig)} cells from alternative data.")
#     print(f"  Example alternative barcode (original): {alt_obs_names_orig[0] if alt_obs_names_orig else 'N/A'}")

#     # 2. Strip sample suffix (e.g., '-Emx1_Ctrl') from alternative barcodes
#     # Regex: Remove hyphen followed by alphanumeric/underscore characters at the end
#     alt_obs_names_stripped = pd.Series(alt_obs_names_orig).str.replace(r'-[A-Za-z0-9_]+$', '', regex=True).tolist()
#     print(f"  Example alternative barcode (stripped): {alt_obs_names_stripped[0] if alt_obs_names_stripped else 'N/A'}")

#     # Create mapping for UMAP (if present, use in-memory coords)
#     alt_umap_df_stripped = None
#     if alt_umap_present and alt_umap_coords is not None:
#         # Use alt_umap_coords loaded earlier
#         alt_umap_df_stripped = pd.DataFrame(alt_umap_coords, index=alt_obs_names_stripped, columns=['alt_UMAP1', 'alt_UMAP2'])
#         alt_umap_df_stripped = alt_umap_df_stripped[~alt_umap_df_stripped.index.duplicated(keep='first')]
#         print(f"  Created mapping for alternative UMAP.")

#     # Create mapping for Leiden clustering (if present, use in-memory obs)
#     alt_leiden_series_stripped = None
#     new_alt_leiden_key = f"{alt_leiden_key}_alt"
#     if alt_leiden_present:
#         alt_leiden_data = adata_alt_obs_in_memory[alt_leiden_key]
#         alt_leiden_series_stripped = pd.Series(alt_leiden_data.values, index=alt_obs_names_stripped)
#         alt_leiden_series_stripped = alt_leiden_series_stripped[~alt_leiden_series_stripped.index.duplicated(keep='first')]
#         print(f"  Created mapping for alternative Leiden clustering ('{alt_leiden_key}'). Will be stored as '{new_alt_leiden_key}'.")

#     # Create mapping for DG majority voting (if present, use in-memory obs)
#     alt_dg_series_stripped = None
#     new_alt_dg_key = "DG_majority_voting_alt"
#     if alt_dg_present:
#         alt_dg_data = adata_alt_obs_in_memory['DG_majority_voting']
#         alt_dg_series_stripped = pd.Series(alt_dg_data.values, index=alt_obs_names_stripped)
#         alt_dg_series_stripped = alt_dg_series_stripped[~alt_dg_series_stripped.index.duplicated(keep='first')]
#         print(f"  Created mapping for alternative DG majority voting. Will be stored as '{new_alt_dg_key}'.")

#     # Create mapping for ISO majority voting (if present, use in-memory obs)
#     alt_iso_series_stripped = None
#     new_alt_iso_key = "ISO_majority_voting_alt"
#     if alt_iso_present:
#         alt_iso_data = adata_alt_obs_in_memory['ISO_majority_voting']
#         alt_iso_series_stripped = pd.Series(alt_iso_data.values, index=alt_obs_names_stripped)
#         alt_iso_series_stripped = alt_iso_series_stripped[~alt_iso_series_stripped.index.duplicated(keep='first')]
#         print(f"  Created mapping for alternative ISO majority voting. Will be stored as '{new_alt_iso_key}'.")

#     # Create mapping for DG confidence score (if present, use in-memory obs)
#     alt_dg_conf_series_stripped = None
#     new_alt_dg_conf_key = "DG_conf_score_alt"
#     if alt_dg_conf_present:
#         alt_dg_conf_data = adata_alt_obs_in_memory['DG_conf_score']
#         alt_dg_conf_series_stripped = pd.Series(alt_dg_conf_data.values, index=alt_obs_names_stripped)
#         alt_dg_conf_series_stripped = alt_dg_conf_series_stripped[~alt_dg_conf_series_stripped.index.duplicated(keep='first')]
#         print(f"  Created mapping for alternative DG confidence score. Will be stored as '{new_alt_dg_conf_key}'.")

#     # Create mapping for ISO confidence score (if present, use in-memory obs)
#     alt_iso_conf_series_stripped = None
#     new_alt_iso_conf_key = "ISO_conf_score_alt"
#     if alt_iso_conf_present:
#         alt_iso_conf_data = adata_alt_obs_in_memory['ISO_conf_score']
#         alt_iso_conf_series_stripped = pd.Series(alt_iso_conf_data.values, index=alt_obs_names_stripped)
#         alt_iso_conf_series_stripped = alt_iso_conf_series_stripped[~alt_iso_conf_series_stripped.index.duplicated(keep='first')]
#         print(f"  Created mapping for alternative ISO confidence score. Will be stored as '{new_alt_iso_conf_key}'.")

#     # 3. Get barcodes from the current adata_merged
#     current_obs_names_with_suffix = adata_merged.obs_names
#     print(f"  Example current barcode (with suffix): {current_obs_names_with_suffix[0] if not current_obs_names_with_suffix.empty else 'N/A'}")

#     # 4. Strip numeric suffix (e.g., '-0') from current barcodes
#     # Regex: Remove hyphen followed by digits at the end
#     current_obs_names_stripped = current_obs_names_with_suffix.str.replace(r'-\d+$', '', regex=True)
#     print(f"  Example current barcode (stripped): {current_obs_names_stripped[0] if not current_obs_names_stripped.empty else 'N/A'}")

#     # 5. Align UMAP and/or Leiden based on the stripped barcodes
#     # Align all potential data sources first based on stripped barcodes,
#     # then set index to the original adata_merged obs names (with suffix)
#     aligned_alt_umap_df = None
#     if alt_umap_df_stripped is not None:
#         aligned_alt_umap_df = alt_umap_df_stripped.reindex(current_obs_names_stripped)
#         aligned_alt_umap_df.index = current_obs_names_with_suffix

#     aligned_alt_leiden_series = None
#     if alt_leiden_series_stripped is not None:
#         aligned_alt_leiden_series = alt_leiden_series_stripped.reindex(current_obs_names_stripped)
#         aligned_alt_leiden_series.index = current_obs_names_with_suffix

#     aligned_alt_dg_series = None
#     if alt_dg_series_stripped is not None:
#         aligned_alt_dg_series = alt_dg_series_stripped.reindex(current_obs_names_stripped)
#         aligned_alt_dg_series.index = current_obs_names_with_suffix

#     aligned_alt_iso_series = None
#     if alt_iso_series_stripped is not None:
#         aligned_alt_iso_series = alt_iso_series_stripped.reindex(current_obs_names_stripped)
#         aligned_alt_iso_series.index = current_obs_names_with_suffix

#     aligned_alt_dg_conf_series = None
#     if alt_dg_conf_series_stripped is not None:
#         aligned_alt_dg_conf_series = alt_dg_conf_series_stripped.reindex(current_obs_names_stripped)
#         aligned_alt_dg_conf_series.index = current_obs_names_with_suffix

#     aligned_alt_iso_conf_series = None
#     if alt_iso_conf_series_stripped is not None:
#         aligned_alt_iso_conf_series = alt_iso_conf_series_stripped.reindex(current_obs_names_stripped)
#         aligned_alt_iso_conf_series.index = current_obs_names_with_suffix


#     # Identify cells with missing data in ANY aligned source
#     cells_to_remove = pd.Index([])
#     print("Checking for missing values in aligned alternative data...")

#     if aligned_alt_umap_df is not None:
#         missing_cells = aligned_alt_umap_df[aligned_alt_umap_df['alt_UMAP1'].isna()].index
#         if not missing_cells.empty:
#             print(f"  Identified {len(missing_cells)} cells missing in alternative UMAP data.")
#             cells_to_remove = cells_to_remove.union(missing_cells)

#     if aligned_alt_leiden_series is not None:
#         missing_cells = aligned_alt_leiden_series[aligned_alt_leiden_series.isna()].index
#         if not missing_cells.empty:
#                 print(f"  Identified {len(missing_cells)} cells missing in alternative Leiden data ('{alt_leiden_key}').")
#                 cells_to_remove = cells_to_remove.union(missing_cells)

#     if aligned_alt_dg_series is not None:
#         missing_cells = aligned_alt_dg_series[aligned_alt_dg_series.isna()].index
#         if not missing_cells.empty:
#             print(f"  Identified {len(missing_cells)} cells missing in alternative DG majority voting data.")
#             cells_to_remove = cells_to_remove.union(missing_cells)

#     if aligned_alt_iso_series is not None:
#         missing_cells = aligned_alt_iso_series[aligned_alt_iso_series.isna()].index
#         if not missing_cells.empty:
#             print(f"  Identified {len(missing_cells)} cells missing in alternative ISO majority voting data.")
#             cells_to_remove = cells_to_remove.union(missing_cells)

#     if aligned_alt_dg_conf_series is not None:
#         missing_cells = aligned_alt_dg_conf_series[aligned_alt_dg_conf_series.isna()].index
#         if not missing_cells.empty:
#             print(f"  Identified {len(missing_cells)} cells missing in alternative DG confidence score data.")
#             cells_to_remove = cells_to_remove.union(missing_cells)

#     if aligned_alt_iso_conf_series is not None:
#         missing_cells = aligned_alt_iso_conf_series[aligned_alt_iso_conf_series.isna()].index
#         if not missing_cells.empty:
#             print(f"  Identified {len(missing_cells)} cells missing in alternative ISO confidence score data.")
#             cells_to_remove = cells_to_remove.union(missing_cells)

#     # Filter adata_merged if any cells need to be removed
#     if not cells_to_remove.empty:
#         n_removed = len(cells_to_remove)
#         print(f"\nWarning: Removing {n_removed} unique cells from the current dataset due to missing data in at least one alternative source.")
#         cells_to_keep = adata_merged.obs_names.difference(cells_to_remove)
#         adata_merged = adata_merged[cells_to_keep.values, :].copy()
#         print(f"Removed {n_removed} cells. Current cell count: {adata_merged.n_obs}")

#         # Important: Re-index ALL aligned data AFTER filtering adata_merged
#         print("Re-indexing aligned alternative data to match filtered dataset...")
#         if aligned_alt_umap_df is not None:
#             aligned_alt_umap_df = aligned_alt_umap_df.reindex(adata_merged.obs_names)
#         if aligned_alt_leiden_series is not None:
#             aligned_alt_leiden_series = aligned_alt_leiden_series.reindex(adata_merged.obs_names)
#         if aligned_alt_dg_series is not None:
#             aligned_alt_dg_series = aligned_alt_dg_series.reindex(adata_merged.obs_names)
#         if aligned_alt_iso_series is not None:
#             aligned_alt_iso_series = aligned_alt_iso_series.reindex(adata_merged.obs_names)
#         if aligned_alt_dg_conf_series is not None:
#             aligned_alt_dg_conf_series = aligned_alt_dg_conf_series.reindex(adata_merged.obs_names)
#         if aligned_alt_iso_conf_series is not None:
#             aligned_alt_iso_conf_series = aligned_alt_iso_conf_series.reindex(adata_merged.obs_names)
#         print("Re-indexing complete.")
#     else:
#             print("No cells removed based on missing alternative data.")


#     # --- Add Aligned Data to Filtered adata_merged ---
#     print("\nAdding aligned alternative data to the final AnnData object...")

#     # Add Alternative UMAP
#     if aligned_alt_umap_df is not None:
#         if aligned_alt_umap_df.shape[0] != adata_merged.n_obs:
#             print(f"  Warning: UMAP dimension mismatch after filtering/re-indexing. Skipping UMAP addition.")
#         # Check for NaNs that might remain *after* reindexing (shouldn't happen if filter logic is correct)
#         elif aligned_alt_umap_df['alt_UMAP1'].isna().any():
#                 print(f"  Warning: NaNs found in UMAP data after filtering/re-indexing. Skipping UMAP addition.")
#         else:
#             adata_merged.obsm[alt_umap_key] = aligned_alt_umap_df[['alt_UMAP1', 'alt_UMAP2']].values
#             print(f"  Added alternative UMAP coordinates to adata_merged.obsm['{alt_umap_key}']")
#     elif alt_umap_present:
#             print(f"  Skipping UMAP addition (likely due to initial absence or previous error).")

#     # Add Alternative Leiden Clustering
#     if aligned_alt_leiden_series is not None:
#         if aligned_alt_leiden_series.shape[0] != adata_merged.n_obs:
#                 print(f"  Warning: Leiden dimension mismatch after filtering/re-indexing. Skipping Leiden addition.")
#         elif aligned_alt_leiden_series.isna().any():
#             print(f"  Warning: NaNs found in Leiden data after filtering/re-indexing. Filling with 'Unknown'.")
#             adata_merged.obs[new_alt_leiden_key] = pd.Categorical(aligned_alt_leiden_series.fillna('Unknown'))
#             print(f"  Added alternative Leiden clustering to adata_merged.obs['{new_alt_leiden_key}'] (NaNs filled)")
#         else:
#                 adata_merged.obs[new_alt_leiden_key] = pd.Categorical(aligned_alt_leiden_series)
#                 print(f"  Added alternative Leiden clustering to adata_merged.obs['{new_alt_leiden_key}']")
#     elif alt_leiden_present:
#             print(f"  Skipping Leiden addition (likely due to initial absence or previous error).")

#     # Add Alternative DG majority voting
#     if aligned_alt_dg_series is not None:
#         if aligned_alt_dg_series.shape[0] != adata_merged.n_obs:
#             print(f"  Warning: DG majority voting dimension mismatch after filtering/re-indexing. Skipping DG MV addition.")
#         else:
#             try:
#                 # Fill NaNs with 'Unassigned' (which should be an existing category)
#                 filled_dg_series = aligned_alt_dg_series.fillna('Unassigned')

#                 # Get original categories from IN-MEMORY data
#                 original_dg_categories = adata_alt_obs_in_memory['DG_majority_voting'].cat.categories
#                 if 'Unassigned' not in original_dg_categories:
#                      print(f"  Warning: 'Unassigned' not found in original DG categories. Adding it.")
#                      # This ideally shouldn't happen based on user data, but as a fallback:
#                      original_dg_categories = original_dg_categories.tolist() + ['Unassigned']

#                 # Create the final categorical series with explicit original categories
#                 final_dg_categorical = pd.Categorical(
#                     filled_dg_series,
#                     categories=original_dg_categories,
#                     ordered=False
#                 )

#                 # Assign the correctly typed Series
#                 adata_merged.obs[new_alt_dg_key] = final_dg_categorical
#                 print(f"  Added alternative DG majority voting to adata_merged.obs['{new_alt_dg_key}'] (NaNs filled with 'Unassigned').")

#             except Exception as e:
#                  print(f"  Error processing or assigning DG majority voting: {e}")
#     elif alt_dg_present:
#             print(f"  Skipping DG MV addition (likely due to initial absence or previous error).")

#     # Add Alternative ISO majority voting
#     if aligned_alt_iso_series is not None:
#         if aligned_alt_iso_series.shape[0] != adata_merged.n_obs:
#             print(f"  Warning: ISO majority voting dimension mismatch after filtering/re-indexing. Skipping ISO MV addition.")
#         else:
#             try:
#                 # Fill NaNs with 'Unassigned' (which should be an existing category)
#                 filled_iso_series = aligned_alt_iso_series.fillna('Unassigned')

#                 # Get original categories from IN-MEMORY data
#                 original_iso_categories = adata_alt_obs_in_memory['ISO_majority_voting'].cat.categories
#                 if 'Unassigned' not in original_iso_categories:
#                      print(f"  Warning: 'Unassigned' not found in original ISO categories. Adding it.")
#                      original_iso_categories = original_iso_categories.tolist() + ['Unassigned']

#                 # Create the final categorical series with explicit original categories
#                 final_iso_categorical = pd.Categorical(
#                     filled_iso_series,
#                     categories=original_iso_categories,
#                     ordered=False
#                 )

#                 # Assign the correctly typed Series
#                 adata_merged.obs[new_alt_iso_key] = final_iso_categorical
#                 print(f"  Added alternative ISO majority voting to adata_merged.obs['{new_alt_iso_key}'] (NaNs filled with 'Unassigned').")

#             except Exception as e:
#                  print(f"  Error processing or assigning ISO majority voting: {e}")
#     elif alt_iso_present:
#             print(f"  Skipping ISO MV addition (likely due to initial absence or previous error).")

#     # Add Alternative DG confidence score (using aligned_alt_dg_conf_series)
#     if aligned_alt_dg_conf_series is not None:
#         if aligned_alt_dg_conf_series.shape[0] != adata_merged.n_obs:
#             print(f"  Warning: DG confidence score dimension mismatch after filtering/re-indexing. Skipping DG Conf addition.")
#         else:
#             # Fill potential remaining NaNs with 0
#             adata_merged.obs[new_alt_dg_conf_key] = aligned_alt_dg_conf_series.fillna(0)
#             print(f"  Added alternative DG confidence score to adata_merged.obs['{new_alt_dg_conf_key}'] (NaNs filled with 0)")
#     elif alt_dg_conf_present:
#             print(f"  Skipping DG Conf addition (likely due to initial absence or previous error).")

#     # Add Alternative ISO confidence score (using aligned_alt_iso_conf_series)
#     if aligned_alt_iso_conf_series is not None:
#         if aligned_alt_iso_conf_series.shape[0] != adata_merged.n_obs:
#             print(f"  Warning: ISO confidence score dimension mismatch after filtering/re-indexing. Skipping ISO Conf addition.")
#         else:
#             adata_merged.obs[new_alt_iso_conf_key] = aligned_alt_iso_conf_series.fillna(0)
#             print(f"  Added alternative ISO confidence score to adata_merged.obs['{new_alt_iso_conf_key}'] (NaNs filled with 0)")
#     elif alt_iso_conf_present:
#             print(f"  Skipping ISO Conf addition (likely due to initial absence or previous error).")

#     # Plot the alternative UMAP, colored by the *alternative* Leiden clustering if available
#     if alt_umap_present and alt_umap_key in adata_merged.obsm:
#         print(f"Plotting alternative UMAP ({alt_umap_key})...")
#         # Use the newly added alternative leiden key for coloring, if it exists and was added
#         plot_color_key = new_alt_leiden_key if new_alt_leiden_key in adata_merged.obs else None
#         if plot_color_key:
#                 sc.pl.embedding(adata_merged, basis=alt_umap_key, color=plot_color_key,
#                                 legend_loc='on data', save=f"_umap_alt_colored_by_{plot_color_key}.png", show=True,
#                                 title=f'Alternative UMAP colored by {plot_color_key}')
#         else:
#                 # Plot without color if alt leiden key is missing or wasn't added
#                 sc.pl.embedding(adata_merged, basis=alt_umap_key,
#                                 save=f"_umap_alt_basic.png", show=True,
#                                 title='Alternative UMAP')
#     elif alt_umap_present:
#             print(f"Skipping alternative UMAP plot because '{alt_umap_key}' was not successfully added to adata_merged.obsm.")

# # %%
# adata_merged

# # %%
# # Save the processed merged dataset (now potentially with alt UMAP)
# output_file = os.path.join(OUTPUT_DIR, 'merged_raw_processed.h5ad')
# print(f"\nSaving processed merged dataset to {output_file}")
# try:
#     # Ensure the output directory exists
#     os.makedirs(os.path.dirname(output_file), exist_ok=True)
#     adata_merged.write(Path(output_file))
#     print("Successfully saved the AnnData object.")
# except Exception as e:
#     print(f"Error saving AnnData object: {e}")

# print("\nScript finished!")

# # %%
# sc.pl.umap(adata_merged, color='leiden_0.8', legend_loc='on data', show=True)