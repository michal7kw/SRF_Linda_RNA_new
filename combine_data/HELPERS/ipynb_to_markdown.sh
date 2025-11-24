#!/bin/bash

# Create the markdown directory if it doesn't exist
mkdir -p markdown

jupyter nbconvert --to markdown raw_0_merge_raw_and_process.ipynb --output-dir=markdown
jupyter nbconvert --to markdown raw_1_annotate_merged_raw.ipynb --output-dir=markdown
jupyter nbconvert --to markdown raw_5_clean_data.ipynb --output-dir=markdown
jupyter nbconvert --to markdown raw_6_add_mapmycells_annotations.ipynb --output-dir=markdown
jupyter nbconvert --to markdown raw_9_finalize_annotation.ipynb --output-dir=markdown
jupyter nbconvert --to markdown raw_10_clusters_similiarty.ipynb --output-dir=markdown
jupyter nbconvert --to markdown raw_10_clusters_similiarty.ipynb --output-dir=markdown
jupyter nbconvert --to markdown raw_16_dge_analysis_final_cell_types_level_one.ipynb --output-dir=markdown
jupyter nbconvert --to markdown raw_16_dge_analysis_merged_GC_level_two.ipynb --output-dir=markdown