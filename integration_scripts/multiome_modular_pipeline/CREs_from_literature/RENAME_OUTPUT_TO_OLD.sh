#!/bin/bash

# Script to rename "output" folders to "output_old" in all CRE pipeline directories

# Base directory
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature"

# List of directories to process
DIRECTORIES=(
    "CREs_encode_all"
    "CREs_encode_GABA_specific"
    "CREs_encode_paper_intersection"
    "CREs_GABA_DEGs_encode"
    "CREs_GABA_DEGs_paper"
    "CREs_paper_exclusive"
    "CREs_paper_specific"
    "CREs_splicing_differential_all_celltypes"
    "CREs_splicing_directional"
    "CREs_splicing_genes_encode"
    "CREs_splicing_genes_paper"
)

echo "=============================================="
echo "Rename output -> output_old"
echo "Started at: $(date)"
echo "=============================================="
echo ""

for DIR in "${DIRECTORIES[@]}"; do
    OUTPUT_PATH="${BASE_DIR}/${DIR}/output"
    OUTPUT_OLD_PATH="${BASE_DIR}/${DIR}/output_old"

    echo -n "${DIR}: "

    if [[ ! -d "${OUTPUT_PATH}" ]]; then
        echo "SKIPPED (no output folder)"
    elif [[ -d "${OUTPUT_OLD_PATH}" ]]; then
        echo "SKIPPED (output_old already exists)"
    else
        mv "${OUTPUT_PATH}" "${OUTPUT_OLD_PATH}"
        echo "RENAMED"
    fi
done

echo ""
echo "Completed at: $(date)"
echo "=============================================="
