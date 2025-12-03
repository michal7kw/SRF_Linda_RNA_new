#!/bin/bash

# Master script to run all CRE analysis pipelines IN PARALLEL
# All jobs are submitted simultaneously via sbatch (non-blocking)

# Base directory
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/integration_scripts/multiome_modular_pipeline/CREs_from_literature"

# List of directories to process (all run in parallel)
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
echo "Master CRE Pipeline Launcher (PARALLEL)"
echo "Started at: $(date)"
echo "=============================================="
echo ""

# Submit all jobs in parallel
for DIR in "${DIRECTORIES[@]}"; do
    FULL_PATH="${BASE_DIR}/${DIR}"

    if [[ -d "${FULL_PATH}" ]] && [[ -f "${FULL_PATH}/0_RUN_ALL.sh" ]]; then
        cd "${FULL_PATH}"
        echo -n "${DIR}: "
        sbatch 0_RUN_ALL.sh
    else
        echo "${DIR}: SKIPPED (missing directory or script)"
    fi
done

cd "${BASE_DIR}"
echo ""
echo "All jobs submitted in parallel at: $(date)"
echo "=============================================="
