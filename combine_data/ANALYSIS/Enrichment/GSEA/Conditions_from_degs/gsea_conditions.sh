#!/bin/bash

# This script assumes that 'mamba' is available and the shell is configured
# to use it (e.g., by running 'conda init bash' or 'mamba init bash').

# Fail on any error
set -e

# Initialize mamba for the current shell session
eval "$(mamba.exe shell hook --shell bash)"

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Activate conda environment
echo "Activating conda environment 'bioinf'..."
if ! mamba activate bioinf; then
    echo "Could not activate conda environment 'bioinf'"
    echo "Please make sure it exists and your shell is configured for conda/mamba."
    exit 1
fi

# Log directory will be relative to the script's location
LOG_DIR="${SCRIPT_DIR}/../logs/gsea_runs_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$LOG_DIR"

# Since the scripts are in the same directory as this .sh file, we'll cd there.
cd "$SCRIPT_DIR"

echo "Running GSEA analysis scripts in parallel. Logs will be saved to ${LOG_DIR}"

run_script() {
    local genotype=$1
    local merge_clusters=$2
    local cluster_name=$3
    local cluster_column=$4

    if [ -n "$cluster_name" ] && [ -n "$cluster_column" ]; then
        local sanitized_cluster_name=$(echo "$cluster_name" | sed 's/ /_/g' | sed 's/\//_/g')
        local log_suffix="${genotype}_cluster_${sanitized_cluster_name}"
        echo "  - Starting gsea_conditions_unified.py for genotype=${genotype}, cluster=${cluster_name}"
        python gsea_conditions_unified.py --genotype "${genotype}" --merge_clusters "False" --cluster_name "${cluster_name}" --cluster_column "${cluster_column}" > "${LOG_DIR}/gsea_conditions_${log_suffix}.log" 2>&1 &
    else
        local log_suffix="${genotype}_merge_${merge_clusters}"
        echo "  - Starting gsea_conditions_unified.py for genotype=${genotype}, merge_clusters=${merge_clusters}"
        python gsea_conditions_unified.py --genotype "${genotype}" --merge_clusters "${merge_clusters}" > "${LOG_DIR}/gsea_conditions_${log_suffix}.log" 2>&1 &
    fi
}

# Run for all 6 cases:
# Genotypes: Emx1, Nestin, both
# Merge Clusters: True, False
run_script "Emx1" "False"
run_script "Nestin" "False"
run_script "both" "False"
run_script "Emx1" "True"
run_script "Nestin" "True"
run_script "both" "True"

# Run for GABA cluster from cell_type_L1 for all genotypes
echo "Starting GSEA analysis for GABA cluster..."
run_script "Emx1" "False" "GABA" "cell_type_L1"
run_script "Nestin" "False" "GABA" "cell_type_L1"
run_script "both" "False" "GABA" "cell_type_L1"

echo "Waiting for all scripts to complete..."
wait

echo "All GSEA analysis scripts have completed." 