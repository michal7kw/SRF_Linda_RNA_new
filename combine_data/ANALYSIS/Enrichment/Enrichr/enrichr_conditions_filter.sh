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
LOG_DIR="${SCRIPT_DIR}/../logs/enrichr_filter_runs_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$LOG_DIR"

# Since the scripts are in the same directory as this .sh file, we'll cd there.
cd "$SCRIPT_DIR"

echo "Running Enrichr filter scripts. Logs will be saved to ${LOG_DIR}"

run_script() {
    local genotype=$1
    local merge_clusters=$2
    local focus=$3
    local cluster_name=$4
    local cluster_column=$5

    local params="--genotype ${genotype} --merge_clusters ${merge_clusters} --focus ${focus}"
    local log_suffix="${genotype}_focus_${focus}"

    if [ -n "$cluster_name" ] && [ -n "$cluster_column" ]; then
        local sanitized_cluster_name=$(echo "$cluster_name" | sed 's/ /_/g' | sed 's/\//_/g')
        params+=" --cluster_name \"${cluster_name}\" --cluster_column \"${cluster_column}\""
        log_suffix+="_cluster_${sanitized_cluster_name}"
        echo "  - Starting enrichr_conditions_filter.py for genotype=${genotype}, focus=${focus}, cluster=${cluster_name}"
    else
        log_suffix+="_merge_${merge_clusters}"
        echo "  - Starting enrichr_conditions_filter.py for genotype=${genotype}, focus=${focus}, merge_clusters=${merge_clusters}"
    fi

    # Use eval to correctly handle parameters with spaces
    eval "python enrichr_conditions_filter.py ${params}" > "${LOG_DIR}/enrichr_conditions_filter_${log_suffix}.log" 2>&1 &
}

# Run for all 6 original cases with 'neuro' focus
# echo "Running focused analysis for Neurodegeneration & Cell Death..."
run_script "Emx1" "False" "neuro"
# run_script "Nestin" "False" "neuro"
# run_script "both" "False" "neuro"
# run_script "Emx1" "True" "neuro"
# run_script "Nestin" "True" "neuro"
# run_script "both" "True" "neuro"

# Run for GABA cluster with 'splicing' focus
# echo "Running focused analysis for RNA Splicing on GABA cluster..."
# run_script "Emx1" "False" "splicing" "GABA" "cell_type_L1"
# run_script "Nestin" "False" "splicing" "GABA" "cell_type_L1"
# run_script "both" "False" "splicing" "GABA" "cell_type_L1"

echo "Waiting for all scripts to complete..."
wait

echo "All Enrichr filter scripts have completed."