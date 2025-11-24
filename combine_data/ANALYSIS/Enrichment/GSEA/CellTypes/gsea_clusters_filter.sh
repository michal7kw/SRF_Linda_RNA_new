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
LOG_DIR="${SCRIPT_DIR}/../logs/gsea_focus_runs_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$LOG_DIR"

# Since the scripts are in the same directory as this .sh file, we'll cd there.
cd "$SCRIPT_DIR"

echo "Running GSEA focus analysis scripts in parallel. Logs will be saved to ${LOG_DIR}"

run_script() {
    local condition=$1
    local genotype=$2
    local merge_clusters=$3
    local log_suffix="${condition}_${genotype}_merge_${merge_clusters}"
    echo "  - Starting gsea_clusters_filter_unified.py for condition=${condition}, genotype=${genotype}, merge_clusters=${merge_clusters}"
    python gsea_clusters_filter_unified.py --condition "${condition}" --genotype "${genotype}" --merge_clusters "${merge_clusters}" > "${LOG_DIR}/gsea_clusters_filter_${log_suffix}.log" 2>&1 &
}

# Run for all 12 cases:
run_script "Control" "both" "True"
run_script "Control" "Emx1" "True"
run_script "Control" "Nestin" "True"
run_script "Mutant" "both" "True"
run_script "Mutant" "Emx1" "True"
run_script "Mutant" "Nestin" "True"
run_script "Control" "both" "False"
run_script "Control" "Emx1" "False"
run_script "Control" "Nestin" "False"
run_script "Mutant" "both" "False"
run_script "Mutant" "Emx1" "False"
run_script "Mutant" "Nestin" "False"

echo "Waiting for all scripts to complete..."
wait

echo "All GSEA focus analysis scripts have completed." 