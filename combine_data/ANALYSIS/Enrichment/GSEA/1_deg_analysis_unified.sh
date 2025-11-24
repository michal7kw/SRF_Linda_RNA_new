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
LOG_DIR="${SCRIPT_DIR}/../logs/deg_runs_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$LOG_DIR"

# Since the scripts are in the same directory as this .sh file, we'll cd there.
cd "$SCRIPT_DIR"

echo "Running DEG analysis scripts in parallel. Logs will be saved to ${LOG_DIR}"

run_script() {
    local genotype=$1
    local merge_clusters=$2
    local log_suffix="${genotype}_merge_${merge_clusters}"
    echo "  - Starting deg_analysis_unified.py for genotype=${genotype}, merge_clusters=${merge_clusters}"
    python deg_analysis_unified.py --genotype "${genotype}" --merge_clusters "${merge_clusters}" > "${LOG_DIR}/deg_analysis_${log_suffix}.log" 2>&1 &
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

echo "Waiting for all scripts to complete..."
wait

echo "All DEG analysis scripts have completed."