#!/bin/bash
# Master script to run Striatum_EP analysis pipeline
# Submits all steps as dependent SLURM jobs

echo "========================================================================"
echo "RUNNING STRIATUM_EP CRE ANALYSIS PIPELINE"
echo "========================================================================"
echo "Started: $(date)"

# 1. Submit Extraction Job
echo "Step 1: Submitting Extraction Job..."
JOB_ID_EXTRACT=$(sbatch --parsable 1_extract_Striatum_EP.sh)
echo "Submitted Extraction Job: $JOB_ID_EXTRACT"

# 2. Submit Heatmap Job (dependent on Extraction)
echo "Step 2: Submitting Heatmap Job..."
JOB_ID_HEATMAP=$(sbatch --parsable --dependency=afterok:$JOB_ID_EXTRACT 4_create_heatmaps.sh)
echo "Submitted Heatmap Job: $JOB_ID_HEATMAP"

# 3. Submit Visualization Job (dependent on Heatmap)
echo "Step 3: Submitting Visualization Job..."
JOB_ID_VIZ=$(sbatch --parsable --dependency=afterok:$JOB_ID_HEATMAP 5_visualize_comparisons.sh)
echo "Submitted Visualization Job: $JOB_ID_VIZ"

echo "========================================================================"
echo "Pipeline submitted successfully!"
echo "Monitor jobs with: squeue -u $USER"
echo "========================================================================"
