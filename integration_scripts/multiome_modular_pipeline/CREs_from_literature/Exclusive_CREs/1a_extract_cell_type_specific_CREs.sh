#!/bin/bash
#SBATCH --job-name=1a_extract_cell_type_specific_CREs
#SBATCH --output=logs/1a_extract_cell_type_specific_CREs.log
#SBATCH --error=logs/1a_extract_cell_type_specific_CREs.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=workq

################################################################################
# Extract Cell-Type-SPECIFIC CREs (Mutually Exclusive)
#
# This script extracts CREs that are EXCLUSIVE to each cell type:
# - GABA-specific: Active in GABA neurons but NOT in Excitatory neurons
# - Excitatory-specific: Active in Excitatory neurons but NOT in GABA neurons
#
# Problem with previous approach:
# - 60% of CREs overlapped between GABA and Excitatory sets
# - This caused identical signal patterns in both
# - Not useful as negative control
#
# New approach:
# - Extract mutually exclusive CRE sets
# - Ensures true positive/negative control comparison
# - Expected: Clear visual difference in ATAC signal
#
# INPUT FILES:
# - data/table_1.xlsx: ENCODE data table 1
# - data/table_2.tsv: ENCODE data table 2
# - data/table_8.txt: ENCODE data table 8
# - extract_cell_type_specific_CREs.py: Python script to extract mutually exclusive CRE sets
#
# OUTPUT FILES:
# - output/GABA_specific_CREs.bed: BED file with GABA-specific CREs (mutually exclusive)
# - output/Excitatory_specific_CREs.bed: BED file with excitatory-specific CREs (mutually exclusive)
# - output/cell_type_specific_CREs_summary.txt: Summary report with overlap analysis
# - logs/0a_extract_cell_type_specific_CREs.log: Execution log file
# - logs/0a_extract_cell_type_specific_CREs.err: Error log file
#
################################################################################

echo "========================================================================"
echo "EXTRACT CELL-TYPE-SPECIFIC CREs (MUTUALLY EXCLUSIVE)"
echo "========================================================================"
echo "Started: $(date)"
echo ""

# Create logs and output directories
mkdir -p logs
mkdir -p output

# Activate conda environment with pandas
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate bigwig

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# Check if required data files exist
echo "Checking prerequisites..."
echo ""

MISSING=0

if [ ! -f "data/table_1.xlsx" ]; then
    echo "✗ Missing: data/table_1.xlsx"
    MISSING=$((MISSING + 1))
else
    echo "✓ Found: data/table_1.xlsx"
fi

if [ ! -f "data/table_2.tsv" ]; then
    echo "✗ Missing: data/table_2.tsv"
    MISSING=$((MISSING + 1))
else
    echo "✓ Found: data/table_2.tsv"
fi

if [ ! -f "data/table_8.txt" ]; then
    echo "✗ Missing: data/table_8.txt"
    MISSING=$((MISSING + 1))
else
    FILE_SIZE=$(du -h data/table_8.txt | cut -f1)
    echo "✓ Found: data/table_8.txt ($FILE_SIZE)"
fi

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "ERROR: Missing $MISSING required data files!"
    exit 1
fi

echo ""
echo "========================================================================"
echo "Running Python extraction script..."
echo "========================================================================"
echo ""

# Run the Python script
python 1a_extract_cell_type_specific_CREs.py

EXIT_CODE=$?

echo ""
echo "========================================================================"

if [ $EXIT_CODE -eq 0 ]; then
    echo "EXTRACTION COMPLETE!"
    echo "========================================================================"
    echo ""

    # Show results
    if [ -f "output/GABA_specific_CREs.bed" ] && [ -f "output/Excitatory_specific_CREs.bed" ]; then
        N_GABA=$(wc -l < output/GABA_specific_CREs.bed)
        N_EXCIT=$(wc -l < output/Excitatory_specific_CREs.bed)

        echo "Output files created:"
        echo "  ★ GABA-specific CREs: $N_GABA"
        echo "     → output/GABA_specific_CREs.bed"
        echo ""
        echo "  ★ Excitatory-specific CREs: $N_EXCIT"
        echo "     → output/Excitatory_specific_CREs.bed"
        echo ""
        echo "  ★ Summary report:"
        echo "     → output/cell_type_specific_CREs_summary.txt"
        echo ""

        # Verify no overlap
        echo "Verifying mutual exclusivity..."
        OVERLAP=$(comm -12 \
            <(awk '{print $1":"$2"-"$3}' output/GABA_specific_CREs.bed | sort) \
            <(awk '{print $1":"$2"-"$3}' output/Excitatory_specific_CREs.bed | sort) \
            | wc -l)

        if [ $OVERLAP -eq 0 ]; then
            echo "  ✓ VERIFIED: 0 overlapping coordinates (mutually exclusive)"
        else
            echo "  ✗ WARNING: $OVERLAP overlapping coordinates found!"
        fi
        echo ""

        # Show first few lines of each BED file
        echo "First 3 GABA-specific CREs:"
        head -3 output/GABA_specific_CREs.bed | column -t
        echo ""

        echo "First 3 Excitatory-specific CREs:"
        head -3 output/Excitatory_specific_CREs.bed | column -t
        echo ""
    fi

    echo "========================================================================"
    echo "NEXT STEPS"
    echo "========================================================================"
    echo ""
    echo "1. Review the summary report:"
    echo "   cat output/cell_type_specific_CREs_summary.txt"
    echo ""
    echo "2. Update your analysis scripts to use the new BED files:"
    echo ""
    echo "   For deepTools (1b_create_heatmaps_specific_CREs_deeptools.sh):"
    echo "   BED_GABA=\"output/GABA_specific_CREs.bed\""
    echo "   BED_EXCITATORY=\"output/Excitatory_specific_CREs.bed\""
    echo ""
    echo "   For Python (3_compare_cell_type_specific_CREs.py):"
    echo "   BED_GABA = \"output/GABA_specific_CREs.bed\""
    echo "   BED_EXCITATORY = \"output/Excitatory_specific_CREs.bed\""
    echo ""
    echo "3. Re-run your heatmap analysis:"
    echo "   sbatch 1b_create_heatmaps_specific_CREs_deeptools.sh"
    echo ""
    echo "4. Expected results:"
    echo "   - GABA-specific CREs: HIGH signal in GABA ATAC samples"
    echo "   - Excitatory-specific CREs: LOW signal in GABA ATAC samples"
    echo "   - Clear visual difference demonstrating cell-type specificity!"
    echo ""
    echo "Completed: $(date)"
    echo "========================================================================"
else
    echo "EXTRACTION FAILED!"
    echo "========================================================================"
    echo ""
    echo "Check the error log for details:"
    echo "  logs/extract_specific_CREs_${SLURM_JOB_ID}.err"
    echo ""
    exit 1
fi
