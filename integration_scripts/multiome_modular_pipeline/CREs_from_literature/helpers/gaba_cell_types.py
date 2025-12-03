#!/usr/bin/env python3
"""
Hippocampal GABAergic Cell Types - Shared Module

This module provides the correct list of hippocampal GABAergic cell types
for filtering Table 16/Table 8 SubType columns.

METHODOLOGY (from helpers/explore_data/identify_hippocampal_gaba_types.py):
1. Hippocampal samples identified from Table 1 (Major Region == 'HPF')
2. Cells from those samples filtered by Class == 'GABA' in Table 2
3. Extract unique CellType values -> 46 GABAergic cell types

IMPORTANT: This EXCLUDES glutamatergic cell types that were incorrectly
matched by previous keyword-based approaches:
- CA1GL1-3, CA3GL1-6: Glutamatergic pyramidal neurons
- DGGR: Dentate gyrus granule cells (glutamatergic)
- PVM: Perivascular macrophages (non-neuronal)
- RGDG: Radial glia

Usage:
    from helpers.gaba_cell_types import HIPPOCAMPAL_GABA_CELLTYPES, is_gaba_subtype

    # Filter dataframe
    df['is_gaba'] = df['SubType'].apply(is_gaba_subtype)
    gaba_df = df[df['is_gaba']]
"""

import pandas as pd

# Explicit list of 46 GABAergic cell types found in hippocampal samples
# Derived from: helpers/explore_data/identify_hippocampal_gaba_types.py
HIPPOCAMPAL_GABA_CELLTYPES = {
    # Parvalbumin interneurons (PV+)
    'PVGA1', 'PVGA2', 'PVGA3', 'PVGA4', 'PVGA5', 'PVGA6', 'PVGA7',
    # Somatostatin interneurons (SST+)
    'SSTGA1', 'SSTGA2', 'SSTGA3', 'SSTGA4', 'SSTGA5', 'SSTGA6',
    'SSTGA7', 'SSTGA8', 'SSTGA9', 'SSTGA10',
    # VIP interneurons
    'VIPGA1', 'VIPGA2', 'VIPGA3', 'VIPGA4',
    # Lamp5 interneurons
    'LAMGA1', 'LAMGA2', 'LAMGA3', 'LAMGA4',
    # Dentate gyrus neuroblasts (GABAergic immature neurons)
    'DGNBL1', 'DGNBL2',
    # Lateral septal complex GABAergic
    'LSXGA3', 'LSXGA4', 'LSXGA5', 'LSXGA7',
    # Medial septal GABAergic
    'MSGA1', 'MSGA2', 'MSGA4', 'MSGA6', 'MSGA7', 'MSGA8', 'MSGA9',
    'MSGA11', 'MSGA12',
    # Other minor populations
    'CNUGA',      # Cerebellar nuclei GABAergic (sparse in hippocampus)
    'OBNBL',      # Olfactory bulb neuroblasts
    'STRGA2', 'STRGA3',  # Striatal GABAergic
    'D2MSN2', 'D2MSN3',  # D2 medium spiny neurons (very sparse)
}

# Excluded cell types (for reference/documentation)
EXCLUDED_GLUTAMATERGIC = {
    # CA1 glutamatergic pyramidal neurons
    'CA1GL1', 'CA1GL2', 'CA1GL3',
    # CA3 glutamatergic pyramidal neurons
    'CA3GL1', 'CA3GL2', 'CA3GL3', 'CA3GL4', 'CA3GL5', 'CA3GL6',
    # Dentate gyrus granule cells (glutamatergic)
    'DGGR',
    # Radial glia
    'RGDG',
    # Non-neuronal cells
    'PVM',  # Perivascular macrophages
}


def is_gaba_subtype(subtype):
    """
    Check if SubType is a hippocampal GABAergic cell type.
    Uses EXACT matching against the curated list of 46 cell types.

    This excludes:
    - Glutamatergic neurons (CA1GL, CA3GL, DGGR, etc.)
    - Non-neuronal cells (PVM, etc.)

    Args:
        subtype: SubType value from Table 16 or similar

    Returns:
        bool: True if subtype is a hippocampal GABAergic cell type
    """
    if pd.isna(subtype):
        return False
    return str(subtype) in HIPPOCAMPAL_GABA_CELLTYPES


def get_gaba_celltypes_summary():
    """Return a formatted summary of GABA cell types for documentation."""
    categories = {
        'Parvalbumin (PV+)': ['PVGA1', 'PVGA2', 'PVGA3', 'PVGA4', 'PVGA5', 'PVGA6', 'PVGA7'],
        'Somatostatin (SST+)': ['SSTGA1', 'SSTGA2', 'SSTGA3', 'SSTGA4', 'SSTGA5', 'SSTGA6',
                                'SSTGA7', 'SSTGA8', 'SSTGA9', 'SSTGA10'],
        'VIP interneurons': ['VIPGA1', 'VIPGA2', 'VIPGA3', 'VIPGA4'],
        'Lamp5 interneurons': ['LAMGA1', 'LAMGA2', 'LAMGA3', 'LAMGA4'],
        'DG neuroblasts': ['DGNBL1', 'DGNBL2'],
        'Lateral septal': ['LSXGA3', 'LSXGA4', 'LSXGA5', 'LSXGA7'],
        'Medial septal': ['MSGA1', 'MSGA2', 'MSGA4', 'MSGA6', 'MSGA7', 'MSGA8', 'MSGA9',
                          'MSGA11', 'MSGA12'],
        'Other (sparse)': ['CNUGA', 'OBNBL', 'STRGA2', 'STRGA3', 'D2MSN2', 'D2MSN3'],
    }

    lines = [
        f"Hippocampal GABAergic cell types ({len(HIPPOCAMPAL_GABA_CELLTYPES)} total):",
        "Method: EXACT cell type matching (not keyword-based)",
        "",
    ]

    for category, types in categories.items():
        lines.append(f"  {category}: {', '.join(types)}")

    lines.extend([
        "",
        "EXCLUDED (glutamatergic/non-neuronal):",
        f"  {', '.join(sorted(EXCLUDED_GLUTAMATERGIC))}",
    ])

    return '\n'.join(lines)


if __name__ == '__main__':
    # Print summary when run directly
    print(get_gaba_celltypes_summary())
    print(f"\nTotal cell types: {len(HIPPOCAMPAL_GABA_CELLTYPES)}")
