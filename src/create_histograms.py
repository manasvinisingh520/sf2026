"""
Load region-specific matrix/gene/cell files, filter cells by mitochondrial percentage
from metadata, then bin filtered cells by Path..Group. and output min-max normalized
gene expression matrices per bin.
"""

import numpy as np
import pandas as pd
import argparse
from pathlib import Path
from utils import read_mtx_file, read_excel_columns


def filter_cells(matrix, cell_names, metadata, mito_max=15.0):
    filtered_metadata = metadata[metadata['cell_annotation'].isin(cell_names)].copy()
    filtered_metadata = filtered_metadata.set_index('cell_annotation').reindex(cell_names)
    n_before_filter = len(filtered_metadata)

    # Filter by mitochondrial percentage only
    mito_mask = (filtered_metadata['percent.mito'] < mito_max).fillna(False)
    combined_mask = mito_mask
    filtered_metadata = filtered_metadata[combined_mask].copy()
    
    print(f"  Before: {n_before_filter:,} cells")
    print(f"  After: {len(filtered_metadata):,} cells")
    
    # Get cell indices to keep (cells that passed filters)
    cells_to_keep = filtered_metadata.index.dropna().tolist()
    cell_indices_to_keep = [i for i, name in enumerate(cell_names) if name in cells_to_keep]
    
    # Filter matrix and cell names
    filtered_matrix = matrix[:, cell_indices_to_keep]
    filtered_cell_names = [cell_names[i] for i in cell_indices_to_keep]
    
    return filtered_matrix, filtered_cell_names, filtered_metadata


def min_max_normalize_rows_sparse(matrix):
    """
    Min-max normalize each gene (row) across cells.
    Returns a CSR sparse matrix where each row is scaled to [0,1].
    """
    from scipy.sparse import coo_matrix

    mat = matrix.tocoo(copy=True).astype(float)
    row_max = np.array(matrix.max(axis=1)).ravel()
    row_min = np.array(matrix.min(axis=1)).ravel()
    denom = row_max - row_min
    denom[denom == 0] = 1.0  # avoid divide-by-zero; rows with zero range become all zeros

    rows = mat.row
    mat.data = (mat.data - row_min[rows]) / denom[rows]

    # entries that were zero but above row_min stay zero implicitly; rows with zero range stay zero
    return mat.tocsr()


def bin_by_path_group(normalized_matrix, cell_names, filtered_metadata):
    """
    Split normalized matrix into bins (Cond1-Cond4) based on Path..Group.
    Returns dict of bin_name -> (matrix, cell_names_in_bin, group_label).
    """
    bins = {}
    groups = filtered_metadata['Path..Group.'].fillna('NA')
    unique_groups = list(dict.fromkeys(groups))  # preserve order
    cond_labels = ['Cond1', 'Cond2', 'Cond3', 'Cond4']
    for idx, group in enumerate(unique_groups[:4]):
        bin_cells_idx = [i for i, name in enumerate(cell_names) if groups.iloc[i] == group]
        bins[cond_labels[idx]] = (
            normalized_matrix[:, bin_cells_idx],
            [cell_names[i] for i in bin_cells_idx],
            group,
        )
    return bins


def main():
    parser = argparse.ArgumentParser(
        description='Load region-specific files, filter cells by mito %, then bin and normalize',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '-region',
        type=str,
        default='EC',
        choices=['EC', 'ITG', 'PFC', 'V1', 'V2'],
        help='Region name to select correct tab in metadata (default: EC)'
    )
    parser.add_argument(
        '-mito_max',
        type=float,
        default=15.0,
        help='Maximum mitochondrial percentage threshold (default: 15.0)'
    )
    parser.add_argument(
        '-data_dir',
        type=str,
        default='data',
        help='Base directory containing region files (default: data)'
    )
    parser.add_argument(
        '-metadata',
        type=str,
        default='data/2025-11-16_Astrocytes_metadata.xlsx',
        help='Path to metadata Excel file (default: shared workbook)'
    )
    
    args = parser.parse_args()
    
    region = args.region.upper()
    data_dir = Path(args.data_dir)
    metadata_path = Path(args.metadata).resolve()

    # Derive file paths from region
    base_prefix = f"2025-10-22_Astrocytes_{region}"
    mtx_path = (data_dir / f"{base_prefix}_matrix.mtx").resolve()
    row_path = (data_dir / f"{base_prefix}_row_annotation.txt").resolve()
    col_path = (data_dir / f"{base_prefix}_cell_annotation.txt").resolve()

    # Map region to tab index (0-based)
    region_to_tab = {
        "EC": 0,
        "ITG": 1,
        "PFC": 2,
        "V2": 3,
        "V1": 4,
    }
    tab_index = region_to_tab.get(region)

    # Load metadata (need Path..Group. for binning)
    metadata = read_excel_columns(
        str(metadata_path),
        columns=['cell_annotation', 'Median.UMI.Counts.per.Cell', 'percent.mito', 'Path..Group.'],
        sheet_name=tab_index
    )
    print(f"Loaded metadata: {metadata.shape}")
    
    # Read the matrix file
    matrix, gene_names, cell_names = read_mtx_file(
        mtx_path=str(mtx_path),
        row_annotation_path=str(row_path),
        col_annotation_path=str(col_path),
        transpose=True  # Transpose so rows=genes, cols=cells
    )
    
    print(f"Loaded: {len(gene_names):,} genes × {len(cell_names):,} cells")
    
    # Filter cells using metadata
    print("\nFiltering cells...")
    filtered_matrix, filtered_cell_names, filtered_metadata = filter_cells(
        matrix,
        cell_names,
        metadata,
        mito_max=args.mito_max
    )

    # Min-max normalize per gene across filtered cells
    print("\nNormalizing (min-max per gene)...")
    norm_matrix = min_max_normalize_rows_sparse(filtered_matrix)

    # Bin by Path..Group. into Cond1-Cond4
    bins = bin_by_path_group(norm_matrix, filtered_cell_names, filtered_metadata)

    # Report bin sizes
    for cond, (mat_bin, cells_bin, group_label) in bins.items():
        print(f"{cond} ({group_label}): genes x cells = {mat_bin.shape}")
    
    return norm_matrix, gene_names, filtered_cell_names, bins


if __name__ == "__main__":
    main()

