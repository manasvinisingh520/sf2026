"""
Utility functions for reading and processing single-cell data.

This module provides helper functions for:
- Reading .mtx files (Matrix Market format)
- Creating AnnData objects from matrices
- Reading Excel files into pandas DataFrames
- Filtering cells based on metadata criteria
- Getting file paths for region-specific data
"""

import scipy.io
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, List, Dict, Union, Tuple


def read_mtx_file(mtx_path, row_annotation_path=None, col_annotation_path=None, 
                  transpose=True, make_gene_names_unique=True):
    """
    Read a .mtx file and optionally load row/column annotations.
    
    Parameters:
    -----------
    mtx_path : str
        Path to the .mtx file
    row_annotation_path : str, optional
        Path to file with row names (genes/features). 
        Each line should contain one name.
    col_annotation_path : str, optional
        Path to file with column names (cell barcodes).
        Each line should contain one name.
    transpose : bool, default=True
        If True, transpose the matrix so rows=genes, cols=cells.
        MTX format typically stores cells as rows and genes as columns.
    make_gene_names_unique : bool, default=True
        If True, ensure gene names are unique by appending _1, _2, etc.
    
    Returns:
    --------
    matrix : scipy.sparse.csr_matrix or scipy.sparse.csc_matrix
        Sparse matrix representation
    gene_names : array-like, optional
        Gene/feature names if row_annotation_path provided
    cell_names : array-like, optional
        Cell barcode names if col_annotation_path provided
    """
    
    # Read the .mtx file using scipy
    # This returns a COO (Coordinate) format sparse matrix
    matrix = scipy.io.mmread(mtx_path)
    print(f"Loaded matrix shape: {matrix.shape}")
    print(f"Number of non-zero entries: {matrix.nnz:,}")
    print(f"Sparsity: {(1 - matrix.nnz / (matrix.shape[0] * matrix.shape[1])) * 100:.2f}%")
    
    # Transpose if needed (MTX format often has cells as rows, genes as columns)
    if transpose:
        matrix = matrix.T
        print(f"Transposed matrix shape: {matrix.shape}")
    
    # Convert to CSR format for efficient row operations (or CSC for column operations)
    # CSR is typically preferred for most operations
    matrix = csr_matrix(matrix)
    
    # Load row annotations (genes/features) if provided
    gene_names = None
    if row_annotation_path:
        # Read gene names from file
        # Handle different formats (quoted strings, CSV, TSV, etc.)
        try:
            with open(row_annotation_path, 'r') as f:
                gene_names = [line.strip().strip('"') for line in f]
            
            # Handle unique gene names (common issue in scRNA-seq)
            if make_gene_names_unique:
                unique_names = []
                name_counts = {}
                for name in gene_names:
                    if name in name_counts:
                        name_counts[name] += 1
                        unique_names.append(f"{name}_{name_counts[name]}")
                    else:
                        name_counts[name] = 0
                        unique_names.append(name)
                gene_names = unique_names
            
            print(f"Loaded {len(gene_names)} gene names")
            
        except Exception as e:
            print(f"Warning: Could not load row annotations: {e}")
    
    # Load column annotations (cell barcodes) if provided
    cell_names = None
    if col_annotation_path:
        try:
            with open(col_annotation_path, 'r') as f:
                cell_names = [line.strip().strip('"') for line in f]
            print(f"Loaded {len(cell_names)} cell barcodes")
        except Exception as e:
            print(f"Warning: Could not load column annotations: {e}")
    
    return matrix, gene_names, cell_names


def read_excel_columns(
    file_path: str,
    columns: Optional[List[str]] = None,
    sheet_name: Optional[Union[str, int]] = 0,
    dtype: Optional[Dict[str, str]] = None,
    engine: Optional[str] = None,
    use_progress: bool = True
) -> pd.DataFrame:
    """
    Read a Microsoft Excel file and return a pandas DataFrame with selected columns.

    Parameters:
    -----------
    file_path : str
        Path to the Excel file (.xlsx, .xls, .xlsm).
    columns : list of str, optional
        Column names to select. If None, returns all columns.
    sheet_name : str or int, optional (default=0)
        Name or index of the sheet to read.
    dtype : dict, optional
        Optional mapping of column name -> dtype to enforce after reading.
    engine : str, optional
        Excel engine to use (e.g., 'openpyxl'). If None, pandas will auto-select.
    use_progress : bool, optional (default=True)
        If True, show a progress bar while reading using openpyxl+tqdm when available.

    Returns:
    --------
    df : pandas.DataFrame
        DataFrame containing the selected columns.
    """
    # If requested, try progress-enabled path using openpyxl + tqdm
    if use_progress:
        try:
            import openpyxl  # type: ignore
            try:
                from tqdm import tqdm  # type: ignore
            except Exception:
                tqdm = None  # fallback without visible progress

            wb = openpyxl.load_workbook(file_path, read_only=True, data_only=True)
            ws = wb.worksheets[sheet_name] if isinstance(sheet_name, int) else wb[sheet_name]

            header_row = next(ws.iter_rows(min_row=1, max_row=1, values_only=True))
            headers = list(header_row)

            # Determine which columns to keep
            if columns is not None:
                name_to_idx = {name: idx for idx, name in enumerate(headers)}
                missing = [c for c in columns if c not in name_to_idx]
                if missing:
                    print(f"Warning: Missing columns not found in Excel: {missing}")
                selected = [c for c in columns if c in name_to_idx]
                selected_indices = [name_to_idx[c] for c in selected]
                out_headers = selected
            else:
                selected_indices = list(range(len(headers)))
                out_headers = headers

            total_rows = max(ws.max_row - 1, 0)
            progress_iter = tqdm(total=total_rows, desc="Reading Excel", unit="row") if tqdm else None

            records = []
            # Iterate data rows
            for row in ws.iter_rows(min_row=2, values_only=True):
                if row is None:
                    if progress_iter:
                        progress_iter.update(1)
                    continue
                # Ensure row has enough length
                row_values = list(row)
                if len(row_values) < len(headers):
                    row_values += [None] * (len(headers) - len(row_values))
                # Select indices
                record = [row_values[i] for i in selected_indices]
                records.append(record)
                if progress_iter:
                    progress_iter.update(1)

            if progress_iter:
                progress_iter.close()
            wb.close()

            df = pd.DataFrame.from_records(records, columns=out_headers)
        except Exception as e:
            print(f"Note: Falling back to pandas.read_excel (reason: {e})")
            df = pd.read_excel(file_path, sheet_name=sheet_name, engine=engine)
    else:
        # Read the sheet via pandas directly
        df = pd.read_excel(file_path, sheet_name=sheet_name, engine=engine)

    # Select requested columns if provided
    if columns is not None:
        missing = [c for c in columns if c not in df.columns]
        if missing:
            print(f"Warning: Missing columns not found in Excel: {missing}")
        present = [c for c in columns if c in df.columns]
        df = df[present]

    # Apply dtype conversions if requested
    if dtype:
        for col, target_dtype in dtype.items():
            if col in df.columns:
                try:
                    df[col] = df[col].astype(target_dtype)
                except Exception as e:
                    print(f"Warning: Could not cast column '{col}' to {target_dtype}: {e}")

    return df


def create_anndata_object(matrix, gene_names=None, cell_names=None, obs=None, transpose=False):
    """
    Create an AnnData object (standard format for single-cell analysis).
    Requires scanpy/anndata package.
    
    Parameters:
    -----------
    matrix : scipy.sparse matrix
        Expression matrix (assumed to be genes × cells if transpose=False)
    gene_names : array-like, optional
        Gene names
    cell_names : array-like, optional
        Cell barcodes
    obs : pandas.DataFrame, optional
        Cell metadata (must have same number of rows as cells)
    transpose : bool, default=False
        If True, transpose the matrix before creating AnnData
        (AnnData expects cells as rows, genes as columns)
    
    Returns:
    --------
    adata : anndata.AnnData
        AnnData object ready for analysis
    """
    try:
        import anndata as ad
        
        # Create AnnData object
        # AnnData expects cells as rows (obs), genes as columns (var)
        if transpose:
            adata = ad.AnnData(matrix.T)
        else:
            # Matrix is already in cells × genes format
            adata = ad.AnnData(matrix)
        
        if gene_names is not None:
            adata.var_names = gene_names
        if cell_names is not None:
            adata.obs_names = cell_names
        if obs is not None:
            adata.obs = obs
        
        return adata
    
    except ImportError:
        print("anndata package not installed. Install with: pip install anndata scanpy")
        return None

def filter_anndata_object(adata, min_genes=200, min_cells=3, min_counts=None, 
                          max_counts=None, mito_max=None):
    """
    Filter an AnnData object based on the number of genes and cells.
    
    Parameters:
    -----------
    adata : anndata.AnnData
        AnnData object to filter
    min_genes : int, default=200
        Minimum number of genes detected per cell
    min_cells : int, default=3
        Minimum number of cells expressing each gene
    min_counts : int, optional
        Minimum total counts (UMIs) per cell
    max_counts : int, optional
        Maximum total counts (UMIs) per cell    
    Returns:
    --------
    adata : anndata.AnnData
        Filtered AnnData object
    """
    import scanpy as sc
    sc.pp.filter_cells(adata, min_genes=min_genes, min_counts=min_counts, max_counts=max_counts)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    if mito_max is not None:
        mask = (adata.obs['percent.mito'] < mito_max).fillna(False)
        adata = adata[mask].copy()
    else:
        print("Warning: mito_max is not set. Skipping mito filter.")

    return adata


def get_region_file_paths(region, data_dir="data", base_prefix="2025-10-22_Astrocytes_{region}"):
    """
    Get file paths for region-specific matrix, row, and column annotation files.
    
    Parameters:
    -----------
    region : str
        Region name (e.g., 'EC', 'ITG', 'PFC', 'V1', 'V2')
    data_dir : str or Path, default="data"
        Base directory containing region files
    base_prefix : str, default="2025-10-22_Astrocytes_{region}"
        File naming pattern with {region} placeholder
    
    Returns:
    --------
    mtx_path : Path
        Path to matrix.mtx file
    row_path : Path
        Path to row_annotation.txt file
    col_path : Path
        Path to cell_annotation.txt file
    """
    data_dir = Path(data_dir)
    prefix = base_prefix.format(region=region)
    mtx_path = (data_dir / f"{prefix}_matrix.mtx").resolve()
    row_path = (data_dir / f"{prefix}_row_annotation.txt").resolve()
    col_path = (data_dir / f"{prefix}_cell_annotation.txt").resolve()
    return mtx_path, row_path, col_path


def get_dge_results_dir(version: int, base_dir: str = None) -> Path:
    """
    Get the path to DGE results directory for a given version.
    
    Parameters:
    -----------
    version : int
        DGE version (1, 2, 3, or 4)
    base_dir : str, optional
        Base directory for results. If None, uses current working directory.
    
    Returns:
    --------
    results_dir : Path
        Path to the DGE results directory (e.g., results/dge4)
    """
    if base_dir is None:
        from pathlib import Path
        base_dir = Path.cwd()
    else:
        base_dir = Path(base_dir)
    
    return base_dir / "results" / f"dge{version}"


def filter_cells(matrix, cell_names, metadata, mito_max=15.0):
    """
    Filter cells based on mitochondrial percentage threshold.
    
    Parameters:
    -----------
    matrix : scipy.sparse matrix
        Expression matrix (genes × cells)
    cell_names : list
        List of cell names/barcodes
    metadata : pandas.DataFrame
        Metadata DataFrame with 'cell_annotation' and 'percent.mito' columns
    mito_max : float, default=15.0
        Maximum mitochondrial percentage threshold
    
    Returns:
    --------
    filtered_matrix : scipy.sparse matrix
        Filtered expression matrix
    filtered_cell_names : list
        List of cell names that passed filters
    filtered_metadata : pandas.DataFrame
        Filtered metadata DataFrame
    """
    # Ensure cell_names matches matrix columns
    n_cols = matrix.shape[1]
    if len(cell_names) != n_cols:
        print(f"  Warning: cell_names length ({len(cell_names)}) != matrix columns ({n_cols})")
        cell_names = cell_names[:n_cols]  # Truncate to match matrix
    
    filtered_metadata = metadata[metadata['cell_annotation'].isin(cell_names)].copy()
    filtered_metadata = filtered_metadata.set_index('cell_annotation').reindex(cell_names)
    n_before_filter = len(filtered_metadata)

    # Debug: Check percent.mito values
    mito_values = filtered_metadata['percent.mito']
    print(f"  Mitochondrial percentage stats:")
    print(f"    Min: {mito_values.min():.2f}, Max: {mito_values.max():.2f}, Mean: {mito_values.mean():.2f}")
    print(f"    NaN count: {mito_values.isna().sum()}")
    print(f"    Threshold: {mito_max}")

    # Filter by mitochondrial percentage only
    mito_mask = (filtered_metadata['percent.mito'] < mito_max).fillna(False)
    combined_mask = mito_mask
    filtered_metadata = filtered_metadata[combined_mask].copy()
    
    print(f"  Before: {n_before_filter:,} cells")
    print(f"  After: {len(filtered_metadata):,} cells")
    print(f"  Filtered out: {n_before_filter - len(filtered_metadata):,} cells")
    
    # Get cell indices to keep (cells that passed filters)
    cells_to_keep = set(filtered_metadata.index.dropna().tolist())
    cell_indices_to_keep = [i for i, name in enumerate(cell_names) if name in cells_to_keep]
    
    # Ensure indices are within matrix bounds
    cell_indices_to_keep = [i for i in cell_indices_to_keep if 0 <= i < n_cols]
    
    # Filter matrix and cell names
    filtered_matrix = matrix[:, cell_indices_to_keep]
    filtered_cell_names = [cell_names[i] for i in cell_indices_to_keep]
    
    return filtered_matrix, filtered_cell_names, filtered_metadata