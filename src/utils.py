"""
Utility functions for reading and processing single-cell data.

This module provides helper functions for:
- Reading .mtx files (Matrix Market format)
- Creating AnnData objects from matrices
- Reading Excel files into pandas DataFrames
"""

import scipy.io
from scipy.sparse import csr_matrix
import pandas as pd
from typing import Optional, List, Dict, Union


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

