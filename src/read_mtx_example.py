"""
Example script to read .mtx files (Matrix Market format) in Python.

This script demonstrates how to load sparse matrices from .mtx files,
commonly used in single-cell RNA sequencing data.

The typical 10X Genomics format includes:
- matrix.mtx: The sparse count matrix
- genes.tsv or features.tsv: Gene/feature names (rows)
- barcodes.tsv: Cell barcodes (columns)
"""

import scipy.io
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
import random


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
    print(f"Matrix format: {matrix.format}")
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


def convert_to_dense(matrix, subset_genes=None, subset_cells=None):
    """
    Convert sparse matrix to dense format (use with caution - can be memory intensive).
    
    Parameters:
    -----------
    matrix : scipy.sparse matrix
        Sparse matrix to convert
    subset_genes : array-like, optional
        Indices of genes to subset
    subset_cells : array-like, optional
        Indices of cells to subset
    
    Returns:
    --------
    dense_matrix : numpy.ndarray
        Dense matrix representation
    """
    if subset_genes is not None or subset_cells is not None:
        if subset_genes is not None:
            matrix = matrix[subset_genes, :]
        if subset_cells is not None:
            matrix = matrix[:, subset_cells]
    
    return matrix.toarray()


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


# Example usage
if __name__ == "__main__":
    # Paths to your files
    mtx_path = r"i:\sf2026\data\2025-10-22_Astrocytes_EC_matrix.mtx"
    row_annotation_path = r"i:\sf2026\data\2025-10-22_Astrocytes_EC_row_annotation.txt"
    col_annotation_path = r"i:\sf2026\data\2025-10-22_Astrocytes_EC_cell_annotation.txt"
    metadata_path = r"i:\sf2026\data\2025-10-18_Astrocytes_metadata.csv"
    metadata = pd.read_csv(metadata_path, low_memory=False)

    # Print all column names and values from the first row (index 0) of the metadata DataFrame
    print("\n" + "=" * 60)
    print("First Row of Metadata:")
    print("=" * 60)
    first_row = metadata.iloc[0]
    for col, val in first_row.items():
        print(f"{col}: {val}")
    
    print("=" * 60)
    print("Reading MTX file with annotations")
    print("=" * 60)
    
    # Read the MTX file
    matrix, gene_names, cell_names = read_mtx_file(
        mtx_path=mtx_path,
        row_annotation_path=row_annotation_path,
        col_annotation_path=col_annotation_path,
        transpose=False  # Typically need to transpose for scRNA-seq
    )
    
    print("\n" + "=" * 60)
    print("Matrix Summary")
    print("=" * 60)
    print(f"Final shape: {matrix.shape}")
    print(f"Rows (genes): {matrix.shape[0]:,}")
    print(f"Columns (cells): {matrix.shape[1]:,}")
    
    if gene_names:
        print(f"\nFirst 5 genes: {gene_names[:5]}")
    if cell_names:
        print(f"First 5 cells: {cell_names[:5]}")
    
    # Example: Access specific gene or cell
    if gene_names and cell_names:
        # Get expression of first gene across all cells
        for i in range(5):
            random_gene_index = random.randint(0, len(gene_names) - 1)  
            random_gene_expression = matrix[random_gene_index, :].toarray().flatten()
            print(f"Expression of '{gene_names[random_gene_index]}' across cells:")
            print(f"  Mean: {np.mean(random_gene_expression):.2f}")
            print(f"  Max: {np.max(random_gene_expression)}")
            print(f"  Non-zero cells: {np.sum(random_gene_expression > 0):,}")

        # Get expression profile of first cell
        first_cell_expression = matrix[:, 0].toarray().flatten()
        print(f"\nExpression profile of '{cell_names[0]}':")
        print(f"  Total UMIs: {np.sum(first_cell_expression):,}")
        print(f"  Genes detected: {np.sum(first_cell_expression > 0):,}")
    
    # Example: Convert to pandas DataFrame (WARNING: can be memory intensive!)
    # Only do this for small subsets!
    print("\n" + "=" * 60)
    print("Converting small subset to DataFrame (for demonstration)")
    print("=" * 60)
    
    if gene_names and cell_names:
        # Just first 100 genes and 10 cells as example
        subset_matrix = matrix[:100, :10].toarray()
        df = pd.DataFrame(
            subset_matrix,
            index=gene_names[:100],
            columns=cell_names[:10]
        )
        print(f"DataFrame shape: {df.shape}")
        print("\nFirst few rows and columns:")
        print(df.iloc[:5, :5])
    
    # Example: Create AnnData object (if scanpy/anndata is installed)
    print("\n" + "=" * 60)
    print("Creating AnnData object (optional)")
    print("=" * 60)
    
    adata = create_anndata_object(
        matrix=matrix,
        gene_names=gene_names,
        cell_names=cell_names,
        transpose=True  # Matrix is genes × cells, need to transpose for AnnData
    )
    
    if adata is not None:
        print(f"AnnData object created successfully!")
        print(f"  Shape: {adata.shape}")
        print(f"  Variables (genes): {adata.n_vars:,}")
        print(f"  Observations (cells): {adata.n_obs:,}")

