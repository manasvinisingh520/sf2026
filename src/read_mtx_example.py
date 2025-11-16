"""
Example script to read .mtx files (Matrix Market format) in Python.

This script demonstrates how to load sparse matrices from .mtx files,
commonly used in single-cell RNA sequencing data.

The typical 10X Genomics format includes:
- matrix.mtx: The sparse count matrix
- genes.tsv or features.tsv: Gene/feature names (rows)
- barcodes.tsv: Cell barcodes (columns)
"""

import pandas as pd
import numpy as np
import random

from utils import read_mtx_file, create_anndata_object, read_excel_columns


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


# Example usage
if __name__ == "__main__":
    # Paths to your files
    mtx_path = r"i:\sf2026\data\2025-10-22_Astrocytes_EC_matrix.mtx"
    row_annotation_path = r"i:\sf2026\data\2025-10-22_Astrocytes_EC_row_annotation.txt" # genes.tsv
    col_annotation_path = r"i:\sf2026\data\2025-10-22_Astrocytes_EC_cell_annotation.txt" # cells.tsv
    metadata_path = r"i:\sf2026\data\2025-11-16_Astrocytes_metadata.xlsx" # metadata.xlsx
    #metadata = pd.read_csv(metadata_path, low_memory=False)
    metadata = read_excel_columns(metadata_path, columns=['cell_annotation', "RIN", "Path..Group."])

    # Print all column names and values from the first row (index 0) of the metadata DataFrame
    print(f"Number of rows in metadata: {len(metadata.index)}")
    print("\n" + "=" * 60)
    print("First Row of Metadata:")
    print("=" * 60)
    first_row = metadata.iloc[0]
    for col, val in first_row.items():
        #print(f"{col}: {val}")
        pass
    
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
    
    # Filter metadata to include only rows where 'cell_annotation' matches the loaded cell_names
    if cell_names and 'cell_annotation' in metadata.columns:
        filtered_metadata = metadata[metadata['cell_annotation'].isin(cell_names)].copy()
        print(f"Filtered metadata shape: {filtered_metadata.shape}")
    else:
        filtered_metadata = metadata
        print("Note: 'cell_annotation' column not found in metadata or cell_names unavailable.")
    # Example: Access specific gene or cell
    print (f"Number of rows in filtered metadata: {len(filtered_metadata.index)}")

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
        transpose=True,  # Matrix is genes × cells, need to transpose for AnnData
        obs=filtered_metadata,
    )
    
    if adata is not None:
        print(f"AnnData object created successfully!")
        print(f"  Shape: {adata.shape}")
        print(f"  Variables (genes): {adata.n_vars:,}")
        print(f"  Observations (cells): {adata.n_obs:,}")

