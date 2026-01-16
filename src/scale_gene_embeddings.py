"""
Script to scale gene embeddings according to expression levels in specific cell(s).

Takes a gene embedding file (TSV format) and scales each gene's embedding vector
by its expression level in the specified cell ID(s).
"""

import pandas as pd
import numpy as np
import argparse
import os
import anndata as ad
from tqdm import tqdm


def load_gene_embeddings(embeddings_path):
    """
    Load gene embeddings from TSV file.
    
    Parameters:
    -----------
    embeddings_path : str
        Path to gene embeddings TSV file
        
    Returns:
    --------
    embeddings_df : pd.DataFrame
        DataFrame with gene IDs as index and embedding dimensions as columns
    """
    print(f"Loading gene embeddings from: {embeddings_path}")
    
    # Read TSV file
    df = pd.read_csv(embeddings_path, sep='\t', header=0)
    
    # First column should be gene IDs, rest are embedding dimensions
    df = df.iloc[1:].copy()
    
    # Set first column as index (gene IDs)
    gene_ids = df.iloc[:, 0].values
    embedding_data = df.iloc[:, 1:].values.astype(float)
    
    # Create DataFrame with gene IDs as index
    embeddings_df = pd.DataFrame(
        embedding_data,
        index=gene_ids,
        columns=range(embedding_data.shape[1])
    )
    
    print(f"Loaded {len(embeddings_df)} genes with {embedding_data.shape[1]} embedding dimensions")
    
    return embeddings_df


def get_expression_levels_for_cell(adata, cell_id):
    """
    Get expression levels for genes in a specific cell.
    
    Parameters:
    -----------
    adata : anndata.AnnData
        AnnData object with expression data
    cell_id : str
        Cell ID to get expression from
        
    Returns:
    --------
    expression_dict : dict
        Dictionary mapping gene ID to expression level
    """
    # Find cell index
    if cell_id not in adata.obs_names:
        raise ValueError(f"Cell ID '{cell_id}' not found in AnnData. Available cells: {len(adata.obs_names)}")
    
    cell_index = adata.obs_names.get_loc(cell_id)
    
    # Get expression for this cell (1 × genes)
    # Note: adata.X contains raw UMI counts
    expression_subset = adata.X[cell_index:cell_index+1, :]
    
    # Convert to dense if sparse
    if hasattr(expression_subset, 'toarray'):
        expression_subset = expression_subset.toarray()
    
    # Apply CPM normalization (Counts Per Million) to account for sequencing depth
    # Calculate total counts for this cell
    total_counts = np.sum(expression_subset)
    if total_counts == 0:
        total_counts = 1  # Avoid division by zero
    
    # Normalize to counts per million
    cpm_subset = (expression_subset / total_counts) * 1e6
    
    # Apply log1p transformation (log(CPM + 1))
    log1p_subset = np.log1p(cpm_subset)
    
    # Flatten to 1D array
    expression_values = log1p_subset.flatten()
    
    # Create dictionary mapping gene IDs to expression
    gene_ids = adata.var_names.values
    expression_dict = dict(zip(gene_ids, expression_values))
    
    print(f"Retrieved expression for {len(expression_dict)} genes")
    print(f"  Expression range: {np.min(expression_values):.2f} - {np.max(expression_values):.2f}")
    print(f"  Mean expression: {np.mean(expression_values):.2f}")
    
    return expression_dict


def scale_embeddings(embeddings_df, expression_dict):
    """
    Scale gene embeddings by expression levels using weight method (embedding * sqrt(log1p(CPM))).
    
    Parameters:
    -----------
    embeddings_df : pd.DataFrame
        DataFrame with gene embeddings (genes × dimensions)
    expression_dict : dict
        Dictionary mapping gene ID to expression level (log1p(CPM))
        
    Returns:
    --------
    scaled_embeddings_df : pd.DataFrame
        Scaled embeddings DataFrame
    """
    print("Scaling embeddings using weight method (embedding * sqrt(log1p(CPM)))")
    print("Using log1p(CPM) expression values directly (no additional normalization)")
    
    # Get expression values for genes in embeddings
    expression_values = []
    
    for gene_id in embeddings_df.index:
        if gene_id in expression_dict:
            expression_values.append(expression_dict[gene_id])
        else:
            # Gene not in expression data, use expression = 0
            expression_values.append(0.0)
    
    expression_array = np.array(expression_values)
    
    # Apply scaling using weight method: embedding * sqrt(expression)
    scaled_embeddings = embeddings_df.values.copy()
    
    # Scale by square root of expression (less aggressive than direct multiplication)
    for i, expr in enumerate(expression_array):
        scaled_embeddings[i, :] *= np.sqrt(expr + 1e-8)  # Add small epsilon to avoid sqrt(0)
    
    # Create new DataFrame
    scaled_embeddings_df = pd.DataFrame(
        scaled_embeddings,
        index=embeddings_df.index,
        columns=embeddings_df.columns
    )
    
    # Count genes with non-zero expression
    n_expressed = np.sum(expression_array > 0)
    print(f"Scaled {len(scaled_embeddings_df)} gene embeddings")
    print(f"  {n_expressed} genes had non-zero expression")
    print(f"  {len(scaled_embeddings_df) - n_expressed} genes had zero expression (embeddings set to zero)")
    
    return scaled_embeddings_df


def save_scaled_embeddings(scaled_embeddings_df, output_path):
    """
    Save scaled embeddings to TSV file in the same format as input.
    
    Parameters:
    -----------
    scaled_embeddings_df : pd.DataFrame
        Scaled embeddings DataFrame
    output_path : str
        Output file path
    """
    print(f"Saving scaled embeddings to: {output_path}")
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
    
    # Prepare data for writing (same format as input: gene_id, dim0, dim1, ...)
    with open(output_path, 'w') as f:
        # Write header row with column indices
        n_dims = len(scaled_embeddings_df.columns)
        header = ',' + ','.join([str(i) for i in range(n_dims)])
        f.write(header + '\n')
        
        # Write each gene's embedding
        for gene_id in tqdm(scaled_embeddings_df.index, desc="Writing embeddings", unit="gene"):
            embedding_values = scaled_embeddings_df.loc[gene_id].values
            row = f"{gene_id}," + ','.join([str(val) for val in embedding_values])
            f.write(row + '\n')
    
    print(f"Successfully saved {len(scaled_embeddings_df)} scaled embeddings")


def main():
    parser = argparse.ArgumentParser(
        description='Scale gene embeddings by expression levels in specific cell(s)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '-path',
        type=str,
        required=True,
        help='Path to gene embeddings TSV file'
    )
    parser.add_argument(
        '-anndata',
        type=str,
        required=True,
        help='Path to AnnData file (.h5ad) containing expression data'
    )
    parser.add_argument(
        '-cell_ids',
        type=str,
        nargs='+',
        required=True,
        help='Cell ID(s) to get expression from (can specify multiple)'
    )
    parser.add_argument(
        '-output',
        type=str,
        default=None,
        help='Output file path pattern (default: adds "_scaled_{cell_id}" before extension for each cell)'
    )
    
    args = parser.parse_args()
    
    print("="*60)
    print("Scaling Gene Embeddings by Expression Levels")
    print("="*60)
    
    # Load gene embeddings
    embeddings_df = load_gene_embeddings(args.path)
    
    # Load AnnData
    print(f"\nLoading AnnData from: {args.anndata}")
    adata = ad.read_h5ad(args.anndata)
    print(f"AnnData shape: {adata.shape} (cells × genes)")
    
    # Process each cell separately
    base, ext = os.path.splitext(args.path)
    
    for cell_id in args.cell_ids:
        print("\n" + "="*60)
        print(f"Processing cell: {cell_id}")
        print("="*60)
        
        # Get expression levels for this cell (log1p(CPM))
        expression_dict = get_expression_levels_for_cell(adata, cell_id)
        
        # Scale embeddings
        scaled_embeddings_df = scale_embeddings(
            embeddings_df,
            expression_dict
        )
        
        # Determine output path
        if args.output is None:
            # Create output path with cell ID
            output_path = f"{base}_scaled_{cell_id}{ext}"
        else:
            # Use provided pattern, insert cell_id before extension
            output_base, output_ext = os.path.splitext(args.output)
            output_path = f"{output_base}_{cell_id}{output_ext}"
        
        # Save scaled embeddings
        save_scaled_embeddings(scaled_embeddings_df, output_path)
    
    print("\n" + "="*60)
    print(f"Complete! Created {len(args.cell_ids)} scaled embedding file(s)")
    print("="*60)


if __name__ == "__main__":
    main()
