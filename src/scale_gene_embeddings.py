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
import matplotlib.pyplot as plt
import umap
HAS_UMAP = True


def load_gene_embeddings(embeddings_path):
    """
    Load gene embeddings from TSV/CSV file.
    
    Format: First row is header (starts with comma, then dimension indices: ,0,1,2,3,...)
    Subsequent rows: gene ID in first column, embedding values in remaining columns.
    File is comma-separated.
    
    Parameters:
    -----------
    embeddings_path : str
        Path to gene embeddings file
        
    Returns:
    --------
    embeddings_df : pd.DataFrame
        DataFrame with gene IDs as index and embedding dimensions as columns
    """
    print(f"Loading gene embeddings from: {embeddings_path}")
    
    # Read CSV file with header (first row contains dimension indices)
    # First column will be unnamed/empty in header, then 0,1,2,3... for dimensions
    try:
        df = pd.read_csv(embeddings_path, sep=',', header=0)
    except:
        # Try tab-separated as fallback
        df = pd.read_csv(embeddings_path, sep='\t', header=0)
    
    # First column contains gene IDs (may be named 'Unnamed: 0' or similar)
    # Get the first column name (should be the gene ID column)
    first_col_name = df.columns[0]
    gene_ids = df[first_col_name].astype(str).values
    
    # Remaining columns are embedding dimensions
    embedding_cols = [col for col in df.columns if col != first_col_name]
    embedding_data = df[embedding_cols].values.astype(float)
    
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
    
    # Debug: Check gene ID formats
    print(f"\nDebug: Checking gene ID matching...")
    print(f"  Sample embedding gene IDs (first 5): {list(embeddings_df.index[:5])}")
    print(f"  Sample AnnData gene IDs (first 5): {list(expression_dict.keys())[:5]}")
    
    # Check how many genes match
    matching_genes = [gid for gid in embeddings_df.index if gid in expression_dict]
    print(f"  Genes in embeddings: {len(embeddings_df)}")
    print(f"  Genes in AnnData: {len(expression_dict)}")
    print(f"  Matching genes: {len(matching_genes)}")

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


def plot_before_after(embeddings_df, scaled_embeddings_df, expression_dict, cell_id, output_dir=None):
    """
    Plot before and after visualization of gene embeddings using UMAP or PCA.
    
    Parameters:
    -----------
    embeddings_df : pd.DataFrame
        Original embeddings DataFrame
    scaled_embeddings_df : pd.DataFrame
        Scaled embeddings DataFrame
    expression_dict : dict
        Dictionary mapping gene ID to expression level
    cell_id : str
        Cell ID (for title)
    output_dir : str, optional
        Directory to save plot (if None, displays plot)
    """
    print("\nGenerating before/after visualization...")
    
    # Get expression values for coloring
    expression_values = []
    for gene_id in embeddings_df.index:
        if gene_id in expression_dict:
            expression_values.append(expression_dict[gene_id])
        else:
            expression_values.append(0.0)
    expression_array = np.array(expression_values)
    
    # Reduce dimensionality for visualization
    original_embeddings = embeddings_df.values
    scaled_embeddings = scaled_embeddings_df.values
    
    print("  Computing UMAP for visualization...")
    reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.1)
    original_2d = reducer.fit_transform(original_embeddings)
    scaled_2d = reducer.fit_transform(scaled_embeddings)
    method_name = "UMAP"
    
    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # Plot original embeddings
    ax1 = axes[0]
    scatter1 = ax1.scatter(
        original_2d[:, 0], original_2d[:, 1],
        c=expression_array, cmap='viridis', alpha=0.6, s=20
    )
    ax1.set_xlabel(f'{method_name} Dimension 1', fontsize=12)
    ax1.set_ylabel(f'{method_name} Dimension 2', fontsize=12)
    ax1.set_title(f'Before Scaling\n(Cell: {cell_id})', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    plt.colorbar(scatter1, ax=ax1, label='log1p(CPM) Expression')
    
    # Plot scaled embeddings
    ax2 = axes[1]
    scatter2 = ax2.scatter(
        scaled_2d[:, 0], scaled_2d[:, 1],
        c=expression_array, cmap='viridis', alpha=0.6, s=20
    )
    ax2.set_xlabel(f'{method_name} Dimension 1', fontsize=12)
    ax2.set_ylabel(f'{method_name} Dimension 2', fontsize=12)
    ax2.set_title(f'After Scaling (×√expression)\n(Cell: {cell_id})', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    plt.colorbar(scatter2, ax=ax2, label='log1p(CPM) Expression')
    
    plt.tight_layout()
    plt.show()
    
    # Print some statistics
    print(f"\n  Visualization statistics:")
    print(f"    Original embedding range: [{np.min(original_embeddings):.4f}, {np.max(original_embeddings):.4f}]")
    print(f"    Scaled embedding range: [{np.min(scaled_embeddings):.4f}, {np.max(scaled_embeddings):.4f}]")
    print(f"    Expression range: [{np.min(expression_array):.4f}, {np.max(expression_array):.4f}]")
    print(f"    Genes with expression > 0: {np.sum(expression_array > 0)} / {len(expression_array)}")


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
        required=False,
        default='data\model_data\EC_AnnData_perCell.h5ad',
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
    parser.add_argument(
        '-plot',
        action='store_true',
        help='Generate before/after visualization plots'
    )
    parser.add_argument(
        '-plot_dir',
        type=str,
        default=None,
        help='Directory to save plots (if not provided and -plot is set, displays plots)'
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
    
    # Print example cell IDs
    print(f"\nAvailable cell IDs (showing first 10):")
    for i, cell_id in enumerate(adata.obs_names[:10]):
        print(f"  {i+1}. {cell_id}")
    if len(adata.obs_names) > 10:
        print(f"  ... and {len(adata.obs_names) - 10} more cells")
    print(f"Total cells: {len(adata.obs_names)}")
    
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
        
        # Generate visualization if requested
        if args.plot:
            plot_before_after(
                embeddings_df,
                scaled_embeddings_df,
                expression_dict,
                cell_id,
                output_dir=args.plot_dir
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
        #save_scaled_embeddings(scaled_embeddings_df, output_path)
    
    print("\n" + "="*60)
    print(f"Complete! Created {len(args.cell_ids)} scaled embedding file(s)")
    print("="*60)


if __name__ == "__main__":
    main()
