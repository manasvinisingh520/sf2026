"""
Script to create gene-by-cell expression matrices for each region.

For each region, this script:
1. Loads the expression data (MTX files)
2. Creates a matrix with genes on y-axis and cells on x-axis
3. Orders cells by condition (1, 2, 3, 4)
4. Optionally visualizes or saves the matrix
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.sparse import csr_matrix
import argparse
from utils import (
    read_mtx_file, 
    create_anndata_object, 
    read_excel_columns, 
    get_region_file_paths,
    filter_anndata_object,
    min_max_normalize,
    filter_cells,
    get_DEGs
)
from config import REGIONS, REGION_TO_TAB, DATA_DIR, METADATA_DATE_PREFIX


def load_region_data(region, cell_type="Astrocytes", base_dir=DATA_DIR, date_prefix=None):
    """Load expression data for a given region."""
    # Set date prefix based on cell type
    if date_prefix is None:
        if cell_type == "Astrocytes":
            date_prefix = "2025-11-16"
        elif cell_type == "Microglia":
            date_prefix = "2025-12-27"
        else:
            raise ValueError(f"Cell type '{cell_type}' not supported. Available cell types: Astrocytes, Microglia")
    
    # Get file paths using utility function
    if cell_type == "Astrocytes":
        base_prefix = f"{date_prefix}_Astrocytes_{region}"
    elif cell_type == "Microglia":
        base_prefix = f"{date_prefix}_Microglia_{region}"
    else:
        raise ValueError(f"Cell type '{cell_type}' not supported. Available cell types: Astrocytes, Microglia")
    
    mtx_path, row_annotation_path, col_annotation_path = get_region_file_paths(
        region, 
        data_dir=base_dir, 
        base_prefix=base_prefix
    )
    
    # Get metadata path based on cell type
    if cell_type == "Astrocytes":
        metadata_path = os.path.join(base_dir, f"{date_prefix}_Astrocytes_metadata.xlsx")
    elif cell_type == "Microglia":
        metadata_path = os.path.join(base_dir, f"{base_prefix}_metadata.csv")
    
    print(f"\nLoading data for region: {region}, cell type: {cell_type}")
    
    # Load metadata based on cell type
    if cell_type == "Astrocytes":
        # Get tab index for this region
        tab_index = REGION_TO_TAB.get(region)
        if tab_index is None:
            raise ValueError(f"Region '{region}' not found in REGION_TO_TAB. Available regions: {list(REGION_TO_TAB.keys())}")
        
        # Load metadata (include RIN and Total.Genes.Detected for filtering)
        metadata = read_excel_columns(
            metadata_path,
            columns=['cell_annotation', 'Path..Group.', 'SampleName', 'percent.mito', 'RIN', 'Total.Genes.Detected'],
            sheet_name=tab_index
        )
    else:  # Microglia
        metadata = pd.read_csv(metadata_path)
        # Rename columns with spaces to use underscores to avoid formula syntax errors
        rename_dict = {}
        if 'Pathology Stage' in metadata.columns:
            rename_dict['Pathology Stage'] = 'Pathology_Stage'
        if 'Donor ID' in metadata.columns:
            rename_dict['Donor ID'] = 'Donor_ID'
        if rename_dict:
            metadata = metadata.rename(columns=rename_dict)
        # Convert condition column to strings
        if 'Pathology_Stage' in metadata.columns:
            metadata['Pathology_Stage'] = metadata['Pathology_Stage'].astype(str)
    
    print(f"Loaded metadata: {metadata.shape}")
    
    # Read the MTX file
    matrix, gene_names, cell_names = read_mtx_file(
        mtx_path=str(mtx_path),
        row_annotation_path=str(row_annotation_path),
        col_annotation_path=str(col_annotation_path),
        transpose=False  # Matrix will be genes × cells
    )
    print(f"Loaded matrix shape: {matrix.shape} (genes × cells)")
    
    # Filter and align metadata to match cell_names
    if cell_names and 'cell_annotation' in metadata.columns:
        # Astrocytes case: use cell_annotation column
        filtered_metadata = metadata[metadata['cell_annotation'].isin(cell_names)].copy()
        filtered_metadata = filtered_metadata.set_index('cell_annotation')
        filtered_metadata = filtered_metadata.reindex(cell_names)
    elif cell_names:
        # Microglia case: cell names are typically in the first column (e.g., 'Unnamed: 0')
        first_col = metadata.columns[0]
        filtered_metadata = metadata.set_index(first_col)
        filtered_metadata = filtered_metadata.reindex(cell_names)
    else:
        # No cell_names available
        filtered_metadata = pd.DataFrame(index=range(matrix.shape[1]))
    
    # Create AnnData object
    # The filtered_metadata index should already match cell_names
    adata = create_anndata_object(
        matrix=matrix,
        gene_names=gene_names,
        cell_names=cell_names,
        transpose=True,  # Matrix is genes × cells, transpose for AnnData (cells × genes)
        obs=filtered_metadata,
    )
    
    if adata is None:
        raise ImportError("Failed to create AnnData object. Make sure anndata is installed.")
    
    print(f"AnnData object created: {adata.shape} (cells × genes)")
    
    # Apply basic quality filters
    adata = filter_anndata_object(adata, min_genes=5800, min_cells=0, min_counts=None, max_counts=None, mito_max=(0.15 if cell_type == "Astrocytes" else None))
    
    print(f"After filtering: {adata.shape} (cells × genes)")
    
    return adata


def order_cells_by_condition(adata, condition_col='Path..Group.', cell_type="Astrocytes"):
    """
    Order cells by condition (1, 2, 3, 4).
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with cells and metadata
    condition_col : str
        Name of the condition column
    cell_type : str
        Cell type (Astrocytes or Microglia) - used to determine condition column if not provided
    
    Returns:
    --------
    adata_ordered : AnnData
        AnnData object with cells ordered by condition
    condition_order : list
        List of conditions in order
    """
    
    if condition_col not in adata.obs.columns:
        raise ValueError(f"Condition column '{condition_col}' not found in adata.obs.")
    
    # Convert condition to string and handle numeric values
    conditions = adata.obs[condition_col].astype(str)
    
    # Define the desired order (1, 2, 3, 4)
    desired_order = ['1', '2', '3', '4']
    
    # Get unique conditions and sort according to desired order
    unique_conditions = conditions.unique()
    condition_order = [c for c in desired_order if c in unique_conditions]
    
    print(f"Condition counts:")
    for cond in condition_order:
        count = (conditions == cond).sum()
        print(f"  Condition {cond}: {count} cells")
    
    # Create a categorical column with ordered categories
    conditions_cat = pd.Categorical(conditions, categories=condition_order, ordered=True)
    
    # Sort by condition
    sort_indices = np.argsort(conditions_cat.codes)
    adata_ordered = adata[sort_indices].copy()
    
    return adata_ordered, condition_order


def create_expression_matrix(adata, genes=None):
    """
    Create a gene-by-cell expression matrix.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object (cells × genes)
    genes : list, optional
        List of gene names to include. If None, includes all genes.
    
    Returns:
    --------
    matrix : np.ndarray or csr_matrix
        Expression matrix (genes × cells)
    gene_names : list
        List of gene names
    cell_names : list
        List of cell names
    """
    # Get expression matrix (currently cells × genes)
    # Keep as sparse if possible for efficiency
    expression = adata.X
    
    expression = expression.T  # Transpose dense array
    # Convert to sparse for min_max_normalize
    expression = csr_matrix(expression)
    
    # Get gene and cell names
    gene_names = list(adata.var_names)
    cell_names = list(adata.obs_names)
    
    # Subset genes if specified
    if genes is not None:
        # Find indices of requested genes
        gene_indices = [i for i, g in enumerate(gene_names) if g in genes]
        expression = expression[gene_indices, :]  # Use 2D indexing for sparse
        gene_names = [gene_names[i] for i in gene_indices]
        print(f"Subset to {len(gene_names)} genes")
    
    # Normalize (expects and returns sparse)
    expression = min_max_normalize(expression)
    
    # Convert to dense numpy array for plotting
    if isinstance(expression, csr_matrix):
        expression = expression.toarray()
    elif not isinstance(expression, np.ndarray):
        expression = np.array(expression)
    
    return expression, gene_names, cell_names


def plot_expression_heatmap(expression_matrix, gene_names, cell_names, condition_col='Path..Group.', 
                            condition_labels=None, region=None):
    """
    Create and display a heatmap of the expression matrix.
    
    Parameters:
    -----------
    expression_matrix : np.ndarray
        Expression matrix (genes × cells)
    gene_names : list
        List of gene names
    cell_names : list
        List of cell names
    condition_col : str
        Name of condition column (for getting condition info if needed)
    condition_labels : list, optional
        List of condition labels for each cell
    region : str, optional
        Region name for title
    """
    n_genes, n_cells = expression_matrix.shape
    
    # Create figure and plot the expression matrix as an image (genes x cells)
    fig, ax = plt.subplots(figsize=(min(20, n_cells * 0.05), min(15, n_genes * 0.3)))
    im = ax.imshow(expression_matrix, aspect='auto', cmap='viridis', interpolation='nearest')
    
    # Set y-axis to show all gene names
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(gene_names, fontsize=8)
    
    # Add condition boundaries if condition_labels provided
    if condition_labels:
        unique_conditions = []
        condition_positions = []
        current_condition = None
        for i, cond in enumerate(condition_labels):
            if cond != current_condition:
                if current_condition is not None:
                    condition_positions.append(i)
                unique_conditions.append(cond)
                current_condition = cond
        condition_positions.append(len(condition_labels))
        
        # Debug: Print condition boundaries
        print(f"Found {len(unique_conditions)} conditions: {unique_conditions}")
        print(f"Condition boundaries at positions: {condition_positions}")
        print(f"Drawing {len(condition_positions) - 1} lines at positions: {condition_positions[:-1]}")
        
        # Draw vertical lines between conditions (all boundaries except the last one which is just the end)
        for pos in condition_positions[:-1]:  # Draw at all boundaries except the final end position
            ax.axvline(x=pos, color='black', linewidth=1, zorder=10)
    
    # Add colorbar, labels, and title
    plt.colorbar(im, ax=ax, label="Expression")
    ax.set_xlabel('Cells (ordered by condition)', fontsize=12)
    ax.set_ylabel('Genes', fontsize=12)
    title = 'Expression Matrix (Genes x Cells)'
    if region:
        title += f' - {region}'
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.show()

def main():
    parser = argparse.ArgumentParser(
        description='Create gene-by-cell expression matrices for each region, ordered by condition'
    )
    parser.add_argument(
        '--regions',
        type=str,
        nargs='+',
        default=REGIONS,
        help=f'Regions to process (default: all regions: {REGIONS})'
    )
    parser.add_argument(
        '-use_DEGs',
        action='store_true',
        help='Use DEGs instead of all genes'
    )
    parser.add_argument(
        '--cell_type',
        type=str,
        default='Astrocytes',
        choices=['Astrocytes', 'Microglia'],
        help='Cell type to process (default: Astrocytes)'
    )
    
    args = parser.parse_args()
    
    # Process each region
    for region in args.regions:
        try:
            # Load data
            adata = load_region_data(region, cell_type=args.cell_type)
            print(f"After loading: {adata.shape} (cells × genes)")
            
            # Determine condition column based on cell type
            if args.cell_type == "Microglia":
                condition_col = 'Pathology_Stage'
            else:
                condition_col = 'Path..Group.'
            
            # Order cells by condition
            adata_ordered, condition_order = order_cells_by_condition(adata, condition_col=condition_col, cell_type=args.cell_type)
            print(f"After ordering cells: {adata_ordered.shape} (cells × genes)")

            if args.use_DEGs:
                print("Using DEGs...")
                gene_list = get_DEGs(region, top_k=5000, cell_type=args.cell_type)
                print(f"Found {len(gene_list)} DEGs")
                if len(gene_list) == 0:
                    print("Warning: get_DEGs returned empty list! Using all genes instead.")
                else:
                    mask = adata_ordered.var_names.isin(gene_list)
                    print(f"Matching genes in adata: {mask.sum()} out of {len(adata_ordered.var_names)}")
                    if mask.sum() == 0:
                        print("Warning: No genes matched! Using all genes instead.")
                    else:
                        adata_ordered = adata_ordered[:, mask].copy()
                        print(f"After DEG filtering: {adata_ordered.shape} (cells × genes)")
            
            # Create expression matrix
            print(f"Creating expression matrix from adata with shape: {adata_ordered.shape}")
            expression_matrix, gene_names, cell_names = create_expression_matrix(
                adata_ordered
            )
            print(f"Expression matrix created: {expression_matrix.shape} (genes × cells)")

            # Get metadata from AnnData object and ensure cell_annotation column exists
            metadata = adata_ordered.obs.copy()
            metadata['cell_annotation'] = adata_ordered.obs_names  # Add cell_annotation column
            print("Filtering...")
            if args.cell_type == "Astrocytes":
                expression_matrix, cell_names, metadata = filter_cells(expression_matrix, cell_names, metadata, extra_filter=False)
            
            # Update condition_labels after filtering to match filtered cells
            # filter_cells returns metadata with cell_annotation as index, aligned with filtered cell_names
            # Reindex to ensure order matches (in case of any misalignment)
            metadata = metadata.reindex(cell_names)
            condition_labels = metadata[condition_col].astype(str).tolist()
            
            print(f"\nExpression matrix shape: {expression_matrix.shape} (genes × cells)")
            plot_expression_heatmap(
                expression_matrix,
                gene_names,
                cell_names,
                condition_labels=condition_labels,
                region=region
            )
            exit()
            
        except Exception as e:
            print(f"Error processing region {region}: {e}\n")
            import traceback
            traceback.print_exc()


if __name__ == "__main__":
    main()

