"""
Load region-specific matrix/gene/cell files, filter cells by mitochondrial percentage
from metadata, then bin filtered cells by Path..Group. and output min-max normalized
gene expression matrices per bin.
"""

import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path
from scipy.sparse import hstack
from utils import read_mtx_file, read_excel_columns, filter_cells, get_region_file_paths
from config import REGION_TO_TAB, DEFAULT_METADATA_PATH


def min_max_normalize(matrix):
    mat = matrix.tocoo(copy=True).astype(float) #coordinate format
    # Convert to CSR for efficient row operations, then get max/min as dense arrays
    matrix_csr = matrix.tocsr()
    row_max = np.array(matrix_csr.max(axis=1).toarray()).ravel()
    row_min = np.array(matrix_csr.min(axis=1).toarray()).ravel()
    denom = row_max - row_min
    denom[denom == 0] = 1.0  # Avoid divide by zero
    mat.data = (mat.data - row_min[mat.row]) / denom[mat.row]
    return mat.tocsr()


def plot_all_genes_expression_histogram(bins, gene_names, region, max_genes=None, min_expression=None):
    cond_labels = ['Cond1', 'Cond2', 'Cond3', 'Cond4']
    combined_matrices = [] #separate matrices for each condition
    boundaries = [0] #tracks where each condition starts and ends
    labels = [] 
    
    for cond in cond_labels:
        mat_bin, _, group_label = bins[cond]
        combined_matrices.append(mat_bin)
        boundaries.append(boundaries[-1] + mat_bin.shape[1])
        labels.append(f'{cond} ({group_label})')
    
    combined_matrix = hstack(combined_matrices).tocsr() if len(combined_matrices) > 1 else combined_matrices[0]
    expression_data = combined_matrix.toarray()
    n_genes, n_cells = expression_data.shape
    
    # Filter genes based on expression thresholds
    gene_indices_to_plot = np.arange(n_genes)
    
    if min_expression is not None:
        # Filter by mean expression across all cells
        mean_expression = np.mean(expression_data, axis=1)
        mask = mean_expression >= min_expression
        gene_indices_to_plot = gene_indices_to_plot[mask]
        print(f"  Filtering by mean expression >= {min_expression}: {len(gene_indices_to_plot):,} genes pass")
    
    # Limit number of genes to plot if specified
    if max_genes is not None and len(gene_indices_to_plot) > max_genes:
        gene_indices_to_plot = gene_indices_to_plot[:max_genes]
        print(f"  Limiting to first {max_genes:,} genes: {len(gene_indices_to_plot):,} genes to plot")
    
    # Apply filters to get final data
    expression_data = expression_data[gene_indices_to_plot, :]
    gene_names_to_plot = [gene_names[i] for i in gene_indices_to_plot]
    n_genes_to_plot = len(gene_indices_to_plot)
    
    fig, ax = plt.subplots(figsize=(16, 8))
    cmap = cm.get_cmap('tab20')
    colors = [cmap((i % 20) / 20) for i in range(n_genes_to_plot)]
    
    print(f"  Plotting {n_genes_to_plot:,} genes (out of {n_genes:,} total)...")
    cell_indices = np.arange(n_cells)
    for gene_idx in range(n_genes_to_plot):
        ax.scatter(cell_indices, expression_data[gene_idx, :], color=colors[gene_idx], 
                  alpha=0.3, s=1, label=gene_names_to_plot[gene_idx] if gene_idx < 20 else None)
    
    for boundary in boundaries[1:-1]:
        ax.axvline(x=boundary, color='black', linestyle='--', linewidth=2, alpha=0.7)
    
    ax.set_xlabel('Cell Index (grouped by condition)', fontsize=12)
    ax.set_ylabel('Normalized Gene Expression', fontsize=12)
    ax.set_title(f'Gene Expression Across All Cells - {region} Region\n({n_genes_to_plot:,} genes, {n_cells:,} cells)', 
                 fontsize=14, fontweight='bold')
    ax.set_xticks([(boundaries[i] + boundaries[i+1]) / 2 for i in range(len(boundaries)-1)])
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.grid(True, alpha=0.3)
    if n_genes_to_plot <= 20:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=1)
    
    plt.tight_layout()
    plt.show()
    return fig


def bin_by_path_group(normalized_matrix, cell_names, filtered_metadata):
    """Split normalized matrix into bins (Cond1-Cond4) based on Path..Group."""
    bins = {}
    groups = filtered_metadata['Path..Group.']
    unique_groups = sorted(list(dict.fromkeys(groups)))
    cond_labels = ['Cond1', 'Cond2', 'Cond3', 'Cond4']
    for idx, group in enumerate(unique_groups):
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
        default=0.15,
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
        default=DEFAULT_METADATA_PATH,
        help='Path to metadata Excel file (default: shared workbook)'
    )
    parser.add_argument(
        '-max_genes',
        type=int,
        default=None,
        help='Maximum number of genes to plot (default: None, plots all genes)'
    )
    parser.add_argument(
        '-min_expression',
        type=float,
        default=0.1,
        help='Minimum mean expression threshold (0-1, since data is normalized). Only plot genes with mean expression >= this value (default: None, no filtering)'
    )
    args = parser.parse_args()
    
    region = args.region.upper()
    mtx_path, row_path, col_path = get_region_file_paths(region, args.data_dir)
    metadata_path = Path(args.metadata).resolve()

    metadata = read_excel_columns(
        str(metadata_path),
        columns=['cell_annotation', 'percent.mito', 'Path..Group.'],
        sheet_name=REGION_TO_TAB[region]
    )
    print(f"Loaded metadata: {metadata.shape}")
    
    matrix, gene_names, cell_names = read_mtx_file(
        str(mtx_path), str(row_path), str(col_path), transpose=False
    )
    
    print("\nFiltering cells...")
    filtered_matrix, filtered_cell_names, filtered_metadata = filter_cells(
        matrix, cell_names, metadata, mito_max=args.mito_max
    )

    print("\nNormalizing (min-max per gene)...")
    bins = bin_by_path_group(
        min_max_normalize(filtered_matrix), 
        filtered_cell_names, 
        filtered_metadata
    )

    for cond, (mat_bin, _, group_label) in bins.items():
        print(f"{cond} ({group_label}): genes x cells = {mat_bin.shape}")
    
    print("\nGenerating expression plot for all genes...")
    plot_all_genes_expression_histogram(
        bins, gene_names, region, 
        max_genes=args.max_genes,
        min_expression=args.min_expression
    )
    
    return bins


if __name__ == "__main__":
    main()

