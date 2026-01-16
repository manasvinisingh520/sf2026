"""
Plot gene embeddings using dimensionality reduction (PCA, t-SNE, UMAP)
Colored by top 5000 DEG union genes or a provided gene list
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
import pickle
from convert_to_ensemble_IDs import create_symbol_to_ensembl_mapping
from find_TF_score import find_top_tf_for_genes


def load_gene_embeddings(embeddings_path):
    """Load gene embeddings from TSV/CSV file.
    Automatically skips the first row (header row with column indices).
    Returns embeddings array and gene IDs.
    """
    # Try comma-separated first (CSV), then tab-separated (TSV)
    try:
        df = pd.read_csv(embeddings_path, sep=',', header=None)
    except:
        df = pd.read_csv(embeddings_path, sep='\t', header=None)
    
    # Extract embeddings: each row has an index in first column, then embedding values
    # We'll take all columns after the first one as embedding dimensions
    embeddings_list = []
    gene_ids = []
    
    # Always skip the first row (row 0) - it's the header
    for idx in range(1, len(df)):
        row = df.iloc[idx]
        # First column is gene ID
        gene_id = str(row.iloc[0])
        # Skip first column (index) and take all remaining columns as embedding values
        embedding_values = row.iloc[1:].values
        # Convert to float, filtering out NaN
        embedding_values = [float(x) for x in embedding_values if pd.notna(x)]
        if len(embedding_values) > 0:
            embeddings_list.append(embedding_values)
            gene_ids.append(gene_id)
    
    if not embeddings_list:
        raise ValueError(f"No embeddings found in file. DataFrame shape: {df.shape}")
    
    embeddings_array = np.array(embeddings_list)
    
    return embeddings_array, gene_ids


def convert_genes_to_ensembl(gene_ids, mapping_file="data/GRN/Priors/TFs_utoronto.csv"):
    """
    Convert gene symbols to Ensembl IDs.
    
    Parameters:
    -----------
    gene_ids : list or set
        List of gene IDs (could be symbols or Ensembl IDs)
    mapping_file : str
        Path to mapping file (TFs_utoronto.csv or similar)
        
    Returns:
    --------
    ensembl_ids : set
        Set of Ensembl IDs (original if already Ensembl, converted if symbols)
    """
    print(f"  convert_genes_to_ensembl called with {len(gene_ids)} genes")
    gene_ids = list(gene_ids)
    
    # Check if genes are already Ensembl IDs
    sample_size = min(100, len(gene_ids))
    ensembl_count = sum(1 for gid in gene_ids[:sample_size] if str(gid).startswith('ENSG'))
    print(f"  Sample check: {ensembl_count}/{sample_size} genes are Ensembl IDs")
    
    if ensembl_count > sample_size * 0.5:  # If more than half are Ensembl IDs
        print(f"  Genes appear to already be Ensembl IDs - skipping conversion")
        return set(gene_ids)
    
    print(f"Converting {len(gene_ids)} gene symbols to Ensembl IDs...")
    
    # Load mapping
    try:
        symbol_to_ensembl = create_symbol_to_ensembl_mapping(mapping_file)
    except FileNotFoundError:
        # Try alternative path
        alt_path = "data/TFs_utoronto.csv"
        if os.path.exists(alt_path):
            symbol_to_ensembl = create_symbol_to_ensembl_mapping(alt_path)
        else:
            print(f"Warning: Could not find mapping file. Trying gene_annotations.pkl...")
            # Try gene_annotations.pkl
            try:
                gene_annotations_file = os.path.join("data", "AtlasData", "gene_annotations.pkl")
                with open(gene_annotations_file, 'rb') as f:
                    gene_annotations = pickle.load(f)
                
                # Check if it has Ensembl ID column
                if 'Ensembl ID' in gene_annotations.columns:
                    mapping_df = gene_annotations[['Gene name', 'Ensembl ID']].dropna()
                    symbol_to_ensembl = dict(zip(mapping_df['Gene name'], mapping_df['Ensembl ID']))
                    print(f"Loaded mapping from gene_annotations.pkl: {len(symbol_to_ensembl)} genes")
                else:
                    print("Warning: gene_annotations.pkl does not have Ensembl ID column")
                    return set(gene_ids)  # Return original if no mapping available
            except Exception as e:
                print(f"Warning: Could not load gene mapping: {e}")
                return set(gene_ids)  # Return original if no mapping available
    
    # Convert genes
    ensembl_ids = []
    converted_count = 0
    for gene_id in gene_ids:
        gene_str = str(gene_id).strip()
        # If already Ensembl ID, keep it
        if gene_str.startswith('ENSG'):
            ensembl_ids.append(gene_str)
        # Otherwise, try to convert
        elif gene_str in symbol_to_ensembl:
            ensembl_ids.append(symbol_to_ensembl[gene_str])
            converted_count += 1
        else:
            # Keep original if no mapping found
            ensembl_ids.append(gene_str)
    
    print(f"  Converted {converted_count} gene symbols to Ensembl IDs")
    print(f"  {len(ensembl_ids) - converted_count} genes unchanged (already Ensembl or no mapping)")
    
    return set(ensembl_ids)


def load_top_5000_deg_union(region, version='final', bins=100, padj_max=0.05, log2fc_min=0.5, cell_type=None):
    """
    Load the union of top 5000 DEG genes for a region.
    
    Parameters:
    -----------
    region : str
        Region name (e.g., 'EC', 'ITG')
    version : str
        DGE version (default: 'final')
    bins : int
        Bin number (default: 10)
    padj_max : float
        padj threshold (default: 0.05)
    log2fc_min : float
        log2fc threshold (default: 0.5)
    cell_type : str, optional
        Cell type (e.g., 'Microglia', 'Astrocytes')
        
    Returns:
    --------
    top_5000_genes : set
        Set of gene IDs (Ensembl IDs) in top 5000 DEG union
    """
    # Get results directory - use results/dge_final_protein_coding
    base_dir = Path("results/dge_final_protein_coding")
    
    if not base_dir.exists():
        print(f"Warning: DGE results directory not found: {base_dir}")
        return set()
    
    # Find all files matching region and bin size
    # File pattern: dge_results_PC_{region}_bins{bins}_seed*_*.csv
    pattern = f"dge_results_PC_{region}_bins{bins}*.csv"
    matching_files = sorted(list(base_dir.glob(pattern)))
    
    print(f"Found {len(matching_files)} DEG files for region {region}")
    
    # If no files found, return empty set
    if len(matching_files) == 0:
        print(f"Warning: No DEG files found for region {region}, version {version}, bins {bins}")
        return set()
    
    # Load all files and take union
    all_results = []
    for file_path in matching_files:
        df = pd.read_csv(file_path, index_col=0)
        all_results.append(df)
    
    # Combine all results - take union of genes, keep best padj for each gene
    combined_results = pd.concat(all_results, axis=0)
    # For genes that appear in multiple files, keep the one with lowest padj
    combined_results = combined_results.sort_values('padj', ascending=True, na_position='last')
    combined_results = combined_results[~combined_results.index.duplicated(keep='first')]
    
    print(f"Union of {len(matching_files)} files: {len(combined_results)} unique genes")
    
    # Filter by thresholds
    filtered_results = combined_results.copy()
    if padj_max is not None:
        filtered_results = filtered_results[filtered_results['padj'] < padj_max]
    if log2fc_min is not None:
        filtered_results = filtered_results[abs(filtered_results['log2FoldChange']) > log2fc_min]
    
    # Sort by padj and get top 5000
    filtered_results = filtered_results.sort_values('padj', ascending=True)
    top_5000 = filtered_results.head(5000)
    
    # Get gene IDs (may be symbols or Ensembl IDs) as set
    top_5000_genes = set(top_5000.index.tolist())
    print(f"Top 5000 DEG genes (after filtering): {len(top_5000_genes)} genes")
    
    # Convert to Ensembl IDs if needed
    print("\nConverting DEG genes to Ensembl IDs...")
    top_5000_genes_ensembl = convert_genes_to_ensembl(top_5000_genes)
    print(f"After conversion: {len(top_5000_genes_ensembl)} Ensembl IDs")
    
    return top_5000_genes_ensembl


def create_tf_color_mapping(gene_tf_dict, gene_ids):
    # Get unique TFs
    unique_tfs = sorted(set(tf for tf in gene_tf_dict.values() if tf is not None))
    
    # Create color mapping using tab20 colormap
    n_tfs = len(unique_tfs)
    if n_tfs > 0:
        cmap = plt.cm.get_cmap('tab20')
        colors_list = [cmap(i / max(n_tfs, 1)) for i in range(n_tfs)]
        tf_to_color = {tf: colors_list[i % len(colors_list)] for i, tf in enumerate(unique_tfs)}
        tf_to_color[None] = 'lightgray'  # For genes without TF
    else:
        tf_to_color = {None: 'lightgray'}
    
    # Create point colors
    point_colors = [gene_tf_dict.get(gene_id, None) for gene_id in gene_ids]
    point_colors = [tf_to_color.get(tf, 'lightgray') for tf in point_colors]
    
    return unique_tfs, tf_to_color, point_colors


def create_tf_legend(unique_tfs, tf_to_color, max_legend_items=20):
    # Limit number of legend items
    tfs_to_show = unique_tfs[:max_legend_items]
    legend_elements = [Patch(facecolor=tf_to_color[tf], label=tf) for tf in tfs_to_show]
    
    if len(unique_tfs) > max_legend_items:
        legend_elements.append(Patch(facecolor='gray', label=f'... and {len(unique_tfs) - max_legend_items} more'))
    
    return legend_elements


def plot_embedding_by_tf(ax, coords, gene_ids, gene_tf_dict, title, xlabel, ylabel):
    """Plot 2D embedding with genes colored by their top TF."""
    # Create color mapping
    unique_tfs, tf_to_color, point_colors = create_tf_color_mapping(gene_tf_dict, gene_ids)
    
    # Count genes per TF
    tf_counts = {}
    for gene_id in gene_ids:
        tf = gene_tf_dict.get(gene_id, None)
        tf_counts[tf] = tf_counts.get(tf, 0) + 1
    
    print(f"  Plotting {len(gene_ids)} genes colored by {len(unique_tfs)} unique TFs")
    
    # Plot all genes
    ax.scatter(coords[:, 0], coords[:, 1], c=point_colors, alpha=0.6, s=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    
    # Add legend
    legend_elements = create_tf_legend(unique_tfs, tf_to_color, max_legend_items=20)
    ax.legend(handles=legend_elements, loc='best', fontsize=6, ncol=2)


def main():
    # Configuration - hardcoded defaults
    embeddings_path = os.path.join("GREmLN", "embeddings", "EC_gene_embeddings.tsv")
    
    # Extract region and cell_type from filename
    filename = os.path.basename(embeddings_path)
    parts = filename.split('_')
    region = parts[0]  # Get first part (region)
    cell_type = parts[1] if len(parts) > 1 and parts[1] in ['Microglia', 'Astrocytes'] else None
    
    version = 'final'
    bins = 100
    
    # Load gene embeddings
    print("Loading gene embeddings...")
    embeddings_matrix, gene_ids = load_gene_embeddings(embeddings_path)
    print(f"Loaded {len(gene_ids)} genes with {embeddings_matrix.shape[1]} dimensions")
    
    # Load top 5000 DEG union
    print("\nLoading top 5000 DEG union...")
    print(f"Region: {region}, Cell type: {cell_type}, Version: {version}, Bins: {bins}")
    top_5000_deg_genes = load_top_5000_deg_union(region, version=version, bins=bins, cell_type=cell_type)
    print(f"Found {len(top_5000_deg_genes)} genes in top 5000 DEG union")
    
    # Convert gene_ids to Ensembl format for comparison (assuming format mismatch)
    # Embeddings are likely in symbol format, DEG results are in Ensembl format
    gene_ids_ensembl = convert_genes_to_ensembl(gene_ids)
    
    # Count overlap using Ensembl IDs
    overlap = len([gid for gid in gene_ids_ensembl if gid in top_5000_deg_genes])
    print(f"Overlap with embeddings (after format conversion): {overlap} genes")
    
    # Load TF information for genes
    print("\nLoading TF information for genes...")
    grn_path = f"data/GRN/Astrocytes/{region}_GRN_ensembl.tsv"
    if not os.path.exists(grn_path):
        # Try alternative path
        grn_path = f"data/GRN/{region}_GRN.tsv"
    
    if os.path.exists(grn_path):
        # Get top TF for all genes in embeddings
        gene_tf_df = find_top_tf_for_genes(grn_path, gene_list=gene_ids)
        
        # Create dictionary mapping gene ID -> top TF
        gene_tf_dict = dict(zip(gene_tf_df['gene'], gene_tf_df['top_tf']))
        
        # Fill in genes that don't have TF info
        for gene_id in gene_ids:
            if gene_id not in gene_tf_dict:
                gene_tf_dict[gene_id] = None
        
        print(f"Found TF information for {len([g for g, tf in gene_tf_dict.items() if tf is not None])} genes")
    else:
        print(f"Warning: GRN file not found at {grn_path}. Coloring genes by DEG instead.")
        gene_tf_dict = None
    
    # Compute UMAP
    reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.1)
    embeddings_2d_umap = reducer.fit_transform(embeddings_matrix)
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    if gene_tf_dict is not None:
        plot_embedding_by_tf(ax, embeddings_2d_umap, gene_ids, gene_tf_dict,
                            f'{region} Gene Embeddings - UMAP (colored by top TF)',
                            'UMAP Dimension 1', 'UMAP Dimension 2')
    else:
        # Fallback to DEG coloring if TF info not available
        genes_to_color = top_5000_deg_genes
        plot_label = "Top 5000 DEG"
        colors = ['red' if gene_id in genes_to_color else 'blue' for gene_id in gene_ids]
        n_red = sum(1 for c in colors if c == 'red')
        ax.scatter(embeddings_2d_umap[:, 0], embeddings_2d_umap[:, 1], c=colors, alpha=0.5, s=10)
        ax.set_xlabel('UMAP Dimension 1')
        ax.set_ylabel('UMAP Dimension 2')
        ax.set_title(f'{region} Gene Embeddings - UMAP ({plot_label} colored)')
        ax.grid(True, alpha=0.3)
        legend_elements = [
            Patch(facecolor='red', label=f'{plot_label} ({n_red} found)'),
            Patch(facecolor='blue', label='Other genes')
        ]
        ax.legend(handles=legend_elements, loc='best', fontsize=8)
    
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
