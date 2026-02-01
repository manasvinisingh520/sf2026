"""
Plot gene embeddings using UMAP visualization.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import umap
import argparse
from utils import create_color_mapping


def load_gene_embeddings(embeddings_path):
    """Load gene embeddings from TSV/CSV file."""
    try:
        df = pd.read_csv(embeddings_path, sep=',', header=None)
    except:
        df = pd.read_csv(embeddings_path, sep='\t', header=None)
    
    embeddings_list = []
    gene_ids = []
    
    # Skip the first row (header)
    for idx in range(1, len(df)):
        row = df.iloc[idx]
        gene_id = str(row.iloc[0])
        embedding_values = row.iloc[1:].values
        embedding_values = [float(x) for x in embedding_values if pd.notna(x)]
        if len(embedding_values) > 0:
            embeddings_list.append(embedding_values)
            gene_ids.append(gene_id)
    
    return np.array(embeddings_list), gene_ids


def concatenate_all_gene_embeddings():
    """
    Automatically find and combine all gene embeddings from all regions for both Astrocytes and Microglia.
    
    Returns:
    --------
    concatenated_embeddings : np.ndarray
        Concatenated embedding matrix (genes × dimensions)
    all_gene_ids : list
        List of gene IDs for each gene
    gene_labels : list
        List of region labels for each gene
    cell_types : list
        List of cell types (Astrocytes/Microglia) for each gene
    """
    embeddings_dir = os.path.join("GREmLN", "embeddings")
    if not os.path.exists(embeddings_dir):
        raise ValueError(f"Embeddings directory not found: {embeddings_dir}")
    
    # Find all gene embedding files
    import glob
    embedding_files = glob.glob(os.path.join(embeddings_dir, "*_gene_embeddings.tsv"))
    embedding_files.extend(glob.glob(os.path.join(embeddings_dir, "*_gene_embeddings.csv")))
    
    print(f"\nFound {len(embedding_files)} gene embedding file(s)...")
    all_embeddings = []
    all_gene_ids = []
    all_gene_labels = []
    all_cell_types = []
    
    for emb_path in sorted(embedding_files):
        filename = os.path.basename(emb_path)
        # Extract region and cell type from filename
        # e.g., "EC_gene_embeddings.tsv" (Astrocytes) or "EC_Microglia_gene_embeddings.tsv" (Microglia)
        base_name = filename.replace('_gene_embeddings.tsv', '').replace('_gene_embeddings.csv', '')
        
        # Determine cell type
        if 'Microglia' in filename:
            cell_type = 'Microglia'
            # Extract region (remove _Microglia from base_name)
            region = base_name.replace('_Microglia', '').replace('Microglia', '').strip('_')
        else:
            cell_type = 'Astrocytes'
            region = base_name.replace('_TF', '')  # Remove _TF if present
        
        print(f"  Loading: {filename} (region: {region}, cell_type: {cell_type})")
        embeddings, gene_ids = load_gene_embeddings(emb_path)
        
        # Skip empty embeddings
        if embeddings.shape[0] == 0:
            print(f"    Warning: Skipping empty embedding file")
            continue
        
        # Check embedding dimensions
        if len(all_embeddings) > 0:
            expected_dim = all_embeddings[0].shape[1]
            if embeddings.shape[1] != expected_dim:
                print(f"    Warning: Skipping {filename} - dimension mismatch ({embeddings.shape[1]} vs {expected_dim})")
                continue
        
        all_embeddings.append(embeddings)
        n_genes = embeddings.shape[0]
        
        # Create gene labels (region_celltype for each gene)
        gene_labels = [f"{region}_{cell_type}"] * n_genes
        all_gene_labels.extend(gene_labels)
        all_cell_types.extend([cell_type] * n_genes)
        all_gene_ids.extend(gene_ids)
    
    if len(all_embeddings) == 0:
        raise ValueError("No valid gene embeddings found to concatenate")
    
    # Concatenate embeddings
    concatenated_embeddings = np.vstack(all_embeddings)
    print(f"\nConcatenated embeddings shape: {concatenated_embeddings.shape}")
    print(f"  Total genes: {concatenated_embeddings.shape[0]}")
    print(f"  Embedding dimensions: {concatenated_embeddings.shape[1]}")
    
    # Count genes per region and cell type
    from collections import Counter
    label_counts = Counter(all_gene_labels)
    print(f"  Genes per region_celltype: {dict(label_counts)}")
    
    cell_type_counts = Counter(all_cell_types)
    print(f"  Genes per cell type: {dict(cell_type_counts)}")
    
    return concatenated_embeddings, all_gene_ids, all_gene_labels, all_cell_types


parser = argparse.ArgumentParser(description='Plot gene embeddings using UMAP visualization.')
parser.add_argument('--region', type=str, default='EC', help='Region to process (default: EC). Can specify multiple regions separated by commas.')
parser.add_argument('--embeddings_path', type=str, default=None, help='Path to embeddings file (default: auto-detect)')
parser.add_argument('--combine_all', action='store_true', help='Combine all gene embeddings from all regions (Astrocytes and Microglia)')
parser.add_argument('--cell_type', type=str, default=None, choices=['Astrocytes', 'Microglia'],
                    help='Filter by cell type (Astrocytes or Microglia). If not specified, includes all cell types.')
parser.add_argument('--color_by', type=str, default='auto', choices=['cell_type', 'region', 'none', 'auto'],
                    help='What to color by: cell_type, region, none (single color), or auto (detect from data)')
args = parser.parse_args()

# Determine if we're loading multiple embeddings
gene_labels = None
cell_types = None

if args.combine_all:
    # Combine all gene embeddings from all regions (Astrocytes and Microglia)
    embeddings_matrix, gene_ids, gene_labels, cell_types = concatenate_all_gene_embeddings()
    region_label = "All Regions (Astrocytes and Microglia)"
else:
    # Load single region
    region = args.region
    if args.embeddings_path is None:
        embeddings_path = os.path.join("GREmLN", "embeddings", f"{region}_gene_embeddings.tsv")
    else:
        embeddings_path = args.embeddings_path
    
    embeddings_matrix, gene_ids = load_gene_embeddings(embeddings_path)
    region_label = region
    gene_labels = None
    cell_types = None

# Filter by cell_type if specified
if args.cell_type is not None:
    if cell_types is None:
        print(f"Warning: cell_type filter specified but no cell type information available. Ignoring filter.")
    else:
        # Filter embeddings, gene_ids, gene_labels, and cell_types
        mask = np.array([ct == args.cell_type for ct in cell_types])
        embeddings_matrix = embeddings_matrix[mask]
        gene_ids = [gene_ids[i] for i in range(len(gene_ids)) if mask[i]]
        if gene_labels is not None:
            gene_labels = [gene_labels[i] for i in range(len(gene_labels)) if mask[i]]
        cell_types = [cell_types[i] for i in range(len(cell_types)) if mask[i]]
        print(f"Filtered to {args.cell_type}: {len(gene_ids)} genes")

print(f"Embedding matrix shape: {embeddings_matrix.shape}")
print(f"  - Number of genes: {embeddings_matrix.shape[0]}")
print(f"  - Embedding dimensions: {embeddings_matrix.shape[1]}")

# Compute UMAP for visualization (2D)
print("\nComputing UMAP for visualization...")
reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.1)
embeddings_2d = reducer.fit_transform(embeddings_matrix)

# Determine coloring scheme based on color_by argument
from matplotlib.patches import Patch

if args.color_by == 'none':
    point_colors = 'gray'
    title = f'{region_label} Gene Embeddings - UMAP'
    color_map = None
    unique_values = None
elif args.color_by == 'cell_type' and cell_types is not None:
    # If cell_types are already simple strings (not labels with underscores), use 'generic' type with custom colors
    # Otherwise use 'cell_type' to extract from labels
    if cell_types and not any('_' in str(ct) for ct in cell_types[:10]):
        # Already extracted cell types, use generic with custom colors
        custom_colors = {
            'Astrocytes': '#1f77b4',      # Blue
            'Microglia': '#ff7f0e',       # Orange
        }
        point_colors, color_map, unique_values, _ = create_color_mapping(cell_types, attribute_type='generic', custom_colors=custom_colors)
    else:
        point_colors, color_map, unique_values, _ = create_color_mapping(cell_types, attribute_type='cell_type')
    title = f'{region_label} Gene Embeddings - UMAP (Colored by Cell Type)'
elif args.color_by == 'region' and gene_labels is not None:
    point_colors, color_map, unique_values, _ = create_color_mapping(gene_labels, attribute_type='region')
    title = f'{region_label} Gene Embeddings - UMAP (Colored by Region)'
elif args.color_by == 'auto':
    # Auto-detect: use cell_type if available, otherwise region, otherwise single color
    if cell_types is not None:
        # If cell_types are already simple strings (not labels with underscores), use 'generic' type with custom colors
        if cell_types and not any('_' in str(ct) for ct in cell_types[:10]):
            custom_colors = {
                'Astrocytes': '#1f77b4',      # Blue
                'Microglia': '#ff7f0e',       # Orange
            }
            point_colors, color_map, unique_values, _ = create_color_mapping(cell_types, attribute_type='generic', custom_colors=custom_colors)
        else:
            point_colors, color_map, unique_values, _ = create_color_mapping(cell_types, attribute_type='cell_type')
        title = f'{region_label} Gene Embeddings - UMAP (Colored by Cell Type)'
    elif gene_labels is not None:
        point_colors, color_map, unique_values, _ = create_color_mapping(gene_labels, attribute_type='region')
        title = f'{region_label} Gene Embeddings - UMAP (Colored by Region)'
    else:
        point_colors = 'gray'
        title = f'{region_label} Gene Embeddings - UMAP'
        color_map = None
        unique_values = None
else:
    # Fallback to single color
    point_colors = 'gray'
    title = f'{region_label} Gene Embeddings - UMAP'
    color_map = None
    unique_values = None

# Plot embedding
fig, ax = plt.subplots(figsize=(12, 10))
ax.scatter(embeddings_2d[:, 0], embeddings_2d[:, 1], c=point_colors, alpha=0.6, s=20)
ax.set_xlabel('UMAP Dimension 1', fontsize=12)
ax.set_ylabel('UMAP Dimension 2', fontsize=12)
ax.set_title(title, fontsize=14, fontweight='bold')
ax.set_xlim(-5, 25)
ax.set_ylim(-15, 15)
ax.grid(True, alpha=0.3)

# Add legend if we have a color mapping
if color_map is not None and unique_values is not None:
    legend_elements = [Patch(facecolor=color_map[val], label=val) for val in unique_values]
    ax.legend(handles=legend_elements, loc='best', fontsize=10)

plt.tight_layout()
plt.show()
