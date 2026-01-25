"""
Plot gene embeddings using UMAP visualization.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import umap
import argparse


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

print(f"Embedding matrix shape: {embeddings_matrix.shape}")
print(f"  - Number of genes: {embeddings_matrix.shape[0]}")
print(f"  - Embedding dimensions: {embeddings_matrix.shape[1]}")

# Compute UMAP for visualization (2D)
print("\nComputing UMAP for visualization...")
reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.1)
embeddings_2d = reducer.fit_transform(embeddings_matrix)

# Color by cell type if combining all, otherwise use single color
from matplotlib.patches import Patch

if args.combine_all and cell_types:
    # Create color mapping: Astrocytes = blue, Microglia = red
    cell_type_colors = {
        'Astrocytes': '#1f77b4',  # Blue
        'Microglia': '#d62728'     # Red
    }
    
    point_colors = [cell_type_colors[ct] for ct in cell_types]
    title = f'{region_label} Gene Embeddings - UMAP (Colored by Cell Type)'
else:
    point_colors = 'gray'
    title = f'{region_label} Gene Embeddings - UMAP'

# Plot embedding
fig, ax = plt.subplots(figsize=(12, 10))
ax.scatter(embeddings_2d[:, 0], embeddings_2d[:, 1], c=point_colors, alpha=0.6, s=20)
ax.set_xlabel('UMAP Dimension 1', fontsize=12)
ax.set_ylabel('UMAP Dimension 2', fontsize=12)
ax.set_title(title, fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)

# Add legend if combining all regions
if args.combine_all and cell_types:
    cell_type_colors = {
        'Astrocytes': '#1f77b4',  # Blue
        'Microglia': '#d62728'     # Red
    }
    legend_elements = [Patch(facecolor=cell_type_colors[ct], label=ct) for ct in ['Astrocytes', 'Microglia'] if ct in set(cell_types)]
    ax.legend(handles=legend_elements, loc='best', fontsize=10)

plt.tight_layout()
plt.show()
