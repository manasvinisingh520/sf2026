"""
Cluster gene embeddings using HDBSCAN.
Clusters both regular gene embeddings and TF gene embeddings.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os
import hdbscan
from sklearn.metrics import silhouette_score
import umap


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


def cluster_with_hdbscan(embeddings, min_cluster_size=10, min_samples=5):
    """
    Cluster embeddings using HDBSCAN.
    
    Parameters:
    -----------
    embeddings : np.ndarray
        Embeddings matrix (n_samples, n_features)
    min_cluster_size : int
        Minimum cluster size for HDBSCAN (default: 10)
    min_samples : int
        Minimum samples for HDBSCAN (default: 5)
        
    Returns:
    --------
    cluster_labels : np.ndarray
        Cluster labels (-1 for noise points)
    clusterer : hdbscan.HDBSCAN
        Fitted HDBSCAN clusterer
    """
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        metric='euclidean'
    )
    cluster_labels = clusterer.fit_predict(embeddings)
    
    n_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
    n_noise = list(cluster_labels).count(-1)
    
    print(f"  Found {n_clusters} clusters")
    print(f"  Noise points: {n_noise}")
    
    # Calculate silhouette score (only if we have clusters)
    if n_clusters > 1:
        # Filter out noise points for silhouette score
        mask = cluster_labels != -1
        if mask.sum() > 1:
            try:
                sil_score = silhouette_score(embeddings[mask], cluster_labels[mask])
                print(f"  Silhouette score: {sil_score:.3f}")
            except:
                print(f"  Silhouette score: N/A (too few points)")
    
    return cluster_labels, clusterer


def plot_clustered_embeddings(embeddings_2d, cluster_labels, gene_ids, title, output_path=None):
    """
    Plot 2D embeddings colored by cluster.
    
    Parameters:
    -----------
    embeddings_2d : np.ndarray
        2D embeddings (n_samples, 2)
    cluster_labels : np.ndarray
        Cluster labels
    gene_ids : list
        List of gene IDs
    title : str
        Plot title
    output_path : str, optional
        Path to save the plot
    """
    unique_clusters = sorted(set(cluster_labels))
    n_clusters = len([c for c in unique_clusters if c != -1])
    
    # Create color mapping
    colors = plt.cm.get_cmap('tab20')
    cluster_colors = {}
    for i, cluster in enumerate(unique_clusters):
        if cluster == -1:
            cluster_colors[cluster] = 'lightgray'
        else:
            cluster_colors[cluster] = colors(i / max(n_clusters, 1))
    
    point_colors = [cluster_colors[label] for label in cluster_labels]
    
    fig, ax = plt.subplots(figsize=(12, 10))
    scatter = ax.scatter(embeddings_2d[:, 0], embeddings_2d[:, 1], 
                        c=point_colors, alpha=0.6, s=10)
    ax.set_xlabel('UMAP Dimension 1')
    ax.set_ylabel('UMAP Dimension 2')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    
    # Create legend
    legend_elements = []
    for cluster in unique_clusters:
        if cluster == -1:
            label = f'Noise ({list(cluster_labels).count(-1)} points)'
        else:
            count = list(cluster_labels).count(cluster)
            label = f'Cluster {cluster} ({count} genes)'
        legend_elements.append(Patch(facecolor=cluster_colors[cluster], label=label))
    
    ax.legend(handles=legend_elements, loc='best', fontsize=8, ncol=2)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"  Saved plot to {output_path}")
    
    plt.show()


def main():
    # Configuration
    region = "EC"
    embeddings_path = os.path.join("GREmLN", "embeddings", f"{region}_gene_embeddings.tsv")
    tf_embeddings_path = os.path.join("GREmLN", "embeddings", f"{region}_TF_gene_embeddings.tsv")
    
    # HDBSCAN parameters
    min_cluster_size = 10
    min_samples = 5
    
    # 1) Cluster regular gene embeddings
    print("=" * 60)
    print(f"Clustering {region} gene embeddings")
    print("=" * 60)
    
    print(f"\nLoading gene embeddings from {embeddings_path}...")
    embeddings_matrix, gene_ids = load_gene_embeddings(embeddings_path)
    print(f"Loaded {len(gene_ids)} genes with {embeddings_matrix.shape[1]} dimensions")
    
    print(f"\nClustering with HDBSCAN (min_cluster_size={min_cluster_size}, min_samples={min_samples})...")
    cluster_labels, clusterer = cluster_with_hdbscan(embeddings_matrix, 
                                                      min_cluster_size=min_cluster_size,
                                                      min_samples=min_samples)
    
    # Compute UMAP for visualization
    print("\nComputing UMAP for visualization...")
    reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.1)
    embeddings_2d = reducer.fit_transform(embeddings_matrix)
    
    # Plot
    plot_clustered_embeddings(embeddings_2d, cluster_labels, gene_ids,
                            f'{region} Gene Embeddings - HDBSCAN Clusters',
                            output_path=f'{region}_gene_embeddings_clusters.png')
    
    # Print cluster distribution
    print(f"\nCluster distribution:")
    cluster_counts = pd.Series(cluster_labels).value_counts().sort_index()
    for cluster, count in cluster_counts.items():
        if cluster == -1:
            print(f"  Noise: {count} genes")
        else:
            print(f"  Cluster {cluster}: {count} genes")
    
    # 2) Cluster TF gene embeddings (if file exists)
    if os.path.exists(tf_embeddings_path):
        print("\n" + "=" * 60)
        print(f"Clustering {region} TF gene embeddings")
        print("=" * 60)
        
        print(f"\nLoading TF gene embeddings from {tf_embeddings_path}...")
        tf_embeddings_matrix, tf_gene_ids = load_gene_embeddings(tf_embeddings_path)
        print(f"Loaded {len(tf_gene_ids)} genes with {tf_embeddings_matrix.shape[1]} dimensions")
        
        print(f"\nClustering with HDBSCAN (min_cluster_size={min_cluster_size}, min_samples={min_samples})...")
        tf_cluster_labels, tf_clusterer = cluster_with_hdbscan(tf_embeddings_matrix,
                                                                min_cluster_size=min_cluster_size,
                                                                min_samples=min_samples)
        
        # Compute UMAP for visualization
        print("\nComputing UMAP for visualization...")
        tf_reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.1)
        tf_embeddings_2d = tf_reducer.fit_transform(tf_embeddings_matrix)
        
        # Plot
        plot_clustered_embeddings(tf_embeddings_2d, tf_cluster_labels, tf_gene_ids,
                                f'{region} TF Gene Embeddings - HDBSCAN Clusters',
                                output_path=f'{region}_TF_gene_embeddings_clusters.png')
        
        # Print cluster distribution
        print(f"\nCluster distribution:")
        tf_cluster_counts = pd.Series(tf_cluster_labels).value_counts().sort_index()
        for cluster, count in tf_cluster_counts.items():
            if cluster == -1:
                print(f"  Noise: {count} genes")
            else:
                print(f"  Cluster {cluster}: {count} genes")
    else:
        print(f"\nWarning: TF gene embeddings file not found at {tf_embeddings_path}")
        print("Skipping TF gene embeddings clustering.")


if __name__ == "__main__":
    main()
