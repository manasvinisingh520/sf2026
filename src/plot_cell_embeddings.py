"""
Plot cell embeddings with clustering and stage distribution analysis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import BoundaryNorm
import os
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import umap
import scanpy as sc


def load_cell_embeddings(embeddings_path):
    """Load cell embeddings from TSV/CSV file."""
    try:
        df = pd.read_csv(embeddings_path, sep=',', header=None)
    except:
        df = pd.read_csv(embeddings_path, sep='\t', header=None)
    
    embeddings_list = []
    for idx in range(1, len(df)):
        row = df.iloc[idx]
        embedding_values = row.iloc[1:].values
        embedding_values = [float(x) for x in embedding_values if pd.notna(x)]
        if len(embedding_values) > 0:
            embeddings_list.append(embedding_values)
        
    return np.array(embeddings_list)


def load_conditions_from_mapping(mapping_path):
    """Load conditions from patient mapping CSV file."""
    mapping_df = pd.read_csv(mapping_path)
    cell_conditions = mapping_df['condition'].values.tolist()
    cell_conditions = [str(c) if pd.notna(c) else None for c in cell_conditions]
    return cell_conditions


def create_color_mapping(cell_conditions):
    """Create color mapping for conditions."""
    unique_conditions = sorted([c for c in set(cell_conditions) if c is not None])
    
    custom_colors = {
        '1': 'green',
        '2': 'lightgreen',
        '3': 'orange',
        '4': 'red'
    }
    
    condition_to_color = {c: custom_colors.get(c, 'gray') for c in unique_conditions}
    point_colors = [condition_to_color.get(cond, 'gray') if cond is not None else 'gray' 
                    for cond in cell_conditions]
    
    return unique_conditions, condition_to_color, point_colors


def create_tf_color_mapping(cell_tfs):
    """Create color mapping for TFs using a colormap."""
    unique_tfs = sorted([tf for tf in set(cell_tfs) if tf is not None])
    n_tfs = len(unique_tfs)
    
    # Use a colormap to assign colors (cycling through tab20 if more than 20 TFs)
    cmap = plt.cm.tab20
    tf_to_color = {}
    for idx, tf in enumerate(unique_tfs):
        tf_to_color[tf] = cmap(idx % 20)
    
    point_colors = [tf_to_color.get(tf, 'gray') if tf is not None else 'gray' 
                    for tf in cell_tfs]
    
    return unique_tfs, tf_to_color, point_colors


def create_legend(unique_conditions, condition_to_color, has_unmapped=False):
    """Create legend elements for conditions that are actually present."""
    elements = [Patch(facecolor=condition_to_color[c], label=c) for c in unique_conditions]
    if has_unmapped:
        elements.append(Patch(facecolor='gray', label='Unmapped'))
    return elements


def create_tf_legend(unique_tfs, tf_to_color, has_unmapped=False, max_legend_items=20):
    """Create legend elements for TFs that are actually present."""
    # Limit legend items to avoid overcrowding
    tfs_to_show = unique_tfs[:max_legend_items]
    elements = [Patch(facecolor=tf_to_color[tf], label=tf) for tf in tfs_to_show]
    if len(unique_tfs) > max_legend_items:
        elements.append(Patch(facecolor='gray', label=f'... and {len(unique_tfs) - max_legend_items} more'))
    if has_unmapped:
        elements.append(Patch(facecolor='gray', label='Unmapped'))
    return elements


def cluster_embeddings(embeddings, n_clusters=6):
    """Cluster embeddings using KMeans."""
    clusterer = KMeans(n_clusters=n_clusters, random_state=42)
    return clusterer.fit_predict(embeddings)


def calculate_stage_percentages(cluster_labels, cell_conditions):
    """Calculate percentage of each stage within each cluster."""
    df = pd.DataFrame({
        'cluster': cluster_labels,
        'stage': cell_conditions
    })
    
    stage_totals = df['stage'].value_counts().to_dict()
    cluster_stage_counts = df.groupby(['cluster', 'stage']).size().unstack(fill_value=0)
    
    percentage_df = pd.DataFrame(index=cluster_stage_counts.index, columns=cluster_stage_counts.columns)
    for stage in cluster_stage_counts.columns:
        stage_total = stage_totals.get(stage, 0)
        if stage_total > 0:
            percentage_df[stage] = (cluster_stage_counts[stage] / stage_total) * 100
        else:
            percentage_df[stage] = 0
    
    return percentage_df.sort_index(), cluster_stage_counts


def calculate_tf_distribution(cluster_labels, cell_tfs):
    """Calculate TF distribution within each cluster."""
    df = pd.DataFrame({
        'cluster': cluster_labels,
        'tf': cell_tfs
    })
    
    # Count TFs per cluster
    cluster_tf_counts = df.groupby(['cluster', 'tf']).size().unstack(fill_value=0)
    
    # Calculate percentages: for each cluster, what % of cells have each TF
    cluster_tf_percentages = cluster_tf_counts.div(cluster_tf_counts.sum(axis=1), axis=0) * 100
    
    # Find top TF per cluster (TF with highest count)
    top_tf_per_cluster = cluster_tf_counts.idxmax(axis=1)
    top_tf_counts = cluster_tf_counts.max(axis=1)
    top_tf_percentages = cluster_tf_percentages.max(axis=1)
    
    return cluster_tf_counts, top_tf_per_cluster, top_tf_counts, top_tf_percentages


def plot_stage_distribution_line(percentage_df, title="Stage Distribution Across Clusters"):
    """Plot line graph showing distribution of each stage across clusters."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    clusters = sorted(percentage_df.index.tolist())
    stages = sorted(percentage_df.columns.tolist())
    
    stage_colors = {
        '1': 'green',
        '2': 'lightgreen',
        '3': 'orange',
        '4': 'red'
    }
    
    for stage in stages:
        percentages = [percentage_df.loc[cluster, stage] for cluster in clusters]
        color = stage_colors.get(stage, 'gray')
        ax.plot(clusters, percentages, marker='o', label=f'Stage {stage}', 
                color=color, linewidth=2, markersize=8)
    
    ax.set_xlabel('Cluster', fontsize=12)
    ax.set_ylabel('Percentage of Stage in Cluster (%)', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(clusters)
    plt.tight_layout()
    return fig


# Load data
embeddings_path = os.path.join("GREmLN", "embeddings", "EC_cell_embeddings.tsv")
embeddings_matrix = load_cell_embeddings(embeddings_path)

# Load AnnData with TF scores
tf_adata_path = os.path.join("data", "GRN", "EC_AnnData_tf_scores.h5ad")
print(f"Loading AnnData with TF scores from {tf_adata_path}...")
tf_adata = sc.read_h5ad(tf_adata_path)

# Extract top TF for each cell (should be in same order as embeddings)
if 'top_tf' not in tf_adata.obs.columns:
    raise ValueError("'top_tf' column not found in AnnData. Run find_TF_score.py first.")
    
cell_tfs = tf_adata.obs['top_tf'].values.tolist()
cell_tfs = [str(tf) if pd.notna(tf) else None for tf in cell_tfs]

unique_tfs, tf_to_color, point_colors = create_tf_color_mapping(cell_tfs)
has_unmapped = any(tf is None for tf in cell_tfs)

# Still load conditions for stage distribution analysis
mapping_path = os.path.join("data", "GRN", "Astrocytes", "EC_Astrocytes_patient_mapping.csv")
cell_conditions = load_conditions_from_mapping(mapping_path)
unique_conditions, condition_to_color, _ = create_color_mapping(cell_conditions)

if len(cell_tfs) != embeddings_matrix.shape[0]:
    raise ValueError(f"Cell count mismatch: embeddings has {embeddings_matrix.shape[0]} cells, "
                    f"but TF data has {len(cell_tfs)} cells")

# Cluster embeddings
print("Clustering embeddings...")
cluster_labels = cluster_embeddings(embeddings_matrix, n_clusters=5)
n_clusters_found = len(np.unique(cluster_labels))
print(f"Found {n_clusters_found} clusters")

# Calculate silhouette score
silhouette_avg = silhouette_score(embeddings_matrix, cluster_labels)
print(f"Silhouette score: {silhouette_avg:.4f}")

# Calculate and print TF distribution
print("\nCalculating TF distribution within clusters...")
cluster_tf_counts, top_tf_per_cluster, top_tf_counts, top_tf_percentages = calculate_tf_distribution(cluster_labels, cell_tfs)

print("\nTop TF in Each Cluster:")
for cluster in sorted(cluster_tf_counts.index):
    top_tf = top_tf_per_cluster[cluster]
    count = top_tf_counts[cluster]
    percentage = top_tf_percentages[cluster]
    print(f"  Cluster {cluster}: {top_tf} ({count} cells, {percentage:.1f}%)")

print("\nTF Distribution Across Clusters (counts):")
print(cluster_tf_counts)

# Compute UMAP
print("\nComputing UMAP...")
reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.1)
embeddings_2d_umap = reducer.fit_transform(embeddings_matrix)

# Plot embedding colored by top TF
fig, ax = plt.subplots(figsize=(10, 8))
ax.scatter(embeddings_2d_umap[:, 0], embeddings_2d_umap[:, 1], c=point_colors, alpha=0.5, s=10)
ax.set_xlabel('UMAP Dimension 1')
ax.set_ylabel('UMAP Dimension 2')
ax.set_title('EC Cell Embeddings - UMAP (Colored by Top TF)')
ax.grid(True, alpha=0.3)
ax.legend(handles=create_tf_legend(unique_tfs, tf_to_color, has_unmapped), 
          loc='best', fontsize=6, ncol=2)
plt.tight_layout()
plt.show()

# Plot embedding colored by cluster
fig, ax = plt.subplots(figsize=(10, 8))
unique_clusters = sorted(np.unique(cluster_labels))
n_clusters = len(unique_clusters)

# Create discrete boundaries for the colormap
boundaries = np.array(unique_clusters + [unique_clusters[-1] + 1]) - 0.5
norm = BoundaryNorm(boundaries, n_clusters)

scatter = ax.scatter(embeddings_2d_umap[:, 0], embeddings_2d_umap[:, 1], 
                    c=cluster_labels, cmap='tab10', norm=norm, alpha=0.6, s=20)
ax.set_xlabel('UMAP Dimension 1')
ax.set_ylabel('UMAP Dimension 2')
ax.set_title('EC Cell Embeddings - UMAP (Colored by Cluster)')
ax.grid(True, alpha=0.3)

# Create discrete colorbar
cbar = plt.colorbar(scatter, ax=ax, label='Cluster', ticks=unique_clusters, boundaries=boundaries)
cbar.set_ticklabels(unique_clusters)
plt.tight_layout()
plt.show()
