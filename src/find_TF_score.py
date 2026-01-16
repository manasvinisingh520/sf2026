import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import issparse
import os
import pickle
from convert_to_ensemble_IDs import create_symbol_to_ensembl_mapping


def load_grn_network(grn_path):
    """
    Load GRN network from ARACNe output (alternating format).
    Format: TF on one line, target+MI on next line (no headers).
    
    Parameters:
    -----------
    grn_path : str
        Path to GRN TSV file
        
    Returns:
    --------
    edges : pd.DataFrame
        DataFrame with columns 'tf', 'target', 'mi'
    """
    edges_list = []
    with open(grn_path, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]  # Remove empty lines
    
    i = 0
    while i < len(lines):
        tf = lines[i]
        if i + 1 < len(lines):
            target_line = lines[i + 1].split()
            if len(target_line) >= 2:
                target = target_line[0]
                mi = float(target_line[1])
                edges_list.append({'tf': tf, 'target': target, 'mi': mi})
            i += 2
        else:
            break
    
    edges = pd.DataFrame(edges_list)
    print(f"Loaded {len(edges)} edges from network file")
    return edges


def build_regulons(edges, genes, topk=200):
    """
    Build regulons from edges.
    
    Parameters:
    -----------
    edges : pd.DataFrame
        DataFrame with columns 'tf', 'target', 'mi'
    genes : set
        Set of genes present in the data
    topk : int
        Keep top-k targets per TF (default: 200)
        
    Returns:
    --------
    regulons : dict
        Dictionary mapping TF -> list of (target, weight) tuples
    """
    # Filter network - keep only targets that are actually present
    edges = edges[edges["tf"].isin(genes) & edges["target"].isin(genes)].copy()
    
    # Keep top-k targets per TF
    edges = edges.sort_values(["tf", "mi"], ascending=[True, False])
    edges = edges.groupby("tf").head(topk)
    
    # Convert MI to weight (scaled)
    edges["weight"] = edges["mi"] / edges["mi"].max()
    
    # Build regulons
    regulons = {}
    for _, row in edges.iterrows():
        tf = row['tf']
        target = row['target']
        weight = row['weight']
        
        if tf not in regulons:
            regulons[tf] = []
        regulons[tf].append((target, weight))
    
    print(f"Built regulons for {len(regulons)} TFs")
    return regulons


def calculate_tf_scores_for_cells(adata, regulons):
    """
    Calculate TF activity scores for each cell.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with expression data
    regulons : dict
        Dictionary mapping TF -> list of (target, weight) tuples
        
    Returns:
    --------
    tf_scores_df : pd.DataFrame
        DataFrame with TFs as columns, cells as rows
    """
    # Create gene index mapping
    gene_to_idx = {gene: idx for idx, gene in enumerate(adata.var_names)}
    
    # Convert expression matrix to dense if sparse
    if issparse(adata.X):
        expression_matrix = adata.X.toarray()
    else:
        expression_matrix = adata.X
    
    n_cells = expression_matrix.shape[0]
    tf_scores = {}
    
    for tf, targets in regulons.items():
        scores = np.zeros(n_cells)
        total_weight = 0.0
        
        for target_gene, weight in targets:
            if target_gene in gene_to_idx:
                gene_idx = gene_to_idx[target_gene]
                scores += expression_matrix[:, gene_idx] * weight
                total_weight += weight
        
        if total_weight > 0:
            scores = scores / total_weight
        
        tf_scores[tf] = scores
    
    # Convert to DataFrame
    tf_scores_df = pd.DataFrame(tf_scores, index=adata.obs_names)
    return tf_scores_df


def find_top_tf_for_cells(adata, grn_path, topk=200, save_tf_scores_path=None):
    """
    Find the strongest TF for each cell.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with expression data
    grn_path : str
        Path to GRN network file
    topk : int
        Keep top-k targets per TF (default: 200)
    save_tf_scores_path : str, optional
        Path to save TF scores CSV file. If None, doesn't save.
        
    Returns:
    --------
    adata : AnnData
        AnnData object with TF scores added:
        - adata.obs['top_tf']: TF with highest score per cell
        - adata.obs['top_tf_score']: Score of top TF per cell
        - adata.obsm['X_tf_scores']: All TF scores (cells x TFs)
        - adata.uns['tf_regulons']: Regulon info
    """
    # Load GRN network
    edges = load_grn_network(grn_path)
    
    # Build regulons
    genes = set(adata.var_names)
    regulons = build_regulons(edges, genes, topk=topk)
    
    # Calculate TF scores for cells
    tf_scores_df = calculate_tf_scores_for_cells(adata, regulons)
    
    # Save TF scores to CSV if path provided
    if save_tf_scores_path:
        tf_scores_df.to_csv(save_tf_scores_path)
        print(f"Saved TF activity scores to {save_tf_scores_path}")
    
    # Find the TF with highest score for each cell
    top_tf_per_cell = tf_scores_df.idxmax(axis=1)
    top_tf_score_per_cell = tf_scores_df.max(axis=1)
    
    # Store in AnnData
    adata.obs['top_tf'] = top_tf_per_cell.values
    adata.obs['top_tf_score'] = top_tf_score_per_cell.values
    adata.obsm['X_tf_scores'] = tf_scores_df.values
    adata.uns['tf_regulons'] = {tf: [t[0] for t in targets] for tf, targets in regulons.items()}
    
    return adata


def convert_targets_to_ensembl(edges):
    """
    Convert target gene symbols to Ensembl IDs.
    
    Parameters:
    -----------
    edges : pd.DataFrame
        DataFrame with 'target' column containing gene symbols
        
    Returns:
    --------
    edges : pd.DataFrame
        DataFrame with 'target' column converted to Ensembl IDs
    """
    # Check if targets are already Ensembl IDs
    sample_size = min(100, len(edges))
    ensembl_count = sum(1 for target in edges['target'].head(sample_size) if str(target).startswith('ENSG'))
    
    if ensembl_count > sample_size * 0.5:
        print("Targets appear to already be Ensembl IDs - skipping conversion")
        return edges
    
    print(f"Converting {len(edges)} target gene symbols to Ensembl IDs...")
    
    # Load mapping
    try:
        mapping_file = "data/GRN/Priors/TFs_utoronto.csv"
        symbol_to_ensembl = create_symbol_to_ensembl_mapping(mapping_file)
    except FileNotFoundError:
        # Try alternative paths
        alt_paths = ["data/TFs_utoronto.csv", "data/AtlasData/gene_annotations.pkl"]
        symbol_to_ensembl = None
        
        for alt_path in alt_paths:
            if os.path.exists(alt_path):
                if alt_path.endswith('.pkl'):
                    try:
                        with open(alt_path, 'rb') as f:
                            gene_annotations = pickle.load(f)
                        if 'Ensembl ID' in gene_annotations.columns:
                            mapping_df = gene_annotations[['Gene name', 'Ensembl ID']].dropna()
                            symbol_to_ensembl = dict(zip(mapping_df['Gene name'], mapping_df['Ensembl ID']))
                            print(f"Loaded mapping from {alt_path}: {len(symbol_to_ensembl)} genes")
                            break
                    except Exception as e:
                        print(f"Warning: Could not load {alt_path}: {e}")
                        continue
                else:
                    try:
                        symbol_to_ensembl = create_symbol_to_ensembl_mapping(alt_path)
                        break
                    except Exception as e:
                        continue
        
        if symbol_to_ensembl is None:
            print("Warning: Could not load gene mapping. Targets will not be converted.")
            return edges
    
    # Convert targets
    converted_targets = []
    converted_count = 0
    for target in edges['target']:
        target_str = str(target).strip()
        if target_str.startswith('ENSG'):
            converted_targets.append(target_str)
        elif target_str in symbol_to_ensembl:
            converted_targets.append(symbol_to_ensembl[target_str])
            converted_count += 1
        else:
            # Keep original if no mapping found
            converted_targets.append(target_str)
    
    edges = edges.copy()
    edges['target'] = converted_targets
    print(f"Converted {converted_count} target symbols to Ensembl IDs")
    
    return edges


def find_top_tf_for_genes(grn_path, gene_list=None):
    """
    Find the strongest TF for each gene (TF with highest MI/weight that regulates it).
    
    Parameters:
    -----------
    grn_path : str
        Path to GRN network file
    gene_list : list or set, optional
        List of genes to analyze (should be Ensembl IDs). If None, analyzes all target genes in the network.
        
    Returns:
    --------
    gene_tf_df : pd.DataFrame
        DataFrame with columns:
        - 'gene': Gene name (Ensembl ID)
        - 'top_tf': TF with strongest connection (highest MI)
        - 'top_tf_mi': MI value of the strongest connection
        - 'top_tf_weight': Weight of the strongest connection
        - 'n_regulating_tfs': Number of TFs that regulate this gene
    """
    # Load GRN network
    edges = load_grn_network(grn_path)
    
    # Check if edges is empty or doesn't have required columns BEFORE filtering
    if len(edges) == 0 or 'target' not in edges.columns:
        print("Warning: No edges found or invalid GRN file format. Returning empty DataFrame.")
        return pd.DataFrame(columns=['gene', 'top_tf', 'top_tf_mi', 'top_tf_weight', 'n_regulating_tfs'])
    
    # Convert target gene symbols to Ensembl IDs
    edges = convert_targets_to_ensembl(edges)
    
    # Filter to requested genes if provided
    if gene_list is not None:
        gene_set = set(gene_list)
        edges = edges[edges["target"].isin(gene_set)].copy()
    
    # Check if edges is empty after filtering
    if len(edges) == 0:
        print("Warning: No edges found after filtering. Returning empty DataFrame.")
        return pd.DataFrame(columns=['gene', 'top_tf', 'top_tf_mi', 'top_tf_weight', 'n_regulating_tfs'])
    
    # Normalize MI to weight (scaled)
    edges["weight"] = edges["mi"] / edges["mi"].max()
    
    # Group by target gene and find TF with highest MI/weight
    gene_tf_results = []
    
    for target_gene, group in edges.groupby("target"):
        # Find TF with highest MI
        max_idx = group["mi"].idxmax()
        top_tf = group.loc[max_idx, "tf"]
        top_tf_mi = group.loc[max_idx, "mi"]
        top_tf_weight = group.loc[max_idx, "weight"]
        n_tfs = len(group)
        
        gene_tf_results.append({
            'gene': target_gene,
            'top_tf': top_tf,
            'top_tf_mi': top_tf_mi,
            'top_tf_weight': top_tf_weight,
            'n_regulating_tfs': n_tfs
        })
    
    if len(gene_tf_results) == 0:
        print("Warning: No gene-TF pairs found. Returning empty DataFrame.")
        return pd.DataFrame(columns=['gene', 'top_tf', 'top_tf_mi', 'top_tf_weight', 'n_regulating_tfs'])
    
    gene_tf_df = pd.DataFrame(gene_tf_results)
    
    if len(gene_tf_df) > 0 and 'top_tf_mi' in gene_tf_df.columns:
        gene_tf_df = gene_tf_df.sort_values('top_tf_mi', ascending=False)
    
    print(f"Found strongest TF for {len(gene_tf_df)} genes")
    
    return gene_tf_df


if __name__ == "__main__":
    # Example: Calculate TF scores for cells
    print("=" * 60)
    print("Calculating TF scores for CELLS")
    print("=" * 60)
    
    # Read AnnData
    adata = sc.read_h5ad("data/GRN/Astrocytes/EC_AnnData.h5ad")
    
    # Normalize data: CPM (Counts Per Million) + log1p + z-score
    print("\nNormalizing expression data (CPM + log1p + z-score)...")
    sc.pp.normalize_total(adata, target_sum=1e6)  # CPM normalization (1 million)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)  # Z-score normalization (mean=0, std=1 per gene)
    print("Normalization complete")
    
    # Calculate TF scores for cells
    grn_path = "data/GRN/Astrocytes/EC_GRN.tsv"
    tf_scores_csv_path = "data/GRN/EC_cell_tf_scores.csv"
    adata = find_top_tf_for_cells(adata, grn_path, topk=200, save_tf_scores_path=tf_scores_csv_path)
    
    # Save the AnnData object
    output_path = "data/GRN/EC_cells_tf_scores.h5ad"
    print(f"\nSaving AnnData with TF scores to {output_path}...")
    adata.write_h5ad(output_path)
    print("Done!")
    
    print(f"\nSummary:")
    print(f"  Cells analyzed: {adata.n_obs}")
    print(f"  Top TF stored in: adata.obs['top_tf']")
    print(f"  Top TF scores stored in: adata.obs['top_tf_score']")
    print(f"  All TF scores stored in: adata.obsm['X_tf_scores']")
    print(f"\nTop 10 most frequent enriched TFs:")
    print(adata.obs['top_tf'].value_counts().head(10))
    