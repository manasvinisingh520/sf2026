"""Analyze intersections between GREmLN gene vocab, AnnData, top k DGE genes, and GRN.
All data sources use Ensembl ID format.
"""

import pandas as pd
import os
import argparse
import pickle
import anndata as ad
from pathlib import Path


def load_gremln_vocab(vocab_path='GREmLN/scGraphLLM/resources/gene_vocab.csv'):
    """Load GREmLN gene vocabulary (Ensembl IDs)."""
    vocab_df = pd.read_csv(vocab_path)
    vocab_genes = set(vocab_df['gene_name'].values)
    # Remove special tokens
    special_tokens = {'<PAD>', '<MASK>', '<CLS>'}
    vocab_genes = vocab_genes - special_tokens
    return vocab_genes


def load_anndata(region, cell_type='Astrocytes'):
    """Load AnnData object for specified region (Ensembl IDs)."""
    anndata_path = os.path.join('data', 'GRN', cell_type, f'{region}_AnnData_ensembl.h5ad')
    if not os.path.exists(anndata_path):
        raise FileNotFoundError(f"AnnData file not found: {anndata_path}")
    
    adata = ad.read_h5ad(anndata_path)
    adata_genes = set(adata.var_names)
    return adata, adata_genes


def load_top_k_genes(region, top_k, cell_type='Astrocytes'):
    """Load top k genes from DGE results and convert to Ensembl IDs."""
    base_dir = Path("results/dge_final")
    dge_files = sorted(base_dir.glob(f"dge_results_{region}_*.csv")) if base_dir.exists() else []
    print(f"Found {len(dge_files)} DGE files")
    
    if not dge_files:
        print("WARNING: No DGE results files found!")
        return set()
    
    all_results = [pd.read_csv(f, index_col=0) for f in dge_files]
    combined_results = pd.concat(all_results, axis=0).sort_values('padj', ascending=True, na_position='last')
    combined_results = combined_results[~combined_results.index.duplicated(keep='first')]
    print(f"  Union: {len(combined_results)} unique genes")
    
    # Get top k genes
    top_k_genes = combined_results[combined_results['padj'].notna()].sort_values('padj', ascending=True).head(top_k)
    top_k_gene_symbols = set(top_k_genes.index)
    print(f"  Top {top_k} genes selected (from {len(combined_results[combined_results['padj'].notna()])} with valid padj)")
    
    # Convert gene symbols to Ensembl IDs using pickle dictionary
    mapping_file = os.path.join('data', 'genes_2_ensembl_ids.pkl')
    print(f"Loading Ensembl ID mapping from: {mapping_file}")
    with open(mapping_file, 'rb') as f:
        symbol_to_ensembl = pickle.load(f)
    print(f"  Loaded {len(symbol_to_ensembl)} gene mappings from pickle file")
    
    top_k_ensembl_ids = {s if s.startswith('ENSG') else symbol_to_ensembl.get(s) 
                          for s in top_k_gene_symbols if s.startswith('ENSG') or symbol_to_ensembl.get(s)}
    unmapped = [s for s in top_k_gene_symbols if not s.startswith('ENSG') and s not in symbol_to_ensembl]
    if unmapped:
        print(f"  Unmapped ({len(unmapped)}): {unmapped[:10]}")
    
    print(f"  Converted to {len(top_k_ensembl_ids)} Ensembl IDs")
    return top_k_ensembl_ids


def load_gene_embeddings(region):
    """Load gene embeddings and extract Ensembl IDs."""
    embeddings_path = os.path.join("GREmLN", "embeddings", f"{region}_gene_embeddings.tsv")
    if not os.path.exists(embeddings_path):
        # Try CSV format
        embeddings_path = os.path.join("GREmLN", "embeddings", f"{region}_gene_embeddings.csv")
        if not os.path.exists(embeddings_path):
            raise FileNotFoundError(f"Gene embeddings file not found for region: {region}")
    
    df_embeddings = pd.read_csv(embeddings_path, sep=',', header=0)
    first_col = df_embeddings.columns[0]
    emb_genes = {str(gene_id).strip() for gene_id in df_embeddings[first_col].values 
                 if str(gene_id).startswith('ENSG')}
    return emb_genes


def load_grn(region, cell_type='Astrocytes'):
    """Load GRN file and extract all genes (Ensembl IDs)."""
    grn_path = os.path.join('data', 'GRN', cell_type, f'{region}_GRN_ensembl.tsv')
    print(f"Loading GRN from: {grn_path}")
    if not os.path.exists(grn_path):
        raise FileNotFoundError(f"GRN file not found: {grn_path}")
    
    df_grn = pd.read_csv(grn_path, sep='\t', header=None, names=['TF', 'Target', 'Score'])
    grn_genes = set(df_grn['TF'].values) | set(df_grn['Target'].values)
    grn_tfs = set(df_grn['TF'].values)
    grn_targets = set(df_grn['Target'].values)
    
    print(f"  Loaded {len(df_grn):,} edges")
    print(f"  Total unique genes: {len(grn_genes)}")
    print(f"  TFs: {len(grn_tfs)}, Targets: {len(grn_targets)}")
    
    return grn_genes, grn_tfs, grn_targets


def main():
    parser = argparse.ArgumentParser(
        description='Analyze intersections between GREmLN vocab, AnnData, top k DGE genes, and GRN',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--region',
        type=str,
        default='EC',
        help='Region to process (default: EC)'
    )
    parser.add_argument(
        '--top_k',
        type=int,
        default=500,
        help='Number of top genes from DGE union (default: 500)'
    )
    parser.add_argument(
        '--cell_type',
        type=str,
        default='Astrocytes',
        help='Cell type (default: Astrocytes)'
    )
    parser.add_argument(
        '--vocab_path',
        type=str,
        default='GREmLN/scGraphLLM/resources/gene_vocab.csv',
        help='Path to gene vocabulary CSV file'
    )
    
    args = parser.parse_args()
    
    # Load data sources
    vocab_genes = load_gremln_vocab(args.vocab_path)
    adata, adata_genes = load_anndata(args.region, args.cell_type)
    top_k_ensembl_ids = load_top_k_genes(args.region, args.top_k, args.cell_type)
    emb_genes = load_gene_embeddings(args.region)
    
    # Calculate intersection
    vocab_adata_intersection = vocab_genes & adata_genes
    
    print(f"GREmLN Vocab genes: {len(vocab_genes)}")
    print(f"AnnData genes: {len(adata_genes)}")
    print(f"Top {args.top_k} DGE genes: {len(top_k_ensembl_ids)}")
    print(f"Intersection (Vocab & AnnData): {len(vocab_adata_intersection)}")
    print(f"  Percentage of Vocab: {len(vocab_adata_intersection)/len(vocab_genes)*100:.1f}%")
    print(f"  Percentage of AnnData: {len(vocab_adata_intersection)/len(adata_genes)*100:.1f}%")
    
    # Intersection between top k genes and Vocab & AnnData intersection
    print(f"\nIntersection: Top {args.top_k} DGE genes & Vocab & AnnData ({len(vocab_adata_intersection)} genes)")
    print("-" * 80)
    topk_vocab_adata = top_k_ensembl_ids & vocab_adata_intersection
    print(f"Intersection size: {len(topk_vocab_adata)} genes")
    print(f"  Percentage of Top {args.top_k}: {len(topk_vocab_adata)/len(top_k_ensembl_ids)*100:.1f}%")
    print(f"  Percentage of Vocab & AnnData ({len(vocab_adata_intersection)}): {len(topk_vocab_adata)/len(vocab_adata_intersection)*100:.1f}%")
    
    # Intersection between top k genes and gene embeddings
    print(f"\nIntersection: Top {args.top_k} DGE genes & Gene Embeddings ({len(emb_genes)} genes)")
    print("-" * 80)
    topk_emb = top_k_ensembl_ids & emb_genes
    print(f"Intersection size: {len(topk_emb)} genes")
    print(f"  Percentage of Top {args.top_k}: {len(topk_emb)/len(top_k_ensembl_ids)*100:.1f}%")
    print(f"  Percentage of Gene Embeddings ({len(emb_genes)}): {len(topk_emb)/len(emb_genes)*100:.1f}%")
    
    # Show what's missing
    vocab_only = vocab_genes - adata_genes
    adata_only = adata_genes - vocab_genes
    
    print(f"\nVocab only (not in AnnData): {len(vocab_only)}")
    print(f"AnnData only (not in Vocab): {len(adata_only)}")


if __name__ == "__main__":
    main()
