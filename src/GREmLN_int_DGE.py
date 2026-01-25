"""Check intersection between gene embeddings Ensembl IDs and top k genes from DGE union."""

import pandas as pd
import os
import argparse
import pickle
from pathlib import Path

parser = argparse.ArgumentParser(description='Check intersection between gene embeddings and top k DGE genes')
parser.add_argument('--region', type=str, default='EC', help='Region to process (default: EC)')
parser.add_argument('--top_k', type=int, default=500, help='Number of top genes from DGE union (default: 500)')
args = parser.parse_args()

# 1. Extract Ensembl IDs from embeddings
embeddings_path = os.path.join("GREmLN", "embeddings", f"{args.region}_gene_embeddings.tsv")
print(f"Reading embeddings: {embeddings_path}")
df_embeddings = pd.read_csv(embeddings_path, sep=',', header=0)
# First column contains gene IDs (Ensembl IDs)
first_col = df_embeddings.columns[0]
ensembl_ids_embeddings = {str(gene_id).strip() for gene_id in df_embeddings[first_col].values 
                          if str(gene_id).startswith('ENSG')}
print(f"Found {len(ensembl_ids_embeddings)} Ensembl IDs in embeddings")

# 2. Load and union DGE results
base_dir = Path("results/dge_final")
dge_files = sorted(base_dir.glob(f"dge_results_{args.region}_*.csv")) if base_dir.exists() else []
print(f"Found {len(dge_files)} DGE files")

all_results = [pd.read_csv(f, index_col=0) for f in dge_files]
if not all_results:
    print("ERROR: No DGE results files loaded!")
    exit(1)

combined_results = pd.concat(all_results, axis=0).sort_values('padj', ascending=True, na_position='last')
combined_results = combined_results[~combined_results.index.duplicated(keep='first')]
print(f"Union: {len(combined_results)} unique genes")

# 3. Get top k genes
top_k_genes = combined_results[combined_results['padj'].notna()].sort_values('padj', ascending=True).head(args.top_k)
top_k_gene_symbols = set(top_k_genes.index)
print(f"Top {args.top_k} genes selected (from {len(combined_results[combined_results['padj'].notna()])} with valid padj)")

# 4. Convert gene symbols to Ensembl IDs using pickle dictionary
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

# 5. Check intersection
intersection = ensembl_ids_embeddings.intersection(top_k_ensembl_ids)
ensembl_to_symbol = {v: k for k, v in symbol_to_ensembl.items()}

print(f"\n{'='*60}\nINTERSECTION ANALYSIS\n{'='*60}")
print(f"Ensembl IDs in embeddings: {len(ensembl_ids_embeddings)}")
print(f"Ensembl IDs in top {args.top_k} DGE genes: {len(top_k_ensembl_ids)}")
print(f"Intersection: {len(intersection)} genes")
if top_k_ensembl_ids:
    print(f"Intersection percentage: {len(intersection)/len(top_k_ensembl_ids)*100:.1f}%")
