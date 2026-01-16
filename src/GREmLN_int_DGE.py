"""
Check intersection between V2 gene embeddings Ensembl IDs and top k genes from V2 DGE union
"""

import pandas as pd
import os
import argparse
from pathlib import Path

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Check intersection between V2 gene embeddings Ensembl IDs and top k genes from V2 DGE union')
parser.add_argument('--top_k', type=int, default=100, help='Number of top genes to select from DGE union (default: 100)')
args = parser.parse_args()
top_k = args.top_k

# 1. Read embeddings file and extract Ensembl IDs
embeddings_path = os.path.join("GREmLN", "embeddings", "V2_gene_embeddings.tsv")
print(f"Reading embeddings file: {embeddings_path}")

# Read the file - first column contains Ensembl ID,embedding_values
df_embeddings = pd.read_csv(embeddings_path, sep='\t', header=None)

# Extract Ensembl IDs from first column (skip first row which is header)
# Each row has format: ENSG00000127603,0.11866271,-1.6478044,...
ensembl_ids_embeddings = []
for idx, row in df_embeddings.iterrows():
    if idx == 0:  # Skip header row
        continue
    first_col = str(row.iloc[0])
    if first_col and first_col.startswith('ENSG'):
        # Split by comma and take first part (Ensembl ID)
        ensembl_id = first_col.split(',')[0]
        ensembl_ids_embeddings.append(ensembl_id)

ensembl_ids_embeddings = set(ensembl_ids_embeddings)
print(f"Found {len(ensembl_ids_embeddings)} Ensembl IDs in embeddings file")

# 2. Find V2 DGE results files
# Check common results directories
base_dir = "results/dge_final"

dge_files = []
if os.path.exists(base_dir):
    dge_dir = Path(base_dir)
    # Pattern for V2 DGE results (multiple comparisons)
    pattern = "dge_results_V2_*.csv"
    matching_files = sorted(list(dge_dir.glob(pattern)))
    if matching_files:
        print(f"Found {len(matching_files)} V2 DGE files in {base_dir}")
        dge_files.extend(matching_files)

print(f"\nUsing {len(dge_files)} DGE result files:")

# 3. Load all DGE files
all_results = []
for file_path in dge_files:
    try:
        df = pd.read_csv(file_path, index_col=0)
        all_results.append(df)
    except Exception as e:
        print(f"Warning: Could not read {file_path.name}: {e}")

# 4. Create union of all results (after loading all files)
if all_results:
    print(f"\nCreating union from {len(all_results)} DGE result files...")
    # Combine all results - take union of genes, keep best padj for each gene
    combined_results = pd.concat(all_results, axis=0)
    # For genes that appear in multiple files, keep the one with lowest padj
    combined_results = combined_results.sort_values('padj', ascending=True, na_position='last')
    combined_results = combined_results[~combined_results.index.duplicated(keep='first')]

    print(f"Union of {len(all_results)} files: {len(combined_results)} unique genes")

    # 5. Get top k genes from union (sorted by padj, with valid padj)
    # Filter to genes with valid padj
    valid_padj = combined_results[combined_results['padj'].notna()].copy()
    valid_padj = valid_padj.sort_values('padj', ascending=True)
    top_k_genes = valid_padj.head(top_k)

    print(f"Top {top_k} genes (by padj) from union:")
    print(f"  Genes with valid padj: {len(valid_padj)}")
    print(f"  Top {top_k} genes selected")

    # Get gene symbols from top k genes (DGE results use gene symbols)
    top_k_gene_symbols = set(top_k_genes.index)

    # 6. Convert gene symbols to Ensembl IDs using Ensembl REST API
    print(f"\nConverting gene symbols to Ensembl IDs using Ensembl REST API...")
    import requests
    import time
    
    gene_symbols_list = list(top_k_gene_symbols)
    symbol_to_ensembl = {}
    
    # Ensembl REST API endpoint for ID conversion
    server = "https://rest.ensembl.org"
    ext = "/xrefs/symbol/homo_sapiens/{symbol}?"
    headers = {"Content-Type": "application/json"}
    
    print(f"  Converting {len(gene_symbols_list)} gene symbols...")
    converted = 0
    failed = 0
    
    for i, symbol in enumerate(gene_symbols_list):
        if i > 0 and i % 10 == 0:
            print(f"    Progress: {i}/{len(gene_symbols_list)} ({converted} converted, {failed} failed)")
            time.sleep(0.1)  # Rate limiting
        
        try:
            url = server + ext.format(symbol=symbol)
            r = requests.get(url, headers=headers)
            
            if r.ok:
                data = r.json()
                # Find Ensembl gene ID (type == "gene")
                for item in data:
                    if item.get('type') == 'gene' and 'id' in item:
                        ensembl_id = item['id']
                        if ensembl_id.startswith('ENSG'):
                            symbol_to_ensembl[symbol] = ensembl_id
                            converted += 1
                            break
                else:
                    failed += 1
            else:
                failed += 1
        except Exception as e:
            failed += 1
            if i < 3:  # Show first few errors
                print(f"    Warning: Could not convert {symbol}: {e}")
    
    print(f"  Successfully mapped {converted}/{len(gene_symbols_list)} gene symbols")
    if failed > 0:
        print(f"  Failed to map {failed} gene symbols")
    
    # Convert top k gene symbols to Ensembl IDs
    top_k_ensembl_ids = set()
    unmapped_symbols = []
    for symbol in top_k_gene_symbols:
        if symbol in symbol_to_ensembl:
            top_k_ensembl_ids.add(symbol_to_ensembl[symbol])
        elif symbol.startswith('ENSG'):
            # Already an Ensembl ID
            top_k_ensembl_ids.add(symbol)
        else:
            unmapped_symbols.append(symbol)
    
    if unmapped_symbols:
        print(f"  Note: {len(unmapped_symbols)} symbols could not be converted to Ensembl IDs")
        if len(unmapped_symbols) <= 10:
            print(f"    Unmapped: {unmapped_symbols}")

    # 7. Check intersection
    intersection = ensembl_ids_embeddings.intersection(top_k_ensembl_ids)
    
    print(f"\n{'='*60}")
    print("INTERSECTION ANALYSIS")
    print(f"{'='*60}")
    print(f"Ensembl IDs in embeddings: {len(ensembl_ids_embeddings)}")
    print(f"Ensembl IDs in top {top_k} DGE genes (after conversion): {len(top_k_ensembl_ids)}")
    print(f"Intersection: {len(intersection)} genes")
    if len(top_k_ensembl_ids) > 0:
        print(f"Intersection percentage: {len(intersection)/len(top_k_ensembl_ids)*100:.1f}% of top {top_k} genes")
    
    if intersection:
        print(f"\nIntersecting genes:")
        # Need to reverse map Ensembl IDs back to gene symbols for display
        ensembl_to_symbol = {v: k for k, v in symbol_to_ensembl.items()} if 'symbol_to_ensembl' in locals() else {}
        for ensembl_id in sorted(intersection):
            # Find which gene symbol this corresponds to
            gene_symbol = None
            for sym, ens_id in symbol_to_ensembl.items():
                if ens_id == ensembl_id:
                    gene_symbol = sym
                    break
            if gene_symbol and gene_symbol in top_k_genes.index:
                gene_info = top_k_genes.loc[gene_symbol]
                print(f"  {ensembl_id} ({gene_symbol}): padj={gene_info['padj']:.4e}, log2FC={gene_info.get('log2FoldChange', 'N/A'):.3f}")
            else:
                print(f"  {ensembl_id}")
    else:
        print("\nNo intersection found.")
        print(f"\nNote: This could mean:")
        print(f"  1. Gene symbols could not be converted to Ensembl IDs")
        print(f"  2. The top {top_k} DGE genes are not in the embeddings file")
        print(f"  3. There's a mismatch in gene ID formats")
else:
    print("\nERROR: No DGE results files were successfully loaded!")
