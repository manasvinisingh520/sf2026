"""
Read all *_row_annotation.txt files from data/AtlasData and create a union of all gene symbols.
"""

import os
import glob
import pandas as pd
from pathlib import Path
from utils import gene_symbols_to_ensembl
import pickle


def read_gene_symbols_from_annotation_file(file_path):
    """
    Read gene symbols from a row annotation file.
    Each line contains one gene symbol (may have quotes around it).
    
    Parameters:
    -----------
    file_path : str
        Path to the row annotation file
        
    Returns:
    --------
    set
        Set of unique gene symbols from the file
    """
    gene_symbols = set()
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                # Strip whitespace and quotes from each line
                gene_symbol = line.strip().strip('"').strip("'")
                if gene_symbol:  # Only add non-empty symbols
                    gene_symbols.add(gene_symbol)
    except Exception as e:
        print(f"Warning: Could not read {file_path}: {e}, skipping...")
        return set()
    
    return gene_symbols


def create_union_of_gene_symbols(data_dir="data/AtlasData", pattern="*_row_annotation.txt"):
    """
    Read all row annotation files and create a union of all gene symbols.
    
    Parameters:
    -----------
    data_dir : str
        Directory containing the annotation files
    pattern : str
        Glob pattern to match annotation files
        
    Returns:
    --------
    set
        Union of all unique gene symbols across all files
    """
    # Find all matching files
    search_pattern = os.path.join(data_dir, pattern)
    annotation_files = glob.glob(search_pattern)
    
    if not annotation_files:
        print(f"No files found matching pattern: {search_pattern}")
        return set()
    
    print(f"Found {len(annotation_files)} annotation files")
    print("-" * 60)
    
    all_gene_symbols = set()
    file_counts = {}
    
    # Read each file and collect gene symbols
    for file_path in sorted(annotation_files):
        file_name = os.path.basename(file_path)
        print(f"Reading: {file_name}...", end=" ")
        
        gene_symbols = read_gene_symbols_from_annotation_file(file_path)
        file_counts[file_name] = len(gene_symbols)
        
        # Add to union
        all_gene_symbols.update(gene_symbols)
        
        print(f"{len(gene_symbols)} unique genes (union now: {len(all_gene_symbols)})")
    
    print("-" * 60)
    print(f"\nTotal unique gene symbols across all files: {len(all_gene_symbols)}")
    print(f"\nGenes per file:")
    for file_name, count in sorted(file_counts.items()):
        print(f"  {file_name}: {count} genes")
    
    return all_gene_symbols


def main():
    """Main function to create union of gene symbols."""
    data_dir = "data/AtlasData"
    
    print("Creating union of gene symbols from row annotation files...")
    print("=" * 60)
    
    # Create union
    union_genes = create_union_of_gene_symbols(data_dir)
    
    if union_genes:
        # Convert to sorted list for easier viewing
        union_genes_list = sorted(list(union_genes))
        
        # Print sample
        print(f"\nSample of gene symbols (first 10): {union_genes_list[:10]}")
        print(f"Sample of gene symbols (last 10): {union_genes_list[-10:]}")

        # Save union genes to a pickle file
        union_genes_file = os.path.join("data", "union_genes.pkl")
        os.makedirs(os.path.dirname(union_genes_file), exist_ok=True)
        with open(union_genes_file, 'wb') as f:
            pickle.dump(union_genes_list, f)
        print(f"\nSaving union genes to: {union_genes_file}")

        # Create ensemble IDs from union of gene symbols
        ensemble_ids = gene_symbols_to_ensembl(union_genes_list, \
            csv_file="data/GRN/Priors/TFs_utoronto.csv", verbose=True)

        # Print a sample of the mapping from gene symbol to Ensembl ID
        print("\nSample mapping of gene symbols to Ensembl IDs (first 10):")
        for i, (gene, ensembl_id) in enumerate(list(ensemble_ids.items())[:10]):
            print(f"  {gene}: {ensembl_id}")

        # same ensemble IDs to a pickle file
        # Optionally save to a file
        ensemble_ids_file = os.path.join("data", "genes_2_ensembl_ids.pkl")
        os.makedirs(os.path.dirname(ensemble_ids_file), exist_ok=True)

        print(f"\nSaving ensemble IDs to: {ensemble_ids_file}")

        with open(ensemble_ids_file, 'wb') as f:
            pickle.dump(ensemble_ids, f)

        print(f"Total genes: {len(ensemble_ids)}")


if __name__ == "__main__":
    main()
