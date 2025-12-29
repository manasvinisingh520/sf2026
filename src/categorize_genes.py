# load astrocyte data for all regions
# find union of all?
# categorize into protein coding, TFs, etc

from typing import Any
from pybiomart import Dataset
from utils import read_mtx_file, get_region_file_paths
from config import REGIONS, DATA_DIR, METADATA_DATE_PREFIX
from pathlib import Path
import pickle
import hashlib
import os
import pandas as pd
import mygene

gene_annotations_file = DATA_DIR + "/gene_annotations.pkl"
cache_file_union = DATA_DIR + "/union_genes.pkl"
unnamed_genes_file = DATA_DIR + "/unnamed_genes.pkl"

def save_gene_annotations():
    # 1. Load astrocyte data for all regions and find union of all genes (with caching)
    # Use the cache_file_union defined at the top of the file
    if os.path.exists(cache_file_union):
        print("Loading union of genes from cache...")
        with open(cache_file_union, 'rb') as f:
            gene_names = pickle.load(f)
        print(f"Loaded {len(gene_names)} genes from cache")
    else:
        print("Loading astrocyte data for all regions...")
        all_gene_sets = []

        for region in REGIONS:
            print(f"  Loading region: {region}")
            # Get file paths using utility function
            mtx_path, row_annotation_path, col_annotation_path = get_region_file_paths(
                region,
                data_dir=DATA_DIR,
                base_prefix=f"{METADATA_DATE_PREFIX}_Astrocytes_{region}"
            )
            
            # Read the MTX file to get gene names
            matrix, gene_names_region, cell_names = read_mtx_file(
                mtx_path=str(mtx_path),
                row_annotation_path=str(row_annotation_path),
                col_annotation_path=str(col_annotation_path),
                transpose=False  # Matrix will be genes × cells
            )
            
            gene_set = set(gene_names_region) if gene_names_region else set()
            all_gene_sets.append(gene_set)
            print(f"    Found {len(gene_set)} genes in {region}")

        # Find union of all genes across all regions
        union_genes = set.union(*all_gene_sets) if all_gene_sets else set()
        print(f"\nUnion of all genes across all regions: {len(union_genes)} unique genes")

        # Convert to sorted list for consistent ordering
        gene_names = sorted(list(union_genes))
        
        # Save to cache (using the cache_file_union defined at the top)
        print("Saving union of genes to cache...")
        with open(cache_file_union, 'wb') as f:
            pickle.dump(gene_names, f)

    # 2. Connect to the Ensembl Human Genes dataset
    dataset = Dataset(name='hsapiens_gene_ensembl', 
                    host='http://www.ensembl.org')

    # Query everything (no filters) to get gene biotypes
    print("\nQuerying Ensembl for gene biotypes...")
    all_annotations = dataset.query(
        attributes=['external_gene_name', 'gene_biotype']
    )

    # Filter the results to only include genes from our union
    gene_annotations = all_annotations[all_annotations['Gene name'].isin(gene_names)]
    gene_annotations = gene_annotations.drop_duplicates(subset=['Gene name'], keep='first')
    print(f"Found annotations for {len(gene_annotations)} genes")

    with open(gene_annotations_file, 'wb') as f:
        pickle.dump(gene_annotations, f)

def save_unnamed_genes(gene_annotations, union_genes):
    union_genes_set = set(union_genes)
    annotated_genes = set(gene_annotations['Gene name'].dropna().unique())
    unnamed_genes = union_genes_set - annotated_genes
    print(f"Found {len(unnamed_genes)} unnamed genes")
    with open(unnamed_genes_file, 'wb') as f:
        pickle.dump(unnamed_genes, f)

def save_mygene_results(df):
    with open(unnamed_genes_file, 'rb') as f:
        unnamed_genes = pickle.load(f)

    # Query mygene for unnamed genes
    print(f"\nQuerying mygene for {len(unnamed_genes)} unnamed genes...")
    mg = mygene.MyGeneInfo()
    unnamed_genes_list = list(unnamed_genes)
    results = mg.querymany(unnamed_genes_list, 
                           scopes='symbol,alias', 
                           fields='type_of_gene,symbol', 
                           species='human')
    
    df = pd.DataFrame(results)
    print(f"Total results returned: {len(df)}")
    
    # Check for notfound entries
    if 'notfound' in df.columns:
        # Handle NaN values - convert to boolean, treating NaN as False (found)
        notfound_mask = df['notfound'].fillna(False).astype(bool)
        notfound_count = notfound_mask.sum()
        found_count = len(df) - notfound_count
        print(f"  Found: {found_count} genes")
        print(f"  Not found: {notfound_count} genes")
        
        # Filter to only found genes (notfound == False or NaN)
        df_found = df[~notfound_mask].copy()
    else:
        found_count = len(df)
        print(f"  (notfound column not present, assuming all found)")
        df_found = df.copy()

    # Add found genes to gene_annotations
    print(f"\nAdding {len(df_found)} found genes to gene_annotations...")
    with open(gene_annotations_file, 'rb') as f:
        gene_annotations = pickle.load(f)
        
    # Create new entries for found genes
    # Deduplicate by query term - keep first occurrence
    df_found_dedup = df_found.drop_duplicates(subset='query', keep='first')
    print(f"  After deduplication: {len(df_found_dedup)} unique genes (was {len(df_found)})")
    
    new_entries = []
    found_gene_names = set()
    for _, row in df_found_dedup.iterrows():
        query_term = row.get('query', '')
        if pd.isna(query_term) or query_term == '':
            continue
        
        gene_type = row.get('type_of_gene', 'unknown')
        new_entry = {'Gene name': query_term, 'Gene type': gene_type}
        new_entries.append(new_entry)
        found_gene_names.add(query_term)
    
    if new_entries:
        new_df = pd.DataFrame(new_entries)
        gene_annotations = pd.concat([gene_annotations, new_df], ignore_index=True)
        print(f"  Added {len(new_entries)} new genes to gene_annotations")
        
        # Save updated gene_annotations
        with open(gene_annotations_file, 'wb') as f:
            pickle.dump(gene_annotations, f)
        print(f"  Saved updated gene_annotations to {gene_annotations_file}")
    else:
        print("  No new genes to add (all queries were invalid)")
    
    # Remove found genes from unnamed_genes
    print(f"\nRemoving {len(found_gene_names)} found genes from unnamed_genes...")
    unnamed_genes_updated = unnamed_genes - found_gene_names
    # Save updated unnamed_genes
    with open(unnamed_genes_file, 'wb') as f:
        pickle.dump(unnamed_genes_updated, f)
    print(f"  Saved updated unnamed_genes to {unnamed_genes_file}")

def save_gene_type(gene_prefix, gene_type, start=True):
    with open(unnamed_genes_file, 'rb') as f:
        unnamed_genes = pickle.load(f)
    
    # Find genes based on start or end match
    if start:
        new_genes = [gene for gene in unnamed_genes if str(gene).startswith(gene_prefix)]
        match_type = "starting with"
    else:
        new_genes = [gene for gene in unnamed_genes if str(gene).endswith(gene_prefix)]
        match_type = "ending with"
    
    print(f"Total unnamed genes: {len(unnamed_genes)}")
    print(f"Genes {match_type} {gene_prefix}: {len(new_genes)}")
    
    if len(new_genes) > 0:
        # Load existing gene_annotations
        with open(gene_annotations_file, 'rb') as f:
            gene_annotations = pickle.load(f)
        
        # Convert to set for deduplication check
        existing_genes = set(gene_annotations['Gene name'].values) if 'Gene name' in gene_annotations.columns else set()
        
        # Create new entries (only if not already present)
        new_entries = []
        new_genes_set = set()
        for gene in new_genes:
            gene_str = str(gene)
            if gene_str not in existing_genes:
                new_entry = {'Gene name': gene_str, 'Gene type': gene_type}
                new_entries.append(new_entry)
                new_genes_set.add(gene_str)
        
        if new_entries:
            new_df = pd.DataFrame(new_entries)
            gene_annotations = pd.concat([gene_annotations, new_df], ignore_index=True)
            print(f"\nAdded {len(new_entries)} genes {match_type} {gene_prefix} to gene_annotations (type: {gene_type})")
            
            # Save updated gene_annotations
            with open(gene_annotations_file, 'wb') as f:
                pickle.dump(gene_annotations, f)
            print(f"Saved updated gene_annotations to {gene_annotations_file}")
        else:
            print(f"\nNo new genes {match_type} {gene_prefix} to add (all already exist in gene_annotations)")
        
        # Remove genes from unnamed_genes
        unnamed_genes_updated = unnamed_genes - new_genes_set
        removed_count = len(unnamed_genes) - len(unnamed_genes_updated)
        print(f"\nRemoved {removed_count} genes {match_type} {gene_prefix} from unnamed_genes")
        print(f"Updated unnamed_genes: {len(unnamed_genes_updated)} genes (was {len(unnamed_genes)})")
        
        # Save updated unnamed_genes
        with open(unnamed_genes_file, 'wb') as f:
            pickle.dump(unnamed_genes_updated, f)
        print(f"Saved updated unnamed_genes to {unnamed_genes_file}")
    else:
        print(f"No genes {match_type} {gene_prefix} found in unnamed_genes")

def main():
    with open(gene_annotations_file, 'rb') as f:
        gene_annotations = pickle.load(f)
    
    print(f"Total genes in gene_annotations: {len(gene_annotations)}")
    print(f"\nGene type value counts:")
    if 'Gene type' in gene_annotations.columns:
        print(gene_annotations['Gene type'].value_counts())
    else:
        print("Available columns:", list(gene_annotations.columns))

if __name__ == "__main__":
    main()