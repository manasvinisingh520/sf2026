"""Load astrocyte data for all regions, find union of genes, and categorize into types."""

from pybiomart import Dataset
from utils import read_mtx_file, get_region_file_paths
from config import REGIONS, DATA_DIR, METADATA_DATE_PREFIX
import pickle
import os
import pandas as pd
import mygene

gene_annotations_file = f"{DATA_DIR}/gene_annotations.pkl"
cache_file_union = f"{DATA_DIR}/union_genes.pkl"
unnamed_genes_file = f"{DATA_DIR}/unnamed_genes.pkl"

def load_pickle(path):
    """Helper to load pickle file."""
    with open(path, 'rb') as f:
        return pickle.load(f)

def save_pickle(obj, path):
    """Helper to save pickle file."""
    with open(path, 'wb') as f:
        pickle.dump(obj, f)

def save_gene_annotations():
    """Load astrocyte data for all regions, find union of genes, and get Ensembl annotations."""
    # Load union of genes (with caching)
    if os.path.exists(cache_file_union):
        print(f"Loading union of genes from cache...")
        gene_names = load_pickle(cache_file_union)
    else:
        print("Loading astrocyte data for all regions...")
        all_gene_sets = []
        for region in REGIONS:
            print(f"  Loading region: {region}")
            mtx_path, row_annotation_path, col_annotation_path = get_region_file_paths(
                region, data_dir=DATA_DIR, base_prefix=f"{METADATA_DATE_PREFIX}_Astrocytes_{region}"
            )
            _, gene_names_region, _ = read_mtx_file(
                mtx_path=str(mtx_path), row_annotation_path=str(row_annotation_path),
                col_annotation_path=str(col_annotation_path), transpose=False
            )
            all_gene_sets.append(set(gene_names_region or []))
        
        gene_names = sorted(set.union(*all_gene_sets) if all_gene_sets else set())
        save_pickle(gene_names, cache_file_union)
        print(f"Saved union of {len(gene_names)} genes to cache")

    # Get Ensembl annotations
    print("\nQuerying Ensembl for gene biotypes...")
    dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    all_annotations = dataset.query(attributes=['external_gene_name', 'gene_biotype'])
    gene_annotations = all_annotations[all_annotations['Gene name'].isin(gene_names)]
    gene_annotations = gene_annotations.drop_duplicates(subset=['Gene name'], keep='first')
    print(f"Found annotations for {len(gene_annotations)} genes")
    save_pickle(gene_annotations, gene_annotations_file)

def save_unnamed_genes(gene_annotations, union_genes):
    """Find genes without annotations and save to file."""
    unnamed_genes = set(union_genes) - set(gene_annotations['Gene name'].dropna().unique())
    print(f"Found {len(unnamed_genes)} unnamed genes")
    save_pickle(unnamed_genes, unnamed_genes_file)

def save_mygene_results(df):
    """Query mygene for unnamed genes and update annotations."""
    unnamed_genes = load_pickle(unnamed_genes_file)
    print(f"\nQuerying mygene for {len(unnamed_genes)} unnamed genes...")
    
    mg = mygene.MyGeneInfo()
    results = mg.querymany(list(unnamed_genes), scopes='symbol,alias', 
                          fields='type_of_gene,symbol', species='human')
    df = pd.DataFrame(results)
    
    # Filter found genes
    notfound_mask = df['notfound'].fillna(False).astype(bool) if 'notfound' in df.columns else pd.Series([False] * len(df))
    df_found = df[~notfound_mask].drop_duplicates(subset='query', keep='first')
    print(f"Found {len(df_found)} genes (not found: {notfound_mask.sum()})")
    
    # Add found genes to annotations
    if len(df_found) > 0:
        gene_annotations = load_pickle(gene_annotations_file)
        new_entries = [{'Gene name': row['query'], 'Gene type': row.get('type_of_gene', 'unknown')} 
                      for _, row in df_found.iterrows() 
                      if pd.notna(row.get('query', '')) and row.get('query', '') != '']
        if new_entries:
            gene_annotations = pd.concat([gene_annotations, pd.DataFrame(new_entries)], ignore_index=True)
            save_pickle(gene_annotations, gene_annotations_file)
            found_gene_names = {e['Gene name'] for e in new_entries}
            unnamed_genes = unnamed_genes - found_gene_names
            save_pickle(unnamed_genes, unnamed_genes_file)
            print(f"Added {len(new_entries)} genes, {len(unnamed_genes)} remaining unnamed")

def save_gene_type(gene_prefix, gene_type, start=True):
    """Categorize genes by prefix/suffix and update annotations."""
    unnamed_genes = load_pickle(unnamed_genes_file)
    match_func = str.startswith if start else str.endswith
    match_type = "starting with" if start else "ending with"
    new_genes = [str(g) for g in unnamed_genes if match_func(str(g), gene_prefix)]
    
    print(f"Unnamed genes: {len(unnamed_genes)}, {match_type} '{gene_prefix}': {len(new_genes)}")
    
    if new_genes:
        gene_annotations = load_pickle(gene_annotations_file)
        existing_genes = set(gene_annotations.get('Gene name', []))
        new_genes_set = {g for g in new_genes if g not in existing_genes}
        
        if new_genes_set:
            new_df = pd.DataFrame([{'Gene name': g, 'Gene type': gene_type} for g in new_genes_set])
            gene_annotations = pd.concat([gene_annotations, new_df], ignore_index=True)
            save_pickle(gene_annotations, gene_annotations_file)
            unnamed_genes = unnamed_genes - new_genes_set
            save_pickle(unnamed_genes, unnamed_genes_file)
            print(f"Added {len(new_genes_set)} genes ({match_type} '{gene_prefix}'), {len(unnamed_genes)} remaining")
        else:
            print(f"All genes {match_type} '{gene_prefix}' already in annotations")

def main():
    """Print summary of gene annotations."""
    gene_annotations = load_pickle(gene_annotations_file)
    print(f"Total genes in gene_annotations: {len(gene_annotations)}")
    if 'Gene type' in gene_annotations.columns:
        print("\nGene type value counts:")
        print(gene_annotations['Gene type'].value_counts())
    else:
        print(f"\nAvailable columns: {list(gene_annotations.columns)}")

if __name__ == "__main__":
    main()