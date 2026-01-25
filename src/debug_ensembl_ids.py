"""
Debug script to check Ensembl ID mapping for EC region Astrocytes.
Returns lists of mapped and unmapped genes.
"""

import os
import pandas as pd
import pickle
import json
from pathlib import Path
from utils import (
    read_mtx_file,
    get_region_file_paths,
    gene_symbols_to_ensembl,
    get_top_k_genes
)
from config import DATA_DIR, REGION_TO_TAB, DGE_RESULTS_BASE_DIR, GENE_ANNOTATIONS_FILE

try:
    from pybiomart import Dataset
    BIOMART_AVAILABLE = True
except ImportError:
    BIOMART_AVAILABLE = False

def get_ec_astrocyte_ensembl_mapping():
    """
    Get Ensembl ID mapping for EC region Astrocytes.
    
    Returns:
    --------
    mapped_genes : list
        List of gene symbols (original format) that were successfully mapped to Ensembl IDs
    unmapped_genes : list
        List of gene symbols that were NOT mapped to Ensembl IDs
    mapping : dict
        Dictionary mapping gene symbol -> Ensembl ID (for mapped genes only)
    """
    region = "EC"
    date_prefix = "2025-11-16"
    
    
    # Get file paths
    base_prefix = f"{date_prefix}_Astrocytes_{region}"
    mtx_path, row_annotation_path, col_annotation_path = get_region_file_paths(
        region,
        data_dir=DATA_DIR,
        base_prefix=base_prefix
    )
    
    # Read the MTX file to get gene names
    matrix, gene_names, cell_names = read_mtx_file(
        mtx_path=str(mtx_path),
        row_annotation_path=str(row_annotation_path),
        col_annotation_path=str(col_annotation_path),
        transpose=False  # Keep as genes × cells
    )
    
    # Strip .### suffixes from gene names (e.g., AL392086.3 -> AL392086)
    import re
    gene_names_stripped = [re.sub(r'\.\d+$', '', gene) for gene in gene_names]
    
    # Create mapping from original names to stripped names
    original_to_stripped = dict(zip(gene_names, gene_names_stripped))
    
    # Get CSV mapping file
    csv_file = os.path.join(DATA_DIR, "GRN", "Priors", "TFs_utoronto.csv")
    if not os.path.exists(csv_file):
        csv_file = os.path.join("data", "GRN", "Priors", "TFs_utoronto.csv")
    
    # Convert gene symbols to Ensembl IDs (using stripped names)
    if os.path.exists(csv_file):
        mapping_stripped = gene_symbols_to_ensembl(gene_names_stripped, csv_file=csv_file, verbose=False)
    else:
        mapping_stripped = gene_symbols_to_ensembl(gene_names_stripped, csv_file=None, verbose=False)
    
    # Map back to original gene names
    mapping = {}
    for original_name, stripped_name in original_to_stripped.items():
        if stripped_name in mapping_stripped:
            mapping[original_name] = mapping_stripped[stripped_name]
    
    # Separate mapped and unmapped genes
    # Keep original gene symbols in mapped_genes (not Ensembl IDs)
    mapped_genes = []
    unmapped_genes = []
    
    for gene_symbol in gene_names:
        if gene_symbol in mapping:
            # Successfully mapped - keep original symbol
            mapped_genes.append(gene_symbol)
        else:
            # Not mapped
            unmapped_genes.append(gene_symbol)
    
    # Results only

    # Return mapped genes (original symbols), unmapped genes, and the mapping dictionary
    return mapped_genes, unmapped_genes, mapping


def extra_mapping(gene_symbols, verbose=True):
    """
    Use biomart (pybiomart) library to map gene symbols to Ensembl IDs.
    
    Parameters:
    -----------
    gene_symbols : list
        List of gene symbols to map
    verbose : bool
        Whether to print progress messages (default: True)
        
    Returns:
    --------
    dict
        Dictionary mapping gene symbol -> Ensembl ID
    """
    if not BIOMART_AVAILABLE:
        if verbose:
            print("Warning: pybiomart library not available. Install with: pip install pybiomart")
        return {}
    
    mapping = {}
    
    try:
        # Connect to Ensembl biomart
        dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
        
        # Try with hgnc_symbol filter first
        try:
            result = dataset.query(
                attributes=['external_gene_name', 'ensembl_gene_id'],
                filters={'hgnc_symbol': gene_symbols}
            )
        except Exception:
            # If hgnc_symbol filter doesn't work, query all and filter in pandas
            result = dataset.query(
                attributes=['external_gene_name', 'ensembl_gene_id']
            )
            if result is not None and len(result) > 0:
                gene_symbols_set = set(gene_symbols)
                result = result[result['Gene name'].isin(gene_symbols_set)]
        
        if result is not None and len(result) > 0:
            for _, row in result.iterrows():
                gene_symbol = row.get('Gene name', '')
                ensembl_id = row.get('Gene stable ID', '')
                if gene_symbol and ensembl_id:
                    if gene_symbol not in mapping:
                        mapping[gene_symbol] = ensembl_id
                
    except Exception:
        return {}
    
    return mapping


def check_dge_overlap(mapped_genes, mapping, region="EC", top_k=500, dge_results_dir=None):
    """
    Check overlap between top DGE genes and mapped genes.
    
    Parameters:
    -----------
    mapped_genes : list
        List of mapped gene symbols
    mapping : dict
        Dictionary mapping gene symbol -> Ensembl ID
    region : str
        Region name (default: "EC")
    top_k : int
        Number of top genes to get from DGE results (default: 500)
    dge_results_dir : str, optional
        Directory containing DGE results. If None, uses default from config.
        
    Returns:
    --------
    overlap : set
        Set of genes that are in both top DGE genes and mapped genes
    only_in_top_dge : set
        Set of top DGE genes that are NOT mapped
    top_dge_genes : set
        Set of top K DGE genes
    combined_dge : pd.DataFrame
        Combined DGE results DataFrame
    """
    if dge_results_dir is None:
        dge_results_dir = Path(DGE_RESULTS_BASE_DIR) / "dge_final"
    else:
        dge_results_dir = Path(dge_results_dir)
    
    # Load all DGE results for the region
    pattern = f"dge_results_{region}_*.csv"
    dge_files = sorted(list(dge_results_dir.glob(pattern)))
    
    if len(dge_files) == 0:
        print(f"  No DGE files found matching pattern: {pattern}")
        print(f"  Searched in: {dge_results_dir}")
        return set(), set(), set(), None
        
    # Collect all genes from all DGE files (union)
    all_dge_data = []
    
    for dge_file in dge_files:
        try:
            df = pd.read_csv(dge_file, index_col=0)
            # Get genes with valid padj
            df_valid = df[df['padj'].notna()].copy()
            all_dge_data.append(df_valid)
        except Exception as e:
            print(f"    Warning: Could not load {dge_file.name}: {e}")
    
    if not all_dge_data:
        print("  No valid DGE data loaded")
        return set(), set(), set(), None
    
    # Combine all DGE data (union)
    combined_dge = pd.concat(all_dge_data)
    # For genes that appear in multiple files, keep the one with best (lowest) padj
    combined_dge = combined_dge.sort_values('padj').groupby(combined_dge.index).first()
        
    # Get top K genes by padj (no filtering by padj threshold or log2FC)
    top_k_list, top_k_set = get_top_k_genes(combined_dge, k=top_k, sort_by='padj', 
                                             padj_threshold=1.0, log2fc_threshold=0)
        
    # Check overlap with mapped genes
    mapped_set = set(mapped_genes)
    overlap = top_k_set & mapped_set
    only_in_top_dge = top_k_set - mapped_set
    
    if len(overlap) > 0:
        for gene in sorted(list(overlap))[:20]:
            ensembl_id = mapping.get(gene, 'N/A')
            padj = combined_dge.loc[gene, 'padj'] if gene in combined_dge.index else 'N/A'
    
    if len(only_in_top_dge) > 0:
        for gene in sorted(list(only_in_top_dge))[:20]:
            padj = combined_dge.loc[gene, 'padj'] if gene in combined_dge.index else 'N/A'
    
    return overlap, only_in_top_dge, top_k_set, combined_dge


def print_gene_type_distribution(genes, gene_annotations_file=None):
    """
    Print gene type distribution for a set of genes.
    
    Parameters:
    -----------
    genes : list or set
        List/set of gene symbols
    gene_annotations_file : str, optional
        Path to gene_annotations.pkl file. If None, uses "data/gene_annotations.pkl".
    """
    if gene_annotations_file is None:
        # Use data/ directory (not DATA_DIR which is data/AtlasData)
        gene_annotations_file = GENE_ANNOTATIONS_FILE
    
    if not os.path.exists(gene_annotations_file):
        print(f"  Warning: Gene annotations file not found: {gene_annotations_file}")
        return
    
    with open(gene_annotations_file, 'rb') as f:
        gene_annotations = pickle.load(f)
    
    # Convert to set for faster lookup
    genes_set = set(genes)
    
    # Get gene types for genes in the set
    gene_types = []
    for gene in genes:
        gene_info = gene_annotations[gene_annotations['Gene name'] == gene]
        if len(gene_info) > 0:
            gene_type = gene_info['Gene type'].iloc[0]
            gene_types.append(gene_type)
        else:
            gene_types.append('unknown')
    
    # Count gene types
    from collections import Counter
    type_counts = Counter(gene_types)
    
    print(f"  Gene type distribution ({len(genes)} genes):")
    for gene_type, count in type_counts.most_common():
        percentage = 100 * count / len(genes)
        print(f"    {gene_type}: {count} ({percentage:.1f}%)")


if __name__ == "__main__":
    mapped, unmapped, mapping = get_ec_astrocyte_ensembl_mapping()
    
    # Check overlap with DGE results
    overlap, only_in_top500, top_500_genes, combined_dge = check_dge_overlap(
        mapped, mapping, region="EC", top_k=500
    )

    # Try extra mapping with biomart
    biomart_mapping = extra_mapping(unmapped, verbose=False)
    
    # Combine original mapping and biomart mapping into one unified mapping
    combined_mapping = mapping.copy()
    combined_mapping.update(biomart_mapping)
    
    # Update overlap with combined mapping
    mapped_set_combined = set(mapped) | set(biomart_mapping.keys())
    overlap_combined = top_500_genes & mapped_set_combined
    
    # Use combined_mapping for the rest of the analysis
    mapping = combined_mapping
    mapped = list(set(mapped) | set(biomart_mapping.keys()))
    
    # Results
    print(f"\n{'='*60}")
    print(f"RESULTS:")
    print(f"{'='*60}")
    print(f"Top 500 DGE genes mapped (original): {len(overlap)} ({100*len(overlap)/500:.1f}%)")
    print(f"Top 500 DGE genes mapped (with biomart): {len(overlap_combined)} ({100*len(overlap_combined)/500:.1f}%)")
    print(f"Improvement: +{len(overlap_combined) - len(overlap)} genes")
    
    # Load gene annotations to filter for protein coding
    gene_annotations_file = os.path.join("data", "gene_annotations.pkl")
    if os.path.exists(gene_annotations_file):
        with open(gene_annotations_file, 'rb') as f:
            gene_annotations = pickle.load(f)
        
        # Create a set of protein coding genes
        protein_coding_genes = set(
            gene_annotations[gene_annotations['Gene type'] == 'protein_coding']['Gene name'].values
        )
        
        # Filter top 500 DGE genes to only protein coding
        top_500_protein_coding = top_500_genes & protein_coding_genes
    
        # Check overlap with mapped genes
        mapped_set = set(mapped)
        protein_coding_mapped = top_500_protein_coding & mapped_set
        protein_coding_not_mapped = top_500_protein_coding - mapped_set
        
        print(f"\nProtein coding genes in top 500 that are mapped: {len(protein_coding_mapped)} ({100*len(protein_coding_mapped)/len(top_500_protein_coding) if len(top_500_protein_coding) > 0 else 0:.1f}%)")
        print(f"Protein coding genes in top 500 that are NOT mapped: {len(protein_coding_not_mapped)} ({100*len(protein_coding_not_mapped)/len(top_500_protein_coding) if len(top_500_protein_coding) > 0 else 0:.1f}%)")
    
    # Load GREmLN gene vocabulary
    gremln_vocab_file = os.path.join("GREmLN", "data", "default_gene_vocab.json")
    if os.path.exists(gremln_vocab_file):
        with open(gremln_vocab_file, 'r') as f:
            gremln_vocab = json.load(f)
        
        # Print first 10 entries
        print(f"\nFirst 10 entries in GREmLN vocabulary:")
        if isinstance(gremln_vocab, dict):
            for i, (gene, node_id) in enumerate(list(gremln_vocab.items())[:10]):
                print(f"  {gene}: {node_id}")
        elif isinstance(gremln_vocab, list):
            for i, entry in enumerate(gremln_vocab[:10]):
                print(f"  {i}: {entry}")
        
        # Extract gene names from vocab (could be keys or values depending on format)
        if isinstance(gremln_vocab, dict):
            # Check if it's gene_name -> node_id or node_id -> gene_name
            if all(isinstance(k, str) and not k.isdigit() for k in list(gremln_vocab.keys())[:10]):
                gremln_genes = set(gremln_vocab.keys())
            else:
                gremln_genes = set(gremln_vocab.values())
        elif isinstance(gremln_vocab, list):
            gremln_genes = set(gremln_vocab)
        else:
            gremln_genes = set()
        
        # Check overlap with mapped top 500 DGE genes (374 genes)
        mapped_top_500 = overlap_combined
        overlap_with_vocab = mapped_top_500 & gremln_genes
        not_in_vocab = mapped_top_500 - gremln_genes
        
        print(f"\nGREmLN Gene Vocabulary:")
        print(f"  Total genes in vocab: {len(gremln_genes)}")
        print(f"  Mapped top 500 DGE genes ({len(mapped_top_500)}): {len(overlap_with_vocab)} in vocab ({100*len(overlap_with_vocab)/len(mapped_top_500) if len(mapped_top_500) > 0 else 0:.1f}%)")
        print(f"  Mapped top 500 DGE genes NOT in vocab: {len(not_in_vocab)} ({100*len(not_in_vocab)/len(mapped_top_500) if len(mapped_top_500) > 0 else 0:.1f}%)")
    else:
        print(f"\nWarning: GREmLN gene vocabulary file not found: {gremln_vocab_file}")