"""
Script to map gene symbols to Ensembl IDs using TFs_utoronto.csv and update GRN.tsv, matrix.csv, and AnnData.h5ad files.
Supports V1, V1, V2, and other regions.
"""

import pandas as pd
import os
import argparse
import anndata as ad
from tqdm import tqdm


def create_symbol_to_ensembl_mapping(csv_file="data/TFs_utoronto.csv"):
    """
    Create a mapping dictionary from HGNC symbol to Ensembl ID.
    
    Parameters:
    -----------
    csv_file : str
        Path to the TFs_utoronto.csv file
        
    Returns:
    --------
    dict
        Dictionary mapping HGNC symbol to Ensembl ID
    """
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"File not found: {csv_file}")
    
    print(f"Reading mapping file: {csv_file}")
    df = pd.read_csv(csv_file)
    
    # Create mapping from HGNC symbol to Ensembl ID
    # Drop rows where either column is NaN
    mapping_df = df[['HGNC symbol', 'Ensembl ID']].dropna()
    
    # Create dictionary (HGNC symbol -> Ensembl ID)
    symbol_to_ensembl = dict(zip(mapping_df['HGNC symbol'], mapping_df['Ensembl ID']))
    
    print(f"Created mapping for {len(symbol_to_ensembl)} symbols")
    
    return symbol_to_ensembl


def update_grn_file(grn_file="data/GRN/Astrocytes/V1_GRN.tsv", mapping=None, output_file="data/GRN/Astrocytes/V1_GRN_ensembl.tsv"):
    """
    Update V1_GRN.tsv by mapping first column (TF symbols) to Ensembl IDs and saving as V1_GRN_ensembl.tsv.
    
    Parameters:
    -----------
    grn_file : str
        Path to V1_GRN.tsv file
    mapping : dict
        Dictionary mapping symbol to Ensembl ID
    output_file : str
        Output file path
    """
    if mapping is None:
        mapping = create_symbol_to_ensembl_mapping()
    
    print(f"\nReading GRN file: {grn_file}")
    
    # First, get file size for progress estimation
    file_size = os.path.getsize(grn_file) / (1024 * 1024)  # Size in MB
    print(f"File size: {file_size:.2f} MB")
    
    # For large files (>50 MB), process in chunks to save memory and show progress
    if file_size > 50:
        print("Processing large file in chunks...")
        chunk_size = 100000  # Process 100k rows at a time (smaller to avoid buffer overflow)
        mapped_count = 0
        unmapped_symbols = set()
        total_rows = 0
        first_chunk = True
        
        # Open output file for writing
        with open(output_file, 'w', encoding='utf-8') as out_f:
            # Read and process in chunks with error handling
            try:
                # Try with on_bad_lines (pandas >= 1.3.0)
                try:
                    chunk_iter = pd.read_csv(
                        grn_file, 
                        sep='\t', 
                        header=None, 
                        names=['TF', 'Target', 'Score'], 
                        chunksize=chunk_size,
                        on_bad_lines='skip',  # Skip malformed lines
                        engine='python'  # Use Python engine for better error handling
                    )
                except TypeError:
                    # Fallback for older pandas versions
                    chunk_iter = pd.read_csv(
                        grn_file, 
                        sep='\t', 
                        header=None, 
                        names=['TF', 'Target', 'Score'], 
                        chunksize=chunk_size,
                        error_bad_lines=False,  # Skip malformed lines (pandas < 1.3.0)
                        warn_bad_lines=False,
                        engine='python'
                    )
                
                # Estimate total chunks for progress
                estimated_rows = int((file_size * 1024 * 1024) / 70)  # ~70 bytes per row estimate
                estimated_chunks = max(1, (estimated_rows // chunk_size) + 1)
                
                with tqdm(total=estimated_chunks, desc="Processing chunks", unit="chunk") as pbar:
                    for chunk in chunk_iter:
                        if len(chunk) == 0:
                            continue
                            
                        total_rows += len(chunk)
                        
                        # Map TF column using vV1torized operation
                        original_tf = chunk['TF'].copy()
                        chunk['TF'] = chunk['TF'].map(mapping).fillna(chunk['TF'])
                        
                        # Count mappings
                        mapped_mask = chunk['TF'] != original_tf
                        mapped_count += mapped_mask.sum()
                        unmapped_symbols.update(original_tf[~mapped_mask].unique())
                        
                        # Write chunk to output file
                        chunk.to_csv(out_f, sep='\t', header=False, index=False, mode='w' if first_chunk else 'a')
                        first_chunk = False
                        
                        pbar.update(1)
            except Exception as e:
                print(f"Error reading chunks: {e}")
                print("Falling back to line-by-line processing...")
                # Fallback: process line by line
                with open(grn_file, 'r', encoding='utf-8') as in_f:
                    for line in tqdm(in_f, desc="Processing lines", unit="line"):
                        try:
                            parts = line.strip().split('\t')
                            if len(parts) >= 3:
                                tf_symbol = parts[0]
                                target = parts[1]
                                score = parts[2]
                                
                                # Map TF symbol
                                mapped_tf = mapping.get(tf_symbol, tf_symbol)
                                if mapped_tf != tf_symbol:
                                    mapped_count += 1
                                else:
                                    unmapped_symbols.add(tf_symbol)
                                
                                # Write to output
                                out_f.write(f"{mapped_tf}\t{target}\t{score}\n")
                                total_rows += 1
                        except Exception as line_error:
                            continue  # Skip malformed lines
        
        print(f"Processed {total_rows} rows")
        print(f"  Mapped {mapped_count} symbols to Ensembl IDs")
        if len(unmapped_symbols) > 0:
            print(f"  {len(unmapped_symbols)} symbols not found in mapping (keeping original)")
            if len(unmapped_symbols) <= 20:
                print(f"  Unmapped symbols: {sorted(unmapped_symbols)}")
        
        # Return a small sample for compatibility
        df = pd.read_csv(grn_file, sep='\t', header=None, names=['TF', 'Target', 'Score'], nrows=1000)
        df['TF'] = df['TF'].map(mapping).fillna(df['TF'])
    else:
        # For smaller files, process normally
        df = pd.read_csv(grn_file, sep='\t', header=None, names=['TF', 'Target', 'Score'])
        print(f"Loaded {len(df)} rows")
        
        # Map first column (TF symbols) to Ensembl IDs
        print("Mapping TF symbols to Ensembl IDs...")
        original_tf = df['TF'].copy()
        df['TF'] = df['TF'].map(mapping).fillna(df['TF'])
        
        # Count mapped symbols
        mapped_mask = df['TF'] != original_tf
        mapped_count = mapped_mask.sum()
        unmapped_symbols = set(original_tf[~mapped_mask].unique())
        
        print(f"  Mapped {mapped_count} symbols to Ensembl IDs")
        if len(unmapped_symbols) > 0:
            print(f"  {len(unmapped_symbols)} symbols not found in mapping (keeping original)")
            if len(unmapped_symbols) <= 20:
                print(f"  Unmapped symbols: {sorted(unmapped_symbols)}")
        
        # Save to file
        print(f"Writing to: {output_file}")
        df.to_csv(output_file, sep='\t', header=False, index=False)
    
    print(f"Successfully updated {output_file}")
    
    return df


def update_matrix_file(matrix_file="data/GRN/Astrocytes/V1_matrix.csv", mapping=None, output_file="data/GRN/Astrocytes/V1_matrix_ensembl.csv"):
    """
    Update V1_matrix.csv by mapping first column (gene symbols) to Ensembl IDs and saving as V1_matrix_ensembl.csv.
    
    Parameters:
    -----------
    matrix_file : str
        Path to V1_matrix.csv file
    mapping : dict
        Dictionary mapping symbol to Ensembl ID
    output_file : str
        Output file path
    """
    if mapping is None:
        mapping = create_symbol_to_ensembl_mapping()
    
    print(f"\nReading matrix file: {matrix_file}")
    print("Reading file (this may take a while for large files)...")
    
    # Read CSV file (no header, first column is gene names)
    # Read in chunks to handle large files efficiently
    chunk_size = 10000
    chunks = []
    mapped_count = 0
    unmapped_symbols = set()
    
    def map_symbol(symbol):
        nonlocal mapped_count, unmapped_symbols
        if symbol in mapping:
            mapped_count += 1
            return mapping[symbol]
        else:
            unmapped_symbols.add(symbol)
            return symbol  # Keep original if not in mapping
    
    # Estimate total chunks for progress bar (approximate)
    file_size_mb = os.path.getsize(matrix_file) / (1024 * 1024)
    estimated_chunks = max(1, int((file_size_mb * 1024) / chunk_size))  # Rough estimate
    
    print("Processing file in chunks...")
    chunk_iter = pd.read_csv(matrix_file, header=None, chunksize=chunk_size)
    with tqdm(desc="Processing chunks", unit="chunk") as pbar:
        for chunk in chunk_iter:
            # Replace first column with mapped values
            chunk.iloc[:, 0] = chunk.iloc[:, 0].apply(map_symbol)
            chunks.append(chunk)
            pbar.update(1)
    
    # Concatenate all chunks
    df = pd.concat(chunks, ignore_index=True)
    
    print(f"Loaded {len(df)} genes")
    print("Mapping gene symbols to Ensembl IDs...")
    print(f"  Mapped {mapped_count} symbols to Ensembl IDs")
    if len(unmapped_symbols) > 0:
        print(f"  {len(unmapped_symbols)} symbols not found in mapping (keeping original)")
        if len(unmapped_symbols) <= 20:
            print(f"  First 20 unmapped symbols: {sorted(list(unmapped_symbols))[:20]}")
    
    # Save to file
    print(f"Writing to: {output_file}")
    df.to_csv(output_file, header=False, index=False)
    
    print(f"Successfully updated {output_file}")
    
    return df


def update_anndata_file(anndata_file="data/GRN/Astrocytes/V1_AnnData.h5ad", mapping=None, output_file="data/GRN/Astrocytes/V1_AnnData_ensembl.h5ad"):
    """
    Update V1_AnnData.h5ad by mapping gene names (var_names) from symbols to Ensembl IDs and saving as V1_AnnData_ensembl.h5ad.
    
    Parameters:
    -----------
    anndata_file : str
        Path to V1_AnnData.h5ad file
    mapping : dict
        Dictionary mapping symbol to Ensembl ID
    output_file : str
        Output file path
    """
    if mapping is None:
        mapping = create_symbol_to_ensembl_mapping()
    
    print(f"\nReading AnnData file: {anndata_file}")
    
    try:
        adata = ad.read_h5ad(anndata_file)
    except Exception as e:
        raise FileNotFoundError(f"Could not read AnnData file: {e}")
    
    print(f"Loaded AnnData: {adata.shape} (cells × genes)")
    print(f"Gene names (var_names) shape: {len(adata.var_names)}")
    
    # Map gene names (var_names) from symbols to Ensembl IDs
    print("Mapping gene symbols to Ensembl IDs...")
    mapped_count = 0
    unmapped_symbols = set()
    
    # Store original gene names in var metadata if not already there
    if 'gene_symbol' not in adata.var.columns:
        adata.var['gene_symbol'] = adata.var_names.copy()
    
    # Create new var_names with Ensembl IDs (with progress bar)
    new_var_names = []
    for symbol in tqdm(adata.var_names, desc="Mapping genes", unit="gene"):
        if symbol in mapping:
            new_var_names.append(mapping[symbol])
            mapped_count += 1
        else:
            new_var_names.append(symbol)  # Keep original if not in mapping
            unmapped_symbols.add(symbol)
    
    # Update var_names
    adata.var_names = new_var_names
    
    print(f"  Mapped {mapped_count} symbols to Ensembl IDs")
    if len(unmapped_symbols) > 0:
        print(f"  {len(unmapped_symbols)} symbols not found in mapping (keeping original)")
        if len(unmapped_symbols) <= 20:
            print(f"  Unmapped symbols: {sorted(list(unmapped_symbols))[:20]}")
    
    # Save to file
    print(f"Writing to: {output_file}")
    adata.write_h5ad(output_file)
    
    print(f"Successfully updated {output_file}")
    
    return adata


def main():
    """Main function to update GRN, matrix, and AnnData files."""
    parser = argparse.ArgumentParser(
        description='Convert gene symbols to Ensembl IDs in GRN, matrix, and AnnData files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '-region',
        type=str,
        default='V1',
        help='Region to process (e.g., V1, V1, V2, ITG, PFC). Default: V1'
    )
    parser.add_argument(
        '-cell_type',
        type=str,
        default='Astrocytes',
        help='Cell type dirV1tory (e.g., Astrocytes, Microglia). Default: Astrocytes'
    )
    parser.add_argument(
        '-grn_only',
        action='store_true',
        help='Only update GRN file'
    )
    parser.add_argument(
        '-matrix_only',
        action='store_true',
        help='Only update matrix file'
    )
    parser.add_argument(
        '-anndata_only',
        action='store_true',
        help='Only update AnnData file'
    )
    
    args = parser.parse_args()
    
    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
    
    csv_file = os.path.join(data_dir, 'GRN', 'Priors', 'TFs_utoronto.csv')
    grn_dir = os.path.join(data_dir, 'GRN', args.cell_type)
    
    # Construct filenames based on cell_type (Microglia has "_Microglia_" in filename)
    if args.cell_type == "Microglia":
        grn_file = os.path.join(grn_dir, f'{args.region}_Microglia_GRN.tsv')
        matrix_file = os.path.join(grn_dir, f'{args.region}_Microglia_matrix.csv')
        anndata_file = os.path.join(grn_dir, f'{args.region}_Microglia_AnnData.h5ad')
        # Output files
        grn_output = os.path.join(grn_dir, f'{args.region}_Microglia_GRN_ensembl.tsv')
        matrix_output = os.path.join(grn_dir, f'{args.region}_Microglia_matrix_ensembl.csv')
        anndata_output = os.path.join(grn_dir, f'{args.region}_Microglia_AnnData_ensembl.h5ad')
    else:
        grn_file = os.path.join(grn_dir, f'{args.region}_GRN.tsv')
        matrix_file = os.path.join(grn_dir, f'{args.region}_matrix.csv')
        anndata_file = os.path.join(grn_dir, f'{args.region}_AnnData.h5ad')
        # Output files
        grn_output = os.path.join(grn_dir, f'{args.region}_GRN_ensembl.tsv')
        matrix_output = os.path.join(grn_dir, f'{args.region}_matrix_ensembl.csv')
        anndata_output = os.path.join(grn_dir, f'{args.region}_AnnData_ensembl.h5ad')
    
    try:
        # Create mapping
        mapping = create_symbol_to_ensembl_mapping(csv_file)
        
        # Determine what to update
        update_grn = args.grn_only or (not args.matrix_only and not args.anndata_only)
        update_matrix = args.matrix_only or (not args.grn_only and not args.anndata_only)
        update_anndata = args.anndata_only or (not args.grn_only and not args.matrix_only)
        
        # Update GRN file
        if update_grn:
            if os.path.exists(grn_file):
                print("\n" + "="*60)
                print(f"Updating {args.region}_GRN.tsv")
                print("="*60)
                update_grn_file(grn_file, mapping, output_file=grn_output)
            else:
                print(f"\nNote: GRN file not found at {grn_file}, skipping...")
        
        # Update matrix file
        if update_matrix:
            if os.path.exists(matrix_file):
                print("\n" + "="*60)
                print(f"Updating {args.region}_matrix.csv")
                print("="*60)
                update_matrix_file(matrix_file, mapping, output_file=matrix_output)
            else:
                print(f"\nNote: Matrix file not found at {matrix_file}, skipping...")
        
        # Update AnnData file
        if update_anndata:
            if os.path.exists(anndata_file):
                print("\n" + "="*60)
                print(f"Updating {args.region}_AnnData.h5ad")
                print("="*60)
                update_anndata_file(anndata_file, mapping, output_file=anndata_output)
            else:
                print(f"\nNote: AnnData file not found at {anndata_file}, skipping...")
        
        print("\n" + "="*60)
        print("Complete!")
        print("="*60)
        
    except Exception as e:
        print(f"\nError: {e}")
        raise


if __name__ == "__main__":
    main()
