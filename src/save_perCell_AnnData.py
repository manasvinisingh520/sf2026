"""
Script to create, filter, and save AnnData object for a given region.
"""

import argparse
import os
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from utils import (
    read_mtx_file,
    create_anndata_object,
    read_excel_columns,
    get_region_file_paths,
    gene_symbols_to_ensembl
)
from config import REGIONS, DATA_DIR, REGION_TO_TAB, METADATA_DATE_PREFIX


def create_and_save_anndata(region, output_dir="data\model_data", date_prefix="2025-11-16"):

    print(f"\n{'='*60}")
    print(f"Processing region: {region}")
    print(f"{'='*60}")
    
    # Get file paths
    base_prefix = f"{date_prefix}_Astrocytes_{region}"
    
    mtx_path, row_annotation_path, col_annotation_path = get_region_file_paths(
        region,
        data_dir=DATA_DIR,
        base_prefix=base_prefix
    )
    # Get metadata path
    metadata_path = os.path.join(DATA_DIR, f"{date_prefix}_Astrocytes_metadata.xlsx")
    
    # Get tab index for this region
    tab_index = REGION_TO_TAB.get(region)
    
    # Load metadata
    metadata = read_excel_columns(
        metadata_path,
        columns=['cell_annotation', 'Path..Group.', 'SampleName', 'percent.mito', 'RIN', 'Total.Genes.Detected', 'Thal', 'Ptau.Total.Tau..A.U..'],
        sheet_name=tab_index
    )
    print(f"Loaded metadata: {metadata.shape}")
    
    # Read the MTX file (genes × cells format)
    matrix, gene_names, cell_names = read_mtx_file(
        mtx_path=str(mtx_path),
        row_annotation_path=str(row_annotation_path),
        col_annotation_path=str(col_annotation_path),
        transpose=False  # Keep as genes × cells
    )
    print(f"Loaded matrix shape: {matrix.shape} (genes × cells)")
    
    # Convert gene symbols to Ensembl IDs
    print(f"\nConverting gene symbols to Ensembl IDs...")
    csv_file = os.path.join(DATA_DIR, "GRN", "Priors", "TFs_utoronto.csv")
    if not os.path.exists(csv_file):
        # Try alternative path
        csv_file = os.path.join("data", "GRN", "Priors", "TFs_utoronto.csv")
    
    # Store original gene names for metadata
    original_gene_names = gene_names.copy()
    
    # Use gene_symbols_to_ensembl from utils (handles CSV + REST API fallback)
    if os.path.exists(csv_file):
        print(f"Using mapping file: {csv_file}")
        mapping = gene_symbols_to_ensembl(gene_names, csv_file=csv_file, verbose=True)
    else:
        print(f"Warning: Mapping file not found at {csv_file}")
        print("  Using Ensembl REST API directly...")
        mapping = gene_symbols_to_ensembl(gene_names, csv_file=None, verbose=True)
    
    # Map gene names to Ensembl IDs
    mapped_count = 0
    unmapped_symbols = set()
    ensembl_gene_names = []
    
    for symbol in tqdm(gene_names, desc="Applying mapping", unit="gene"):
        if symbol in mapping:
            ensembl_gene_names.append(mapping[symbol])
            mapped_count += 1
        elif symbol.startswith('ENSG'):
            # Already an Ensembl ID
            ensembl_gene_names.append(symbol)
            mapped_count += 1
        else:
            # Keep original if not in mapping
            ensembl_gene_names.append(symbol)
            unmapped_symbols.add(symbol)
    
    print(f"\nMapping results:")
    print(f"  Mapped {mapped_count}/{len(gene_names)} genes to Ensembl IDs ({100*mapped_count/len(gene_names):.1f}%)")
    if len(unmapped_symbols) > 0:
        print(f"  {len(unmapped_symbols)} genes not found in mapping (keeping original)")
        if len(unmapped_symbols) <= 20:
            print(f"  Sample unmapped genes: {sorted(list(unmapped_symbols))[:20]}")
    
    # Use Ensembl IDs as gene names
    gene_names = ensembl_gene_names
    
    # Filter and align metadata to match cell_names
    filtered_metadata = metadata[metadata['cell_annotation'].isin(cell_names)].copy()
    filtered_metadata = filtered_metadata.set_index('cell_annotation')
    filtered_metadata = filtered_metadata.reindex(cell_names)
    
    # Create AnnData object
    print(f"\nCreating AnnData object...")
    adata = create_anndata_object(
        matrix=matrix,
        gene_names=gene_names,
        cell_names=cell_names,
        transpose=True,  # Matrix is genes × cells, transpose for AnnData (cells × genes)
        obs=filtered_metadata,
    )
    
    if adata is None:
        raise ImportError("Failed to create AnnData object. Make sure anndata is installed.")
    
    # Store original gene symbols in var metadata if they were converted
    if 'gene_symbol' not in adata.var.columns and len(original_gene_names) == len(adata.var_names):
        adata.var['gene_symbol'] = original_gene_names
    
    print(f"AnnData object created: {adata.shape} (cells × genes)")
    
    # Create output filename
    output_filename = f"{region}_AnnData_perCell.h5ad"
    output_path = os.path.join(output_dir, output_filename)
    
    # Save as h5ad
    adata.write_h5ad(output_path)
    print(f"Successfully saved AnnData object to {output_path}")
    
    return output_path


def main():
    parser = argparse.ArgumentParser(
        description='Create, filter, and save AnnData object for a region',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '-region',
        type=str,
        choices=REGIONS,
        required=True,
        help=f'Region name. Choices: {", ".join(REGIONS)}'
    )
    args = parser.parse_args()

    try:
        create_and_save_anndata(
            region=args.region
        )
    except Exception as e:
        print(f"\nError: {e}")
        raise


if __name__ == "__main__":
    main()
