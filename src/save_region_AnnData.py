"""
Script to create, filter, and save AnnData object for a given region.
"""

import argparse
import os
import pandas as pd
from pathlib import Path
from utils import (
    read_mtx_file,
    create_anndata_object,
    read_excel_columns,
    get_region_file_paths,
    filter_anndata_object
)
from config import REGIONS, DATA_DIR, REGION_TO_TAB, METADATA_DATE_PREFIX


def create_and_save_anndata(region, output_dir="data\GRN", date_prefix="2025-11-16"):

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
        columns=['cell_annotation', 'Path..Group.', 'SampleName', 'percent.mito', 'RIN', 'Total.Genes.Detected'],
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
    
    print(f"AnnData object created: {adata.shape} (cells × genes)")
    
    # Filter AnnData object
    adata = filter_anndata_object(adata, mito_max=0.15)
    print(f"After filtering: {adata.shape} (cells × genes)")
    
    # Create output filename
    output_filename = f"{region}_AnnData.h5ad"
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
        output_path = create_and_save_anndata(
            region=args.region
        )
    except Exception as e:
        print(f"\nError: {e}")
        raise


if __name__ == "__main__":
    main()
