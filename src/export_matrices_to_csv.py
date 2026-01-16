"""
Script to export expression matrices for each astrocyte region to CSV files.

For each region, creates a CSV file using pseudobulk aggregation (100 cells per bin):
- Genes as rows
- Pseudobulk bins as columns
- Comma-separated values
- Gene names in first column
- No header row
"""

import os
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.sparse import csr_matrix
from tqdm import tqdm
from utils import (
    read_mtx_file, 
    get_region_file_paths, 
    create_anndata_object,
    read_excel_columns,
    aggregate_cells_into_pseudobulk,
    filter_anndata_object
)
from config import REGIONS, DATA_DIR, REGION_TO_TAB


def export_region_matrix(region, output_dir=None, date_prefix="2025-11-16", 
                        target_cells_per_bin=100, seed=100, cell_type="Astrocytes"):
    """
    Export expression matrix for a region to a text file using pseudobulk aggregation.
    
    Parameters:
    -----------
    region : str
        Region name (EC, ITG, PFC, V1, V2)
    output_dir : str
        Output directory for text files
    date_prefix : str
        Date prefix for file naming
    target_cells_per_bin : int
        Target number of cells per bin for pseudobulk aggregation (default: 100)
    seed : int
        Random seed for cell shuffling (default: 100)
    cell_type : str
        Cell type (e.g., 'Astrocytes', 'Microglia'). Default: 'Astrocytes'
    """
    print(f"\n{'='*60}")
    print(f"Processing region: {region}")
    print(f"Cell type: {cell_type}")
    print(f"{'='*60}")
    
    # Set output directory based on cell_type if not provided
    if output_dir is None:
        output_dir = os.path.join("data", "GRN", cell_type)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get file paths based on cell_type
    if cell_type == "Astrocytes":
        base_prefix = f"{date_prefix}_Astrocytes_{region}"
    elif cell_type == "Microglia":
        # Update date prefix for Microglia if needed
        base_prefix = f"2025-12-27_Microglia_{region}"
    else:
        base_prefix = f"{date_prefix}_{cell_type}_{region}"
    
    mtx_path, row_annotation_path, col_annotation_path = get_region_file_paths(
        region,
        data_dir=DATA_DIR,
        base_prefix=base_prefix
    )
    
    # Get metadata path based on cell_type
    if cell_type == "Astrocytes":
        metadata_path = os.path.join(DATA_DIR, f"{date_prefix}_Astrocytes_metadata.xlsx")
    elif cell_type == "Microglia":
        metadata_path = os.path.join(DATA_DIR, f"2025-12-27_Microglia_{region}_metadata.csv")
    else:
        metadata_path = os.path.join(DATA_DIR, f"{date_prefix}_{cell_type}_metadata.xlsx")
    
    # Get tab index for this region (only for Astrocytes)
    if cell_type == "Microglia":
        tab_index = None
    else:
        tab_index = REGION_TO_TAB.get(region)
    # Load metadata
    print(f"Loading metadata from: {metadata_path}")
    if cell_type == "Microglia":
        # Microglia uses CSV format
        metadata = pd.read_csv(metadata_path)
        # Rename columns with spaces to use underscores
        rename_dict = {}
        if 'Pathology Stage' in metadata.columns:
            rename_dict['Pathology Stage'] = 'Path..Group.'
        if 'Donor ID' in metadata.columns:
            rename_dict['Donor ID'] = 'Donor_ID'  # Use Donor_ID for consistency
        if rename_dict:
            metadata = metadata.rename(columns=rename_dict)
    else:
        # Astrocytes uses Excel format
        metadata = read_excel_columns(
            metadata_path,
            columns=['cell_annotation', 'Path..Group.', 'SampleName', 'percent.mito', 'RIN', 'Total.Genes.Detected'],
            sheet_name=tab_index
        )
    print(f"Loaded metadata: {metadata.shape}")
    
    # Read the MTX file (genes × cells format)
    print(f"Loading matrix from: {mtx_path}")
    matrix, gene_names, cell_names = read_mtx_file(
        mtx_path=str(mtx_path),
        row_annotation_path=str(row_annotation_path),
        col_annotation_path=str(col_annotation_path),
        transpose=False  # Keep as genes × cells
    )
    print(f"Loaded matrix shape: {matrix.shape} (genes × cells)")
    
    # Filter and align metadata to match cell_names
    if cell_type == "Microglia":
        # For Microglia, align by index if cell_names match metadata index
        # Or create a simple DataFrame with the metadata columns
        if cell_names and len(cell_names) == len(metadata):
            # Assume cell_names align with metadata rows
            filtered_metadata = metadata.copy()
            filtered_metadata.index = cell_names
        else:
            # Create empty DataFrame with same number of cells
            filtered_metadata = pd.DataFrame(index=cell_names if cell_names else range(matrix.shape[1]))
            # Try to add available columns from metadata
            for col in metadata.columns:
                if len(metadata) == len(filtered_metadata):
                    filtered_metadata[col] = metadata[col].values
    elif cell_names and 'cell_annotation' in metadata.columns:
        filtered_metadata = metadata[metadata['cell_annotation'].isin(cell_names)].copy()
        filtered_metadata = filtered_metadata.set_index('cell_annotation')
        filtered_metadata = filtered_metadata.reindex(cell_names)
    else:
        filtered_metadata = pd.DataFrame(index=range(matrix.shape[1]))
    
    # Create AnnData object
    print("Creating AnnData object...")
    adata = create_anndata_object(
        matrix=matrix,
        gene_names=gene_names,
        cell_names=cell_names,
        transpose=True,  # Matrix is genes × cells, transpose for AnnData (cells × genes)
        obs=filtered_metadata,
    )

    if cell_type == "Microglia":
        pass
    else:
        adata = filter_anndata_object(adata, mito_max=0.15)
    
    if adata is None:
        raise ImportError("Failed to create AnnData object. Make sure anndata is installed.")

    
    print(f"AnnData object created: {adata.shape} (cells × genes)")
    
    # Aggregate cells into pseudobulk bins
    print(f"\nAggregating cells into pseudobulk bins ({target_cells_per_bin} cells per bin)...")
    adata_pseudobulk = aggregate_cells_into_pseudobulk(
        adata, 
        target_cells_per_bin=target_cells_per_bin,
        filter_patients_cell_threshold=80,
        seed=seed,
        patient_col='SampleName' if cell_type != "Microglia" else "Donor_ID"
    )
    
    print(f"After aggregation: {adata_pseudobulk.shape} (bins × genes)")
    
    # Save pseudobulk AnnData object
    anndata_output_filename = f"{region}_{cell_type}_AnnData.h5ad" if cell_type == "Microglia" else f"{region}_AnnData.h5ad"
    anndata_output_path = os.path.join(output_dir, anndata_output_filename)
    print(f"Saving pseudobulk AnnData object to: {anndata_output_path}")
    adata_pseudobulk.write_h5ad(anndata_output_path)
    print(f"Successfully saved AnnData object")
    
    # Convert to genes × bins format for export
    expression = adata_pseudobulk.X
    expression = expression.toarray()
    expression = expression.T  # Transpose to genes × bins
    
    # Get gene and bin names
    gene_names = list(adata_pseudobulk.var_names)
    bin_names = list(adata_pseudobulk.obs_names)
    
    # Extract patient names from bin names or from obs metadata
    # Bin names are in format: "{patient_name}_bin{bin_idx+1}"
    if cell_type == "Microglia":
        # For Microglia, use Donor_ID
        if 'Donor_ID' in adata_pseudobulk.obs.columns:
            patient_names = adata_pseudobulk.obs['Donor_ID'].values.tolist()
            print(f"Using patient names from adata_pseudobulk.obs['Donor_ID']")
        else:
            # Fallback: extract from bin names
            patient_names = []
            for bin_name in bin_names:
                if '_bin' in bin_name:
                    patient_name = bin_name.rsplit('_bin', 1)[0]
                    patient_names.append(patient_name)
                else:
                    patient_names.append(bin_name)
            print(f"Extracted patient names from bin names")
    elif 'SampleName' in adata_pseudobulk.obs.columns:
        # Use SampleName from metadata if available (more reliable)
        patient_names = adata_pseudobulk.obs['SampleName'].values.tolist()
        print(f"Using patient names from adata_pseudobulk.obs['SampleName']")
    else:
        # Extract from bin names (format: "patient_bin1")
        patient_names = []
        for bin_name in bin_names:
            if '_bin' in bin_name:
                patient_name = bin_name.rsplit('_bin', 1)[0]  # Split from right, take first part
                patient_names.append(patient_name)
            else:
                patient_names.append(bin_name)  # Fallback if format is different
        print(f"Extracted patient names from bin names")
    
    # Create DataFrame
    df = pd.DataFrame(
        expression,
        index=gene_names,
        columns=bin_names
    )
    
    # Create output filename (include cell_type)
    if cell_type == "Microglia":
        output_filename = f"{region}_Microglia_matrix.csv"
    else:
        output_filename = f"{region}_matrix.csv"
    patient_mapping_filename = f"{region}_{cell_type}_patient_mapping.csv"
    output_path = os.path.join(output_dir, output_filename)
    
    # Create patient mapping file
    patient_mapping_path = os.path.join(output_dir, patient_mapping_filename)
    
    # Save patient to bin mapping
    patient_mapping_df = pd.DataFrame({
        'bin_name': bin_names,
        'patient_name': patient_names
    })
    # Also add condition if available
    if 'Path..Group.' in adata_pseudobulk.obs.columns:
        patient_mapping_df['condition'] = adata_pseudobulk.obs['Path..Group.'].values
        print(f"  Added condition information to mapping")
    patient_mapping_df.to_csv(patient_mapping_path, index=False)
    print(f"Saved patient mapping to: {patient_mapping_path}")
    print(f"  {len(set(patient_names))} unique patients across {len(bin_names)} bins")
        
    # Save to CSV file (comma-separated) with progress bar
    print(f"Saving to: {output_path}")
    n_genes, n_bins = df.shape
    
    # Write file with progress bar (NO header row)
    with open(output_path, 'w') as f:
        # Write data rows with progress bar (no header)
        matrix_values = df.values
        gene_names_list = df.index.tolist()
        
        for i in tqdm(range(n_genes), desc="Writing rows", unit="genes"):
            gene_name = gene_names_list[i]
            row_data = ','.join(matrix_values[i, :].astype(str))
            f.write(f"{gene_name},{row_data}\n")
    
    print(f"Successfully exported {df.shape[0]} genes × {df.shape[1]} bins")
    
    return output_path


def main():
    """Export matrices for astrocyte regions using pseudobulk aggregation."""
    parser = argparse.ArgumentParser(
        description='Export expression matrices for astrocyte regions to CSV files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '-regions',
        type=str,
        nargs='+',
        choices=REGIONS + ['CrossRegion'],
        default=REGIONS + ['CrossRegion'],
        help=f'Region(s) to process (default: all regions). Choices: {", ".join(REGIONS)}'
    )
    parser.add_argument(
        '-cell_type',
        type=str,
        default="Astrocytes",
        help='Cell type to filter by (e.g., "Astrocytes", "Microglia"). If not specified, no filtering is applied.'
    )
    
    args = parser.parse_args()
    
    # Get regions to process
    regions_to_process = args.regions if isinstance(args.regions, list) else [args.regions]
    
    print("="*60)
    print("Using pseudobulk aggregation (100 cells per bin)")
    print(f"Processing regions: {', '.join(regions_to_process)}")
    print("="*60)
        
    # Process each region
    for region in regions_to_process:
        try:
            export_region_matrix(region, target_cells_per_bin=100, seed=100, cell_type=args.cell_type)
        except Exception as e:
            print(f"\nError processing region {region}: {e}")
            exit()
    
    print(f"\n{'='*60}")
    print("Export complete!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

