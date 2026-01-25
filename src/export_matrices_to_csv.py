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
import pandas as pd
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
    print(f"\n{'='*60}\nProcessing region: {region} | Cell type: {cell_type}\n{'='*60}")
    
    # Set output directory
    output_dir = output_dir or os.path.join("data", "GRN", cell_type)
    os.makedirs(output_dir, exist_ok=True)
    
    # Configure paths based on cell_type
    is_microglia = cell_type == "Microglia"
    is_astrocytes = cell_type == "Astrocytes"
    microglia_date = "2025-12-27"
    
    base_prefix = f"{microglia_date if is_microglia else date_prefix}_{cell_type}_{region}"
    if is_microglia:
        metadata_path = os.path.join(DATA_DIR, f"{microglia_date}_Microglia_{region}_metadata.csv")
    else:
        metadata_path = os.path.join(DATA_DIR, f"{date_prefix}_Astrocytes_metadata.xlsx")
    
    tab_index = None if is_microglia else REGION_TO_TAB.get(region)
    patient_col = 'Donor_ID' if is_microglia else 'SampleName'
    
    mtx_path, row_annotation_path, col_annotation_path = get_region_file_paths(
        region, data_dir=DATA_DIR, base_prefix=base_prefix
    )
    # Load metadata
    print(f"Loading metadata from: {metadata_path}")
    if is_microglia:
        metadata = pd.read_csv(metadata_path)
        metadata = metadata.rename(columns={'Pathology Stage': 'Path..Group.', 'Donor ID': 'Donor_ID'})
    else:
        metadata = read_excel_columns(
            metadata_path, columns=['cell_annotation', 'Path..Group.', 'SampleName', 'percent.mito', 'RIN', 'Total.Genes.Detected'],
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
    
    # Align metadata to cell_names
    if is_microglia and cell_names and len(cell_names) == len(metadata):
        filtered_metadata = metadata.copy()
        filtered_metadata.index = cell_names
    elif not is_microglia and cell_names and 'cell_annotation' in metadata.columns:
        filtered_metadata = metadata[metadata['cell_annotation'].isin(cell_names)].set_index('cell_annotation').reindex(cell_names)
    else:
        filtered_metadata = pd.DataFrame(index=cell_names if cell_names else range(matrix.shape[1]))
        if len(metadata) == len(filtered_metadata):
            filtered_metadata = pd.concat([filtered_metadata, metadata.reset_index(drop=True)], axis=1)
    
    # Create and filter AnnData object
    adata = create_anndata_object(matrix=matrix, gene_names=gene_names, cell_names=cell_names, transpose=True, obs=filtered_metadata)
    if adata is None:
        raise ImportError("Failed to create AnnData object. Make sure anndata is installed.")
    if is_astrocytes:
        adata = filter_anndata_object(adata, mito_max=0.15)
    print(f"AnnData created: {adata.shape} (cells × genes)")
    
    # Aggregate and save pseudobulk
    print(f"\nAggregating cells into pseudobulk ({target_cells_per_bin} cells/bin)...")
    adata_pseudobulk = aggregate_cells_into_pseudobulk(
        adata, target_cells_per_bin=target_cells_per_bin, filter_patients_cell_threshold=80, seed=seed, patient_col=patient_col
    )
    print(f"After aggregation: {adata_pseudobulk.shape} (bins × genes)")
    
    anndata_path = os.path.join(output_dir, f"{region}_{cell_type}_AnnData.h5ad" if is_microglia else f"{region}_AnnData.h5ad")
    adata_pseudobulk.write_h5ad(anndata_path)
    print(f"Saved AnnData to: {anndata_path}")
    
    # Convert to genes × bins format
    expression = adata_pseudobulk.X.toarray().T
    gene_names = list(adata_pseudobulk.var_names)
    bin_names = list(adata_pseudobulk.obs_names)
    
    # Extract patient names (prefer metadata, fallback to bin names)
    patient_names = (adata_pseudobulk.obs[patient_col].values.tolist() 
                     if patient_col in adata_pseudobulk.obs.columns
                     else [bn.rsplit('_bin', 1)[0] if '_bin' in bn else bn for bn in bin_names])
    
    # Create and save patient mapping
    df = pd.DataFrame(expression, index=gene_names, columns=bin_names)
    patient_mapping_df = pd.DataFrame({'bin_name': bin_names, 'patient_name': patient_names})
    if 'Path..Group.' in adata_pseudobulk.obs.columns:
        patient_mapping_df['condition'] = adata_pseudobulk.obs['Path..Group.'].values
    
    mapping_path = os.path.join(output_dir, f"{region}_{cell_type}_patient_mapping.csv")
    patient_mapping_df.to_csv(mapping_path, index=False)
    print(f"Saved patient mapping: {len(set(patient_names))} patients, {len(bin_names)} bins")
        
    # Save CSV (no header, gene name in first column)
    output_path = os.path.join(output_dir, f"{region}_Microglia_matrix.csv" if is_microglia else f"{region}_matrix.csv")
    print(f"Saving to: {output_path}")
    with open(output_path, 'w') as f:
        for gene_name, row in tqdm(zip(df.index, df.values), total=len(df), desc="Writing genes", unit="genes"):
            f.write(f"{gene_name},{','.join(row.astype(str))}\n")
    print(f"Exported {df.shape[0]} genes × {df.shape[1]} bins")
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

