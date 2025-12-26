"""
Script to create tables showing cell counts by condition and patient for each region.

For each of the five regions (EC, ITG, PFC, V2, V1), this script creates a table with:
- Rows: Conditions (from Path..Group. column)
- Columns: Patients (from SampleName column)
- Values: Number of cells
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from utils import read_mtx_file, create_anndata_object, read_excel_columns, filter_cells, get_region_file_paths
from config import REGIONS, REGION_TO_TAB, DATA_DIR, METADATA_DATE_PREFIX

# Set style for better-looking plots
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300


def load_region_data(region, base_dir="data", date_prefix="2025-11-16"):
    # Get file paths using utility function
    mtx_path, row_annotation_path, col_annotation_path = get_region_file_paths(
        region, 
        data_dir=base_dir, 
        base_prefix=f"{date_prefix}_Astrocytes_{{region}}"
    )
    metadata_path = os.path.join(base_dir, f"{date_prefix}_Astrocytes_metadata.xlsx")
    
    # Get tab index for this region
    tab_index = REGION_TO_TAB.get(region)
    if tab_index is None:
        raise ValueError(f"Region '{region}' not found in REGION_TO_TAB. Available regions: {list(REGION_TO_TAB.keys())}")
    
    # Load metadata (include UMI and mito columns for filtering)
    metadata = read_excel_columns(
        metadata_path,
        columns=['cell_annotation', 'Path..Group.', 'SampleName', 'Median.UMI.Counts.per.Cell', 'percent.mito'],
        sheet_name=tab_index
    )
    print(f"Loaded metadata: {metadata.shape}")
    
    # Read the MTX file
    matrix, gene_names, cell_names = read_mtx_file(
        mtx_path=str(mtx_path),
        row_annotation_path=str(row_annotation_path),
        col_annotation_path=str(col_annotation_path),
        transpose=False  # Matrix will be genes × cells
    )
    print(f"Loaded matrix shape: {matrix.shape} (genes × cells)")
    
    # Use filter_cells from utils to filter by mitochondrial percentage
    filtered_matrix, filtered_cell_names, filtered_metadata = filter_cells(
        matrix, cell_names, metadata, mito_max=0.15
    )
    # Update matrix and cell_names for AnnData creation
    matrix = filtered_matrix
    cell_names = filtered_cell_names
    
    # Create AnnData object
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
    
    return adata


def create_cell_count_table(adata, condition_col='Path..Group.', patient_col='SampleName'):
    # Check that required columns exist
    if condition_col not in adata.obs.columns:
        raise ValueError(f"Condition column '{condition_col}' not found in adata.obs. Available columns: {list(adata.obs.columns)}")
    if patient_col not in adata.obs.columns:
        raise ValueError(f"Patient column '{patient_col}' not found in adata.obs. Available columns: {list(adata.obs.columns)}")
    
    # Create a DataFrame with condition and patient columns
    df = pd.DataFrame({
        'condition': adata.obs[condition_col],
        'patient': adata.obs[patient_col]
    })
    
    # Remove rows with missing values
    df = df.dropna()
    # Create pivot table
    table = pd.crosstab(df['condition'], df['patient'], margins=False)
    # Sort rows and columns for better readability
    table = table.sort_index(axis=0)  # Sort conditions
    table = table.sort_index(axis=1)  # Sort patients
    
    return table


def visualize_table(table, region):
    # Create figure with appropriate size
    fig, ax = plt.subplots(figsize=(max(12, len(table.columns) * 0.8), max(8, len(table.index) * 0.6)))
    
    # Create mask for zeros (they will appear white/neutral)
    zero_mask = (table == 0)
    
    # Create heatmap with reversed colormap (red for low values, blue for high values)
    sns.heatmap(
        table,
        annot=True,  # Show cell counts in each box
        fmt='.0f',    # Format as integers (no decimals) but allow float type
        cmap='RdBu_r',  # Reversed Red-Blue: red for low, blue for high
        cbar_kws={'label': 'Number of Cells'},
        ax=ax,
        linewidths=0.5,
        linecolor='gray',
        mask=zero_mask  # Mask zeros so they appear white
    )
    
    # Set labels and title
    ax.set_xlabel('Patient', fontsize=12, fontweight='bold')
    ax.set_ylabel('Condition', fontsize=12, fontweight='bold')
    ax.set_title(f'Cell Counts by Condition and Patient - Region: {region}', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Rotate labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Show the plot
    plt.show()


def main():
    # Configuration
    from config import DATA_DIR, METADATA_DATE_PREFIX
    base_dir = DATA_DIR
    date_prefix = METADATA_DATE_PREFIX
    
    # Check if base directory exists
    if not os.path.exists(base_dir):
        print(f"Warning: Base directory '{base_dir}' not found. Trying relative path 'data'")
        base_dir = "data"
    
    print("=" * 60)
    print("Loading and processing all regions...")
    print("=" * 60)
    
    # First pass: Load all data and create all tables
    all_tables = {}
    
    for region in REGIONS:
        print(f"\nProcessing region: {region}")
        # Load data for this region
        adata = load_region_data(region, base_dir=base_dir, date_prefix=date_prefix)
        # Create cell count table
        table = create_cell_count_table(adata)
        all_tables[region] = table
    
    print("\n" + "=" * 60)
    print("All data loaded and processed. Displaying visualizations...")
    print("=" * 60)
    
    # Second pass: Display all visualizations
    for region in REGIONS:
        if region in all_tables:
            visualize_table(all_tables[region], region)

    return all_tables


if __name__ == "__main__":
    tables = main()

