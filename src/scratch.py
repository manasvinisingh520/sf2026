"""Print the number of patients for each Astrocytes region."""

import pandas as pd
import os
import anndata as ad
from config import REGION_TO_TAB, DATA_DIR
from utils import read_excel_columns

def print_patient_counts():
    """Print the number of unique patients for each Astrocytes region."""
    regions = ["EC", "ITG", "PFC", "V2", "V1"]
    
    print("="*60)
    print("Number of patients per Astrocytes region")
    print("="*60)
    
    for region in regions:
        # Try to load from perCell AnnData file first
        anndata_path = os.path.join('data', 'model_data', f'{region}_AnnData_perCell.h5ad')
        
        if os.path.exists(anndata_path):
            try:
                adata = ad.read_h5ad(anndata_path)
                n_patients = adata.obs['SampleName'].nunique()
                print(f"{region}: {n_patients} patients")
            except Exception as e:
                print(f"{region}: Error loading AnnData - {e}")
        else:
            # Fallback to metadata file
            metadata_path = os.path.join(DATA_DIR, "2025-11-16_Astrocytes_metadata.xlsx")
            tab_index = REGION_TO_TAB.get(region)
            
            if tab_index is not None and os.path.exists(metadata_path):
                try:
                    metadata = read_excel_columns(
                        metadata_path,
                        columns=['cell_annotation', 'SampleName'],
                        sheet_name=tab_index
                    )
                    n_patients = metadata['SampleName'].nunique()
                    print(f"{region}: {n_patients} patients (from metadata)")
                except Exception as e:
                    print(f"{region}: Error loading metadata - {e}")
            else:
                print(f"{region}: No data file found")
    
    print("="*60)

if __name__ == "__main__":
    print_patient_counts()
