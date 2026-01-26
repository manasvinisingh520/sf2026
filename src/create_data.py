"""Scale all gene embeddings according to cell expression and save in pickle file.
For Astrocytes only, with 10-fold cross-validation split and target variables.
"""

import pandas as pd
import numpy as np
import os
import pickle
import gc
import anndata as ad
import glob
from tqdm import tqdm
from sklearn.model_selection import KFold
from utils import get_region_file_paths, read_excel_columns, read_mtx_file, scale_embeddings, create_anndata_object
from config import REGION_TO_TAB, DATA_DIR


def load_gene_embeddings(embeddings_path, gene_list: list[str]):
    """Load gene embeddings from TSV/CSV file."""
    """ gene_list: is list of all genes in the region data"""
    # Read CSV file with header (first row contains dimension indices)
    # File is comma-separated: first column is gene IDs, remaining columns are embedding dimensions
    # Try comma-separated first, then tab-separated as fallback
    df = pd.read_csv(embeddings_path, sep=',', header=0)
    
    # Create a dictionary with gene_ids as key and embedding_data as value
    gene_ids = df.iloc[:, 0].astype(str).values  # All rows, first column = gene IDs
    embedding_data = df.iloc[:, 1:].values.astype(float)  # All rows, columns 1 onwards = embeddings
    gene_embedding_dict = {gene_id: embedding for gene_id, embedding in zip(gene_ids, embedding_data)}

    # Convert ensembl ids to gene symbols (reverse the mapping)
    with open("data/genes_2_ensembl_ids.pkl", "rb") as f:
        symbol_to_ensembl = pickle.load(f)

    
    # for each gene in gene_list, if the gene is in gene_embedding_dict, keep the embedding
    not_found_genes = []
    embedding_data = {}
    for gene in gene_list:
        gene_ensembl_id = symbol_to_ensembl.get(gene)
        if gene_ensembl_id in gene_embedding_dict:
            embedding_data[gene] = gene_embedding_dict[gene_ensembl_id]
        else:
            not_found_genes.append(gene)
    print(f"Not found genes: {len(not_found_genes)}")
    print(f"Not found genes: {not_found_genes[:10]}")
    print(f"Found genes: {len(embedding_data)}")
    return embedding_data

def find_anndata_file(region, cell_type='Astrocytes'):
    """Find AnnData file for a given region and cell type (perCell only)."""
    # Use perCell AnnData file specifically
    if cell_type == "Astrocytes":
        path = os.path.join('data', 'model_data', f'{region}_AnnData_perCell.h5ad')
    elif cell_type == "Microglia":
        path = os.path.join('data', 'model_data', f'{region}_Microglia_AnnData_perCell.h5ad')
    else:
        raise ValueError(f"Invalid cell type: {cell_type}")
    
    if os.path.exists(path):
        return path
    
    return None

def create_anndata(region):
    
    base_prefix = f"2025-11-16_Astrocytes_{region}"
    
    mtx_path, row_annotation_path, col_annotation_path = get_region_file_paths(
        region,
        data_dir=DATA_DIR,
        base_prefix=base_prefix
    )
    # Get metadata path
    metadata_path = os.path.join(DATA_DIR, f"2025-11-16_Astrocytes_metadata.xlsx")
    
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
    
    # Align metadata with cell_names (set cell_annotation as index)
    metadata = metadata.set_index('cell_annotation')
    metadata = metadata.reindex(cell_names)  # Align with cell_names order
    
    # Filter cells with unknown Thal category BEFORE creating AnnData
    valid_cell_mask = (metadata['Thal'] != 'unk').fillna(False)
    valid_cell_names = metadata[valid_cell_mask].index.tolist()
    valid_metadata = metadata[valid_cell_mask].copy()
    
    # Filter cell_names and matrix to match valid cells
    cell_indices_to_keep = [i for i, name in enumerate(cell_names) if name in valid_cell_names]
    filtered_cell_names = [cell_names[i] for i in cell_indices_to_keep]
    filtered_matrix = matrix[:, cell_indices_to_keep]  # genes × cells
    
    # Create AnnData object (transpose=True because matrix is genes × cells)
    adata = create_anndata_object(
        matrix=filtered_matrix,
        gene_names=gene_names,
        cell_names=filtered_cell_names,
        obs=valid_metadata,
        transpose=True  # Matrix is genes × cells, transpose for AnnData
    )
    
    return adata

def main():
    print("="*80)
    print("Scaling All Gene Embeddings by Cell Expression (Astrocytes Only)")
    print("10-Fold Cross-Validation Split with Target Variables")
    print("="*80)

    # Target dimensions for PCA reduction
    target_dims = 16
    target_cell_type = "Astrocytes"
    n_folds = 8
    
    # Find all gene embedding files
    embeddings_dir = os.path.join("GREmLN", "embeddings")
    
    # Filter to only Astrocyte files
    if target_cell_type == "Astrocytes":
        regions = ["EC", "ITG", "PFC", "V2", "V1"]
    elif target_cell_type == "Microglia":
        regions = ["CrossRegion", "EC", "ITG", "PFC", "V2", "V1"]
    else:
        raise ValueError(f"Invalid cell type: {target_cell_type}")
    
    for region in regions:
        print(f"Processing {region}")
        if target_cell_type == "Astrocytes":
            emb_path = os.path.join(embeddings_dir, f"{region}_gene_embeddings.tsv")
        elif target_cell_type == "Microglia":
            emb_path = os.path.join(embeddings_dir, f"{region}_Microglia_gene_embeddings.tsv")

        # load the region data
        adata = create_anndata(region)
        # Load gene embeddings
        embeddings_dict = load_gene_embeddings(emb_path, gene_list=adata.var_names)

        # Filter the adata for which not found genes are not in the adata.var_names
        adata = adata[:, adata.var_names.isin(embeddings_dict.keys())].copy()
        print(f"Filtered adata shape: {adata.shape}")

         # Get expression matrix with cells as rows and genes as columns
        expression_matrix = np.array(adata.X.toarray(), dtype=np.float16)

        # create embeddings matrix from embeddings_dict
        embeddings_matrix = np.array([embeddings_dict[gene] for gene in adata.var_names], dtype=np.float16)
        print(f"Embeddings matrix shape: {embeddings_matrix.shape}")
        
        # Verify shapes match (should always match since we filtered to common genes)
        if embeddings_matrix.shape[0] != expression_matrix.shape[1]:
            print(f"ERROR: Shape mismatch! Embeddings: {embeddings_matrix.shape[0]} genes, Expression: {expression_matrix.shape[1]} genes.")
            quit()
        
        # Extract ptau_info and thal_info from obs
        ptau_info = adata.obs['Ptau.Total.Tau..A.U..']
        thal_info = adata.obs['Thal']
        
        # expression_matrix = d, embedding = A
        print("Scaling embeddings...")
        region_scaled_embeddings = scale_embeddings(embeddings_matrix, expression_matrix, n_components=target_dims)
        del embeddings_matrix, expression_matrix
        # region_scaled_embeddings shape: (n_cells, n_genes, n_dims)

        # Collect all patient embeddings into lists for KFold splitting
        all_patient_embeddings = []
        all_cell_names = []
        all_ptau = []
        all_thal = []
        unique_patients = adata.obs['SampleName'].unique()
        print(f"Number of unique patients: {len(unique_patients)}")

        # verify that the number of cells in the embeddings is the same as the number of cells in the adata.obs
        if region_scaled_embeddings.shape[0] != len(adata.obs):
            print(f"ERROR: Mismatch! Embeddings: {region_scaled_embeddings.shape[0]} cells, AnnData obs: {len(adata.obs)} cells")
            quit()

        for patient_name in unique_patients:
            patient_mask = adata.obs['SampleName'] == patient_name
            patient_embeddings = region_scaled_embeddings[patient_mask].astype(np.float16)
            patient_cell_names = adata.obs.index[patient_mask].tolist()
            all_patient_embeddings.append(patient_embeddings)
            all_cell_names.append(patient_cell_names)
            # ptau and thal have only one value per patient, so take the first (or unique) value
            all_ptau.append(ptau_info[patient_mask].iloc[0])
            all_thal.append(thal_info[patient_mask].iloc[0])
        
        # Split patients into 8 folds (not cells, but patients)
        n_patients = len(unique_patients)
        kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
        
        embeddings_folds = []
        patient_info_folds = []
        ptau_folds = []
        thal_folds = []
        
        print(f"Splitting {n_patients} patients into {n_folds} folds...")
        for fold_idx, (_, fold_patient_indices) in enumerate(kf.split(range(n_patients))):
            # Concatenate embeddings from patients in this fold
            fold_embeddings_list = [all_patient_embeddings[i] for i in fold_patient_indices]
            fold_embeddings = np.concatenate(fold_embeddings_list, axis=0).astype(np.float16)
            
            # Get ptau and thal for patients in this fold (one value per patient)
            fold_ptau = [all_ptau[i] for i in fold_patient_indices]
            fold_thal = [all_thal[i] for i in fold_patient_indices]
            fold_patient_info = {unique_patients[i]: all_cell_names[i] for i in fold_patient_indices}
            
            embeddings_folds.append(fold_embeddings)
            patient_info_folds.append(fold_patient_info)
            ptau_folds.append(fold_ptau)
            thal_folds.append(fold_thal)
            
            n_fold_cells = fold_embeddings.shape[0]
            n_fold_patients = len(fold_patient_indices)
            print(f"  Fold {fold_idx + 1}: {n_fold_patients} patients, {n_fold_cells} cells")

        # Save folds
        emb_output_dir = os.path.join("data", "model_data", f"{region}_data_for_model_training.pkl")
        with open(emb_output_dir, "wb") as f:
            pickle.dump((embeddings_folds, patient_info_folds, ptau_folds, thal_folds), f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Saved {n_folds} folds of rescaled embeddings to {emb_output_dir}")
        print(f"  Each fold contains: embeddings (3D array), patient_info (dict: patient_name -> cell_names), ptau, thal")
        del embeddings_folds, patient_info_folds, ptau_folds, thal_folds, region_scaled_embeddings
        del adata, fold_embeddings, fold_patient_info, fold_ptau, fold_thal
        del all_patient_embeddings, all_cell_names, all_ptau, all_thal
        del unique_patients
        gc.collect()
        
        # Try to get memory info if psutil is available
        try:
            import psutil
            process = psutil.Process()
            memory_mb = process.memory_info().rss / (1024 * 1024)
            print(f"Memory after saving: {memory_mb:.2f} MB")
        except ImportError:
            print("Memory info not available (psutil not installed)")
  
if __name__ == "__main__":
    main()
