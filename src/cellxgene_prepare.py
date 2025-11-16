import scanpy as sc
import pandas as pd

def prepare_cellxgene_data(raw_data_file: str, output_file: str, cell_annotation_file: str, gene_annotation_file: str):
    """
    Prepare data for CellXGene
    """
    # Read the raw data file
    print ("Reading raw data file...")  
    adata = sc.read_mtx(raw_data_file)
    adata = adata.T # transpose the matrix so that cells are rows and genes are columns
    print ("Raw data file read successfully")
    print ("Reading gene annotation file...")
    adata.var_names = pd.read_csv(gene_annotation_file, header=None)[0]
    print ("Gene annotation file read successfully")
    print ("Reading cell annotation file...")
    adata.obs_names = pd.read_csv(cell_annotation_file, header=None)[0]
    print ("Cell annotation file read successfully")
    print ("Writing CellXGene data to H5AD file...")
    adata.write_h5ad(output_file)
    print ("CellXGene data written to H5AD file successfully")

    # Step 1: Calculate QC metrics and store in the anndata object
    print ("Calculating QC metrics...")
    sc.pp.calculate_qc_metrics(adata, inplace=True) 
    print ("QC metrics calculated successfully")

    # Step 2: Normalize with a very vanilla recipe
    print ("Normalizing data...")
    normalized_data = sc.pp.recipe_seurat(adata, copy=True)
    print ("Data normalized successfully")

    # Step 3: Do some basic preprocessing to run PCA and compute the neighbor graph  
    print ("Running PCA...")
    sc.pp.pca(normalized_data)
    print ("PCA run successfully")
    print ("Computing neighbor graph...")
    sc.pp.neighbors(normalized_data)
    print ("Neighbor graph computed successfully")
    print ("Infer clusters with the Louvain algorithm...")
    sc.tl.louvain(normalized_data)
    print ("Clusters inferred successfully")
    print ("Computing tsne and umap embeddings...")
    sc.tl.umap(normalized_data)
    print ("Tsne and umap embeddings computed successfully")
    print ("Writing to output file...")
    sc.write(output_file, normalized_data)
    print ("Output file written successfully")


if __name__ == "__main__":
    prepare_cellxgene_data(raw_data_file=r"i:\sf2026\data\2025-10-22_Astrocytes_EC_matrix.mtx",
         output_file=r"i:\sf2026\data\2025-10-22_Astrocytes_EC_cellxgene_processed.h5ad",
         cell_annotation_file=r"i:\sf2026\data\2025-10-22_Astrocytes_EC_cell_annotation.txt",
         gene_annotation_file=r"i:\sf2026\data\2025-10-22_Astrocytes_EC_row_annotation.txt")
    