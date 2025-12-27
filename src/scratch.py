from utils import get_DEGs

gene_list = get_DEGs("V2", top_k=500, cell_type="Microglia")
print(len(gene_list))