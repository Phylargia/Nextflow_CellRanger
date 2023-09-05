# This is a simple python script that utilises memory_profiler to capture the memory usage of Scanpy Pre-processing (scRNAseq Analysis)
import numpy as np
import pandas as pd
import scanpy as sc
from memory_profiler import profile

@profile
def load_and_process_data():
    data_path = 'C:/Users/Kenod/Pictures/count/filtered_feature_bc_matrix'
    adata = sc.read_10x_mtx(data_path, var_names='gene_symbols', cache=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.X[0,:].sum()
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
    sc.tl.umap(adata)
    sc.tl.tsne(adata)
    sc.tl.leiden(adata, resolution=0.3)
    
if __name__ == "__main__":
    load_and_process_data()
       
