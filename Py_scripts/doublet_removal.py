# https://doubletdetection.readthedocs.io/en/latest/tutorial.html
# https://github.com/JonathanShor/DoubletDetection

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scvi
import seaborn as sns
import doubletdetection

# Define path directory to working directory containing input files
data_path = 'C:/Users/Kenod/Pictures/count/filtered_feature_bc_matrix'

# Change working directory and list files
os.chdir(data_path)
files = os.listdir()
print(files)

# Read in data
adata = sc.read_10x_mtx(data_path, var_names='gene_symbols', cache=True)

# Summary
adata.shape
adata

# Doublet Removal
adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_cells=1)

clf = doubletdetection.BoostClassifier(
    n_iters=10,
    clustering_algorithm="louvain",
    standard_scaling=True,
    pseudocount=0.1,
    n_jobs=-1,
)
doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
doublet_score = clf.doublet_score()

adata.obs["doublet"].value_counts()

# Pre-process
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Plot
sc.pl.umap(adata, color=["doublet", "doublet_score"])
