#!env python3
import sys
import numpy as np
import pandas as pd
import scanpy as sc

data = sc.read(sys.argv[1],cache=False)
result_file = sys.argv[2]
adata = data.T
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=400)
sc.pp.filter_genes(adata, min_cells=3)
mito_genes = adata.var_names.str.startswith('mt-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
adata.obs['n_counts'] = adata.X.sum(axis=1)
adata = adata[adata.obs['n_genes'] > 400, :]
#adata = adata[adata.obs['percent_mito'] < 0.3, :]
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata)
#sc.pl.highly_variable_genes(adata)
sc.pp.regress_out(adata, ['n_counts', 'n_genes'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.louvain(adata)
sc.tl.umap(adata)
sc.tl.tsne(adata)
adata.write(result_file)
