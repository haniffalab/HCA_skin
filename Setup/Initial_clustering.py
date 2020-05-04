#!/usr/bin/env python
# coding: utf-8

#%% Import
import os.path
import numpy as np
import pandas as pd
import scanpy as sc # v1.4.3
import scrublet as scr
import sys
import bbknn
from statsmodels import robust
import sys
import matplotlib.pyplot as plt
import forceatlas2


# In[2]:


import os.path
os.chdir('/home/ngr18/hcaskin/raw_data')

adata=sc.read_h5ad('adata_scrubletremoved.h5ad')

sc.pp.filter_cells(adata, min_genes=50)
sc.pp.filter_genes(adata, min_cells=3)
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1


sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True)

adata = adata[adata.obs['n_genes'] < 6000, :]
adata = adata[adata.obs['n_genes'] > 400, :]
adata = adata[adata.obs['n_counts'] > 1000, :]
adata = adata[adata.obs['percent_mito'] < 0.2, :]

adata.raw = adata
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, adata.X, min_mean=0.0125, max_mean=5, min_disp=0.25)
sc.pl.highly_variable_genes(adata)


filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=0.0125, max_mean=5, min_disp=0.25)
adata_filtered = adata[:, filter_result.gene_subset]

sc.tl.pca(adata_filtered, n_comps=50)

import bbknn
bbknn.bbknn(adata_filtered, batch_key='donor_id', copy=False)

sc.tl.leiden(adata_filtered, resolution=1.6, key_added='leiden')
sc.tl.umap(adata_filtered)


sc.tl.rank_genes_groups(adata_filtered, 'leiden', method='wilcoxon', use_raw=False, n_genes=500)

result = adata_filtered.uns['rank_genes_groups']
groups = result['names'].dtype.names
adata_DEGs=pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(500)


initial_clustering_genes=['GATA3', 'KRTDAP', 'KRT1', 'DMKN', 'ITGA6', 'KRT14', 'KRT5', 'S100A2', 'CRYAB', 'TYRP1', 
'PTGDS', 'DCT', 'PMEL', 'FGL2', 'PMP22', 'LUM', 'APOD', 'NT5E', 'TNFAIP6', 'DCN', 
'CFD', 'RGS5'cytof_genes=['GATA3', 'KRTDAP', 'KRT1', 'DMKN', 'ITGA6', 'KRT14', 'KRT5', 'S100A2', 'CRYAB', 'TYRP1', 
'PTGDS', 'DCT', 'PMEL', 'FGL2', 'PMP22', 'LUM', 'APOD', 'NT5E', 'TNFAIP6', 'DCN', 
'CFD', 'RGS5', 'MGP', 'TFPI', 'CCL21', 'TFF3', 'FABP4', 'CD34', 'PECAM1', 'TM4SF1',
'SPARCL1', 'SERPINE1', 'KIT', 'TPSB2', 'TPSAB1', 'CTSG', 'GNLY', 'XCL2', 'XCL1', 
'KLRB1', 'IL7R', 'CXCR4', 'CD3E', 'CD8A', 'CD8B', 'CD7', 'CD3D', 'TRAC', 
'CD5', 'CD4', 'FOXP3', 'TIGIT', 'CD247', 'IGKC', 'IGLC2', 'IGHA1', 'GZMB', 
'JCHAIN', 'PLEK', 'ACOT7', 'HLA-DQB2', 'HLA-DQA1', 'CCR7', 'HLA-DRA', 'CXCL8', 'RNASE1', 
'CCL3','CD207','CD79A']

sc.tl.dendrogram(adata_filtered, groupby='leiden')

