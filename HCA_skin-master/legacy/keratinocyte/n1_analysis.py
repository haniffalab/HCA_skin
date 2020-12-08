#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
#%% Import
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os.path
import pandas as pd
import scanpy as sc # v1.4.3
import sys
###############################################################################
#%% Settings
sc.settings.set_figure_params(dpi=100)
###############################################################################
#%% Scanpy workflow
adata = sc.read_h5ad('KCs.h5ad')
adata.shape
adata.raw = adata
###############################################################################
#%% Normalize
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
###############################################################################
#%% Logarithmize the data
sc.pp.log1p(adata)
###############################################################################
#%% Variable genes
sc.pp.highly_variable_genes(adata, max_mean=4)
sum(adata.var['highly_variable'])
sc.pl.highly_variable_genes(adata, save='.png')
###############################################################################
#%% PCA
sc.tl.pca(adata, n_comps=50)
sc.pl.pca_variance_ratio(adata, log=True)
###############################################################################
#%% BBKNN
import bbknn # v1.3.2
bbknn.bbknn(adata, batch_key='donor_id', copy=False, n_pcs=20)
###############################################################################
#%% UMAP
sc.tl.umap(adata)
sc.pl.umap(adata, color='donor_id')
sc.pl.umap(adata, color='Status')
adata.write('kc_preprocessed.h5ad')
###############################################################################
#%% Clustering
sc.tl.leiden(adata, resolution=2, key_added='leiden_res2')
sc.pl.umap(adata, color='leiden_res2', save='_leiden_res2.png', )
# Cells were subset according to disease state and sampling site, and were annotated separately.
