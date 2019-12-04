#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
#%% Import
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import os.path
import pandas as pd
import scanpy as sc # v1.4.3
import sys
###############################################################################
#%% Settings
sc.settings.set_figure_params(dpi=100)
###############################################################################
#%% Subsetting
adata = sc.read_h5ad('KCs.h5ad')
adata.raw = adata
adata = adata[adata.obs['Status'].isin(['Healthy'])]
adata.raw
save_file = 'kc_healthy.h5ad'
adata.write(save_file)

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
###############################################################################
#%% Variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.1, max_mean=4)
sc.pl.highly_variable_genes(adata)
###############################################################################
#%% PCA
sc.tl.pca(adata, n_comps=50)
sc.pl.pca_variance_ratio(adata, log=True)
###############################################################################
#%% BBKNN
import bbknn # v1.3.2
bbknn.bbknn(adata, batch_key='donor_id', copy=False, n_pcs=50)
###############################################################################
#%% UMAP
sc.tl.umap(adata)
sc.pl.umap(adata, color='donor_id')
sc.pl.umap(adata, color='Status')
###############################################################################
#%% Clustering
res = 2 # chosen resolution
leiden_res = 'leiden_res' + str(res)
sc.tl.leiden(adata, resolution=res, key_added=leiden_res)
sc.pl.umap(adata, color=leiden_res, save='_' + leiden_res + '.png', )
###############################################################################
#%% FDG
sc.tl.draw_graph(adata, layout='fa')
sc.pl.draw_graph(adata, layout='fa')
save_file = 'kc_healthy_preprocessed_clustered.h5ad'
adata.write(save_file)

###############################################################################
#%% Choose clustering:
sc.tl.dendrogram(adata, groupby='leiden_res2', n_pcs=None, use_rep=None, var_names=None, use_raw=None, cor_method='pearson', linkage_method='complete', key_added=None)
###############################################################################
#%% Annotation
adata.obs['Annotation_general'] = adata.obs['leiden_res2']
clusters = list(map(str, range(0, 18+1)))
celltypes = [
    'Pre_proliferation_KC',
    'Pre_proliferation_KC',
    'Pre_proliferation_KC',
    'Pre_proliferation_KC',
    'Post_proliferation_KC',
    'Post_proliferation_KC',
    'Pre_proliferation_KC',
    'Post_proliferation_KC',
    'Post_proliferation_KC',
    'Post_proliferation_KC',
    'Post_proliferation_KC',
    'Post_proliferation_KC',
    'CD83_KC',
    'Proliferating_KC',
    'Post_proliferation_KC',
    'Pre_proliferation_KC',
    'CD83_KC',
    'CD83_KC',
    'Post_proliferation_KC',
    ]
colours_general = [
    '#E87D72',
    '#0e6c8b',
    '#b8bae5',
    '#87AC34',
    ]
cell_dict = dict(zip(clusters, celltypes))
adata.obs['Annotation_general'] = adata.obs['Annotation_general'].map(cell_dict)

sc.settings.set_figure_params(dpi=150)
sc.pl.draw_graph(adata, color='Annotation_general', save='_Annotation_general.png', palette=colours_general)
sc.pl.umap(adata, color='Annotation_general', save='_Annotation_general.png', palette=colours_general)
###############################################################################
#%% Feature plots for figure B
import matplotlib.pyplot as plt
gr_col = mpl.colors.LinearSegmentedColormap.from_list('rg', ['silver', 'red'])

sc.pl.draw_graph(adata, color='CNFN', color_map=gr_col, use_raw=False)
sc.pl.draw_graph(adata, color='CNFN', save='_CNFN_red.png', color_map=gr_col, use_raw=False)
sc.pl.draw_graph(adata, color='KRT10', save='_KRT10_red.png', color_map=gr_col, use_raw=False)
sc.pl.draw_graph(adata, color='CDK1', save='_CDK1_red.png', color_map=gr_col, use_raw=False)
sc.pl.draw_graph(adata, color='KRT15', save='_KRT15_red.png', color_map=gr_col, use_raw=False)
sc.pl.draw_graph(adata, color='HMMR', save='_HMMR_red.png', color_map=gr_col, use_raw=False)
sc.pl.draw_graph(adata, color='HILPDA', save='_HILPDA_red.png', color_map=gr_col, use_raw=False)
###############################################################################
#%% Combine clusters:
adata.obs['Annotation'] = adata.obs['leiden_res2']
clusters = list(map(str, range(0, 18+1)))
celltypes = [
    'c1',
    'c1',
    'c1',
    'c5',
    'c7',
    'c7',
    'c5',
    'c3',
    'c6',
    'c3',
    'c4',
    'c3',
    'c8',
    'c2',
    'c3',
    'c1',
    'c8',
    'c8',
    'c3',
    ]
cell_dict = dict(zip(clusters, celltypes))
adata.obs['Annotation'] = adata.obs['Annotation'].map(cell_dict)

kc_colours2 = [
    '#00A9FF',
    '#7CAE00',
    '#00BFC4',
    '#218b5e',
    '#ffac7f',
    '#CD9600',
    'grey',
    '#F8766D']

sc.pl.draw_graph(adata, color='Annotation', palette=kc_colours2)
###############################################################################
#%% Figure D - FDG
rcParams['figure.figsize'] = 5, 5
sc.pl.draw_graph(adata, color='Annotation', save='_Annotation.png', palette=kc_colours2)
sc.pl.draw_graph(adata, color='Annotation', save='_Annotation_labels.png', palette=kc_colours2, legend_loc='on data',)
#%% Figure D - PAGA
sc.tl.paga(adata, groups='Annotation')
sc.pl.paga(adata, layout='fa', threshold=0.1, node_size_scale=3, edge_width_scale=1, frameon=False, node_size_power=0, random_state=8)
sc.pl.paga(adata, layout='fa', threshold=0.1, node_size_scale=3, edge_width_scale=1, frameon=False, node_size_power=0)
#%% Figure F - heatmap
dotplot2_genes = list(reversed(['CST6', 'KRT2', 'SDR9C7', 'CYP4F22', 'TGM1', 'PRSS8', 'SULT2B1', 'CERS3', 'ABCA12', 'KRTDAP', 'CKAP4', 'SPINK5', 'KLK7', 'CLIP1', 'TJP1', 'ELOVL4', 'OCLN', 'PRDM1', 'RELA', 'SLC27A4', 'NIPAL4', 'SPTSSB']))
sc.pl.dotplot(adata, dotplot2_genes, groupby='Annotation', use_raw=False, save='_FIGURE_F.png', color_map=dot_col)
###############################################################################
#%% Markers: Annotation_general
sc.tl.rank_genes_groups(adata, 'Annotation_general', method='wilcoxon', use_raw=False)
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='.pdf')
adata.uns['rank_genes_groups_Annotation_general'] = adata.uns['rank_genes_groups']
#%% Markers: Annotation
sc.tl.rank_genes_groups(adata, 'Annotation', n_genes=500, method='wilcoxon', use_raw=False)
adata.uns['rank_genes_groups_Annotation'] = adata.uns['rank_genes_groups']
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='.pdf')
result = adata.uns['rank_genes_groups_Annotation']
res_genes = pd.DataFrame(result['names'])
groups = result['names'].dtype.names
res_df = pd.DataFrame( {
    group + '_' + key: result[key][group]
    for group in groups for key in ['names', 'pvals', 'pvals_adj', 'logfoldchanges'] })
res_df.to_csv('rank_genes_groups_Annotation.csv')
###############################################################################
#%% DPT
sc.tl.diffmap(adata, n_comps=15)
sc.pl.scatter(adata, basis='diffmap', color='Annotation', save='.png')
adata.uns['iroot'] = np.flatnonzero(adata.obs['Annotation'] == 'c1')[0]
sc.tl.dpt(adata, n_dcs=10, n_branchings=0, min_group_size=0.01, allow_kendall_tau_shift=True)
sc.pl.scatter(adata, basis='diffmap', color='dpt_pseudotime', save='_pdt.png')
sc.pl.draw_graph(adata, color='dpt_pseudotime', save=('_dpt_pseudotime.png'), use_raw=False)
sc.pl.paga(adata, layout='fa', threshold=0.1, node_size_scale=3, edge_width_scale=1, frameon=False, node_size_power=0, color='dpt_pseudotime', save='_kc_paga_pdt.pdf')
#save_file = 'kc_healthy_preprocessed_clustered.h5ad'
#adata.write(save_file) # saved
#adata = sc.read_h5ad('kc_healthy_preprocessed_clustered.h5ad')
###############################################################################
#%% Save annotation
replacement = adata.obs['Annotation']  # cells of the subset
replacement.to_pickle('replacement_kc_healthy_Annotation.pkl')
replacement = adata.obs['Annotation_general']  # cells of the subset
replacement.to_pickle('replacement_kc_healthy_Annotation_general.pkl')
