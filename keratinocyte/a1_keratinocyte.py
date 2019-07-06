import numpy as np
import pandas as pd
import scipy
from scipy.sparse import csr_matrix, csc_matrix
import anndata
import scanpy
import scanpy.api as sc
import matplotlib.pyplot as plt
#%%
sc.settings.verbosity = 4
sc.logging.print_versions()
#%%
# Section 1: prepare AnnData
adata = scanpy.api.read_h5ad('kc.h5ad')
# Calculate FA (force atlas) force-directed graph (FDG):
sc.pp.neighbors(adata, n_pcs=20, use_rep='X_harmony')
adata.uns['neighbors']
sc.settings.set_figure_params(dpi=150)
scanpy.api.tl.paga(adata, groups='Clustering_res_0_35_harmony_merged')
#%%
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color='Clustering_res_0_35_harmony_merged', layout='fa')
adata.write("./kc_scanpy.h5ad")
np.savetxt('X_FA_FDG.txt', adata.obsm['X_draw_graph_fa'])
adata.obs.index.to_series().to_csv('cell_names.txt', index=False)
#%% Plot with publication colours:
colourpalette = ["#00BFC4", "#00A9FF", "#A9A9A9", "#ffac7f", "#CD9600", "#F8766D", "#7CAE00", "#218b5e"]
sc.pl.draw_graph(adata, color='Clustering_res_0_35_harmony_merged', layout='fa', palette=colourpalette)
sc.pl.paga(adata, layout='fr', threshold=0.04, node_size_scale=3, edge_width_scale=1.5, labels=['', '', '', '', '', '', '', ''], frameon=False, node_size_power=0, save='_kc_agathr0.04.png')
sc.pl.paga(adata, layout='fr', threshold=0.04, node_size_scale=3, edge_width_scale=1.5, frameon=False, node_size_power=0, save='_kc_agathr0.04_labelled.png')
