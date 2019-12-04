#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import scanpy as sc # v1.4.3
from matplotlib import rcParams
adata.obs['no_annot'] = 'none'  # so that plot orders by pdt
adata.obs['no_annot'] = adata.obs['no_annot'].astype('category')
adata.obs['no_annot']
nodes = ['none']
rcParams['figure.figsize'] = 4, 6
################################################################
with open('genes.pdt.std.modified.txt', 'r') as f:  # one gene per line, compiled from the results of Monocle
    genes = f.readlines()
f.closed
genes = [x.strip() for x in genes]

import matplotlib as mpl
colours = mpl.colors.LinearSegmentedColormap.from_list('rg', ['#3aa4c7', '#FF3030'])

plot = sc.pl.paga_path(
    adata,
    nodes,
    groups_key='no_annot',
    keys=genes,
    use_raw=True,
    normalize_to_zero_one=True,
    show_node_names=True,
    n_avg=50,
    show_yticks=True,
    show_colorbar=False,
    color_maps_annotations={'distance': 'plasma'},
    color_map=colours,
    return_data=False,
    show=False,
    save='pdt.std.heatmap.png')
