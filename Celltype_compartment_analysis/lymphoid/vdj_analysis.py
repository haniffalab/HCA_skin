#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc  # v1.4.3
import pyvdj
import pickle
import upsetplot
###############################################################################
#%% Data loaded
sc.settings.set_figure_params(dpi=150)
adata = sc.read_h5ad('lymphoid_preprocessed.h5ad')
adata.shape
sc.pl.umap(adata)
###############################################################################
#%% Add final annotation prepared separately
sc.pl.umap(adata, color='Annotation')

replacement = pd.read_csv('annotation_cellnames_corrected/Skin_lymphoid_EZ_all_vdj_annot.csv', index_col=0)
adata.obs['Annotation'] = adata.obs['Annotation'].astype(str)
adata.obs['Annotation'].loc[replacement.index] = replacement['Annotation']

replacement = pd.read_csv('annotation_cellnames_corrected/Skin_lymphoid_healthy_all_vdj_annot.csv', index_col=0)
adata.obs['Annotation'].loc[replacement.index] = replacement['Annotation']

replacement = pd.read_csv('annotation_cellnames_corrected/Skin_lymphoid_psor_all_vdj_annot.csv', index_col=0)
adata.obs['Annotation'].loc[replacement.index] = replacement['Annotation']

colours_general = [
    '#d0e1f2',
    '#4ba2cf',
    '#329a22',
    '#00407c',
    '#ec139a',
    '#a3dc82',
    'black',
    'red',
    '#a02b9d',
    '#B482D6',
    ]
sc.pl.umap(adata, color='Annotation', palette=colours_general)
sc.pl.umap(adata, color='Annotation', legend_loc='on data', )

sc.pl.umap(adata, color='Annotation', save='_Annotation.png', )
sc.pl.umap(adata, color='Annotation', save='_Annotation_label.png', legend_loc='on data', )

colours_no_nk = [
    '#c5c5c5',
    '#c5c5c5',
    '#c5c5c5',
    '#c5c5c5',
    '#c5c5c5',
    '#a3dc82',
    'red',
    'black',
    '#a02b9d',
    '#B482D6',
    ]
sc.pl.umap(adata, color='Annotation', palette=colours_no_nk, save='_Annotation_colours_no_nk.png', )
sc.pl.umap(adata, color='Annotation', palette=colours_no_nk, legend_loc='on data', save='_Annotation_colours_no_nk_label.png',)
###############################################################################
# VDJ data
manifest = pd.read_csv('../../sample_manifest.csv')

manifest = manifest[manifest['sample_id'].isin(adata.obs['sample_id'].unique())]  # keep the ones in adata
manifest = manifest[manifest['vdj'].notnull()]

paths = '../../data/mtx/' + manifest['vdj'] + '/filtered_contig_annotations.csv'
samples = dict(zip(paths, manifest['sample_id']))
samples
cellnames = adata.obs_names
cellbarcode = cellnames.str.split("-").str[:2].str.join("-") # cell barcode part + '-1'
adata.obs['vdj_obs'] = cellbarcode.astype(str) + "_" + adata.obs['sample_id'].astype(str)
#%% Load data
adata = pyvdj.load_vdj(samples, adata, obs_col='vdj_obs', cellranger=3)
# VDJ annotation
adata = pyvdj.add_obs(adata, obs=['clonotype', 'is_clone', 'any_productive', 'all_productive'])
###############################################################################
#%% Figure 4F: we prepare an annotation column for cells part of an expanded cells
# and use this to make a df of proportion of expanded within each celltype / status-site
n = 2 # clonotypes with at least 2 clones Are considered expanded
adata.obs['vdj_expanded_n_2'] = adata.obs['vdj_clone_count']
adata.obs['vdj_expanded_n_2'][adata.obs['vdj_expanded_n_2'] < n] = 0
adata.obs['vdj_expanded_n_2'][adata.obs['vdj_expanded_n_2'] >= n] = 1

obs_vdj_df = adata.obs[adata.obs['vdj_has_vdjdata'] == True]
obs_vdj_df.keys().to_list()

t_celltypes = ['Tc', 'Treg', 'Th', 'Tc IL13 IL22', 'Tc17_Th17']

obs_vdj_df = obs_vdj_df[obs_vdj_df['Annotation'].isin(t_celltypes)]obs_vdj_df['Status_Site'] = obs_vdj_df['Status'].astype(str) + '_' + obs_vdj_df['Site'].astype(str)

t_exp_counts = obs_vdj_df.groupby(['Annotation', 'Status_Site']).vdj_expanded_n_2.value_counts()
t_exp_counts.to_csv('t_cell_expanded_proportions.csv', header=True)

t_exp_counts_df = pd.read_csv('t_cell_expanded_proportions.csv', header=0)

celltypes = ['Tc', 'Tc IL13 IL22', 'Tc17_Th17', 'Th', 'Treg',]
# Save proportions as percent of total:
Proportion = []
Status_site = []
Celltype = []

for celltype in celltypes:
    c_df = t_exp_counts_df[t_exp_counts_df['Annotation'] == celltype]
    expanded_raw = c_df[c_df['vdj_expanded_n_2'] == 1]['vdj_expanded_n_2.1'].to_list()
    notexpan_raw = c_df[c_df['vdj_expanded_n_2'] == 0]['vdj_expanded_n_2.1'].to_list()

    c_df.Status_Site
    total_raw = [x+y for x, y in zip(expanded_raw, notexpan_raw)]
    expanded = [x/y*100 for x, y in zip(expanded_raw, total_raw)]

    Proportion = Proportion + expanded
    Status_site = Status_site + c_df.Status_Site.unique().tolist()
    Celltype = Celltype + [celltype] * len(c_df.Status_Site.unique())


pct_exp_only_df = pd.DataFrame({
    'Celltype': pd.Series(Celltype),
    'Status_site': pd.Series(Status_site),
    'Proportion': pd.Series(Proportion),
})

pct_exp_only_df
pct_exp_only_df.to_csv('pct_exp_only_df.csv', index=False)
# This file is used in R/A1_barplot_percent_expanded_T_cell_types_side.R to make the figure
###############################################################################
#%% Get counts for testing difference in cell proportions, using R/A2_cell_counts.R script:
# New annotation columns:
obs_vdj_df['Site_2'] = obs_vdj_df['Site']
clusters = ['lesion', 'non_lesion']
celltypes = ['L', 'NL']
cell_dict = dict(zip(clusters, celltypes))
obs_vdj_df['Site_2'] = obs_vdj_df['Site_2'].map(cell_dict)

obs_vdj_df['donor_id_Site_2_Annotation_vdj_expanded_n_2'] = (
    obs_vdj_df['donor_id'].astype(str) +
    '_' +
    obs_vdj_df['Site_2'].astype(str) +
    '_' +
    obs_vdj_df['Annotation'].astype(str) +
    '_' +
    obs_vdj_df['vdj_expanded_n_2'].astype(str)
    )

cell_counts = obs_vdj_df['donor_id_Site_2_Annotation_vdj_expanded_n_2'].value_counts()
batch = []  # donor + site_annotation
variable = []  # expanded / not expanded (i.e. vdj_expanded_n_2)
value = []  # cell count
timepoint = []  # status (i.e. donor[0]) + site_annotation

for donor in obs_vdj_df['donor_id'].unique():
    print(donor)
    for site in obs_vdj_df['Site_2'].unique():
        print(site)
        for annotation in obs_vdj_df['Annotation'].unique():
            for celltype in obs_vdj_df['vdj_expanded_n_2'].unique():
                label = donor + '_' + site + '_' + annotation + '_' + str(celltype)
                try:
                    value.append(cell_counts[label])
                except:
                    pass
                else:
                    print('Writing %s' % label)
                    batch_value = donor + '_' + site + '_' + annotation
                    batch.append(batch_value)

                    variable.append(str(celltype))

                    timepoint.append(donor[0] + '_' + site + '_' + annotation)

d = {'batch': batch,
    'variable': variable,
    'value': value,
    'timepoint': timepoint,
    }
counts_df = pd.DataFrame(d)
counts_df.to_csv('t_cell_counts_df_site.csv')
# The csv file is used in R/A2_cell_counts.R
###############################################################################
#%% Figure 4G: heatmap of shared clonotypes between cell types
t_celltypes = ['Tc', 'Treg', 'Th', 'Tc IL13 IL22', 'Tc17_Th17']

cld_set = obs[obs['vdj_has_vdjdata']]
cld_set = cld_set[cld_set.Annotation.isin(t_celltypes)]

cld = dict()
for s_s in cld_set['Annotation'].unique():
    print(s_s)
    cld_set_s = cld_set[cld_set.Annotation == s_s]
    cld_set_s = set(cld_set_s.vdj_donor_clonotype)
    cld[s_s] = cld_set_s
    print(len(cld_set_s))
######################################
# Populate a matrix by comparing pairs
groups = list(cld.keys())
dat_list = []
for k, v in cld.items():
    print(k + '  vs')

    values = []
    for l, w in cld.items():
        if l in groups:
            # Change this section if you want to exclude/include self-comparison
            if l == k:
                values.append(np.nan)
            else:
                print(l)
                print(len(cld[k] & cld[l]))
                values.append( len(cld[k] & cld[l]) )
        else:
            print("-")
            values.append(np.nan)

    print()
    dat_list.append(values)
    groups.pop(0)

dat_list

# Plot results
##############
dat = np.array(dat_list)
dat
x_labels = list(cld.keys())
y_labels = list(cld.keys())

sc.settings.set_figure_params(dpi=100)

dat.T  # is used so that labels are next to squares

# The below functions are in heatmap_functions.py
# Adapted from https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
fig, ax = plt.subplots()
im, cbar = heatmap(dat.T, x_labels, y_labels, ax=ax,
                   cmap='magma_r', cbarlabel='Shared clonotypes')
texts = annotate_heatmap(im, data=dat.T, valfmt='{x:.0f}', textcolors=['black', 'white'], threshold=7000)  # threshold is for switching colours
ax.grid(False)
figure_shared = plt.gcf()
figure_shared.savefig('Figure_4G_shared_clonotypes_celltype.pdf')
figure_shared.savefig('Figure_4G_shared_clonotypes_celltype.ps')

###############################################################################
#%% Supplementary for Figure 4: we export the Shannon indices, then plot in R
# The 'downsampling()' function used below is defined in downsampling.py
meta = 'donor_id_Status_Site'
category = 'vdj_clonotype'

filtering_cells = adata.obs
filtering_cells['donor_id_Status_Site'] = filtering_cells['donor_id'].astype(str) + '_' + filtering_cells['Status'].astype(str) + '_' + filtering_cells['Site'].astype(str)
filtering_cells = filtering_cells[filtering_cells['vdj_has_vdjdata'] == True]

t_celltypes = ['Tc', 'Tc17_Th17', 'Tc IL13 IL22', 'Th', 'Treg',]
cutoff_values = [240, 20, 85, 494, 196]  # replicates with low cell count were excluded using these cutoffs
filter_cell_dict = dict(zip(t_celltypes, cutoff_values))

shi_df_all = pd.DataFrame()
for t, c in filter_cell_dict.items():
    collection_df = filtering_cells
    collection_df = collection_df[collection_df.Annotation == t]
    collection_df.shape
    shannon_index_dict = downsampling(collection_df, meta, category, cutoff=c)

    for s, i in shannon_index_dict.items():
        print(s)
        print(np.median(i))

    shi_dict = {}
    for s, i in shannon_index_dict.items():
        shi_dict[s] = np.median(i)

    shi_df = pd.DataFrame({
        'name': list(shi_dict),
        'shi': list(shi_dict.values()),
        'celltype': t
        })
    shi_df_all = pd.concat([shi_df_all, shi_df], ignore_index=True)

shi_df_all

# Additional columns to help plotting:
repl_dict = dict(zip(filtering_cells.donor_id_Status_Site, filtering_cells.Site))
shi_df_all['Site'] = shi_df_all.name
shi_df_all['Site'].replace(to_replace=repl_dict, inplace=True)

repl_dict = dict(zip(filtering_cells.donor_id_Status_Site, filtering_cells.donor_id))
shi_df_all['donor_id'] = shi_df_all.name
shi_df_all['donor_id'].replace(to_replace=repl_dict, inplace=True)

repl_dict = dict(zip(filtering_cells.donor_id_Status_Site, filtering_cells.Status))
shi_df_all['Status'] = shi_df_all.name
shi_df_all['Status'].replace(to_replace=repl_dict, inplace=True)

shi_df_all['Status_Site'] = shi_df_all.Status.astype(str) + "_" + shi_df_all['Site'].astype(str)

colours_t = [
    '#a3dc82',  # Tc
    'red',
    'black',  # Tc17_Th17
    '#a02b9d',
    '#B482D6',
    ]
celltypes = ['Tc', 'Tc IL13 IL22', 'Tc17_Th17', 'Th', 'Treg',]
celltype_colour_dict = dict(zip(celltypes, colours_t))

shi_df_all['colour'] = shi_df_all['celltype']
shi_df_all['colour'].replace(to_replace=celltype_colour_dict, inplace=True)

Status_Site_list = list(shi_df_all['order'].unique())
orders = [2, 1, 4, 3, 0]
order_dict = dict(zip(Status_Site_list, orders))
shi_df_all['order'] = shi_df_all['Status_Site']
shi_df_all['order'].replace(to_replace=order_dict, inplace=True)

shi_df_all.to_csv('shi_df_all_t_cells.csv')
# The csv file is used in R/A3_shi_dotplot.R
###############################################################################
