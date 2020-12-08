#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os.path
import pandas as pd
from matplotlib import rcParams
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
# Add annotation that was prepared separately
adata.obs['Annotation'] = 'Not_annotated'
# healthy
replacement = pd.read_pickle('replacement_kc_healthy_Annotation_general.pkl')
adata.obs['Annotation'].loc[replacement.index] = replacement

adata.obs['Annotation'] = adata.obs['Annotation'].astype(str)
# psor lesion
replacement = pd.read_csv('data_kc/kc_psor_lesion_Annot.csv', index_col=0)
adata.obs['Annotation'].loc[replacement.index] = replacement['Annotation']

# psor nonlesion
replacement = pd.read_csv('data_kc/kc_psor_nonlesion_Annot.csv', index_col=0)
adata.obs['Annotation'].loc[replacement.index] = replacement['Annotation']

# ecz lesion
replacement = pd.read_csv('data_kc/kc_eczema_lesion_Annot.csv', index_col=0)
adata.obs['Annotation'].loc[replacement.index] = replacement['Annotation']

# ecz nonlesion
replacement = pd.read_csv('data_kc/kc_eczema_nonlesion_Annot.csv', index_col=0)
adata.obs['Annotation'].loc[replacement.index] = replacement['Annotation']

adata.obs['Annotation'].unique()

# Removed low-quality ecz cells which had no annotation:
adata = adata[~adata.obs['Annotation'].isin(['Not_annotated'])]
#adata.write('kc_all_annot.h5ad')
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
#
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
sc.pl.umap(adata, color='donor_id', save='_donor_id.png')
sc.pl.umap(adata, color='Status', save='_Status.png')
sc.pl.umap(adata, color='Annotation', save='_Annotation.png')
###############################################################################
#%% Psoriasis
adata.obs['Status_Site_Annotation'] = adata.obs['Status'].astype(str) + '_' + adata.obs['Site'].astype(str) + '_' + adata.obs['Annotation'].astype(str)
sc.pl.umap(adata, color='Status_Site_Annotation', save='_Status_Site_Annotation.png')
adata.obs['Status_Site_Annotation'].cat.categories

# Dictionary of tested: reference annotation
comparisons = {
    'Psoriasis_lesion_Pre_proliferation_KC': 'Healthy_non_lesion_Pre_proliferation_KC',
    'Psoriasis_lesion_Post_proliferation_KC': 'Healthy_non_lesion_Post_proliferation_KC',
    'Eczema_lesion_Pre_proliferation_KC': 'Healthy_non_lesion_Pre_proliferation_KC',
    'Eczema_lesion_Post_proliferation_KC': 'Healthy_non_lesion_Post_proliferation_KC',
    'Psoriasis_lesion_Pre_proliferation_KC': 'Eczema_lesion_Pre_proliferation_KC',
    'Psoriasis_lesion_Post_proliferation_KC': 'Eczema_lesion_Post_proliferation_KC',
    }

res_genes_dict = dict()
for t, r in comparisons.items():
    print(t)
    print(r)
    rank_groups = [t, r]

    sc.tl.rank_genes_groups(adata, 'Status_Site_Annotation', n_genes=100, method='wilcoxon', use_raw=False, groups=rank_groups, reference=r)

    uns_name = 'rank_genes_groups_' + t
    print(uns_name)
    adata.uns[uns_name] = adata.uns['rank_genes_groups']

    pdf_filename = '_' + t + '_vs_' + r + '.pdf'
    print(pdf_filename)
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, save=pdf_filename)

    result = adata.uns[uns_name]
    groups = result['names'].dtype.names
    res_df = pd.DataFrame( {
        group + '_' + key: result[key][group]
        for group in groups for key in ['names', 'pvals', 'pvals_adj', 'logfoldchanges'] })
    csv_filename = 'rank_genes_groups_' + t  + '_vs_' + r + '.csv'
    res_df.to_csv(csv_filename)
    print(csv_filename)

    res_genes_dict[t + '_vs_' + r] = res_df
    print()
    print()

#%% Preproliferation KC both in P and E:
P_pre = 'rank_genes_groups_Psoriasis_lesion_Pre_proliferation_KC'
P_post = 'rank_genes_groups_Psoriasis_lesion_Post_proliferation_KC'
E_pre = 'rank_genes_groups_Eczema_lesion_Pre_proliferation_KC'
E_post = 'rank_genes_groups_Eczema_lesion_Post_proliferation_KC'

res_genes_dict[P_pre].columns = ['names', 'pvals', 'pvals_adj', 'logfoldchanges']
P_pre_genes = set(res_genes_dict[P_pre].names)
res_genes_dict[E_pre].columns = ['names', 'pvals', 'pvals_adj', 'logfoldchanges']
E_pre_genes = set(res_genes_dict[E_pre].names)
len(P_pre_genes & E_pre_genes)
P_E_pre = P_pre_genes & E_pre_genes

# Psoriasis and eczema pre-proliferation DE genes vs healthy:
P_E_pre_df = res_genes_dict[P_pre][res_genes_dict[P_pre].names.isin(P_E_pre)]
P_only_pre_df = res_genes_dict[P_pre][~res_genes_dict[P_pre].names.isin(P_E_pre)]
E_only_pre_df = res_genes_dict[E_pre][~res_genes_dict[E_pre].names.isin(P_E_pre)]

res_genes_dict[P_post].columns = ['names', 'pvals', 'pvals_adj', 'logfoldchanges']
P_post_genes = set(res_genes_dict[P_post].names)

res_genes_dict[E_post].columns = ['names', 'pvals', 'pvals_adj', 'logfoldchanges']
E_post_genes = set(res_genes_dict[E_post].names)

len(P_post_genes & E_post_genes)
P_E_post = P_post_genes & E_post_genes

# Psoriasis and eczema post-proliferation DE genes vs healthy:
P_E_post_df = res_genes_dict[P_post][res_genes_dict[P_post].names.isin(P_E_post)]
P_only_post_df = res_genes_dict[P_post][~res_genes_dict[P_post].names.isin(P_E_post)]
E_only_post_df = res_genes_dict[E_post][~res_genes_dict[E_post].names.isin(P_E_post)]
###############################################################################
#%% Make new annotation columns:
adata.obs['Status_2'] = adata.obs['Status']
clusters = ['Healthy', 'Psoriasis', 'Eczema']
celltypes = ['HS', 'P', 'E']
cell_dict = dict(zip(clusters, celltypes))
adata.obs['Status_2'] = adata.obs['Status_2'].map(cell_dict)
sc.pl.umap(adata, color='Status_2', save='Status_2.png', )

adata.obs['Site_2'] = adata.obs['Site']
clusters = ['lesion', 'non_lesion']
celltypes = ['L', 'NL']
cell_dict = dict(zip(clusters, celltypes))
adata.obs['Site_2'] = adata.obs['Site_2'].map(cell_dict)
sc.pl.umap(adata, color='Site_2', save='Site_2.png', )
###############################################################################
#%% Combine annotation to get desired labels:
adata.obs['Status_Annotation_2_Site'] = (
    adata.obs['Status_2'].astype(str) +
    '_' +
    adata.obs['Annotation_2'].astype(str) +
    '_' +
    adata.obs['Site_2'].astype(str)
    )

# This section removes the 'NL' label from healthy samples:
adata.obs['annotation_multi'] = adata.obs['Status_Annotation_2_Site']
clusters = ['E_post_NL', 'E_pre_NL', 'E_pro_NL', 'E_post_inf_NL', 'E_post_L', 'E_pre_L', 'E_post_inf_L', 'E_pro_L', 'P_pre_L', 'P_post_L', 'P_post_inf_L', 'P_pro_L', 'P_pre_NL', 'P_post_NL', 'P_post_inf_NL', 'P_pro_NL', 'HS_post_NL', 'HS_post_inf_NL', 'HS_pre_NL', 'HS_pro_NL']
celltypes = ['E_post_NL', 'E_pre_NL', 'E_pro_NL', 'E_post_inf_NL', 'E_post_L', 'E_pre_L', 'E_post_inf_L', 'E_pro_L', 'P_pre_L', 'P_post_L', 'P_post_inf_L', 'P_pro_L', 'P_pre_NL', 'P_post_NL', 'P_post_inf_NL', 'P_pro_NL', 'HS_post', 'HS_post_inf', 'HS_pre', 'HS_pro']
cell_dict = dict(zip(clusters, celltypes))
adata.obs['annotation_multi'] = adata.obs['annotation_multi'].map(cell_dict)
adata.obs['annotation_multi'].unique()
sc.pl.umap(adata, color='annotation_multi', save='_annotation_multi.png', )
###############################################################################
#%% Prepare % tables
status_site_val_counts = adata.obs['Status_2_Site_2'].value_counts()  # annotation previously prepared
inf_val = adata.obs['Status_Annotation_2_Site'].value_counts()

celltype_list = []
status_site_list = []
pct_list = []
for status in adata.obs['Status_2'].unique():
    print(status)
    for site in adata.obs['Site_2'].unique():
        print(site)
        status_site = status + '_' + site
        try:
            total = status_site_val_counts[status_site]
        except:
            pass
        else:
            print('Writing %s' % status_site)
            for label in ['pre', 'pro', 'post', 'post_inf']:
                lab = status + '_' + label + '_' + site
                print(lab)

                pct = inf_val[lab] / total * 100

                celltype_list.append(label)
                status_site_list.append(status_site)
                pct_list.append(pct)

d = {'Celltype': celltype_list,
     'Status_site': status_site_list,
     'Proportion': pct_list,
     }
pct_inf_df = pd.DataFrame(d)

pct_inf_df.to_csv('pct_inf_df_proportion.csv', index=False)
# The barplot in the figure was prepared using the R/barplot.R script
###############################################################################
# Get counts for R script glm method:
# Get counts for each condition_donor_celltype:
cell_counts = adata.obs['donor_id_Annotation_2_Site_2'].value_counts()
batch = []
variable = []
value = []
timepoint = []

for donor in adata.obs['donor_id'].unique():
    print(donor)
    for site in adata.obs['Site_2'].unique():
        print(site)
        for celltype in adata.obs['Annotation_2'].unique():
            label = donor + '_' + celltype + '_' + site
            try:
                value.append(cell_counts[label])
            except:
                pass
            else:
                print('Writing %s' % label)
                batch_value = donor + '_' + site
                batch.append(batch_value)

                variable.append(celltype)

                timepoint.append(donor[0] + '_' + site)

d = {'batch': batch,
    'variable': variable,
    'value': value,
    'timepoint': timepoint,
    }
counts_df = pd.DataFrame(d)
counts_df.to_csv('cell_counts_df.csv')
###############################################################################
#%% Plot results
cat_ordered = ['HS_post_inf_NL', 'HS_post_NL', 'HS_pro_NL', 'HS_pre_NL', 'E_post_inf_NL', 'E_post_NL', 'E_pro_NL', 'E_pre_NL', 'E_post_inf_L', 'E_post_L', 'E_pro_L', 'E_pre_L',  'P_post_inf_NL', 'P_post_NL', 'P_pro_NL', 'P_pre_NL', 'P_post_inf_L', 'P_post_L', 'P_pro_L', 'P_pre_L',]
adata.obs['Status_2_Annotation_final_Site_2'].cat.reorder_categories(pd.Index(cat_ordered), inplace=True)

rcParams['figure.figsize'] = 7, 6
disease_genes = [
# Psoriasis and eczema pre-proliferation
    'S100A7',
    'S100A8',
    'S100A9',
    'TMSB10',
    'ALDOA',
    'IFITM3',
    'IFI27',
    'KRT6A',
    'KRT6B',
    'KRT6C',
    'KRT16',
    'IFITM1',
    'MIF',
# Psoriasis and eczema post-proliferation
    'SERPINB4',
    'CTSC',
    'CST3',
    'C10orf99',
    'GJB2',
# Psoriasis pre only
    'HSPA8',
    'HSPA5',
    'FABP5',
    'HSPH1',
    'GJB6',
# Psoriasis post only
    'CRABP2',
    'AKR1B10',
    'PI3',
    'SPRR2A',
    'CALML3',
# Eczema pre only
    'RPS4Y1',
    'ARPC1B',
# Eczema post only
    'SERPINB13',
    'ANGPTL4',
    'PLSCR3',
   ]

#%% Figure 3H
sc.pl.dotplot(
    adata,
    var_names=disease_genes,
    groupby='Status_2_Annotation_final_Site_2',
    use_raw=False ,
    save='_FIGURE_disease_genes.png')
disease_genes = pd.Series(disease_genes)[pd.Series(disease_genes).isin(data_genes)].to_list()
len(disease_genes)

rcParams['figure.figsize'] = 8, 6
sc.pl.dotplot(
    adata,
    var_names=list(set(disease_genes)),
    groupby='Status_2_Annotation_final_Site_2',
    use_raw=False,
    save='_FIGURE_eczema_frompaper2.png')

sc.pl.dotplot(
    adata,
    var_names=disease_genes,
    groupby='Status_2_Annotation_final_Site_2',
    use_raw=False,)
    save='_FIGURE_eczema_disease_genes.png')
