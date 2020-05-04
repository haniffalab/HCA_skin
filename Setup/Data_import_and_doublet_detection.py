#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import scanpy as sc # v1.4.3\n",
import scrublet as scr
import sys
import bbknn
from statsmodels import robust
import sys
import matplotlib.pyplot as plt
import forceatlas2 as fa
import os.path

###Set working directory
os.chdir('/home/ngr18/hcaskin/raw_data')


###############################################################################
#%% Settings
sc.settings.set_figure_params(dpi=150)
###############################################################################
#%% Load -- csv prepared from the manifest xlsx file
manifest = pd.read_csv('sample_manifest.csv')
#manifest['sample_id']
#manifest['vdj']
###############################################################################
#%% Load
path = './'
names = manifest['sample_id'].tolist()
filenames = [os.path.join(path, n) for n in names]



#%% Concatenate and save
adatas = [sc.read_10x_mtx(filename, cache=True) for filename in filenames]

adata = adatas[0].concatenate(adatas[1:], join='inner', batch_key='sample_id', batch_categories=names, index_unique='-')

save_file = 'adata_raw.h5ad'
adata.write(save_file)



###############################################################################
#%% Metadata
for metadata in manifest.columns[1:-1]: # sample_id and vdj not needed
    print(metadata)
    adata.obs[metadata] = adata.obs['sample_id']
    replacement = manifest[metadata].tolist()
    adata.obs[metadata].replace(to_replace=names, value=replacement, inplace=True)

#adata.obs['donor_id']
#adata.obs_keys
###############################################################################


meta_10x_channels = 'sample_id'

RUNs, DSs, CELLs, THRs, MEDs, MADs, CUTs, no_thr = [], [], [], [], [], [], [], []

# Loop through channels in anndata object:
orig_stdout = sys.stdout
sys.stdout = open('scrublet_output_file_mad.txt', 'w')

for run in adata.obs[meta_10x_channels].unique():
    print(run)
    ad = adata[adata.obs[meta_10x_channels] == run, :]
    x = ad.X
    scrub = scr.Scrublet(x)
    ds, prd = scrub.scrub_doublets()
    RUNs.append(run)
    DSs.append(ds)
    CELLs.append(ad.obs_names)
    # MAD calculation of threshold:
    MED = np.median(ds)
    MAD = robust.mad(ds)
    CUT = (MED + (4 * MAD))
    MEDs.append(MED)
    MADs.append(MAD)
    CUTs.append(CUT)

    try:  # not always can calculate automatic threshold
        THRs.append(scrub.threshold_)
        print('Threshold found by scrublet')
    except:
        THRs.append(0.4)
        no_thr.append(run)
        print('No threshold found, assigning 0.4 to', run)
        scrub.call_doublets(threshold=0.4) # so that it can make the plot
    fig = scrub.plot_histogram()
    fig[0].savefig(run + '.png')

    # Alternative histogram for MAD-based cutoff
    scrub.call_doublets(threshold=CUT)
    fig = scrub.plot_histogram()
    fig[0].savefig(run + '_mad_' + '.png')
    plt.close('all')
    print()
    print()

print()
print('The following sample(s) do not have automatic threshold:')
print(no_thr)

sys.stdout.close()
sys.stdout = orig_stdout

ns = np.array(list(map(len, DSs)))

tbl = pd.DataFrame({
    'run': np.repeat(RUNs, ns),
    'ds': np.concatenate(DSs),
    'thr': np.repeat(THRs, ns),
    'mad_MED': np.repeat(MEDs, ns),
    'mad_MAD': np.repeat(MADs, ns),
    'mad_thr': np.repeat(CUTs, ns),
    }, index=np.concatenate(CELLs))

tbl['auto_prd'] = tbl['ds'] > tbl['thr']
tbl['mad_prd'] = tbl['ds'] > tbl['mad_thr']

tbl.to_csv('doublets_score_mad.csv', header=True, index=True)



adata.obs['mad_prd']=tbl['mad_prd']
adata = adata[tbl['mad_prd'] != True]


save_file = 'adata_scrubletremoved.h5ad'
adata.write(save_file)

