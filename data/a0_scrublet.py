import numpy as np
import pandas as pd
import scipy
import scipy.io
from scipy.sparse import csr_matrix, csc_matrix
import anndata
import scanpy as sc
import scrublet as scr
import sys
path = '../data/megamatrix_nodoublets/'
adata = sc.read_mtx(path + 'export_skin_all_data.txt')
adata.var_names = np.genfromtxt(path + 'export_skin_all_rownames.txt', dtype=str)
metadata = pd.read_table(path + 'export_skin_all_metadata.txt', index_col=0, header=0)
adata.obs = metadata
# Scrublet section:
RUNs, DSs, PRDs, CELLs = [], [], [], []
sys.stdout = open('scrublet_output_file.txt', 'w')

for run in adata.obs['Sanger_sample_ID'].unique():
    print(run)
    ad = adata[adata.obs['Sanger_sample_ID'] == run, :]
    x = ad.X
    scrub = scr.Scrublet(x)
    ds, prd = scrub.scrub_doublets()
    RUNs.append(run)
    DSs.append(ds)
    CELLs.append(ad.obs_names)
    PRDs.append(prd)
    print()
    print()

ns = np.array(list(map(len, DSs)))

tbl = pd.DataFrame({
    'run': np.repeat(RUNs, ns),
    'ds': np.concatenate(DSs),
    'prd': np.concatenate(PRDs)
    }, index=np.concatenate(CELLs))

tbl.to_csv('doublets_score.csv', header=True, index=True)
