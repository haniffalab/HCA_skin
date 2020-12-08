annot = 'Annotation'  # obs column name
clusters = ['c1', 'c2', 'c3', 'c4']  # clusters or cell types along the path
###############################################################################
import scanpy as sc # v1.4.3
from scipy import sparse, io
# Subset data to include only cells in your path
adata = sc.read_h5ad('kc_healthy_preprocessed_clustered.h5ad')
adata = adata[adata.obs[annot].isin(clusters)]
adata.raw
###############################################################################
# Export genes
genes = adata.var.copy()
genes['gene_short_name'] = genes.index  # otherwise monocle throws: "Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions"
genes[['gene_short_name', 'highly_variable']].to_csv('pdt_genes.csv', header=True)
# Export annotation
export_columns = ['dpt_pseudotime', annot, ]  # may specify additional columns
metadata = adata.obs[export_columns]
metadata.to_csv('pdt_metadata.csv')
# Export counts as sparse matrix from raw (adata.raw.X must contain raw count values)
c_sparse = adata.raw.X.transpose()
c_sparse = c_sparse.astype(int)

io.mmwrite('pdt_counts.mtx', c_sparse)

# This produces 3 files that will be used by R/run_monocle.R
