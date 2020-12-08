#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified on 031019 

@author: Issac
@mod: Issac
"""

# gene expression should be at data.X
# DRs should be at data.obsm
# categoeis should be columns at data.obs

# gene expression should be at data.X
# DRs should be at data.obsm
# categoeis should be columns at data.obs

import sys
file_name     = sys.argv[1]
category      = sys.argv[2]
output_folder = sys.argv[3]

import matplotlib; matplotlib.use('Agg');
import scanpy.api as sc
from os.path import join
import os.path
from os import path
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
import pandas as pd

data = sc.read(file_name)

#Caution, un-commenting below three lines can cause errors, do so with caution
#data.X = data.raw.X
#sc.pp.normalize_per_cell(data, counts_per_cell_after=1e4)
#sc.pp.log1p(data)

#replace col names
data.obs.columns = data.obs.columns.str.replace(" ", "")

expression_path = join(output_folder, 'expression.mtx')

if path.exists(expression_path):
    print("expression.mtx is present in target folder, to save on overhead, assuming save_gene_names is present too")
    
else:
    print("Prior exported expression.mtx not found procedding write new file")
    # save expression data
    expression_data_filename = join(output_folder, 'expression')
    gene_expression = csr_matrix(data.X)
    mmwrite(expression_data_filename, gene_expression)

# save gene names
    gene_names = "\n".join(data.var_names.tolist())
    with open(join(output_folder, 'gene_names.txt'), "w") as gene_names_fobj:
        gene_names_fobj.write(gene_names)

#for aggregate matrix.py
variables = pd.read_csv(output_folder + "/variables.csv")
variable_list = list(variables['x'])
for i in variable_list:
    print ("categorical variable exported: " + i)
    # get categories
    cat_concat = ':' + data.obs[i].astype(str)
    categories = "\n".join(cat_concat.tolist())
    with open(join(output_folder, str(i)+".txt"), "w") as categories_file:
        categories_file.write(categories)





