import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvdj
import random
from scipy import stats

def downsampling(collection_df, meta, category, cutoff=0):
    # collection_df: df of annotation (celltype, species) and sample (donor, region) for each meaurement (cell, specimen)
    # meta: group by this column of collection_df
    # category: count by this column of collection_df
    collection_dict = dict()
    subset_lengths = []
    for s in collection_df[meta].unique():
        subset_s = collection_df.loc[collection_df[meta] == s][category]

        s_length = len(subset_s)
        if s_length > cutoff:
            subset_lengths.append(s_length)
            collection_dict[s] = subset_s.astype(str)
        else:
            print('Sample "%s" skipped (count too low: %s)' % (s, s_length))

    shannon_index_dict = dict()
    x = min(subset_lengths)

    for s, subset_s in collection_dict.items():
        shannon_index_dict[s] = []

        l = len(subset_s)
        print()
        print('Subsampling sample %s' % s)
        print()
        for r in range(0, 100):
            print('Iteration %d ' % r, end='\r')
            x_rand_index = random.sample(range(l), x)
            subset_s_x = subset_s.iloc[x_rand_index]
            sd = dict(subset_s_x.value_counts())
            shannon_index = pyvdj.shannon(sd)
            shannon_index_dict[s].append(shannon_index)

        plt.hist(shannon_index_dict[s])
        plt.title(s)
        plt.xlabel('Shannon index')
        plt.show()
        plt.close()

    return shannon_index_dict
