from collections import OrderedDict
from itertools import chain
import harmony
import pandas as pd
import numpy as np


def concat_for_harmony(adatas, sample_names=None, min_cell_count=1):
    counts_dict = OrderedDict()
    for adata, sample in zip(adatas, sample_names):
        print(sample)
        counts = adata_to_df(adata)
        # Update sample names
        counts.index = sample + '_' + counts.index.astype(str)

        # Remove zero count genes
        counts = counts.loc[:, counts.sum() > 0]
        counts = counts.astype(np.int16)

        # Update dictionary
        counts_dict[sample] = counts

    # Concatenate cells and genes
    print('Concatenating data..')
    all_cells = list(chain(*[list(counts_dict[sample].index)
                             for sample in sample_names]))
    all_genes = list(
        chain(*[list(counts_dict[sample].columns) for sample in sample_names]))
    all_genes = list(set(all_genes))

    # Counts matrix
    counts = pd.DataFrame(0, index=all_cells,
                          columns=all_genes, dtype=np.int16)
    for sample in sample_names:
        sample_counts = counts_dict[sample]
        counts.loc[sample_counts.index, sample_counts.columns] = sample_counts

    # Filter out low detection genes
    gs = counts.sum()
    counts = counts.loc[:, counts.columns[gs > min_cell_count]]

    return counts

def adata_to_df(adata):
    return pd.DataFrame(data = adata.X.todense(), 
                 index = adata.obs_names,
                 columns = adata.var_names)


    