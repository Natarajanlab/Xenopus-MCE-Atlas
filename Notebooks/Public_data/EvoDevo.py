import pandas as pd
import numpy as np
import os
import scanpy as sc
import harmony

import seaborn as sns


from scipy.stats import spearmanr
import pandas as pd

def main(datasets, metas):
    celltype_db = {}
    rows = []

    
    for i, iset in enumerate(datasets):
        row = []
        row_dataset = metas[i].Dataset[0]
        for j, jset in enumerate(datasets):

            col_dataset = metas[j].Dataset[0]


            merge = pd.concat([datasets[i], datasets[j]], axis = 1)
            merge.dropna(inplace = True)

            merge_anno = pd.concat([metas[i], metas[j]], axis = 0)

            # DE
            counts1 = datasets[i].loc[merge.index]
            counts2 = datasets[j].loc[merge.index]

            print('Doing DE round {}.{}'.format(i,j))

            df1 = do_de(counts1, metas[i])
            df2 = do_de(counts2, metas[j])
            
            #Save DE
            df1['cell-type'] = metas[i]['Dataset'][0] + '_' + df1['cell-type'].astype(str)
            for c in df1['cell-type'].unique():
                celltype_db[c] = df1[df1['cell-type'] == c].gene.tolist()
            
            

            DE_genes = list(set(df1.gene).intersection(set(df2.gene)))


            merge = merge.T

            merge['Dataset'] = merge_anno['Dataset']
            merge['Cell-type'] = merge_anno['Dataset'].astype(str) + '_' + merge_anno['Cell-type'].astype(str)
            merge['Cell-type'] = merge['Cell-type']
            mean = merge.groupby(['Dataset', 'Cell-type']).mean().T

            glob_min = merge.groupby(['Dataset']).mean().T

            for c in mean.columns:
                mean[c] = mean[c].div(glob_min[c[0]], axis = 0)


            mean_DE = mean.loc[DE_genes]
            spearman, pvalues = calculate_spearman(mean_DE)

            row.append(spearman.loc[row_dataset, col_dataset].astype(float))
        rows.append(pd.concat(row, axis = 1))
            
    mat = pd.concat(rows, axis = 0)   

    return mat, celltype_db

    
    
def do_de(data, meta):
    
    adata = sc.AnnData(data.T)
    adata.obs = meta

    
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata   
    sc.tl.rank_genes_groups(adata, 'Cell-type', method='t-test')
    res = adata.uns['rank_genes_groups']

    groups = res['names'].dtype.names
    test = pd.DataFrame(
            {group + '_' + key[:1]: res[key][group]
            for group in groups for key in ['names', 'pvals', 'logfoldchanges']})

    n = test.loc[:,test.columns.str.endswith('n')]
    p = test.loc[:,test.columns.str.endswith('p')]
    l = test.loc[:,test.columns.str.endswith('l')]

    n.columns = [x[:-2] for x in n.columns]
    p.columns = [x[:-2] for x in p.columns]
    l.columns = [x[:-2] for x in l.columns]

    n = n.melt()
    n.columns = ['cell-type', 'gene']

    p = p.melt()
    p.columns = ['cell-type', 'p-value']

    l = l.melt()
    l.columns = ['cell-type', 'logfoldchanges']
    df = pd.concat([n, p['p-value'], l['logfoldchanges']], axis = 1)

    df = df[(df['p-value']<0.05) & (df['logfoldchanges']>1)]

    df = df[(~df.gene.str.startswith('LOC') ) & 
       (~df.gene.str.startswith('Xetrov')) ]
    
    return df


def calculate_spearman(df):
    df = df.dropna()._get_numeric_data()
    dfcols = df.columns
    pearson = pd.DataFrame(index = dfcols, columns = dfcols)
    pvalues = pd.DataFrame(index = dfcols, columns = dfcols)
    
    for r in dfcols:
        for c in dfcols:
            p = spearmanr(df[r], df[c])
            pearson.loc[r,c] = float(p[0])
            pvalues.loc[r,c] = float(p[1])
    return pearson, pvalues