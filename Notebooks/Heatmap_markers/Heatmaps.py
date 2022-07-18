import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.preprocessing import StandardScaler

def two_step_order(binar, mean_exp_mat, stages, sort_scaled=True):
    if sort_scaled:
        

        scaler = StandardScaler()
        scaled = scaler.fit_transform(mean_exp_mat.T)
        mean_exp_mat = pd.DataFrame(data = scaled.T, index = mean_exp_mat.index, columns = mean_exp_mat.columns)
    order = []
    stage_genes = {}
    b = binar.sort_values(by=stages, ascending = False)
    for st in stages:
        b[b[st] == 1]
        genes = b[b[st] == 1].index.tolist()
        order.append(mean_exp_mat.loc[genes, st].sort_values().index.tolist())
        stage_genes[st] = pd.Series(mean_exp_mat.loc[genes, st].sort_values().index.tolist())
    return [x for y in order for x in y], stage_genes

def map_to_human(gene_list):
    df = pd.read_csv('/work/sduknn/Andreas/notebooks/xenopus/Final_notebooks/XENLA_XenBase20190115_prot.txt', sep = '\t', index_col = 0)
    mapper = {x:h for x,h in zip(df.index, df['Unnamed: 5'])}
    out = []
    for i in gene_list:
        try:
            a = mapper[i]
            if 'NoName' not in a:
                out.append(a)
        except:
            pass
    print('Mapped {} % to human orthologs'.format((len(out) / len(gene_list) * 100)))
    return out

def map_gene_dict(gene_dict, mapper):
    if mapper is not None:
        for group in gene_dict.keys():
            h_genes = []
            for gene in gene_dict[group]:
                try:
                    h_genes.append(mapper[gene])
                except:
                    pass
            gene_dict[group]= [x for y in h_genes for x in y]
        return gene_dict
    else:
        print("Please supply xenopus to human dictionary for mapping")
        
def internal_order(mean_exp, gene_groups, sort_order):
    #Binarize matrix
    binar = mean_exp.apply(np.argmax, axis = 1)
    binar = pd.crosstab(binar.index, binar.values)
    binar.columns = mean_exp.columns
    #Sort matrix
    binar_sorted = binar.sort_values(by=sort_order, ascending = False)#.index.tolist()
    store = {}
    #Sort gene groups by binarized matrix
    
    for group in gene_groups.keys():
        store[group] = []
        order = []
        for st in sort_order:
            genes = binar_sorted[binar_sorted[st] == 1].index.tolist()
            g = mean_exp.loc[genes, st].sort_values().index.tolist()
            order.append(g)
            store[group].append(pd.DataFrame([g,
                                             [st]*len(g)]))
        order = [x for y in order for x in y]

        indices = {c: i for i, c in enumerate(order)}
        gene_groups[group] = sorted(gene_groups[group], key=indices.get)
    return gene_groups, store
        
    
def group_heatmap(mean_exp, gene_groups, order = None, gene_group_order = None, plot_scaled=False, human_genes=False, mapper = None, figure_scale = True, figsize = (5,10), sort_rows = True, label_size = 8, col_colors=False, minmax = None, center = 0):
    if gene_group_order is None:
        gene_group_order = list(gene_groups.keys())
        
    if order is None:
        order = mean_exp.columns
        
    if plot_scaled:
        from sklearn.preprocessing import StandardScaler

        scaler = StandardScaler()
        scaled = scaler.fit_transform(mean_exp.T)
        mean_exp = pd.DataFrame(data = scaled.T, index = mean_exp.index, columns = mean_exp.columns)
    elif plot_scaled == 'minmax':
        from sklearn.preprocessing import MinMaxScaler
        scaler = MinMaxScaler()
        scaled = scaler.fit_transform(mean_exp.T)
        mean_exp = pd.DataFrame(data = scaled.T, index = mean_exp.index, columns = mean_exp.columns)
     
    if minmax == None:
        _min, _max = np.min(mean_exp.values), np.max(mean_exp.values)
    else: 
        _min, _max = minmax
    
    if human_genes:
        gene_groups = map_gene_dict(gene_groups, mapper)
    
    for gdx, group in enumerate(gene_group_order):
        idx = np.where([ge in mean_exp.index for ge in gene_groups[group]])[0]
        gene_groups[group] = [gene_groups[group][i] for i in idx]
    
    if figure_scale:
        fig_scale = [len(gene_groups[x]) for x in gene_group_order]
        fig_scale = np.array(fig_scale) / max(fig_scale)
    else:
        fig_scale = [1] * len(gene_group_order)
        
        
    if sort_rows:
        gene_groups, store = internal_order(mean_exp, gene_groups, list(order))
        
    fig, axs = plt.subplots( len(gene_group_order), 2, 
                            figsize = figsize,
                           gridspec_kw={'height_ratios':fig_scale, 'width_ratios':[1,0.02]})
    axs = np.ravel(axs)    

    cbar = True  
    cbar_ax = axs[1]
    order_dict = {}
    for group_idx, gdx in enumerate(range(0,len(gene_group_order)*2, 2)):
        
        group = gene_group_order[group_idx]
        sns.heatmap(mean_exp.loc[gene_groups[group], order], ax = axs[gdx], 
                  cmap = 'RdBu_r',
                      cbar=cbar, 
                    cbar_ax=cbar_ax,
                   vmin = _min, vmax = _max,
                   yticklabels=True, 
                   #col_colors =col_colors,
                    center = center
                   )
        order_dict[group] = gene_groups[group]
        axs[gdx ].set_title(group)
        cbar = False
        cbar_ax = None
        axs[gdx].yaxis.set_tick_params(labelsize=label_size)
        if gdx != len(gene_group_order*2) -2:
            axs[gdx].set_ylabel('')
            axs[gdx].set_xlabel('')
            axs[gdx].set_xticks([])
            axs[gdx].set_xticklabels([])
    for i in range(3,len(gene_group_order)*2, 2):
        axs[i].set_axis_off()
    return order_dict, store