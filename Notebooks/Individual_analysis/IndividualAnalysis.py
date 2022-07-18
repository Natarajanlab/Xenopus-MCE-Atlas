import warnings
import os
import numpy as np
import pandas as pd
from itertools import chain
from sklearn.preprocessing import StandardScaler
from scipy.stats import gaussian_kde

import matplotlib
from matplotlib import font_manager
import matplotlib.pyplot as plt
import palantir

stages = ['st08','st10.5','st12','st13','st16','st18','st20','st22','st24','st27']

def plot_ind_cell_clusters(adatas, clusters, stages, cluster_colors):
    """Plot cell clusters on the tSNE map
    :param tsne: tSNE map
    :param clusters: Results of the determine_cell_clusters function
    """

    # Cluster colors
    n_clusters = 15

    # Set up figure
    n_cols = 10
    n_rows = 1
    fig = plt.figure(figsize=[ n_cols*3, 1*3])
    gs = plt.GridSpec(n_rows + 2, n_cols+2,
                      height_ratios=np.append([0.75, 0.75], np.repeat(1, n_rows)))

    for idx, i in enumerate(stages):
        
        layout = pd.read_csv('Ind_coord/coordinates_{}.txt'.format(i), header = None, index_col = 0)
        layout.index = adatas[idx].obs_names
        layout.columns = ['x', 'y']
        row = int(np.floor(idx / n_cols))
        ax = plt.subplot(gs[0:2, idx])
        ax.scatter(layout['x'], layout['y'], s=3,
               c=cluster_colors[clusters[layout.index]])
    
        ax.set_axis_off()
        ax.set_title(i + '\n{} Cells'.format(len(adatas[idx].obs_names)), fontsize=10)
        ax.set_rasterized(True)
    plt.show()
        
    
def plot_gene_expression(adatas, data, stages, genes, plot_scale=False,
                         n_cols=5, percentile=0, cmap=matplotlib.cm.Spectral_r):
    """ Plot gene expression on tSNE maps
    :param genes: Iterable of strings to plot on tSNE
    """
    fig = palantir.plot.FigureGrid(10, len(stages))

    

    for idx, st, ax in zip(range(0,len(stages)), stages, fig):
            #data_sub = counts.loc[counts.index.str.startswith(st), :]

        layout = pd.read_csv('Ind_coord/coordinates_{}.txt'.format(st), header = None, index_col = 0)
        layout.index = [st + '_' + x for x in list(adatas[idx].obs_names)]
        layout.columns = ['x', 'y']


        # Data
        c = data.loc[list(layout.index), genes]
        vmin = np.percentile(c[~np.isnan(c)], percentile)
        vmax = np.percentile(c[~np.isnan(c)], 100 - percentile)

        ax.scatter(layout['x'], layout['y'], s=3, color='lightgrey')
        ax.scatter(layout.loc[list(layout.index), 'x'], layout.loc[list(layout.index), 'y'], s=3,
                       c=c, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_axis_off()
        ax.set_title(genes)

        if plot_scale:
            normalize = matplotlib.colors.Normalize(
                    vmin=vmin, vmax=vmax)
            cax, _ = matplotlib.colorbar.make_axes(ax)
            matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)

        ax.set_rasterized(True)
    
    plt.show()
    
    
    
def plot_gene_expression_umap(adatas, data, stages, genes, plot_scale=False,
                         n_cols=5, percentile=0, cmap=matplotlib.cm.Spectral_r):
    """ Plot gene expression on tSNE maps
    :param genes: Iterable of strings to plot on tSNE
    """
    fig = palantir.plot.FigureGrid(10, len(stages))

    

    for idx, st, ax in zip(range(0,len(stages)), stages, fig):
            #data_sub = counts.loc[counts.index.str.startswith(st), :]

        layout = pd.read_csv('Ind_coord/umap_coord_{}.txt'.format(st), header = None, index_col = 0)
        layout.reset_index(inplace = True)
        layout.index = [st + '_' + x for x in list(adatas[idx].obs_names)]
        layout.columns = ['x', 'y']


        # Data
        c = data.loc[list(layout.index), genes]
        vmin = np.percentile(c[~np.isnan(c)], percentile)
        vmax = np.percentile(c[~np.isnan(c)], 100 - percentile)

        ax.scatter(layout['x'], layout['y'], s=3, color='lightgrey')
        ax.scatter(layout.loc[list(layout.index), 'x'], layout.loc[list(layout.index), 'y'], s=3,
                       c=c, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_axis_off()
        ax.set_title(genes)

        if plot_scale:
            normalize = matplotlib.colors.Normalize(
                    vmin=vmin, vmax=vmax)
            cax, _ = matplotlib.colorbar.make_axes(ax)
            matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)

        ax.set_rasterized(True)
    
    plt.show()
    