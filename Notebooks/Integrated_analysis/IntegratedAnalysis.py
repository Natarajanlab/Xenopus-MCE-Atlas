import warnings
import os
import numpy as np
import pandas as pd
from itertools import chain
from sklearn.preprocessing import StandardScaler
from scipy.stats import gaussian_kde
import palantir
import matplotlib.pyplot as plt

import matplotlib
from matplotlib import font_manager



def plot_cell_clusters(tsne, clusters, cluster_colors, cluster_order):
    """Plot cell clusters on the tSNE map
    :param tsne: tSNE map
    :param clusters: Results of the determine_cell_clusters function
    
    #Modified palantir method!
    """

    # Cluster number of clusters
    n_clusters = len(clusters.unique())

    # Set up figure
    n_cols = 20
    n_rows = int(np.ceil(n_clusters / n_cols))
    fig = plt.figure(figsize=[2 * n_cols, 2 * (n_rows + 2)])
    gs = plt.GridSpec(n_rows + 2, n_cols,
                      height_ratios=np.append([0.75, 0.75], np.repeat(1, n_rows)))

    # Clusters
    ax = plt.subplot(gs[0:2, 2:4])
    ax.scatter(tsne['x'], tsne['y'], s=3,
               c=cluster_colors[clusters[tsne.index]])
    ax.set_axis_off()
    ax.set_rasterized(True)

    # Branch probabilities
    for i, cluster in enumerate(cluster_order):
        row = int(np.floor(i / n_cols))
        ax = plt.subplot(gs[row + 2, i % n_cols])
        ax.scatter(tsne.loc[:, 'x'], tsne.loc[:, 'y'], s=3, color='lightgrey')
        cells = clusters.index[clusters == cluster]
        ax.scatter(tsne.loc[cells, 'x'], tsne.loc[cells, 'y'],
                   s=3, color=cluster_colors[cluster])
        ax.set_axis_off()
        ax.set_title(cluster, fontsize=10)
        ax.set_rasterized(True)
        
        
def plot_comp_clusters(tsne, clusters, timepoint_colors=False, clust_set = [1,2]):
    """Plot cell clusters on the tSNE map
    :param tsne: tSNE map
    :param clusters: Results of the determine_cell_clusters function
    
    #Modified palantir method!
    """

    # Cluster colors
    n_clusters = 15
    if timepoint_colors == False:
        cluster_colors = pd.Series(sns.color_palette(
            'hls', n_clusters), index=[str(i) for i in range(0,15)])
    else:
        cluster_colors = pd.Series(sns.color_palette(
            'rainbow', 10), index=['st08','st10.5','st12','st13','st16','st18','st20','st24','st27'])
        

    # Set up figure
    n_cols = 6
    n_rows = int(np.ceil(n_clusters / n_cols))
    fig = plt.figure(figsize=[2 * n_cols, 2 * (n_rows + 2)])
    gs = plt.GridSpec(n_rows + 2, n_cols,
                      height_ratios=np.append([0.75, 0.75], np.repeat(1, n_rows)))

    # Clusters
    ax = plt.subplot(gs[0:2, 2:4])
    boo = [c in clust_set for c in clusters]
    cells = clusters.index[boo]
    ax.scatter(tsne.loc[:, 'x'], tsne.loc[:, 'y'], s=3, color='lightgrey')

    ax.scatter(tsne.loc[cells, 'x'], tsne.loc[cells, 'y'],
               s=3, c=cluster_colors[clusters[cells]])
    ax.set_axis_off()
    ax.set_title('clusters: {}'.format(clust_set), fontsize=10)
    ax.set_rasterized(True)

    

#Modified palantir methods!
def plot_palantir_results(pr_res, tsne):
    """ Plot Palantir results on tSNE
    """

    # Set up figure
    n_branches = pr_res.branch_probs.shape[1]
    n_cols = 6
    n_rows = int(np.ceil(n_branches / n_cols))
    fig = plt.figure(figsize=[2 * n_cols, 2 * (n_rows + 2)])
    gs = plt.GridSpec(n_rows + 2, n_cols,
                      height_ratios=np.append([0.75, 0.75], np.repeat(1, n_rows)))
    cmap = matplotlib.cm.plasma
    # Pseudotime
    ax = plt.subplot(gs[0:2, 1:3])
    c = pr_res.pseudotime[tsne.index]
    ax.scatter(tsne.loc[:, 'x'], tsne.loc[:, 'y'], s=3,
               cmap=matplotlib.cm.plasma, c=c)
    normalize = matplotlib.colors.Normalize(
        vmin=np.min(c), vmax=np.max(c))
    cax, _ = matplotlib.colorbar.make_axes(ax)
    cbar = matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
    ax.set_axis_off()
    ax.set_title('Pseudotime')
    ax.set_rasterized(True)

    # Entropy
    ax = plt.subplot(gs[0:2, 3:5])
    c = pr_res.entropy[tsne.index]
    ax.scatter(tsne.loc[:, 'x'], tsne.loc[:, 'y'], s=3,
               cmap=matplotlib.cm.plasma, c=c)
    normalize = matplotlib.colors.Normalize(
        vmin=np.min(c), vmax=np.max(c))
    cax, _ = matplotlib.colorbar.make_axes(ax)
    cbar = matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
    ax.set_axis_off()
    ax.set_title('Differentiation potential')
    ax.set_rasterized(True)

    # Branch probabilities
    order = [2, 3, 1, 4, 0, 5]
    row = 2
    for i, branch in enumerate(pr_res.branch_probs.columns):
        row = int(np.floor(i / n_cols))
        ax = plt.subplot(gs[row + 2, order[i]])
        c = pr_res.branch_probs.loc[tsne.index, branch]
        ax.scatter(tsne.loc[:, 'x'], tsne.loc[:, 'y'], s=3,
                   cmap=matplotlib.cm.plasma, c=c)
        normalize = matplotlib.colors.Normalize(
            vmin=np.min(c), vmax=np.max(c))
        cax, _ = matplotlib.colorbar.make_axes(ax)
        cbar = matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
        ax.set_axis_off()
        ax.set_title(branch, fontsize=10)
        ax.set_rasterized(True)
        
    

def plot_gene_expression(data, tsne, genes, plot_scale=False,
                         n_cols=5, percentile=0, cmap=matplotlib.cm.Spectral_r):
    """ Plot gene expression on tSNE maps
    :param genes: Iterable of strings to plot on tSNE
    """

    not_in_dataframe = set(genes).difference(data.columns)
    if not_in_dataframe:
        if len(not_in_dataframe) < len(genes):
            print('The following genes were either not observed in the experiment, '
                  'or the wrong gene symbol was used: {!r}'.format(not_in_dataframe))
        else:
            print('None of the listed genes were observed in the experiment, or the '
                  'wrong symbols were used.')
            return

    # remove genes missing from experiment
    genes = pd.Series(genes)[pd.Series(genes).isin(data.columns)]

    # Plot
    cells = data.index.intersection(tsne.index)
    fig = palantir.plot.FigureGrid(len(genes), n_cols)

    for g, ax in zip(genes, fig):
        # Data
        c = data.loc[cells, g]
        vmin = np.percentile(c[~np.isnan(c)], percentile)
        vmax = np.percentile(c[~np.isnan(c)], 100 - percentile)

        ax.scatter(tsne['x'], tsne['y'], s=3, color='lightgrey')
        ax.scatter(tsne.loc[cells, 'x'], tsne.loc[cells, 'y'], s=3,
                   c=c, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_axis_off()
        ax.set_title(g)

        if plot_scale:
            normalize = matplotlib.colors.Normalize(
                vmin=vmin, vmax=vmax)
            cax, _ = matplotlib.colorbar.make_axes(ax)
            matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
            
        ax.set_rasterized(True)