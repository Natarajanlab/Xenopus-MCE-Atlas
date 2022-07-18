import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt

def global_enrich(Hvgs, allGenes, go_db):
    
    x = []
    anno = []
    grouping = []
    fisher = []
    for go in go_db.keys():
        
        go_subset = [x for x in go_db[go] if x in allGenes]
        
        #Quantify
        N_background = sum([x in go_subset for x in allGenes])
        total_background =  len(allGenes)
    
        n_genes = sum([x in go_subset for x in Hvgs])
        total_genes = len(Hvgs)

        _, p = stats.fisher_exact([[n_genes, N_background], 
                                   [total_genes - n_genes ,  total_background - N_background]])

        #Adding it all
        x.append(n_genes / total_genes)
        anno.append(go)
        grouping.append('Global')
        fisher.append(p)
    df = pd.DataFrame({'GO':anno, 
                 'Group':grouping, 
                 'Enrichment':x,
                'p':fisher})
    return df   


def go_enrich(go_db, gene_sets, background, allgenes):

    x = []
    anno = []
    stage = []
    fisher = []
    for go in go_db.keys():
        
        go_subset = [x for x in go_db[go] if x in allgenes]
        N_background = sum([x in go_subset for x in background])
        total_background =  len(background)

        for st in gene_sets.keys():

            #Quantify
            n_genes = sum([x in go_subset for x in gene_sets[st]])
            total_genes = len(gene_sets[st])

            #Fishers
            _, p = stats.fisher_exact([[n_genes, N_background], [total_genes - n_genes ,  total_background - N_background]])

            #Adding it all
            x.append(n_genes / total_genes)
            anno.append(go)
            stage.append(st)
            fisher.append(p)
    df = pd.DataFrame({'GO':anno, 
                 'Group':stage, 
                 'Enrichment':x,
                'p':fisher
                      }
                     )
    return df

def plot_dots(enrichDF, alpha = 0.05, trim_labels = False, x_order = None):
    if trim_labels:
        enrichDF.GO = [' '.join(x[1:]).title() for x in enrichDF.GO.str.split('_')]
    

    # Mapping from column names to integer coordinates
    if x_order is not None:
        x_labels = x_order
    else:
        x_labels = [v for v in sorted(enrichDF.Group.unique())]
    y_labels = [v for v in sorted(enrichDF.GO.unique())]
    x_to_num = {p[1]:p[0] for p in enumerate(x_labels)} 
    y_to_num = {p[1]:p[0] for p in enumerate(y_labels)} 
    
    print(len(x_labels), len(y_labels))
    fig, ax = plt.subplots(figsize = (len(x_labels)*0.50 , len(y_labels)*0.40))
    size_scale = 1000
    #Plot background
    ax.scatter(
        x=enrichDF.Group.map(x_to_num), # Use mapping for x
        y=enrichDF.GO.map(y_to_num), # Use mapping for y
        s=enrichDF.Enrichment * size_scale, # Vector of sizes, proportional to size parameter
        marker='o',
        c = 'lightgrey' )    
    
    #Plot significant
    sub = enrichDF[enrichDF.p< alpha]
    scatter = ax.scatter(
        x=sub.Group.map(x_to_num), # Use mapping for x
        y=sub.GO.map(y_to_num), # Use mapping for y
        s=sub.Enrichment * size_scale, # Vector of sizes, proportional to size parameter
        marker='o',
        c = sub.p,
        cmap = 'Reds_r', 
        vmin=0, vmax = alpha, 
        )    
    
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    
    
    for s in [0.05, 0.1, 0.2, 0.3]:
        ax.scatter([], [], c='lightgray', alpha=0.3, s=s * size_scale,
                label=str(s*100) + ' %')
        plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Set Enrichment', bbox_to_anchor=(1, 0.7))
    
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    cmap = plt.get_cmap('Reds_r')
    cNorm  = colors.Normalize(vmin=0, vmax=alpha)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    scalarMap.set_array([])
    
    cax = fig.add_axes([0.8, 0.75, 0.05, 0.1])
    plt.colorbar(scalarMap, cax = cax, ticks=np.linspace(0,alpha, 3),
                                            orientation='vertical',label = 'Fishers p')
    
    ax.set_xticks([x_to_num[v] for v in x_labels])
    ax.set_xticklabels(x_labels, rotation=90, horizontalalignment='right')
    ax.set_yticks([y_to_num[v] for v in y_labels])
    ax.set_yticklabels(y_labels)
    
    ax.grid(False, 'major')
    ax.grid(True, 'minor')
    ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
    ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)
    
    ax.set_xlim([-0.5, max([v for v in x_to_num.values()]) + 0.5]) 
    ax.set_ylim([-0.5, max([v for v in y_to_num.values()]) + 0.5])
    return ax 