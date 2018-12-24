# -*- coding: utf-8 -*-
"""
This module allows for plotting of annotated graphs using matplotlib

"""

import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import numpy as np
Path=mpath.Path

def _draw_beta(ax,start,stop,label="",linewidth=1.2):
    '''
    This is a private function which draws a named
    arrow on a matplotlib axis frim start to stop 
    at 0.5 height
    '''
    ax.add_patch(mpatches.FancyArrow(start-0.5, 0.5, stop-start, 0,width=0.5, head_width=1, head_length=1,edgecolor='k',facecolor='k'))
    ax.text(start+(stop-start)/2,-1,label,style='italic',horizontalalignment= 'center', verticalalignment= 'baseline')



def _draw_helix(ax,start,stop,label="",linewidth=1.2):
    '''
    This is a private function which draws a named
    helix on a matplotlib axis frim start to stop 
    at 0.5 height
    '''
    ax.add_patch(mpatches.PathPatch(Path([(start-0.5,0),(start,-0.15),(start+0.5,0)],
                                          [Path.MOVETO,Path.CURVE3,Path.CURVE3]),
                 fc="none", edgecolor='k',linewidth=linewidth))
    for pos in np.arange(start,stop):
        ax.add_patch(mpatches.Ellipse((pos+0.5, 0.35), 0.4, 0.7, angle=-20, linewidth=linewidth,edgecolor='k', fill=False, zorder=2))
        ax.add_patch(mpatches.PathPatch(Path([(pos+0.5,0),(pos+1,-0.15),(pos+1.5,0)],
                                             [Path.MOVETO,Path.CURVE3,Path.CURVE3]),
                     fc="none", edgecolor='k',linewidth=linewidth))
    ax.text(start+(stop-start)/2,-1,label,style='italic',horizontalalignment= 'center', verticalalignment= 'baseline')
    
def _annotate_sequence(dataax, annot_ax,resids,features):
    for feature in features:
        start=feature['sel'][0]+resids[0]
        stop=feature['sel'][1]+resids[0]
        if feature['style']=='helix':
            _draw_helix(annot_ax,start,stop,feature['text'])
        elif feature['style']=='beta':
            _draw_beta(annot_ax,start,stop,feature['text'])
        elif feature['style']=='block':
            label = dataax.xaxis.get_ticklabels(minor=True)[feature['sel'][0]]
            label.set_bbox(dict(facecolor='none', edgecolor='red',boxstyle='square,pad=0.1'))
            
def plot_on_seq(data,seq,**kwargs):
    '''
    this function creates a new figure with sequence annotated bar chart.
    data - 1d numpy array with per residue values
    seq  - string with sequence
    
    available keyword arguments:
    filename - name to save file (png, svg, pdf and other formats available)
    features - dictionary with sequence features from 
               seq_tools.hist_ss.get_hist_ss_in_aln_for_shade
    resids   - np array with residue numbers of same length as data and seq
    y_axis_label - string label of Y axis
        figsize  - tupple with figure size in inches (width,heigth)
    dpi      - int with DPI value
    '''
    
    # Trying to predict 'optimal' width for the plot
    if not 'figsize' in kwargs:
        width=len(seq)/5
        heigth=width/10
        figsize=(width,heigth)
    else:
        figsize=kwargs['figsize']
    if not 'dpi' in kwargs:
        dpi=300
    else:
        dpi=kwargs['dpi']
        
    # Creating subplot grid (upper for data, lower for annotation)
    gs  = gridspec.GridSpec(2, 1, height_ratios=[1, 0.5])    
    fig = plt.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1],sharex=ax1)    
    
    # Populating values if abscent
    if not 'resids' in kwargs:
        resids=np.arange(data.size)
    else:
        resids=kwargs['resids']
        
    # Plotting the data
    ax1.bar(resids,data)    
    ax1.xaxis.set_ticks(resids,minor=True)
    ax1.xaxis.set_ticklabels(seq,minor=True)
    if 'y_axis_label' in kwargs:
        ax1.set_ylabel(kwargs['y_axis_label'])
    ax1.grid(True,'major')
    ax1.tick_params(axis=u'both', which=u'major',length=0,labeltop=True)  
    ax1.tick_params(axis=u'both', which=u'major',length=0,labelbottom=False) 
    ax1.tick_params(axis=u'both', which=u'minor',length=0) 
    
    # Preparing 2nd axis for annotating
    ax2.set_aspect('equal')
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_ylim(-0.5,1.1)
    
    # Hiding all elements in the axes, except plot itself
    ax2.axis('off')
    if 'features' in kwargs:
        _annotate_sequence(ax1,ax2,resids,kwargs['features'])
    
    if 'filename' in kwargs:
        fig.savefig(kwargs['filename'])
    plt.show()
    return(fig)

def heatplot_on_seq(data,seq,**kwargs):
    ''' 
    this function creates a new figure with sequence annotated bar chart.
    data - 2d numpy array with per residue values
           cols - resids, rows - values
    seq  - string with sequence
    
    available keyword arguments:
    filename - name to save file (png, svg, pdf and other formats available)
    features - dictionary with sequence features from 
               seq_tools.hist_ss.get_hist_ss_in_aln_for_shade
    resids   - np array with residue numbers of same length as data and seq
    y_axis_label - string label of Y axis
    y_axis_values  - 1d numpy array with row names
    colorbar_label - string label of colorbar axis
    cmap - matplotlib colormap
    figsize  - tupple with figure size in inches (width,heigth)
    dpi      - int with DPI value
    '''
    
    # Try to predict 'optimal' width for the plot
    if not 'figsize' in kwargs:
        width=len(seq)/5
        heigth=width/5
        figsize=(width,heigth)
    else:
        figsize=kwargs['figsize']
    if not 'dpi' in kwargs:
        dpi=300
    else:
        dpi=kwargs['dpi']
    factor=(figsize[0]/figsize[0])/15
    
    # Creating subplot grid (upper for data, lower for annotation)
    gs  = gridspec.GridSpec(2, 1, height_ratios=[1, factor])
    fig = plt.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0],sharex=ax1)
    
    # Populating values if abscent
    if not 'resids' in kwargs:
        resids=np.arange(data.size)
    else:
        resids=kwargs['resids']
    if not 'y_axis_values' in kwargs:
        y_names=np.arange(data.shape[0])
    else:
        y_names=kwargs['y_axis_values']
    if not 'cmap' in kwargs:
        cmap='viridis'
    else:
        cmap=kwargs['cmap']
    
    # Plotting ang adjusting the heatmap
    im=ax1.imshow(data,extent=(resids[0]-0.5,resids[-1]+0.5,y_names[0],y_names[-1]),aspect = 'auto',cmap=cmap)
    
    # Adding all ticks and their labels
    ax1.xaxis.set_ticks(resids,minor=True)
    ax1.xaxis.set_ticklabels(seq,minor=True)
    
    if 'y_axis_label' in kwargs:
        ax1.set_ylabel(kwargs['y_axis_label'])
    
    # Tinkering with the grid and ticks
    ax1.grid(True,'major')
    ax1.tick_params(axis=u'both', which=u'major',length=0,labeltop=True)  
    ax1.tick_params(axis=u'both', which=u'major',length=0,labelbottom=False) 
    ax1.tick_params(axis=u'both', which=u'minor',length=0) 
    
    # Preparing 2nd axis for annotating
    ax2.set_aspect('equal')
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_ylim(-0.5,1.1)
    
    # Hiding all elements in the axes, except plot itself
    ax2.axis('off')
    
    # Trying to tighten the graph, although I do not think it will work with gridspec subplots
    fig.tight_layout()
    
    # 'Uninvasive' addition of the colorbar
    box = ax1.get_position()
    cb1=plt.colorbar(im,ax=ax1,orientation="vertical")    
    ax1.set_position(box)    
    cb1.ax.set_position([box.x0*1.02 + box.width * 1.02, box.y0, 0.02, box.height])
    if 'colorbar_label' in kwargs:
        cb1.ax.set_ylabel(kwargs['colorbar_label'])    
    
    # Annotating the sequence at the 2nd axes
    if 'features' in kwargs:
        _annotate_sequence(ax1,ax2,resids,kwargs['features'])
    
    # Optional saving
    if 'filename' in kwargs:
        fig.savefig(kwargs['filename'])
    plt.show()
    return(fig)
