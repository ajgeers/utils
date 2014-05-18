"""Module for plots."""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def colorbar(path='colorbar.pdf', vmin=0, vmax=1, unit='[-]',
             colors='RdBu_r', nlevels=None, orientation='horizontal',
             extend='neither'):
    """Create colorbar as separate figure and write to disk.

    Args:
        vmin: Lower limit of colorbar range.
        vmax: Upper limit of colorbar range.
        unit: String to label colorbar.
        colors: Colormap (see matplotlib website for options).
        nlevel: Number of colors for discrete colorbar; 'None' gives continuous
            colorbar.
        orientation ('horizontal', 'vertical'): Orientation of colorbar.
        extend ('neither', 'both', 'min', 'max'): Extend colorbar to show
            out-of-range values.

    """

    if orientation == 'horizontal':
        figsize = (9, .5)
    else:
        figsize = (.5, 9)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    # Initiate colormap
    cmap = mpl.cm.get_cmap(colors)

    # Define norm to correspond to the data for which the colorbar will be
    # used. If nlevels is given, the colorbar will be discrete.
    if nlevels == None:
        norm = mpl.colors.Normalize(vmin, vmax)
    else:
        bounds = np.linspace(vmin, vmax, nlevels)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    # ColorBarBase derives from ScalarMappable and puts a colorbar in a
    # specified axes
    cbar = mpl.colorbar.ColorbarBase(ax,
                                     cmap=cmap,
                                     norm=norm,
                                     orientation=orientation,
                                     ticks=[vmin, (vmax+vmin) / 2., vmax],
                                     extend=extend)

    # Set font properties and align ticks within colorbar edges
    if orientation == 'horizontal':
        ax.set_xlabel(unit, fontsize=36)
        ax.xaxis.set_tick_params(labelsize=32)
        ax.get_xticklabels()[0].set_horizontalalignment('left')
        ax.get_xticklabels()[2].set_horizontalalignment('right')
    else:
        ax.set_xlabel(unit, fontsize=36)
        ax.yaxis.set_tick_params(labelsize=32)
        ax.get_yticklabels()[0].set_verticalalignment('bottom')
        ax.get_yticklabels()[2].set_verticalalignment('top')

    # Hack to remove white gaps in Preview
    cbar.solids.set_edgecolor("face")

    fig.savefig(path, bbox_inches='tight', pad_inches=0.05)
