import matplotlib.pyplot as plt
import colorcet as cc
import matplotlib as mpl
import numpy as np

mpl.rcParams.update({
    "figure.figsize": (9, 6),
    "figure.dpi": 200,
    "font.family": "serif",
    "font.size": 20,
    "text.usetex": True,
    'text.latex.preamble': "\n".join([
        r"\usepackage{amsmath}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage{bm}"
#        r"\usepackage{siunitx}",
        ])
})

mpl.rcParams['font.size'] = 10

def add_timestamp(plt, time, xy=(0.75, 1.05), fmt="%.2f"):
    bbox = dict(boxstyle="round", facecolor="wheat", edgecolor='none')
    label = r"t = " + fmt % (time)
    plt.annotate(label, xy=xy, xycoords="axes fraction", bbox=bbox)

def add_annotation(ax, label, xy, **kwargs):
    bbox = dict(boxstyle="round", facecolor="wheat", edgecolor='none')
    ax.annotate(label, xy=xy, xycoords="axes fraction", bbox=bbox)

def remove_xticks(ax):
    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

def remove_yticks(ax):
    ax.tick_params(axis='y', which='both', right=False, left=False)

# assumes fdata ordering (x, y, z)
def make_imshow(ax, plane, loc, fdata, ncr, cmap='rainbow4', colorbar=True):
    origin = ncr.get_box_origin()
    extent = ncr.get_box_extent()

    # plane = 'xy':
    imin = origin[0]
    imax = origin[0] + extent[0]
    jmin = origin[1]
    jmax = origin[1] + extent[1]


    if plane == 'xy':
        pl = fdata[:, :, loc]
        xlab = r'$x$'
        ylab = r'$y$'
    if plane == 'xz':
        jmin = origin[2]
        jmax = origin[2] + extent[2]
        pl = fdata[:, loc, :]
        xlab = r'$x$'
        ylab = r'$z$'
    elif plane == 'yz':
        imin = jmin
        imax = jmax
        jmin = origin[2]
        jmax = origin[2] + extent[2]
        pl = fdata[loc, :, :]
        xlab = r'$y$'
        ylab = r'$z$'

    im = ax.imshow(X=pl.transpose(),
                   cmap=cc.cm[cmap],
                   interpolation='bilinear',
                   origin='lower',
                   extent=[imin, imax, jmin, jmax])

    ticks   = np.pi * np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
    ticklab = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']
    ax.set_xticks(ticks, ticklab)
    ax.set_yticks(ticks, ticklab)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    cbar = None
    if colorbar:
        cbar = ax.cax.colorbar(im)
        cbar.formatter.set_powerlimits((0, 0))
        # 18 July 2022
        # https://stackoverflow.com/questions/34039396/matplotlib-colorbar-scientific-notation-offset
        cbar.ax.yaxis.set_offset_position('left')
    return im, cbar
