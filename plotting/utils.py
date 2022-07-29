import matplotlib.pyplot as plt
import colorcet as cc
import matplotlib as mpl
import numpy as np
import os

mpl.rcParams.update({
    "figure.figsize": (9, 6),
    "figure.dpi": 200,
    "font.family": "serif",
    "font.size": 11,
    "text.usetex": True,
    'legend.framealpha': 1.0,
    'lines.linewidth': 0.75,
    'grid.color':     'b0b0b0',
    'grid.linestyle': 'dotted',
    'grid.linewidth': 0.25,
    'grid.alpha':     1.0,
    'text.latex.preamble': "\n".join([
        r"\usepackage{amsmath}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage{bm}"
#        r"\usepackage{siunitx}",
        ])
})

# 28 July 2022
# https://stackoverflow.com/questions/42086276/get-default-line-colour-cycle
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

def add_timestamp(plt, time, xy=(0.75, 1.05), fmt="%.2f", **kwargs):
    bbox = dict(boxstyle="round", facecolor="wheat", edgecolor='none')
    label = r"t = " + fmt % (time)
    plt.annotate(label, xy=xy, xycoords="axes fraction", bbox=bbox, **kwargs)

def add_annotation(ax, label, xy, **kwargs):
    bbox = dict(boxstyle="round", facecolor="wheat", edgecolor='none')
    ax.annotate(label, xy=xy, xycoords="axes fraction", bbox=bbox, **kwargs)

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

def make_mean_profiles(ax, ncr, step, fields, labels):
    z = ncr.get_all('z')
    n = len(z)
    zticks   = np.pi * np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
    zticklab = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']
    markersize = 3
    markers = ['o', 'x', '+']
    colors = ['blue', 'red', 'green']
    linewidth = 0.75

    for i, field in enumerate(fields):
        bar = np.zeros(n)
        data = ncr.get_dataset(step=step, name=field)
        bar = data.mean(axis=(1, 2))
        ax.plot(bar, z,
                color=colors[i],
                marker=markers[i],
                markersize=markersize,
                linewidth=linewidth,
                label=labels[i])
        ax.set_aspect(1)
        ax.set_yticks(ticks=zticks, labels=zticklab)
        ax.set_xticks([-1, 0, 1])
        ax.set_xlim([-1.1, 1.1])
    return ax

def make_rms_profiles(ax, ncr, step, fields, labels):
    z = ncr.get_all('z')
    n = len(z)
    zticks   = np.pi * np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
    zticklab = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']
    markersize = 3
    markers = ['o', 'x', '+']
    colors = ['blue', 'red', 'green']
    linewidth = 0.75

    for i, field in enumerate(fields):
        rms = np.zeros(n)
        data = ncr.get_dataset(step=step, name=field)
        rms = np.sqrt((data ** 2).mean(axis=(1, 2)))

        ax.plot(rms, z,
                color=colors[i],
                marker=markers[i],
                markersize=markersize,
                linewidth=linewidth,
                label=labels[i])
        ax.set_aspect(1)
        ax.set_yticks(ticks=zticks, labels=zticklab)
        ax.set_xticks([0, 1, 2, 3])
        ax.set_xlim([-0.1, 3.1])
    return ax


def save_figure(plt, figpath, fignum=1, overwrite=False):
    figname = 'fig' + str(fignum) + '.eps'
    fname = os.path.join(figpath, figname)

    if os.path.exists(fname) and not overwrite:
        print("Figure '" + fname + "' already exists.")
        plt.close()
        exit()

    print("Save figure as:", fname)
    plt.savefig(fname=fname, format='eps')
    plt.close()

# 7 July 2022
# https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest(t, tref):
    idx = np.argmin(np.abs(tref - t))
    return idx
