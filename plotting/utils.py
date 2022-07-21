import matplotlib.pyplot as plt
import colorcet as cc
import matplotlib as mpl
import numpy as np
import os
import plotly.graph_objects as go

mpl.rcParams.update({
    "figure.figsize": (9, 6),
    "figure.dpi": 200,
    "font.family": "serif",
    "font.size": 20,
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

mpl.rcParams['font.size'] = 10

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

# 13 July 2022
# https://github.com/plotly/plotly.py/issues/2189
def mpl_to_plotly(cmap, pl_entries=11, rdigits=2):
    # cmap - colormap
    # pl_entries - int = number of Plotly colorscale entries
    # rdigits - int -=number of digits for rounding scale values
    scale = np.linspace(0, 1, pl_entries)
    colors = (cmap(scale)[:, :3]*255).astype(np.uint8)
    pl_colorscale = [[round(s, rdigits), f'rgb{tuple(color)}'] for s, color in zip(scale, colors)]
    return pl_colorscale

def make_volume_rendering(ncr, step, field, **kwargs):
    show_axes = kwargs.pop('show_axes', True)
    show_colorbar = kwargs.pop('show_colorbar', True)
    fmt = kwargs.pop('fmt', 'png')
    width = kwargs.pop('width', 1024)
    height = kwargs.pop('height', 1024)
    opacity = kwargs.pop('opacity', 0.01)
    surface_count = kwargs.pop('surface_count', 150)
    cmap_name = kwargs.pop('cmap_name', 'rainbow4')
    scale = kwargs.pop('scale', 1.5)
    margin = kwargs.pop('margin', dict(r=10, b=5, l=10, t=5))

    X, Y, Z = ncr.get_meshgrid()
    vor = ncr.get_dataset(step=step, name=field)

    rainbow4_cmap = cc.cm[cmap_name]
    rainbow4 = mpl_to_plotly(rainbow4_cmap, rainbow4_cmap.N)

    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=vor.flatten(),
        opacity=opacity, # needs to be small to see through all surfaces
        surface_count=surface_count, # needs to be a large number for good volume rendering
        colorscale=rainbow4,
        colorbar=dict(
                #title=dict(text='$|\omega|$', side='top'),
                thickness = 20,
                len=0.75,
                xpad=5,
                orientation='v',
                tickfont=dict(
                    size=22,
                    color="black"
                ),
            ),
        showscale=show_colorbar
        ))


    tickvals = np.array([-1.5, -0.75, 0.0, 0.75, 1.5])
    ticktext = [' -3/2 ', ' -3/4 ', ' 0 ', ' 3/4 ', ' 3/2 ']

    # 12 July 2022
    # https://plotly.com/python/3d-camera-controls/
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=-0.1),
        eye=dict(x=1.5, y=1.5, z=1.2)
    )

    fig.update_layout(
        scene_camera=camera,
        scene = dict(
                        xaxis = dict(
                            tickmode = 'array',
                            tickvals = tickvals,
                            ticktext = ticktext,
                            visible  = show_axes
                        ),
                        yaxis = dict(
                            tickmode = 'array',
                            tickvals = tickvals,
                            ticktext = ticktext,
                            visible  = show_axes
                        ),
                        zaxis = dict(
                            tickmode = 'array',
                            tickvals = tickvals,
                            ticktext = ticktext,
                            visible  = show_axes
                        ),
                        xaxis_title=dict(text=r'x', font=dict(size=24, color='black')),
                        yaxis_title=dict(text=r'y', font=dict(size=24, color='black')),
                        zaxis_title=dict(text=r'z', font=dict(size=24, color='black'))),
        margin=margin,
        font=dict(
            family="Arial", # Times New Roman
            size=16,
            color="black"
        ),
    )

    fig.write_image("temp_figure." + fmt, scale=scale, width=width, height=height)
    image = plt.imread("temp_figure." + fmt, format=fmt)
    return image
