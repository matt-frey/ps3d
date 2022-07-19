import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import plotly.graph_objects as go
from nc_reader import nc_reader
import colorcet as cc
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create KE evolution including volume rendering.')
parser.add_argument('--filepath',
                    type=str,
                    help='file path')

parser.add_argument('--grid',
                    type=int,
                    help='grid')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--overwrite',
                    help='overwrite figures',
                    action='store_true')

parser.add_argument('--fignum',
                    type=int,
                    help='figure number')

args = parser.parse_args()
filepath = args.filepath
grid = args.grid
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

print()
print("\tFilepath:  ", filepath)
print("\tGrid:      ", grid)
print("\tSave path: ", save_path)
print("\tOverwrite: ", overwrite)
print("\tFignum:    ", fignum)
print()

mpl.rcParams['font.size'] = 16

# 12 July 2022
# https://plotly.com/python/v3/matplotlib-colorscales/#formatting-the-colormap
def matplotlib_to_plotly(cmap, pl_entries):
    h = 1.0/(pl_entries-1)
    pl_colorscale = []

    for k in range(pl_entries):
        #C = map(np.uint8, np.array(cmap(k*h)[:3])*255)
        C = np.array(cmap(k*h)[:3])*255
        pl_colorscale.append([k*h, 'rgb'+str((C[0], C[1], C[2]))])

    return pl_colorscale

rainbow4 = matplotlib_to_plotly(cc.cm['rainbow4'], 256)

def get_screenshot(step):
    ncreader = nc_reader()
    ncreader.open(os.path.join(filepath, 'beltrami_' + str(grid) + '_fields.nc'))
    t = ncreader.get_all('t')
    X, Y, Z = ncreader.get_meshgrid()
    vor = ncreader.get_dataset(step=step, name='vorticity_magnitude')
    ncreader.close()

    opacity = 0.01
    surface_count = 150


    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=vor.flatten(),
        #isomin=0.01,
        #isomax=0.99,
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
            )
        ))


    tickvals = np.array([-1.5, -0.75, 0.0, 0.75, 1.5])
    ticktext = [' -3/2 ', ' -3/4 ', ' 0 ', ' 3/4 ', ' 3/2 ']

    # 12 July 2022
    # https://plotly.com/python/3d-camera-controls/
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=-0.1),
        eye=dict(x=1.4, y=1.4, z=1.4)
    )


    fig.update_layout(
        scene_camera=camera,
        scene = dict(
                        xaxis = dict(
                            tickmode = 'array',
                            tickvals = tickvals,
                            ticktext = ticktext
                        ),
                        yaxis = dict(
                            tickmode = 'array',
                            tickvals = tickvals,
                            ticktext = ticktext
                        ),
                        zaxis = dict(
                            tickmode = 'array',
                            tickvals = tickvals,
                            ticktext = ticktext
                        ),
                        xaxis_title=dict(text=r'x', font=dict(size=24, color='black')),
                        yaxis_title=dict(text=r'y', font=dict(size=24, color='black')),
                        zaxis_title=dict(text=r'z', font=dict(size=24, color='black'))),
                        margin=dict(r=10, b=5, l=10, t=5),
        font=dict(
            family="Arial", # Times New Roman
            size=16,
            color="black"
        ),
    )


    fig.write_image("temp_fig.png", scale=1.5, width=1024, height=1024)

    image = plt.imread("temp_fig.png", format='png')

    return t[step], image

# 7 July 2022
# https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest(t, tref):
    idx = (np.abs(tref - t)).argmin()
    print('Closest time to', tref, 'is', t[idx])
    return idx



def add_inset(ax, bounds, t, ke, step):

    tref, arr = get_screenshot(step)

    idx = find_nearest(t=t, tref=tref)

    axins = ax.inset_axes(bounds=bounds, zorder=-1)
    ax.indicate_inset(bounds=[t[idx], ke[idx], 0.0, 0.0],
                    inset_ax=axins, zorder=1,
                    edgecolor='darkgrey', alpha=1.0)
    axins.imshow(arr)
    add_timestamp(axins, t[idx], xy=(0.05, 0.9), fmt="%.2f")

    # 7 July 2022
    # https://stackoverflow.com/questions/37039685/hide-tick-label-values-but-keep-axis-labels
    axins.set_xticks([])
    axins.set_yticks([])

    # 7 July 2022
    # https://stackoverflow.com/questions/1982770/matplotlib-changing-the-color-of-an-axis
    axins.spines['bottom'].set_color('darkgrey')
    axins.spines['top'].set_color('darkgrey')
    axins.spines['right'].set_color('darkgrey')
    axins.spines['left'].set_color('darkgrey')


plt.figure(figsize=(13, 6.5), dpi=200)
ax = plt.gca()

t, ke, en = np.loadtxt(os.path.join(filepath, 'beltrami_' + str(grid) + '_ecomp.asc'),
                       skiprows=1, unpack=True)

ncelli = 1.0 / grid ** 3

# calculate mean KE and mean EN
ke *= ncelli
ax.plot(t, ke, color='blue', linewidth=1.0)
ax.set_xlabel(r'time, $t$')
ax.set_ylabel(r'average kinetic energy, $\langle\mathcal{K}\rangle$')

w = 0.45
h = w
add_inset(ax=ax, bounds=[-0.1, 0.47, w, h], t=t, ke=ke, step=0)
add_inset(ax=ax, bounds=[ 0.17, 0.47, w, h], t=t, ke=ke, step=6)
add_inset(ax=ax, bounds=[ 0.17, 0.01, w, h], t=t, ke=ke, step=7)
add_inset(ax=ax, bounds=[ 0.6, 0.47, w, h], t=t, ke=ke, step=10)

plt.tight_layout()

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
