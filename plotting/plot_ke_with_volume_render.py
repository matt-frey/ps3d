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
cmap_name = 'rainbow4'

rainbow4_cmap = cc.cm[cmap_name]

rainbow4 = mpl_to_plotly(rainbow4_cmap, rainbow4_cmap.N)

def get_screenshot(step):
    ncreader = nc_reader()
    ncreader.open(os.path.join(filepath, 'beltrami_' + str(grid) + '_fields.nc'))
    t = ncreader.get_all('t')
    image = make_volume_rendering(ncr=ncreader, step=step, field='vorticity_magnitude')
    ncreader.close()
    return t[step], image


def add_inset(ax, bounds, t, ke, step):

    tref, arr = get_screenshot(step)

    idx = find_nearest(t=t, tref=tref)

    axins = ax.inset_axes(bounds=bounds, zorder=-1)
    ax.indicate_inset(bounds=[t[idx], ke[idx], 0.0, 0.0],
                    inset_ax=axins, zorder=1,
                    edgecolor='darkgrey', alpha=1.0)
    axins.imshow(arr)
    add_timestamp(axins, t[idx], xy=(0.05, 0.9), fmt="%.1f", fontsize=12)

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
