import numpy as np
import matplotlib.pyplot as plt
from tools.nc_reader import nc_reader
import matplotlib as mpl
import colorcet as cc
from mpl_toolkits.axes_grid1 import ImageGrid
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create kinetic energy and enstrophy height profile figures.')
parser.add_argument('--filename',
                    type=str,
                    help='output file')
parser.add_argument('--steps',
                    type=int,
                    nargs=10,
                    default=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                    help='add 6 steps to plot')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--overwrite',
                    help='overwrite figures',
                    action='store_true')

parser.add_argument('--fignum',
                    type=int,
                    help='figure number',
                    default=12)


args = parser.parse_args()
fname = args.filename
steps = np.asarray(args.steps)
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

print()
print("\tFilename:  ", fname)
print("\tSteps:     ", steps)
print("\tSave path: ", save_path)
print("\tOverwrite: ", overwrite)
print("\tFignum:   ", fignum)
print()

ncreader = nc_reader()
ncreader.open(fname)
t = ncreader.get_all('t')

fig = plt.figure(figsize=(len(steps)*1.0, 5), dpi=400)
grid = ImageGrid(fig, 111,
                 nrows_ncols=(2, 5),
                 aspect=True,
                 axes_pad=(0.1, 0.3),
                 direction='row',
                 share_all=True,
                 cbar_location="bottom",
                 cbar_mode=None,
                 cbar_size="4%",
                 cbar_pad=0.0)

for i, step in enumerate(steps):
    ax = grid[i]

    labels = [None] * 2
    if i == 2:
        labels = [
            r'$\bar{|\bm{u}|^2}/2$',
            r'$\bar{|\bm{\omega}|^2}/2$']
    if i < 5:
        remove_xticks(ax)

    if not i == 0 and not i == 5:
        remove_yticks(ax)

    make_mean_profiles(ax=ax,
                       ncr=ncreader,
                       step=step,
                       fields=['kinetic_energy', 'enstrophy'],
                       labels=labels)

    add_timestamp(ax, t[step], xy=(0.03, 1.06), fmt="%.2f")

grid[0].set_ylabel(r'$z$')
grid[5].set_ylabel(r'$z$')

grid[2].legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.4))

#plt.tight_layout()

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()

ncreader.close()
