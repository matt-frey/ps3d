import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
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
                    nargs='+',
                    help='add steps to plot')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--fignum',
                    type=int,
                    help='figure number',
                    default=12)

parser.add_argument('--overwrite',
                    help='overwrite figures',
                    action='store_true')


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

n = len(steps)

fig = plt.figure(figsize=(n, 3), dpi=400)
grid = ImageGrid(fig, 111,
                 nrows_ncols=(1, n),
                 aspect=True,
                 axes_pad=(0.2, 0.3),
                 direction='row',
                 share_all=True,
                 cbar_location="bottom",
                 cbar_mode=None,
                 cbar_size="4%",
                 cbar_pad=0.0)

for i, step in enumerate(steps):
    ax = grid[i]

    labels = [None] * 2
    if i == int(n / 2):
        labels = [
            r'$\langle|\bm{u}|^2\rangle/2$',
            r'$\langle|\bm{\omega}|^2\rangle/2$']

    if not i == 0:
        remove_yticks(ax)

    make_mean_profiles(ax=ax,
                       ncr=ncreader,
                       step=step,
                       fields=['kinetic_energy', 'enstrophy'],
                       labels=labels,
                       normalise=True,
                       xticks=[0, 0.5, 1],
                       xlim=[-0.1, 1.1])

    add_timestamp(ax, t[step], xy=(0.03, 1.06), fmt="%.2f")

grid[0].set_ylabel(r'$z$')
grid[5].set_ylabel(r'$z$')

grid[int(n / 2)].legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.4))

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()

ncreader.close()
