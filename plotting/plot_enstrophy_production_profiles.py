import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib as mpl
import colorcet as cc
from mpl_toolkits.axes_grid1 import ImageGrid
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create enstrophy production height profile figures.')
parser.add_argument('--filenames',
                    type=str,
                    nargs=7,
                    help='output files')

parser.add_argument('--file_path',
                    type=str,
                    help='path to output files')

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
fnames = args.filenames
file_path = args.file_path
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

if fignum is None:
    print("No figure numbers provided.")
    exit()

print()
print("\tSave path:   ", save_path)
print("\tOverwrite:   ", overwrite)
print("\tFignum:     ", fignum)
print()

ncreader = nc_reader()

fig, axs = plt.subplots(1, 7, figsize=(9, 3), dpi=400, sharey=True)
grid = axs.flatten()

for i, fname in enumerate(fnames):

    ax = grid[i]
    
    labels = [None] * 3
    if i == 3:
        labels = [
            r'$\upsilon_1$',
            r'$\upsilon_2$',
            r'$\upsilon_3$']

    if not i == 0:
        remove_yticks(ax)
        
    ff = os.path.join(file_path, fname)
    with open(ff) as f:
        line = f.readline()
    t = float(line.split('# time')[1])
        
    z, up1, up2, up3 = np.loadtxt(ff, skiprows=2, unpack=True)

    ax.plot(up1, z, color=colors[0], label=labels[0])
    ax.plot(up2, z, color=colors[1], label=labels[1])
    ax.plot(up3, z, color=colors[2], label=labels[2])
    ax.grid(zorder=-1)
    
    add_timestamp(ax, t, xy=(0.035, 1.06), fmt="%.2f")

grid[3].legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.35))
grid[0].set_ylabel(r'$z$')

zticks   = np.pi * np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
zticklab = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']
ax.set_yticks(ticks=zticks, labels=zticklab)

plt.subplots_adjust(hspace=0.3, wspace=0.4)

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
