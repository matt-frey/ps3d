import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib as mpl
import colorcet as cc
from mpl_toolkits.axes_grid1 import ImageGrid
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create vorticity height profile figures.')
parser.add_argument('--filenames',
                    type=str,
                    nargs='+',
                    help='output file')

parser.add_argument('--steps',
                    type=int,
                    nargs=9,
                    default=[0, 1, 2, 3, 4, 5, 6, 7, 8],
                    help='9 steps to plot')

parser.add_argument('--file_numbers',
                    type=int,
                    nargs=9,
                    default=[0, 0, 0, 0, 0, 0, 0, 0, 0],
                    help='file numbers')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--overwrite',
                    help='overwrite figures',
                    action='store_true')

parser.add_argument('--fignums',
                    type=int,
                    nargs=2,
                    help='figure numbers')

args = parser.parse_args()
fnames = args.filenames
file_numbers = args.file_numbers
steps = np.asarray(args.steps)
save_path = args.save_path
overwrite = args.overwrite
fignums = args.fignums

if fignums is None:
    print("No figure numbers provided.")
    exit()

print()
print("\tFilename:    ", fnames)
print("\tFile numbers:", file_numbers)
print("\tSteps:       ", steps)
print("\tSave path:   ", save_path)
print("\tOverwrite:   ", overwrite)
print("\tFignums:     ", fignums)
print()

ncreader = nc_reader()

#
# Vorticity mean profile
#
fig = plt.figure(figsize=(9, 3), dpi=400)
grid = ImageGrid(fig, 111,
                 nrows_ncols=(1, 9),
                 aspect=True,
                 axes_pad=(0.1, 0.3),
                 direction='row',
                 share_all=True,
                 cbar_location="bottom",
                 cbar_mode=None,
                 cbar_size="4%",
                 cbar_pad=0.0)

for i, step in enumerate(steps):

    ncreader.open(fnames[file_numbers[i]])
    t = ncreader.get_all('t')
    
    ax = grid[i]

    labels = [None] * 3
    if i == 4:
        labels = [
            r'$\langle\xi\rangle$',
            r'$\langle\eta\rangle$',
            r'$\langle\zeta\rangle$']

    if not i == 0:
        remove_yticks(ax)

    make_mean_profiles(ax=ax,
                       ncr=ncreader,
                       step=step,
                       fields=['x_vorticity', 'y_vorticity', 'z_vorticity'],
                       labels=labels)

    ncreader.close()

    add_timestamp(ax, t[step], xy=(0.03, 1.06), fmt="%.2f")

grid[0].set_ylabel(r'$z$')
grid[5].set_ylabel(r'$z$')

grid[4].legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.4))

save_figure(plt=plt, figpath=save_path, fignum=fignums[0], overwrite=overwrite)
plt.close()

#
# Vorticity rms profile
#
#fig = plt.figure(figsize=(9, 3), dpi=400)
fig, axs = plt.subplots(2, 9, figsize=(9, 6), dpi=400, sharey=True)
grid = axs.flatten()
#grid = ImageGrid(fig, 111,
#                 nrows_ncols=(1, 9),
#                 aspect=True,
#                 axes_pad=(0.1, 0.3),
#                 direction='row',
#                 share_all=True,
#                 cbar_location="bottom",
#                 cbar_mode=None,
#                 cbar_size="4%",
#                 cbar_pad=0.0)

for i, step in enumerate(steps):
    ncreader.open(fnames[file_numbers[i]])
    t = ncreader.get_all('t')

    vel_ax = grid[i]
    vor_ax = grid[i+9]

    vel_labels = [None] * 3
    vor_labels = [None] * 3
    if i == 4:
        vel_labels = [
            r'$u_{\mathrm{rms}}$',
            r'$v_{\mathrm{rms}}$',
            r'$w_{\mathrm{rms}}$']
        vor_labels = [
            r'$\xi_{\mathrm{rms}}$',
            r'$\eta_{\mathrm{rms}}$',
            r'$\zeta_{\mathrm{rms}}$']

    if not i == 0 and not i == 9:
        remove_yticks(vel_ax)
        remove_yticks(vor_ax)

    make_rms_profiles(ax=vel_ax,
                      ncr=ncreader,
                      step=step,
                      fields=['x_velocity', 'y_velocity', 'z_velocity'],
                      labels=vel_labels)

    make_rms_profiles(ax=vor_ax,
                      ncr=ncreader,
                      step=step,
                      fields=['x_vorticity', 'y_vorticity', 'z_vorticity'],
                      labels=vor_labels)

    ncreader.close()
    
    add_timestamp(vel_ax, t[step], xy=(0.03, 1.06), fmt="%.2f")

grid[0].set_ylabel(r'$z$')
grid[9].set_ylabel(r'$z$')

grid[4].legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.35))
grid[9+4].legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.2))

plt.subplots_adjust(hspace=0.3)

save_figure(plt=plt, figpath=save_path, fignum=fignums[1], overwrite=overwrite)
plt.close()
