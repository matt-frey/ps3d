import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib as mpl
import colorcet as cc
from mpl_toolkits.axes_grid1 import ImageGrid
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create vorticity and velocity height profile figures.')
parser.add_argument('--filenames',
                    type=str,
                    nargs='+',
                    help='output file')

parser.add_argument('--steps',
                    type=int,
                    nargs=7,
                    default=[0, 1, 2, 3, 4, 5, 6],
                    help='7 steps to plot')

parser.add_argument('--file_numbers',
                    type=int,
                    nargs=7,
                    default=[0, 0, 0, 0, 0, 0, 0],
                    help='file numbers')

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
file_numbers = args.file_numbers
steps = np.asarray(args.steps)
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

if fignum is None:
    print("No figure numbers provided.")
    exit()

print()
print("\tFilename:    ", fnames)
print("\tFile numbers:", file_numbers)
print("\tSteps:       ", steps)
print("\tSave path:   ", save_path)
print("\tOverwrite:   ", overwrite)
print("\tFignum:     ", fignum)
print()

ncreader = nc_reader()

#
# Vorticity and velocity rms profile
#
fig, axs = plt.subplots(2, 7, figsize=(9, 6), dpi=400, sharey=True)
grid = axs.flatten()

for i, step in enumerate(steps):
    ncreader.open(fnames[file_numbers[i]])
    t = ncreader.get_all('t')

    vel_ax = grid[i]
    vor_ax = grid[i+7]

    vel_labels = [None] * 3
    vor_labels = [None] * 3
    if i == 3:
        vel_labels = [
            r'$u_{\mathrm{rms}}$',
            r'$v_{\mathrm{rms}}$',
            r'$w_{\mathrm{rms}}$']
        vor_labels = [
            r'$\xi_{\mathrm{rms}}$',
            r'$\eta_{\mathrm{rms}}$',
            r'$\zeta_{\mathrm{rms}}$']

    if not i == 0 and not i == 7:
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
    
    add_timestamp(vel_ax, t[step], xy=(0.035, 1.06), fmt="%.2f")

grid[0].set_ylabel(r'$z$')
grid[7].set_ylabel(r'$z$')

grid[3].legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.35))
grid[7+3].legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.2))

plt.subplots_adjust(hspace=0.3, wspace=0.4)

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
