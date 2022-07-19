import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create figure xy')
parser.add_argument('--filename',
                    type=str,
                    help='output file')
parser.add_argument('--steps',
                    type=int,
                    nargs=6,
                    help='add 6 steps to plot')

parser.add_argument('--plane',
                    type=str,
                    help="'xy', 'xz' or 'yz'")

parser.add_argument('--loc',
                    type=int,
                    help='index location to extract plane')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--auto_num',
                    help='auto select figure number',
                    action='store_true')

args = parser.parse_args()
fname = args.filename
steps = np.asarray(args.steps)
plane = args.plane
loc = args.loc
save_path = args.save_path
auto_num = args.auto_num

print("Filename: ", fname)
print("Steps:    ", steps)
print("Plane:    ", plane)
print("Location: ", loc)
print("Save path:", save_path)
print("Auto num: ", auto_num)

ncreader = nc_reader()
ncreader.open(fname)

t = ncreader.get_all('t')

#
# Vorticity magnitude
#
fig = plt.figure(figsize=(8, 5), dpi=200)
grid = ImageGrid(fig, 111,
                 nrows_ncols=(2, 3),
                 aspect=True,
                 axes_pad=(0.4, 0.3),
                 direction='row',
                 share_all=True,
                 cbar_location="right",
                 cbar_mode='each',
                 cbar_size="4%",
                 cbar_pad=0.1)

for i, step in enumerate(steps):
    vor = ncreader.get_dataset(step=step, name='vorticity_magnitude')

    ax = grid[i]
    im, cbar = make_imshow(ax=ax,
                           plane=plane,
                           loc=loc,
                           fdata=vor,
                           ncr=ncreader,
                           cmap='rainbow4',
                           colorbar=True)

    if i < 3:
        remove_xticks(ax)

    if i == 0 and i == 3:
        pass
    else:
        remove_yticks(ax)

    add_timestamp(ax, t[step], xy=(0.03, 1.06), fmt="%.2f")

add_annotation(grid[2], r'$z = -\pi/2$', xy=(0.6, 1.2), fontsize=12)

save_figure(plt=plt, figpath=save_path, fignum=1, auto=auto_num)
plt.close()

#
# Pressure
#
fig = plt.figure(figsize=(8, 5), dpi=200)
grid = ImageGrid(fig, 111,
                 nrows_ncols=(2, 3),
                 aspect=True,
                 axes_pad=(0.4, 0.3),
                 direction='row',
                 share_all=True,
                 cbar_location="right",
                 cbar_mode='each',
                 cbar_size="4%",
                 cbar_pad=0.1)

for i, step in enumerate(steps):
    pres = ncreader.get_dataset(step=step, name='pressure')

    ax = grid[i]

    im, cbar = make_imshow(ax=ax,
                           plane='xy',
                           loc=loc,
                           fdata=pres,
                           ncr=ncreader,
                           cmap='rainbow4',
                           colorbar=True)
    if i < -1:
        remove_xticks(ax)

    if i == 0 and i == 3:
        pass
    else:
        remove_yticks(ax)

    add_timestamp(ax, t[step], xy=(0.03, 1.06), fmt="%.2f")

add_annotation(grid[2], r'$z = -\pi/2$', xy=(0.6, 1.2), fontsize=12)

save_figure(plt=plt, figpath=save_path, fignum=2, auto=auto_num)

ncreader.close()
