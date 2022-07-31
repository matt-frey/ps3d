import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create cross section difference figures.')
parser.add_argument('--filename',
                    type=str,
                    help='output file')
parser.add_argument('--step',
                    type=int,
                    default=5,
                    help='step to consider')

parser.add_argument('--plane',
                    type=str,
                    help="'xy', 'xz' or 'yz'")

parser.add_argument('--loc',
                    type=int,
                    default=0,
                    help='index location to extract plane')

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

parser.add_argument('--zlabel',
                    default=r'$z = -\pi/2$',
                    help='z-label for location')

args = parser.parse_args()
fname = args.filename
step = args.step
plane = args.plane
loc = args.loc
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum
zlabel = args.zlabel

print()
print("\tFilename:  ", fname)
print("\tStep:      ", step)
print("\tPlane:     ", plane)
print("\tLocation:  ", loc)
print("\tSave path: ", save_path)
print("\tOverwrite: ", overwrite)
print("\tFignum:   ", fignum)
print()

ncreader = nc_reader()
ncreader.open(fname)

t = ncreader.get_all('t')

fields = ['x_velocity',
          'y_velocity',
          'z_velocity',
          'x_vorticity',
          'y_vorticity',
          'z_vorticity']

labels = [r'$u(t) - u(0)$',
          r'$v(t) - v(0)$',
          r'$w(t) - w(0)$',
          r'$\xi(t) - \xi(0)$',
          r'$\eta(t) - \eta(0)$',
          r'$\zeta(t) - \zeta(0)$']

fig = plt.figure(figsize=(8, 5), dpi=200)
grid = ImageGrid(fig, 111,
                nrows_ncols=(2, 3),
                aspect=True,
                axes_pad=(0.45, 0.3),
                direction='row',
                share_all=True,
                cbar_location="right",
                cbar_mode='each',
                cbar_size="4%",
                cbar_pad=0.05)

for i, field in enumerate(fields):
    fdata = ncreader.get_dataset(step=step, name=field)
    idata = ncreader.get_dataset(step=0, name=field)
    fdata = fdata - idata

    ax = grid[i]
    im, cbar = make_imshow(ax=ax,
                        plane=plane,
                        loc=loc,
                        fdata=fdata,
                        ncr=ncreader,
                        cmap='rainbow4',
                        colorbar=True)

    if i < 3:
        remove_xticks(ax)

    if i == 0 or i == 3:
        pass
    else:
        remove_yticks(ax)

    add_annotation(ax, labels[i], xy=(0.03, 1.06), fmt="%.2f")

add_timestamp(grid[2], t[step], xy=(0.03, 1.2), fontsize=12)

if not zlabel is None:
    add_annotation(grid[2], zlabel, xy=(0.6, 1.2), fontsize=12)

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()

ncreader.close()
