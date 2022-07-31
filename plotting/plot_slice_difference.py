import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib.colors as mpl_colors
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

# r'$(v(t) - v(0)) / |\bm{u}(t) - \bm{u}(0)|_{\mathrm{\max}}$',
labels = [r'$\Delta u(t) / |\Delta\bm{u}(t))|_{\mathrm{\max}}$',
          r'$\Delta v(t) / |\Delta\bm{u}(t)|_{\mathrm{\max}}$',
          r'$\Delta w(t) / |\Delta\bm{u}(t)|_{\mathrm{\max}}$',
          r'$\Delta\xi(t) / |\Delta\bm{u}(t)|_{\mathrm{\max}}$',
          r'$\Delta\eta(t) / |\Delta\bm{u}(t)|_{\mathrm{\max}}$',
          r'$\Delta\zeta(t) / |\Delta\bm{u}(t)|_{\mathrm{\max}}$']

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

# find normalization: we normalize the differences by the maximum
# magnitude of the velocity field difference:
udiff = ncreader.get_dataset(step=step, name='x_velocity') - \
    ncreader.get_dataset(step=0, name='x_velocity')

vdiff =	ncreader.get_dataset(step=step, name='y_velocity') - \
    ncreader.get_dataset(step=0, name='y_velocity')

wdiff =	ncreader.get_dataset(step=step, name='z_velocity') - \
    ncreader.get_dataset(step=0, name='z_velocity')

max_magn = get_max_magnitude(udiff, vdiff, wdiff, plane, loc)

for i, field in enumerate(fields):
    fdata = ncreader.get_dataset(step=step, name=field)
    idata = ncreader.get_dataset(step=0, name=field)

    fdata = (fdata - idata) / max_magn

    ax = grid[i]
    im, cbar = make_imshow(ax=ax,
                           plane=plane,
                           loc=loc,
                           fdata=fdata,
                           ncr=ncreader,
                           cmap='coolwarm',
                           norm=mpl_colors.CenteredNorm(vcenter=0.0),
                           colorbar=True)

    if i < 3:
        remove_xticks(ax)

    if i == 0 or i == 3:
        pass
    else:
        remove_yticks(ax)

    add_annotation(ax, labels[i], xy=(0.03, 1.06))

add_timestamp(grid[0], t[step], xy=(-0.18, 1.25), fontsize=12)

if not zlabel is None:
    add_annotation(grid[0], zlabel, xy=(-0.5, 1.25), fontsize=12)

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()

ncreader.close()
