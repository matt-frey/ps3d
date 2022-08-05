import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
from mpl_toolkits.axes_grid1 import ImageGrid
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create cross section figures.')
parser.add_argument('--filenames',
                    type=str,
                    nargs='+',
                    help='output file')
parser.add_argument('--steps',
                    type=int,
                    nargs=6,
                    default=[0, 1, 2, 3, 4, 5],
                    help='add 9 steps to plot')

parser.add_argument('--file_numbers',
                    type=int,
                    nargs=6,
                    default=[0, 0, 0, 0, 0, 0],
                    help='the file by index to read data from')

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

parser.add_argument('--field',
                    help='field to plot')

parser.add_argument('--norms',
                    nargs=6,
                    help='color map norms',
                    default=[None]*6)

parser.add_argument('--colormaps',
                    nargs=6,
                    help='color maps',
                    default=['rainbow4']*6)

parser.add_argument('--zlabel',
                    default=r'$z = -\pi/2$',
                    help='z-label for location')

args = parser.parse_args()
fnames = args.filenames
steps = np.asarray(args.steps)
cmaps = args.colormaps
norms = args.norms
file_numbers = np.asarray(args.file_numbers)
plane = args.plane
loc = args.loc
save_path = args.save_path
overwrite = args.overwrite
field = args.field
fignum = args.fignum
zlabel = args.zlabel

if field is None:
    print("No field provided.")
    exit()

print()
print("\tFiles:       ", fnames)
print("\tField:       ", field)
print("\tSteps:       ", steps)
print("\tColormaps:   ", cmaps)
print("\tNorms:       ", norms)
print("\tFile numbers:", file_numbers)
print("\tPlane:       ", plane)
print("\tLocation:    ", loc)
print("\tSave path:   ", save_path)
print("\tOverwrite:   ", overwrite)
print("\tFignum:      ", fignum)
print()
    
ncreader = nc_reader()

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

for i, step in enumerate(steps):
    ncreader.open(fnames[file_numbers[i]])
    t = ncreader.get_all('t')

    fdata = ncreader.get_dataset(step=step, name=field)

    ax = grid[i]
    im, cbar = make_imshow(ax=ax,
                           plane=plane,
                           loc=loc,
                           fdata=fdata,
                           ncr=ncreader,
                           cmap=cmaps[i],
                           cmap_norm=norms[i],
                           colorbar=True)

    ncreader.close()
    
    if i < 3:
        remove_xticks(ax)

    if i == 0 or i == 3:
        pass
    else:
        remove_yticks(ax)

    add_timestamp(grid[i], t[step], xy=(0.03, 1.06), fontsize=12)

if not zlabel is None:
    add_annotation(grid[0], zlabel, xy=(-0.4, 1.25), fontsize=12)

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
