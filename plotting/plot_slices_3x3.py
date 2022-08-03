import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib as mpl
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
                    nargs=9,
                    default=[0, 1, 2, 3, 4, 5, 6, 7, 8],
                    help='add 9 steps to plot')

parser.add_argument('--file_numbers',
                    type=int,
                    nargs=9,
                    default=[0, 0, 0, 0, 0, 0, 0, 0, 0],
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

parser.add_argument('--fignums',
                    type=int,
                    nargs='+',
                    help='figure numbers')

parser.add_argument('--fields',
                    nargs='+',
                    help='fields to plot')

parser.add_argument('--norms',
                    nargs=9,
                    help='color map norms',
                    default=[None]*9)

parser.add_argument('--zlabel',
                    default=r'$z = -\pi/2$',
                    help='z-label for location')

args = parser.parse_args()
fnames = args.filenames
steps = np.asarray(args.steps)
file_numbers = np.asarray(args.file_numbers)
plane = args.plane
loc = args.loc
save_path = args.save_path
overwrite = args.overwrite
fields = args.fields
norms = args.norms
fignums = args.fignums
zlabel = args.zlabel


if fields is None:
    print("No fields provided.")
    exit()

if fignums is None:
    fignums = np.arange(1, len(fields) + 1, dtype=int)

if not len(fields) == len(fignums):
    print("Number of fields and number of figures do not match.")
    exit()

print()
print("\tFiles:       ", fnames)
print("\tFields:      ", fields)
print("\tSteps:       ", steps)
print("\tFile numbers:", file_numbers)
print("\tNorms:       ", norms)
print("\tPlane:       ", plane)
print("\tLocation:    ", loc)
print("\tSave path:   ", save_path)
print("\tOverwrite:   ", overwrite)
print("\tFignums:     ", fignums)
print()

ncreader = nc_reader()

for j, field in enumerate(fields):
    print('Plotting', field)
    fig = plt.figure(figsize=(8, 7.5), dpi=200)
    grid = ImageGrid(fig, 111,
                    nrows_ncols=(3, 3),
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
                               cmap='rainbow4',
                               cmap_norm=norms[i],
                               colorbar=True)

        ncreader.close()
        
        if i < 6:
            remove_xticks(ax)

        if i == 0 or i == 3 or i == 6:
            pass
        else:
            remove_yticks(ax)

        add_timestamp(ax, t[step], xy=(0.03, 1.06), fmt="%.2f")

    if not zlabel is None:
        add_annotation(grid[0], zlabel, xy=(-0.5, 1.25), fontsize=12)

    save_figure(plt=plt, figpath=save_path, fignum=fignums[j], overwrite=overwrite)
    plt.close()
