import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
from iso_surface import iso_surface
from mpl_toolkits.axes_grid1 import ImageGrid
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create iso-surface plots.')
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

args = parser.parse_args()
fnames = args.filenames
steps = np.asarray(args.steps)
file_numbers = np.asarray(args.file_numbers)
save_path = args.save_path
overwrite = args.overwrite
fields = args.fields
fignums = args.fignums

if fields is None:
    print("No fields provided.")
    exit()

if fignums is None:
    fignums = np.arange(1, len(fields) + 1, dtype=int)

if not len(fields) == len(fignums):
    print("Number of fields and number of figures do not match.")
    exit()

print()
print("\tFiles:           ", fnames)
print("\tFields:          ", fields)
print("\tSteps:           ", steps)
print("\tFile numbers:    ", file_numbers)
print("\tSave path:       ", save_path)
print("\tOverwrite:       ", overwrite)
print("\tFignums:         ", fignums)
print()

ncreader = nc_reader()

for j, field in enumerate(fields):
    print('Plotting', field)
    fig = plt.figure(figsize=(7.5, 7), dpi=500)
    grid = ImageGrid(fig, 111,
                    nrows_ncols=(3, 3),
                    aspect=True,
                    axes_pad=(0.01, 0.01),
                    direction='row',
                    share_all=True,
                    cbar_location="right",
                    cbar_mode='none',
                    cbar_size="4%",
                    cbar_pad=0.05)

    for i, step in enumerate(steps):
        ncreader.open(fnames[file_numbers[i]])
        t = ncreader.get_all('t')
        ncreader.close()
        im = plt.imread(os.path.join(save_path, 'temp_fig' + str(fignums[j]) + '_' + str(i) + '.png'))
        ax = grid[i]
        ax.imshow(im)
        ax.axis('off')
        add_timestamp(ax, t[step], xy=(0.2, 0.95), fmt="%.2f", fontsize=8)

    save_figure(plt=plt, figpath=save_path, fignum=fignums[j], overwrite=overwrite)
    plt.close()
