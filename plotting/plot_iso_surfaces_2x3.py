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
                    nargs=6,
                    default=[0, 1, 2, 3, 4, 5],
                    help='add 3 steps to plot')

parser.add_argument('--file_numbers',
                    type=int,
                    nargs=6,
                    default=[0, 0, 0, 0, 0, 0],
                    help='the file by index to read data from')

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
steps = np.asarray(args.steps)
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum
file_numbers = np.asarray(args.file_numbers)

print()
print("\tFiles:           ", fnames)
print("\tSteps:           ", steps)
print("\tFile numbers:    ", file_numbers)
print("\tSave path:       ", save_path)
print("\tOverwrite:       ", overwrite)
print("\tFignum:          ", fignum)
print()

ncreader = nc_reader()

fig = plt.figure(figsize=(7.5, 4.5), dpi=550)
grid = ImageGrid(fig, 111,
                 nrows_ncols=(2, 3),
                 aspect=True,
                 axes_pad=(0.1, 0.1),
                 direction='row',
                 share_all=True,
                 cbar_location="right",
                 cbar_mode='none',
                 cbar_size="4%",
                 cbar_pad=0.05)

for i in range(6):
    ncreader.open(fnames[file_numbers[i]])
    t = ncreader.get_all('t')
    ncreader.close()
    step = steps[i]

    im = plt.imread(os.path.join(save_path, 'temp_fig' + str(fignum) + '_' + str(i) + '.png'))
    ax = grid[i]
    ax.imshow(im)
    ax.axis('off')
    add_timestamp(ax, t[step], xy=(0.03, 1.0), fmt="%.2f", fontsize=8)

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()

#for i in range(6):
#    os.remove(os.path.join(save_path, 'temp_fig' + str(fignum) + '_' + str(i) + '.png'))
